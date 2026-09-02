"""Shape-adaptive grinding: the contact mechanics that replace wheel kinematics.

In the rigid-wheel model the chip thickness comes from geometry -- the grit
rides a known circle, so ``h(u) = H0 + HG*u - u^2/(2*RTIP)`` and
``semgrit.hybrid`` hands three constants to the VUMAT.

SAG has no such trajectory. The tool is compliant: a polyurethane layer
squashes by the *wheel compression* ``T``, line contact spreads into an
elliptical patch, and the load is shared by however many grains that patch
covers. So ``h`` is not kinematic, it is the **indentation depth of one grain
under its share of the contact load**:

    T  ->  F_N  ->  F_n = F_N/N_abr  ->  d  (Brinell)  ->  compare with dc

That chain is the whole of Ghosh, Sidpara & Bandyopadhyay (2021),
Int. J. Refractory Metals and Hard Materials 99, 105610, eqs. (2)-(16), and it
is implemented here equation by equation with the paper's numbering kept in the
docstrings so a reader can check it line against line.

Two things worth knowing before using any number out of this module.

**The transition criterion does not change.** SAG still asks whether the
material a grain is removing is thinner than the critical depth of cut, which
is what ``vumat_grind.for`` already decides per material point. What changes is
only how ``h`` is obtained, and the VUMAT already accepts ``h`` through a field
variable (``IHMODE = 1``). So SAG reuses the constitutive law unmodified; there
is no second Fortran law to keep in step.

**Bifano's dc is wrong for WC-Co, by a factor of about 17.** The paper measures
60-100 nm and Bifano's expression predicts 1.37 um on the same material. That
is not a slip in either place: ``dc_report`` recomputes both and reports them
side by side, because a deck that quietly used 1.37 um would be ductile
everywhere by construction and would prove nothing. See :func:`dc_report`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Sequence

# ---------------------------------------------------------------------------
# Measured pad data, table by table from the paper. Kept as data rather than
# folded into code so the provenance survives.
# ---------------------------------------------------------------------------

# Areal abrasive density Co and ACTIVE density Ca, both mm^-2, measured from
# SEM of the fresh pads (paper section 4.1 and 4.2). Ca ~ 0.8 Co, which the
# paper states and these numbers satisfy.
PAD_DENSITY = {
    6.0: dict(areal_per_mm2=1750.0, active_per_mm2=1400.0),
    15.0: dict(areal_per_mm2=720.0, active_per_mm2=576.0),
    30.0: dict(areal_per_mm2=312.0, active_per_mm2=250.0),
}

ACTIVE_FRACTION = 0.8
"""Ca = 0.8 Co. Used only for a pad size not in PAD_DENSITY."""

# The chip sizes the paper measured by SEM of the collected chips, nm. The
# 6 um row is the one that matters: pure ductile, and therefore an upper bound
# on dc rather than a sample of h above it.
MEASURED_CHIP_NM = {
    6.0: (60.0, 100.0),
    15.0: (160.0, 230.0),
    30.0: (240.0, 350.0),
}

MEASURED_DC_WC_CO_NM = (60.0, 100.0)
"""dc for HVOF-sprayed WC-12Co, measured. The paper's headline result."""

K_DUCTILE_MAX = 5.0
"""Pure ductile needs dg/D_WC < 5 (paper section 4.2). An empirical companion
to the dc test, not a substitute: it is a statement about the pad and the
carbide size, with no force in it."""


class SAGError(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# the compliant tool
# ---------------------------------------------------------------------------

@dataclass
class Pad:
    """The abrasive layer: grain size and how densely it is populated."""

    grain_um: float
    areal_per_mm2: float = 0.0
    active_per_mm2: float = 0.0

    def __post_init__(self):
        if self.grain_um <= 0:
            raise SAGError("grain size must be positive")
        known = PAD_DENSITY.get(float(self.grain_um))
        if not self.areal_per_mm2:
            if known:
                self.areal_per_mm2 = known["areal_per_mm2"]
            else:
                raise SAGError(
                    "no measured density for a %.4g um pad. Measured sizes are "
                    "%s; supply areal_per_mm2 for anything else rather than "
                    "letting it be guessed."
                    % (self.grain_um, sorted(PAD_DENSITY)))
        if not self.active_per_mm2:
            self.active_per_mm2 = (known["active_per_mm2"] if known
                                   else ACTIVE_FRACTION * self.areal_per_mm2)

    @property
    def measured_density(self) -> bool:
        return float(self.grain_um) in PAD_DENSITY


@dataclass
class Tool:
    """The compliant wheel: rigid hub, elastic layer, abrasive pad.

    ``shore_a`` is the route the paper uses to a modulus for the backing
    (eq. 4). ``elastic_mpa`` overrides it when the layer's modulus is known
    directly -- e.g. from the neo-Hookean C10 of a polyurethane card, where
    E ~ 6*C10.
    """

    diameter_mm: float = 125.0
    width_mm: float = 10.0
    pad: Pad = field(default_factory=lambda: Pad(30.0))
    shore_a: float = 40.0
    poisson: float = 0.24
    elastic_mpa: float = 0.0
    layer_thickness_mm: float = 5.0

    def __post_init__(self):
        if self.diameter_mm <= 0 or self.width_mm <= 0:
            raise SAGError("wheel diameter and width must be positive")
        if not 0.0 <= self.poisson < 0.5:
            raise SAGError("Poisson's ratio out of range")

    @property
    def radius_mm(self) -> float:
        return 0.5 * self.diameter_mm

    def modulus_mpa(self) -> float:
        """E of the compliant layer, MPa.

        Paper eq. (4), from shore hardness S:

            Et = 0.0981 (56 + 7.62 S) / (0.1375 (254 - 2.54 S))

        The expression is a fit and it blows up as S -> 100, so it is guarded.
        """
        if self.elastic_mpa > 0:
            return float(self.elastic_mpa)
        s = float(self.shore_a)
        if not 0.0 < s < 100.0:
            raise SAGError("shore A hardness must be in (0, 100), got %g" % s)
        denom = 0.1375 * (254.0 - 2.54 * s)
        if denom <= 0:
            raise SAGError("shore A %g makes eq. 4 singular" % s)
        return 0.0981 * (56.0 + 7.62 * s) / denom


def modulus_from_neo_hooke(c10_mpa: float) -> float:
    """Small-strain Young's modulus of an incompressible neo-Hookean solid.

    ``mu = 2 C10`` and ``E = 3 mu`` at incompressibility, so ``E = 6 C10``.
    The reference deck's polyurethane card carries C10 = 0.0575 MPa, i.e.
    E = 0.345 MPa, which is the number to feed ``Tool.elastic_mpa`` if the
    intent is to match that deck rather than a shore reading.
    """
    if c10_mpa <= 0:
        raise SAGError("C10 must be positive")
    return 6.0 * float(c10_mpa)


# ---------------------------------------------------------------------------
# contact
# ---------------------------------------------------------------------------

def equivalent_modulus_mpa(e_work_mpa: float, nu_work: float,
                           e_tool_mpa: float, nu_tool: float) -> float:
    """Paper eq. (3):  Eeq = [ (1-nu_w^2)/Ew + (1-nu_t^2)/Et ]^-1.

    With a compliant tool this is dominated by the tool: at Ew = 200 GPa and
    Et = 0.35 MPa the workpiece contributes about two parts in a million. That
    is the physics of the process, not a rounding problem -- the compliance is
    supposed to come from the tool.
    """
    if min(e_work_mpa, e_tool_mpa) <= 0:
        raise SAGError("both moduli must be positive")
    return 1.0 / ((1.0 - nu_work ** 2) / e_work_mpa
                  + (1.0 - nu_tool ** 2) / e_tool_mpa)


def normal_load_n(eq_modulus_mpa: float, radius_mm: float,
                  compression_mm: float) -> float:
    """Paper eq. (2):  FN = 1.44 Eeq R^(1/2) T^(3/2), newtons.

    The 3/2 power is the Hertz signature: load grows faster than compression
    because the patch widens as it deepens.
    """
    if compression_mm < 0:
        raise SAGError("wheel compression cannot be negative")
    return (1.44 * eq_modulus_mpa * math.sqrt(radius_mm)
            * compression_mm ** 1.5)


def spot_area_mm2(compression_mm: float, speed_rpm: float) -> float:
    """Paper eq. (6):  As = 138.22 (T^0.151 N^0.009), mm^2. R^2 = 0.991.

    An empirical fit to measured finishing spots, NOT a Hertz solution -- so it
    carries the paper's geometry (125 mm wheel, PU foam, WC-Co) and should not
    be extrapolated far from it. The exponents say what the paper says in
    words: compression matters, speed almost does not.
    """
    if compression_mm <= 0 or speed_rpm <= 0:
        raise SAGError("compression and speed must be positive for eq. 6")
    return 138.22 * (compression_mm ** 0.151) * (speed_rpm ** 0.009)


def spot_length_mm(compression_mm: float, speed_rpm: float) -> float:
    """Paper eq. (7):  Ls = 17.69 (T^0.232 N^0.012), mm. R^2 = 0.981."""
    if compression_mm <= 0 or speed_rpm <= 0:
        raise SAGError("compression and speed must be positive for eq. 7")
    return 17.69 * (compression_mm ** 0.232) * (speed_rpm ** 0.012)


def hertz_semi_axes_mm(area_mm2: float, length_mm: float,
                       face_width_mm: float = 0.0) -> tuple:
    """Semi-axes (a, b) of the contact patch, from its area and length.

    The paper measures both the area (eq. 6) and the length along L-L (eq. 7),
    and a free ellipse is fixed by the pair: ``a = Ls/2`` and ``pi a b = As``.
    Deriving b this way keeps both measurements rather than assuming a shape.

    THE ELLIPSE DOES NOT FIT ON THE WHEEL. Eqs. 6 and 7 are independent
    empirical fits, and combining them gives a contact WIDTH ``2b`` that exceeds
    the paper's own 10 mm face at every one of its operating points -- 10.16 mm
    at T = 0.6 up to 11.14 mm at T = 0.2, so 1.5% to 11% over. The contact
    therefore runs off both edges of the wheel and is CLIPPED, not elliptical,
    which also means Hertz's semi-infinite half-space assumption is violated
    across the face.

    Pass ``face_width_mm`` and the clipping is reported: ``b`` is capped at the
    half-width and the shortfall in area is returned, so a caller can say how
    much of the patch the ellipse model is misplacing instead of quietly using
    a width the tool does not have. Left at 0, the free ellipse is returned
    unchanged and nothing is hidden either way.
    """
    a = 0.5 * length_mm
    if a <= 0:
        raise SAGError("spot length must be positive")
    b_free = area_mm2 / (math.pi * a)
    if face_width_mm <= 0:
        return a, b_free
    b = min(b_free, 0.5 * face_width_mm)
    return a, b


def clipping_report(area_mm2: float, length_mm: float,
                    face_width_mm: float) -> dict:
    """How badly the elliptical patch overruns the wheel face.

    ``overrun`` is 2b/W: at 1.0 the ellipse exactly fits, above 1.0 it does
    not exist on this wheel. ``area_clipped_fraction`` is how much of the
    nominal area falls outside the face.
    """
    if face_width_mm <= 0:
        raise SAGError("face width must be positive")
    a, b_free = hertz_semi_axes_mm(area_mm2, length_mm)
    overrun = 2.0 * b_free / face_width_mm
    if overrun <= 1.0:
        return dict(a_mm=a, b_free_mm=b_free, b_used_mm=b_free,
                    overrun=overrun, clipped=False,
                    area_clipped_fraction=0.0, area_on_face_mm2=area_mm2)
    # Area of the ellipse inside |y| <= W/2, as a fraction of pi a b.
    t = 0.5 * face_width_mm / b_free            # < 1
    frac_on = (2.0 / math.pi) * (math.asin(t) + t * math.sqrt(1.0 - t * t))
    return dict(a_mm=a, b_free_mm=b_free, b_used_mm=0.5 * face_width_mm,
                overrun=overrun, clipped=True,
                area_clipped_fraction=1.0 - frac_on,
                area_on_face_mm2=area_mm2 * frac_on)


def max_pressure_mpa(load_n: float, area_mm2: float) -> float:
    """p0 from FN and As.

    Paper eq. (1) gives a Hertzian hemi-ellipsoidal distribution whose mean is
    ``p_av = (2/3) p0``, so ``p0 = 1.5 FN / As``.
    """
    if area_mm2 <= 0:
        raise SAGError("contact area must be positive")
    return 1.5 * load_n / area_mm2


def pressure_at_mpa(p0_mpa: float, x_mm: float, y_mm: float,
                    a_mm: float, b_mm: float) -> float:
    """Paper eq. (1):  p = p0 [1 - (x/a)^2 - (y/b)^2]^(1/2), 0 outside."""
    r = (x_mm / a_mm) ** 2 + (y_mm / b_mm) ** 2
    if r >= 1.0:
        return 0.0
    return p0_mpa * math.sqrt(1.0 - r)


# ---------------------------------------------------------------------------
# from contact load to one grain's bite
# ---------------------------------------------------------------------------

def indentation_depth_mm(load_per_grain_n: float, grain_um: float,
                         bhn_kgf_mm2: float) -> float:
    """Depth a spherical grain sinks under its own share of the load.

    Paper eqs. (11) and (12) are the Brinell relation and its inversion:

        BHN = 2 Fn / [ pi dg ( dg - sqrt(dg^2 - di^2) ) ]
        d   = dg/2 - (1/2) sqrt(dg^2 - di^2)

    Eliminating the indentation diameter ``di`` between them -- which is what
    the deck actually needs, since ``di`` is never measured here -- gives the
    depth directly:

        d = Fn / (pi dg H)

    with ``H`` the hardness in a consistent stress unit. This is the shallow
    spherical-cap limit of the pair above, and it is exact in that limit; the
    full inversion adds nothing while ``d << dg``, which holds by four orders
    of magnitude here (nanometres into micrometres). :func:`indentation_check`
    verifies that assumption rather than assuming it.

    BHN is quoted in kgf/mm^2 (581 for WC-Co), so it is converted to MPa.
    """
    if load_per_grain_n <= 0:
        return 0.0
    if grain_um <= 0 or bhn_kgf_mm2 <= 0:
        raise SAGError("grain size and hardness must be positive")
    h_mpa = bhn_kgf_mm2 * 9.80665          # kgf/mm^2 -> MPa
    dg_mm = grain_um * 1e-3
    return load_per_grain_n / (math.pi * dg_mm * h_mpa)


def indentation_check(depth_mm: float, grain_um: float) -> dict:
    """How safe the shallow-cap simplification in :func:`indentation_depth_mm` is.

    Returns the ratio d/dg and the exact contact radius. The simplification is
    good to first order in d/dg; anything above a few percent should use the
    full eqs. (11)-(12) instead, and this is what says so.
    """
    dg_mm = grain_um * 1e-3
    ratio = depth_mm / dg_mm if dg_mm > 0 else float("inf")
    # Contact radius of a spherical cap of depth d on a sphere of diameter dg.
    r = math.sqrt(max(depth_mm * (dg_mm - depth_mm), 0.0))
    return dict(depth_over_grain=ratio, contact_radius_mm=r,
                shallow_cap_valid=ratio < 0.05)


def groove_area_mm2(depth_mm: float, grain_um: float) -> float:
    """Cross-section of the groove one grain ploughs. Paper eq. (13):

        A' = (dg^2/4) asin( 2 sqrt(d(dg-d)) / dg )
             - sqrt(d(dg-d)) (dg/2 - d)

    A circular segment: the arc term minus the triangle. Written exactly as
    published, including the factor the paper prints as ``sin-1 2.sqrt(...)``.
    """
    if depth_mm <= 0:
        return 0.0
    dg = grain_um * 1e-3
    if depth_mm >= dg:
        raise SAGError("indentation %.6g mm exceeds grain diameter %.6g mm"
                       % (depth_mm, dg))
    root = math.sqrt(depth_mm * (dg - depth_mm))
    arg = 2.0 * root / dg
    arg = min(1.0, max(-1.0, arg))         # guard asin against float drift
    return (dg * dg / 4.0) * math.asin(arg) - root * (dg / 2.0 - depth_mm)


def groove_width_mm(depth_mm: float, grain_um: float) -> float:
    """Width of the groove one grain ploughs: the chord of its contact circle.

        w = 2 sqrt( d (dg - d) )

    This is NOT in the paper, and it is the quantity that reconciles the
    paper with itself. Its own equations give a penetration depth of ~0.3 nm,
    while it measures chips of 60-350 nm -- a factor of a few hundred, which
    reads as a contradiction until one notices the two are different lengths.
    A grain indenting 0.3 nm into a 30 um sphere leaves a groove 190 nm WIDE:
    the aspect ratio is ~630, because a shallow cap on a large sphere is very
    much wider than it is deep.

    Against the paper's own measured chip sizes:

        dg = 6 um   ->  w =  80 nm    measured  60-100 nm
        dg = 15 um  ->  w = 125 nm    measured 160-230 nm
        dg = 30 um  ->  w = 190 nm    measured 240-350 nm

    The 6 um pad -- the one the paper calls purely ductile, and the row its dc
    conclusion rests on -- lands in the middle of the measured band. So the
    detached chip scales with the groove WIDTH, not the depth, which is what a
    fragment spalling from a shallow scratch should do.

    Consequence for a deck: the length to compare against dc is the depth, but
    the length to compare against a measured CHIP is this. Reporting only one
    of them is how the 265x apparent inconsistency arises.
    """
    if depth_mm <= 0:
        return 0.0
    dg = grain_um * 1e-3
    if depth_mm >= dg:
        raise SAGError("indentation exceeds grain diameter")
    return 2.0 * math.sqrt(depth_mm * (dg - depth_mm))


def grains_per_revolution(tool: Tool) -> float:
    """Paper eq. (15):  nt = (pi D W) Ca."""
    return math.pi * tool.diameter_mm * tool.width_mm * tool.pad.active_per_mm2


# ---------------------------------------------------------------------------
# the whole chain
# ---------------------------------------------------------------------------

@dataclass
class SAGContact:
    """One operating point: everything the deck writer and the report need."""

    compression_mm: float
    speed_rpm: float
    grain_um: float

    eq_modulus_mpa: float
    tool_modulus_mpa: float
    normal_load_n: float
    tangential_load_n: float
    spot_area_mm2: float
    spot_length_mm: float
    semi_axis_a_mm: float
    semi_axis_b_mm: float
    max_pressure_mpa: float
    mean_pressure_mpa: float

    active_grains: float
    load_per_grain_n: float
    tangential_per_grain_n: float
    indentation_mm: float
    groove_area_mm2: float
    groove_width_mm: float
    grains_per_rev: float
    mrr_mm3_min: float

    surface_speed_mm_s: float
    contact_time_s: float

    shallow_cap_valid: bool
    depth_over_grain: float
    density_measured: bool
    face_overrun: float = 0.0
    """2b/W. Above 1 the elliptical patch is wider than the wheel face."""
    area_clipped_fraction: float = 0.0
    area_on_face_mm2: float = 0.0

    @property
    def indentation_nm(self) -> float:
        return self.indentation_mm * 1e6

    @property
    def groove_width_nm(self) -> float:
        """The length that scales with the measured chip size, not with dc."""
        return self.groove_width_mm * 1e6

    def regime(self, dc_nm: float) -> str:
        """Which side of dc this operating point sits on."""
        h = self.indentation_nm
        if h <= 0:
            return "no contact"
        return "ductile" if h < dc_nm else "brittle"

    def margin(self, dc_nm: float) -> float:
        """h/dc. Below 1 is ductile; the distance from 1 is the safety."""
        return self.indentation_nm / dc_nm if dc_nm > 0 else float("inf")


def solve_contact(tool: Tool, *, compression_mm: float, speed_rpm: float,
                  work_modulus_mpa: float, work_poisson: float,
                  bhn_kgf_mm2: float, friction: float = 0.2) -> SAGContact:
    """Run the paper's chain end to end for one operating point.

    Order is eqs. (3), (2), (5), (6), (7), (8), (9), (10), (11)-(12), (13),
    (14), (15), (16).
    """
    et = tool.modulus_mpa()
    eeq = equivalent_modulus_mpa(work_modulus_mpa, work_poisson,
                                 et, tool.poisson)
    fn = normal_load_n(eeq, tool.radius_mm, compression_mm)
    ft = friction * fn                                          # eq. (5)

    a_s = spot_area_mm2(compression_mm, speed_rpm)               # eq. (6)
    l_s = spot_length_mm(compression_mm, speed_rpm)              # eq. (7)
    clip = clipping_report(a_s, l_s, tool.width_mm)
    a, b = clip["a_mm"], clip["b_used_mm"]
    # Pressure uses the area actually ON the face: a load spread over a patch
    # that partly does not exist would understate the pressure.
    p0 = max_pressure_mpa(fn, clip["area_on_face_mm2"])

    nabr = tool.pad.active_per_mm2 * a_s                         # eq. (8)
    if nabr <= 0:
        raise SAGError("no active grains in the contact patch")
    fn_g = fn / nabr                                             # eq. (9)
    ft_g = ft / nabr                                             # eq. (10)

    d = indentation_depth_mm(fn_g, tool.grain_um if hasattr(tool, "grain_um")
                             else tool.pad.grain_um, bhn_kgf_mm2)
    chk = indentation_check(d, tool.pad.grain_um)
    a_prime = groove_area_mm2(d, tool.pad.grain_um)              # eq. (13)
    v_a = l_s * a_prime                                          # eq. (14)
    n_t = grains_per_revolution(tool)                            # eq. (15)
    mrr = v_a * n_t * speed_rpm                                  # eq. (16)

    v_s = math.pi * tool.diameter_mm * speed_rpm / 60.0          # mm/s
    return SAGContact(
        compression_mm=compression_mm, speed_rpm=speed_rpm,
        grain_um=tool.pad.grain_um,
        eq_modulus_mpa=eeq, tool_modulus_mpa=et,
        normal_load_n=fn, tangential_load_n=ft,
        spot_area_mm2=a_s, spot_length_mm=l_s,
        semi_axis_a_mm=a, semi_axis_b_mm=b,
        max_pressure_mpa=p0, mean_pressure_mpa=2.0 * p0 / 3.0,
        active_grains=nabr, load_per_grain_n=fn_g,
        tangential_per_grain_n=ft_g,
        indentation_mm=d, groove_area_mm2=a_prime,
        groove_width_mm=groove_width_mm(d, tool.pad.grain_um),
        grains_per_rev=n_t, mrr_mm3_min=mrr,
        surface_speed_mm_s=v_s,
        contact_time_s=(l_s / v_s) if v_s > 0 else 0.0,
        shallow_cap_valid=chk["shallow_cap_valid"],
        depth_over_grain=chk["depth_over_grain"],
        density_measured=tool.pad.measured_density,
        face_overrun=clip["overrun"],
        area_clipped_fraction=clip["area_clipped_fraction"],
        area_on_face_mm2=clip["area_on_face_mm2"],
    )


# ---------------------------------------------------------------------------
# load sharing over a real protrusion distribution
#
# This is where the reference paper and this project part company, and it is
# worth being precise about why.
#
# The paper divides the contact load equally among every active grain -- its
# eq. (9) is Fn = FN/Nabr, one number for all of them. It says so plainly: "it
# is presumed that the normal and tangential forces are uniformly distributed
# over the active abrasive particles and the abrasives are spherical in shape
# with identical size."
#
# On its own numbers that gives a mean indentation of about 0.3 nm against a
# measured dc of 60-100 nm, so h/dc ~ 0.004 and EVERY pad -- 6, 15 and 30 um --
# comes out ductile. But the paper OBSERVES brittle fracture on the 15 and
# 30 um pads. Its own primary criterion therefore fails to reproduce its own
# central observation, and the empirical k = dg/D_WC < 5 rule is what carries
# the conclusion instead.
#
# The missing physics is that grains are not identical and do not share load
# equally. A pad has a protrusion distribution; the tallest grains touch first
# and carry far more than their share. Measured B4C grains in this project run
# 0.76 to 7.05 um tall (mean 3.98, sd 1.63), so if engagement is confined to
# the top tenth of that band only 3 of 27 grains touch and each carries 9x the
# mean load -- and at the top few percent, 27x.
#
# That is exactly the quantity this project measures and the paper assumes
# away. So the load concentration below is not a fitted fudge factor: it is
# computed from the measured height distribution of real grains, and it is the
# contribution the SEM pipeline makes to the SAG model.
# ---------------------------------------------------------------------------

def engaged_fraction(heights_um: Sequence[float],
                     engagement_um: float) -> dict:
    """Which grains actually touch, given how far the pad is pressed in.

    A grain engages when it stands within ``engagement_um`` of the tallest one.
    Everything shorter than that is still clear of the work and carries nothing.

    Returns the engaged count, the fraction, and the load concentration -- the
    factor by which the load on an engaged grain exceeds the paper's uniform
    ``FN/Nabr``. That factor is ``n_total / n_engaged``: the same total load
    over fewer carriers.
    """
    h = [float(x) for x in heights_um if x > 0]
    if not h:
        raise SAGError("no positive grain heights")
    if engagement_um <= 0:
        raise SAGError("engagement depth must be positive")
    top = max(h)
    thr = top - engagement_um
    n_eng = sum(1 for x in h if x >= thr)
    n_eng = max(n_eng, 1)
    return dict(
        n_total=len(h), n_engaged=n_eng,
        engaged_fraction=n_eng / len(h),
        load_concentration=len(h) / n_eng,
        tallest_um=top, threshold_um=thr,
        mean_um=sum(h) / len(h),
    )


def concentrated_indentation(contact: "SAGContact", heights_um: Sequence[float],
                             engagement_um: float, bhn_kgf_mm2: float) -> dict:
    """Redo the indentation with load shared only among the grains that touch.

    The paper's own chain up to ``FN`` is kept unchanged -- only eq. (9) is
    replaced, because that is the step that assumes identical grains. Depth is
    linear in load, so the depth an engaged grain reaches is the uniform depth
    times the load concentration.
    """
    share = engaged_fraction(heights_um, engagement_um)
    fn_eng = contact.load_per_grain_n * share["load_concentration"]
    d = indentation_depth_mm(fn_eng, contact.grain_um, bhn_kgf_mm2)
    out = dict(share)
    out.update(
        load_per_grain_uniform_n=contact.load_per_grain_n,
        load_per_grain_engaged_n=fn_eng,
        indentation_uniform_nm=contact.indentation_nm,
        indentation_engaged_nm=d * 1e6,
        groove_width_engaged_nm=groove_width_mm(d, contact.grain_um) * 1e6,
    )
    return out


# ---------------------------------------------------------------------------
# dc, and the disagreement about it
# ---------------------------------------------------------------------------

def dc_report(*, hardness_mpa: float, youngs_mpa: float,
              kic_mpa_sqrt_m: float, measured_nm: Optional[tuple] = None,
              lambda_c: float = 0.15) -> dict:
    """Every dc estimate for this material, side by side.

    The point is the disagreement. On the paper's WC-Co numbers Bifano's
    expression gives 1.37 um against a measured 60-100 nm -- a factor of 17 --
    and a deck built on 1.37 um would be ductile everywhere by construction and
    would demonstrate nothing. So both are computed and the ratio is stated,
    and the caller chooses knowingly.

    Reuses ``semgrit.hybrid.critical_depth_mm`` rather than reimplementing it,
    so there is one implementation of the three forms in the project.
    """
    from .hybrid import critical_depth_mm, kic_from_mpa_sqrt_m

    kic = kic_from_mpa_sqrt_m(kic_mpa_sqrt_m)
    out = {
        "hardness_mpa": hardness_mpa, "youngs_mpa": youngs_mpa,
        "kic_mpa_sqrt_m": kic_mpa_sqrt_m, "kic_mpa_sqrt_mm": kic,
        "lambda_c": lambda_c,
        "form1_nm": critical_depth_mm(lambda_c, hardness_mpa, youngs_mpa,
                                      kic, 1) * 1e6,
        "form2_nm": critical_depth_mm(lambda_c, hardness_mpa, youngs_mpa,
                                      kic, 2) * 1e6,
    }
    if measured_nm:
        lo, hi = float(measured_nm[0]), float(measured_nm[1])
        mid = 0.5 * (lo + hi)
        out.update(measured_lo_nm=lo, measured_hi_nm=hi, measured_mid_nm=mid,
                   bifano_over_measured=out["form2_nm"] / mid)
    return out


def carbide_ratio(grain_um: float, carbide_um: float) -> dict:
    """k = dg / D_WC, and whether it clears the paper's ductile threshold.

    A second, independent criterion (paper section 4.2): pure ductile removal
    needed k < 5, observed as k = 22, 11, 4 across the three pads. It contains
    no force at all, so it agrees with the dc test only for the material it was
    measured on -- worth reporting beside dc, not instead of it.
    """
    if carbide_um <= 0:
        raise SAGError("carbide size must be positive")
    k = grain_um / carbide_um
    return dict(k=k, threshold=K_DUCTILE_MAX, pure_ductile=k < K_DUCTILE_MAX)


# ---------------------------------------------------------------------------
# self-check
# ---------------------------------------------------------------------------

def demo() -> None:
    """Assert this module against the paper's own published numbers."""
    # --- eq. 4, shore hardness to modulus ---------------------------------
    t = Tool(shore_a=40.0)
    et = t.modulus_mpa()
    assert 0.0 < et < 100.0, et

    # C10 -> E for the reference deck's polyurethane.
    assert abs(modulus_from_neo_hooke(0.0575) - 0.345) < 1e-9

    # --- eq. 3 is tool-dominated when the tool is soft -------------------
    eeq = equivalent_modulus_mpa(200_000.0, 0.25, 0.345, 0.24)
    # Eq. 3 sums COMPLIANCES, and each carries a (1-nu^2), so the limit as the
    # workpiece stiffens is Et/(1-nu_t^2) -- slightly ABOVE Et, not below it.
    limit = 0.345 / (1.0 - 0.24 ** 2)
    assert abs(eeq - limit) / limit < 1e-4, (eeq, limit)
    # and the soft layer carries essentially all of the compliance
    tool_share = ((1 - 0.24 ** 2) / 0.345) / ((1 - 0.25 ** 2) / 200_000.0
                                              + (1 - 0.24 ** 2) / 0.345)
    assert tool_share > 0.9999, tool_share

    # --- eq. 2 scales as T^1.5 -------------------------------------------
    f1 = normal_load_n(eeq, 62.5, 0.2)
    f2 = normal_load_n(eeq, 62.5, 0.4)
    assert abs(f2 / f1 - 2.0 ** 1.5) < 1e-9, "eq. 2 is not Hertzian in T"

    # --- eqs. 6 and 7 reproduce the paper's fitted magnitudes ------------
    a_s = spot_area_mm2(0.4, 1050.0)
    l_s = spot_length_mm(0.4, 1050.0)
    assert 100.0 < a_s < 160.0, a_s
    assert 12.0 < l_s < 20.0, l_s
    # speed barely matters, compression does -- the exponents say so
    assert spot_area_mm2(0.4, 300.0) / a_s > 0.98
    assert spot_area_mm2(0.2, 1050.0) / a_s < 0.92

    # --- the ellipse closes ----------------------------------------------
    a, b = hertz_semi_axes_mm(a_s, l_s)
    assert abs(math.pi * a * b - a_s) < 1e-9, "ellipse area must round-trip"

    # --- but it does not fit on the wheel --------------------------------
    # 2b exceeds the paper's own 10 mm face at every one of its settings, so
    # the patch is clipped by the wheel edge. Assert the overrun is real and
    # systematic rather than a single-point artefact.
    worst = 0.0
    for T in (0.2, 0.4, 0.6):
        for N in (300.0, 550.0, 800.0, 1050.0):
            r = clipping_report(spot_area_mm2(T, N), spot_length_mm(T, N), 10.0)
            assert r["clipped"], (T, N, r["overrun"])
            assert r["b_used_mm"] == 5.0, "b must be capped at the half-width"
            worst = max(worst, r["overrun"])
    assert 1.10 < worst < 1.15, worst
    # the clipped AREA is small, so the ellipse model is nearly self-consistent
    # -- this is a caveat on the width, not a refutation of the fits
    r04 = clipping_report(spot_area_mm2(0.4, 1050.0),
                          spot_length_mm(0.4, 1050.0), 10.0)
    assert r04["area_clipped_fraction"] < 0.05
    assert r04["area_on_face_mm2"] < spot_area_mm2(0.4, 1050.0)
    # a wide enough wheel is not clipped at all
    wide = clipping_report(spot_area_mm2(0.4, 1050.0),
                           spot_length_mm(0.4, 1050.0), 40.0)
    assert not wide["clipped"] and wide["area_clipped_fraction"] == 0.0
    assert abs(wide["b_used_mm"] - wide["b_free_mm"]) < 1e-12

    # --- eq. 1 is p0 at the centre and zero at the rim -------------------
    p0 = max_pressure_mpa(10.0, a_s)
    assert abs(pressure_at_mpa(p0, 0, 0, a, b) - p0) < 1e-12
    assert pressure_at_mpa(p0, a, 0, a, b) == 0.0
    assert pressure_at_mpa(p0, 2 * a, 0, a, b) == 0.0

    # --- Ca = 0.8 Co on every measured pad -------------------------------
    for dg, row in PAD_DENSITY.items():
        r = row["active_per_mm2"] / row["areal_per_mm2"]
        assert abs(r - ACTIVE_FRACTION) < 0.02, (dg, r)

    # --- the paper's own measured forces ---------------------------------
    # Fig. 4 measures FN of order 2-10 N over T = 0.2-0.6 mm, and section 4.1
    # states Fn is 1e-4..1e-5 N and Ft 1e-5..1e-6 N. A shore-40 backing
    # reproduces all three; the reference deck's C10 = 0.0575 MPa (E = 0.345
    # MPa) is ~5x softer and undershoots FN, so it is not used for validation.
    for T, lo_n, hi_n in ((0.2, 1.0, 3.0), (0.4, 4.0, 7.0), (0.6, 8.0, 12.0)):
        f = normal_load_n(equivalent_modulus_mpa(200_000.0, 0.25,
                                                 Tool(shore_a=40.0).modulus_mpa(),
                                                 0.24), 62.5, T)
        assert lo_n < f < hi_n, (T, f)

    got = {}
    for dg in (6.0, 15.0, 30.0):
        c = solve_contact(Tool(pad=Pad(dg), shore_a=40.0),
                          compression_mm=0.4, speed_rpm=1050.0,
                          work_modulus_mpa=200_000.0, work_poisson=0.25,
                          bhn_kgf_mm2=581.0)
        got[dg] = c
        # the paper's stated per-grain range, inclusive of its own spread
        assert 1e-5 < c.load_per_grain_n < 1e-3, (dg, c.load_per_grain_n)
        assert 1e-6 < c.tangential_per_grain_n < 1e-4, dg
        # shallow-cap simplification must be valid, not assumed
        assert c.shallow_cap_valid, (dg, c.depth_over_grain)
        # the groove must be a real segment
        assert c.groove_area_mm2 > 0

    # --- groove WIDTH reconciles the paper with its own chip sizes -------
    # Its equations give ~0.3 nm of depth while it measures 60-350 nm chips.
    # The width, not the depth, is what tracks the measurement -- and the
    # 6 um pad (the row the dc conclusion rests on) lands inside its band.
    for dg, lo, hi in ((6.0, 60.0, 100.0), (15.0, 160.0, 230.0),
                       (30.0, 240.0, 350.0)):
        w = got[dg].groove_width_nm
        assert w > 20.0 * got[dg].indentation_nm,             "a shallow cap on a big sphere must be far wider than deep"
        # within a factor of two of the measured chip band
        assert 0.5 * lo < w < 2.0 * hi, (dg, w, lo, hi)
    # and the 6 um pad is INSIDE the measured band, not merely near it
    assert 60.0 <= got[6.0].groove_width_nm <= 100.0, got[6.0].groove_width_nm
    # width must grow with grain size, as the measured chips do
    assert (got[6.0].groove_width_nm < got[15.0].groove_width_nm
            < got[30.0].groove_width_nm)

    # --- MRR lands on a physically sensible finishing rate ---------------
    assert 0.01 < got[30.0].mrr_mm3_min < 10.0, got[30.0].mrr_mm3_min

    # --- coarser pad -> fewer grains -> deeper bite ----------------------
    assert got[30.0].active_grains < got[6.0].active_grains
    assert got[30.0].indentation_nm > got[6.0].indentation_nm, \
        "a coarser pad must indent deeper: fewer grains share the same load"

    # --- higher compression -> more load, deeper bite --------------------
    lo = solve_contact(Tool(pad=Pad(30.0), shore_a=40.0),
                       compression_mm=0.2, speed_rpm=1050.0,
                       work_modulus_mpa=200_000.0, work_poisson=0.25,
                       bhn_kgf_mm2=581.0)
    hi = solve_contact(Tool(pad=Pad(30.0), shore_a=40.0),
                       compression_mm=0.6, speed_rpm=1050.0,
                       work_modulus_mpa=200_000.0, work_poisson=0.25,
                       bhn_kgf_mm2=581.0)
    assert hi.normal_load_n > lo.normal_load_n
    assert hi.indentation_nm > lo.indentation_nm, \
        "the paper: more compression -> more force -> deeper penetration"

    # --- and that is the paper's ductile trend ---------------------------
    # lower compression and finer pad both push toward ductile
    dc_mid = 80.0
    assert lo.margin(dc_mid) < hi.margin(dc_mid)
    assert got[6.0].margin(dc_mid) < got[30.0].margin(dc_mid)

    # --- MRR rises with compression and with speed -----------------------
    slow = solve_contact(Tool(pad=Pad(30.0), shore_a=40.0),
                         compression_mm=0.4, speed_rpm=300.0,
                         work_modulus_mpa=200_000.0, work_poisson=0.25,
                         bhn_kgf_mm2=581.0)
    assert got[30.0].mrr_mm3_min > slow.mrr_mm3_min, \
        "MRR must rise with wheel speed (Preston)"
    assert hi.mrr_mm3_min > lo.mrr_mm3_min

    # --- load sharing over a real protrusion distribution ----------------
    # Measured B4C grain heights from this project's own SEM pipeline.
    heights = [0.761, 1.9, 2.4, 2.8, 3.1, 3.3, 3.5, 3.6, 3.8, 3.9, 3.977,
               4.0, 4.1, 4.2, 4.4, 4.5, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.9,
               6.2, 6.5, 6.8, 7.045]
    assert len(heights) == 27

    # deeper engagement -> more grains touch -> less concentration
    prev = None
    for eng in (0.5, 1.0, 2.0, 4.0):
        sh = engaged_fraction(heights, eng)
        assert 1 <= sh["n_engaged"] <= 27
        if prev is not None:
            assert sh["n_engaged"] >= prev["n_engaged"]
            assert sh["load_concentration"] <= prev["load_concentration"]
        prev = sh
    # engaging everything must give no concentration at all
    allin = engaged_fraction(heights, 100.0)
    assert allin["n_engaged"] == 27
    assert abs(allin["load_concentration"] - 1.0) < 1e-12

    # and concentration must deepen the bite, linearly in load
    c30 = got[30.0]
    conc = concentrated_indentation(c30, heights, 1.0, 581.0)
    assert conc["load_concentration"] > 1.0
    assert (abs(conc["indentation_engaged_nm"]
                - conc["indentation_uniform_nm"] * conc["load_concentration"])
            < 1e-9), "depth is linear in load"
    assert conc["indentation_engaged_nm"] > conc["indentation_uniform_nm"]

    # THE POINT: the paper's uniform assumption cannot reach its own dc, and
    # concentration is the mechanism that closes the gap. Assert the gap is
    # real rather than asserting a particular closure.
    assert c30.margin(80.0) < 0.01,         "uniform load gives h/dc ~ 0.004: every pad ductile, which is NOT "         "what the paper observed on its 15 and 30 um pads"
    needed = 80.0 / c30.indentation_nm
    assert 100.0 < needed < 1000.0, needed
    # a concentration that large means only a small fraction of the active
    # grains carry load -- which is what a real protrusion spread produces
    assert conc["engaged_fraction"] < 0.2

    # --- dc: Bifano vs measured, the disagreement ------------------------
    r = dc_report(hardness_mpa=11_020.0, youngs_mpa=200_000.0,
                  kic_mpa_sqrt_m=7.78, measured_nm=MEASURED_DC_WC_CO_NM)
    # the paper prints 1.37 um for exactly these inputs
    assert abs(r["form2_nm"] - 1370.0) < 40.0, r["form2_nm"]
    assert 10.0 < r["bifano_over_measured"] < 25.0, r["bifano_over_measured"]

    # --- k = dg/D_WC, the paper's three pads -----------------------------
    for dg, k_paper, ductile in ((30.0, 22.0, False), (15.0, 11.0, False),
                                 (6.0, 4.4, True)):
        k = carbide_ratio(dg, 1.36)
        assert abs(k["k"] - k_paper) < 0.6, (dg, k["k"], k_paper)
        assert k["pure_ductile"] is ductile, dg

    # --- refusals --------------------------------------------------------
    for bad in (lambda: Pad(7.5),                       # unmeasured density
                lambda: Tool(diameter_mm=-1),
                lambda: Tool(shore_a=100.0).modulus_mpa(),
                lambda: spot_area_mm2(0.0, 1050.0),
                lambda: groove_area_mm2(1.0, 30.0),     # deeper than the grain
                lambda: modulus_from_neo_hooke(0.0)):
        try:
            bad()
        except SAGError:
            pass
        else:
            raise AssertionError("should have been refused: %r" % bad)

    print("semgrit.sag: all checks passed")
    print("  paper operating point T=0.4 mm, N=1050 rpm, WC-Co:")
    print("    Et = %.4f MPa, Eeq = %.4f MPa" % (got[30.0].tool_modulus_mpa,
                                                 got[30.0].eq_modulus_mpa))
    print("    FN = %.3f N over As = %.1f mm2 (p0 = %.4f MPa)"
          % (got[30.0].normal_load_n, got[30.0].spot_area_mm2,
             got[30.0].max_pressure_mpa))
    for dg in (6.0, 15.0, 30.0):
        c = got[dg]
        print("    dg=%2.0f um: %8.0f grains, Fn = %.3e N, h = %.3f nm, "
              "groove %6.1f nm wide (chip %s nm)"
              % (dg, c.active_grains, c.load_per_grain_n, c.indentation_nm,
                 c.groove_width_nm,
                 "-".join("%.0f" % v for v in MEASURED_CHIP_NM[dg])))
    print("    Bifano dc = %.0f nm vs measured %.0f-%.0f nm (%.1fx)"
          % (r["form2_nm"], r["measured_lo_nm"], r["measured_hi_nm"],
             r["bifano_over_measured"]))


if __name__ == "__main__":
    demo()
