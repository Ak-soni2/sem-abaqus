"""The ductile/brittle switch that ``vumat_grind.for`` runs on.

The subroutine decides, per material point, whether to follow Johnson-Cook
with strain-gradient enhancement or Johnson-Holmquist II. It makes that
decision from two numbers it cannot work out for itself:

``dc``
    the critical depth of cut for the ductile-brittle transition, a material
    property;
``h(u)``
    the undeformed chip thickness the grit takes at the material point's
    station along the scratch, a kinematic property of this particular deck.

Both are computed here, in Python, where they can be checked against the
placement the writer actually used, and handed over as three numbers in the
material card (H0, HG, RTIP). The Fortran therefore carries no knowledge of
wheels, infeeds or rotation senses -- which is what keeps it verifiable on a
single material point.

Why h is a straight line plus a parabola
----------------------------------------
With ONE grit the trajectory is exact. The tip rides a circle of radius
``r_tip`` about an axis that translates outward at the infeed speed ``v_r``,
and it crosses the block at surface speed ``v_s``. Writing ``u`` for the
tangential station of a point:

    time the grit arrives there   t(u) = (u0 - u) / v_s
    tip reach along e_r at t      r_tip*cos(u/r_tip) + v_r*t
    chip thickness                h(u) = that, minus the ground radius

    h(u) = H0 + HG*u - u^2/(2*r_tip),      HG = -v_r/v_s

The linear term is the wedge in every textbook figure of a grit trajectory --
rubbing, then ploughing, then shearing -- and here it is produced by the
radial infeed rather than by a table feed. The quadratic term is the sagitta
of the grit's circular path: 15 nm over a 48 um block on a 50 mm wheel, which
would be ignorable except that dc is of that same order.

The classical traverse-grinding form ``h(theta) = L_g (v_w/v_s) sin(theta)``
is the same straight line over a block much shorter than the contact arc, so
a traverse case can be fed in through the same three numbers.

H0 is not evaluated from r_tip and r_ground independently. It is pinned to the
deck's own tangency: at t = 0 the governing grit vertex sits exactly the
standoff clear of the ground face, so h at that vertex's station is known to
the picometre, and H0 follows. Deriving it any other way would let a
sub-micron disagreement between "tallest tip" and "tallest tip inside the
footprint" leak into a quantity being compared against a 5 nm threshold.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

# ``vumat_grind.for`` reads exactly this many constants, and SDV12 is its
# deletion flag. Changing either here without changing the Fortran is the one
# way to make a deck that runs and is wrong, so both are asserted on write.
N_HYBRID_PROPS = 56
HYBRID_DEPVAR = 20
HYBRID_DELETE_SDV = 12

# ---------------------------------------------------------------------------
# Placeholder ductile-regime constants for the sandstone the JH-2 card
# describes. PLACEHOLDERS, and the deck says so in a comment block.
#
# A is not invented: 90.0 MPa is the quasi-static uniaxial compressive
# strength of that same JH-2 card (the intersection of the elastic path with
# the intact strength surface, which vumat_jh2.for's header quotes and
# verify_vumat_grind.py re-derives). Starting the ductile branch's yield there
# makes the two laws agree at the transition in uniaxial compression instead
# of stepping across it.
#
# Everything else is a defensible order of magnitude for a quartz-bonded rock
# and nothing more. Calibrate B, n, C, m and D1..D5 against your own
# nanoindentation or scratch data before quoting a force from this model.
# ---------------------------------------------------------------------------
JC_SANDSTONE_PLACEHOLDER = dict(
    a_mpa=90.0, b_mpa=50.0, n=0.50, c=0.020, m=1.0, edot0=1.0,
    youngs_mpa=6500.0, poisson=0.21,
    density_kg_m3=2350.0, specific_heat_j_kgk=800.0, taylor_quinney=0.9,
    t0_k=293.15, tmelt_k=1473.15,
    burgers_mm=5.0e-7, taylor_factor=3.0, alpha=0.3, sge_exponent=1.0,
    r_prime=2.0, sge_shear_mpa=0.0,
    d1=0.0, d2=0.15, d3=-1.5, d4=0.0, d5=0.0, dcrit=1.0,
)

# Hardness and toughness for the same rock, also placeholders.
DC_SANDSTONE_PLACEHOLDER = dict(hardness_mpa=1000.0, kic_mpa_sqrt_mm=0.3,
                                lambda_c=0.15, dc_form=1)


class HybridError(RuntimeError):
    pass


@dataclass
class HybridParams:
    """Everything ``vumat_grind.for`` needs beyond the JH-2 card.

    Units follow the deck: mm, tonne, s, MPa, N. Two of those bite:

    * ``kic`` is a stress intensity, so in this system it is MPa*sqrt(mm),
      which is 31.623 times the MPa*sqrt(m) that toughness is normally
      quoted in. Use :func:`kic_from_mpa_sqrt_m`.
    * ``specific_heat_j_kgk`` is given in the familiar J/(kg K) and converted
      on the way into the card.
    """

    enabled: bool = False

    # ---- Johnson-Cook flow ---------------------------------------------
    a_mpa: float = 90.0
    b_mpa: float = 50.0
    n: float = 0.50
    c: float = 0.020
    m: float = 1.0
    edot0: float = 1.0

    youngs_mpa: float = 6500.0
    poisson: float = 0.21
    density_kg_m3: float = 2350.0
    specific_heat_j_kgk: float = 800.0
    taylor_quinney: float = 0.9
    t0_k: float = 293.15
    tmelt_k: float = 1473.15

    # ---- strain-gradient enhancement -----------------------------------
    burgers_mm: float = 5.0e-7
    """Burgers vector. Also the floor on the gradient length: eta = 4 ep / h
    is a dislocation density divided by b, so a length below b would describe
    a lattice curvature no lattice can hold, and h -> 0 in the rubbing zone
    would otherwise return an infinite flow stress."""
    taylor_factor: float = 3.0
    alpha: float = 0.3
    sge_exponent: float = 1.0
    """Lambda. 1.0 is the blanking paper's eq. 7; the micro-milling paper
    leaves it as a fitted exponent delta."""
    r_prime: float = 2.0
    """r'. 2.0 makes the blanking form identical to the peening and
    micro-milling forms."""
    sge_shear_mpa: float = 0.0
    """Shear modulus used inside the SGE term. 0 = derive from E and nu."""

    # ---- Johnson-Cook damage --------------------------------------------
    d1: float = 0.0
    d2: float = 0.15
    d3: float = -1.5
    d4: float = 0.0
    d5: float = 0.0
    dcrit: float = 1.0

    # ---- the transition --------------------------------------------------
    dc_mm: float = 0.0
    """Critical depth of cut. 0 = compute it from hardness and toughness."""
    lambda_c: float = 0.15
    hardness_mpa: float = 1000.0
    kic: float = 0.3
    """Fracture toughness in MPa*sqrt(mm)."""
    dc_form: int = 1
    """1: dc = lambda_c (H/E)^0.5 (Kc/H)^2, the form on the project's slide.
    2: dc = lambda_c (E/H) (Kc/H)^2, Bifano, Dow & Scattergood (1991), whose
    calibrated lambda_c is 0.15.

    They differ by (E/H)^1.5 -- a factor of 17 on this rock -- so lambda_c is
    NOT transferable between them. Whichever you pick, say which in the paper.
    """

    h_source: int = 0
    """0 h from the coordinates (single grit), 1 from field variable 1,
    2 force ductile everywhere, 3 force brittle everywhere. 2 and 3 exist so a
    hybrid deck can be run as a pure-JC or pure-JH-2 deck without rebuilding
    it, which is how you find out how much of a result the switch caused."""

    # ---- prescribing the trajectory directly ---------------------------
    # Normally H0 and HG are DERIVED from where the grit was seated and how fast
    # the wheel feeds in, and that derivation is the thing that keeps the card
    # honest -- it cannot disagree with the geometry the deck writes.
    #
    # These two exist for the case where the trajectory is the specification
    # rather than the consequence: "cut this exact profile", from a measured
    # groove or a figure in a paper. Setting either one overrides that term of
    # h(u) = H0 + HG*u - u^2/(2*RTIP) and the deck says so in its own header, so
    # a reader can see the number was imposed rather than derived.
    #
    # RTIP is deliberately NOT overridable: it is the grit's real radius on the
    # real wheel, and letting it be typed would allow a curvature no wheel in
    # the model has.
    h0_override_mm: Optional[float] = None
    """Peak/offset term of h(u). None leaves it derived from the seating."""
    hg_override: Optional[float] = None
    """Wedge slope of h(u). None leaves it derived from the infeed. Zero is a
    real value here and means a pure arc: no infeed, depth from curvature."""

    placeholder: bool = True
    """True means the ductile constants have not been calibrated for this
    material and the deck should say so, loudly, in its own header."""

    def validate(self) -> None:
        if self.youngs_mpa <= 0:
            raise HybridError("youngs_mpa must be positive")
        if not -1.0 < self.poisson < 0.5:
            raise HybridError("poisson must be in (-1, 0.5)")
        if self.a_mpa <= 0:
            raise HybridError(
                "a_mpa is the ductile yield stress and must be positive; a "
                "zero here would make the ductile branch flow at no load")
        if self.burgers_mm <= 0:
            raise HybridError("burgers_mm must be positive: it is the floor "
                              "on the strain-gradient length")
        if self.dc_form not in (1, 2):
            raise HybridError("dc_form must be 1 or 2")
        if self.h_source not in (0, 1, 2, 3):
            raise HybridError("h_source must be 0, 1, 2 or 3")
        if self.dcrit <= 0 or self.dcrit > 1.0:
            raise HybridError("dcrit must be in (0, 1]")
        if self.dc_mm <= 0 and not (self.hardness_mpa > 0 and self.kic > 0):
            raise HybridError(
                "give dc_mm directly, or a positive hardness and toughness "
                "to compute it from")
        if self.critical_depth_mm() <= 0:
            raise HybridError("the critical depth of cut came out non-positive")

    # -- derived ----------------------------------------------------------
    @property
    def shear_mpa(self) -> float:
        return self.youngs_mpa / (2.0 * (1.0 + self.poisson))

    def critical_depth_mm(self) -> float:
        if self.dc_mm > 0:
            return float(self.dc_mm)
        return critical_depth_mm(self.lambda_c, self.hardness_mpa,
                                 self.youngs_mpa, self.kic, self.dc_form)


def kic_from_mpa_sqrt_m(kic_mpa_sqrt_m: float) -> float:
    """MPa*sqrt(m) -> MPa*sqrt(mm), which is what a mm-MPa deck needs.

    Getting this wrong scales dc by 1000, which is the difference between a
    transition at 5 nm and one at 5 um -- i.e. between a model that is mostly
    brittle and one that is mostly ductile. It is a conversion worth having a
    function for.
    """
    return float(kic_mpa_sqrt_m) * math.sqrt(1000.0)


def critical_depth_mm(lambda_c: float, hardness_mpa: float, youngs_mpa: float,
                      kic: float, form: int = 1) -> float:
    """The critical depth of cut, in the same length unit ``kic`` implies.

    ``form`` 1 is ``lambda_c (H/E)^0.5 (Kc/H)^2`` and 2 is Bifano's
    ``lambda_c (E/H) (Kc/H)^2``. Both are dimensionally sound: ``(Kc/H)^2``
    carries the length and the modulus ratio is dimensionless. They come from
    the same indentation-fracture family and differ only in how the crack
    initiation is referenced, which is why the constant in front is not
    shared.
    """
    if hardness_mpa <= 0 or youngs_mpa <= 0:
        raise HybridError("hardness and Young's modulus must be positive")
    ratio = (kic / hardness_mpa) ** 2
    if form == 2:
        return lambda_c * (youngs_mpa / hardness_mpa) * ratio
    return lambda_c * math.sqrt(hardness_mpa / youngs_mpa) * ratio


# ---------------------------------------------------------------------------
# the chip-thickness field
# ---------------------------------------------------------------------------

@dataclass
class ChipField:
    """h(u) = h0 + hg*u - u^2/(2*rtip), with u the tangential station."""

    theta_c: float          # rad, the frame the field is written in
    h0_mm: float
    hg: float               # dimensionless slope, dh/du
    rtip_mm: float
    u_gov_mm: float         # station of the vertex the block was seated on
    h_entry_mm: float       # chip thickness where the grit meets the block
    h_exit_mm: float        # and where it leaves
    transition_u_mm: Optional[float] = None
    """FIRST station where h crosses dc, if any crossing lies inside the block.
    This is the number the whole model exists to produce: the point along the
    scratch where removal stops being ductile."""
    transition_all_mm: tuple = ()
    """Every crossing, not just the first. A prescribed arc that rises above dc
    and comes back down crosses twice, and a header that mentioned only one of
    them would describe half the experiment."""

    def h_at(self, u: float) -> float:
        h = self.h0_mm + self.hg * u
        if self.rtip_mm > 0:
            h -= u * u / (2.0 * self.rtip_mm)
        return max(h, 0.0)


def chip_field(place: dict, motion: Optional[dict], wp,
               *, rotation_reversed: bool = False,
               dc_mm: float = 0.0,
               h0_override: Optional[float] = None,
               hg_override: Optional[float] = None) -> ChipField:
    """Work out H0, HG and RTIP for the grit this deck actually placed.

    ``place`` is :func:`semgrit.rigid_wheel.place_workpiece`'s output and
    ``motion`` is :func:`semgrit.analysis.wheel_motion`'s, so the field is
    derived from the same placement and the same infeed the deck writes -- not
    from a re-derivation that could disagree with it.
    """
    frames = place["frames"]
    gov = place["gov"]
    if gov is None or not frames:
        raise HybridError("the workpiece was never seated on a grit, so there "
                          "is no scratch to define a chip thickness along")
    if wp is None:
        raise HybridError("a hybrid deck needs a workpiece")

    # The governing vertex: the one inside the footprint with the largest
    # radial reach. That is the vertex the ground face was placed tangent to,
    # so its chip thickness at t = 0 is exactly the standoff, and everything
    # else follows from it.
    v = frames[gov]
    hb, hz = wp.length_mm / 2.0, wp.width_mm / 2.0
    inside = (np.abs(v[:, 1]) <= hb) & (np.abs(v[:, 2]) <= hz)
    if not inside.any():
        inside = np.ones(len(v), dtype=bool)
    idx = int(np.argmax(np.where(inside, v[:, 0], -np.inf)))
    a_gov, u_gov = float(v[idx, 0]), float(v[idx, 1])
    rtip = math.hypot(a_gov, u_gov)
    if rtip <= 0:
        raise HybridError("the governing grit vertex sits on the wheel axis")

    # Stations are crossed at the speed of the GRIT TIP, omega*r_tip, not at
    # the speed of the ground face. The two differ by the protrusion, which is
    # parts per million of the radius -- but so is nothing else in this file,
    # and the whole point of pinning H0 to the tangency was to stop that class
    # of small inconsistency reaching a nanometre-scale threshold.
    if motion is None:
        v_r, v_s = 0.0, 1.0
    else:
        v_r = float(motion["radial_speed_mm_s"])
        v_s = float(motion["omega_rad_s"]) * rtip
        if v_s <= 0:
            raise HybridError("the wheel is not turning, so no grit ever "
                              "reaches a station: h(u) is undefined")

    r_ground = float(place["r_ground"])
    # h at the governing station at t = 0 is minus the standoff: the ground
    # face was put that far above the vertex.
    h_gov = a_gov - r_ground

    way = -1.0 if rotation_reversed else 1.0
    hg = -way * v_r / v_s
    h0 = h_gov - hg * u_gov + (u_gov * u_gov / (2.0 * rtip) if rtip > 0 else 0.0)

    if h0_override is not None:
        h0 = float(h0_override)
    if hg_override is not None:
        hg = float(hg_override)

    fld = ChipField(theta_c=float(place["theta_c"]), h0_mm=h0, hg=hg,
                    rtip_mm=rtip, u_gov_mm=u_gov,
                    h_entry_mm=0.0, h_exit_mm=0.0)
    # The entry edge is the end the grains arrive from, which follows the
    # rotation sense exactly as the placement does.
    u_entry = way * hb
    fld.h_entry_mm = fld.h_at(u_entry)
    fld.h_exit_mm = fld.h_at(-u_entry)
    if dc_mm > 0:
        st = transition_stations(fld, -hb, hb, dc_mm)
        fld.transition_all_mm = tuple(st)
        fld.transition_u_mm = st[0] if st else None
    return fld


def transition_stations(fld: ChipField, u_lo: float, u_hi: float,
                        dc: float, samples: int = 2001) -> list:
    """EVERY station in [u_lo, u_hi] where h(u) crosses dc, in order.

    h(u) = H0 + HG*u - u^2/(2*RTIP) is a downward parabola, so it can cross a
    level twice. For a deck whose depth comes from a radial infeed the linear
    term dominates over a block far shorter than the wheel radius, h is monotone
    in practice, and a single bisection between the two ends is right -- which
    is what this function used to do, and why it only ever returned one station.

    That assumption fails the moment the profile is prescribed rather than fed
    in. An arc that starts at zero depth, peaks above dc and returns to zero
    crosses dc twice and has EQUAL h at both ends, so an endpoint test sees no
    sign change at all and reports "never crosses" -- the confident wrong answer
    a deck header should never give. Scanning for sign changes and refining each
    one costs a few thousand evaluations of a quadratic and cannot miss a root
    the sampling resolves.
    """
    if samples < 3:
        samples = 3
    span = u_hi - u_lo
    us = [u_lo + span * i / (samples - 1) for i in range(samples)]
    fs = [fld.h_at(u) - dc for u in us]
    out = []
    for i in range(samples):
        if fs[i] == 0.0:
            out.append(us[i])
            continue
        if i == 0:
            continue
        if fs[i - 1] * fs[i] < 0.0:
            lo, hi, f_lo = us[i - 1], us[i], fs[i - 1]
            for _ in range(120):
                mid = 0.5 * (lo + hi)
                if (fld.h_at(mid) - dc > 0.0) == (f_lo > 0.0):
                    lo = mid
                else:
                    hi = mid
            out.append(0.5 * (lo + hi))
    return out


def _transition_station(fld: ChipField, u_lo: float, u_hi: float,
                        dc: float) -> Optional[float]:
    """The FIRST station where h(u) crosses dc, or None.

    Kept for the callers that want a single number -- notably the one field on
    :class:`ChipField`. Use :func:`transition_stations` when the profile may be
    non-monotone; this returns only the earliest of however many there are.
    """
    st = transition_stations(fld, u_lo, u_hi, dc)
    return st[0] if st else None


# ---------------------------------------------------------------------------
# the material card
# ---------------------------------------------------------------------------

KGM3_TO_TONNE_MM3 = 1.0e-12
# J/(kg K) -> mJ/(tonne K), the consistent unit in mm-MPa-tonne-s.
JKGK_TO_MJ_TONNEK = 1.0e6


def hybrid_props(jh2_constants: Sequence[float], p: HybridParams,
                 fld: Optional[ChipField],
                 *, edot0: float = 1.0, edmin: float = 1.0,
                 itcut: int = 1, fsmax: float = -1.0) -> list[float]:
    """The 56 constants ``vumat_grind.for`` reads, in its order.

    ``jh2_constants`` is the existing 17-value JH-2 card, unchanged, so a deck
    that was running the brittle model keeps exactly the brittle model it had
    wherever h >= dc.
    """
    j = list(jh2_constants)
    if len(j) != 17:
        raise HybridError("the JH-2 card is 17 constants; got %d" % len(j))
    out = j + [edot0, edmin, float(itcut), fsmax]
    out += [p.a_mpa, p.b_mpa, p.n, p.c, p.m, p.edot0,
            p.youngs_mpa, p.poisson,
            p.density_kg_m3 * KGM3_TO_TONNE_MM3,
            p.specific_heat_j_kgk * JKGK_TO_MJ_TONNEK,
            p.taylor_quinney, p.t0_k, p.tmelt_k,
            p.burgers_mm, p.taylor_factor, p.alpha, p.sge_exponent,
            p.r_prime, p.sge_shear_mpa]
    out += [p.d1, p.d2, p.d3, p.d4, p.d5, p.dcrit]
    out += [p.critical_depth_mm(), p.lambda_c, p.hardness_mpa, p.kic,
            float(p.dc_form)]
    if fld is None:
        out += [0.0, 0.0, 0.0, 0.0, float(p.h_source)]
    else:
        out += [fld.theta_c, fld.h0_mm, fld.hg, fld.rtip_mm,
                float(p.h_source)]
    if len(out) != N_HYBRID_PROPS:
        raise HybridError("built %d constants, vumat_grind.for reads %d"
                          % (len(out), N_HYBRID_PROPS))
    return out


def plan_hybrid(plan: dict, p: Optional[HybridParams] = None):
    """The chip-thickness field a hybrid deck *would* carry. Writes nothing.

    Takes a :func:`semgrit.build_deck.plan_deck` result, so the preview and the
    build read the same placement and the same infeed. Returns
    ``(ChipField, dc_mm)``.
    """
    import dataclasses

    from .analysis import wheel_motion

    pr = plan["_params"]
    an = pr.analysis
    if p is None:
        p = getattr(an, "hybrid", None)
    if p is None:
        raise HybridError("no HybridParams: set analysis.hybrid")
    dc = p.critical_depth_mm()
    step_time = float((plan.get("cost") or {}).get("step_time_s") or 0.0)
    if step_time <= 0:
        raise HybridError("the step has no duration, so there is no infeed "
                          "rate and no chip-thickness ramp")
    # plan_deck has already resolved depth_of_cut_um = 0 into a real number.
    an2 = dataclasses.replace(an,
                              depth_of_cut_um=float(plan["depth_of_cut_um"]))
    motion = wheel_motion(an2, plan["_place"]["theta_c"],
                          pr.surface_speed_mm_s, pr.outer_radius_mm, step_time)
    fld = chip_field(plan["_place"], motion, plan["_wp"],
                     rotation_reversed=bool(an.rotation_reversed), dc_mm=dc,
                     h0_override=getattr(p, "h0_override_mm", None),
                     hg_override=getattr(p, "hg_override", None))
    return fld, dc


def summary_text(fld: ChipField, dc_mm: float, wp, p: HybridParams) -> str:
    """What the switch will do to this deck, as a block of text."""
    hb = wp.length_mm / 2.0
    L: list[str] = []
    a = L.append
    a("DUCTILE / BRITTLE TRANSITION")
    a("  critical depth of cut dc : %.4f nm   (%.6e mm)"
      % (dc_mm * 1e6, dc_mm))
    a("    %s" % ("dc = lambda_c (H/E)^0.5 (Kc/H)^2" if p.dc_form == 1
                  else "dc = lambda_c (E/H) (Kc/H)^2   [Bifano 1991]"))
    a("    lambda_c %.4g, H %.4g MPa, E %.4g MPa, Kc %.4g MPa*sqrt(mm)"
      % (p.lambda_c, p.hardness_mpa, p.youngs_mpa, p.kic))
    a("")
    a("CHIP THICKNESS ALONG THE SCRATCH   h(u) = H0 + HG u - u^2/(2 RTIP)")
    a("  H0 %.6e mm, HG %.6e, RTIP %.4f mm" % (fld.h0_mm, fld.hg, fld.rtip_mm))
    a("  at the entry edge (u = %+.4f mm) : %8.4f nm" % (hb, fld.h_entry_mm * 1e6))
    a("  at the exit edge  (u = %+.4f mm) : %8.4f nm" % (-hb, fld.h_exit_mm * 1e6))
    a("")
    if fld.transition_u_mm is None:
        if fld.h_entry_mm >= dc_mm and fld.h_exit_mm >= dc_mm:
            a("  h >= dc EVERYWHERE: the whole scratch is BRITTLE (JH-2).")
            a("  To see a ductile regime, cut shallower, raise dc, or use the")
            a("  Bifano form -- it gives a dc about (E/H)^1.5 larger.")
        else:
            a("  h < dc EVERYWHERE: the whole scratch is DUCTILE (JC+SGE).")
            a("  Nothing will fracture. Cut deeper or lower dc.")
    else:
        u_t = fld.transition_u_mm
        ductile_mm = abs(hb - u_t) if fld.hg < 0 else abs(u_t + hb)
        a("  TRANSITION at u = %+.6f mm from the block centre." % u_t)
        a("  ductile for the first %.4f mm of the %.4f mm scratch (%.1f%%),"
          % (ductile_mm, 2 * hb, 100.0 * ductile_mm / (2 * hb)))
        a("  brittle for the rest. Plot SDV13 after the run to see it.")
    a("")
    a("  SDV13 = branch (1 ductile, 2 brittle)   SDV14 = h at that point")
    a("  SDV15 = dc                              SDV19 = SGE amplification")
    if p.h_source in (2, 3):
        a("")
        a("  NOTE: the switch is OVERRIDDEN -- %s"
          % _H_SOURCE_TEXT[p.h_source])
    if p.placeholder:
        a("")
        a("  NOTE: the Johnson-Cook constants are PLACEHOLDERS. A is this")
        a("  material's own JH-2 quasi-static compressive strength so the two")
        a("  branches meet at the transition; the rest are order-of-magnitude.")
    return "\n".join(L)


_H_SOURCE_TEXT = {
    0: "h from the material point's own coordinates (single grit)",
    1: "h from field variable 1",
    2: "FORCED DUCTILE everywhere -- the switch is disabled",
    3: "FORCED BRITTLE everywhere -- this is plain JH-2",
}


def write_hybrid_material(w, p: HybridParams, wp, jh2_constants,
                          fld: Optional[ChipField], n_depvar: int,
                          element_deletion: bool) -> dict:
    """Write the ``*Material`` block for the hybrid law. Returns a summary."""
    p.validate()
    props = hybrid_props(jh2_constants, p, fld)
    dc = p.critical_depth_mm()

    w("**\n")
    w("** ---------------- HYBRID DUCTILE / BRITTLE MATERIAL --------------\n")
    w("** Johnson-Cook + strain-gradient enhancement where the grit takes a\n")
    w("** cut thinner than dc, Johnson-Holmquist II where it takes a thicker\n")
    w("** one. Requires vumat_grind.for:\n")
    w("**   abaqus job=... input=... user=vumat_grind.for double=both cpus=N\n")
    w("**\n")
    w("** critical depth of cut dc  : %.6e mm  (%.4f nm)\n" % (dc, dc * 1e6))
    w("**   from   %s\n" % ("lambda_c*(H/E)^0.5*(Kc/H)^2" if p.dc_form == 1
                            else "lambda_c*(E/H)*(Kc/H)^2  [Bifano 1991]"))
    w("**   lambda_c %.6g, H %.6g MPa, E %.6g MPa, Kc %.6g MPa*sqrt(mm)\n"
      % (p.lambda_c, p.hardness_mpa, p.youngs_mpa, p.kic))
    w("** h source                  : %s\n" % _H_SOURCE_TEXT[p.h_source])
    if fld is not None and p.h_source == 0:
        w("** h(u) = H0 + HG*u - u^2/(2*RTIP), u = tangential station in mm\n")
        w("**   theta_c %.9f rad, H0 %.9e mm, HG %.9e, RTIP %.6f mm\n"
          % (fld.theta_c, fld.h0_mm, fld.hg, fld.rtip_mm))
        # Say which terms were imposed. A derived H0 is a consequence of
        # the seating and can be checked against it; an imposed one is a
        # specification and cannot. A reader is entitled to know which.
        imposed = []
        if getattr(p, "h0_override_mm", None) is not None:
            imposed.append("H0")
        if getattr(p, "hg_override", None) is not None:
            imposed.append("HG")
        if imposed:
            w("**   %s PRESCRIBED, not derived from seating and infeed:\n"
              % " and ".join(imposed))
            w("**   this deck cuts a SPECIFIED trajectory. The geometry\n")
            w("**   is still the real one; only the depth profile was\n")
            w("**   imposed.\n")
        w("**   h at the entry edge   : %.6e mm  (%.4f nm)\n"
          % (fld.h_entry_mm, fld.h_entry_mm * 1e6))
        w("**   h at the exit edge    : %.6e mm  (%.4f nm)\n"
          % (fld.h_exit_mm, fld.h_exit_mm * 1e6))
        # Report EVERY crossing. A prescribed arc rises above dc and comes back
        # down, so it transitions twice, and the endpoints are equal -- naming
        # only the first, or testing only the ends, describes half the run.
        st = list(getattr(fld, "transition_all_mm", ())
                  or ([fld.transition_u_mm]
                      if fld.transition_u_mm is not None else []))
        if not st:
            # No crossing at all: decide which regime from a point that is
            # actually inside the cut, not from the ends, which on an arc are
            # both at zero depth and would always read ductile.
            mid = fld.h_at(0.0)
            side = "ENTIRELY BRITTLE" if mid >= dc else "ENTIRELY DUCTILE"
            w("**   the whole scratch is %s: h never crosses dc\n"
              % side)
        elif len(st) == 1:
            w("**   ductile-brittle transition at u = %+.6f mm from the block\n"
              % st[0])
            w("**   centre -- plot SDV13 to see it.\n")
        else:
            w("**   %d ductile-brittle transitions, at u =\n" % len(st))
            for u in st:
                w("**     %+.6f mm from the block centre\n" % u)
            w("**   the cut rises through dc and falls back, so the scratch\n")
            w("**   runs ductile - brittle - ductile. Plot SDV13 to see it.\n")
    if p.placeholder:
        w("**\n")
        w("** WARNING: the Johnson-Cook constants below are PLACEHOLDERS for\n")
        w("** this material. A is the JH-2 card's own quasi-static uniaxial\n")
        w("** compressive strength so the two branches meet at the\n")
        w("** transition, but B, n, C, m and D1..D5 are order-of-magnitude\n")
        w("** values. Calibrate them before quoting a force from this run.\n")
    w("**\n")
    w("*Material, name=%s\n" % wp.material)
    w("*Density\n%.8e,\n" % (p.density_kg_m3 * KGM3_TO_TONNE_MM3))
    if element_deletion:
        w("*Depvar, delete=%d\n%d,\n" % (HYBRID_DELETE_SDV, n_depvar))
    else:
        w("*Depvar\n%d,\n" % n_depvar)
    # Full round-trip precision, four to a line. repr() gives the shortest
    # decimal string that reads back as the same double, so the card is exact
    # rather than merely close.
    #
    # This is not fussiness. dc here is 5.3e-6 mm and H0 is 9.2e-4 mm, so a
    # card written at the %g default would hand the subroutine a chip
    # thickness whose last digits are noise, against a threshold six orders
    # smaller than the block. The deck already writes coordinates at %.12e for
    # the same reason: this project has been bitten once by a rounding that
    # moved a tangency, and the switch is more sensitive than the tangency was.
    # EIGHT per line. Not a style choice: Abaqus reads *User Material
    # constants eight to a line, and four to a line is rejected outright with
    #   ***ERROR: THERE ARE INVALID DATA ASSOCIATED WITH THIS USER DEFINED
    #             MATERIAL DEFINITION
    # which is what every deck did on its first real submission. Eight values
    # at 17 significant digits is about 190 characters, inside Abaqus' 256.
    w("*User Material, constants=%d\n" % len(props))
    for i in range(0, len(props), 8):
        w(", ".join(repr(float(v)) for v in props[i:i + 8]) + "\n")

    return {
        "dc_mm": dc,
        "dc_nm": dc * 1e6,
        "dc_form": p.dc_form,
        "n_props": len(props),
        "n_depvar": n_depvar,
        "h_source": p.h_source,
        "props": props,
        "chip_field": (None if fld is None else dict(
            theta_c_rad=fld.theta_c, h0_mm=fld.h0_mm, hg=fld.hg,
            rtip_mm=fld.rtip_mm, u_gov_mm=fld.u_gov_mm,
            h_entry_mm=fld.h_entry_mm, h_exit_mm=fld.h_exit_mm,
            transition_u_mm=fld.transition_u_mm,
            transition_all_mm=list(getattr(fld, "transition_all_mm", ())),
            # Which terms were imposed rather than derived. The geometry
            # verifier checks H0 and HG against the seating and the infeed,
            # which is exactly the right check for a deck whose profile is a
            # CONSEQUENCE -- and exactly the wrong one for a deck whose profile
            # is the SPECIFICATION. It has to be told which it is looking at.
            h0_prescribed=getattr(p, "h0_override_mm", None) is not None,
            hg_prescribed=getattr(p, "hg_override", None) is not None)),
        "placeholder_constants": bool(p.placeholder),
    }
