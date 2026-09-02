"""Two Abaqus decks for shape-adaptive grinding, at the two scales it needs.

SAG cannot be put in one deck, and the reason is arithmetic rather than taste.
At the paper's own operating point the contact patch is 128 mm^2 and holds
32,000 grains on a 30 um pad or 179,000 on a 6 um pad. Resolving dc = 80 nm
across that area at the project standard of five elements per dc needs about
5e11 elements. For scale, the finest deck this project has ever run
(``ARC_80NM``) is 172,000 elements and 74 hours on eight cores -- for ONE grain.

So the process is modelled at both ends and each deck answers what it can:

``MACRO`` -- the whole compliant wheel
    Rigid hub, hyperelastic/viscoelastic polyurethane ring, every grain on the
    pad present as measured geometry, pressed into the workpiece by the wheel
    compression T and then rotated. Its job is the *contact*: how big the patch
    really is, how the pressure is distributed over it, how many grains are
    actually engaged and what load each one carries. The workpiece is meshed for
    contact, not for dc, so this deck deliberately CANNOT show a transition --
    and says so in its own header rather than leaving a reader to assume.

``MICRO`` -- one patch of that contact, resolved
    A small block with a handful of measured grains, meshed at dc/5, driven by
    the force per grain the MACRO deck computed. Its job is the *transition*:
    SDV13, ductile against brittle, at a resolution that can carry it.

The link between them is a single number -- the load one grain carries -- and it
is written into both decks' headers so the pair cannot be quoted out of step.
This is the standard two-scale argument, and it is honest as long as the
coupling is stated: MACRO gives MICRO its boundary condition, MICRO gives back a
removal mode that MACRO cannot resolve.

Everything is reused rather than reimplemented: the grain solids come from the
SEM pipeline, ``semgrit.sag`` supplies the contact mechanics, and the MICRO deck
is a normal hybrid deck written by ``build_deck`` with h delivered through the
VUMAT's existing field-variable path. There is no new constitutive law here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Sequence

from . import sag

# The polyurethane in the reference deck, verbatim, with its density corrected.
#
# The reference carries rho = 2.2e-13 tonne/mm^3, which is 0.22 kg/m^3 -- about
# 5000x too light, lighter than air, and lighter than any real foam (30-100
# kg/m^3 open-cell, ~1100 solid). In Abaqus/Explicit density sets the stable
# increment, so an absurdly light ring runs faster; it also has no inertia. That
# is a speed hack, and it is not defensible in a paper, so the default here is
# real polyurethane and mass scaling is applied openly and reported.
PU_SOLID_DENSITY_KG_M3 = 1100.0
PU_FOAM_DENSITY_KG_M3 = 60.0
FACETS_PER_GRAIN = 121.0
"""Measured, not assumed: RUN_ME/2_multi_abrasive writes 1,452 grit facets for
12 grains. General-contact cost scales with facet count, so this is what turns a
grain count into a cost."""

PU_REFERENCE_DECK_DENSITY_KG_M3 = 0.22
"""What the reference deck actually used. Kept so ``--match-reference`` can
reproduce it exactly, clearly labelled, rather than the number being lost."""


PRESS_MACH = 0.005
"""Press-in speed as a fraction of the compliant layer's dilatational wave
speed. The press has to be quasi-static or the contact pressure is inertial
rather than Hertzian: a ramp fast compared with the wave speed launches a stress
wave through the layer and the patch never reaches equilibrium. The reference
deck presses ~0.8 mm in 0.0005 s, which is v/c = 0.09 and only 1.8 wave transits
of a 5 mm layer -- so its contact pressure is not the steady one its own Hertz
comparison assumes."""

HOLD_TAUS = 3.0
"""Dwell, in Prony relaxation times, between pressing and grinding.

The polyurethane is viscoelastic, so its stiffness depends on how long the load
has been there. The paper measures a STEADY force with a load cell, and the
Hertz model it compares against uses ``moduli=LONG TERM`` -- the RELAXED
modulus. A deck that presses and immediately grinds is still in the glassy
state and would carry a stiffer layer, a smaller patch and a higher force per
grain than the experiment. Three time constants leaves ~5% of the relaxation
unfinished, which is below the spread on the measured forces."""


class SAGDeckError(RuntimeError):
    pass


@dataclass
class Polyurethane:
    """The compliant layer's material card.

    Neo-Hookean plus a one-term Prony series, which is what the reference deck
    uses and what a polyurethane at these rates needs: hyperelastic because the
    strains are large, viscoelastic because the contact is transient and a
    purely elastic layer would give back all of its energy.
    """

    c10_mpa: float = 0.0575
    d1: float = 1.0
    prony_g: float = 0.11
    prony_k: float = 0.05
    prony_tau_s: float = 0.01
    density_kg_m3: float = PU_SOLID_DENSITY_KG_M3
    thickness_mm: float = 5.0

    def __post_init__(self):
        if self.c10_mpa <= 0:
            raise SAGDeckError("C10 must be positive")
        if self.thickness_mm <= 0:
            raise SAGDeckError("the compliant layer needs a positive thickness")
        if not 0.0 <= self.prony_g < 1.0:
            raise SAGDeckError(
                "Prony g1 = %g. It is the RELAXED fraction of the shear "
                "modulus and must be in [0, 1); at 1 the layer has no "
                "long-term stiffness at all." % self.prony_g)
        if self.density_kg_m3 <= 0:
            raise SAGDeckError("density must be positive")

    @property
    def modulus_mpa(self) -> float:
        """Small-strain E, for the Hertz model. E = 6 C10 at incompressibility."""
        return sag.modulus_from_neo_hooke(self.c10_mpa)

    @property
    def shear_mpa(self) -> float:
        return 2.0 * self.c10_mpa

    def density_tonne_mm3(self) -> float:
        return self.density_kg_m3 * 1e-12

    def long_term_fraction(self) -> float:
        """Share of the instantaneous shear modulus still carried at t -> inf."""
        return 1.0 - self.prony_g

    def cards(self) -> list:
        """The ``*Material`` block, 8 values per line as Abaqus requires."""
        return [
            "*Density",
            " %.7g," % self.density_tonne_mm3(),
            "*Hyperelastic, neo hooke, moduli=LONG TERM",
            " %.7g, %.7g" % (self.c10_mpa, self.d1),
            "*Viscoelastic, time=PRONY",
            " %.7g, %.7g, %.7g" % (self.prony_g, self.prony_k,
                                   self.prony_tau_s),
        ]


@dataclass
class SAGParams:
    """One SAG operating point, plus how to turn it into decks."""

    # --- the tool ---------------------------------------------------------
    diameter_mm: float = 125.0
    width_mm: float = 10.0
    grain_um: float = 6.0
    pad_areal_per_mm2: float = 0.0
    hub_diameter_mm: float = 0.0
    """Rigid hub OD. 0 -> diameter - 2*layer thickness."""
    polyurethane: Polyurethane = field(default_factory=Polyurethane)
    shore_a: float = 40.0
    use_shore_modulus: bool = True
    """True -> Hertz uses eq. 4 on ``shore_a`` (which reproduces the paper's
    measured forces). False -> uses 6*C10 from the PU card (which matches the
    reference deck, and is ~5x softer)."""

    # --- the process ------------------------------------------------------
    compression_mm: float = 0.4
    speed_rpm: float = 1050.0
    friction: float = 0.2
    press_time_s: float = 0.0
    """Press-in duration. 0 -> derived so the press stays quasi-static
    (v/c = PRESS_MACH of the compliant layer's own wave speed)."""
    hold_time_s: float = 0.0
    """Dwell after pressing. 0 -> HOLD_TAUS * the Prony relaxation time, so the
    layer relaxes to its long-term modulus before grinding starts."""
    grind_time_s: float = 0.02

    # --- the workpiece ----------------------------------------------------
    material: str = "wc_co"
    carbide_um: float = 1.36
    """Mean carbide size, for the paper's k = dg/D_WC < 5 criterion."""
    bhn_kgf_mm2: float = 581.0

    # --- discretisation ---------------------------------------------------
    elements_per_dc: float = 5.0
    """Elements across dc, THROUGH THE DEPTH. The project standard is 5."""
    elements_per_grain: float = 20.0
    """Elements across one grain, IN-PLANE. dc/5 in-plane would spend the
    whole budget on the empty space between grains; what the in-plane mesh has
    to carry is the shape of the cutting edge."""
    micro_patch_mm: float = 0.0
    """Side of the resolved MICRO patch. 0 -> sized from the grain count."""
    micro_grains: int = 3
    macro_grain_cap: int = 400_000
    """Most grains the MACRO deck will place, and it is deliberately large.

    A cap is honoured by NARROWING THE SECTOR, never by thinning the pad, so a
    small cap does not give a cheap approximation of the real contact -- it
    gives a narrow sector whose own curvature cannot span the indent, i.e. a
    flat punch. At T = 0.4 mm on a 62.5 mm radius the sector must be at least
    12.97 deg, and a 6 um pad at 1,400 grains/mm2 puts ~262,000 grains in the
    17 deg the contact arc wants. At 121 facets each that is ~32 M
    general-contact facets, which is expensive and is what fidelity costs here.

    Lower it only if a narrow, non-contacting sector is genuinely what is
    wanted -- ``plan`` reports ``spans_indent`` so that choice is visible."""
    macro_sector_mode: str = "contact"
    """``contact`` sizes the sector to the arc that actually touches, which is
    the physically correct choice and the default. ``cap`` sizes it to whatever
    ``macro_grain_cap`` allows, which is cheaper and can be too flat to make
    contact at all -- ``plan`` reports ``spans_indent`` either way."""
    macro_sector_deg: float = 0.0
    """Sector of wheel to model. 0 -> the smallest sector that holds the cap.
    A full wheel is never right here: at 1,400 grains/mm2 a 125 mm wheel with a
    10 mm face carries 5.5 million grains."""
    mass_scale_factor: float = 0.0
    """Fixed mass scaling. 0 -> none. Reported either way; the reference used
    50 while also under-densifying the ring, which double-counts the speedup."""

    name: str = "sag"
    cores: int = 8

    def __post_init__(self):
        if self.compression_mm <= 0:
            raise SAGDeckError(
                "wheel compression must be positive: it is what creates the "
                "contact patch, and at T = 0 there is no SAG process")
        if self.speed_rpm <= 0:
            raise SAGDeckError("wheel speed must be positive")
        if self.elements_per_dc < 3.0:
            raise SAGDeckError(
                "elements_per_dc = %g cannot resolve a transition across dc. "
                "The project standard is 5." % self.elements_per_dc)
        if self.micro_grains < 1:
            raise SAGDeckError("the MICRO deck needs at least one grain")
        if self.macro_sector_mode not in ("contact", "cap"):
            raise SAGDeckError(
                "macro_sector_mode must be 'contact' or 'cap', not %r"
                % (self.macro_sector_mode,))
        if self.elements_per_grain < 6.0:
            raise SAGDeckError(
                "elements_per_grain = %g cannot carry a grain's shape in-plane"
                % self.elements_per_grain)
        if self.press_time_s > 0 and self.press_mach() > 0.05:
            raise SAGDeckError(
                "press_time_s = %g gives v/c = %.3f in the compliant layer. "
                "Above ~0.01 the contact pressure is inertial rather than "
                "Hertzian and the patch never reaches equilibrium; leave it 0 "
                "to have a quasi-static time derived."
                % (self.press_time_s, self.press_mach()))
        if self.hold_time_s > 0 and                 self.hold_time_s < self.polyurethane.prony_tau_s:
            raise SAGDeckError(
                "hold_time_s = %g is shorter than the Prony relaxation time "
                "%g, so the layer is still glassy when grinding starts and "
                "will carry more load than the relaxed experiment."
                % (self.hold_time_s, self.polyurethane.prony_tau_s))
        if self.compression_mm >= self.polyurethane.thickness_mm:
            raise SAGDeckError(
                "compression %.3f mm meets or exceeds the %.3f mm compliant "
                "layer, so the hub would contact the work and the tool is no "
                "longer compliant"
                % (self.compression_mm, self.polyurethane.thickness_mm))

    def wave_speed_mm_s(self) -> float:
        """Dilatational wave speed of the compliant layer, mm/s.

        Uses the small-strain modulus of the neo-Hookean card (E = 6*C10) at
        the layer's own density. This is what sets both the quasi-static press
        rate and the layer's share of the stable time increment.
        """
        e = self.polyurethane.modulus_mpa
        rho = self.polyurethane.density_tonne_mm3()
        return math.sqrt(e / rho)

    def press_time(self) -> float:
        if self.press_time_s > 0:
            return self.press_time_s
        return self.compression_mm / (PRESS_MACH * self.wave_speed_mm_s())

    def hold_time(self) -> float:
        if self.hold_time_s > 0:
            return self.hold_time_s
        return HOLD_TAUS * self.polyurethane.prony_tau_s

    def press_mach(self) -> float:
        """v/c of the press. Above ~0.01 the contact is inertial."""
        t = self.press_time()
        return (self.compression_mm / t) / self.wave_speed_mm_s() if t > 0             else float("inf")

    def wave_transits(self) -> float:
        """How many times a wave crosses the layer during the press."""
        c = self.wave_speed_mm_s()
        return self.press_time() * c / self.polyurethane.thickness_mm

    def tool(self) -> sag.Tool:
        pad = sag.Pad(self.grain_um, areal_per_mm2=self.pad_areal_per_mm2)
        return sag.Tool(
            diameter_mm=self.diameter_mm, width_mm=self.width_mm, pad=pad,
            shore_a=self.shore_a, poisson=0.24,
            elastic_mpa=(0.0 if self.use_shore_modulus
                         else self.polyurethane.modulus_mpa),
            layer_thickness_mm=self.polyurethane.thickness_mm)

    def hub_diameter(self) -> float:
        if self.hub_diameter_mm > 0:
            return self.hub_diameter_mm
        return self.diameter_mm - 2.0 * self.polyurethane.thickness_mm


def plan(p: SAGParams) -> dict:
    """Everything both decks will contain, and the cost, without writing a byte.

    Deliberately mirrors ``build_deck.plan_deck``: the notebook previews and the
    figures read this, so a preview cannot disagree with the deck that follows.
    """
    from . import materials

    w = materials.get(p.material)
    hp = w.hybrid_params()
    dc_mm = hp.critical_depth_mm()
    dc_nm = dc_mm * 1e6

    contact = sag.solve_contact(
        p.tool(), compression_mm=p.compression_mm, speed_rpm=p.speed_rpm,
        work_modulus_mpa=hp.youngs_mpa, work_poisson=hp.poisson,
        bhn_kgf_mm2=p.bhn_kgf_mm2, friction=p.friction)

    dc_forms = sag.dc_report(
        hardness_mpa=hp.hardness_mpa, youngs_mpa=hp.youngs_mpa,
        kic_mpa_sqrt_m=w.dc["kic_mpa_sqrt_m"],
        measured_nm=(sag.MEASURED_DC_WC_CO_NM if p.material == "wc_co"
                     else None),
        lambda_c=hp.lambda_c)

    k = sag.carbide_ratio(p.grain_um, p.carbide_um)

    # --- what the MACRO deck cannot do, quantified -----------------------
    el_dc_mm = dc_mm / p.elements_per_dc
    full_patch_elements = contact.spot_area_mm2 / (el_dc_mm ** 2)

    # --- MICRO patch sizing ----------------------------------------------
    # Side of a square holding the requested number of grains at the pad's own
    # active density, so the patch is a real sample of the contact rather than
    # an arbitrary box.
    side = p.micro_patch_mm or math.sqrt(
        p.micro_grains / p.tool().pad.active_per_mm2)
    depth = max(20.0 * dc_mm, 4.0 * contact.indentation_mm)
    n_micro = int(round(p.tool().pad.active_per_mm2 * side * side))

    # dc/5 is needed THROUGH THE DEPTH, not in-plane. The transition test is
    # h < dc, and both are depths -- so the direction that has to resolve dc is
    # the one h is measured along. In-plane the length that matters is the
    # grain, not dc: the mesh has to carry the shape of the cutting edge.
    #
    # This is not a shortcut, it is what the existing decks already do.
    # ARC_80NM resolves dc = 52.9 nm with 10.6 nm depth elements and 350 nm
    # in-plane -- a 33:1 anisotropy -- and it is the deck the whole project's
    # transition result rests on. Meshing in-plane at dc/5 as well would put
    # 5,282 elements across an 84 um patch whose grains are 27 um apart, i.e.
    # spend almost all of them on empty space between grains.
    el_depth_mm = el_dc_mm
    el_inplane_mm = (p.grain_um * 1e-3) / p.elements_per_grain
    surf_band = min(depth, 10.0 * dc_mm)
    nx = max(int(round(side / el_inplane_mm)), 1)
    ny = nx
    nz_fine = max(int(round(surf_band / el_depth_mm)), 1)
    nz_coarse = 12
    micro_elements = nx * ny * (nz_fine + nz_coarse)

    # --- MACRO sizing -----------------------------------------------------
    # A full wheel is out of the question: 5.5 M grains on a 6 um pad. So the
    # deck carries the smallest SECTOR that still holds the grain cap, which is
    # the same argument rim_depth_mm makes in the rigid pipeline -- model the
    # part that touches, not the part that exists.
    pad_area_full = math.pi * p.diameter_mm * p.width_mm
    dens = p.tool().pad.active_per_mm2
    grains_full = dens * pad_area_full
    # The sector has to be wide enough that its OWN CURVATURE spans the
    # indent. A sector of angle s has sagitta R(1 - cos(s/2)); if that is less
    # than the wheel compression T, the modelled arc is effectively flat and
    # the deck describes a flat punch pressed into the work, not a wheel. On a
    # 62.5 mm radius at T = 0.4 mm that needs s >= 12.97 deg.
    #
    # Worth noting: the paper's MEASURED spot length gives a contact arc of
    # 14.29 deg at the same setting, so the geometric requirement and the
    # experiment agree to 10% by two unrelated routes.
    r_out = 0.5 * p.diameter_mm
    ratio = 1.0 - p.compression_mm / r_out
    sector_min_deg = (2.0 * math.degrees(math.acos(max(-1.0, ratio)))
                      if ratio > -1.0 else 360.0)
    contact_arc_deg = 2.0 * math.degrees(math.asin(
        min(contact.spot_length_mm / p.diameter_mm, 1.0)))

    if p.macro_sector_deg > 0:
        sector_deg = min(p.macro_sector_deg, 360.0)
    elif p.macro_sector_mode == "contact":
        # The physically right answer: model the arc that actually touches,
        # with a little margin so the patch is not clipped by the cut faces.
        sector_deg = min(360.0, max(contact_arc_deg * 1.2, sector_min_deg))
    else:
        need_area = p.macro_grain_cap / dens if dens > 0 else pad_area_full
        sector_deg = min(360.0, 360.0 * need_area / pad_area_full)
    sector_area = pad_area_full * sector_deg / 360.0
    sector_grains = int(round(dens * sector_area))
    sector_arc_mm = math.pi * p.diameter_mm * sector_deg / 360.0

    # Under-populating the pad is NOT a harmless economy. The contact load is
    # shared among the grains that are actually there, so placing 5,000 of the
    # 261,877 a 17 deg sector really holds makes the force per grain 52x too
    # high -- and that force is the single number coupling MACRO to MICRO. So
    # a cap is only honoured by NARROWING THE SECTOR to match it, never by
    # thinning the pad within a sector.
    capped = sector_grains > p.macro_grain_cap
    if capped and p.macro_sector_deg <= 0:
        sector_deg = min(sector_deg,
                         360.0 * (p.macro_grain_cap / dens) / pad_area_full)
        sector_area = pad_area_full * sector_deg / 360.0
        sector_grains = int(round(dens * sector_area))
        sector_arc_mm = math.pi * p.diameter_mm * sector_deg / 360.0
    grains_placed = sector_grains
    density_honoured = abs(
        grains_placed / max(sector_area, 1e-30) - dens) / dens < 0.02

    return dict(
        params=p,
        material=dict(key=p.material, label=w.label,
                      dc_nm=dc_nm, dc_measured=w.dc_measured,
                      youngs_mpa=hp.youngs_mpa, hardness_mpa=hp.hardness_mpa),
        contact=contact,
        dc_forms=dc_forms,
        carbide=k,
        timing=dict(
            wave_speed_mm_s=p.wave_speed_mm_s(),
            press_time_s=p.press_time(), press_mach=p.press_mach(),
            wave_transits=p.wave_transits(),
            hold_time_s=p.hold_time(),
            hold_taus=p.hold_time() / p.polyurethane.prony_tau_s,
            grind_time_s=p.grind_time_s,
            total_time_s=p.press_time() + p.hold_time() + p.grind_time_s,
            press_velocity_mm_s=p.compression_mm / p.press_time(),
        ),
        macro=dict(
            hub_diameter_mm=p.hub_diameter(),
            layer_thickness_mm=p.polyurethane.thickness_mm,
            pu_density_kg_m3=p.polyurethane.density_kg_m3,
            grains_on_full_pad=int(round(grains_full)),
            grains_in_contact=int(round(contact.active_grains)),
            sector_deg=sector_deg, sector_area_mm2=sector_area,
            sector_arc_mm=sector_arc_mm,
            grains=grains_placed,
            grain_facets=int(round(grains_placed * FACETS_PER_GRAIN)),
            density_honoured=density_honoured,
            grains_per_mm2=grains_placed / max(sector_area, 1e-30),
            pad_density_per_mm2=dens,
            capped=capped, grain_cap=p.macro_grain_cap,
            contact_arc_deg=contact_arc_deg,
            sector_min_deg=sector_min_deg,
            sector_sagitta_mm=r_out * (1.0 - math.cos(
                math.radians(sector_deg) / 2.0)),
            spans_indent=(r_out * (1.0 - math.cos(math.radians(sector_deg)
                                                  / 2.0))
                          >= p.compression_mm),
            resolves_dc=False,
        ),
        micro=dict(
            side_mm=side, depth_mm=depth, grains=max(n_micro, 1),
            element_mm=el_depth_mm, element_inplane_mm=el_inplane_mm,
            elements=micro_elements, nx=nx, nz_fine=nz_fine,
            elements_per_dc=p.elements_per_dc,
            elements_per_grain=p.elements_per_grain,
            anisotropy=el_inplane_mm / el_depth_mm,
            load_per_grain_n=contact.load_per_grain_n,
            indentation_nm=contact.indentation_nm,
            resolves_dc=True,
        ),
        infeasible=dict(
            full_patch_elements=full_patch_elements,
            note=("resolving dc across the whole %.1f mm2 patch would need "
                  "%.2e elements" % (contact.spot_area_mm2,
                                     full_patch_elements)),
        ),
    )


def header_lines(pl: dict) -> list:
    """The block written at the top of both decks.

    Long, on purpose. Every number a reader would otherwise have to trust is
    stated with where it came from, and the two things that would silently
    invalidate a conclusion -- that MACRO cannot resolve dc, and that Bifano
    disagrees with the measurement on this material -- are stated as facts
    rather than left in a report file nobody opens.
    """
    p: SAGParams = pl["params"]
    c: sag.SAGContact = pl["contact"]
    m = pl["material"]
    d = pl["dc_forms"]
    k = pl["carbide"]
    L = []
    a = L.append
    a("** " + "=" * 74)
    a("** SHAPE ADAPTIVE GRINDING")
    a("**")
    a("** A compliant tool: rigid hub, %.1f mm polyurethane layer, %.0f um"
      % (p.polyurethane.thickness_mm, p.grain_um))
    a("** abrasive pad. The layer squashes by the wheel compression, line")
    a("** contact spreads into an area, and the load is shared by every grain")
    a("** that area covers -- which is why the force per grain collapses and a")
    a("** brittle material can be removed ductilely.")
    a("**")
    a("** THE PROCESS")
    a("**   wheel                 %.1f mm dia x %.1f mm wide"
      % (p.diameter_mm, p.width_mm))
    a("**   compression T         %.4f mm" % p.compression_mm)
    a("**   speed N               %.1f rpm  (%.1f mm/s surface)"
      % (p.speed_rpm, c.surface_speed_mm_s))
    a("**   pad                   %.0f um grain, %.0f active grains/mm2%s"
      % (p.grain_um, p.tool().pad.active_per_mm2,
         "" if c.density_measured else "  (density NOT measured, supplied)"))
    a("**   workpiece             %s" % m["label"])
    a("**")
    a("** THE CONTACT  (Ghosh, Sidpara & Bandyopadhyay 2021, eqs. 2-16)")
    a("**   tool modulus Et       %.4f MPa   (%s)"
      % (c.tool_modulus_mpa,
         "eq. 4 from shore A %.0f" % p.shore_a if p.use_shore_modulus
         else "6*C10 from the PU card"))
    a("**   equivalent Eeq        %.4f MPa   (eq. 3)" % c.eq_modulus_mpa)
    a("**   normal load FN        %.4f N      (eq. 2, ~T^1.5)" % c.normal_load_n)
    a("**   tangential FT         %.4f N      (eq. 5, mu = %.2f)"
      % (c.tangential_load_n, p.friction))
    a("**   spot area As          %.2f mm2    (eq. 6, empirical fit)"
      % c.spot_area_mm2)
    a("**   spot length Ls        %.3f mm     (eq. 7, empirical fit)"
      % c.spot_length_mm)
    a("**   ellipse semi-axes     a = %.3f mm, b = %.3f mm"
      % (c.semi_axis_a_mm, c.semi_axis_b_mm))
    a("**   peak pressure p0      %.5f MPa   (eq. 1)" % c.max_pressure_mpa)
    a("**   active grains Nabr    %.0f          (eq. 8)" % c.active_grains)
    a("**   load per grain Fn     %.4e N   (eq. 9)  <-- COUPLES THE TWO DECKS"
      % c.load_per_grain_n)
    a("**   indentation d         %.4f nm     (eqs. 11-12)" % c.indentation_nm)
    a("**   groove width          %.2f nm     (2*sqrt(d(dg-d)))"
      % c.groove_width_nm)
    a("**   MRR                   %.4f mm3/min (eq. 16)" % c.mrr_mm3_min)
    a("**")
    a("** THE TRANSITION")
    a("**   dc used               %.2f nm  (%s)"
      % (m["dc_nm"], "MEASURED" if m["dc_measured"] else "computed"))
    a("**   h / dc                %.4f  ->  %s"
      % (c.margin(m["dc_nm"]), c.regime(m["dc_nm"]).upper()))
    a("**   Bifano form 2         %.2f nm" % d["form2_nm"])
    a("**   form 1                %.2f nm" % d["form1_nm"])
    if "bifano_over_measured" in d:
        a("**")
        a("**   Bifano is %.1fx the measured value on this material, and that"
          % d["bifano_over_measured"])
        a("**   is the reference paper's own conclusion, not a slip here: a")
        a("**   thermally sprayed multiphase coating is not sintered bulk. The")
        a("**   deck uses the MEASUREMENT. Had it used Bifano's %.2f nm the"
          % d["form2_nm"])
        a("**   whole scratch would be ductile by construction and would")
        a("**   demonstrate nothing.")
    a("**   k = dg/D_WC           %.2f  (pure ductile needs < %.0f) -> %s"
      % (k["k"], k["threshold"], "yes" if k["pure_ductile"] else "no"))
    a("** " + "=" * 74)
    return L


def macro_header(pl: dict) -> list:
    L = header_lines(pl)
    mac = pl["macro"]
    inf = pl["infeasible"]
    t = pl["timing"]
    a = L.append
    a("**")
    a("** THIS IS THE MACRO DECK: THE CONTACT, NOT THE TRANSITION")
    a("**")
    a("** It carries the compliant tool -- rigid hub, %.1f mm polyurethane"
      % mac["layer_thickness_mm"])
    a("** ring, and %s measured grains." % format(mac["grains"], ","))
    a("**")
    a("** THREE STEPS, and the first two are timed by the layer's own physics:")
    a("**   1 PRESS  %.6f s   %.1f mm/s = %.4f of the layer's wave speed"
      % (t["press_time_s"], t["press_velocity_mm_s"], t["press_mach"]))
    a("**                        (%.1f wave transits of the %.1f mm layer, so"
      % (t["wave_transits"], mac["layer_thickness_mm"]))
    a("**                        the patch reaches equilibrium rather than")
    a("**                        being loaded inertially)")
    a("**   2 HOLD   %.6f s   %.1f Prony relaxation times, so the layer"
      % (t["hold_time_s"], t["hold_taus"]))
    a("**                        relaxes to its LONG-TERM modulus -- which is")
    a("**                        the state the measured force and the Hertz")
    a("**                        comparison both correspond to")
    a("**   3 GRIND  %.6f s   at %.1f rpm" % (t["grind_time_s"],
                                              pl["params"].speed_rpm))
    a("**   total    %.6f s" % t["total_time_s"])
    a("**")
    a("** A SECTOR, not the whole wheel: %.3f deg, %.3f mm of arc, %.3f mm2"
      % (mac["sector_deg"], mac["sector_arc_mm"], mac["sector_area_mm2"]))
    a("** of pad. The full %.0f mm wheel would carry %s grains, which is not"
      % (pl["params"].diameter_mm, format(mac["grains_on_full_pad"], ",")))
    a("** a deck. This is the same argument rim_depth_mm makes in the rigid")
    a("** pipeline: model the part that touches, not the part that exists.")
    a("**")
    a("** %s grains x %.0f facets = %s general-contact facets. For scale, the"
      % (format(mac["grains"], ","), FACETS_PER_GRAIN,
         format(mac["grain_facets"], ",")))
    a("** existing multi-abrasive deck runs 12 grains and 1,452 facets.")
    if mac["capped"]:
        a("**")
        a("** CAPPED at %s grains. The sector's own density would put more"
          % format(mac["grain_cap"], ","))
        a("** there; raise macro_grain_cap knowingly if the cost is acceptable.")
    a("**")
    a("** The real contact arc is %.2f deg (spot length %.3f mm), so the"
      % (mac["contact_arc_deg"], pl["contact"].spot_length_mm))
    a("** sector models %.1f%% of it -- a REPRESENTATIVE STRIP of the contact,"
      % min(100.0, 100.0 * mac["sector_deg"] / max(mac["contact_arc_deg"],
                                                   1e-9)))
    a("** which is what makes the pressure and per-grain load meaningful")
    a("** without carrying the whole patch.")
    a("**")
    a("** It CANNOT show a ductile-brittle transition, and no plot from it")
    a("** should be read as one. %s" % inf["note"])
    a("** For scale, the finest deck in this project is 172,000 elements and")
    a("** 74 h on 8 cores, for a single grain. What this deck gives is the")
    a("** contact: patch size, pressure distribution, engaged grain count and")
    a("** the load per grain -- which is what drives the MICRO deck.")
    a("**")
    a("** Polyurethane density %.1f kg/m3." % mac["pu_density_kg_m3"])
    if abs(mac["pu_density_kg_m3"] - PU_REFERENCE_DECK_DENSITY_KG_M3) < 0.1:
        a("** WARNING: that is the REFERENCE DECK's value, which is about")
        a("** 5000x too light -- lighter than air. It shortens the run by")
        a("** inflating the stable increment and removes the ring's inertia.")
        a("** Kept only to reproduce that deck; do not quote a force from it.")
    a("** " + "=" * 74)
    return L


def micro_header(pl: dict) -> list:
    L = header_lines(pl)
    mic = pl["micro"]
    a = L.append
    a("**")
    a("** THIS IS THE MICRO DECK: THE TRANSITION, RESOLVED")
    a("**")
    a("** A %.5f x %.5f mm patch of the contact above, %.5f mm deep,"
      % (mic["side_mm"], mic["side_mm"], mic["depth_mm"]))
    a("** holding %d measured grain(s) -- the number the pad's own active"
      % mic["grains"])
    a("** density puts in an area that size, so the patch is a sample of the")
    a("** real contact and not an arbitrary box.")
    a("**")
    a("** Each grain is driven by Fn = %.4e N, the load per grain the"
      % mic["load_per_grain_n"])
    a("** MACRO deck's contact solution gives. That single number is the")
    a("** coupling between the two decks; quoting one without the other is")
    a("** quoting half the model.")
    a("**")
    a("** Mesh %.3f nm through the depth = dc/%.1f, and %.0f nm in-plane"
      % (mic["element_mm"] * 1e6, mic["elements_per_dc"],
         mic["element_inplane_mm"] * 1e6))
    a("** = grain/%.0f -- a %.0f:1 anisotropy, deliberate. h and dc are both"
      % (mic["elements_per_grain"], mic["anisotropy"]))
    a("** DEPTHS, so the depth is what must resolve dc; in-plane the length")
    a("** that matters is the grain. ARC_80NM in this project does the same at")
    a("** 33:1 and is where the existing transition result comes from.")
    a("** About %s elements. SDV13 is the result: 1 ductile, 2 brittle."
      % format(mic["elements"], ","))
    a("** " + "=" * 74)
    return L


def summary_text(pl: dict) -> str:
    """Human-readable plan, for a notebook cell or a log."""
    return "\n".join(ln[3:] if ln.startswith("** ") else ln[2:]
                     for ln in macro_header(pl))


def demo() -> None:
    """Self-check: the plan must be arithmetically consistent and refuse abuse."""
    p = SAGParams(grain_um=6.0, compression_mm=0.4, speed_rpm=1050.0,
                  material="wc_co")
    pl = plan(p)

    c = pl["contact"]
    m = pl["material"]

    # the material is the paper's, with its MEASURED dc
    assert m["dc_measured"], "wc_co must use the measured dc"
    assert abs(m["dc_nm"] - 80.0) < 1e-9, m["dc_nm"]

    # and the 6 um pad must come out ductile, which is the paper's result
    assert c.regime(m["dc_nm"]) == "ductile"
    assert c.margin(m["dc_nm"]) < 1.0

    # Bifano must be reported and must disagree
    assert pl["dc_forms"]["bifano_over_measured"] > 10.0

    # k < 5 for the 6 um pad on 1.36 um carbide
    assert pl["carbide"]["pure_ductile"]

    # the MACRO deck must declare it cannot resolve dc, and the number that
    # justifies that must be genuinely enormous
    assert pl["macro"]["resolves_dc"] is False
    mac = pl["macro"]
    assert mac["sector_deg"] < 360.0, "MACRO must be a sector, not the wheel"
    assert mac["grains_on_full_pad"] > 1e6, mac["grains_on_full_pad"]
    assert mac["grains"] <= p.macro_grain_cap
    assert abs(mac["grain_facets"] - mac["grains"] * FACETS_PER_GRAIN) < 1.0

    # THE SECTOR MUST BE ABLE TO MAKE CONTACT. A sector of angle s has its own
    # sagitta R(1 - cos(s/2)); below the wheel compression the modelled arc is
    # effectively flat, and the deck describes a flat punch, not a wheel.
    assert mac["spans_indent"], (mac["sector_deg"], mac["sector_sagitta_mm"])
    assert mac["sector_sagitta_mm"] >= p.compression_mm
    assert mac["sector_deg"] >= mac["sector_min_deg"]
    # The geometric minimum must agree with the paper's MEASURED contact arc.
    # Two unrelated routes to the same angle, so this is a real cross-check on
    # the empirical fits rather than a restatement of them.
    assert abs(mac["sector_min_deg"] - mac["contact_arc_deg"]) \
        / mac["contact_arc_deg"] < 0.15, (mac["sector_min_deg"],
                                          mac["contact_arc_deg"])

    # THE PAD DENSITY MUST NEVER BE THINNED. A cap is honoured by narrowing the
    # sector; placing fewer grains than the sector holds would share the same
    # contact load among fewer carriers and make the per-grain force -- the one
    # number coupling MACRO to MICRO -- wrong in proportion. A 5,000-grain cap
    # on a 17 deg sector would have made it 52x too high.
    assert mac["density_honoured"], (mac["grains_per_mm2"],
                                     mac["pad_density_per_mm2"])
    assert abs(mac["grains_per_mm2"] - mac["pad_density_per_mm2"]) \
        / mac["pad_density_per_mm2"] < 0.02

    tight = plan(SAGParams(grain_um=6.0, macro_grain_cap=5000))
    assert tight["macro"]["grains"] <= 5000
    assert tight["macro"]["sector_deg"] < mac["sector_deg"], \
        "a cap must narrow the sector, not thin the pad"
    assert tight["macro"]["density_honoured"], "the pad was thinned"
    assert not tight["macro"]["spans_indent"], \
        "and a 5,000-grain sector is honestly reported as too flat to contact"

    # an explicit sector is honoured even when it cannot contact -- but said so
    fixed = plan(SAGParams(grain_um=6.0, macro_sector_deg=1.0,
                           macro_grain_cap=10 ** 9))
    assert abs(fixed["macro"]["sector_deg"] - 1.0) < 1e-9
    assert not fixed["macro"]["spans_indent"]

    # every pad size must reach a contacting sector at the default cap
    for dg in (6.0, 15.0, 30.0):
        q = plan(SAGParams(grain_um=dg))["macro"]
        assert q["spans_indent"], dg
        assert q["density_honoured"], dg

    try:
        SAGParams(macro_sector_mode="whatever")
    except SAGDeckError:
        pass
    else:
        raise AssertionError("an unknown sector mode must be refused")
    assert pl["infeasible"]["full_patch_elements"] > 1e9

    # the MICRO patch must be a real sample: the grains it claims must be what
    # the pad's density puts in that area
    mic = pl["micro"]
    dens = p.tool().pad.active_per_mm2
    assert abs(mic["grains"] - dens * mic["side_mm"] ** 2) < 1.0
    assert mic["resolves_dc"] is True
    # and it must actually resolve dc THROUGH THE DEPTH
    assert abs(mic["element_mm"] * 1e6 - m["dc_nm"] / p.elements_per_dc) < 1e-9
    # in-plane it resolves the GRAIN, and that is a much coarser length
    assert abs(mic["element_inplane_mm"] * 1e6
               - p.grain_um * 1000.0 / p.elements_per_grain) < 1e-6
    assert mic["anisotropy"] > 5.0, mic["anisotropy"]
    # tractable, unlike the full patch
    assert mic["elements"] < 5e6, mic["elements"]
    assert pl["infeasible"]["full_patch_elements"] / mic["elements"] > 1e4

    # coarser pad -> deeper bite -> closer to brittle (the paper's trend)
    coarse = plan(SAGParams(grain_um=30.0, compression_mm=0.4,
                            speed_rpm=1050.0, material="wc_co"))
    assert (coarse["contact"].indentation_nm
            > pl["contact"].indentation_nm)
    assert not coarse["carbide"]["pure_ductile"], "30 um on 1.36 um is k=22"

    # more compression -> more load (eq. 2)
    hard = plan(SAGParams(grain_um=6.0, compression_mm=0.6, speed_rpm=1050.0,
                          material="wc_co"))
    assert hard["contact"].normal_load_n > c.normal_load_n

    # the PU card writes 8-per-line-safe blocks and real density by default
    pu = Polyurethane()
    assert pu.density_kg_m3 == PU_SOLID_DENSITY_KG_M3
    cards = pu.cards()
    assert any("neo hooke" in x for x in cards)
    assert any("PRONY" in x for x in cards)
    assert abs(pu.modulus_mpa - 0.345) < 1e-9
    assert abs(pu.long_term_fraction() - 0.89) < 1e-9

    # headers must state the two things that would invalidate a conclusion
    mh = "\n".join(macro_header(pl))
    assert "CANNOT show a ductile-brittle transition" in mh
    assert "A SECTOR, not the whole wheel" in mh
    assert "Bifano is" in mh
    ih = "\n".join(micro_header(pl))
    assert "SDV13" in ih
    assert "COUPLES THE TWO DECKS" in "\n".join(header_lines(pl))

    # the reference deck's density must be flagged when it is used
    ref = plan(SAGParams(
        grain_um=6.0,
        polyurethane=Polyurethane(
            density_kg_m3=PU_REFERENCE_DECK_DENSITY_KG_M3)))
    assert "5000x too light" in "\n".join(macro_header(ref))

    # --- refusals --------------------------------------------------------
    for bad, why in (
        (lambda: SAGParams(compression_mm=0.0), "zero compression"),
        (lambda: SAGParams(speed_rpm=0.0), "zero speed"),
        (lambda: SAGParams(elements_per_dc=2.0), "cannot resolve dc"),
        (lambda: SAGParams(compression_mm=6.0), "compression exceeds layer"),
        (lambda: SAGParams(micro_grains=0), "no grains"),
        (lambda: SAGParams(elements_per_grain=3.0), "grain unresolved"),
        (lambda: Polyurethane(c10_mpa=0.0), "zero C10"),
        (lambda: Polyurethane(prony_g=1.0), "fully relaxing Prony"),
        (lambda: Polyurethane(density_kg_m3=0.0), "zero density"),
        (lambda: SAGParams(grain_um=7.5).tool(), "unmeasured pad density"),
    ):
        try:
            bad()
        except (SAGDeckError, sag.SAGError):
            pass
        else:
            raise AssertionError("should have been refused: %s" % why)

    print("semgrit.sagdeck: all checks passed")
    print(summary_text(pl))


if __name__ == "__main__":
    demo()
