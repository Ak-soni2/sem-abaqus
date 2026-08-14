"""Turn the geometry deck into one you can submit straight from the terminal.

Everything that was being added by hand in Abaqus/CAE -- step, mass scaling, the wheel
motion, contact, the JH-2 user material, section controls, restart and output -- is
written here instead, so the ``.inp`` the notebook produces needs no GUI at all:

    abaqus job=grind input=wheel.inp user=vumat_jh2_3.for double=both cpus=8 interactive

Two details this module exists to get right, both of which cost real time when they
were done by hand:

**The infeed.** A wheel with only a rotation BC spins on its own axis and never touches
anything: at t = 0 one grit is tangent, it turns away, and every grit behind it sits a
micron or two below the surface for ever. Grinding needs the wheel fed *into* the work.
The radial direction depends on where the block ended up, so the infeed components can
only be computed here, where ``theta_workpiece`` is known.

**The deletion flag.** Brittle fracture needs three things switched on together --
``*Depvar, delete=N``, ``ELEMENT DELETION=YES`` on the section controls, *and* a VUMAT
that actually drives the flag to zero. Any one missing and damaged material stays in
the mesh carrying residual strength, and the result looks ductile.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Sequence

# Baranowski et al., Sci. Rep. 14, 3339 (2024), in mm-MPa-tonne-s.
# Order is the one both supplied VUMATs read:
#   K1 G HEL PHEL T A B C N M beta D1 D2 K2 K3 SFMAX SIGHEL
JH2_SANDSTONE = (3735.6, 2686.0, 1982.0, 1374.0, 8.0, 0.71, 0.30, 0.022,
                 0.55, 0.40, 1.0, 0.002, 1.20, 9000.0, 22000.0, 0.25, 912.0)


@dataclass
class AnalysisParams:
    """Everything needed to make the deck submittable. All of it user-facing."""

    enabled: bool = False

    # ---- step -----------------------------------------------------------
    step_time_s: float = 0.0        # 0 = derive from travel / surface speed
    nlgeom: bool = True
    mass_scaling_factor: float = 10.0
    """Explicit runs faster with scaled mass but the physics degrades: the factor
    multiplies density, so 30 makes the workpiece 30x heavier and badly distorts
    inertia at grinding speeds. 10 is a reasonable compromise, 1 disables it."""
    bulk_viscosity: tuple = (0.06, 1.2)

    # ---- wheel motion ---------------------------------------------------
    depth_of_cut_um: float = 3.0
    """Radial infeed over the step. Without this nothing touches. Must stay under the
    bond-rim clearance reported by the build, or the bond hits the workpiece."""
    # Note: the infeed is a constant velocity, so the depth always ramps linearly
    # from zero to depth_of_cut_um across the step. Engagement therefore builds up,
    # which both animates better and avoids a step change in overclosure. Holding a
    # constant depth instead would need an *Amplitude; it is deliberately not offered
    # rather than offered as a switch that does nothing.

    # ---- workpiece material ---------------------------------------------
    material_model: str = "jh2"      # 'jh2' | 'elastic' | 'hybrid'
    """``hybrid`` is the ductile/brittle law of vumat_grind.for: Johnson-Cook
    with strain-gradient enhancement where the grit takes a cut thinner than
    the critical depth, JH-2 where it takes a thicker one. It keeps the same
    ``jh2_constants`` for the brittle half, so a hybrid deck is the brittle
    deck plus a ductile regime, not a different model."""
    jh2_constants: Sequence[float] = JH2_SANDSTONE
    jh2_density_kg_m3: float = 2350.0
    hybrid: Optional["object"] = None
    """A :class:`semgrit.hybrid.HybridParams`. Required by 'hybrid'."""
    n_depvar: int = 12
    element_deletion: bool = True
    hourglass: str = "ENHANCED"

    # ---- contact --------------------------------------------------------
    contact_scope: str = "engaging"  # 'engaging' | 'all exterior' | 'none'
    friction: float = 0.2

    # ---- how the block is held ------------------------------------------
    fix_back_face: bool = True
    fix_ends: bool = False
    fix_sides: bool = False

    # ---- output ---------------------------------------------------------
    field_frames: int = 60
    restart_intervals: int = 10
    """More than 1 matters: with 1 the only restart state is written at the end of the
    step, so an interrupted run cannot be resumed at all."""
    element_output: str = "S, PEEQ, SDV, STATUS"
    node_output: str = "U, V"
    history_preselect: bool = True
    history_reference_node: bool = True
    """Write RF and RM at the wheel's reference node to the history output.

    PRESELECT gives the energies but not this, and without it the .odb holds no
    grinding force at all -- the post-processing script would have nothing to read."""
    history_intervals: int = 2000
    """Samples of the reaction force over the step. 2000, not 200, and this is the
    one output setting that cannot be recovered after the run.

    A grit crosses one 0.3 um element in 0.3e-3/30000 = 1.0e-8 s. At 200 samples
    over an 1.8e-6 s step the interval is 9.0e-9 s -- 1.11 samples per element
    passage, i.e. exactly Nyquist on the dominant excitation of the whole model.
    The peak force would be a random draw somewhere below the true peak and the
    Ft/Fn ratio would inherit the aliasing noise. The solver increment is
    5.4e-11 s, so 2000 samples is still 186 increments apart, and one node at
    2000 samples is a few hundred kB."""
    rotation_reversed: bool = False
    """Turn the wheel the other way: the surface then travels toward *increasing* theta
    and grains arrive from the block's low-theta end.

    The default sense is the one the placement code assumes -- 'first grit at entry'
    takes the highest theta because grains arrive from there. Reversing this without
    reversing that would put the block at the exit, so the deck header and the entry
    edge are both derived from this one flag."""

    def effective_contact_scope(self) -> str:
        """``contact_scope``, promoted to 'all exterior' when elements can die.

        'engaging' pairs the grits against ``A_WP_GROUND_SURF``, which
        ``rigid_wheel`` builds as ``ES_WP_GROUND, S1`` -- face S1 of the FIRST
        element in each depth column, i.e. the top layer and nothing else. That
        is correct for a deck with no deletion and silently fatal for one with
        it: the cut is 4-5 of the 0.03 um layers deep, so once the top layer
        erodes the grit is cutting through a hole with no contact resistance at
        all. The job still completes, the .sta still says success, and every
        force is the response of a 0.03 um membrane.

        Abaqus cannot re-expose interior faces of a face-based surface as
        elements die; general contact with ALL EXTERIOR can, because it tracks
        the exterior of the mesh rather than a fixed face list. 24,960 hexes
        plus 4,168 rigid facets is affordable for that.

        The promotion is unconditional rather than a user choice because there
        is no combination of the two in which 'engaging' is right.
        """
        if self.contact_scope == "engaging" and self.element_deletion:
            return "all exterior"
        return self.contact_scope

    def validate(self) -> None:
        if self.material_model not in ("jh2", "elastic", "hybrid"):
            raise ValueError(
                "material_model must be 'jh2', 'elastic' or 'hybrid'")
        if self.contact_scope not in ("engaging", "all exterior", "none"):
            raise ValueError("contact_scope must be 'engaging', 'all exterior' or 'none'")
        if self.material_model in ("jh2", "hybrid") and len(self.jh2_constants) != 17:
            raise ValueError("the JH-2 VUMATs read props 1..17; got %d constants"
                             % len(self.jh2_constants))
        if self.material_model == "hybrid":
            from .hybrid import HYBRID_DEPVAR, HybridParams
            if not isinstance(self.hybrid, HybridParams):
                raise ValueError(
                    "material_model='hybrid' needs analysis.hybrid to be a "
                    "HybridParams; see semgrit.hybrid")
            self.hybrid.validate()
            if self.n_depvar < HYBRID_DEPVAR:
                raise ValueError(
                    "vumat_grind.for writes %d state variables (SDV13 is the "
                    "branch, SDV14 the chip thickness, SDV19 the strain-"
                    "gradient factor); n_depvar is only %d, so the ones that "
                    "say what the model did would be silently dropped"
                    % (HYBRID_DEPVAR, self.n_depvar))
        if self.element_deletion and self.n_depvar < 12:
            raise ValueError("element deletion uses SDV12 as the flag, so n_depvar "
                             "must be at least 12")
        if self.mass_scaling_factor < 1.0:
            raise ValueError("mass_scaling_factor must be >= 1")
        if self.field_frames < 1 or self.restart_intervals < 1:
            raise ValueError("field_frames and restart_intervals must be >= 1")


def wheel_motion(p: AnalysisParams, theta_c_rad: float, surface_speed_mm_s: float,
                 radius_mm: float, step_time_s: float) -> dict:
    """Rotation and infeed for the reference node, in the global frame."""
    omega = surface_speed_mm_s / radius_mm
    ae = p.depth_of_cut_um / 1000.0
    v_r = ae / step_time_s if step_time_s > 0 else 0.0
    # e_r points from the wheel axis out through the workpiece, and the block sits
    # *outside* the rim with its ground face looking back at the axis. So feeding in
    # means translating the wheel along +e_r: that is what closes the gap.
    #
    # This was the other way round, and it is why nothing ever cut. The geometry was
    # right -- the tallest grain under the block exactly tangent at t=0, verified to a
    # picometre -- and then the step commanded the wheel to retract by precisely the
    # depth of cut, so the one grain that touched let go on the first increment. The
    # job ran to completion, reported no error, and ground nothing. Nothing checked the
    # *direction* of the commanded motion until now; verify_rigid_deck.py does.
    return dict(
        omega_rad_s=omega, rpm=omega * 30.0 / math.pi,
        # Negative VR3 turns the surface toward decreasing theta. A positive rotation
        # about +Z carries +X toward +Y, so the sign and the sentence must be read
        # together -- the deck header claimed the opposite for months.
        vr3=(omega if p.rotation_reversed else -omega),
        v1=math.cos(theta_c_rad) * v_r,
        v2=math.sin(theta_c_rad) * v_r,
        radial_speed_mm_s=v_r, depth_of_cut_mm=ae,
        sweep_mm=omega * step_time_s * radius_mm,
    )


def write_section_and_material(w, p: AnalysisParams, wp, chip_field=None):
    """Section controls plus the workpiece material, written after *End Assembly.

    ``chip_field`` is a :class:`semgrit.hybrid.ChipField` and is only used by
    the hybrid law, which needs to know where along the scratch the grit takes
    a thick cut and where it takes a thin one. Returns the hybrid summary dict
    when that law is written, otherwise None.
    """
    if p.material_model == "elastic":
        w("*Material, name=%s\n" % wp.material)
        w("*Density\n%.8e,\n" % (wp.density_kg_m3 * 1e-12))
        w("*Elastic\n%.8e, %.4f\n" % (wp.youngs_modulus_mpa, wp.poisson_ratio))
        return None
    if p.material_model == "hybrid":
        from .hybrid import write_hybrid_material
        return write_hybrid_material(w, p.hybrid, wp, list(p.jh2_constants),
                                     chip_field, p.n_depvar,
                                     p.element_deletion)
    w("** Johnson-Holmquist II, driven by a VUMAT. Supply it on the command line:\n")
    w("**   abaqus job=... input=... user=vumat_jh2_3.for double=both cpus=N\n")
    w("** The VUMAT must set the deletion flag SDV%d to 0 once D reaches 1, or\n"
      % p.n_depvar)
    w("** nothing is ever deleted and the result looks ductile rather than brittle.\n")
    w("*Material, name=%s\n" % wp.material)
    w("*Density\n%.8e,\n" % (p.jh2_density_kg_m3 * 1e-12))
    if p.element_deletion:
        w("*Depvar, delete=%d\n%d,\n" % (p.n_depvar, p.n_depvar))
    else:
        w("*Depvar\n%d,\n" % p.n_depvar)
    c = list(p.jh2_constants)
    w("*User Material, constants=%d\n" % len(c))
    for i in range(0, len(c), 8):
        w(", ".join("%g" % v for v in c[i:i + 8]) + "\n")
    return None


def write_step(w, p: AnalysisParams, wp, motion: dict, step_time_s: float,
               engage_surface: bool) -> None:
    """The whole history definition, ready to submit."""
    w("**\n** ---------------- ANALYSIS ----------------\n**\n")
    w("*Section Controls, name=EC-1, hourglass=%s, ELEMENT DELETION=%s\n"
      % (p.hourglass, "YES" if p.element_deletion else "NO"))
    w("1., 1., 1.\n")
    w("*Surface Interaction, name=IntProp-1\n")
    if p.friction > 0:
        w("*Friction\n%g,\n" % p.friction)
    else:
        w("*Friction\n0.,\n")

    # How the block is held. Applied outside the step so it holds from t = 0.
    holds = [("A_WP_BACK_FACE", p.fix_back_face),
             ("A_WP_END_A", p.fix_ends), ("A_WP_END_B", p.fix_ends),
             ("A_WP_SIDE_A", p.fix_sides), ("A_WP_SIDE_B", p.fix_sides)]
    if any(on for _, on in holds):
        w("**\n*Boundary\n")
        for nm, on in holds:
            if on:
                w("%s, ENCASTRE\n" % nm)

    w("**\n*Step, name=Step-1, nlgeom=%s\n" % ("YES" if p.nlgeom else "NO"))
    w("*Dynamic, Explicit\n, %g\n" % step_time_s)
    w("*Bulk Viscosity\n%g, %g\n" % p.bulk_viscosity)
    if p.mass_scaling_factor > 1.0:
        w("*Fixed Mass Scaling, factor=%g\n" % p.mass_scaling_factor)

    w("**\n** The wheel: one velocity BC on the single rigid-body reference node.\n")
    w("** DOF 1 and 2 are the radial infeed -- the depth of cut. Without them the\n")
    w("** wheel spins on the spot and never touches the workpiece.\n")
    w("*Boundary, type=VELOCITY\n")
    w("A_WHEEL_REF, 1, 1, %.6f\n" % motion["v1"])
    w("A_WHEEL_REF, 2, 2, %.6f\n" % motion["v2"])
    w("A_WHEEL_REF, 3, 3\n")
    w("A_WHEEL_REF, 4, 4\n")
    w("A_WHEEL_REF, 5, 5\n")
    w("A_WHEEL_REF, 6, 6, %.6f\n" % motion["vr3"])

    scope = p.effective_contact_scope()
    if scope != "none":
        w("**\n*Contact, op=NEW\n")
        if scope != p.contact_scope:
            w("** contact scope promoted %r -> %r: element deletion is on, and a\n"
              % (p.contact_scope, scope))
            w("** face-based surface cannot re-expose interior faces as elements\n")
            w("** die. See AnalysisParams.effective_contact_scope.\n")
        if scope == "all exterior":
            w("*Contact Inclusions, ALL EXTERIOR\n")
        else:
            # Only the grits that can reach the block. On a wheel with half a million
            # facets this is far cheaper than tracking all of them.
            w("*Contact Inclusions\n")
            w("%s, A_WP_GROUND_SURF\n"
              % ("A_GRITS_ENGAGE_SURF" if engage_surface else "A_GRITS_SURF"))
        w("*Contact Property Assignment\n ,  , IntProp-1\n")

    w("**\n*Restart, write, number interval=%d, time marks=NO\n" % p.restart_intervals)
    w("*Output, field, number interval=%d\n" % p.field_frames)
    w("*Node Output\n%s\n" % p.node_output)
    # Scoped to the workpiece: it is the only deformable part, and an unscoped request
    # asks the rigid facets and the mass elements for stress and damage they cannot
    # have. Abaqus only warns, but twelve warnings hide the one that matters.
    w("*Element Output, elset=A_WP_ALL, directions=YES\n%s\n" % p.element_output)
    # One sampling interval for BOTH history blocks. The energy balance is what
    # decides whether the run is trustworthy (ALLAE/ALLIE, ALLKE/ALLIE), so
    # sampling it at Abaqus' 20-interval default while the force is sampled 2000
    # times means the two cannot be read against each other.
    n_hist = max(int(p.history_intervals), 1)
    dt_out = (step_time_s / n_hist) if step_time_s > 0 else 0.0
    hist_suffix = ("" if dt_out <= 0 else ", time interval=%.9e" % dt_out)
    if p.history_preselect:
        w("*Output, history, variable=PRESELECT%s\n" % hist_suffix)
    if p.history_reference_node:
        # The reaction at the driven reference node is the grinding force: the wheel
        # is one rigid body, so whatever the work pushes back with is reacted here.
        # PRESELECT does not include it, and without it the .odb holds no force at
        # all -- the job finishes and there is nothing to plot.
        #
        # NUMBER INTERVAL is a *field*-output parameter; on history output Abaqus
        # rejects it outright ("THE PARAMETER HISTORY CANNOT BE USED WITH THE
        # PARAMETER NUMBER INTERVAL"), which is fatal at the pre-processing stage. The
        # history equivalent is TIME INTERVAL, so convert the count into one.
        w("*Output, history%s\n" % hist_suffix)
        w("*Node Output, nset=A_WHEEL_REF\n")
        w("RF1, RF2, RF3, RM1, RM2, RM3, U1, U2, UR3\n")
    w("*End Step\n")
