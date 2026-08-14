"""A runnable Abaqus/Explicit multi-grit scratch simulation on a stone workpiece.

Units: **mm, tonne, s, MPa, N** throughout (density in tonne/mm^3).

Why micro-scale
---------------
Resolving a 3.4 um grit cutting into stone needs elements smaller than the grit,
~1.5 um. The real contact arc at 10 um depth of cut is 1.41 mm long, which at that
element size is order 10^8 elements -- far out of reach. So this builds the standard
alternative from the grinding literature: a **multi-grit scratch model** over a
patch of the wheel a few tens of microns across, with the grits swept across a
small stone block.

Simplifications, each stated so they can be judged:

* Grits are **rigid** (R3D3) and driven by a constant translational velocity rather
  than rotated about the wheel axis. Over the few microseconds simulated the wheel
  turns well under a tenth of a degree, so rotation and translation are
  indistinguishable at this scale.
* The **bond is omitted**. With rigid grits it carries no load and only adds
  elements; grit pull-out is outside the scope of a cutting check.
* The workpiece is a flat block. Wheel curvature over a 0.1 mm chord on a 100 mm
  radius is 0.0125 um, far below the element size.
* Side and bottom faces are **fixed**, so stress waves reflect. Over 4 us a wave
  travels ~17 mm in stone against a 0.1 mm model, so reflections are present. Fine
  for verifying mechanics; use infinite elements or a larger domain for quantitative
  force predictions.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Optional, TextIO

import numpy as np

from .abaqus import KGM3_TO_TONNE_MM3, MATERIALS, _write_int_set
from .wheel import UM_PER_MM, WheelModel, _rotation_matrix


@dataclass
class StoneMaterial:
    """Workpiece material.

    Carries both a Johnson-Holmquist II parameter set and a simpler
    pressure-dependent plasticity model. JH-2 is the intended constitutive law for
    brittle rock, but in Abaqus/Explicit it is normally supplied through a **VUMAT**
    user subroutine, so the deck cannot activate it on its own. The simpler model is
    written active so the deck runs immediately and the mechanics can be checked;
    the JH-2 constants are written alongside, ready to switch on.
    """

    name: str = "STONE_GRANITE"
    density_kg_m3: float = 2650.0
    youngs_modulus_mpa: float = 50_000.0
    poisson_ratio: float = 0.25

    # --- Drucker-Prager, used for the runnable version -------------------
    friction_angle_deg: float = 51.0
    dilation_angle_deg: float = 15.0
    compressive_strength_mpa: float = 150.0
    failure_strain: float = 0.10

    # --- Johnson-Holmquist II constants (representative granite) ---------
    # Ordered as a VUMAT would normally expect them. These are literature-typical
    # values for granite and MUST be checked against your own source before any
    # quantitative claim.
    jh2: dict = field(default_factory=lambda: {
        "G_shear_modulus_MPa": 20_000.0,
        "A_intact_strength": 0.95,
        "N_intact_exponent": 0.62,
        "B_fractured_strength": 0.35,
        "M_fractured_exponent": 0.62,
        "C_strain_rate": 0.005,
        "D1_damage": 0.005,
        "D2_damage": 0.70,
        "K1_bulk_MPa": 25_700.0,
        "K2_MPa": -4_500.0,
        "K3_MPa": 300_000.0,
        "T_tensile_strength_MPa": 15.0,
        "HEL_MPa": 4_500.0,
        "P_HEL_MPa": 2_930.0,
        "SIGMA_HEL_MPa": 2_350.0,
        "BETA_bulking": 1.0,
        "EPS_ref_strain_rate": 1.0,
    })

    @property
    def density_tonne_mm3(self) -> float:
        return self.density_kg_m3 * KGM3_TO_TONNE_MM3

    @property
    def wave_speed_mm_s(self) -> float:
        return math.sqrt(self.youngs_modulus_mpa / self.density_tonne_mm3)


@dataclass
class WorkpieceSpec:
    """Stone block geometry and mesh."""

    length_mm: float = 0.10      # along the grinding (grit travel) direction
    width_mm: float = 0.06       # along the wheel axis
    depth_mm: float = 0.025      # into the workpiece
    element_size_mm: float = 0.0015

    def divisions(self) -> tuple[int, int, int]:
        h = self.element_size_mm
        return (
            max(int(round(self.length_mm / h)), 1),
            max(int(round(self.width_mm / h)), 1),
            max(int(round(self.depth_mm / h)), 1),
        )

    def n_elements(self) -> int:
        a, b, c = self.divisions()
        return a * b * c


@dataclass
class ProcessSpec:
    """Grinding process parameters."""

    wheel_speed_m_s: float = 30.0
    engagement_um: float = 2.0
    """How deep the tallest grit tips sit below the workpiece top surface."""
    approach_gap_um: float = 5.0
    """Clearance between the leading grit and the workpiece at t = 0, so nothing
    starts already interpenetrating."""
    duration_s: Optional[float] = None
    """Defaults to the time for the leading grit to cross the whole workpiece."""
    mass_scaling_target_dt: Optional[float] = None

    @property
    def wheel_speed_mm_s(self) -> float:
        return self.wheel_speed_m_s * 1000.0


# --------------------------------------------------------------------------
# Workpiece mesh
# --------------------------------------------------------------------------

def build_workpiece_mesh(
    wp: WorkpieceSpec,
    origin: np.ndarray,
    e_travel: np.ndarray,
    e_axial: np.ndarray,
    e_up: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """Structured C3D8R box.

    ``origin`` is the centre of the top face; ``e_up`` points out of that face
    (radially outward, toward the wheel). The block extends ``depth_mm`` in ``-e_up``.
    """
    nl, nw, nd = wp.divisions()
    u = np.linspace(-wp.length_mm / 2.0, wp.length_mm / 2.0, nl + 1)
    v = np.linspace(-wp.width_mm / 2.0, wp.width_mm / 2.0, nw + 1)
    w = np.linspace(-wp.depth_mm, 0.0, nd + 1)

    def nid(i: int, j: int, k: int) -> int:
        return (i * (nw + 1) + j) * (nd + 1) + k

    nodes = np.empty(((nl + 1) * (nw + 1) * (nd + 1), 3), dtype=np.float64)
    for i, ui in enumerate(u):
        for j, vj in enumerate(v):
            for k, wk in enumerate(w):
                nodes[nid(i, j, k)] = origin + e_travel * ui + e_axial * vj + e_up * wk

    hexes: list[tuple[int, ...]] = []
    for i in range(nl):
        for j in range(nw):
            for k in range(nd):
                hexes.append((
                    nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                    nid(i, j, k + 1), nid(i + 1, j, k + 1),
                    nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1),
                ))

    sets = {
        "WP_BOTTOM": np.array([nid(i, j, 0) for i in range(nl + 1) for j in range(nw + 1)]),
        "WP_TOP": np.array([nid(i, j, nd) for i in range(nl + 1) for j in range(nw + 1)]),
        "WP_ENTRY": np.array([nid(0, j, k) for j in range(nw + 1) for k in range(nd + 1)]),
        "WP_EXIT": np.array([nid(nl, j, k) for j in range(nw + 1) for k in range(nd + 1)]),
        "WP_SIDE_A": np.array([nid(i, 0, k) for i in range(nl + 1) for k in range(nd + 1)]),
        "WP_SIDE_B": np.array([nid(i, nw, k) for i in range(nl + 1) for k in range(nd + 1)]),
    }
    return nodes, np.asarray(hexes, dtype=np.int64), sets


# --------------------------------------------------------------------------
# Deck writer
# --------------------------------------------------------------------------

def write_grinding_sim(
    path: str,
    model: WheelModel,
    workpiece: Optional[WorkpieceSpec] = None,
    process: Optional[ProcessSpec] = None,
    stone: Optional[StoneMaterial] = None,
    grain_material: str = "diamond",
) -> dict:
    """Write a complete, runnable Abaqus/Explicit grinding deck."""
    wp = workpiece or WorkpieceSpec()
    pr = process or ProcessSpec()
    st = stone or StoneMaterial()
    if not model.placements:
        raise ValueError("the wheel model has no grains placed")

    R = model.spec.outer_radius_mm

    # Local frame at the middle of the grit patch.
    thetas = [p.theta_deg for p in model.placements]
    arc_lo = math.radians(min(thetas)) * R
    arc_hi = math.radians(max(thetas)) * R
    grit_reach = max(
        float(np.hypot(
            *(( (s.vertices - s.centroid_um) / UM_PER_MM
                @ _rotation_matrix(p.rotation_axis, math.radians(p.rotation_angle_deg)).T
                + p.translation_mm)[:, :2].T)).max())
        for p in model.placements
        for s in (model.shapes[p.shape_index],)
    )
    max_protrusion_mm = grit_reach - R

    gap = pr.approach_gap_um / 1000.0
    engage = pr.engagement_um / 1000.0

    # Workpiece sits just ahead of the leading grit, along the travel direction.
    s_leading = arc_hi
    s_centre = s_leading + gap + wp.length_mm / 2.0
    theta_w = s_centre / R

    e_up = np.array([math.cos(theta_w), math.sin(theta_w), 0.0])
    e_travel = np.array([-math.sin(theta_w), math.cos(theta_w), 0.0])
    e_axial = np.array([0.0, 0.0, 1.0])

    # Top face placed so the tallest grit tips bite `engage` deep.
    r_top = R + max_protrusion_mm - engage
    origin = e_up * r_top

    wp_nodes, wp_hexes, wp_sets = build_workpiece_mesh(
        wp, origin, e_travel, e_axial, e_up
    )

    travel = gap + wp.length_mm
    duration = pr.duration_s if pr.duration_s else travel / pr.wheel_speed_mm_s
    v_vec = e_travel * pr.wheel_speed_mm_s

    dt_stable = wp.element_size_mm / st.wave_speed_mm_s
    n_inc = duration / dt_stable

    used_shapes = sorted({p.shape_index for p in model.placements})
    shape_verts = [(s.vertices - s.centroid_um) / UM_PER_MM for s in model.shapes]

    with open(path, "w", encoding="ascii", newline="\n") as fh:
        _header(fh, model, wp, pr, st, max_protrusion_mm, r_top, duration,
                dt_stable, n_inc, theta_w)

        # ---- grit parts (rigid surfaces) ----
        for idx in used_shapes:
            s = model.shapes[idx]
            fh.write(f"*Part, name=GRAIN-{idx + 1}\n*Node\n")
            for i, v in enumerate(shape_verts[idx], start=1):
                fh.write(f"{i}, {v[0]: .8e}, {v[1]: .8e}, {v[2]: .8e}\n")
            fh.write("*Element, type=R3D3\n")
            for e, tri in enumerate(s.faces, start=1):
                fh.write(f"{e}, {tri[0] + 1}, {tri[1] + 1}, {tri[2] + 1}\n")
            fh.write(f"*Elset, elset=GRAIN_ALL, generate\n1, {len(s.faces)}, 1\n")
            ref = len(shape_verts[idx]) + 1
            fh.write(f"*Node\n{ref}, 0., 0., 0.\n")
            fh.write(f"*Nset, nset=GRAIN_REF\n{ref},\n")
            fh.write("*Rigid Body, ref node=GRAIN_REF, elset=GRAIN_ALL\n")
            fh.write("*Surface, type=ELEMENT, name=GRAIN_SURF\nGRAIN_ALL, SPOS\n")
            fh.write("*End Part\n")

        # ---- workpiece part ----
        fh.write("*Part, name=WORKPIECE\n*Node\n")
        for i, v in enumerate(wp_nodes, start=1):
            fh.write(f"{i}, {v[0]: .8e}, {v[1]: .8e}, {v[2]: .8e}\n")
        fh.write("*Element, type=C3D8R\n")
        for e, h in enumerate(wp_hexes, start=1):
            fh.write(f"{e}, " + ", ".join(str(int(n) + 1) for n in h) + "\n")
        fh.write(f"*Elset, elset=WP_ALL, generate\n1, {len(wp_hexes)}, 1\n")
        for nm, ids in wp_sets.items():
            fh.write(f"*Nset, nset={nm}\n")
            _write_int_set(fh, [int(i) + 1 for i in ids])
        fh.write("*Surface, type=ELEMENT, name=WP_SURF\nWP_ALL,\n")
        # Element deletion needs a section control with distortion control off.
        fh.write("*Solid Section, elset=WP_ALL, material="
                 f"{st.name}, controls=EC-1\n,\n")
        fh.write("*End Part\n")

        # ---- assembly ----
        fh.write("*Assembly, name=ASSEMBLY\n")
        fh.write("*Instance, name=WP-1, part=WORKPIECE\n*End Instance\n")
        for p in model.placements:
            fh.write(f"*Instance, name=G-{p.placement_id}, "
                     f"part=GRAIN-{p.shape_index + 1}\n")
            c = p.translation_mm
            fh.write(f"{c[0]: .8e}, {c[1]: .8e}, {c[2]: .8e}\n")
            if abs(p.rotation_angle_deg) > 1e-12:
                a = p.rotation_axis
                q = c + a
                fh.write(f"{c[0]: .8e}, {c[1]: .8e}, {c[2]: .8e}, "
                         f"{q[0]: .8e}, {q[1]: .8e}, {q[2]: .8e}, "
                         f"{p.rotation_angle_deg: .8e}\n")
            fh.write("*End Instance\n")

        names = [f"G-{p.placement_id}.GRAIN_REF" for p in model.placements]
        fh.write("*Nset, nset=ALL_GRAIN_REF\n")
        for i in range(0, len(names), 8):
            fh.write(", ".join(names[i : i + 8]) + "\n")
        for nm in wp_sets:
            fh.write(f"*Nset, nset=A_{nm}, instance=WP-1\n{nm},\n")
        fh.write("*End Assembly\n")

        # ---- section controls (element deletion) ----
        fh.write("** Element deletion lets cut material leave the mesh, which is what\n")
        fh.write("** forms a chip. Without it the workpiece only deforms.\n")
        fh.write("*Section Controls, name=EC-1, ELEMENT DELETION=YES\n1., 1., 1.\n")

        # ---- materials ----
        _write_stone_material(fh, st)
        gm = MATERIALS[grain_material]
        fh.write(f"** {gm.note}\n*Material, name={gm.name}\n")
        fh.write(f"*Density\n{gm.density_kg_m3 * KGM3_TO_TONNE_MM3: .8e},\n")
        fh.write(f"*Elastic\n{gm.youngs_modulus_mpa: .8e}, {gm.poisson_ratio: .4f}\n")

        # ---- contact ----
        fh.write("*Surface Interaction, name=GRIT_STONE\n")
        fh.write("*Friction\n0.15,\n")
        fh.write("*Surface Behavior, pressure-overclosure=HARD\n")
        fh.write("*Contact, op=NEW\n")
        fh.write("*Contact Inclusions, ALL EXTERIOR\n")
        fh.write("*Contact Property Assignment\n ,  , GRIT_STONE\n")

        # ---- initial conditions and BCs ----
        fh.write("** The workpiece is held; the grits sweep across it.\n")
        fh.write("*Boundary\n")
        fh.write("A_WP_BOTTOM, ENCASTRE\n")
        fh.write("A_WP_ENTRY, 1, 3\n")
        fh.write("A_WP_EXIT, 1, 3\n")
        fh.write("A_WP_SIDE_A, 3, 3\n")
        fh.write("A_WP_SIDE_B, 3, 3\n")

        # ---- step ----
        fh.write("*Step, name=Grind, nlgeom=YES\n")
        fh.write(f"*Dynamic, Explicit\n, {duration: .8e}\n")
        if pr.mass_scaling_target_dt:
            fh.write(f"*Fixed Mass Scaling, dt={pr.mass_scaling_target_dt:.3e}, "
                     f"type=below min\n")
        fh.write("** Grits driven at the wheel surface speed. Translation rather than\n")
        fh.write(f"** rotation: over {duration * 1e6:.2f} us the wheel turns only "
                 f"{math.degrees(pr.wheel_speed_mm_s * duration / R):.4f} deg.\n")
        fh.write("*Boundary, type=VELOCITY\n")
        fh.write(f"ALL_GRAIN_REF, 1, 1, {v_vec[0]: .8e}\n")
        fh.write(f"ALL_GRAIN_REF, 2, 2, {v_vec[1]: .8e}\n")
        fh.write("ALL_GRAIN_REF, 3, 3, 0.\n")
        fh.write("ALL_GRAIN_REF, 4, 6, 0.\n")
        fh.write("*Output, field, number interval=20\n")
        fh.write("*Node Output\nU, V, RF\n")
        fh.write("*Element Output, directions=YES\nS, PEEQ, STATUS\n")
        fh.write("*Output, history, time interval="
                 f"{duration / 200: .6e}\n")
        fh.write("*Node Output, nset=ALL_GRAIN_REF\nRF1, RF2, RF3\n")
        fh.write("*End Step\n")

    return {
        "path": path,
        "size_bytes": os.path.getsize(path),
        "n_grits": len(model.placements),
        "n_grit_parts": len(used_shapes),
        "n_grit_elements": sum(len(model.shapes[p.shape_index].faces)
                               for p in model.placements),
        "n_workpiece_elements": int(len(wp_hexes)),
        "n_workpiece_nodes": int(len(wp_nodes)),
        "duration_s": duration,
        "dt_stable_s": dt_stable,
        "n_increments": n_inc,
        "engagement_um": pr.engagement_um,
        "max_protrusion_um": max_protrusion_mm * 1000.0,
        "wheel_speed_m_s": pr.wheel_speed_m_s,
        "workpiece_top_radius_mm": r_top,
        "theta_workpiece_deg": math.degrees(theta_w),
    }


def _header(fh: TextIO, model, wp, pr, st, max_prot, r_top, duration,
            dt_stable, n_inc, theta_w) -> None:
    fh.write("*Heading\n")
    fh.write("** Multi-grit scratch simulation: SEM-measured abrasive grits on stone\n")
    fh.write("** Units: mm, tonne, s, MPa, N. Density in tonne/mm^3.\n")
    fh.write("** Wheel axis = Z; grits translate tangentially.\n")
    fh.write(f"** wheel diameter (mm)      : {model.spec.diameter_mm:g}\n")
    fh.write(f"** grits                    : {len(model.placements)}\n")
    fh.write(f"** max grit protrusion (um) : {max_prot * 1000:.3f}\n")
    fh.write(f"** engagement depth (um)     : {pr.engagement_um:g}\n")
    fh.write(f"** workpiece top radius (mm): {r_top:.6f}\n")
    fh.write(f"** workpiece (mm)           : {wp.length_mm} x {wp.width_mm} x {wp.depth_mm}"
             f", {wp.element_size_mm * 1000:g} um elements\n")
    fh.write(f"** wheel speed (m/s)        : {pr.wheel_speed_m_s:g}\n")
    fh.write(f"** step time (s)            : {duration:.6e}\n")
    fh.write(f"** stable increment (s)     : {dt_stable:.6e}  -> ~{n_inc:,.0f} increments\n")
    fh.write("**\n")


def _write_stone_material(fh: TextIO, st: StoneMaterial) -> None:
    fh.write("** ==================================================================\n")
    fh.write("** WORKPIECE: stone\n")
    fh.write("**\n")
    fh.write("** Johnson-Holmquist II (JH-2) is the intended law for brittle rock, but\n")
    fh.write("** Abaqus/Explicit has no native JH-2 keyword -- it is supplied through a\n")
    fh.write("** VUMAT. The constants are listed below in the usual VUMAT order; add\n")
    fh.write("** *User Material and your VUMAT, then comment out the Drucker-Prager\n")
    fh.write("** block that follows.\n")
    fh.write("**\n")
    fh.write("** These are literature-typical GRANITE values. Verify them against your\n")
    fh.write("** own source before quoting any result.\n")
    fh.write("**\n")
    for k, v in st.jh2.items():
        fh.write(f"**   {k:28s} = {v:g}\n")
    fh.write("**\n")
    fh.write("** *Material, name=%s\n" % st.name)
    fh.write("** *Density\n**  %.8e,\n" % st.density_tonne_mm3)
    fh.write("** *User Material, constants=%d\n" % len(st.jh2))
    vals = list(st.jh2.values())
    for i in range(0, len(vals), 8):
        fh.write("**  " + ", ".join(f"{v:g}" for v in vals[i : i + 8]) + "\n")
    fh.write("** *Depvar, delete=2\n**  2,\n")
    fh.write("** ==================================================================\n")
    fh.write("** ACTIVE model: Mises plasticity with shear failure, so the deck runs and\n")
    fh.write("** the mechanics can be checked before the VUMAT is wired in.\n")
    fh.write("**\n")
    fh.write("** Not Drucker-Prager: *Shear Failure is only accepted as a suboption of\n")
    fh.write("** *Plastic. Pairing it with *Drucker Prager gives\n")
    fh.write("**   'in keyword *SHEARFAILURE ... The keyword is misplaced'\n")
    fh.write("** and the whole import aborts. Drucker-Prager is the better law for rock,\n")
    fh.write("** but then element deletion needs a different route (*Tensile Failure, or\n")
    fh.write("** the JH-2 VUMAT above, which handles damage itself).\n")
    fh.write(f"*Material, name={st.name}\n")
    fh.write(f"*Density\n{st.density_tonne_mm3: .8e},\n")
    fh.write(f"*Elastic\n{st.youngs_modulus_mpa: .8e}, {st.poisson_ratio: .4f}\n")
    fh.write("*Plastic\n")
    fh.write(f"{st.compressive_strength_mpa: .8e}, 0.\n")
    fh.write(f"{st.compressive_strength_mpa * 1.2: .8e}, {st.failure_strain: .4f}\n")
    fh.write("*Shear Failure, element deletion=YES\n")
    fh.write(f"{st.failure_strain: .4f},\n")
