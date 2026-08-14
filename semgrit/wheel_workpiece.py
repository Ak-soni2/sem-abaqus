"""Geometry-only Abaqus deck: grinding wheel plus a workpiece, assembled in contact.

Writes parts, the assembly, node sets and element-based surfaces, and nothing else.
No step, no interaction, no boundary condition, no load -- the analyst adds those.

Element types are Abaqus/Explicit ones throughout: **C3D8R** for the bond and the
workpiece, **R3D3** for the rigid grit facets. Abaqus/Explicit has no
full-integration C3D8, so writing C3D8R means no library conversion is needed after
import.

Positioning
-----------
The workpiece is placed so the **tallest grit tip is exactly tangent** to the surface
that will be ground -- zero initial penetration, one grit touching at one point, every
other grit clear.

Positioning by the *bond* outer radius instead would start the analysis with the grits
buried up to their full protrusion inside the workpiece. Abaqus contact would open with
that much overclosure, which produces a violent first increment and usually an
immediate abort. The bond therefore sits one maximum-protrusion clear of the
workpiece, with the grits bridging the gap, exactly as a real wheel does.

Coordinate system: wheel axis on **Z**, sector spanning theta = 0..sector_deg from
**+X**. The workpiece lies just outside the wheel radius, so the face being ground
points back toward the wheel axis.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from typing import Optional, TextIO

import numpy as np

from .abaqus import KGM3_TO_TONNE_MM3, MATERIALS, _write_int_set
from .wheel import UM_PER_MM, WheelModel, _rotation_matrix


@dataclass
class WorkpieceBlock:
    """Rectangular workpiece, meshed with C3D8R."""

    length_mm: float = 0.10      # along the wheel's tangential direction
    width_mm: float = 0.06       # along the wheel axis
    depth_mm: float = 0.025      # away from the ground surface
    element_size_mm: float = 0.0015

    name: str = "WORKPIECE"
    material: str = "STONE"
    density_kg_m3: float = 2650.0
    youngs_modulus_mpa: float = 50_000.0
    poisson_ratio: float = 0.25

    # Per-direction overrides; ``None`` falls back to the isotropic size above.
    # Worth having because the three directions do not cost the same thing. The
    # stable time increment follows the *smallest* element dimension, so coarsening
    # only the axial direction cuts the element count without lengthening the run,
    # whereas coarsening the cutting direction or the depth blurs the chip.
    element_size_length_mm: Optional[float] = None
    element_size_width_mm: Optional[float] = None
    element_size_depth_mm: Optional[float] = None

    # Graded mesh through the depth. ``surface_layer_mm = 0`` keeps the uniform mesh.
    #
    # Worth having because the depth is the direction the chip is removed *into*, so it
    # is what resolves chip thickness -- and it is also the direction with the most
    # wasted material, since the damage lives in the top few microns of a block that
    # may be hundreds deep. Grading gives a fine surface and a thick body at once, and
    # it is free in time: the stable increment follows the *smallest* element, so as
    # long as the surface layer is no finer than the cutting-direction size, dt does
    # not move at all.
    surface_layer_mm: float = 0.0
    """Depth of the finely meshed zone at the ground face. 0 = uniform mesh."""
    depth_growth: float = 1.3
    """Geometric ratio applied to successive layers below the fine zone."""
    max_depth_element_mm: float = 0.0
    """Cap on layer thickness, to stop the deep elements becoming slivers. 0 = no cap."""

    def requested_sizes(self) -> tuple[float, float, float]:
        h = self.element_size_mm
        return (self.element_size_length_mm or h,
                self.element_size_width_mm or h,
                self.element_size_depth_mm or h)

    def depth_coordinates(self) -> np.ndarray:
        """Node positions through the depth, starting at the ground face."""
        d = self.depth_mm
        h0 = self.requested_sizes()[2]
        if self.surface_layer_mm <= 0:
            nd = max(int(round(d / h0)), 1)
            return np.linspace(0.0, d, nd + 1)

        cap = self.max_depth_element_mm or float("inf")
        g = max(self.depth_growth, 1.0)
        xs = [0.0]
        for _ in range(max(int(round(self.surface_layer_mm / h0)), 1)):
            if xs[-1] + h0 >= d - 1e-12:
                break
            xs.append(xs[-1] + h0)
        h = h0
        while xs[-1] < d - 1e-12:
            h = min(h * g, cap)
            if xs[-1] + h >= d - 1e-12:
                break
            xs.append(xs[-1] + h)
        # Absorb the remainder into the last layer rather than adding a thin one: a
        # sliver at the back would become the smallest element in the model and drag
        # the stable increment down for no benefit at all.
        if d - xs[-1] < h0 and len(xs) > 1:
            xs[-1] = d
        else:
            xs.append(d)
        return np.asarray(xs, dtype=np.float64)

    def depth_layer_range(self) -> tuple[float, float]:
        w = self.depth_coordinates()
        dz = np.diff(w)
        return float(dz.min()), float(dz.max())

    def divisions(self) -> tuple[int, int, int]:
        hl, hw, hd = self.requested_sizes()
        return (max(int(round(self.length_mm / hl)), 1),
                max(int(round(self.width_mm / hw)), 1),
                len(self.depth_coordinates()) - 1)

    def element_sizes(self) -> tuple[float, float, float]:
        """Element size actually achieved, after rounding to a whole element count.

        The block keeps the dimensions the user asked for, so a requested size that
        does not divide them exactly is rounded. Reporting the request instead of
        this would understate the stable increment. The depth entry is the *finest*
        layer, since that is the one the stable increment sees.
        """
        nl, nw, _ = self.divisions()
        return (self.length_mm / nl, self.width_mm / nw, self.depth_layer_range()[0])

    def min_element_size(self) -> float:
        """The dimension that sets the stable time increment.

        For a rectangular brick Abaqus/Explicit takes the characteristic length as
        volume / largest face area, which for edges a <= b <= c is exactly ``a``.
        """
        return min(self.element_sizes())

    def n_elements(self) -> int:
        a, b, c = self.divisions()
        return a * b * c


def rotate_placements_about_z(placements: list, offset_deg: float) -> list:
    """Rigidly rotate placed grits about the wheel axis by ``offset_deg``.

    Lets the bond span a wide, visibly curved arc while the grits sit in a small
    window where they can actually engage the workpiece. Spreading grits over the
    whole arc instead would mean thousands of them, almost all too far from the
    workpiece to ever touch it -- pure contact-search cost.

    Each grit's placement is ``x -> centre + R_axis(v)``. Rotating the whole grit by
    ``Rz`` gives ``Rz*centre + (Rz*R_axis)(v)``, so the new centre is ``Rz*centre``
    and the new orientation is the *product* ``Rz*R_axis``, re-decomposed to
    axis-angle. Rotating the axis alone would be wrong.
    """
    import copy

    from .wheel import matrix_to_axis_angle

    a = math.radians(offset_deg)
    rz = np.array([[math.cos(a), -math.sin(a), 0.0],
                   [math.sin(a), math.cos(a), 0.0],
                   [0.0, 0.0, 1.0]])
    out = []
    for p in placements:
        q = copy.copy(p)
        q.translation_mm = rz @ p.translation_mm
        r_old = _rotation_matrix(p.rotation_axis, math.radians(p.rotation_angle_deg))
        axis, angle = matrix_to_axis_angle(rz @ r_old)
        q.rotation_axis = axis
        q.rotation_angle_deg = angle
        q.theta_deg = p.theta_deg + offset_deg
        out.append(q)
    return out


def grit_tip_radius(model: WheelModel) -> tuple[float, float]:
    """Largest and smallest radial reach over every placed grit, in mm."""
    lo = math.inf
    hi = -math.inf
    for p in model.placements:
        s = model.shapes[p.shape_index]
        rot = _rotation_matrix(p.rotation_axis, math.radians(p.rotation_angle_deg))
        v = (s.vertices - s.centroid_um) / UM_PER_MM @ rot.T + p.translation_mm
        r = np.hypot(v[:, 0], v[:, 1])
        hi = max(hi, float(r.max()))
        lo = min(lo, float(r.min()))
    return hi, lo


def build_block_mesh(
    wp: WorkpieceBlock,
    ground_face_centre: np.ndarray,
    e_len: np.ndarray,
    e_wid: np.ndarray,
    e_depth: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """Structured C3D8R block. ``e_depth`` points away from the ground face."""
    nl, nw, nd = wp.divisions()
    u = np.linspace(-wp.length_mm / 2.0, wp.length_mm / 2.0, nl + 1)
    v = np.linspace(-wp.width_mm / 2.0, wp.width_mm / 2.0, nw + 1)
    # Not a linspace: the depth may be graded, fine at the ground face and coarsening
    # into the body. ``depth_coordinates`` returns the uniform case unchanged.
    w = wp.depth_coordinates()
    if len(w) - 1 != nd:
        raise ValueError("depth coordinates disagree with divisions()")

    def nid(i: int, j: int, k: int) -> int:
        return (i * (nw + 1) + j) * (nd + 1) + k

    nodes = np.empty(((nl + 1) * (nw + 1) * (nd + 1), 3), dtype=np.float64)
    for i, ui in enumerate(u):
        for j, vj in enumerate(v):
            for k, wk in enumerate(w):
                nodes[nid(i, j, k)] = (
                    ground_face_centre + e_len * ui + e_wid * vj + e_depth * wk
                )

    hexes: list[tuple[int, ...]] = []
    for i in range(nl):
        for j in range(nw):
            for k in range(nd):
                hexes.append((
                    nid(i, j, k), nid(i + 1, j, k),
                    nid(i + 1, j + 1, k), nid(i, j + 1, k),
                    nid(i, j, k + 1), nid(i + 1, j, k + 1),
                    nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1),
                ))

    sets = {
        "WP_GROUND_FACE": np.array(
            [nid(i, j, 0) for i in range(nl + 1) for j in range(nw + 1)]),
        "WP_BACK_FACE": np.array(
            [nid(i, j, nd) for i in range(nl + 1) for j in range(nw + 1)]),
        "WP_END_A": np.array(
            [nid(0, j, k) for j in range(nw + 1) for k in range(nd + 1)]),
        "WP_END_B": np.array(
            [nid(nl, j, k) for j in range(nw + 1) for k in range(nd + 1)]),
        "WP_SIDE_A": np.array(
            [nid(i, 0, k) for i in range(nl + 1) for k in range(nd + 1)]),
        "WP_SIDE_B": np.array(
            [nid(i, nw, k) for i in range(nl + 1) for k in range(nd + 1)]),
        "WP_ALL_NODES": np.arange(len(nodes)),
    }
    return nodes, np.asarray(hexes, dtype=np.int64), sets


# Boundary faces of a C3D8R with nodes 0-3 on one face and 4-7 above them.
_HEX_FACE_ID = {
    "S1": (0, 1, 2, 3), "S2": (4, 5, 6, 7), "S3": (0, 1, 5, 4),
    "S4": (1, 2, 6, 5), "S5": (2, 3, 7, 6), "S6": (3, 0, 4, 7),
}


def write_wheel_workpiece_inp(
    path: str,
    model: WheelModel,
    workpiece: Optional[WorkpieceBlock] = None,
    clearance_um: float = 0.0,
    grain_material: str = "diamond",
    bond_material: str = "vitrified_bond",
    model_name: str = "GRINDING_WHEEL_AND_WORKPIECE",
) -> dict:
    """Write the geometry-only assembled deck. Returns a summary dict."""
    wp = workpiece or WorkpieceBlock()
    if not model.placements:
        raise ValueError("the wheel model has no grits placed")

    R = model.spec.outer_radius_mm
    tip_hi, _ = grit_tip_radius(model)
    max_protrusion_mm = tip_hi - R

    # Tangency: ground face exactly at the tallest tip, plus any requested gap.
    r_ground = tip_hi + clearance_um / 1000.0

    # Put the workpiece at the angular centre of the grit patch.
    thetas = [p.theta_deg for p in model.placements]
    theta_c = math.radians(0.5 * (min(thetas) + max(thetas)))
    e_r = np.array([math.cos(theta_c), math.sin(theta_c), 0.0])
    e_t = np.array([-math.sin(theta_c), math.cos(theta_c), 0.0])
    e_z = np.array([0.0, 0.0, 1.0])

    wp_nodes, wp_hexes, wp_sets = build_block_mesh(
        wp, e_r * r_ground, e_t, e_z, e_r)

    used_shapes = sorted({p.shape_index for p in model.placements})
    shape_verts = [(s.vertices - s.centroid_um) / UM_PER_MM for s in model.shapes]
    nl, nw, nd = wp.divisions()

    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("*Heading\n")
        fh.write("** %s\n" % model_name)
        fh.write("** Geometry only: parts, assembly, sets and surfaces. No step, no\n")
        fh.write("** interaction, no boundary condition, no load.\n")
        fh.write("** Units: mm, tonne, s, MPa, N.  Wheel axis = Z.\n")
        fh.write("** Element types are Explicit: C3D8R solids, R3D3 rigid facets.\n")
        fh.write("**\n")
        fh.write("** wheel diameter (mm)        : %g\n" % model.spec.diameter_mm)
        fh.write("** wheel width (mm)           : %g\n" % model.spec.width_mm)
        fh.write("** sector (deg)               : %g\n" % model.spec.sector_deg)
        fh.write("** bond outer radius (mm)     : %.6f\n" % R)
        fh.write("** grits                      : %d\n" % len(model.placements))
        fh.write("** max grit protrusion (um)   : %.4f\n" % (max_protrusion_mm * 1000))
        fh.write("** tallest grit tip radius(mm): %.6f\n" % tip_hi)
        fh.write("** workpiece ground face at   : r = %.6f mm  (tangent to that tip)\n"
                 % r_ground)
        fh.write("** initial clearance (um)     : %g\n" % clearance_um)
        fh.write("** workpiece (mm)             : %g x %g x %g, %g um elements\n" % (
            wp.length_mm, wp.width_mm, wp.depth_mm, wp.element_size_mm * 1000))
        fh.write("** workpiece elements         : %d C3D8R\n" % len(wp_hexes))
        fh.write("** workpiece placed at theta  : %.6f deg\n" % math.degrees(theta_c))
        fh.write("**\n")
        fh.write("** Move the wheel radially inward by your depth of cut to engage.\n")
        fh.write("**\n")

        # ---------------- grit parts ----------------
        for idx in used_shapes:
            s = model.shapes[idx]
            fh.write("*Part, name=GRAIN-%d\n*Node\n" % (idx + 1))
            for i, v in enumerate(shape_verts[idx], start=1):
                fh.write("%d, %.9e, %.9e, %.9e\n" % (i, v[0], v[1], v[2]))
            fh.write("*Element, type=R3D3\n")
            for e, tri in enumerate(s.faces, start=1):
                fh.write("%d, %d, %d, %d\n" % (e, tri[0] + 1, tri[1] + 1, tri[2] + 1))
            fh.write("*Elset, elset=GRAIN_ALL, generate\n1, %d, 1\n" % len(s.faces))
            ref = len(shape_verts[idx]) + 1
            fh.write("*Node\n%d, 0., 0., 0.\n" % ref)
            fh.write("*Nset, nset=GRAIN_REF\n%d,\n" % ref)
            fh.write("*Rigid Body, ref node=GRAIN_REF, elset=GRAIN_ALL\n")
            fh.write("*Surface, type=ELEMENT, name=GRAIN_SURF\nGRAIN_ALL, SPOS\n")
            fh.write("*End Part\n")

        # ---------------- bond ----------------
        fh.write("*Part, name=WHEEL_BOND\n*Node\n")
        for i, v in enumerate(model.body_nodes, start=1):
            fh.write("%d, %.9e, %.9e, %.9e\n" % (i, v[0], v[1], v[2]))
        fh.write("*Element, type=C3D8R\n")
        for e, h in enumerate(model.body_hexes, start=1):
            fh.write("%d, " % e + ", ".join(str(int(n) + 1) for n in h) + "\n")
        fh.write("*Elset, elset=BOND_ALL, generate\n1, %d, 1\n" % len(model.body_hexes))
        fh.write("*Solid Section, elset=BOND_ALL, material=%s\n,\n"
                 % MATERIALS[bond_material].name)
        for name, ids in model.node_sets.items():
            if name.startswith("_"):
                continue
            fh.write("*Nset, nset=%s\n" % name)
            _write_int_set(fh, [int(i) + 1 for i in ids])
        for key, (surf, face) in {
            "_EL_OUTER": ("BOND_OUTER_SURF", "S4"),
            "_EL_BORE": ("BOND_BORE_SURF", "S6"),
            "_EL_ZMIN": ("BOND_ZMIN_SURF", "S1"),
            "_EL_ZMAX": ("BOND_ZMAX_SURF", "S2"),
            "_EL_SECTOR_START": ("BOND_SECTOR_START_SURF", "S3"),
            "_EL_SECTOR_END": ("BOND_SECTOR_END_SURF", "S5"),
        }.items():
            ids = model.node_sets.get(key)
            if ids is None or len(ids) == 0:
                continue
            # Element sets are named without a leading underscore: CAE uses "_"
            # prefixes for its own internal sets, and these have to be referenced
            # by name from the assembly.
            fh.write("*Elset, elset=ES_%s\n" % surf)
            _write_int_set(fh, [int(i) + 1 for i in ids])
            fh.write("*Surface, type=ELEMENT, name=%s\nES_%s, %s\n" % (surf, surf, face))
        fh.write("*End Part\n")

        # ---------------- workpiece ----------------
        fh.write("*Part, name=%s\n*Node\n" % wp.name)
        for i, v in enumerate(wp_nodes, start=1):
            fh.write("%d, %.9e, %.9e, %.9e\n" % (i, v[0], v[1], v[2]))
        fh.write("*Element, type=C3D8R\n")
        for e, h in enumerate(wp_hexes, start=1):
            fh.write("%d, " % e + ", ".join(str(int(n) + 1) for n in h) + "\n")
        fh.write("*Elset, elset=WP_ALL, generate\n1, %d, 1\n" % len(wp_hexes))
        fh.write("*Solid Section, elset=WP_ALL, material=%s\n,\n" % wp.material)
        for name, ids in wp_sets.items():
            fh.write("*Nset, nset=%s\n" % name)
            _write_int_set(fh, [int(i) + 1 for i in ids])
        # The ground face is S1 of the k=0 layer of elements.
        ground_els = []
        for i in range(nl):
            for j in range(nw):
                ground_els.append(((i * nw) + j) * nd + 1)   # 1-based element id
        fh.write("*Elset, elset=ES_WP_GROUND\n")
        _write_int_set(fh, ground_els)
        fh.write("*Surface, type=ELEMENT, name=WP_GROUND_SURF\nES_WP_GROUND, S1\n")
        fh.write("*End Part\n")

        # ---------------- assembly ----------------
        fh.write("*Assembly, name=ASSEMBLY\n")
        fh.write("*Instance, name=BOND-1, part=WHEEL_BOND\n*End Instance\n")
        for p in model.placements:
            fh.write("*Instance, name=G-%d, part=GRAIN-%d\n"
                     % (p.placement_id, p.shape_index + 1))
            c = p.translation_mm
            fh.write("%.9e, %.9e, %.9e\n" % (c[0], c[1], c[2]))
            if abs(p.rotation_angle_deg) > 1e-12:
                a = p.rotation_axis
                q = c + a
                fh.write("%.9e, %.9e, %.9e, %.9e, %.9e, %.9e, %.9e\n"
                         % (c[0], c[1], c[2], q[0], q[1], q[2], p.rotation_angle_deg))
            fh.write("*End Instance\n")
        # The workpiece nodes are already at their final coordinates.
        fh.write("*Instance, name=WP-1, part=%s\n*End Instance\n" % wp.name)

        names = ["G-%d.GRAIN_REF" % p.placement_id for p in model.placements]
        fh.write("*Nset, nset=ALL_GRIT_REF\n")
        for i in range(0, len(names), 8):
            fh.write(", ".join(names[i : i + 8]) + "\n")
        for name in model.node_sets:
            if name.startswith("_"):
                continue
            fh.write("*Nset, nset=A_%s, instance=BOND-1\n%s,\n" % (name, name))
        for name in wp_sets:
            fh.write("*Nset, nset=A_%s, instance=WP-1\n%s,\n" % (name, name))
        # An assembly-level element surface must reference an instance-qualified
        # ELEMENT SET plus a face identifier. It cannot reference a part-level
        # *Surface: doing so gives
        #   "The following sets were not found when generating the surface ..."
        # and aborts the whole import.
        fh.write("*Surface, type=ELEMENT, name=A_WP_GROUND_SURF\n")
        fh.write("WP-1.ES_WP_GROUND, S1\n")
        if model.node_sets.get("_EL_OUTER") is not None:
            fh.write("*Surface, type=ELEMENT, name=A_BOND_OUTER_SURF\n")
            fh.write("BOND-1.ES_BOND_OUTER_SURF, S4\n")
        fh.write("*End Assembly\n")

        # ---------------- materials (placeholders) ----------------
        fh.write("** Placeholder materials so the *Solid Section blocks resolve and the\n")
        fh.write("** deck imports cleanly. Replace with your own, including the JH-2\n")
        fh.write("** VUMAT for the workpiece.\n")
        gm = MATERIALS[grain_material]
        bm = MATERIALS[bond_material]
        for m in (gm, bm):
            fh.write("*Material, name=%s\n" % m.name)
            fh.write("*Density\n%.8e,\n" % (m.density_kg_m3 * KGM3_TO_TONNE_MM3))
            fh.write("*Elastic\n%.8e, %.4f\n" % (m.youngs_modulus_mpa, m.poisson_ratio))
        fh.write("*Material, name=%s\n" % wp.material)
        fh.write("*Density\n%.8e,\n" % (wp.density_kg_m3 * KGM3_TO_TONNE_MM3))
        fh.write("*Elastic\n%.8e, %.4f\n" % (wp.youngs_modulus_mpa, wp.poisson_ratio))

    return {
        "path": path,
        "size_bytes": os.path.getsize(path),
        "n_grits": len(model.placements),
        "n_grit_parts": len(used_shapes),
        "n_grit_facets": sum(len(model.shapes[p.shape_index].faces)
                             for p in model.placements),
        "n_bond_elements": int(len(model.body_hexes)),
        "n_workpiece_elements": int(len(wp_hexes)),
        "n_workpiece_nodes": int(len(wp_nodes)),
        "bond_outer_radius_mm": R,
        "max_grit_protrusion_um": max_protrusion_mm * 1000.0,
        "tallest_grit_tip_radius_mm": tip_hi,
        "workpiece_ground_radius_mm": r_ground,
        "clearance_um": clearance_um,
        "theta_workpiece_deg": math.degrees(theta_c),
        "bond_to_workpiece_gap_um": (r_ground - R) * 1000.0,
    }
