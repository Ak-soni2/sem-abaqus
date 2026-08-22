"""All-rigid grinding wheel plus one deformable workpiece. Geometry only.

Difference from :mod:`semgrit.wheel_workpiece`
----------------------------------------------
There, the bond was a deformable C3D8R solid and every grit was its own rigid body
with its own reference node. That gave one rigid body per grit -- 75 to 800 of them --
so driving the wheel meant constraining hundreds of independent reference nodes, and
the bond's small elements fed into the Explicit stable time increment.

Here the **whole wheel is a single discrete rigid body**: the bond rim is a closed
rigid *surface* shell (R3D4) and every grit's facets (R3D3) are merged into the same
part, tied to **one** reference node on the wheel axis at the origin. Consequences:

* one ``*Rigid Body``, one reference node, so wheel rotation is one BC;
* the bond contributes nothing to the stable time increment, and nothing to the
  deformable element count -- only the workpiece does;
* no ``*Solid Section`` and no material for anything but the workpiece.

Merging rather than instancing is deliberate. A rigid body is defined by one element
set, and an element set cannot span parts, so the grits have to live in the same part
as the bond. The instance transform is therefore baked into the vertices:
``v_world = v_local @ R.T + t``, which is exactly what ``*Instance`` computes from a
translation followed by a rotation about an axis through the translated centre.

Positioning
-----------
The ground face is placed so that the **tallest grit material anywhere inside the
workpiece footprint is exactly tangent to it** -- zero penetration, contact at one
point. The footprint test clips each grit facet to the block's tangential/axial extent
before taking its radial reach, because a facet can cross the footprint edge between
two vertices; testing vertices alone leaves up to one facet-edge of penetration.

Taking the globally tallest grit instead would be wrong here: grits cover the whole
arc, so the tallest one is typically a millimetre away from the workpiece, and
referencing it would leave every grit that can actually engage several microns clear
of the surface -- nothing would cut.

Coordinate system: wheel axis on **Z**, sector spanning theta = 0..sector_deg from
**+X**, workpiece just outside the rim so its ground face looks back at the axis.
Units mm, tonne, s, MPa, N.
"""

from __future__ import annotations

import math
import os
from typing import Optional, Sequence

import numpy as np

from .abaqus import KGM3_TO_TONNE_MM3, MATERIALS, _write_int_set
from .wheel import UM_PER_MM, WheelModel, WheelSpec, _rotation_matrix
from .wheel_workpiece import WorkpieceBlock, build_block_mesh

# Face groups of the rim shell, and the outward direction each one must point.
_SHELL_GROUPS = ("OUTER", "BORE", "ZMIN", "ZMAX", "SECTOR_START", "SECTOR_END")

# Coordinates are written with this many digits, giving a quantum of about
# 2.5e-11 mm at a radius of 25 mm. The earlier %.9e left a quantum of 1e-8 mm,
# which was coarse enough to push grit facets 8e-9 mm through a ground face that
# had been placed exactly tangent in full precision.
_COORD_FMT = "%.12e"

# Tangency is computed from the *quantised* vertices, so the file is self-consistent.
# This guard then covers the workpiece nodes' own quantisation, keeping the initial
# overclosure at exactly zero rather than a few tens of femtometres either way.
_WRITE_GUARD_MM = 1e-9


def _node_line(i: int, v) -> str:
    return "%d, %s, %s, %s\n" % (i, _COORD_FMT % v[0], _COORD_FMT % v[1],
                                 _COORD_FMT % v[2])


def quantise(v: np.ndarray) -> np.ndarray:
    """Round coordinates to the precision ``_COORD_FMT`` will actually emit.

    Geometry has to be checked against what lands in the file, not against what was
    computed. Placing a face tangent to an unrounded vertex and then writing both to
    finite precision leaves the sign of the gap up to the last decimal digit.
    """
    flat = np.asarray(v, dtype=np.float64).ravel()
    out = np.fromiter((float(_COORD_FMT % x) for x in flat), dtype=np.float64,
                      count=flat.size)
    return out.reshape(np.shape(v))


def build_rim_shell(
    spec: WheelSpec,
) -> tuple[np.ndarray, np.ndarray, dict[str, list[int]]]:
    """Closed rigid surface of the rim sector, meshed with 4-node quads.

    A discrete rigid body only needs a surface, so the rim is emitted as the
    bounding faces of the swept annular sector rather than as a solid. Closing it
    (bore, both axial faces, both cut faces) costs a handful of elements and makes
    the sector read as a curved *block* in the viewport instead of a ribbon, which is
    what makes the arc visible at all at short arc lengths.

    A full 360 deg wheel wraps instead: the last column of nodes *is* the first, so
    there is no seam and no pair of cut faces.

    Every quad is ordered so its normal points **out of** the body, so all contact
    surfaces are ``SPOS``.
    """
    r0, r1 = spec.inner_radius_mm, spec.outer_radius_mm
    if r1 - r0 <= 0:
        raise ValueError("rim has zero radial thickness")

    full = spec.is_full_circle
    n_r = max(int(spec.radial_divisions), 1)
    n_z = max(int(spec.axial_divisions), 1)
    n_t = spec.circumferential_divisions()
    n_t_nodes = n_t if full else n_t + 1

    radii = np.linspace(r0, r1, n_r + 1)
    zs = np.linspace(-spec.width_mm / 2.0, spec.width_mm / 2.0, n_z + 1)
    thetas = (np.linspace(0.0, 2 * math.pi, n_t, endpoint=False) if full
              else np.linspace(0.0, spec.sector_rad, n_t + 1))

    def nid(i_r: int, i_t: int, i_z: int) -> int:
        # The modulo wraps a full circle and is a no-op on a sector, where
        # n_t_nodes = n_t + 1 and i_t never exceeds n_t.
        return (i_r * n_t_nodes + (i_t % n_t_nodes)) * (n_z + 1) + i_z

    nodes = np.empty(((n_r + 1) * n_t_nodes * (n_z + 1), 3), dtype=np.float64)
    for i_r, r in enumerate(radii):
        for i_t, tt in enumerate(thetas):
            x, y = r * math.cos(tt), r * math.sin(tt)
            for i_z, z in enumerate(zs):
                nodes[nid(i_r, i_t, i_z)] = (x, y, z)

    quads: list[tuple[int, int, int, int]] = []
    groups: dict[str, list[int]] = {k: [] for k in _SHELL_GROUPS}

    def add(name: str, quad: tuple[int, int, int, int]) -> None:
        groups[name].append(len(quads))
        quads.append(quad)

    # Orderings chosen so the right-hand-rule normal is the outward one; see the
    # numeric check in verify_rigid_deck.py, which recomputes every normal.
    for i_t in range(n_t):
        for i_z in range(n_z):
            add("OUTER", (nid(n_r, i_t, i_z), nid(n_r, i_t + 1, i_z),
                          nid(n_r, i_t + 1, i_z + 1), nid(n_r, i_t, i_z + 1)))
            add("BORE", (nid(0, i_t, i_z), nid(0, i_t, i_z + 1),
                         nid(0, i_t + 1, i_z + 1), nid(0, i_t + 1, i_z)))
    for i_r in range(n_r):
        for i_t in range(n_t):
            add("ZMIN", (nid(i_r, i_t, 0), nid(i_r, i_t + 1, 0),
                         nid(i_r + 1, i_t + 1, 0), nid(i_r + 1, i_t, 0)))
            add("ZMAX", (nid(i_r, i_t, n_z), nid(i_r + 1, i_t, n_z),
                         nid(i_r + 1, i_t + 1, n_z), nid(i_r, i_t + 1, n_z)))
    if not full:
        for i_r in range(n_r):
            for i_z in range(n_z):
                add("SECTOR_START", (nid(i_r, 0, i_z), nid(i_r + 1, 0, i_z),
                                     nid(i_r + 1, 0, i_z + 1), nid(i_r, 0, i_z + 1)))
                add("SECTOR_END", (nid(i_r, n_t, i_z), nid(i_r, n_t, i_z + 1),
                                   nid(i_r + 1, n_t, i_z + 1), nid(i_r + 1, n_t, i_z)))

    quads = np.asarray(quads, dtype=np.int64)

    # Only the boundary of the lattice is meshed, so with more than one radial or axial
    # division the interior nodes belong to no element. Emitting them would leave
    # (n_r-1)(n_t-1)(n_z-1) orphans in the part -- 95,381 of them on a 12x300x30 shell,
    # bloating the deck for nothing. Keep the used ones and renumber.
    used = np.unique(quads)
    if len(used) < len(nodes):
        remap = np.full(len(nodes), -1, dtype=np.int64)
        remap[used] = np.arange(len(used))
        nodes = nodes[used]
        quads = remap[quads]

    return nodes, quads, groups


def rim_mass_properties(
    spec: WheelSpec, density_kg_m3: float
) -> tuple[float, np.ndarray]:
    """Mass and inertia tensor about the origin of the solid rim sector.

    The rigid body is a surface, so it has no mass of its own; without any, Abaqus
    aborts unless every reference-node degree of freedom carries a boundary
    condition. Supplying the inertia of the *solid* sector the surface encloses makes
    the part behave sensibly whatever subset of DOF the analyst prescribes, and --
    unlike ``density``/``thickness`` on ``*Rigid Body`` -- adds no contact offset, so
    the tangency built above survives untouched.

    Closed forms for the annular sector r0..r1, theta 0..dt, z -w/2..w/2:
    ``m = rho (r1^2-r0^2)/2 dt w`` and the second moments follow from
    ``int r^3 dr = (r1^4-r0^4)/4``.
    """
    rho = density_kg_m3 * KGM3_TO_TONNE_MM3
    r0, r1 = spec.inner_radius_mm, spec.outer_radius_mm
    dt = spec.sector_rad
    w = spec.width_mm

    a2 = (r1 ** 2 - r0 ** 2) / 2.0
    a4 = (r1 ** 4 - r0 ** 4) / 4.0
    m = rho * a2 * dt * w

    # Angular integrals over 0..dt.
    c2 = 0.5 * dt + 0.25 * math.sin(2 * dt)          # int cos^2
    s2 = 0.5 * dt - 0.25 * math.sin(2 * dt)          # int sin^2
    sc = 0.5 * math.sin(dt) ** 2                     # int sin cos

    ixx_m = rho * w * a4 * c2          # int x^2 dm
    iyy_m = rho * w * a4 * s2          # int y^2 dm
    ixy_m = rho * w * a4 * sc          # int x y dm
    izz_m = rho * a2 * dt * w ** 3 / 12.0   # int z^2 dm

    tensor = np.array([
        [iyy_m + izz_m, -ixy_m, 0.0],
        [-ixy_m, ixx_m + izz_m, 0.0],
        [0.0, 0.0, ixx_m + iyy_m],
    ])
    return float(m), tensor


def bake_grit(model: WheelModel, placement) -> np.ndarray:
    """World-coordinate vertices of one placed grit, in mm.

    ``*Instance`` applies the translation and then a rotation about an axis through
    the *translated* centre, which reduces to ``v_local @ R.T + t``. Baking with any
    other order silently moves every grit.
    """
    shape = model.shapes[placement.shape_index]
    local = (shape.vertices - shape.centroid_um) / UM_PER_MM
    rot = _rotation_matrix(placement.rotation_axis,
                           math.radians(placement.rotation_angle_deg))
    return local @ rot.T + placement.translation_mm


def _clip_to_slab(poly: np.ndarray, axis: int, lo: float, hi: float) -> np.ndarray:
    """Sutherland-Hodgman clip of a convex polygon to ``lo <= p[axis] <= hi``."""
    for sign, bound in ((1.0, lo), (-1.0, hi)):
        if len(poly) == 0:
            return poly
        out: list[np.ndarray] = []
        n = len(poly)
        for i in range(n):
            a, b = poly[i], poly[(i + 1) % n]
            da = sign * (a[axis] - bound)
            db = sign * (b[axis] - bound)
            if da >= 0:
                out.append(a)
            if (da >= 0) != (db >= 0):
                out.append(a + (b - a) * (da / (da - db)))
        poly = np.asarray(out) if out else np.zeros((0, 3))
    return poly


def ground_radius(
    grit_frames: Sequence[np.ndarray],
    faces: Sequence[np.ndarray],
    half_len: float,
    half_wid: float,
) -> tuple[float, int]:
    """Radius of the ground face that makes the grits exactly tangent to it.

    ``grit_frames[i]`` holds each grit's vertices already expressed in the tangent
    frame as columns ``(a, b, z)`` = (radial, tangential, axial). Returns the largest
    radial reach of any grit *material* lying inside the footprint
    ``|b| <= half_len, |z| <= half_wid``, together with the grit that governs it.

    Facets are clipped to the footprint rather than filtered by vertex, because ``a``
    is affine over a triangle: the maximum over the clipped polygon is the true
    maximum of the grit surface inside the footprint, while a vertex-only test misses
    facets that cross the boundary and leaves them penetrating the block.
    """
    best = -math.inf
    who = -1
    for gi, (v, f) in enumerate(zip(grit_frames, faces)):
        inside = (np.abs(v[:, 1]) <= half_len) & (np.abs(v[:, 2]) <= half_wid)
        if inside.all():
            hi = float(v[:, 0].max())
        else:
            # Only facets touching the footprint can matter.
            hi = -math.inf
            for tri in f:
                p = v[tri]
                if inside[tri].all():
                    hi = max(hi, float(p[:, 0].max()))
                    continue
                if (np.abs(p[:, 1]).min() > half_len
                        or np.abs(p[:, 2]).min() > half_wid):
                    continue
                q = _clip_to_slab(p, 1, -half_len, half_len)
                q = _clip_to_slab(q, 2, -half_wid, half_wid) if len(q) else q
                if len(q):
                    hi = max(hi, float(q[:, 0].max()))
        if hi > best:
            best, who = hi, gi
    if who < 0:
        raise ValueError("no grit material lies inside the workpiece footprint")
    return best, who



def engaging_grits(frames: Sequence[np.ndarray], wp, engage_window_mm: float = 0.0):
    """Indices of the grits that can reach the block, and the window that decided it.

    This *is* the deck's ``ES_GRITS_ENGAGE`` set -- the writer calls this, and so does
    anything that wants to show the contact surface. It used to live inline in the
    writer, which meant a viewer highlighting "the contact" was drawing a lookalike
    computed from a different window. On a project where two implementations of one
    idea have twice agreed on the same wrong answer, that is not a risk worth taking.
    """
    if wp is None:
        return list(range(len(frames))), float("inf")
    win = wp.length_mm / 2.0 + (engage_window_mm or 0.0)
    return ([i for i, v in enumerate(frames)
             if float(np.abs(v[:, 1]).min()) <= win], win)


WP_POSITIONS = ("centred", "first grit at entry", "under the tallest grit",
                "custom angle")
"""Where along the arc the block sits. The wheel turns so its surface travels toward
decreasing theta, so grits arrive from the high-theta end -- that end is the entry."""


def place_workpiece(model: WheelModel, wp, clearance_um: float = 0.0,
                    position: str = "centred", position_deg: float = 0.0,
                    entry_high_theta: bool = True) -> dict:
    """Bake the grits and decide where the workpiece block sits on them.

    One implementation, three callers: the deck writer, the planner behind the preview,
    and the auto depth-of-cut. Duplicating it would let the preview promise a clearance
    the deck does not deliver.

    ``position`` chooses the angular station:

    ``centred``
        Mid-arc, which leaves grit either side however the wheel turns. A full wheel
        has no distinguished angle, so theta = 0 is nominal.
    ``first grit at entry``
        The block's entry edge -- the high-theta end, the one grits reach first --
        sits at the leading grain, so the pass begins with the first abrasive right
        at the edge and every grain downstream then sweeps across.
    ``under the tallest grit``
        Centred on the most protruding grain the block can reach, which is the grain
        that will take the deepest cut.
    ``custom angle``
        ``position_deg``, measured from the global x axis.

    Whatever is asked for, the footprint has to contain grit material or there is no
    tangent radius to sit on. The block subtends a tiny fraction of a wheel (0.048 mm
    of a 157 mm circumference is 0.11 deg), so an empty footprint is common; when it
    happens the block falls back to the tallest grit *it can reach axially*, because
    picking the globally tallest fails whenever the dressed face is wider than the
    block.
    """
    spec = model.spec
    baked = [quantise(bake_grit(model, q)) for q in model.placements]
    faces = [model.shapes[q.shape_index].faces for q in model.placements]
    reach = [float(np.hypot(v[:, 0], v[:, 1]).max()) for v in baked]

    def frame(tc):
        er = np.array([math.cos(tc), math.sin(tc), 0.0])
        et = np.array([-math.sin(tc), math.cos(tc), 0.0])
        ez = np.array([0.0, 0.0, 1.0])
        return er, et, ez, np.column_stack([er, et, ez])

    if position not in WP_POSITIONS:
        raise ValueError("wp_position must be one of %s, not %r"
                         % (", ".join(WP_POSITIONS), position))

    theta_c = 0.0 if spec.is_full_circle else 0.5 * spec.sector_rad
    relocated = False
    requested_deg = None
    cand: list[int] = []
    if wp is not None:
        hb, hz = wp.length_mm / 2.0, wp.width_mm / 2.0
        # Only grits the block can reach across the face are candidates: a grain on a
        # part of the dressed face the block never runs over cannot govern anything.
        cand = [i for i, v in enumerate(baked) if float(np.abs(v[:, 2]).min()) <= hz]

        if position != "centred":
            if not cand:
                raise ValueError(
                    "no grit lies within the workpiece's %.4f mm width; widen the "
                    "block or narrow the dressed face (grit_width_window_mm)"
                    % wp.width_mm)
            if position == "under the tallest grit":
                c = baked[max(cand, key=lambda i: reach[i])].mean(axis=0)
                theta_c = math.atan2(c[1], c[0])
            elif position == "first grit at entry":
                # Whichever end the grains arrive from is the entry. With the default
                # sense (VR3 = -omega, surface toward decreasing theta) that is the
                # high-theta end; reversing the wheel reverses this too, or the block
                # would be placed at the exit and the pass would begin with no grit on
                # it at all.
                if entry_high_theta:
                    lead = max(float(np.arctan2(baked[i][:, 1], baked[i][:, 0]).max())
                               for i in cand)
                else:
                    lead = min(float(np.arctan2(baked[i][:, 1], baked[i][:, 0]).min())
                               for i in cand)
                # The block's entry edge subtends atan(hb / r_ground) at the axis, but
                # r_ground is not known until the block has been seated -- which needs
                # theta_c. Seat it once at the tip radius, then redo the alignment
                # against the radius that came out. The correction is nanometres, but
                # without it the entry edge the plan reports is not the edge that was
                # aligned, and the two disagree in the last digits.
                r_ref = max(reach)
                _way = 1.0 if entry_high_theta else -1.0
                for _ in range(2):
                    theta_c = lead - _way * math.atan2(hb, r_ref)
                    _, _, _, b2 = frame(theta_c)
                    try:
                        r_t2, _ = ground_radius([v @ b2 for v in baked], faces, hb, hz)
                    except ValueError:
                        break
                    r_ref = r_t2 + clearance_um / 1000.0 + _WRITE_GUARD_MM
            else:                                     # custom angle
                theta_c = math.radians(position_deg)
            requested_deg = math.degrees(theta_c)

        _, _, _, basis = frame(theta_c)
        try:
            ground_radius([v @ basis for v in baked], faces, hb, hz)
        except ValueError:
            if not cand:
                raise ValueError(
                    "no grit lies within the workpiece's %.4f mm width; widen the "
                    "block or narrow the dressed face (grit_width_window_mm)"
                    % wp.width_mm)
            c = baked[max(cand, key=lambda i: reach[i])].mean(axis=0)
            theta_c = math.atan2(c[1], c[0])
            relocated = True

    e_r, e_t, e_z, basis = frame(theta_c)
    frames = [v @ basis for v in baked]
    R = spec.outer_radius_mm
    out = dict(baked=baked, faces=faces, frames=frames, theta_c=theta_c,
               relocated=relocated, tip_global=max(reach),
               position=position, requested_theta_deg=requested_deg,
               reachable=(cand if wp is not None else list(range(len(baked)))),
               # how far every grain stands proud of the bond, in microns -- the
               # number you need before choosing a standoff
               protrusion_um=[(r - R) * 1000.0 for r in reach],
               e_r=e_r, e_t=e_t, e_z=e_z, basis=basis,
               r_tangent=None, r_ground=None, gov=None, clearance_um=None)
    if wp is not None:
        r_t, gov = ground_radius(frames, faces, wp.length_mm / 2.0, wp.width_mm / 2.0)
        out.update(r_tangent=r_t, gov=gov,
                   r_ground=r_t + clearance_um / 1000.0 + _WRITE_GUARD_MM,
                   clearance_um=(r_t - R) * 1000.0)
    return out


def write_rigid_wheel_inp(
    path: str,
    model: WheelModel,
    workpiece: Optional[WorkpieceBlock] = None,
    clearance_um: float = 0.0,
    wp_position: str = "centred",
    wp_position_deg: float = 0.0,
    bond_density_kg_m3: float = 2700.0,
    engage_window_mm: Optional[float] = None,
    model_name: str = "RIGID_WHEEL_AND_WORKPIECE",
    analysis=None,
    step_time_s: float = 0.0,
    surface_speed_mm_s: float = 30_000.0,
    include_bond: bool = True,
) -> dict:
    """Write the all-rigid-wheel geometry deck. Returns a summary dict.

    ``include_bond=False`` omits the bond rim shell, leaving the abrasive and the
    workpiece as the only bodies in the deck. On a single-grit deck the rim is
    2,812 of 2,928 rigid facets and never reaches the work -- the grit protrudes
    past it -- so it is 96% of the wheel mesh contributing nothing but general-
    contact facets to search. The rim's mass and inertia are KEPT either way:
    they stand for the spindle the grit is attached to, every reference-node DOF
    is velocity-driven regardless, and the tensor has to stay positive definite.
    """
    if not model.placements:
        raise ValueError("the wheel model has no grits placed")
    wp = workpiece
    spec = model.spec
    R = spec.outer_radius_mm

    if include_bond:
        shell_nodes, shell_quads, shell_groups = build_rim_shell(spec)
        shell_nodes = quantise(shell_nodes)
    else:
        shell_nodes = np.zeros((0, 3), dtype=float)
        shell_quads = []
        shell_groups = {name: [] for name in _SHELL_GROUPS}

    motion = None
    hybrid_info = None
    _entry_high = not bool(getattr(analysis, "rotation_reversed", False))
    pl = place_workpiece(model, wp, clearance_um, wp_position, wp_position_deg,
                         _entry_high)
    baked, faces, frames = pl["baked"], pl["faces"], pl["frames"]
    theta_c, relocated, tip_global = pl["theta_c"], pl["relocated"], pl["tip_global"]
    e_r, e_t, e_z, basis = pl["e_r"], pl["e_t"], pl["e_z"], pl["basis"]
    r_tangent, gov, r_ground = pl["r_tangent"], pl["gov"], pl["r_ground"]

    if wp is None:
        wp_nodes = wp_hexes = None
        wp_sets = {}
        nl = nw = nd = 0
        win = float("inf")
        engage = list(range(len(baked)))
    else:
        wp_nodes, wp_hexes, wp_sets = build_block_mesh(
            wp, e_r * r_ground, e_t, e_z, e_r)
        nl, nw, nd = wp.divisions()
        engage, win = engaging_grits(frames, wp, engage_window_mm)

    min_infeed_um = None
    if wp is not None and engage:
        hz = wp.width_mm / 2.0
        gaps = [(r_ground - float(frames[i][sel, 0].max())) * 1000.0
                for i in engage
                for sel in [(np.abs(frames[i][:, 1]) <= win)
                            & (np.abs(frames[i][:, 2]) <= hz)]
                if sel.any()]
        min_infeed_um = min(gaps) if gaps else None

    mass, inertia = rim_mass_properties(spec, bond_density_kg_m3)
    eig = np.linalg.eigvalsh(inertia)
    if mass <= 0 or eig.min() <= 0:
        raise ValueError("rim inertia tensor is not positive definite")

    # --- node / element numbering in the single WHEEL part -------------------
    n_shell_nodes = len(shell_nodes)
    n_shell_els = len(shell_quads)
    grit_node_base: list[int] = []
    off = n_shell_nodes
    for v in baked:
        grit_node_base.append(off)
        off += len(v)
    ref_node = off + 1                      # 1-based, after every surface node
    n_grit_facets = sum(len(f) for f in faces)

    with open(path, "w", encoding="ascii", newline="\n") as fh:
        w = fh.write
        arc = R * spec.sector_rad
        sag = arc * arc / (8.0 * R)
        rim = R - spec.inner_radius_mm
        w("*Heading\n")
        w("** %s\n" % model_name)
        # Conditional on analysis.enabled: with the analysis on, this deck DOES
        # carry a step, contact, boundary conditions and output -- and the
        # header claimed the opposite in a file whose *Step keyword sits 67,000
        # lines further down. A header that contradicts its own body costs more
        # with a sceptical reader than no header would.
        if analysis is not None and getattr(analysis, "enabled", False):
            w("** Run-ready: %s, the assembly, sets and surfaces, plus the step,\n"
              % ("two parts" if wp is not None else "one part"))
            w("** contact, boundary conditions and output requests.\n")
        else:
            w("** Geometry only: %s, the assembly, sets and surfaces.\n"
              % ("two parts" if wp is not None else "one part"))
            w("** No step, no interaction, no boundary condition, no load, "
              "no output.\n")
        w("** Units: mm, tonne, s, MPa, N.   Wheel axis = Z.\n")
        w("**\n")
        if include_bond:
            w("** WHEEL      : ONE discrete rigid body -- bond rim shell (R3D4) plus every\n")
            w("**              grit facet (R3D3), all tied to reference node %d at the\n"
              % ref_node)
        else:
            w("** ABRASIVE   : ONE discrete rigid body -- grit facets (R3D3) ONLY, no bond\n")
            w("**              rim. The abrasive and the workpiece are the only bodies in\n")
            w("**              this deck. Tied to reference node %d at the\n" % ref_node)
        w("**              origin on the wheel axis. Rotate the wheel with a single BC\n")
        # This line used to read "VR3 = -omega for +X-to-+Y travel", which is backwards:
        # a positive rotation about +Z carries +X toward +Y, so +X-to-+Y travel is
        # VR3 = +omega. The deck writes -omega and every other statement in the code
        # says the surface travels toward DECREASING theta. A header that contradicts
        # the BC beneath it is how the infeed sign survived for months, so it is now
        # stated once, correctly, and cross-checked by verify_rigid_deck.py.
        w("**              on WHEEL-1.WHEEL_REF. VR3 = %somega turns the surface toward\n"
          % ("-" if _entry_high else "+"))
        w("**              %s theta (%s), so grits arrive at\n"
          % ("decreasing" if _entry_high else "increasing",
             "+Y toward +X" if _entry_high else "+X toward +Y"))
        w("**              the workpiece from its %s-theta end.\n"
          % ("high" if _entry_high else "low"))
        if wp is not None:
            w("** WORKPIECE  : the only deformable part, C3D8R, material %s.\n"
              % wp.material)
        else:
            w("** No workpiece in this deck; add your own and position its ground face\n")
            w("** at r = %.6f mm to touch the tallest grit tip.\n" % tip_global)
        w("**\n")
        w("** wheel outer radius (mm)     : %.6f   (diameter %g)\n" % (R, spec.diameter_mm))
        w("** sector (deg)                : %.6f%s\n"
          % (spec.sector_deg, "   (full wheel)" if spec.is_full_circle else ""))
        w("** arc length (mm)             : %.6f%s\n"
          % (arc, "   (circumference)" if spec.is_full_circle else ""))
        w("** rim depth (mm)              : %.6f\n" % rim)
        if not spec.is_full_circle:
            w("** sagitta of the arc (um)     : %.3f  = %.0f%% of the rim depth\n"
              % (sag * 1000.0, 100.0 * sag / rim))
        w("** wheel axial width (mm)      : %g\n" % spec.width_mm)
        w("** grits                       : %d  (%d R3D3 facets)\n"
          % (len(model.placements), n_grit_facets))
        w("** rigid-body mass (tonne)     : %.6e\n" % mass)
        w("** tallest grit tip, whole arc : r = %.6f mm\n" % tip_global)
        if wp is not None:
            w("** grits able to engage        : %d  (within %.4f mm arc of the block)\n"
              % (len(engage), win))
            w("** governing grit in footprint : placement %d, r = %.6f mm\n"
              % (model.placements[gov].placement_id, r_tangent))
            w("** workpiece ground face at    : r = %.6f mm  (tangent, penetration 0)\n"
              % r_ground)
            w("** initial clearance (um)      : %g\n" % clearance_um)
            hl, hw, hd = wp.element_sizes()
            w("** workpiece (mm)              : %g x %g x %g\n"
              % (wp.length_mm, wp.width_mm, wp.depth_mm))
            w("** workpiece elements          : %d C3D8R  (%d x %d x %d)\n"
              % (len(wp_hexes), nl, nw, nd))
            w("** element size (um)           : %.4f cutting x %.4f axial x %.4f depth\n"
              % (hl * 1000.0, hw * 1000.0, hd * 1000.0))
            w("** governs the stable dt (um)  : %.4f  (the smallest of the three)\n"
              % (wp.min_element_size() * 1000.0))
            w("** workpiece centre at theta   : %.6f deg%s\n"
              % (math.degrees(theta_c),
                 "   (moved here: the nominal angle had no grit under it)"
                 if relocated else ""))
            w("**\n")
            w("** To cut, move WHEEL_REF radially inward by the depth of cut before or\n")
            w("** at the start of the step; the grits already graze the surface.\n")
        w("**\n")

        # ---------------- the wheel: one rigid part ----------------
        w("*Part, name=WHEEL\n*Node\n")
        for i, v in enumerate(shell_nodes, start=1):
            w(_node_line(i, v))
        for gi, v in enumerate(baked):
            b = grit_node_base[gi]
            for j, p in enumerate(v, start=1):
                w(_node_line(b + j, p))
        w("%d, 0., 0., 0.\n" % ref_node)

        if len(shell_quads):
            w("*Element, type=R3D4\n")
            for e, q in enumerate(shell_quads, start=1):
                w("%d, %d, %d, %d, %d\n"
                  % (e, q[0] + 1, q[1] + 1, q[2] + 1, q[3] + 1))
        w("*Element, type=R3D3\n")
        eid = n_shell_els
        grit_el_first: list[int] = []
        grit_el_last: list[int] = []
        for gi, f in enumerate(faces):
            b = grit_node_base[gi]
            grit_el_first.append(eid + 1)
            for tri in f:
                eid += 1
                w("%d, %d, %d, %d\n"
                  % (eid, b + tri[0] + 1, b + tri[1] + 1, b + tri[2] + 1))
            grit_el_last.append(eid)

        w("*Elset, elset=ES_WHEEL_ALL, generate\n1, %d, 1\n" % eid)
        w("*Elset, elset=ES_GRITS, generate\n%d, %d, 1\n" % (n_shell_els + 1, eid))
        for name in _SHELL_GROUPS:
            ids = shell_groups[name]
            if not ids:
                continue
            w("*Elset, elset=ES_BOND_%s\n" % name)
            _write_int_set(fh, [i + 1 for i in ids])
        if engage and len(engage) < len(faces):
            w("*Elset, elset=ES_GRITS_ENGAGE\n")
            ids = [e for i in engage
                   for e in range(grit_el_first[i], grit_el_last[i] + 1)]
            _write_int_set(fh, ids)

        # Mass and inertia live on their own elements so no contact offset is
        # introduced; density/thickness on *Rigid Body would shift the surface.
        w("*Element, type=MASS, elset=ES_WHEEL_MASS\n%d, %d\n" % (eid + 1, ref_node))
        w("*Mass, elset=ES_WHEEL_MASS\n%.9e,\n" % mass)
        w("*Element, type=ROTARYI, elset=ES_WHEEL_ROTI\n%d, %d\n" % (eid + 2, ref_node))
        w("*Rotary Inertia, elset=ES_WHEEL_ROTI\n")
        w("%.9e, %.9e, %.9e, %.9e, %.9e, %.9e\n"
          % (inertia[0, 0], inertia[1, 1], inertia[2, 2],
             inertia[0, 1], inertia[0, 2], inertia[1, 2]))

        w("*Nset, nset=WHEEL_REF\n%d,\n" % ref_node)
        w("*Rigid Body, ref node=WHEEL_REF, elset=ES_WHEEL_ALL\n")
        w("*Surface, type=ELEMENT, name=WHEEL_SURF\n")
        w("ES_GRITS, SPOS\n")
        if shell_groups.get("OUTER"):
            w("ES_BOND_OUTER, SPOS\n")
        w("*Surface, type=ELEMENT, name=GRITS_SURF\nES_GRITS, SPOS\n")
        w("*End Part\n")

        # ---------------- the workpiece: the only deformable part ----------------
        if wp is not None:
            w("*Part, name=%s\n*Node\n" % wp.name)
            for i, v in enumerate(wp_nodes, start=1):
                w(_node_line(i, v))
            w("*Element, type=C3D8R\n")
            for e, h in enumerate(wp_hexes, start=1):
                w("%d, " % e + ", ".join(str(int(n) + 1) for n in h) + "\n")
            w("*Elset, elset=WP_ALL, generate\n1, %d, 1\n" % len(wp_hexes))
            if analysis is not None and analysis.enabled:
                # The section controls carry ELEMENT DELETION, without which the
                # VUMAT's deletion flag does nothing and no chip ever separates.
                w("*Solid Section, elset=WP_ALL, controls=EC-1, material=%s\n,\n"
                  % wp.material)
            else:
                w("*Solid Section, elset=WP_ALL, material=%s\n,\n" % wp.material)
            for name, ids in wp_sets.items():
                w("*Nset, nset=%s\n" % name)
                _write_int_set(fh, [int(i) + 1 for i in ids])
            ground_els = [((i * nw) + j) * nd + 1 for i in range(nl) for j in range(nw)]
            w("*Elset, elset=ES_WP_GROUND\n")
            _write_int_set(fh, ground_els)
            w("*Surface, type=ELEMENT, name=WP_GROUND_SURF\nES_WP_GROUND, S1\n")
            w("*End Part\n")

        # ---------------- assembly ----------------
        w("*Assembly, name=ASSEMBLY\n")
        w("*Instance, name=WHEEL-1, part=WHEEL\n*End Instance\n")
        if wp is not None:
            w("*Instance, name=WP-1, part=%s\n*End Instance\n" % wp.name)
        w("*Nset, nset=A_WHEEL_REF, instance=WHEEL-1\nWHEEL_REF,\n")
        for name in wp_sets:
            w("*Nset, nset=A_%s, instance=WP-1\n%s,\n" % (name, name))
        if wp is not None and analysis is not None and analysis.enabled:
            # An assembly-level handle on the workpiece elements, so the step's
            # *Element Output can be scoped to them. Without it the request reaches
            # every element in the model and Abaqus warns twelve times that S, PEEQ,
            # SDV and STATUS mean nothing on R3D3, R3D4, MASS and ROTARYI -- true, and
            # noise that buries a real message.
            #
            # Written only when there is a step to use it: a geometry-only deck has no
            # output request, and adding an unused set would change decks that are
            # already validated against Abaqus.
            w("*Elset, elset=A_WP_ALL, instance=WP-1\nWP_ALL,\n")
        # An assembly *Surface must name an instance-qualified ELEMENT SET plus a
        # face or side identifier. Referring to a part-level *Surface instead gives
        # "The following sets were not found when generating the surface ..." and
        # aborts the import.
        w("*Surface, type=ELEMENT, name=A_WHEEL_SURF\n")
        w("WHEEL-1.ES_GRITS, SPOS\n")
        if shell_groups.get("OUTER"):
            w("WHEEL-1.ES_BOND_OUTER, SPOS\n")
        w("*Surface, type=ELEMENT, name=A_GRITS_SURF\nWHEEL-1.ES_GRITS, SPOS\n")
        if engage and len(engage) < len(faces):
            w("*Surface, type=ELEMENT, name=A_GRITS_ENGAGE_SURF\n")
            w("WHEEL-1.ES_GRITS_ENGAGE, SPOS\n")
        if wp is not None:
            w("*Surface, type=ELEMENT, name=A_WP_GROUND_SURF\n")
            w("WP-1.ES_WP_GROUND, S1\n")
        w("*End Assembly\n")

        # ---------------- material: workpiece only ----------------
        w("** The wheel is rigid, so it needs no section and no material.\n")
        if wp is not None and analysis is not None and analysis.enabled:
            # Run-ready: the JH-2 material, then the whole history definition, so the
            # deck can be submitted from the command line with no CAE step at all.
            from .analysis import write_section_and_material, write_step, wheel_motion
            # The motion is now computed BEFORE the material, because the
            # hybrid law's chip-thickness field is a function of the infeed
            # and the rotation sense. Nothing else about the order matters.
            motion = wheel_motion(analysis, theta_c, surface_speed_mm_s, R, step_time_s)
            chip = None
            if getattr(analysis, "material_model", "") == "hybrid":
                from .hybrid import chip_field
                chip = chip_field(
                    pl, motion, wp,
                    rotation_reversed=bool(
                        getattr(analysis, "rotation_reversed", False)),
                    dc_mm=analysis.hybrid.critical_depth_mm())
            hybrid_info = write_section_and_material(w, analysis, wp, chip)
            write_step(w, analysis, wp, motion, step_time_s,
                       engage_surface=bool(engage) and len(engage) < len(faces))
        elif wp is not None:
            w("** Placeholder elasticity for the workpiece so *Solid Section resolves;\n")
            w("** replace with your JH-2 *User Material / VUMAT block.\n")
            w("*Material, name=%s\n" % wp.material)
            w("*Density\n%.8e,\n" % (wp.density_kg_m3 * KGM3_TO_TONNE_MM3))
            w("*Elastic\n%.8e, %.4f\n" % (wp.youngs_modulus_mpa, wp.poisson_ratio))

    return {
        "path": path,
        "run_ready": bool(analysis is not None and analysis.enabled and wp is not None),
        "motion": motion,
        "size_bytes": os.path.getsize(path),
        "n_grits": len(model.placements),
        "n_grits_engaging": len(engage),
        "n_grit_facets": int(n_grit_facets),
        "n_bond_shell_quads": int(n_shell_els),
        "n_wheel_rigid_elements": int(n_shell_els + n_grit_facets),
        "n_workpiece_elements": int(len(wp_hexes)) if wp is not None else 0,
        "n_workpiece_nodes": int(len(wp_nodes)) if wp is not None else 0,
        "has_workpiece": wp is not None,
        "full_wheel": spec.is_full_circle,
        "wheel_ref_node": int(ref_node),
        "wheel_mass_tonne": mass,
        "wheel_inertia_zz_tonne_mm2": float(inertia[2, 2]),
        "outer_radius_mm": R,
        "sector_deg": spec.sector_deg,
        "arc_length_mm": R * spec.sector_rad,
        "rim_depth_mm": R - spec.inner_radius_mm,
        "sagitta_um": (R * spec.sector_rad) ** 2 / (8.0 * R) * 1000.0,
        "governing_grit_placement_id": (int(model.placements[gov].placement_id)
                                        if gov is not None else None),
        "tangent_radius_mm": r_tangent,
        "tallest_tip_whole_arc_mm": tip_global,
        "workpiece_ground_radius_mm": r_ground,
        "max_engaging_protrusion_um": ((r_tangent - R) * 1000.0
                                       if r_tangent is not None else None),
        # Radial infeed at which the *first* grain anywhere in the swept band reaches
        # the work. With a standoff this is not zero, and a depth of cut below it
        # means the wheel turns for the whole step and never touches anything -- the
        # failure that looks like a completed job with no damage and no contact force.
        "min_engaging_infeed_um": min_infeed_um,
        "hybrid": hybrid_info,
        "clearance_um": clearance_um,
        "theta_workpiece_deg": math.degrees(theta_c),
        "workpiece_relocated_to_tallest_grit": relocated,
    }
