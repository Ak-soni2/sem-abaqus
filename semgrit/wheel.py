"""Assembly of a grinding wheel (or an angular sector of one) from measured grains.

Coordinate system
-----------------
Wheel axis along **Z**. Radius lies in the XY plane, the sector spans
``theta = 0 .. sector_deg`` measured from +X, and the width runs from
``-width/2`` to ``+width/2`` along Z. Putting the axis on Z makes cyclic
symmetry about Z natural for sector models.

Lengths are millimetres, matching the usual Abaqus mm/N/tonne/MPa system. Grain
measurements arrive in microns and are converted on the way in.

Design notes / departures from the original notebook:

* An angular **sector** can be generated (30 deg, 180 deg, 360 deg, ...) instead
  of always a full wheel, and only the outer **rim** need be modelled. Together
  these cut element count by orders of magnitude, which is what makes a
  high-grain-count model tractable.
* Grains are placed by **Poisson-disk sampling** on the unwrapped rim surface, so
  they cannot overlap. The original drew independent uniform angles with no
  spacing constraint, so grains freely interpenetrated each other.
* Grain count is derived from abrasive **concentration** and the measured mean
  grain volume, rather than being an unexplained input.
* Everything is seeded. The original called ``np.random`` throughout with no
  seed, so no run could be reproduced.
* Grains are *not* boolean-unioned into the body. For an Abaqus contact or tie
  formulation they must remain separate bodies; fusing them would destroy the
  interfaces the analysis needs.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

from .grain3d import GrainSolid

UM_PER_MM = 1000.0

# C100 is defined as 4.4 carat/cm^3 of diamond, which is 25% by volume.
VOLUME_FRACTION_AT_C100 = 0.25

# Spacing target for Bridson sampling, as a fraction of the theoretical maximum
# density. The sampler actually achieves ~0.67 on this implementation (525 points
# in 50x10 at min_dist 0.8); asking for 0.50 deliberately overshoots the required
# count so the set can be subsampled down to exactly the requested density.
#
# The sampler must always run to completion: truncating it at N points leaves
# them clustered around the initial seed, which confined 10k grains to 1.4 deg of
# a 30 deg sector.
_POISSON_SPACING_TARGET = 0.50


class WheelBuildError(RuntimeError):
    pass


@dataclass
class WheelSpec:
    """Geometry of the wheel body."""

    diameter_mm: float
    width_mm: float
    sector_deg: float = 360.0
    """Angular extent to generate. 360 = full wheel, 180 = half, 30 = a wedge."""
    rim_depth_mm: Optional[float] = None
    """Radial thickness of the modelled body. ``None`` means solid to the axis
    (or to ``hub_diameter_mm``). Modelling only the rim is normally what you want:
    the grinding action is confined to a shallow surface layer."""
    hub_diameter_mm: float = 0.0

    # rim mesh density
    radial_divisions: int = 4
    axial_divisions: int = 8
    circumferential_divisions_per_deg: float = 1.0

    def __post_init__(self) -> None:
        if self.diameter_mm <= 0:
            raise ValueError("diameter_mm must be positive")
        if self.width_mm <= 0:
            raise ValueError("width_mm must be positive")
        if not (0.0 < self.sector_deg <= 360.0):
            raise ValueError("sector_deg must be in (0, 360]")
        if self.rim_depth_mm is not None and self.rim_depth_mm <= 0:
            raise ValueError("rim_depth_mm must be positive")
        if self.hub_diameter_mm < 0 or self.hub_diameter_mm >= self.diameter_mm:
            raise ValueError("hub_diameter_mm must be in [0, diameter_mm)")

    @property
    def outer_radius_mm(self) -> float:
        return self.diameter_mm / 2.0

    @property
    def inner_radius_mm(self) -> float:
        if self.rim_depth_mm is not None:
            return max(self.outer_radius_mm - self.rim_depth_mm, self.hub_diameter_mm / 2.0)
        return self.hub_diameter_mm / 2.0

    @property
    def is_full_circle(self) -> bool:
        return abs(self.sector_deg - 360.0) < 1e-9

    @property
    def sector_rad(self) -> float:
        return math.radians(self.sector_deg)

    @property
    def surface_area_mm2(self) -> float:
        """Area of the outer cylindrical surface within the sector."""
        return self.outer_radius_mm * self.sector_rad * self.width_mm

    @property
    def arc_length_mm(self) -> float:
        return self.outer_radius_mm * self.sector_rad

    def circumferential_divisions(self) -> int:
        return max(int(round(self.sector_deg * self.circumferential_divisions_per_deg)), 3)


@dataclass
class GrainPopulationSpec:
    """How many grains to place, how they sit, and how they are oriented."""

    areal_density_per_mm2: Optional[float] = None
    """Direct specification. Takes precedence over concentration."""
    concentration: Optional[float] = None
    """Abrasive concentration number; C100 = 25 vol%."""
    volume_fraction: Optional[float] = None
    """Explicit abrasive volume fraction, overriding ``concentration``."""

    max_grains: int = 200_000
    """Hard cap, so a slip of the pen cannot try to emit a billion elements."""

    protrusion_mean: float = 0.55
    protrusion_std: float = 0.12
    protrusion_min: float = 0.25
    protrusion_max: float = 0.85
    """Fraction of grain height standing above the bond surface."""

    max_tilt_deg: float = 35.0
    """Grains are seated tip-outward, then tilted by up to this angle in a random
    direction. Fully uniform orientation would bury many tips in the bond."""

    spacing_factor: float = 1.05
    """Poisson-disk minimum centre spacing, as a multiple of the local grain
    footprint radius sum. >1 leaves a bond gap between neighbours."""

    seed: int = 20260728

    def resolve_volume_fraction(self) -> Optional[float]:
        if self.volume_fraction is not None:
            return self.volume_fraction
        if self.concentration is not None:
            return VOLUME_FRACTION_AT_C100 * self.concentration / 100.0
        return None


@dataclass
class GrainPlacement:
    """One grain instanced onto the wheel."""

    placement_id: int
    shape_index: int
    """Index into the grain library; many placements share one shape, which is
    what keeps the Abaqus deck small."""
    translation_mm: np.ndarray      # (3,)
    rotation_axis: np.ndarray       # (3,) unit
    rotation_angle_deg: float
    theta_deg: float
    axial_mm: float
    radius_mm: float                # radial position of the grain centroid
    protrusion_mm: float
    bounding_radius_mm: float

    def matrix(self) -> np.ndarray:
        """4x4 homogeneous transform: rotate about the centroid, then translate."""
        return _compose(self.rotation_axis, self.rotation_angle_deg, self.translation_mm)


@dataclass
class WheelModel:
    """A meshed wheel body plus the grains placed on it."""

    spec: WheelSpec
    population: GrainPopulationSpec
    body_nodes: np.ndarray          # (N, 3) mm
    body_hexes: np.ndarray          # (M, 8) C3D8 connectivity
    placements: list[GrainPlacement]
    shapes: list[GrainSolid]
    requested_grains: int
    achieved_grains: int
    node_sets: dict[str, np.ndarray] = field(default_factory=dict)
    stats: dict = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    @property
    def n_body_nodes(self) -> int:
        return len(self.body_nodes)

    @property
    def n_body_elements(self) -> int:
        return len(self.body_hexes)

    def total_grain_tets(self) -> int:
        return sum(self.shapes[p.shape_index].n_tets for p in self.placements)

    def total_grain_faces(self) -> int:
        return sum(len(self.shapes[p.shape_index].faces) for p in self.placements)


# --------------------------------------------------------------------------
# Rotations
# --------------------------------------------------------------------------

def _compose(axis: np.ndarray, angle_deg: float, translation: np.ndarray) -> np.ndarray:
    m = np.eye(4)
    m[:3, :3] = _rotation_matrix(axis, math.radians(angle_deg))
    m[:3, 3] = translation
    return m


def _rotation_matrix(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    a = np.asarray(axis, dtype=np.float64)
    n = np.linalg.norm(a)
    if n < 1e-12:
        return np.eye(3)
    a = a / n
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    x, y, z = a
    return np.array(
        [
            [c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
            [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
            [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)],
        ]
    )


def matrix_to_axis_angle(r: np.ndarray) -> tuple[np.ndarray, float]:
    """Decompose a rotation matrix into axis + angle in degrees.

    Abaqus ``*Instance`` accepts exactly one translation and one rotation about an
    arbitrary axis, so every orientation must be expressed this way.
    """
    m = np.asarray(r, dtype=np.float64)[:3, :3]
    cos_t = float(np.clip((np.trace(m) - 1.0) / 2.0, -1.0, 1.0))
    angle = math.acos(cos_t)
    if angle < 1e-12:
        return np.array([0.0, 0.0, 1.0]), 0.0
    if abs(angle - math.pi) < 1e-6:
        # Near 180 deg the skew part vanishes; take the axis from M + I.
        a = m + np.eye(3)
        col = int(np.argmax(np.linalg.norm(a, axis=0)))
        axis = a[:, col]
        axis = axis / max(np.linalg.norm(axis), 1e-15)
        return axis, math.degrees(angle)
    axis = np.array([m[2, 1] - m[1, 2], m[0, 2] - m[2, 0], m[1, 0] - m[0, 1]])
    axis = axis / (2.0 * math.sin(angle))
    axis = axis / max(np.linalg.norm(axis), 1e-15)
    return axis, math.degrees(angle)


# --------------------------------------------------------------------------
# Rim mesh (structured C3D8 hexahedra)
# --------------------------------------------------------------------------

def build_rim_mesh(
    spec: WheelSpec, with_element_sets: bool = False
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """Structured hex mesh of the annular sector.

    An annular sector maps exactly onto a structured grid in
    ``(r, theta, z)``, so a well-shaped C3D8 mesh comes out directly with no
    tetrahedral meshing needed.

    Returns nodes, hex connectivity, and named node sets (outer surface, the two
    sector cut faces, the bore, and the axial faces).
    """
    r0, r1 = spec.inner_radius_mm, spec.outer_radius_mm
    if r1 - r0 <= 0:
        raise WheelBuildError("wheel body has zero radial thickness")

    n_r = max(int(spec.radial_divisions), 1)
    n_z = max(int(spec.axial_divisions), 1)
    n_t = spec.circumferential_divisions()

    radii = np.linspace(r0, r1, n_r + 1)
    zs = np.linspace(-spec.width_mm / 2.0, spec.width_mm / 2.0, n_z + 1)
    if spec.is_full_circle:
        thetas = np.linspace(0.0, 2 * math.pi, n_t, endpoint=False)
        n_t_nodes = n_t
    else:
        thetas = np.linspace(0.0, spec.sector_rad, n_t + 1)
        n_t_nodes = n_t + 1

    # node index = ((i_r * n_t_nodes) + i_t) * (n_z+1) + i_z
    def nid(i_r: int, i_t: int, i_z: int) -> int:
        return (i_r * n_t_nodes + (i_t % n_t_nodes)) * (n_z + 1) + i_z

    nodes = np.empty((( n_r + 1) * n_t_nodes * (n_z + 1), 3), dtype=np.float64)
    for i_r, r in enumerate(radii):
        for i_t, t in enumerate(thetas):
            x, y = r * math.cos(t), r * math.sin(t)
            for i_z, z in enumerate(zs):
                nodes[nid(i_r, i_t, i_z)] = (x, y, z)

    hexes: list[tuple[int, ...]] = []
    # Element indices per boundary, so element-based *Surface blocks can be
    # written. With the node ordering below the outer cylindrical face is S4,
    # the bore S6, the sector cut faces S3/S5, and the axial faces S1/S2.
    el_outer: list[int] = []
    el_bore: list[int] = []
    el_start: list[int] = []
    el_end: list[int] = []
    el_zmin: list[int] = []
    el_zmax: list[int] = []

    for i_r in range(n_r):
        for i_t in range(n_t):
            for i_z in range(n_z):
                e = len(hexes)
                if i_r == n_r - 1:
                    el_outer.append(e)
                if i_r == 0:
                    el_bore.append(e)
                if i_z == 0:
                    el_zmin.append(e)
                if i_z == n_z - 1:
                    el_zmax.append(e)
                if not spec.is_full_circle:
                    if i_t == 0:
                        el_start.append(e)
                    if i_t == n_t - 1:
                        el_end.append(e)
                # Abaqus C3D8: nodes 1-4 around the -z face, 5-8 around the +z
                # face with node 5 sitting above node 1.
                hexes.append(
                    (
                        nid(i_r, i_t, i_z),
                        nid(i_r + 1, i_t, i_z),
                        nid(i_r + 1, i_t + 1, i_z),
                        nid(i_r, i_t + 1, i_z),
                        nid(i_r, i_t, i_z + 1),
                        nid(i_r + 1, i_t, i_z + 1),
                        nid(i_r + 1, i_t + 1, i_z + 1),
                        nid(i_r, i_t + 1, i_z + 1),
                    )
                )

    hex_arr = np.asarray(hexes, dtype=np.int64)

    sets: dict[str, np.ndarray] = {}
    sets["WHEEL_OUTER"] = np.array(
        [nid(n_r, i_t, i_z) for i_t in range(n_t_nodes) for i_z in range(n_z + 1)],
        dtype=np.int64,
    )
    sets["WHEEL_BORE"] = np.array(
        [nid(0, i_t, i_z) for i_t in range(n_t_nodes) for i_z in range(n_z + 1)],
        dtype=np.int64,
    )
    sets["WHEEL_ZMIN"] = np.array(
        [nid(i_r, i_t, 0) for i_r in range(n_r + 1) for i_t in range(n_t_nodes)],
        dtype=np.int64,
    )
    sets["WHEEL_ZMAX"] = np.array(
        [nid(i_r, i_t, n_z) for i_r in range(n_r + 1) for i_t in range(n_t_nodes)],
        dtype=np.int64,
    )
    if not spec.is_full_circle:
        # The two radial cut faces. Apply cyclic-symmetry or symmetry BCs here.
        sets["SECTOR_FACE_START"] = np.array(
            [nid(i_r, 0, i_z) for i_r in range(n_r + 1) for i_z in range(n_z + 1)],
            dtype=np.int64,
        )
        sets["SECTOR_FACE_END"] = np.array(
            [nid(i_r, n_t, i_z) for i_r in range(n_r + 1) for i_z in range(n_z + 1)],
            dtype=np.int64,
        )
    if with_element_sets:
        # Generated in the same (i_r, i_z) order on both cut faces, so entry k of
        # SECTOR_FACE_START pairs with entry k of SECTOR_FACE_END. That makes
        # *Equation-based periodic boundary conditions straightforward.
        sets["_EL_OUTER"] = np.asarray(el_outer, dtype=np.int64)
        sets["_EL_BORE"] = np.asarray(el_bore, dtype=np.int64)
        sets["_EL_ZMIN"] = np.asarray(el_zmin, dtype=np.int64)
        sets["_EL_ZMAX"] = np.asarray(el_zmax, dtype=np.int64)
        if not spec.is_full_circle:
            sets["_EL_SECTOR_START"] = np.asarray(el_start, dtype=np.int64)
            sets["_EL_SECTOR_END"] = np.asarray(el_end, dtype=np.int64)
    return nodes, hex_arr, sets


def hex_volumes(nodes: np.ndarray, hexes: np.ndarray) -> np.ndarray:
    """Volume of each hexahedron by decomposition into 6 tetrahedra.

    Used to confirm the connectivity ordering is right: a negative value means
    the element is inverted and Abaqus would reject it.
    """
    # Node order 0-3 bottom, 4-7 top (matching the C3D8 convention above).
    tets = (
        (0, 1, 3, 4), (1, 2, 3, 6), (1, 3, 4, 6),
        (3, 4, 6, 7), (1, 4, 5, 6),
    )
    total = np.zeros(len(hexes))
    for a, b, c, d in tets:
        p0 = nodes[hexes[:, a]]
        e1 = nodes[hexes[:, b]] - p0
        e2 = nodes[hexes[:, c]] - p0
        e3 = nodes[hexes[:, d]] - p0
        total += np.einsum("ij,ij->i", np.cross(e1, e2), e3) / 6.0
    return total


# --------------------------------------------------------------------------
# Poisson-disk sampling on the unwrapped rim
# --------------------------------------------------------------------------

def jittered_grid_2d(
    width: float,
    height: float,
    count: int,
    min_dist: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, float]:
    """Stratified (jittered-grid) sampling with a guaranteed minimum separation.

    The rectangle is divided into cells of side ``c``; one point goes in each
    chosen cell, offset by a uniform jitter of at most ``(c - min_dist)/2`` in
    each axis. Two points in edge-adjacent cells can therefore approach no closer
    than ``c - 2*jitter = min_dist``, so the spacing constraint holds by
    construction with no rejection testing at all.

    This is O(count) and vectorised, where Bridson sampling is an interpreted
    loop that becomes impractical at the ~10^5 grains a full wheel needs. It
    also wraps exactly for a 360 deg model, since the cells tile the
    circumference.

    Returns ``(points, achieved_min_dist)``. ``achieved_min_dist`` is ``0`` when
    the requested count cannot respect ``min_dist``, in which case the caller
    must fall back to explicit overlap rejection.
    """
    if count <= 0:
        return np.zeros((0, 2)), 0.0
    aspect = width / max(height, 1e-12)
    n_x = max(int(round(math.sqrt(count * aspect))), 1)
    n_y = max(int(math.ceil(count / n_x)), 1)
    while n_x * n_y < count:
        n_x += 1
    cx, cy = width / n_x, height / n_y
    c = min(cx, cy)

    total = n_x * n_y
    chosen = rng.choice(total, size=min(count, total), replace=False)
    chosen = np.sort(chosen)
    ix, iy = chosen % n_x, chosen // n_x

    # When the cell is SMALLER than the required separation there is no room to
    # jitter, and `jitter` comes out exactly 0 -- which does not degrade
    # gracefully, it emits pure cell centres. On the shipped multi-abrasive
    # configuration (c = 0.008252 mm against min_dist = 0.013893 mm) that
    # produced an 8 x 2 lattice: all twelve grits at exactly two axial
    # positions, z = +/-4.1259 um, 374 nm from the free face of a 9 um block.
    # Two of them cut, both in the same lane, so the "multi-abrasive" result was
    # one unconfined edge scratch and no groove-groove interaction was
    # representable at any depth of cut.
    #
    # A lattice is a legitimate fallback, but it has to be reported as one:
    # `achieved` returning min_dist while the sampler ignored min_dist entirely
    # is what let it pass unnoticed.
    jitter = max((c - min_dist) / 2.0, 0.0)
    # Degenerate only when the lattice collapses the AXIAL spread, which is the
    # case that breaks the physics: with n_y rows the grits occupy exactly n_y
    # distinct z, and at n_y = 2 every grit sits on one of two lines. On a full
    # wheel at high density n_y is large, the lattice is a perfectly good
    # sampling of the face, and disturbing it would move a frozen reference deck
    # for no gain -- so the test is on n_y, not on `jitter` alone.
    degenerate = (jitter <= 0.0 and c < min_dist
                  and n_y < max(3, len(chosen) // 4))
    achieved = min_dist if jitter > 0 or c >= min_dist else 0.0
    jx = rng.uniform(-jitter, jitter, size=len(chosen)) if jitter > 0 else np.zeros(len(chosen))
    jy = rng.uniform(-jitter, jitter, size=len(chosen)) if jitter > 0 else np.zeros(len(chosen))

    x = (ix + 0.5) * cx + jx
    y = (iy + 0.5) * cy + jy
    if degenerate:
        # Break the lattice. The separation cannot be honoured whatever we do --
        # the caller asked for more grits than fit -- but a jittered spread over
        # the available cell is strictly better than a lattice, because it does
        # not concentrate every grit onto two coordinates. Half a cell, which is
        # the most that can be applied without leaving the cell.
        jx = rng.uniform(-cx / 2.0, cx / 2.0, size=len(chosen))
        jy = rng.uniform(-cy / 2.0, cy / 2.0, size=len(chosen))
        x = (ix + 0.5) * cx + jx
        y = (iy + 0.5) * cy + jy
    np.clip(y, 0.0, height, out=y)
    x %= width
    return np.column_stack([x, y]), achieved


def poisson_disk_2d(
    width: float,
    height: float,
    min_dist: float,
    rng: np.random.Generator,
    periodic_x: bool = False,
    k: int = 30,
    max_points: Optional[int] = None,
) -> np.ndarray:
    """Bridson Poisson-disk sampling in a rectangle.

    ``periodic_x`` wraps the first axis, which is what makes a full 360 deg wheel
    seamless at theta = 0 instead of leaving a visible seam or overlapping grains
    across it.
    """
    if min_dist <= 0:
        raise ValueError("min_dist must be positive")
    cell = min_dist / math.sqrt(2.0)
    gw = max(int(math.ceil(width / cell)), 1)
    gh = max(int(math.ceil(height / cell)), 1)
    grid = -np.ones((gw, gh), dtype=np.int64)

    pts: list[tuple[float, float]] = []
    active: list[int] = []

    def dx(a: float, b: float) -> float:
        d = abs(a - b)
        return min(d, width - d) if periodic_x else d

    def fits(x: float, y: float) -> bool:
        gx, gy = int(x / cell), int(y / cell)
        for ix in range(gx - 2, gx + 3):
            jx = ix % gw if periodic_x else ix
            if jx < 0 or jx >= gw:
                continue
            for iy in range(gy - 2, gy + 3):
                if iy < 0 or iy >= gh:
                    continue
                idx = grid[jx, iy]
                if idx < 0:
                    continue
                px, py = pts[idx]
                if dx(x, px) ** 2 + (y - py) ** 2 < min_dist * min_dist:
                    return False
        return True

    def add(x: float, y: float) -> None:
        pts.append((x, y))
        grid[min(int(x / cell), gw - 1), min(int(y / cell), gh - 1)] = len(pts) - 1
        active.append(len(pts) - 1)

    add(float(rng.uniform(0, width)), float(rng.uniform(0, height)))
    while active:
        if max_points is not None and len(pts) >= max_points:
            break
        ai = int(rng.integers(len(active)))
        idx = active[ai]
        ox, oy = pts[idx]
        placed = False
        for _ in range(k):
            ang = float(rng.uniform(0, 2 * math.pi))
            rad = float(min_dist * (1.0 + rng.random()))
            nx, ny = ox + rad * math.cos(ang), oy + rad * math.sin(ang)
            if periodic_x:
                nx %= width
            elif not (0 <= nx < width):
                continue
            if not (0 <= ny < height):
                continue
            if fits(nx, ny):
                add(nx, ny)
                placed = True
                break
        if not placed:
            active.pop(ai)
    return np.asarray(pts, dtype=np.float64)


# --------------------------------------------------------------------------
# Assembly
# --------------------------------------------------------------------------

def estimate_grain_count(
    spec: WheelSpec, population: GrainPopulationSpec, shapes: Sequence[GrainSolid]
) -> tuple[int, dict]:
    """Number of grains to place on the sector surface.

    With an explicit areal density this is just density x area. From a
    concentration it is a documented first-order estimate: the volumetric number
    density is the abrasive volume fraction divided by the mean grain volume, and
    the grains that reach the surface are those whose centres lie within roughly
    one mean grain height of it.
    """
    area = spec.surface_area_mm2
    info: dict = {"surface_area_mm2": area}

    if population.areal_density_per_mm2 is not None:
        density = float(population.areal_density_per_mm2)
        info["source"] = "areal_density"
    else:
        vf = population.resolve_volume_fraction()
        if vf is None:
            raise WheelBuildError(
                "specify areal_density_per_mm2, concentration, or volume_fraction"
            )
        if not shapes:
            raise WheelBuildError("no grain shapes available to size the population")
        vols_mm3 = np.array([s.mesh_volume_um3 for s in shapes]) / (UM_PER_MM ** 3)
        heights_mm = np.array([s.height_um for s in shapes]) / UM_PER_MM
        mean_v = float(vols_mm3.mean())
        mean_h = float(heights_mm.mean())
        if mean_v <= 0:
            raise WheelBuildError("mean grain volume is non-positive")
        n_vol = vf / mean_v            # grains per mm^3
        density = n_vol * mean_h       # grains per mm^2 in the surface layer
        info.update(
            {
                "source": "concentration",
                "volume_fraction": vf,
                "mean_grain_volume_mm3": mean_v,
                "mean_grain_height_mm": mean_h,
                "volumetric_density_per_mm3": n_vol,
            }
        )
    info["areal_density_per_mm2"] = density
    uncapped = int(round(density * area))
    n = max(min(uncapped, population.max_grains), 0)
    info["uncapped_grains"] = uncapped
    info["requested_grains"] = n
    info["capped"] = n < uncapped
    if n < uncapped:
        info["cap_note"] = (
            f"concentration implies {uncapped} grains ({density:.0f}/mm2) but "
            f"max_grains={population.max_grains}; the model covers only "
            f"{100.0 * n / uncapped:.1f}% of the implied abrasive. Raise "
            f"--max-grains, shrink --sector, or give --areal-density directly."
        )
    return n, info


def build_wheel(
    spec: WheelSpec,
    shapes: Sequence[GrainSolid],
    population: Optional[GrainPopulationSpec] = None,
) -> WheelModel:
    """Place measured grains on a wheel (or sector) and mesh the body."""
    population = population or GrainPopulationSpec(concentration=100.0)
    shapes = list(shapes)
    if not shapes:
        raise WheelBuildError("no grain shapes supplied")

    warnings: list[str] = []
    nodes, hexes, sets = build_rim_mesh(spec, with_element_sets=True)

    vols = hex_volumes(nodes, hexes)
    if (vols <= 0).any():
        raise WheelBuildError(
            f"{int((vols <= 0).sum())} inverted hexahedra in the rim mesh"
        )

    requested, info = estimate_grain_count(spec, population, shapes)
    if info.get("cap_note"):
        warnings.append(info["cap_note"])

    # Footprint radius of a grain in the tangential/axial plane, in mm.
    foot = np.array(
        [
            0.5 * max(float(s.extent_um()[0]), float(s.extent_um()[1])) / UM_PER_MM
            for s in shapes
        ]
    )
    mean_foot = float(foot.mean())

    # Separation is enforced on the grains' bounding spheres, using the *largest*
    # grain in the library. Spacing on the mean grain size let two large
    # neighbours interpenetrate (8474 overlapping pairs on a 30 deg sector).
    max_radius = float(max(s.bounding_radius_um for s in shapes) / UM_PER_MM)
    required_gap = 2.0 * max_radius * population.spacing_factor

    arc = spec.arc_length_mm
    area = max(spec.surface_area_mm2, 1e-12)
    rng = np.random.default_rng(population.seed)

    samples, achieved_gap = jittered_grid_2d(
        width=arc,
        height=spec.width_mm,
        count=requested,
        min_dist=required_gap,
        rng=rng,
    )
    # Guaranteed separation means the rejection pass below cannot fire; it is kept
    # as a safety net for the saturated case.
    guaranteed = achieved_gap >= required_gap > 0
    spacing_limited = not guaranteed and requested > 0
    if spacing_limited:
        cell_side = math.sqrt(area / max(requested, 1))
        warnings.append(
            f"requested density needs {cell_side:.5f} mm cells but the largest "
            f"grain needs {required_gap:.5f} mm clearance; grains will be "
            f"rejected where they collide"
        )
    if len(samples) < requested:
        warnings.append(
            f"grid provided {len(samples)} of {requested} requested positions"
        )
    min_dist = max(achieved_gap, required_gap)

    # Centroid-relative vertices in mm, computed once per shape rather than once
    # per placement.
    local_cache = [(s.vertices - s.centroid_um) / UM_PER_MM for s in shapes]

    placements: list[GrainPlacement] = []
    # Uniform grid for the bounding-sphere rejection test.
    cell = max(2.0 * max_radius, 1e-9)
    grid: dict[tuple[int, int, int], list[int]] = {}
    accepted_centres: list[np.ndarray] = []
    accepted_radii: list[float] = []
    n_rejected = 0

    def cell_of(p: np.ndarray) -> tuple[int, int, int]:
        return (int(math.floor(p[0] / cell)), int(math.floor(p[1] / cell)),
                int(math.floor(p[2] / cell)))

    def overlaps(centre: np.ndarray, radius: float) -> bool:
        cx, cy, cz = cell_of(centre)
        for i in (cx - 1, cx, cx + 1):
            for j in (cy - 1, cy, cy + 1):
                for k in (cz - 1, cz, cz + 1):
                    for idx in grid.get((i, j, k), ()):
                        d = float(np.linalg.norm(centre - accepted_centres[idx]))
                        if d < radius + accepted_radii[idx]:
                            return True
        return False

    R = spec.outer_radius_mm
    for s_arc, z_ax in samples:
        shape_index = int(rng.integers(len(shapes)))
        shape = shapes[shape_index]

        theta = (s_arc / R) if R > 0 else 0.0
        radial = np.array([math.cos(theta), math.sin(theta), 0.0])

        # Seat the grain tip outward: map local +Z onto the outward radial
        # direction, spin randomly about it, then tilt by a bounded angle.
        spin = float(rng.uniform(0, 2 * math.pi))
        tilt = math.radians(float(rng.uniform(0, population.max_tilt_deg)))
        tilt_dir = float(rng.uniform(0, 2 * math.pi))

        r_align = _align_z_to(radial)
        r_spin = _rotation_matrix(radial, spin)
        tilt_axis = np.cross(radial, _perp(radial, tilt_dir))
        r_tilt = _rotation_matrix(tilt_axis, tilt)
        rot = r_tilt @ r_spin @ r_align

        frac = _truncated_normal(
            rng,
            population.protrusion_mean,
            population.protrusion_std,
            population.protrusion_min,
            population.protrusion_max,
        )
        protrusion = frac * shape.height_um / UM_PER_MM
        target = R + protrusion

        # Radial and tangential components of every grain vertex, expressed
        # without ever forming the rotated point cloud: (R v).e = v.(R^T e), so
        # two matrix-vector products replace an N x 3 matrix product. At 10^5
        # grains that difference dominates the runtime.
        local = local_cache[shape_index]
        a = local @ (rot.T @ radial)          # outward component
        tangent = np.array([-radial[1], radial[0], 0.0])
        b = local @ (rot.T @ tangent)         # tangential component

        # Seat the grain so its furthest point sits exactly `protrusion` above the
        # bond surface. Projecting onto the radial direction alone is not enough:
        # the surface is curved, so a vertex offset tangentially by v_t sits
        # v_t^2 / 2r further out than its radial projection suggests -- about
        # 1 nm here, small but enough to make an exactness check fail.
        r_centre = target - float(a.max())
        for _ in range(8):
            reached = float(np.sqrt((r_centre + a) ** 2 + b ** 2).max())
            delta = target - reached
            if abs(delta) < 1e-12:
                break
            r_centre += delta

        axis, angle = matrix_to_axis_angle(rot)
        centre = radial * r_centre + np.array([0.0, 0.0, z_ax - spec.width_mm / 2.0])

        bound_r = float(shape.bounding_radius_um / UM_PER_MM)
        if not guaranteed:
            if overlaps(centre, bound_r):
                n_rejected += 1
                continue
            grid.setdefault(cell_of(centre), []).append(len(accepted_centres))
            accepted_centres.append(centre)
            accepted_radii.append(bound_r)

        placements.append(
            GrainPlacement(
                placement_id=len(placements) + 1,
                shape_index=shape_index,
                translation_mm=centre,
                rotation_axis=axis,
                rotation_angle_deg=angle,
                theta_deg=math.degrees(theta),
                axial_mm=float(z_ax - spec.width_mm / 2.0),
                radius_mm=float(r_centre),
                protrusion_mm=float(protrusion),
                bounding_radius_mm=bound_r,
            )
        )

    if n_rejected:
        warnings.append(
            f"rejected {n_rejected} of {len(samples)} candidate positions where the "
            f"grain bounding sphere intersected a neighbour (conservative test: "
            f"spheres can touch where the polyhedra do not)"
        )

    stats = dict(info)
    stats.update(
        {
            "achieved_grains": len(placements),
            "achieved_areal_density_per_mm2": (
                len(placements) / spec.surface_area_mm2
                if spec.surface_area_mm2 > 0
                else float("nan")
            ),
            "min_spacing_mm": min_dist,
            "spacing_limited_by_grain_size": bool(spacing_limited),
            "n_rejected_overlaps": int(n_rejected),
            "mean_grain_footprint_mm": 2 * mean_foot,
            "body_nodes": int(len(nodes)),
            "body_elements": int(len(hexes)),
            "min_hex_volume_mm3": float(vols.min()),
            "sector_deg": spec.sector_deg,
            "rim_depth_mm": spec.outer_radius_mm - spec.inner_radius_mm,
        }
    )

    return WheelModel(
        spec=spec,
        population=population,
        body_nodes=nodes,
        body_hexes=hexes,
        placements=placements,
        shapes=shapes,
        requested_grains=requested,
        achieved_grains=len(placements),
        node_sets=sets,
        stats=stats,
        warnings=warnings,
    )


def _align_z_to(target: np.ndarray) -> np.ndarray:
    """Rotation taking +Z onto ``target``."""
    z = np.array([0.0, 0.0, 1.0])
    t = target / max(np.linalg.norm(target), 1e-15)
    c = float(np.dot(z, t))
    if c > 1 - 1e-12:
        return np.eye(3)
    if c < -1 + 1e-12:
        return _rotation_matrix(np.array([1.0, 0.0, 0.0]), math.pi)
    axis = np.cross(z, t)
    return _rotation_matrix(axis, math.acos(max(-1.0, min(1.0, c))))


def _perp(v: np.ndarray, phase: float) -> np.ndarray:
    """A unit vector perpendicular to ``v``, selected by ``phase``."""
    a = np.array([0.0, 0.0, 1.0]) if abs(v[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(v, a)
    e1 /= max(np.linalg.norm(e1), 1e-15)
    e2 = np.cross(v, e1)
    e2 /= max(np.linalg.norm(e2), 1e-15)
    return math.cos(phase) * e1 + math.sin(phase) * e2


def _truncated_normal(
    rng: np.random.Generator, mean: float, std: float, lo: float, hi: float
) -> float:
    for _ in range(64):
        v = float(rng.normal(mean, std))
        if lo <= v <= hi:
            return v
    return float(min(max(mean, lo), hi))


def check_grain_overlaps(model: WheelModel) -> dict:
    """Bounding-sphere overlap audit of the placed grains.

    Poisson-disk sampling constrains centres on the unwrapped surface; this
    verifies the constraint actually holds in 3D, where grains also differ in
    radial position and size.
    """
    if not model.placements:
        return {"n_pairs_checked": 0, "n_overlapping": 0, "worst_overlap_mm": 0.0}

    centres = np.array([p.translation_mm for p in model.placements])
    radii = np.array([p.bounding_radius_mm for p in model.placements])
    max_r = float(radii.max())

    # A KD-tree neighbour query, not a sweep along one axis. Sorting by x prunes
    # nothing on a narrow sector, where every grain shares almost the same x, so
    # that approach degenerated to O(n^2) and took minutes at 30k grains.
    from scipy.spatial import cKDTree

    tree = cKDTree(centres)
    candidates = tree.query_pairs(r=2.0 * max_r, output_type="ndarray")
    pairs = int(len(candidates))
    n_overlap = 0
    worst = 0.0
    if pairs:
        ia, ib = candidates[:, 0], candidates[:, 1]
        d = np.linalg.norm(centres[ia] - centres[ib], axis=1)
        gap = d - (radii[ia] + radii[ib])
        bad = gap < 0
        n_overlap = int(bad.sum())
        if n_overlap:
            worst = float(gap[bad].min())
    return {
        "n_pairs_checked": pairs,
        "n_overlapping": n_overlap,
        "worst_overlap_mm": abs(worst),
        "note": (
            "bounding-sphere test; spheres may intersect where the actual "
            "polyhedra do not"
        ),
    }
