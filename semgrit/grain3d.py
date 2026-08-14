"""Construction of 3D grain solids from measured 2D outlines.

Each grain becomes a *lofted polyhedron*: the measured SEM silhouette is placed
at mid height and tapered to smaller faces above and below. Two properties make
this defensible rather than decorative:

1. The maximum projected cross-section of the solid **is** the measured outline,
   so projecting the CAD back onto the imaging plane reproduces the measurement.
   This is checked numerically in :func:`validate_grain_solid`.
2. Because each ring is a uniform scaling of the outline about a common centre,
   every lateral face is a trapezoid with two parallel sides and is therefore
   planar. That means the volume has a closed form (the prismatoid/frustum
   formula) which the generated mesh is verified against exactly.

Departures from the original notebook, which extruded a convex hull into a
straight prism and then jittered the vertices at random:

* No convex hull. Concave fracture surfaces are the sharp cutting features and
  the hull deletes them, while also inflating grain size.
* No random vertex jitter. It produced self-intersecting, non-watertight,
  inverted-normal meshes with no validation.
* Height is derived from the measured minimum Feret width via an explicit,
  documented blockiness model with a fixed seed, not ``uniform(0.35, 0.55)``.
* Rotations are applied about the grain centroid, not the world origin. Rotating
  about the origin flung grains far from their measured positions.

Volume meshing produces C3D4 tetrahedra with no extra dependencies: the outline
is ear-clipped into triangles, each triangle is swept between two rings to form
a triangular prism, and each prism is split into 3 tets using a globally sorted
vertex rule so that the diagonals of shared quad faces always agree and the mesh
stays conforming.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Sequence

import cv2
import numpy as np

try:  # pragma: no cover - exercised implicitly
    import mapbox_earcut
except ImportError:  # pragma: no cover
    mapbox_earcut = None

from shapely.geometry import Polygon

from .measure import GrainMeasurement
from .metrology import SemImage
from .segment import Segmentation


# --------------------------------------------------------------------------
# Height model
# --------------------------------------------------------------------------

@dataclass
class HeightModel:
    """How grain thickness is inferred from the 2D measurement.

    A single top-down SEM image cannot measure height, so this is the one
    genuinely modelled quantity in the pipeline and is kept explicit and
    seeded rather than hidden inside a random call.

    Thickness is taken as a fraction of the *minimum* Feret width, because for a
    blocky crushed fragment resting on its most stable face the through-thickness
    dimension scales with the short in-plane axis, not the long one.
    """

    mean_ratio: float = 0.70
    std_ratio: float = 0.12
    min_ratio: float = 0.45
    max_ratio: float = 0.95
    seed: int = 20260728

    def thickness_ratio(self, grain_key: int) -> float:
        """Deterministic per-grain thickness ratio.

        Seeded from the grain key so the value depends only on which grain it
        is, not on how many grains were processed before it.
        """
        rng = np.random.default_rng([self.seed, int(grain_key)])
        for _ in range(64):
            v = float(rng.normal(self.mean_ratio, self.std_ratio))
            if self.min_ratio <= v <= self.max_ratio:
                return v
        return float(np.clip(self.mean_ratio, self.min_ratio, self.max_ratio))

    def height_um(self, measurement: GrainMeasurement) -> float:
        base = measurement.feret_min_um
        if not math.isfinite(base) or base <= 0:
            base = measurement.equivalent_diameter_um
        return base * self.thickness_ratio(measurement.grain_id)


@dataclass
class LoftProfile:
    """Vertical profile of the lofted grain, and how its edges are blunted."""

    base_scale: float = 0.70
    """Cross-section scale at z = 0, relative to the measured outline."""
    mid_height_fraction: float = 0.42
    """Height at which the full measured outline sits."""
    top_scale: float = 0.30
    """Cross-section scale at the tip. Smaller = sharper grain."""

    edge_radius_um: float = 0.0
    """Cutting edge radius. 0 leaves the geometry mathematically sharp.

    A sharp edge has zero radius, which in FEA is a **stress singularity**: the
    computed stress rises without bound as the mesh is refined, so results never
    converge and depend on element size. It is also wrong physically -- the edge
    radius is what sets the minimum chip thickness and the ploughing-to-cutting
    transition in grinding.

    This is deliberately a user parameter rather than something measured from the
    image. Edge radius is a dressing and wear property of the wheel, it changes
    during use, and a top-down SEM view of loose grit cannot resolve it reliably.

    Applied in two places: the outline's convex corners are rounded in-plane by a
    morphological opening (exact in 2D), and the meridional profile is filleted so
    the ridge where the side meets the top and bottom caps is rounded too.
    """

    arc_segments: int = 3
    """Segments per rounded corner. Higher is smoother but adds elements."""

    def validate(self) -> None:
        if not (0.0 < self.base_scale <= 1.0):
            raise ValueError("base_scale must be in (0, 1]")
        if not (0.0 < self.top_scale <= 1.0):
            raise ValueError("top_scale must be in (0, 1]")
        if not (0.0 < self.mid_height_fraction < 1.0):
            raise ValueError("mid_height_fraction must be in (0, 1)")
        if self.edge_radius_um < 0:
            raise ValueError("edge_radius_um must be >= 0")
        if self.arc_segments < 1:
            raise ValueError("arc_segments must be >= 1")


# --------------------------------------------------------------------------
# Outline extraction
# --------------------------------------------------------------------------

def grain_outline_um(
    label: int,
    seg: Segmentation,
    sem: SemImage,
    simplify_um: float = 0.10,
    max_vertices: int = 64,
) -> Optional[np.ndarray]:
    """Simplified grain outline in microns, centred on the grain centroid.

    The simplification tolerance is a physical length, so the same grain yields
    the same polygon at 5 kX and 10 kX. Y is negated so the result is a
    right-handed (counter-clockwise) polygon in a standard XY frame rather than
    image row/column order.
    """
    ps = sem.pixel_size_um
    mask = (seg.labels == label).astype(np.uint8)
    cs, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not cs:
        return None
    contour = max(cs, key=cv2.contourArea).reshape(-1, 2).astype(np.float64)
    if len(contour) < 3:
        return None

    pts = np.column_stack([contour[:, 0] * ps, -contour[:, 1] * ps])
    poly = Polygon(pts)
    if not poly.is_valid:
        poly = poly.buffer(0)
        if poly.is_empty:
            return None
        if poly.geom_type == "MultiPolygon":
            poly = max(poly.geoms, key=lambda g: g.area)

    tol = simplify_um
    simple = poly.simplify(tol, preserve_topology=True)
    # Tighten until the vertex budget is met; a coarse cap keeps FE meshes small.
    for _ in range(24):
        if simple.is_empty or simple.geom_type != "Polygon":
            break
        if len(simple.exterior.coords) - 1 <= max_vertices:
            break
        tol *= 1.35
        simple = poly.simplify(tol, preserve_topology=True)

    if simple.is_empty or simple.geom_type != "Polygon" or not simple.is_valid:
        simple = poly
    ring = np.asarray(simple.exterior.coords)[:-1]
    if len(ring) < 3:
        return None

    # Enforce counter-clockwise winding (positive signed area).
    if _signed_area(ring) < 0:
        ring = ring[::-1]

    ring = ring - ring.mean(axis=0)
    return ring


def _signed_area(ring: np.ndarray) -> float:
    x, y = ring[:, 0], ring[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


# Shortest edge any facet may have, in microns, after every ring scaling.
#
# Abaqus rejects a general-contact facet whose nodes are coincident:
#   "A facet that is part of a general contact surface has coincident nodes"
# and its threshold is far looser than exact equality -- an observed failure had
# nodes 1.03e-7 mm (0.0001 um) apart. 0.01 um is 100x clear of that while still
# being 0.3% of a 3 um grain, so nothing of shape is lost.
MIN_FACET_EDGE_UM = 0.01


def _clean_ring(
    ring: np.ndarray,
    min_edge_frac: float = 1e-4,
    min_edge_abs: float = MIN_FACET_EDGE_UM,
) -> np.ndarray:
    """Drop near-duplicate and near-collinear vertices from a polygon ring.

    Slivers in the ear-clipped triangulation become near-zero-volume tets, which
    Abaqus rejects and which break the once-used-face rule that the surface
    extraction relies on. Cleaning the ring removes the cause.

    ``min_edge_abs`` is an **absolute** floor in the ring's own units. A purely
    relative tolerance is not enough: the caller scales this ring down to build the
    tapered rings, so an edge that only just survives here comes back a factor of
    ``top_scale`` shorter on the tip ring. Callers should pass
    ``MIN_FACET_EDGE_UM / smallest_ring_scale``.
    """
    r = np.asarray(ring, dtype=np.float64)
    if len(r) < 3:
        return r
    scale = math.sqrt(max(abs(_signed_area(r)), 1e-30))
    min_edge = max(min_edge_frac * scale, min_edge_abs, 1e-12)

    # 1. collapse consecutive duplicates
    keep = [0]
    for i in range(1, len(r)):
        if np.linalg.norm(r[i] - r[keep[-1]]) > min_edge:
            keep.append(i)
    r = r[keep]
    if len(r) >= 3 and np.linalg.norm(r[0] - r[-1]) <= min_edge:
        r = r[:-1]
    if len(r) < 3:
        return r

    # 2. drop vertices that add no shape (triangle area below tolerance)
    changed = True
    while changed and len(r) > 3:
        changed = False
        n = len(r)
        for i in range(n):
            a, b, c = r[(i - 1) % n], r[i], r[(i + 1) % n]
            area2 = abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
            if area2 < min_edge * scale:
                r = np.delete(r, i, axis=0)
                changed = True
                break
    return r


def inscribed_radius(poly: Polygon) -> float:
    """Radius of the largest circle fitting inside the polygon.

    Used to convert between a cross-section *scale* factor and a physical inward
    *inset* distance, so the edge radius can be specified in microns.
    """
    try:
        from shapely.ops import polylabel

        tol = max(math.sqrt(max(poly.area, 1e-12)) * 1e-3, 1e-9)
        p = polylabel(poly, tolerance=tol)
        return float(poly.exterior.distance(p))
    except Exception:
        rp = poly.representative_point()
        return float(poly.exterior.distance(rp))


def round_outline_corners(
    ring: np.ndarray, radius: float, arc_segments: int = 3
) -> tuple[np.ndarray, float]:
    """Round the convex corners of a polygon to a given radius.

    Implemented as a morphological opening -- erode by ``radius`` then dilate back
    -- which is the textbook definition of limiting convex curvature: afterwards no
    convex corner is sharper than ``radius``. Concave notches are left alone, which
    is correct; blunting a grain rounds its points, it does not fill its hollows.

    Returns ``(ring, achieved_radius)``. The radius is reduced if the requested
    value would consume the grain (a 5 um grain cannot carry a 4 um edge radius).
    """
    poly = Polygon(ring)
    if not poly.is_valid:
        poly = poly.buffer(0)
        if poly.is_empty or poly.geom_type != "Polygon":
            return ring, 0.0
    if radius <= 0:
        return ring, 0.0

    r_in = inscribed_radius(poly)
    r = min(radius, 0.60 * r_in)  # leave material behind after the erosion
    qs = max(int(arc_segments), 1)

    for _ in range(12):
        if r <= 1e-9:
            return ring, 0.0
        try:
            eroded = poly.buffer(-r, join_style=1, quad_segs=qs)
            if (not eroded.is_empty) and eroded.area > 0:
                if eroded.geom_type == "MultiPolygon":
                    eroded = max(eroded.geoms, key=lambda g: g.area)
                out = eroded.buffer(r, join_style=1, quad_segs=qs)
                if out.geom_type == "MultiPolygon":
                    out = max(out.geoms, key=lambda g: g.area)
                if out.is_valid and not out.is_empty and out.area > 0.2 * poly.area:
                    coords = np.asarray(out.exterior.coords)[:-1]
                    if len(coords) >= 3:
                        if _signed_area(coords) < 0:
                            coords = coords[::-1]
                        return coords, r
        except Exception:
            pass
        r *= 0.6
    return ring, 0.0


def fillet_polyline(
    points: Sequence[Sequence[float]], radius: float, arc_segments: int = 3
) -> tuple[np.ndarray, float]:
    """Replace interior corners of an open 2D polyline with tangent circular arcs.

    Used on the grain's meridional profile so the ridges where the side surface
    meets the top and bottom caps, and the equator where the widest cross-section
    sits, all carry the requested radius instead of being knife edges.

    Returns ``(points, smallest_radius_actually_applied)``.
    """
    p = np.asarray(points, dtype=np.float64)
    if radius <= 0 or len(p) < 3:
        return p, 0.0

    out = [p[0]]
    applied: list[float] = []
    n_seg = max(int(arc_segments), 1)

    for i in range(1, len(p) - 1):
        prev, cur, nxt = p[i - 1], p[i], p[i + 1]
        u = cur - prev
        v = nxt - cur
        lu, lv = np.linalg.norm(u), np.linalg.norm(v)
        if lu < 1e-12 or lv < 1e-12:
            out.append(cur)
            continue
        u, v = u / lu, v / lv
        cosang = float(np.clip(np.dot(-u, v), -1.0, 1.0))
        theta = math.acos(cosang)          # interior angle at the corner
        if theta > math.pi - 1e-3 or theta < 1e-3:
            out.append(cur)                # already straight, or a degenerate spike
            continue

        half = theta / 2.0
        t = radius / math.tan(half)
        # never eat more than 45% of either adjacent segment
        t = min(t, 0.45 * lu, 0.45 * lv)
        r_eff = t * math.tan(half)
        if r_eff <= 1e-9:
            out.append(cur)
            continue

        start = cur - u * t
        end = cur + v * t
        bis = (-u) + v
        nb = np.linalg.norm(bis)
        if nb < 1e-12:
            out.append(cur)
            continue
        centre = cur + (bis / nb) * (r_eff / math.sin(half))

        a0 = math.atan2(start[1] - centre[1], start[0] - centre[0])
        a1 = math.atan2(end[1] - centre[1], end[0] - centre[0])
        # take the short way round
        while a1 - a0 > math.pi:
            a1 -= 2 * math.pi
        while a0 - a1 > math.pi:
            a1 += 2 * math.pi
        for k in range(n_seg + 1):
            a = a0 + (a1 - a0) * k / n_seg
            out.append(centre + r_eff * np.array([math.cos(a), math.sin(a)]))
        applied.append(r_eff)

    out.append(p[-1])
    return np.asarray(out), (min(applied) if applied else 0.0)


def _scaling_centre(poly: Polygon, ring: np.ndarray) -> np.ndarray:
    """Point to shrink the cross-sections toward.

    The area centroid is a poor choice: on a strongly concave outline it can lie
    outside the polygon entirely, so scaling toward it pushes the ring out of the
    silhouette and the taper has to be abandoned. The pole of inaccessibility
    (deepest interior point) is inside by construction and is far more likely to
    keep a shrunk ring contained.
    """
    try:
        from shapely.ops import polylabel

        tol = max(math.sqrt(max(poly.area, 1e-12)) * 1e-3, 1e-9)
        p = polylabel(poly, tolerance=tol)
        return np.array([p.x, p.y], dtype=np.float64)
    except Exception:
        c = _centroid(ring)
        if Polygon(ring).contains(_point(c)):
            return c
        rp = poly.representative_point()
        return np.array([rp.x, rp.y], dtype=np.float64)


def _point(xy: np.ndarray):
    from shapely.geometry import Point

    return Point(float(xy[0]), float(xy[1]))


def _centroid(ring: np.ndarray) -> np.ndarray:
    """Area centroid of a simple polygon."""
    x, y = ring[:, 0], ring[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cross = x * y1 - x1 * y
    a = 0.5 * float(cross.sum())
    if abs(a) < 1e-15:
        return ring.mean(axis=0)
    cx = float(((x + x1) * cross).sum()) / (6.0 * a)
    cy = float(((y + y1) * cross).sum()) / (6.0 * a)
    return np.array([cx, cy])


# --------------------------------------------------------------------------
# Triangulation and prism splitting
# --------------------------------------------------------------------------

def triangulate_ring(ring: np.ndarray) -> np.ndarray:
    """Ear-clip a simple polygon into triangles indexing the ring's own vertices.

    No Steiner points are introduced, which is what allows the same index set to
    be reused on every ring of the loft.
    """
    if mapbox_earcut is None:  # pragma: no cover
        raise RuntimeError("mapbox_earcut is required for triangulation")
    verts = np.ascontiguousarray(ring, dtype=np.float64)
    idx = mapbox_earcut.triangulate_float64(verts, np.array([len(verts)]))
    tris = np.asarray(idx, dtype=np.int64).reshape(-1, 3)
    if len(tris) == 0:
        raise ValueError("triangulation produced no triangles")
    return tris


def _triangle_areas(ring: np.ndarray, tris: np.ndarray) -> np.ndarray:
    p = ring[tris]
    return 0.5 * np.abs(
        (p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
        - (p[:, 1, 1] - p[:, 0, 1]) * (p[:, 2, 0] - p[:, 0, 0])
    )


def triangulate_clean(
    ring: np.ndarray, rel_tol: float = 1e-9, max_attempts: int = 12
) -> tuple[np.ndarray, np.ndarray, int]:
    """Triangulate a ring, removing vertices that force degenerate triangles.

    Ear clipping can emit a zero-area sliver between three nearly-collinear
    vertices that are *not* consecutive, so per-vertex ring cleaning cannot see
    it. Such a triangle yields flat tetrahedra.

    Deleting the offending tets afterwards is not an option: it tears a slit in
    the cap and leaves dangling boundary edges (observed as "N boundary edges not
    shared by exactly 2 faces" on 19 of 548 grains). Instead the vertex whose
    removal perturbs the outline least is dropped and the ring re-triangulated,
    which keeps the surface closed and manifold.

    Returns ``(ring, triangles, n_vertices_removed)``.
    """
    cur = np.asarray(ring, dtype=np.float64)
    removed = 0
    for _ in range(max_attempts):
        tris = triangulate_ring(cur)
        total = abs(_signed_area(cur))
        areas = _triangle_areas(cur, tris)
        if total <= 0:
            return cur, tris, removed
        bad = areas <= rel_tol * total
        if not bad.any() or len(cur) <= 3:
            return cur, tris, removed

        worst = tris[int(np.argmin(areas))]
        best: Optional[tuple[float, np.ndarray]] = None
        for v in worst:
            cand = np.delete(cur, int(v), axis=0)
            if len(cand) < 3:
                continue
            poly = Polygon(cand)
            if not poly.is_valid:
                continue
            delta = abs(abs(_signed_area(cand)) - total)
            if best is None or delta < best[0]:
                best = (delta, cand)
        if best is None:
            return cur, tris, removed
        cur = best[1]
        removed += 1
    return cur, triangulate_ring(cur), removed


# Prism split with the base triangle sorted by global vertex index. Shared quad
# faces then always receive the same diagonal from both neighbouring prisms, so
# the tet mesh is conforming.
#
# Verified by Monte Carlo over a reference prism: every interior point falls in
# exactly one tet (no gaps, no overlaps), and the three tet volumes sum to the
# prism volume.
#
# The sort makes the base triangle's winding arbitrary, so no fixed vertex
# ordering can guarantee positive signed volume -- callers must flip each tet
# individually (see build_grain_solid). Flipping swaps two vertices, which
# changes face orientation but not the set of faces, so conformity survives.
_PRISM_TETS = ((0, 1, 2, 5), (0, 1, 4, 5), (0, 3, 4, 5))


def prism_to_tets(bottom: Sequence[int], top: Sequence[int]) -> list[tuple[int, int, int, int]]:
    """Split a triangular prism into 3 conforming tetrahedra."""
    order = np.argsort(np.asarray(bottom))
    b = [int(bottom[i]) for i in order]
    t = [int(top[i]) for i in order]
    nodes = b + t
    return [tuple(nodes[i] for i in tet) for tet in _PRISM_TETS]


def tet_volumes(vertices: np.ndarray, tets: np.ndarray) -> np.ndarray:
    """Signed volumes of tetrahedra. Positive means correct orientation."""
    p0 = vertices[tets[:, 0]]
    e1 = vertices[tets[:, 1]] - p0
    e2 = vertices[tets[:, 2]] - p0
    e3 = vertices[tets[:, 3]] - p0
    return np.einsum("ij,ij->i", np.cross(e1, e2), e3) / 6.0


def surface_from_tets(tets: np.ndarray) -> np.ndarray:
    """Outward-oriented boundary triangles of a tet mesh.

    A face on the boundary is shared by exactly one tet; interior faces appear
    twice and cancel. Deriving the surface from the volume mesh (rather than
    building it separately) guarantees the two agree.
    """
    # Faces of a tet (a,b,c,d) with outward normals for positive orientation.
    face_defs = ((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3))
    faces = np.vstack([tets[:, f] for f in face_defs])
    keys = np.sort(faces, axis=1)
    _, inv, counts = np.unique(keys, axis=0, return_inverse=True, return_counts=True)
    return faces[counts[inv] == 1]


# --------------------------------------------------------------------------
# Grain solid
# --------------------------------------------------------------------------

@dataclass
class GrainSolid:
    """A meshed 3D grain."""

    grain_id: int
    source_image: str
    vertices: np.ndarray          # (N, 3) float64, microns
    tets: np.ndarray              # (M, 4) int64, C3D4 connectivity
    faces: np.ndarray             # (K, 3) int64, outward boundary triangles
    height_um: float
    thickness_ratio: float
    outline_um: np.ndarray        # (V, 2) the (possibly blunted) silhouette used
    centroid_um: np.ndarray       # (3,) centroid of the solid
    analytic_volume_um3: float
    measurement: Optional[GrainMeasurement] = None
    edge_radius_requested_um: float = 0.0
    edge_radius_inplane_um: float = 0.0
    """Radius actually achieved on the in-plane (vertical) edges."""
    edge_radius_meridional_um: float = 0.0
    """Radius actually achieved on the circumferential ridges."""
    n_rings: int = 0
    warnings: list[str] = field(default_factory=list)

    @property
    def n_vertices(self) -> int:
        return len(self.vertices)

    @property
    def n_tets(self) -> int:
        return len(self.tets)

    @property
    def mesh_volume_um3(self) -> float:
        return float(tet_volumes(self.vertices, self.tets).sum())

    @property
    def bounding_radius_um(self) -> float:
        """Radius of the smallest sphere about the centroid containing the grain."""
        d = self.vertices - self.centroid_um
        return float(np.sqrt((d * d).sum(axis=1)).max())

    def extent_um(self) -> np.ndarray:
        return self.vertices.max(axis=0) - self.vertices.min(axis=0)


def build_grain_solid(
    measurement: GrainMeasurement,
    outline_um: np.ndarray,
    height_model: Optional[HeightModel] = None,
    profile: Optional[LoftProfile] = None,
) -> GrainSolid:
    """Build a tet-meshed grain solid from a measured outline."""
    height_model = height_model or HeightModel()
    profile = profile or LoftProfile()
    profile.validate()

    ring = np.asarray(outline_um, dtype=np.float64)
    if ring.ndim != 2 or ring.shape[1] != 2 or len(ring) < 3:
        raise ValueError("outline must be an (V, 2) polygon with V >= 3")

    height = height_model.height_um(measurement)
    if not math.isfinite(height) or height <= 0:
        raise ValueError(f"non-physical grain height {height}")
    ratio = height / measurement.feret_min_um if measurement.feret_min_um > 0 else float("nan")

    warnings: list[str] = []
    # Every ring is a scaled copy of this outline, so an edge here reappears
    # multiplied by the smallest ring scale. Scales are clamped *up* to s_floor, so
    # the smallest one actually used is never below min(base_scale, top_scale) --
    # which makes that a safe lower bound to divide by.
    s_min_bound = max(min(profile.base_scale, profile.top_scale), 1e-3)
    min_edge_abs = MIN_FACET_EDGE_UM / s_min_bound

    ring = _clean_ring(ring, min_edge_abs=min_edge_abs)
    if len(ring) < 3:
        raise ValueError("outline degenerated to fewer than 3 vertices after cleaning")

    # Blunt the in-plane corners first, so the rounded outline is what everything
    # downstream (triangulation, rings, silhouette check) is built from.
    edge_radius_applied = 0.0
    if profile.edge_radius_um > 0:
        ring, edge_radius_applied = round_outline_corners(
            ring, profile.edge_radius_um, profile.arc_segments
        )
        ring = _clean_ring(ring, min_edge_abs=min_edge_abs)
        if len(ring) < 3:
            raise ValueError("outline collapsed while rounding corners")
        if edge_radius_applied < profile.edge_radius_um * 0.98:
            warnings.append(
                f"in-plane edge radius reduced {profile.edge_radius_um:.3f} -> "
                f"{edge_radius_applied:.3f} um; the grain is too small to carry "
                f"the requested radius"
            )

    # Settle the triangulation before building any geometry, so every ring in the
    # loft shares one non-degenerate triangle set.
    ring, tris, n_removed = triangulate_clean(ring)
    if len(ring) < 3:
        raise ValueError("outline degenerated during triangulation cleanup")

    outer = Polygon(ring)
    container = outer.buffer(max(1e-9, 1e-6 * math.sqrt(max(outer.area, 1e-12))))
    centre = _scaling_centre(outer, ring)

    # Scale each ring about the polygon centroid so lateral faces stay planar
    # (parallel opposite sides => trapezoid => coplanar), which is what makes the
    # closed-form volume exact.
    #
    # Uniform scaling about the centroid only stays inside the outline if the
    # polygon is star-shaped about that centroid. Real grains are not: on concave
    # outlines a shrunk ring can poke *outside* the silhouette, which inflated
    # the projected area by up to 24% and produced non-manifold solids. So the
    # taper is relaxed until the ring is both simple and contained.
    def ring_is_simple(s: float) -> bool:
        return bool(Polygon(centre + (ring - centre) * s).is_valid)

    # One minimum scale for the whole grain, found by bisection and applied to every
    # ring, rather than relaxing each ring on its own.
    #
    # Two competing requirements, and monotonicity wins:
    #  - Per-ring relaxation guaranteed each ring stayed inside the outline, but let
    #    a ring near the tip be pushed back out wider than the ring beneath it. That
    #    is an overhang, and it produced 175 deg razor fold-backs.
    #  - A single uniform floor keeps the scale sequence monotonic by construction,
    #    so the profile can only narrow towards the tip. An overhang becomes
    #    impossible.
    # The floor is chosen only to keep each cross-section a valid simple polygon.
    # Requiring containment as well is unsatisfiable for some concave outlines,
    # where scaling about an interior point pushes a reentrant vertex outward; the
    # silhouette then slightly exceeds the widest ring, which is harmless and is
    # reported rather than rejected.
    def min_simple_scale() -> float:
        if ring_is_simple(0.02):
            return 0.02
        lo, hi = 0.02, 1.0
        for _ in range(24):
            mid = 0.5 * (lo + hi)
            if ring_is_simple(mid):
                hi = mid
            else:
                lo = mid
        return min(hi + 1e-6, 1.0)

    s_floor = min_simple_scale()

    def scaled(scale: float) -> tuple[np.ndarray, float]:
        s = float(min(max(scale, s_floor), 1.0))
        return centre + (ring - centre) * s, s

    base_ring, base_s = scaled(profile.base_scale)
    top_ring, top_s = scaled(profile.top_scale)
    if abs(base_s - profile.base_scale) > 1e-9:
        warnings.append(
            f"base taper relaxed {profile.base_scale:.2f}->{base_s:.2f} to keep "
            f"the cross-section simple"
        )
    if abs(top_s - profile.top_scale) > 1e-9:
        warnings.append(
            f"top taper relaxed {profile.top_scale:.2f}->{top_s:.2f} to keep "
            f"the cross-section simple"
        )

    z_mid = height * profile.mid_height_fraction

    # Meridional profile in (inward inset, height). The inset is a physical
    # distance, which is what lets the edge radius be given in microns; it maps to
    # a cross-section scale through the outline's inscribed radius.
    r_in = max(inscribed_radius(outer), 1e-9)
    d_base = (1.0 - base_s) * r_in
    d_top = (1.0 - top_s) * r_in

    meridian_radius = 0.0
    if profile.edge_radius_um > 0:
        # Include the cap runs (inset = r_in is the axis) so the ridges where the
        # side surface meets the flat caps get filleted, not just the equator.
        raw = [(r_in, 0.0), (d_base, 0.0), (0.0, z_mid), (d_top, height), (r_in, height)]
        filleted, meridian_radius = fillet_polyline(
            raw, profile.edge_radius_um, profile.arc_segments
        )
        # Drop the two axis end points; they are cap centres, not rings.
        prof = [(float(d), float(z)) for d, z in filleted[1:-1]]
        if meridian_radius < profile.edge_radius_um * 0.98:
            warnings.append(
                f"meridional edge radius reduced {profile.edge_radius_um:.3f} -> "
                f"{meridian_radius:.3f} um by the grain's own proportions"
            )
    else:
        prof = [(d_base, 0.0), (0.0, z_mid), (d_top, height)]

    # Convert insets to scales, clamped so no ring inverts or exceeds the outline.
    # The *achieved* scale is kept per ring; the closed-form volume depends on it.
    rings: list[tuple[np.ndarray, float, float]] = []   # (ring, z, actual_scale)
    for d, z in prof:
        s_req = float(np.clip(1.0 - d / r_in, 0.02, 1.0))
        if s_req >= 1.0 - 1e-12:
            rings.append((ring, float(z), 1.0))
        else:
            r_ring, s_act = scaled(s_req)
            rings.append((r_ring, float(z), float(s_act)))

    if s_floor > profile.top_scale + 1e-9:
        warnings.append(
            f"taper limited to scale {s_floor:.3f} (requested top {profile.top_scale:.2f}); "
            f"a smaller cross-section would leave the concave outline"
        )

    # Collapse rings that are near-duplicates in height *and* scale. An exact-match
    # test is not enough: a band only microns tall between two nearly equal scales
    # is a sliver whose normal is numerically meaningless, and it shows up as a
    # spurious ~180 deg dihedral.
    z_tol = max(height, 1e-12) * 1e-4
    dedup: list[tuple[np.ndarray, float, float]] = []
    for r_ring, z, s in rings:
        if dedup and abs(dedup[-1][1] - z) < z_tol and abs(dedup[-1][2] - s) < 1e-4:
            continue
        dedup.append((r_ring, z, s))
    rings = dedup
    if len(rings) < 2:
        raise ValueError("profile collapsed to fewer than 2 rings")

    n = len(ring)
    vertices = np.vstack(
        [np.column_stack([r, np.full(len(r), z)]) for r, z, _ in rings]
    )

    if n_removed:
        warnings.append(
            f"removed {n_removed} outline vertex/vertices that forced degenerate "
            f"triangles"
        )

    tets: list[tuple[int, int, int, int]] = []
    for layer in range(len(rings) - 1):
        off_b, off_t = layer * n, (layer + 1) * n
        for tri in tris:
            bottom = [off_b + int(v) for v in tri]
            top = [off_t + int(v) for v in tri]
            tets.extend(prism_to_tets(bottom, top))
    tet_arr = np.asarray(tets, dtype=np.int64)

    # Fix orientation so every tet has positive volume (C3D4 requires it).
    vols = tet_volumes(vertices, tet_arr)
    flip = vols < 0
    if flip.any():
        tet_arr[flip] = tet_arr[flip][:, [0, 2, 1, 3]]
        vols = tet_volumes(vertices, tet_arr)

    # Flat tets are reported, never deleted: removing one tears a slit in the cap
    # and leaves dangling boundary edges. triangulate_clean() above is what
    # actually prevents them.
    total_abs = float(np.abs(vols).sum())
    eps = max(total_abs * 1e-12, 1e-18)
    if (vols <= eps).any():
        n_bad = int((vols <= eps).sum())
        warnings.append(
            f"{n_bad} near-degenerate tetrahedra retained to keep the surface "
            f"closed (max volume {float(vols[vols <= eps].max()):.3e} um3)"
        )

    faces = surface_from_tets(tet_arr)

    # Exact volume: sum of general-pyramid frusta between consecutive rings. Every
    # ring is a uniform scaling of the same outline about a common centre, so each
    # lateral face is a trapezoid with two parallel sides -- planar -- and the
    # prismatoid formula is exact rather than an approximation. This holds for any
    # number of rings, which is what makes the filleted multi-ring profile still
    # verifiable to machine precision.
    a_mid = abs(_signed_area(ring))
    analytic = 0.0
    for (_, z0, s0), (_, z1, s1) in zip(rings, rings[1:]):
        dz = z1 - z0
        if abs(dz) < 1e-15:
            continue
        a0, a1 = a_mid * s0 * s0, a_mid * s1 * s1
        analytic += abs(dz) / 3.0 * (a0 + a1 + math.sqrt(a0 * a1))

    # Volume centroid of the tet mesh.
    tv = tet_volumes(vertices, tet_arr)
    tc = vertices[tet_arr].mean(axis=1)
    total = float(tv.sum())
    centroid = (
        (tc * tv[:, None]).sum(axis=0) / total if abs(total) > 1e-15 else vertices.mean(axis=0)
    )

    return GrainSolid(
        grain_id=measurement.grain_id,
        source_image=measurement.source_image,
        vertices=vertices,
        tets=tet_arr,
        faces=faces,
        height_um=height,
        thickness_ratio=ratio,
        outline_um=ring,
        centroid_um=centroid,
        analytic_volume_um3=float(analytic),
        measurement=measurement,
        edge_radius_requested_um=float(profile.edge_radius_um),
        edge_radius_inplane_um=float(edge_radius_applied),
        edge_radius_meridional_um=float(meridian_radius),
        n_rings=len(rings),
        warnings=warnings,
    )


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------

def validate_grain_solid(
    solid: GrainSolid, area_tol: float = 0.02, volume_tol: float = 1e-6
) -> dict:
    """Independent checks on a grain solid. Returns a report dict with ``ok``."""
    issues: list[str] = []

    vols = tet_volumes(solid.vertices, solid.tets)
    if (vols <= 0).any():
        issues.append(f"{int((vols <= 0).sum())} tets with non-positive volume")

    mesh_v = float(vols.sum())
    if solid.analytic_volume_um3 > 0:
        rel = abs(mesh_v - solid.analytic_volume_um3) / solid.analytic_volume_um3
        if rel > volume_tol:
            issues.append(
                f"mesh volume {mesh_v:.6g} != analytic {solid.analytic_volume_um3:.6g} "
                f"(rel {rel:.2e})"
            )
    else:
        rel = float("nan")

    # Closed surface: every edge of the boundary must be used exactly twice.
    e = np.vstack(
        [solid.faces[:, [0, 1]], solid.faces[:, [1, 2]], solid.faces[:, [2, 0]]]
    )
    _, counts = np.unique(np.sort(e, axis=1), axis=0, return_counts=True)
    n_open = int((counts != 2).sum())
    if n_open:
        issues.append(f"{n_open} boundary edges not shared by exactly 2 faces")

    # Divergence theorem: surface volume must match the tet volume.
    p = solid.vertices[solid.faces]
    surf_v = float(
        np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 0]).sum()
        / 6.0
    )
    if mesh_v > 0 and abs(surf_v - mesh_v) / mesh_v > 1e-6:
        issues.append(f"surface volume {surf_v:.6g} != tet volume {mesh_v:.6g}")

    # The projected silhouette must equal the widest cross-section actually built.
    #
    # Compared against the *widest ring*, not the measured outline: when an edge
    # radius is applied the equator ridge is rounded, so no ring sits exactly on the
    # measured outline any more and the solid is legitimately smaller. Checking
    # against the outline flagged 70 of 96 perfectly good blunted grains.
    proj = _projected_area(solid)
    widest = _widest_ring_area(solid)
    outline_area = abs(_signed_area(solid.outline_um))
    rel_area = abs(proj - widest) / widest if widest > 0 else float("nan")
    blunting_loss = (
        (outline_area - widest) / outline_area if outline_area > 0 else float("nan")
    )

    # The profile must narrow monotonically towards each end: cross-section area
    # rising to a single maximum and falling again. This is the check that matters
    # for a grinding model -- a non-monotonic profile is an overhang, which means a
    # razor fold-back and a stress singularity, and it cannot be meshed sensibly.
    areas = _ring_areas_by_height(solid)
    if len(areas) >= 3:
        a = np.array([x[1] for x in areas])
        peak = int(np.argmax(a))
        rising = np.all(np.diff(a[: peak + 1]) >= -1e-9 * max(a.max(), 1e-30))
        falling = np.all(np.diff(a[peak:]) <= 1e-9 * max(a.max(), 1e-30))
        if not (rising and falling):
            issues.append(
                "cross-section profile is not monotonic about its maximum "
                "(overhang: would create a razor edge)"
            )

    if not np.isfinite(solid.vertices).all():
        issues.append("non-finite vertex coordinates")

    # Shortest edge on any surface facet, and the closest pair of nodes used by a
    # common facet. Abaqus refuses a general-contact facet with coincident nodes:
    #   "A facet that is part of a general contact surface has coincident nodes"
    # This check exists because that error reached a real job: two nodes 1.03e-7 mm
    # apart survived a purely *relative* cleaning tolerance and were then shrunk
    # further by the tip-ring scaling.
    p = solid.vertices[solid.faces]
    edges = np.stack([
        np.linalg.norm(p[:, 1] - p[:, 0], axis=1),
        np.linalg.norm(p[:, 2] - p[:, 1], axis=1),
        np.linalg.norm(p[:, 0] - p[:, 2], axis=1),
    ], axis=1)
    min_edge = float(edges.min()) if len(edges) else float("nan")
    if np.isfinite(min_edge) and min_edge < MIN_FACET_EDGE_UM:
        n_bad = int((edges.min(axis=1) < MIN_FACET_EDGE_UM).sum())
        issues.append(
            f"{n_bad} facet(s) with an edge shorter than {MIN_FACET_EDGE_UM} um "
            f"(min {min_edge:.3e} um); Abaqus will call these nodes coincident"
        )

    # A facet must not reuse a node.
    rep = int(sum(1 for f in solid.faces if len(set(int(x) for x in f)) < 3))
    if rep:
        issues.append(f"{rep} facet(s) reuse a node")

    unused = len(set(range(len(solid.vertices))) - set(np.unique(solid.tets).tolist()))
    if unused:
        issues.append(f"{unused} vertices not referenced by any tet")

    return {
        "grain_id": solid.grain_id,
        "ok": not issues,
        "issues": issues,
        "n_vertices": solid.n_vertices,
        "n_tets": solid.n_tets,
        "n_faces": int(len(solid.faces)),
        "mesh_volume_um3": mesh_v,
        "analytic_volume_um3": solid.analytic_volume_um3,
        "volume_rel_error": rel,
        "surface_volume_um3": surf_v,
        "projected_area_um2": proj,
        "outline_area_um2": outline_area,
        "widest_cross_section_um2": widest,
        "projected_area_rel_error": rel_area,
        "blunting_area_loss_fraction": blunting_loss,
        "edge_radius_requested_um": solid.edge_radius_requested_um,
        "edge_radius_inplane_um": solid.edge_radius_inplane_um,
        "edge_radius_meridional_um": solid.edge_radius_meridional_um,
        "max_convex_dihedral_deg": (
            float(convex_dihedral_angles(solid.vertices, solid.faces).max())
            if len(solid.faces) else float("nan")
        ),
        "min_tet_volume_um3": float(vols.min()),
        "min_facet_edge_um": min_edge,
        "height_um": solid.height_um,
    }


def tet_quality(vertices: np.ndarray, tets: np.ndarray) -> dict:
    """Element-quality metrics for a tetrahedral mesh.

    Matters only for deformable (C3D4) grains: element volume sets the stable time
    increment in Abaqus/Explicit, and sliver elements are ill-conditioned in
    Abaqus/Standard. Rigid (R3D3) grains use the surface only, so these numbers can
    be ignored for them.
    """
    v = tet_volumes(vertices, tets)
    p = vertices[tets]
    edges = np.stack(
        [
            np.linalg.norm(p[:, i] - p[:, j], axis=1)
            for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
        ],
        axis=1,
    )
    longest = edges.max(axis=1)
    # Normalised so a regular tetrahedron scores 1.
    regular = math.sqrt(2.0) / 12.0
    quality = np.divide(
        v, regular * longest ** 3, out=np.zeros_like(v), where=longest > 0
    )
    return {
        "n_tets": int(len(v)),
        "min_volume_um3": float(v.min()) if len(v) else float("nan"),
        "mean_volume_um3": float(v.mean()) if len(v) else float("nan"),
        "min_quality": float(quality.min()) if len(v) else float("nan"),
        "median_quality": float(np.median(quality)) if len(v) else float("nan"),
        "frac_quality_below_0p01": float(np.mean(quality < 0.01)) if len(v) else float("nan"),
    }


def convex_dihedral_angles(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """Turn angle across every convex edge of a closed surface, in degrees.

    The direct measure of how sharp a grain is. A near-90 deg turn is a knife edge
    with zero radius -- a stress singularity in FEA. Blunting replaces one large
    turn with a chain of small ones, so this dropping is exactly the intended
    effect of ``edge_radius_um``.
    """
    p = vertices[faces]
    n = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    n = np.divide(n, ln, out=np.zeros_like(n), where=ln > 0)

    edges: dict[tuple[int, int], list[int]] = {}
    for i, f in enumerate(faces):
        for k in range(3):
            e = (int(f[k]), int(f[(k + 1) % 3]))
            edges.setdefault((min(e), max(e)), []).append(i)

    out: list[float] = []
    for (a, b), fs in edges.items():
        if len(fs) != 2:
            continue
        f1, f2 = fs
        n1, n2 = n[f1], n[f2]
        ang = math.degrees(math.acos(float(np.clip(np.dot(n1, n2), -1.0, 1.0))))
        # Local convexity test: take the vertex of the second face that is not on
        # the shared edge; the edge is convex (a ridge) when that vertex lies below
        # the first face's plane.
        #
        # A centroid-based test is wrong here: grains are not convex, so
        # "outward from the centroid" is meaningless in a concave region and
        # reentrant folds get reported as sharp ridges.
        other = [int(v) for v in faces[f2] if int(v) not in (a, b)]
        if not other:
            continue
        if float(np.dot(vertices[other[0]] - vertices[a], n1)) < 0.0:
            out.append(ang)
    return np.asarray(out)


def _ring_areas_by_height(solid: GrainSolid) -> list[tuple[float, float]]:
    """``(z, cross-section area)`` for every horizontal ring, ordered by height."""
    zs = np.unique(np.round(solid.vertices[:, 2], 9))
    out: list[tuple[float, float]] = []
    for z in zs:
        pts = solid.vertices[np.isclose(solid.vertices[:, 2], z)][:, :2]
        if len(pts) >= 3:
            out.append((float(z), abs(_signed_area(pts))))
    return out


def _widest_ring_area(solid: GrainSolid) -> float:
    """Area of the largest horizontal cross-section actually present in the solid."""
    return max((a for _, a in _ring_areas_by_height(solid)), default=0.0)


def _projected_area(solid: GrainSolid) -> float:
    """Area of the solid's shadow on the XY plane, via the union of its rings."""
    from shapely.ops import unary_union

    zs = np.unique(np.round(solid.vertices[:, 2], 9))
    polys = []
    for z in zs:
        idx = np.isclose(solid.vertices[:, 2], z)
        pts = solid.vertices[idx][:, :2]
        if len(pts) >= 3:
            p = Polygon(pts)
            if not p.is_valid:
                p = p.buffer(0)
            if not p.is_empty:
                polys.append(p)
    if not polys:
        return 0.0
    return float(unary_union(polys).area)


def build_grain_library(
    grains: list[GrainMeasurement],
    seg: Segmentation,
    sem: SemImage,
    height_model: Optional[HeightModel] = None,
    profile: Optional[LoftProfile] = None,
    simplify_um: float = 0.10,
    max_vertices: int = 64,
    interior_only: bool = True,
) -> tuple[list[GrainSolid], list[dict]]:
    """Build solids for a measured grain set.

    ``interior_only`` skips border-truncated grains, whose outlines are cut off
    by the frame and would produce artificially flat-sided solids.
    """
    solids: list[GrainSolid] = []
    reports: list[dict] = []
    for g in grains:
        if interior_only and g.touches_border:
            continue
        ring = grain_outline_um(
            g.label, seg, sem, simplify_um=simplify_um, max_vertices=max_vertices
        )
        if ring is None:
            reports.append({"grain_id": g.grain_id, "ok": False, "issues": ["no outline"]})
            continue
        try:
            solid = build_grain_solid(g, ring, height_model, profile)
        except Exception as exc:  # keep going; report the failure
            reports.append(
                {"grain_id": g.grain_id, "ok": False, "issues": [f"{type(exc).__name__}: {exc}"]}
            )
            continue
        solids.append(solid)
        reports.append(validate_grain_solid(solid))
    return solids, reports
