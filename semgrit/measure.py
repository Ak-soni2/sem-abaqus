"""Per-grain shape measurement in physical units.

Everything is measured on the *actual* segmented outline. The original pipeline
measured a convex hull simplified to 8-10 vertices, which inflates area, erases
the concave fracture surfaces, and turns the reported "facet edges" into
artefacts of the polygon simplification rather than features of the grain.

Grains truncated by the image border are measured but flagged. They must be
excluded from size distributions: a truncated grain looks smaller than it is, so
including them biases the distribution low (38-67% of regions touch the border
in these images).
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict, field
from typing import Iterable, Optional

import cv2
import numpy as np
from skimage.measure import perimeter_crofton

from .metrology import SemImage
from .segment import Segmentation


@dataclass
class GrainMeasurement:
    """Shape descriptors for a single grain, all lengths in microns."""

    grain_id: int
    label: int
    source_image: str

    # --- position ---
    centroid_x_um: float
    centroid_y_um: float

    # --- size ---
    area_um2: float
    perimeter_um: float
    equivalent_diameter_um: float
    """Diameter of the circle with the same area -- the standard particle size."""
    feret_max_um: float
    """Maximum caliper diameter (longest dimension)."""
    feret_min_um: float
    """Minimum caliper diameter (the width that governs sieve passage)."""
    major_axis_um: float
    minor_axis_um: float
    inscribed_radius_um: float

    # --- shape ---
    aspect_ratio: float
    """feret_max / feret_min. 1.0 = equiaxed."""
    elongation: float
    circularity: float
    """4*pi*A/P^2. 1.0 = perfect circle, lower = more angular/irregular."""
    solidity: float
    """area / convex-hull area. Quantifies exactly the concavity a convex hull
    would have thrown away."""
    convexity: float
    """hull perimeter / actual perimeter. Surface roughness measure."""
    orientation_deg: float

    # --- cutting-edge geometry (relevant to grinding) ---
    n_corners: int
    min_corner_angle_deg: float
    """Sharpest corner on the outline: the most aggressive cutting point."""
    mean_corner_angle_deg: float

    # --- provenance / quality ---
    touches_border: bool
    pixel_area: int
    outline_px: int = 0
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def _largest_contour(mask: np.ndarray) -> Optional[np.ndarray]:
    """Outline of a label mask at full resolution (no polygon simplification)."""
    cs, _ = cv2.findContours(
        mask.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE
    )
    if not cs:
        return None
    return max(cs, key=cv2.contourArea)


def _feret_diameters(contour: np.ndarray) -> tuple[float, float]:
    """Max and min caliper diameters, in pixels, by rotating calipers.

    ``cv2.minAreaRect`` gives the minimum-area rectangle, whose short side is a
    good min-Feret estimate, but the max-Feret is the true maximum point-to-point
    distance over the convex hull, which the rectangle can underestimate.
    """
    hull = cv2.convexHull(contour).reshape(-1, 2).astype(np.float64)
    if len(hull) < 2:
        return 0.0, 0.0

    # max Feret: exact maximum pairwise distance on the hull (hulls are small)
    d2 = ((hull[:, None, :] - hull[None, :, :]) ** 2).sum(-1)
    feret_max = float(np.sqrt(d2.max()))

    # min Feret: minimum width over all hull edge directions
    best = math.inf
    n = len(hull)
    for i in range(n):
        p, q = hull[i], hull[(i + 1) % n]
        e = q - p
        L = math.hypot(e[0], e[1])
        if L < 1e-12:
            continue
        nx, ny = -e[1] / L, e[0] / L
        proj = hull[:, 0] * nx + hull[:, 1] * ny
        best = min(best, float(proj.max() - proj.min()))
    feret_min = 0.0 if not math.isfinite(best) else best
    return feret_max, feret_min


def _corner_angles(
    contour: np.ndarray, pixel_size_um: float, min_edge_um: float = 0.25
) -> tuple[int, list[float]]:
    """Interior angles at the dominant corners of the *real* outline.

    The outline is simplified with a tolerance expressed in microns (not as a
    fraction of perimeter, which made the result depend on grain size), and the
    convex hull is deliberately not used, so concave notches survive.
    """
    tol_px = max(min_edge_um / pixel_size_um, 1.0)
    approx = cv2.approxPolyDP(contour, tol_px, True).reshape(-1, 2).astype(np.float64)
    n = len(approx)
    if n < 3:
        return n, []
    angles: list[float] = []
    for i in range(n):
        prev = approx[(i - 1) % n]
        cur = approx[i]
        nxt = approx[(i + 1) % n]
        v1, v2 = prev - cur, nxt - cur
        n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
        if n1 < 1e-9 or n2 < 1e-9:
            continue
        cosang = float(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))
        angles.append(math.degrees(math.acos(cosang)))
    return n, angles


def measure_grain(
    label: int,
    grain_id: int,
    seg: Segmentation,
    sem: SemImage,
) -> Optional[GrainMeasurement]:
    """Measure one segmented grain."""
    ps = sem.pixel_size_um
    mask = seg.labels == label
    pixel_area = int(mask.sum())
    if pixel_area == 0:
        return None

    contour = _largest_contour(mask)
    if contour is None or len(contour) < 3:
        return None

    warnings: list[str] = []

    # Area from the pixel count, not the polygon: it is exact for the
    # segmentation and immune to contour-tracing bias.
    area_um2 = pixel_area * ps * ps

    # Crofton perimeter, not the traced contour length. Tracing a digital
    # boundary follows the pixel staircase and overstates the perimeter.
    #
    # Benchmarked on squares rotated 0-45 deg against the analytic perimeter:
    #   Crofton (4 dir)  mean -0.74%, rms 2.72%, worst 5.14%   <- used
    #   Crofton (2 dir)  mean -0.73%, rms 10.68%, worst 21.07%
    #   traced contour   mean +4.34%, rms 5.15%, worst 7.40%
    # and +0.28% on a digitised circle. The residual is orientation-dependent,
    # so circularity carries roughly +/-10% uncertainty on a single grain.
    perimeter_um = float(perimeter_crofton(mask, directions=4)) * ps

    eq_diam = 2.0 * math.sqrt(area_um2 / math.pi)

    fmax_px, fmin_px = _feret_diameters(contour)
    feret_max = fmax_px * ps
    feret_min = fmin_px * ps

    # Second-moment axes.
    ys, xs = np.nonzero(mask)
    cx, cy = xs.mean(), ys.mean()
    x0, y0 = xs - cx, ys - cy
    cov = np.array(
        [[(x0 * x0).mean(), (x0 * y0).mean()], [(x0 * y0).mean(), (y0 * y0).mean()]]
    )
    evals, evecs = np.linalg.eigh(cov)
    evals = np.clip(evals, 0.0, None)
    minor_axis = 4.0 * math.sqrt(evals[0]) * ps
    major_axis = 4.0 * math.sqrt(evals[1]) * ps
    orientation = math.degrees(math.atan2(evecs[1, 1], evecs[0, 1]))

    hull = cv2.convexHull(contour)
    hull_area_um2 = float(cv2.contourArea(hull)) * ps * ps
    hull_perim_um = float(cv2.arcLength(hull, True)) * ps

    # Solidity compares two *polygon* areas. Mixing the pixel-count area with a
    # polygon hull area yields values above 1 for small or thin grains, because
    # a contour traced through pixel centres encloses less than the pixels it
    # covers (seen on 2 of 116 grains before this was corrected).
    poly_area_um2 = float(cv2.contourArea(contour)) * ps * ps
    solidity = poly_area_um2 / hull_area_um2 if hull_area_um2 > 0 else float("nan")

    # Convexity likewise compares two traced lengths.
    traced_perim_um = float(cv2.arcLength(contour, True)) * ps
    convexity = hull_perim_um / traced_perim_um if traced_perim_um > 0 else float("nan")

    circularity = (
        4.0 * math.pi * area_um2 / (perimeter_um * perimeter_um)
        if perimeter_um > 0
        else float("nan")
    )

    inscribed = float(seg.distance_um[mask].max())

    n_corners, angles = _corner_angles(contour, ps)
    min_angle = float(min(angles)) if angles else float("nan")
    mean_angle = float(np.mean(angles)) if angles else float("nan")

    aspect = feret_max / feret_min if feret_min > 1e-9 else float("nan")
    elongation = 1.0 - (feret_min / feret_max) if feret_max > 1e-9 else float("nan")

    if solidity > 1.02:
        warnings.append(f"solidity {solidity:.3f} > 1; contour/area mismatch")
    if not math.isfinite(aspect):
        warnings.append("degenerate min-Feret")

    return GrainMeasurement(
        grain_id=grain_id,
        label=int(label),
        source_image=sem.path,
        centroid_x_um=float(cx * ps),
        centroid_y_um=float(cy * ps),
        area_um2=float(area_um2),
        perimeter_um=float(perimeter_um),
        equivalent_diameter_um=float(eq_diam),
        feret_max_um=float(feret_max),
        feret_min_um=float(feret_min),
        major_axis_um=float(major_axis),
        minor_axis_um=float(minor_axis),
        inscribed_radius_um=inscribed,
        aspect_ratio=float(aspect),
        elongation=float(elongation),
        circularity=float(circularity),
        solidity=float(solidity),
        convexity=float(convexity),
        orientation_deg=float(orientation),
        n_corners=int(n_corners),
        min_corner_angle_deg=min_angle,
        mean_corner_angle_deg=mean_angle,
        touches_border=bool(label in seg.border_labels),
        pixel_area=pixel_area,
        outline_px=int(len(contour)),
        warnings=warnings,
    )


def measure_all(seg: Segmentation, sem: SemImage) -> list[GrainMeasurement]:
    """Measure every segmented grain, ordered by descending area."""
    out: list[GrainMeasurement] = []
    for i, label in enumerate(seg.label_ids, start=1):
        m = measure_grain(label, i, seg, sem)
        if m is not None:
            out.append(m)
    out.sort(key=lambda g: -g.area_um2)
    for i, g in enumerate(out, start=1):
        g.grain_id = i
    return out


# --------------------------------------------------------------------------
# Population statistics
# --------------------------------------------------------------------------

def _percentiles(v: np.ndarray, qs: Iterable[float]) -> dict[str, float]:
    if v.size == 0:
        return {f"d{int(q)}": float("nan") for q in qs}
    return {f"d{int(q)}": float(np.percentile(v, q)) for q in qs}


def grain_statistics(
    grains: list[GrainMeasurement],
    sem: SemImage,
    interior_only: bool = True,
) -> dict:
    """Population statistics for a grain set.

    ``interior_only`` excludes border-truncated grains, which is required for an
    unbiased size distribution. The areal number density is computed by the
    same rule so it stays consistent with the size statistics.
    """
    used = [g for g in grains if not g.touches_border] if interior_only else list(grains)
    eq = np.array([g.equivalent_diameter_um for g in used], float)
    fmax = np.array([g.feret_max_um for g in used], float)
    fmin = np.array([g.feret_min_um for g in used], float)
    ar = np.array([g.aspect_ratio for g in used], float)
    sol = np.array([g.solidity for g in used], float)
    circ = np.array([g.circularity for g in used], float)
    area = np.array([g.area_um2 for g in used], float)

    def stat(v: np.ndarray) -> dict:
        v = v[np.isfinite(v)]
        if v.size == 0:
            return {"n": 0}
        return {
            "n": int(v.size),
            "mean": float(v.mean()),
            "std": float(v.std(ddof=1)) if v.size > 1 else 0.0,
            "min": float(v.min()),
            "max": float(v.max()),
            **_percentiles(v, (10, 50, 90)),
        }

    total_area = float(area.sum())
    return {
        "source_image": sem.path,
        "pixel_size_um": sem.pixel_size_um,
        "pixel_size_source": sem.pixel_size_source,
        "field_width_um": sem.width_um,
        "field_height_um": sem.height_um,
        "field_area_um2": sem.field_area_um2,
        "n_grains_total": len(grains),
        "n_grains_border": sum(1 for g in grains if g.touches_border),
        "n_grains_used": len(used),
        "interior_only": interior_only,
        "areal_density_per_mm2": (
            len(used) / (sem.field_area_um2 / 1e6) if sem.field_area_um2 > 0 else float("nan")
        ),
        "area_coverage_fraction": (
            total_area / sem.field_area_um2 if sem.field_area_um2 > 0 else float("nan")
        ),
        "equivalent_diameter_um": stat(eq),
        "feret_max_um": stat(fmax),
        "feret_min_um": stat(fmin),
        "aspect_ratio": stat(ar),
        "solidity": stat(sol),
        "circularity": stat(circ),
    }


CSV_COLUMNS = [
    "grain_id", "source_image", "label", "touches_border",
    "centroid_x_um", "centroid_y_um",
    "area_um2", "perimeter_um", "equivalent_diameter_um",
    "feret_max_um", "feret_min_um", "major_axis_um", "minor_axis_um",
    "inscribed_radius_um",
    "aspect_ratio", "elongation", "circularity", "solidity", "convexity",
    "orientation_deg", "n_corners", "min_corner_angle_deg", "mean_corner_angle_deg",
    "pixel_area", "outline_px",
]


def write_csv(grains: list[GrainMeasurement], path: str) -> None:
    import csv

    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        w.writeheader()
        for g in grains:
            w.writerow(g.to_dict())


def write_json(payload: dict, path: str) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, default=str)
