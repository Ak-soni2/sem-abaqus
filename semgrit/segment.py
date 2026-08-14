"""Segmentation of individual abrasive grains from an SEM micrograph.

The grains in these images are bright, angular and heavily in contact, so the
work splits into two problems: separating grain from background, and splitting
grains that touch.

Design notes / departures from the original notebook:

* No histogram equalisation. It is a non-monotonic, global remap that amplifies
  noise and shifts the threshold for no benefit; Otsu already adapts.
* No blanking of a fixed image rectangle. That destroyed ~5% of real grain area.
* Watershed runs on a continuous elevation map built from the distance
  transform and the intensity gradient. The original used a *binary* Canny
  image, which is all plateaus and therefore places boundaries arbitrarily.
* Every size threshold is specified in microns, not pixels, so the same physical
  grain is treated identically at 5 kX and 10 kX.
* Grains touching the image border are kept but flagged, because including
  truncated grains in a size distribution biases it low.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import cv2
import numpy as np
from scipy import ndimage
from skimage.feature import peak_local_max
from skimage.filters import threshold_multiotsu, threshold_otsu
from skimage.morphology import h_maxima, remove_small_holes
from skimage.segmentation import watershed

from .metrology import SemImage


@dataclass
class SegmentationParams:
    """Tunable parameters, all physical units so they are magnification-free."""

    # --- denoising -----------------------------------------------------------
    median_um: float = 0.06
    """Median filter radius. Removes SE shot noise without rounding corners."""

    bilateral_um: float = 0.0
    """Optional edge-preserving smoothing radius; 0 disables."""

    # --- thresholding --------------------------------------------------------
    threshold_method: str = "multiotsu"
    """'otsu' | 'multiotsu'. Multi-Otsu splits background / dull face / bright
    face and keeps the upper two classes, which stops Otsu from slicing through
    mid-grey grain facets."""

    background_percentile: float = 0.0
    """If > 0, flatten illumination by dividing out a heavily blurred estimate."""

    # --- morphological cleanup ----------------------------------------------
    close_um: float = 0.08
    open_um: float = 0.05
    fill_holes_um2: float = 0.5

    # --- splitting touching grains ------------------------------------------
    min_grain_um: float = 0.9
    """Smallest expected grain diameter. Drives seed spacing."""

    h_maxima_um: float = 0.12
    """A distance-transform peak must exceed its surroundings by this much (in
    microns of inscribed radius) to seed a separate grain.

    Set deliberately low so the watershed *over*-segments; false splits are then
    removed by the evidence-based merge below. A single global value cannot do
    both jobs: raising it to suppress fragments also merges genuinely distinct
    grains, and solidity cannot detect the over-segmentation because splitting
    one grain in two leaves both halves convex."""

    gradient_weight: float = 1.0
    """Blend of intensity gradient into the watershed elevation. 0 = pure
    shape-based splitting, higher = boundaries snap harder to visible edges.
    Measured on this dataset: 0 gives a boundary/interior gradient ratio of
    ~1.6-1.9, while >=0.5 gives 3.5-7.5."""

    merge_weak_boundaries: bool = True
    """Remove splits whose dividing line does not lie on a real intensity edge."""

    min_edge_strength: float = 1.5
    """Minimum mean |grad| along a grain-grain boundary, as a multiple of the
    mean |grad| over all grain pixels, for the split to be kept.

    The boundary-strength distribution on this dataset is bimodal -- weak
    artefacts cluster at 1.05-1.46 and genuine edges start around 2.2 -- so this
    sits in the gap. Raising it further changes nothing (identical results at
    1.5 and 2.0 on the B4C images)."""

    min_neck_ratio: float = 0.72
    """A split is also kept, regardless of edge strength, when the grains meet
    at a narrow neck: boundary inscribed radius divided by the smaller region's
    maximum inscribed radius below this value indicates two touching particles
    rather than one grain cut in half."""

    compactness: float = 0.0
    """skimage watershed compactness; small positive values discourage leaking
    through narrow necks."""

    # --- region filtering ----------------------------------------------------
    min_area_um2: float = 0.7
    max_area_fraction: float = 0.25
    """Reject a region covering more than this fraction of the field: it is a
    merge failure, not a grain."""

    def px(self, microns: float, pixel_size_um: float) -> int:
        return int(round(microns / pixel_size_um))

    def odd_px(self, microns: float, pixel_size_um: float, minimum: int = 3) -> int:
        n = self.px(microns, pixel_size_um)
        n = max(n, minimum)
        return n if n % 2 == 1 else n + 1


@dataclass
class Segmentation:
    """Result of segmenting one micrograph."""

    labels: np.ndarray            # int32, 0 = background
    foreground: np.ndarray        # bool
    distance_um: np.ndarray       # float32 distance transform in microns
    n_seeds: int
    pixel_size_um: float
    params: SegmentationParams
    border_labels: set[int] = field(default_factory=set)
    rejected: dict[str, int] = field(default_factory=dict)
    threshold_values: tuple[float, ...] = ()

    @property
    def label_ids(self) -> list[int]:
        ids = np.unique(self.labels)
        return [int(i) for i in ids if i != 0]

    @property
    def n_grains(self) -> int:
        return len(self.label_ids)

    def interior_label_ids(self) -> list[int]:
        return [i for i in self.label_ids if i not in self.border_labels]


def _flatten_background(img: np.ndarray, sigma_px: float) -> np.ndarray:
    """Divide out a large-scale illumination estimate, preserving mean level."""
    if sigma_px <= 0:
        return img
    bg = cv2.GaussianBlur(img.astype(np.float32), (0, 0), sigma_px)
    bg = np.maximum(bg, 1e-3)
    flat = img.astype(np.float32) / bg
    flat *= float(bg.mean())
    return np.clip(flat, 0, 255).astype(np.uint8)


def _threshold(img: np.ndarray, params: SegmentationParams) -> tuple[np.ndarray, tuple[float, ...]]:
    """Separate grain from background."""
    if params.threshold_method == "multiotsu":
        try:
            thresholds = threshold_multiotsu(img, classes=3)
            # Keep everything above the *lower* threshold: dull grain facets and
            # bright grain facets are both grain.
            return img > thresholds[0], tuple(float(t) for t in thresholds)
        except ValueError:
            pass  # degenerate histogram, fall through to Otsu
    t = threshold_otsu(img)
    return img > t, (float(t),)


def _gradient_magnitude(img: np.ndarray) -> np.ndarray:
    """Sobel gradient magnitude of a lightly smoothed image."""
    sm = cv2.GaussianBlur(img, (3, 3), 0).astype(np.float32)
    gx = cv2.Sobel(sm, cv2.CV_32F, 1, 0, ksize=3)
    gy = cv2.Sobel(sm, cv2.CV_32F, 0, 1, ksize=3)
    return np.hypot(gx, gy)


class _UnionFind:
    def __init__(self, n: int) -> None:
        self.parent = list(range(n))

    def find(self, a: int) -> int:
        while self.parent[a] != a:
            self.parent[a] = self.parent[self.parent[a]]
            a = self.parent[a]
        return a

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


def boundary_evidence(
    labels: np.ndarray, gradient: np.ndarray, distance_um: np.ndarray
) -> dict[tuple[int, int], dict[str, float]]:
    """Measure, for every pair of adjacent regions, how well-supported their
    shared boundary is by the image.

    Returns a mapping ``(low_label, high_label) -> {edge_strength, neck_ratio,
    n_px}``. ``edge_strength`` is the mean gradient along the boundary divided by
    the mean gradient over all labelled pixels; ``neck_ratio`` is the largest
    inscribed radius on the boundary divided by the smaller region's maximum
    inscribed radius.
    """
    fg = labels > 0
    if not fg.any():
        return {}
    ref = float(gradient[fg].mean())
    if ref <= 0:
        return {}

    max_dist: dict[int, float] = {}
    for lid in np.unique(labels):
        if lid == 0:
            continue
        max_dist[int(lid)] = float(distance_um[labels == lid].max())

    acc: dict[tuple[int, int], list[float]] = {}
    for dy, dx in ((0, 1), (1, 0), (1, 1), (1, -1)):
        a = labels[
            max(0, -dy) : labels.shape[0] - max(0, dy),
            max(0, -dx) : labels.shape[1] - max(0, dx),
        ]
        b = labels[
            max(0, dy) : labels.shape[0] + min(0, dy),
            max(0, dx) : labels.shape[1] + min(0, dx),
        ]
        ga = gradient[
            max(0, -dy) : gradient.shape[0] - max(0, dy),
            max(0, -dx) : gradient.shape[1] - max(0, dx),
        ]
        da = distance_um[
            max(0, -dy) : distance_um.shape[0] - max(0, dy),
            max(0, -dx) : distance_um.shape[1] - max(0, dx),
        ]
        sel = (a != b) & (a > 0) & (b > 0)
        if not sel.any():
            continue
        lo = np.minimum(a[sel], b[sel])
        hi = np.maximum(a[sel], b[sel])
        gv = ga[sel]
        dv = da[sel]
        for l_, h_, g_, d_ in zip(lo.tolist(), hi.tolist(), gv.tolist(), dv.tolist()):
            key = (l_, h_)
            e = acc.setdefault(key, [0.0, 0.0, 0.0])
            e[0] += g_
            e[1] += 1.0
            e[2] = max(e[2], d_)

    out: dict[tuple[int, int], dict[str, float]] = {}
    for (l_, h_), (gsum, n, dmax) in acc.items():
        if n <= 0:
            continue
        smaller = min(max_dist.get(l_, 0.0), max_dist.get(h_, 0.0))
        out[(l_, h_)] = {
            "edge_strength": (gsum / n) / ref,
            "neck_ratio": (dmax / smaller) if smaller > 0 else 1.0,
            "n_px": n,
        }
    return out


def _merge_unsupported_splits(
    labels: np.ndarray,
    gradient: np.ndarray,
    distance_um: np.ndarray,
    params: SegmentationParams,
) -> tuple[np.ndarray, int]:
    """Undo watershed splits that the image does not support.

    A split survives if the dividing line lies on a real intensity edge, or if
    the two regions meet at a narrow neck (the signature of two particles in
    contact). Everything else is a watershed artefact and is merged back.

    Decisions are per original boundary and applied with union-find, so the
    result is independent of processing order.
    """
    evidence = boundary_evidence(labels, gradient, distance_um)
    if not evidence:
        return labels, 0

    n_max = int(labels.max())
    uf = _UnionFind(n_max + 1)
    merged = 0
    for (a, b), info in sorted(evidence.items()):
        keep = (
            info["edge_strength"] >= params.min_edge_strength
            or info["neck_ratio"] <= params.min_neck_ratio
        )
        if not keep:
            uf.union(a, b)
            merged += 1
    if merged == 0:
        return labels, 0

    remap = np.zeros(n_max + 1, dtype=np.int32)
    next_id = 1
    for lid in range(1, n_max + 1):
        root = uf.find(lid)
        if remap[root] == 0:
            remap[root] = next_id
            next_id += 1
        remap[lid] = remap[root]
    remap[0] = 0
    return remap[labels], merged


def segment_grains(
    sem: SemImage, params: Optional[SegmentationParams] = None
) -> Segmentation:
    """Segment individual grains from a calibrated SEM image."""
    params = params or SegmentationParams()
    ps = sem.pixel_size_um
    img = sem.intensity

    # ---- 1. denoise ------------------------------------------------------
    work = img
    k = params.odd_px(params.median_um, ps)
    if k >= 3:
        work = cv2.medianBlur(work, min(k, 5) if k <= 5 else 5)
    if params.bilateral_um > 0:
        d = params.px(params.bilateral_um, ps)
        if d >= 1:
            work = cv2.bilateralFilter(work, d=2 * d + 1, sigmaColor=25, sigmaSpace=d)

    # ---- 2. optional illumination flattening -----------------------------
    if params.background_percentile > 0:
        work = _flatten_background(work, sigma_px=params.background_percentile / ps)

    # ---- 3. threshold ----------------------------------------------------
    fg, thresholds = _threshold(work, params)

    # ---- 4. morphological cleanup ---------------------------------------
    def ellipse(microns: float) -> Optional[np.ndarray]:
        n = params.px(microns, ps)
        if n < 1:
            return None
        return cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * n + 1, 2 * n + 1))

    fg_u8 = fg.astype(np.uint8)
    se_close = ellipse(params.close_um)
    if se_close is not None:
        fg_u8 = cv2.morphologyEx(fg_u8, cv2.MORPH_CLOSE, se_close)
    se_open = ellipse(params.open_um)
    if se_open is not None:
        fg_u8 = cv2.morphologyEx(fg_u8, cv2.MORPH_OPEN, se_open)
    fg = fg_u8.astype(bool)

    hole_px = max(int(round(params.fill_holes_um2 / (ps * ps))), 1)
    fg = remove_small_holes(fg, hole_px)

    # ---- 5. distance transform ------------------------------------------
    dist_px = cv2.distanceTransform(fg.astype(np.uint8), cv2.DIST_L2, 5)
    dist_um = (dist_px * ps).astype(np.float32)

    # ---- 6. seeds --------------------------------------------------------
    # h-maxima on the distance map: a peak must rise h above its surroundings
    # before it counts as a separate grain. Unlike a bare local-maximum search
    # this merges the ragged plateau of peaks inside one angular grain.
    h_um = max(params.h_maxima_um, 1e-6)
    seeds_mask = h_maxima(dist_um, h=h_um).astype(bool)
    seeds_mask &= fg

    if not seeds_mask.any():
        min_dist_px = max(params.px(params.min_grain_um / 2.0, ps), 1)
        coords = peak_local_max(
            dist_um, min_distance=min_dist_px, labels=fg, exclude_border=False
        )
        seeds_mask = np.zeros(fg.shape, bool)
        if coords.size:
            seeds_mask[tuple(coords.T)] = True

    seeds, n_seeds = ndimage.label(seeds_mask)

    # ---- 7. watershed on a continuous elevation map ----------------------
    # -distance carves basins at grain centres and ridges at the necks between
    # touching grains; the gradient term pulls boundaries onto visible edges.
    grad = _gradient_magnitude(work)
    d_norm = dist_um / max(float(dist_um.max()), 1e-9)
    elevation = -d_norm
    if params.gradient_weight > 0:
        g_norm = grad / max(float(grad.max()), 1e-9)
        elevation = elevation + params.gradient_weight * g_norm

    labels = watershed(
        elevation,
        markers=seeds,
        mask=fg,
        compactness=params.compactness,
        watershed_line=False,
    ).astype(np.int32)

    # ---- 8. drop splits the image does not support -----------------------
    n_merged = 0
    if params.merge_weak_boundaries:
        labels, n_merged = _merge_unsupported_splits(labels, grad, dist_um, params)

    # ---- 9. filter regions ----------------------------------------------
    min_area_px = max(int(round(params.min_area_um2 / (ps * ps))), 1)
    max_area_px = int(params.max_area_fraction * fg.size)
    rejected = {"too_small": 0, "too_large": 0, "merged_splits": n_merged}

    out = np.zeros_like(labels)
    next_id = 1
    for lid in range(1, int(labels.max()) + 1):
        m = labels == lid
        area = int(m.sum())
        if area == 0:
            continue
        if area < min_area_px:
            rejected["too_small"] += 1
            continue
        if area > max_area_px:
            rejected["too_large"] += 1
            continue
        out[m] = next_id
        next_id += 1

    border = set()
    if out.size:
        for edge in (out[0, :], out[-1, :], out[:, 0], out[:, -1]):
            border.update(int(v) for v in np.unique(edge) if v != 0)

    return Segmentation(
        labels=out,
        foreground=fg,
        distance_um=dist_um,
        n_seeds=int(n_seeds),
        pixel_size_um=ps,
        params=params,
        border_labels=border,
        rejected=rejected,
        threshold_values=thresholds,
    )
