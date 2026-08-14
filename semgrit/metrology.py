"""Reading SEM images and establishing the pixel->micron calibration.

The single most important job in this pipeline. Everything downstream (grain
sizes, CAD dimensions, wheel grain counts) is a multiple of ``pixel_size_um``,
so a silent error here scales the whole model.

Primary source is the instrument's own metadata. Zeiss SmartSEM writes a large
key/value block into TIFF tag 34118 which contains ``AP_IMAGE_PIXEL_SIZE`` --
the exact calibrated pixel size. The burnt-in scale bar is measured
independently and used only to *cross-check* that number.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from PIL import Image
from scipy import ndimage

# Zeiss SmartSEM private TIFF tags: 34118 is the ASCII block, 34119 the UTF-16
# copy of the same data.
ZEISS_TAG_ASCII = 34118
ZEISS_TAG_UTF16 = 34119
TIFF_COLORMAP_TAG = 320

# Scale-bar labels only ever take 1-2-5 x decade values, which lets us recover
# the label without OCR by snapping the measured bar length.
_NICE_VALUES_UM = np.array(
    [10 ** e * m for e in range(-4, 5) for m in (1.0, 2.0, 5.0)], dtype=float
)

# A row with more distinct grey levels than this is real scan data, not an
# overlay panel. Sits ~6x clear of both populations on this dataset.
_DATABAR_MAX_LEVELS = 48

_LENGTH_UNITS_TO_UM = {
    "pm": 1e-6,
    "nm": 1e-3,
    "um": 1.0,
    "mm": 1e3,
    "cm": 1e4,
    "m": 1e6,
}


class MetrologyError(RuntimeError):
    """Raised when the pixel size cannot be established or fails validation."""


def _normalise_unit(unit: str) -> Optional[str]:
    """Map a unit string onto a key of ``_LENGTH_UNITS_TO_UM``.

    Zeiss writes the micron sign as a non-ASCII byte which PIL surfaces as
    latin-1 'µ' or a replacement char, so any single non-ASCII character
    followed by 'm' is treated as micrometres.
    """
    u = unit.strip()
    if not u:
        return None
    if u in _LENGTH_UNITS_TO_UM:
        return u
    # strip a leading non-ascii micro sign of any encoding
    if len(u) == 2 and u.endswith("m") and (ord(u[0]) > 127 or u[0] in "µμ"):
        return "um"
    lowered = u.lower()
    if lowered in _LENGTH_UNITS_TO_UM:
        return lowered
    return None


def parse_length_um(text: str) -> Optional[float]:
    """Extract a length in microns from e.g. ``'Image Pixel Size = 29.30 nm'``."""
    if text is None:
        return None
    m = re.search(
        r"=\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*(\S*)", text
    )
    if not m:
        return None
    value = float(m.group(1))
    unit = _normalise_unit(m.group(2))
    if unit is None:
        return None
    return value * _LENGTH_UNITS_TO_UM[unit]


def parse_angle_deg(text: str) -> Optional[float]:
    """Extract an angle in degrees from e.g. ``'Stage at T =   0.0 deg'``."""
    if text is None:
        return None
    m = re.search(r"=\s*([-+]?\d+(?:\.\d+)?)", text)
    return float(m.group(1)) if m else None


def parse_zeiss_metadata(img: Image.Image) -> dict[str, str]:
    """Return the Zeiss ``AP_*``/``DP_*``/``SV_*`` key -> description mapping.

    The block is a flat list of lines where a key line is followed by its
    human-readable value line.
    """
    tags = getattr(img, "tag_v2", None)
    if tags is None:
        return {}

    blob = None
    for tag in (ZEISS_TAG_ASCII, ZEISS_TAG_UTF16):
        if tag not in tags:
            continue
        raw = tags[tag]
        if isinstance(raw, bytes):
            for enc in ("utf-16-le", "latin-1"):
                try:
                    raw = raw.decode(enc)
                    break
                except UnicodeDecodeError:
                    continue
        if isinstance(raw, str):
            # the UTF-16 copy read as latin-1 is riddled with NULs
            raw = raw.replace("\x00", "")
            if "AP_" in raw or "DP_" in raw:
                blob = raw
                break
    if blob is None:
        return {}

    lines = [ln.strip() for ln in re.split(r"\r\n|\r|\n", blob)]
    out: dict[str, str] = {}
    for i, line in enumerate(lines[:-1]):
        if re.fullmatch(r"(?:AP|DP|SV)_[A-Z0-9_]+", line):
            out.setdefault(line, lines[i + 1])
    return out


def _read_intensity(img: Image.Image) -> np.ndarray:
    """Return the true 8-bit intensity plane.

    For palette ("P") SEM TIFFs the palette index *is* the grey level: the
    colormap is an ``i -> i*257`` ramp except for 27 reserved annotation colours
    (every 9th index) which SmartSEM uses to draw the databar. Those reserved
    entries never occur in the micrograph, so taking raw indices both recovers
    the exact grey level and avoids the colour contamination that
    ``cv2.imread`` + ``cvtColor`` would introduce in the databar.
    """
    if img.mode == "P":
        raw = np.asarray(img)
        cmap = img.tag_v2.get(TIFF_COLORMAP_TAG) if hasattr(img, "tag_v2") else None
        if cmap is not None and len(cmap) >= 768:
            cm = np.asarray(cmap, dtype=np.int64)
            ideal = np.arange(256, dtype=np.int64) * 257
            reserved = np.where(cm[:256] != ideal)[0]
            used = np.unique(raw)
            clash = np.intersect1d(used, reserved)
            if clash.size:
                # Only meaningful if it happens inside the micrograph; the
                # caller crops the databar away, so just record it.
                pass
        return raw.astype(np.uint8)

    if img.mode in ("I;16", "I;16B", "I", "F"):
        a = np.asarray(img).astype(np.float64)
        lo, hi = float(a.min()), float(a.max())
        if hi <= lo:
            return np.zeros(a.shape, np.uint8)
        return np.round((a - lo) / (hi - lo) * 255.0).astype(np.uint8)

    if img.mode != "L":
        img = img.convert("L")
    return np.asarray(img).astype(np.uint8)


def find_databar(
    intensity: np.ndarray, search_frac: float = 0.70, min_rows: int = 8
) -> tuple[int, int]:
    """Locate the burnt-in databar as a ``(top, bottom)`` row range.

    The Zeiss databar is a near-white panel spanning the full width, fringed by
    a solid black frame. It is found as the *longest contiguous run* of
    databar-like rows in the lower part of the image.

    Deliberately not anchored to the last row: SmartSEM can leave an orphan row
    of scan data *below* the panel (row 767 of these 768-row images), which
    defeats a bottom-anchored scan. Returns ``(h, h)`` when no databar is found.
    """
    h, w = intensity.shape
    if h == 0 or w == 0:
        return h, h

    # Synthetic overlay graphics are drawn from a tiny palette, so a databar row
    # holds only a handful of distinct grey levels, while any real SEM scan row
    # holds many. Measured on this dataset: databar rows 2-14 distinct levels,
    # micrograph rows 91-225 -- a ~6x margin either side of the threshold.
    #
    # This replaces a brightness-fraction test, which silently ate 12 rows of
    # real data on images whose bottom edge happens to be bright.
    srt = np.sort(intensity, axis=1)
    n_levels = (np.diff(srt, axis=1) != 0).sum(axis=1) + 1
    databar_like = n_levels <= _DATABAR_MAX_LEVELS

    lo = int(search_frac * h)
    best: tuple[int, int] = (h, h)
    best_len = 0
    r = lo
    while r < h:
        if not databar_like[r]:
            r += 1
            continue
        start = r
        while r < h and databar_like[r]:
            r += 1
        if r - start > best_len:
            best_len, best = r - start, (start, r)
    if best_len < min_rows:
        return h, h
    return best


def find_databar_top(intensity: np.ndarray, search_frac: float = 0.70) -> int:
    """Row index where the burnt-in databar begins, or image height if absent."""
    return find_databar(intensity, search_frac)[0]


@dataclass
class ScaleBar:
    """A scale bar measured off the burnt-in databar."""

    length_px: float
    bbox: tuple[int, int, int, int]  # y0, x0, y1, x1 in full-image coords
    implied_um: Optional[float] = None       # length_px * reference pixel size
    snapped_label_um: Optional[float] = None  # nearest 1-2-5 value
    pixel_size_um: Optional[float] = None    # snapped_label_um / length_px


def measure_scale_bar(
    intensity: np.ndarray,
    databar_top: int,
    databar_bottom: Optional[int] = None,
    left_frac: float = 0.34,
) -> Optional[ScaleBar]:
    """Measure the databar scale bar, tick-centre to tick-centre.

    The bar is drawn as an I-beam: a thin horizontal rule capped by two taller
    vertical ticks. The calibrated distance is between the *tick centres*, not
    the bounding-box width -- using the bbox overestimates by the tick
    thickness (2 px of 68 px = 2.5% here).
    """
    h, w = intensity.shape
    bottom = h if databar_bottom is None else min(databar_bottom, h)
    if databar_top >= bottom:
        return None
    panel = intensity[databar_top:bottom, : max(int(left_frac * w), 8)]
    if panel.size == 0:
        return None

    white = float(intensity.max())
    ink = panel < 0.90 * white
    if not ink.any():
        return None

    # Drop the databar frame, which otherwise connects every glyph into one
    # component: any component spanning nearly the full analysed width or the
    # full panel height is frame, not content.
    labels, n = ndimage.label(ink, structure=np.ones((3, 3), bool))
    ph, pw = panel.shape
    best: Optional[ScaleBar] = None
    for idx, sl in enumerate(ndimage.find_objects(labels), start=1):
        ys, xs = sl
        bh, bw = ys.stop - ys.start, xs.stop - xs.start
        if bw >= 0.97 * pw or bh >= 0.97 * ph:
            continue  # frame
        if bw < 15 or bh < 2:
            continue
        comp = labels[sl] == idx

        # The rule must span essentially the whole component width. This is what
        # separates the bar (fill 1.00) from glyph decoys such as the 'm' of
        # "2 um", whose horizontal stroke only reaches ~0.88.
        row_fill = comp.mean(axis=1)
        if row_fill.max() < 0.95:
            continue
        rule_thickness = int((row_fill >= 0.95).sum())
        if rule_thickness < 1 or rule_thickness > max(4, bh // 2):
            continue

        # Ticks are columns markedly taller than the rule. A scale bar has
        # exactly two of them; the three stems of an 'm' give three groups.
        col_height = comp.sum(axis=0)
        tick = col_height >= max(rule_thickness * 3, int(0.5 * bh))
        groups = _contiguous_groups(np.where(tick)[0])

        if len(groups) == 2:
            length_px = float(groups[1].mean() - groups[0].mean())
        elif not groups and bh <= 3:
            # Plain rule with no end ticks: bbox width is the best available
            # estimate (slightly overestimates by one pixel).
            length_px = float(bw - 1)
        else:
            continue

        if length_px < 15:
            continue
        if best is None or length_px > best.length_px:
            best = ScaleBar(
                length_px=length_px,
                bbox=(databar_top + ys.start, xs.start,
                      databar_top + ys.stop, xs.stop),
            )
    return best


def _contiguous_groups(indices: np.ndarray) -> list[np.ndarray]:
    """Split a sorted index array into runs of consecutive values."""
    if indices.size == 0:
        return []
    breaks = np.where(np.diff(indices) > 1)[0]
    return np.split(indices, breaks + 1)


def snap_to_nice(value_um: float, rel_tol: float = 0.08) -> Optional[float]:
    """Snap a length to the nearest 1-2-5 x decade value, if close enough."""
    if not np.isfinite(value_um) or value_um <= 0:
        return None
    i = int(np.argmin(np.abs(np.log(_NICE_VALUES_UM) - np.log(value_um))))
    nice = float(_NICE_VALUES_UM[i])
    return nice if abs(nice - value_um) / nice <= rel_tol else None


@dataclass
class SemImage:
    """A loaded SEM micrograph with a validated pixel->micron calibration."""

    path: str
    intensity: np.ndarray          # micrograph only, databar removed
    full_intensity: np.ndarray     # as stored, including databar
    pixel_size_um: float
    pixel_size_source: str         # 'metadata' | 'scalebar' | 'override'
    databar_top: int
    metadata: dict[str, str] = field(default_factory=dict)
    scale_bar: Optional[ScaleBar] = None
    scalebar_agreement: Optional[float] = None  # relative difference vs used size
    stage_tilt_deg: Optional[float] = None
    magnification: Optional[str] = None
    warnings: list[str] = field(default_factory=list)

    @property
    def shape(self) -> tuple[int, int]:
        return self.intensity.shape

    @property
    def height_um(self) -> float:
        return self.intensity.shape[0] * self.pixel_size_um

    @property
    def width_um(self) -> float:
        return self.intensity.shape[1] * self.pixel_size_um

    @property
    def field_area_um2(self) -> float:
        return self.height_um * self.width_um

    def um_to_px(self, microns: float) -> float:
        return microns / self.pixel_size_um

    def px_to_um(self, pixels: float) -> float:
        return pixels * self.pixel_size_um

    def summary(self) -> str:
        sb = (
            f"{self.scale_bar.length_px:.1f}px->{self.scale_bar.snapped_label_um}um"
            if self.scale_bar and self.scale_bar.snapped_label_um
            else "not found"
        )
        agree = (
            f"{100 * self.scalebar_agreement:+.2f}%"
            if self.scalebar_agreement is not None
            else "n/a"
        )
        return (
            f"{self.path}: {self.intensity.shape[1]}x{self.intensity.shape[0]} px  "
            f"{self.pixel_size_um * 1000:.2f} nm/px ({self.pixel_size_source})  "
            f"field {self.width_um:.2f}x{self.height_um:.2f} um  "
            f"mag {self.magnification}  bar {sb}  agreement {agree}"
        )


def load_sem_image(
    path: str,
    pixel_size_um: Optional[float] = None,
    max_scalebar_disagreement: float = 0.05,
    require_untilted: bool = True,
) -> SemImage:
    """Load an SEM image and establish its calibration.

    Parameters
    ----------
    pixel_size_um
        Explicit override. Use for images with no usable metadata or scale bar.
    max_scalebar_disagreement
        If metadata and the measured scale bar disagree by more than this
        relative amount, raise rather than silently proceed -- this is exactly
        the failure that made the original pipeline 15-30x wrong.
    require_untilted
        Warn when the stage was tilted, because a tilt foreshortens one image
        axis and a single isotropic pixel size no longer applies.
    """
    img = Image.open(path)
    Image.MAX_IMAGE_PIXELS = None
    full = _read_intensity(img)
    meta = parse_zeiss_metadata(img)

    databar_top, databar_bottom = find_databar(full)
    micrograph = full[:databar_top, :]
    if micrograph.size == 0:
        raise MetrologyError(f"{path}: databar detection consumed the whole image")

    warnings: list[str] = []

    meta_px = parse_length_um(meta.get("AP_IMAGE_PIXEL_SIZE")) or parse_length_um(
        meta.get("AP_PIXEL_SIZE")
    )

    bar = measure_scale_bar(full, databar_top, databar_bottom)

    # Decide the pixel size.
    if pixel_size_um is not None:
        used, source = float(pixel_size_um), "override"
    elif meta_px is not None:
        used, source = meta_px, "metadata"
    elif bar is not None:
        raise MetrologyError(
            f"{path}: no pixel size in metadata. A scale bar of "
            f"{bar.length_px:.1f} px was found but its label cannot be read "
            f"without OCR -- pass pixel_size_um explicitly."
        )
    else:
        raise MetrologyError(
            f"{path}: no AP_IMAGE_PIXEL_SIZE metadata and no scale bar found; "
            f"pass pixel_size_um explicitly."
        )

    if used <= 0 or not np.isfinite(used):
        raise MetrologyError(f"{path}: non-physical pixel size {used}")

    # Cross-check against the drawn bar.
    agreement = None
    if bar is not None:
        bar.implied_um = bar.length_px * used
        bar.snapped_label_um = snap_to_nice(bar.implied_um)
        if bar.snapped_label_um:
            bar.pixel_size_um = bar.snapped_label_um / bar.length_px
            agreement = (bar.pixel_size_um - used) / used
            if abs(agreement) > max_scalebar_disagreement:
                raise MetrologyError(
                    f"{path}: scale bar disagrees with {source} pixel size by "
                    f"{100 * agreement:+.1f}% (bar {bar.length_px:.1f} px implies "
                    f"{bar.pixel_size_um * 1000:.2f} nm/px, {source} says "
                    f"{used * 1000:.2f} nm/px). Refusing to guess."
                )
        else:
            warnings.append(
                f"scale bar of {bar.length_px:.1f} px implies "
                f"{bar.implied_um:.3f} um, which is not a 1-2-5 value; "
                f"cross-check skipped"
            )
    else:
        warnings.append("no scale bar found; pixel size not cross-checked")

    tilt = parse_angle_deg(meta.get("AP_STAGE_AT_T"))
    if require_untilted and tilt is not None and abs(tilt) > 1.0:
        warnings.append(
            f"stage tilted {tilt:.1f} deg: one image axis is foreshortened by "
            f"cos(tilt)={np.cos(np.radians(tilt)):.3f}; isotropic pixel size "
            f"is not valid"
        )

    mag = meta.get("AP_MAG")
    if mag:
        mag = mag.split("=", 1)[-1].strip()

    # Sanity: metadata field width should equal columns * pixel size.
    meta_w = parse_length_um(meta.get("AP_WIDTH"))
    if meta_w is not None:
        expected = full.shape[1] * used
        if abs(meta_w - expected) / meta_w > 0.02:
            warnings.append(
                f"AP_WIDTH ({meta_w:.2f} um) disagrees with columns*pixel_size "
                f"({expected:.2f} um)"
            )

    return SemImage(
        path=path,
        intensity=micrograph,
        full_intensity=full,
        pixel_size_um=used,
        pixel_size_source=source,
        databar_top=databar_top,
        metadata=meta,
        scale_bar=bar,
        scalebar_agreement=agreement,
        stage_tilt_deg=tilt,
        magnification=mag,
        warnings=warnings,
    )
