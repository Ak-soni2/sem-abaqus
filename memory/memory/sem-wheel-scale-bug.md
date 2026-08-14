---
name: sem-wheel-scale-bug
description: Zeiss SEM TIFFs in this project carry exact pixel size in tag 34118; the original notebook's scale-bar detection was 15-30x wrong
metadata:
  type: project
---

The SEM `.tif` files in this project are **Zeiss SmartSEM** images that embed full instrument state in TIFF tag **34118** (~1508 lines of `AP_`/`DP_`/`SV_` key/value pairs). `AP_IMAGE_PIXEL_SIZE` gives the exact pixel size — no scale-bar detection or OCR is needed.

Measured facts (verified 2026-07-28 on all 14 images):
- 10 kX images → 29.30 nm/px; 5 kX images → 58.59 nm/px. All 1024x768, databar burned over the bottom 72 rows (micrograph = rows 0..693).
- The original notebook's scale-bar detector thresholded for *bright* objects, but the Zeiss databar is a **white background with a dark bar**, so it locked onto the databar strip (w=1019 px) instead of the bar (68 px). Result: every measurement was **14.9x too small at 10 kX and 29.9x at 5 kX**.
- The true scale bar is an I-beam; the correct length is **tick-centre to tick-centre** (68 px), not the bounding-box width (70 px).
- The TIFF palette reserves 27 non-gray annotation colours (every 9th index). They occupy ~88% of databar pixels but **0% of the micrograph**, so `cv2.imread` does not corrupt grain intensities — but it does make databar/OCR parsing unreliable. Reading raw palette indices gives true gray for the micrograph.

**Why:** trusting the drawn scale bar over the embedded metadata was the single largest error in the project, and it silently propagated into the CAD and wheel dimensions.

**How to apply:** read pixel size from tag 34118 as the primary source; use tick-centre scale-bar measurement only as a cross-check (agrees to ~0.4%). Also check `AP_STAGE_AT_T` for stage tilt before trusting Y-axis distances. See [[sem-grinding-wheel-abaqus-goal]].
