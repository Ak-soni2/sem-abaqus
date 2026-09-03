# EXPECTED — written before the run

This is the prediction for the **15 µm** pad, recorded before any result is
seen. A model interpreted only after the fact can accommodate any outcome.

## The paper's observation for this pad

Section 4.2: *a few traces of brittle fracture are identified... the amounts of brittle fracture are significantly lower than that of 30 um pad.* k = 11.

## What this deck predicts

| quantity | predicted |
|---|---|
| per-grain normal force | 4.125e-05 N |
| groove width | 96.0 nm |
| paper's measured chip | 160–230 nm |
| indentation depth | 0.1537 nm |
| h / dc | 0.00192 |
| passes to reach H·dc | 20 |
| **SDV13 outcome** | **MOSTLY DUCTILE, with some brittle elements late in the passes** |

## What would falsify this

Either extreme falsifies it: entirely ductile through all passes, or a brittle fraction at or above the 30 um deck's. This pad should sit strictly between the other two.

## What CANNOT be checked against the paper

- **Force magnitudes.** The WC-Co Johnson-Cook constants `B, n, C, m` and
  `D1..D5` are placeholders — order-of-magnitude values, not measurements.
  Only `A` is derived, from the JH-2 card's own quasi-static strength. Until
  those are calibrated against nanoindentation or scratch data, no force this
  deck prints is quotable.
- **Surface roughness.** The paper's S_a = 21 nm is a statistical surface from
  roughly 20,000 grain crossings. This deck runs 20 passes of one grain.
  There is no S_a to compare.
- **MRR.** The paper's 0.58 mm³/min comes from its analytical eq. (16) over a
  128 mm² spot. This deck removes material from one patch over microseconds.

**SDV13 — the branch map — is the result.** Everything else is diagnostic.
