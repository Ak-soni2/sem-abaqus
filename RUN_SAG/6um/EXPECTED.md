# EXPECTED — written before the run

This is the prediction for the **6 µm** pad, recorded before any result is
seen. A model interpreted only after the fact can accommodate any outcome.

## The paper's observation for this pad

Section 4.2: *no brittle fracture is found in this case and the material is removed through pure ductile mode.* The ratio k = d_WC/d_g = 4.41 is below the paper's threshold of 5.

## What this deck predicts

| quantity | predicted |
|---|---|
| per-grain normal force | 1.697e-05 N |
| groove width | 61.6 nm |
| paper's measured chip | 60–100 nm |
| indentation depth | 0.1580 nm |
| h / dc | 0.00198 |
| passes to reach H·dc | 20 |
| **SDV13 outcome** | **ENTIRELY DUCTILE — SDV13 = 1 everywhere, all passes** |

## What would falsify this

**Any brittle element.** The paper reports pure ductile removal for this pad and the model predicts h/dc = 0.002, so a single SDV13 = 2 is a real disagreement, not noise.

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
