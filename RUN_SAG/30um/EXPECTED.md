# EXPECTED — written before the run

This is the prediction for the **30 µm** pad, recorded before any result is
seen. A model interpreted only after the fact can accommodate any outcome.

## The paper's observation for this pad

Section 4.2: *although the material is mostly removed through plastic flow, several traces of brittle fracture are observed.* k = 22, far above the threshold of 5.

## What this deck predicts

| quantity | predicted |
|---|---|
| per-grain normal force | 9.505e-05 N |
| groove width | 145.7 nm |
| paper's measured chip | 240–350 nm |
| indentation depth | 0.1770 nm |
| h / dc | 0.00221 |
| passes to reach H·dc | 20 |
| **SDV13 outcome** | **BRITTLE REGIONS PRESENT — SDV13 = 2 somewhere** |

## What would falsify this

**Entirely ductile through all passes.** That would mean the energy criterion accumulates too slowly to reach H·dc under a realistic per-grain load — a finding about the criterion's rate, not a bug in the deck. It is the single most informative negative result available here, which is why this pad is the one to run first.

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
