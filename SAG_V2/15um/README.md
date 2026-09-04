# SAG_V2 · 15 µm pad

**Job name: `sagv2_15um`.** If the Abaqus console prints `Abaqus JOB sagv2_15um`, this is
the deck running. Anything else — `micro_15um` especially — is an older copy
from a previous round.

Self-contained. Copy the folder, run `run.bat` (Windows) or `run.sh`.

## What this deck cuts

| | |
|---|---|
| cut depth | **195 nm** — the paper's *measured* chip thickness for this pad |
| dc | 80 nm (measured, sections 4.2) |
| **chip / dc** | **2.44** |
| elements through the cut | 12.2 at 16.0 nm |
| passes over one track | 16 |
| estimate | ~162 h on 8 cores |

The paper measures chips directly and compares them against dc — that
comparison *is* its result. 6 µm gives 60–100 nm (at dc, pure ductile); 15 µm
gives 160–230; 30 µm gives 240–350 (both above dc, fracture observed).

The static Brinell indentation from eqs. 11–12 is a *different quantity* —
0.18 nm on the 30 µm pad, three orders below the measurement, 1/90th of one
element and smaller than a WC unit cell. A deck built on it runs for days with
a flat energy history. That is what the earlier rounds did.

## First 30 seconds tell you it is working

The packager prints the model mass. It should read **4.64981e-14** — workpiece
plus the grain. If it reads 4.02778e-14, the grain has no mass and the deck
is an old one.

Then watch KINETIC ENERGY. It should leave zero within the first few output
frames, once the grain closes its 9.8 nm standoff — about 5%
of the way through step LOAD. Energy that stays at exactly zero past the first
frame means nothing is touching.

## What to plot

**SDV13 is the result**: 1 = ductile, 2 = brittle. Plot it after every pass,
not only at the end — the criterion accumulates, so *when* a point flips is the
physics.

`postprocess_odb.py` runs automatically and writes `sagv2_15um_sdv13.csv` (brittle
fraction per frame per pass) and `sagv2_15um_summary.json`.

## Expected outcome

**Mostly ductile, some brittle.** chip/dc = 2.44. The paper sees *a few traces of brittle fracture ... significantly lower than the 30 µm pad*. This deck should sit strictly between the other two.

## Caveats

- Johnson-Cook constants for WC-Co are **placeholders** except `A`. SDV13 is
  quotable; force magnitudes are not.
- The energy criterion is regularised by element length, so PSI is calibrated
  **for this mesh** (16.0 nm through the depth).
