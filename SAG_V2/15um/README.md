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

**Check the job name, not the mass.** The console prints `Abaqus JOB sagv2_15um`.
That is the only identifier that distinguishes this build from an earlier one.

The packager reports a model mass of **4.02778e-14** — the workpiece alone.
That is correct and expected: the grain is a rigid body of R3D3 facets, which
carry no volume, and Abaqus permits a massless rigid body when every
translational dof is constrained. All three are, in every step here. (An
earlier README told you to expect a larger number. That was wrong — the
`*Mass` card it referred to was silently ignored by Abaqus, and a second
attempt to add mass properly aborted the input processor. The mass was only
ever wanted as a build identifier, and the job name does that job.)

Then watch KINETIC ENERGY. It is zero while the grain closes its
9.75 nm standoff, which the smooth-step ramp covers by
**13.5%** of step LOAD — output **frame 3 of 20**.
It must be non-zero by then. If frame 4 is still exactly zero,
nothing is touching and the run is not worth continuing.

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
