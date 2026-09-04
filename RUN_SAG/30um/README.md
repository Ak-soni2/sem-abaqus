# SAG 30 µm pad — run folder

Self-contained. Copy the whole folder to your work directory and run
`run.bat` (Windows) or `run.sh` (Linux). Nothing is referenced outside it.

```
run.bat
```

That runs three stages: `abaqus verify -user_exp` (can this machine build
a subroutine at all), a **datacheck** (seconds — catches keyword and material
card errors before the queue), then the solve, then the postprocessor.

**Estimated 96 h on 8 cores.**

## What this deck cuts, and why

The cut depth is the **measured chip thickness** for this pad, not the static
Brinell indentation. The paper measures chips directly (section 4.2, Fig. 17)
and compares them against `dc` — that comparison *is* its result:

| pad | measured chip | chip / dc | paper observed |
|---|---|---|---|
| 6 µm | 60–100 nm | **1.00** | pure ductile |
| 15 µm | 160–230 nm | **2.44** | brittle + plastic |
| 30 µm | 240–350 nm | **3.69** | brittle + plastic |

The Brinell indentation from eqs. 11–12 is a different quantity: 0.18 nm for
the 30 µm pad, three orders below what the paper measured for the same pad. It
is also unmodellable — 1/90th of one element and smaller than a WC unit cell —
so a deck built on it runs for days with the energy history flat.

## Read EXPECTED.md first

It records what this deck predicts, written before the run. The point of the
exercise is a falsifiable test, and that only works if the prediction is on
record beforehand.

## What to look at

`postprocess_odb.py` runs automatically and writes:

- `micro_30um_sdv13.csv` — brittle fraction for **every frame of every pass**
- `micro_30um_summary.json` — the pass at which brittle elements first appear

**SDV13 is the result**: 1 = ductile, 2 = brittle. Plot it in CAE after each
pass, not only at the end — the energy criterion accumulates, so *when* a point
flips is the physics.

## Requirements

- Abaqus with a working Fortran toolchain (`abaqus verify -user_exp`)
- `double=both` — non-negotiable, and the failure is silent if omitted
- `vumat_grind2.for` (58 constants, energy criterion) — **already in this
  folder**, not `vumat_grind.for`

## Honest scope

This deck tests **one** thing: whether SDV13 flips in the right order across
the three pads. It is not a numerical reproduction of the paper — force
magnitudes rest on placeholder Johnson-Cook constants, and surface roughness
needs orders of magnitude more grain passes than are simulated here.
See `EXPECTED.md` for the full list of what cannot be compared.
