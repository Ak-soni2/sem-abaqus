# SAG 30 µm pad — run folder

Self-contained. Copy the whole folder to your work directory and run
`run.bat` (Windows) or `run.sh` (Linux). Nothing is referenced outside it.

```
run.bat
```

That runs three stages: `abaqus verify -user_explicit` (can this machine build
a subroutine at all), a **datacheck** (seconds — catches keyword and material
card errors before the queue), then the solve, then the postprocessor.

**Estimated 96 h on 8 cores.**

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

- Abaqus with a working Fortran toolchain (`abaqus verify -user_explicit`)
- `double=both` — non-negotiable, and the failure is silent if omitted
- `vumat_grind2.for` (58 constants, energy criterion) — **already in this
  folder**, not `vumat_grind.for`

## Honest scope

This deck tests **one** thing: whether SDV13 flips in the right order across
the three pads. It is not a numerical reproduction of the paper — force
magnitudes rest on placeholder Johnson-Cook constants, and surface roughness
needs orders of magnitude more grain passes than are simulated here.
See `EXPECTED.md` for the full list of what cannot be compared.
