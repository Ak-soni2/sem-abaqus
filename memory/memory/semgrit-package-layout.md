---
name: semgrit-package-layout
description: "semgrit module roles, the one deck-building entry point, the two verifiers, and how to regenerate the Colab notebook"
metadata: 
  node_type: memory
  type: project
  originSessionId: 751e6fdc-4044-47a8-b6d3-ca037006409a
  modified: 2026-08-05T10:43:17.814Z
---

As of 2026-07-28 the Colab notebook (`SEM_WHEEL (1).ipynb` / `sem_wheel (1).py`) was replaced by a `semgrit/` package in the same folder. The notebook files were **left in place** for reference, not deleted.

Modules: `metrology.py` (Zeiss TIFF calibration), `segment.py` (watershed + evidence-based merge), `measure.py` (shape descriptors), `grain3d.py` (lofted grain solids + C3D4 tets), `wheel.py` (sector wheel + grain placement), `abaqus.py` (`.inp` writer), `step.py` (STEP/STL CAD export for SOLIDWORKS), `verify.py`, `cli.py`.

Two output targets with very different scale limits: Abaqus `.inp` handles 10^5 grains via `*Part`/`*Instance`, while a CAD kernel is far slower per body — STEP export is capped (`--step-max-grains`, default 200) and samples grains evenly across the sector. Never raise that cap far without warning the user; SOLIDWORKS chokes on thousands of solid bodies.

Entry points:
- `python -m semgrit analyze "*.tif" -o results --verify`
- `python -m semgrit wheel results -o out --diameter 100 --width 10 --sector 30 --rim-depth 2 --areal-density 40 --grain-element R3D3 --verify`
- `python verify_all.py` — 26-check gate (unit + integration + export round-trip + determinism); exits non-zero on failure.

Environment note: this machine's default `python` is a broken msys64 build (`pip` fails to import). Use `C:\Users\MALAV PAREKH\AppData\Local\Programs\Python\Python312\python.exe`, which has all dependencies installed. `manifold3d` is absent so trimesh boolean ops are unavailable — irrelevant, because grains must stay separate bodies for Abaqus contact.

**Why:** the user asked for changes to be verified repeatedly and completely, so the verification suite is the contract; re-run `verify_all.py` after any change rather than spot-checking.

**Deck building has one entry point (added 2026-07-31):** `build_deck.py` —
`build_deck(DeckParams(...), solids, outdir)` writes the `.inp`, CAE loader, placements
CSV, report JSON and optional STEP/STL. It calls `rigid_wheel.write_rigid_wheel_inp`,
which emits the whole wheel as **one** discrete rigid body (see
[[wheel-all-discrete-rigid]]). `wheel_workpiece.py` is the older per-grit-rigid-body
writer, kept only because earlier decks were built with it. `DeckParams` covers sector
mode (`arc`/`angle`/`full`), grit mode (`concentration`/`areal_density`/`count`/`single`),
workpiece on/off, and CAD.

**Verification is two independent implementations:** `verify_rigid_deck.py` (re-derives
geometry from the node coordinates) and `verify_rigid_deck2.py` (cross-checks header and
report against the mesh, integrates mass and inertia numerically). 83 checks together on
a deck with a workpiece; both must exit 0. `verify_all.py` still gates the measurement
half.

**Cheapest regression test in the project:** `python _check_presets.py` — `build_deck.PRESETS`
reproduces the two Abaqus-validated decks (`FINAL_RIGID/`, `SINGLE_GRIT/`) byte-for-byte on
non-comment content. Run it after touching `rigid_wheel.py`. It must be fed
`WHEEL_FIXED/1_measurements/grain_library.pkl` (96 solids, 712 grits); a freshly measured
library has 548 solids, repacks the rim and reports a difference that is **not** a
regression — that mistake cost two wrong "DIFFERS" verdicts, which is why the path is
hard-coded in the script rather than passed in.

**Deliverable:** `SEM_TO_ABAQUS_GRINDING_WHEEL.ipynb`, a self-contained Colab notebook
with the whole package embedded as a base64 blob. **Never edit it by hand** — change the
package, regenerate with `_make_notebook.py`, then re-run `_run_notebook_test.py`, which
executes it headlessly in a scratch dir with the dev copy kept off `sys.path`. Two
nbformat traps that bit once: `source` lines must keep trailing `\n` or the cell
collapses to one line, and the blob must be spliced in as quoted Python literals.

**How to apply:** key tuned constants that were derived empirically and should not be changed casually — segmentation `min_edge_strength=1.5` (sits in a measured bimodal gap), `h_maxima_um=0.12` (deliberate over-segmentation), Crofton perimeter with 4 directions, and jittered-grid sampling (never truncate Bridson). See [[sem-wheel-scale-bug]] and [[sem-grinding-wheel-abaqus-goal]].
