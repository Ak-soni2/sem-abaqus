# sem-abaqus — Full Project Context

Repo: https://github.com/Ak-soni2/sem-abaqus (branch: `full-pipeline-sic`)
Purpose of this file: single-document handoff so any AI/human can understand the entire project without reading the codebase first.

> **See also `IMPROVEMENTS.md`** for the most recent round of work (single-abrasive
> bond removal, mesh refinement, launcher fixes, and the figure kit), including
> before/after numbers and two data-integrity problems found in the archived
> results. Sections 4, 5 and 8 below have been updated to match.

## 1. What this project is

A pipeline that turns SEM (scanning electron microscope) micrographs of abrasive grinding grit into a working Abaqus/Explicit FEA simulation of a grinding wheel cutting a workpiece, using a custom hybrid ductile/brittle material model (VUMAT). End goal: an Abaqus-importable wheel model that predicts material removal mode (ductile plastic flow vs. brittle fracture) based on measured grit geometry and grinding kinematics — not just a CAD/measurement exercise.

Pipeline stages:
1. SEM micrograph (TIFF) → grain segmentation → 3D grain geometry
2. Grains assembled onto a wheel (sector or full) as one rigid body
3. Wheel + workpiece + hybrid material assembled into an Abaqus `.inp` deck
4. Deck run in Abaqus/Explicit with a custom VUMAT subroutine
5. Results post-processed (forces, energy, SDVs) and verified

## 2. Physics / mechanics

- **Ductile-vs-brittle material removal**: below a critical depth of cut `dc`, brittle solids (rock, ceramic) are removed via plastic flow (ductile regime, smooth/low-damage); above `dc`, removal is via fracture/chipping (brittle regime, higher damage/rougher). This is the central hypothesis the whole hybrid VUMAT encodes.
- **`dc` formulas** — three coexist in code, not reconciled to one:
  - `dc = λc(H/E)^0.5 (Kc/H)²`
  - Bifano 1991: `dc = λc(E/H)(Kc/H)²` — differs from the above by `(E/H)^1.5` (~17× on sandstone)
  - New energy-based criterion (from `_derive_grind2.py`/`vumat_grind2.for`): local plastic work per unit area vs. fracture energy: `Wp·Lc ≥ Ψ·Kc²/E`, calibrated by default to trip at the same `dc` as the others.
- **Brittle branch constitutive model**: Johnson-Holmquist II (JH-2) — pressure-dependent strength surface, damage accumulation, bulking, `SFMAX`-capped fractured strength. Reference implementation `vumat_jh2.for` independently verified against the JH94 bulking benchmark to 1.3%.
- **Ductile branch constitutive model**: Johnson-Cook (JC) flow/damage + a Strain-Gradient Enhancement (SGE) size-effect term (Taylor/GND-based), from three Yadav et al. papers (2022/2024/2026):
  `σe = σ_JC · sqrt(1 + (r'·η·b·(M·α·G)² / σ_JC²)^Λ)`, with `η = 4·εp/h`.
- **Damage continuity**: when material flips from ductile to brittle mode (or vice versa) mid-simulation, JC damage is inherited as JH-2 damage (state carried across the transition) rather than reset or double-counted.
- **Materials modeled**:
  - `sandstone` — original placeholder JH-2 card
  - `silicon_carbide` — SiC-N Holmquist/Johnson card (note: this material was mislabeled upstream in old code as "monocrystalline silicon"; corrected in `semgrit/materials.py`)
  - Johnson-Cook constants are explicit **placeholders** for both materials except constant `A`, which is derived to match each material's own JH-2 quasi-static compressive strength (so ductile and brittle branches agree at the transition point). Real calibration against nanoindentation/scratch data is still pending.
- **Grinding kinematics**: rotational surface speed (e.g. −1200 rad/s ⇒ 30 m/s at r=25mm) plus radial infeed for depth of cut `ae`; per-grit chip thickness (`H0`, `HG`, `RTIP`) computed from grit placement + wheel motion and fed into the 56-value VUMAT material card (`semgrit/hybrid.py`). The Fortran subroutine itself carries **no** process/kinematic knowledge — all of that is precomputed in Python and passed in as PROPS.

## 3. Code architecture

### `semgrit/` — core package (SEM measurement → Abaqus deck)
- `metrology.py` — Zeiss SmartSEM TIFF tag 34118 calibration; fixes a documented 15–30× scale bug present in the original notebook.
- `segment.py` — watershed segmentation of grains (continuous elevation map, evidence-based split retention).
- `measure.py` — 25 shape descriptors per grain (Crofton perimeter, Feret diameter, etc.).
- `grain3d.py` — lofts each 2D grain outline into a watertight 3D polyhedron (ear-clip + prism→tet meshing).
- `wheel.py`, `rigid_wheel.py`, `wheel_workpiece.py` — assemble grains onto a wheel. `rigid_wheel.py` is the **current entry point**: emits the whole wheel (bond + all grits) as a single discrete-rigid body, workpiece kept deformable. (Deliberate choice — see memory note below — over per-grit rigid bodies, for stable time-increment and simplicity.)
- `abaqus.py`, `build_deck.py` — `.inp` file writer. `build_deck.py` is the single entry point for deck building (`DeckParams` controls sector mode, grit density mode, hybrid material on/off, workpiece on/off). Has a `PRESETS` dict reproduced byte-for-byte by `_check_presets.py`.
- `hybrid.py` — computes chip-thickness kinematics (`H0`, `HG`, `RTIP`) and the full 56-value hybrid VUMAT material card from grit placement + wheel motion.
- `materials.py` — material registry: JH-2 cards, JC+SGE constants, hardness/toughness values, with unit-consistency self-checks.
- `step.py`, `cadviewer.py` — STEP/STL CAD export (for viewing in SOLIDWORKS etc.).
- `verify.py` — shared verification helpers.

### `semgrit_multi/` — multi-abrasive extension
- `envelope.py` — swept chip-thickness envelope across many grits with time-ordered shadowing.
- `fieldinject.py` — injects field variables into the `.inp` deck.
- `trajectory.py` — replays measured (real) grit trajectories.
- `swmode.py`, `ablate.py` — ablation-study arm switching (forced-ductile / forced-brittle / geometric-hybrid comparison arms).

### Root `_*.py` scripts — one-off build/derivation/test utilities (not part of the importable package)
- `_derive_grind2.py` — derives `vumat_grind2.for` from `vumat_grind.for` via **anchored text surgery** (not hand-copy-paste), so the ~1000-line shared body can't silently drift between the two files. Has a `--check` mode that fails CI-style if the derived file is stale.
- `_make_notebook.py` / `_make_notebook2.py` — regenerate the self-contained Colab notebooks by embedding the `semgrit`/`semgrit_multi` packages as base64 blobs inside the `.ipynb`.
- `_make_run_packages.py`, `_make_final_rigid.py`, `_make_probe.py`, `_make_single_grit.py`, `_build_arc30.py`, `_build_fullwheel.py`, `_build_presentation.py`, `_make_present.py`, `_make_d50.py` — build the various deliverable deck folders (`RUN_ME*/`, `FINAL_RIGID/`, `PROBE_*/`, `SINGLE_GRIT/`, presentation assets).
- `_check_presets.py` — cheapest regression test: reproduces `FINAL_RIGID`/`SINGLE_GRIT` decks byte-for-byte from `build_deck.PRESETS`.
- `_fix_grinding_job.py`, `_fix_grinding_job2.py`, `_fix_grinding_job4.py` — iterative fixups applied to specific job decks during debugging.
- `_cadv_probe.py`, `_add_visuals.py` — CAD viewer probing / visualization helpers.
- `_hybrid_test/` — folder with driver/reference implementations for testing the hybrid VUMAT in isolation.

### Fortran VUMATs (root + copied into `RUN_ME*/` folders)
- `vumat_jh2.for` — standalone JH-2 reference. 17 constants, 12 SDVs (state-dependent variables).
- `vumat_grind.for` — the hybrid ductile/brittle law. 56 PROPS, 20 SDVs. PROPS 1–21 are a **byte-exact prefix** of `vumat_jh2.for`'s material card, so setting `PROPS(56)=3` reproduces plain JH-2 to 0 ULP (used as a correctness anchor).
- `vumat_grind2.for` — derived from `vumat_grind.for` (via `_derive_grind2.py`), adds the local energy-based ductile/brittle criterion. 58 PROPS, 22 SDVs (adds `SWMODE`, `PSI`).
- `vumat_jc_damage (1).for` — earlier, superseded exploration (JC-damage-only, no JH-2 branch). Kept for reference; documented as broken because naive `(1-D)` stress scaling "freezes" damage and material never actually fails.

### Notebooks
- `SEM_TO_ABAQUS_GRINDING_WHEEL.ipynb`, `SEM_TO_ABAQUS_MULTI_ABRASIVE.ipynb` — **current** deliverables. Self-contained Colab notebooks with `semgrit`/`semgrit_multi` embedded as base64. Never hand-edited directly — always regenerated via `_make_notebook*.py`.
- `SEM_WHEEL (1).ipynb` / `sem_wheel (1).py` — original/superseded notebook; this is what the `semgrit/` package replaced. Kept for reference only.
- `SEM_WHEEL_semgrit_COLAB.ipynb` — earlier intermediate Colab variant.

## 4. Abaqus model setup

- **Workpiece**: deformable, C3D8R elements (reduced-integration hex). Mesh is graded in **both** depth (fine at the ground face) and axially (fine down the groove lane, `width_band_mm`). The depth element is **derived from each material's own `dc`** via `ELEMENTS_PER_DC = 5.0` in `_make_run_packages.py`, not hardcoded — sandstone 0.15 × 0.15 × 0.0175 µm (612,480 elements), SiC 0.15 × 0.15 × 0.0106 µm (872,320 elements). Element aspect ratio 8.5:1 and 14.2:1. *(Previously a flat 0.30 × 1.5385 × 0.030 µm for both materials: 51:1 aspect ratio and only 1.76 elements across SiC's `dc`.)*
- **Wheel**: R3D3/R3D4 discrete rigid elements for both grits and bond. **Entire wheel merged into one `*Rigid Body`** with a single reference node on the wheel axis — one rotational velocity BC drives the whole assembly (deliberate simplification over per-grit rigid bodies).
- **`include_bond`** (`DeckParams`): the single-abrasive deck sets this `False`, so only the grit facets and the workpiece exist — no bond rim. The rim's mass/inertia are still written (they stand for the spindle; every reference-node DOF is velocity-driven). The multi-abrasive and ablation decks keep the rim.
- **Contact**: general contact; `A_GRITS_ENGAGE_SURF` set restricts contact-pair computation to grits that actually reach the workpiece block (cost optimization).
- **Boundary conditions**: rotational velocity (`VR3`) for surface speed, plus radial infeed velocity components for depth of cut `ae`; workpiece fixed via encastre/side/end node sets.
- **Material assignment**: `*User Material` (VUMAT) card with 56 or 58 constants, **must be written 8 values per line** — 4-per-line is silently rejected by Abaqus (a real past failure; the exact error message is preserved in `error/*.dat` files for reference). `*Depvar, delete=12` ties element deletion to the JH-2 damage-to-1 SDV flag.
- **Solver**: Abaqus/Explicit, `double=both` precision is **mandatory** — nanometre-scale chip-thickness comparisons on a 25mm-radius wheel need double precision throughout, not just for the packer.
- **`CAE_SCRIPTS/`** — three hand-written Abaqus/CAE Python scripts (not pipeline-generated; only recently un-gitignored — it had been wrongly treated as generated output):
  - `fix_coincident_facets.py` — deletes near-zero-area rigid contact facets from coincident-node meshing artifacts, without touching materials/BCs/VUMAT assignment. Has a `MAX_FRACTION` safety backstop and a `DRY_RUN` mode.
  - `set_explicit_no_deletion.py` — bulk-sets every part to the Explicit element library with deletion off, avoiding ~100 manual "Element Type" dialog operations in CAE. Deliberately never touches step definitions (considered too destructive to auto-fix).
  - `verify_explicit_settings.py` — writes the model out to `.inp` and parses it back to confirm element types/section controls/step procedure actually took effect, rather than trusting the CAE GUI dialog state.

## 5. Deliverable packages (folder-by-folder)

> **The `.inp` files are not in git.** After the mesh refinement they are 88–155 MB each and eleven exceed GitHub's 100 MB limit, so `RUN_ME*/**/*.inp` (and the `.npy` field arrays) are gitignored. Rebuild with `python _make_run_packages.py --all` — deterministic, seeded, and gated. Everything else in those folders is tracked.

- `RUN_ME/` (sandstone) and `RUN_ME_SIC/` (silicon carbide) — each contains 4 sub-decks: `1_single_abrasive`, `2_multi_abrasive`, `3_energy_criterion`, `4_ablation` (3 comparison arms: forced-ductile, forced-brittle [bit-identical to plain JH-2], geometric-hybrid). This 2×4 matrix is the core experiment set. Has `MANIFEST.json` and `README.md`.
- `PROBE_sandstone/`, `PROBE_silicon_carbide/` — minimal 8-element smoke-test decks with pre-computed expected SDV values in `EXPECTED.md`, used to confirm Abaqus field-variable plumbing reaches the VUMAT correctly before spending hours on a full job.
- `FINAL_RIGID/`, `SINGLE_GRIT/`, `WHEEL_FIXED/` — earlier geometry-only rigid-wheel deliverables predating the hybrid material work (placeholder `*Elastic` workpiece material, no VUMAT/BCs).
- `REPOST/plots.py` — **the figure kit**. Reads the archived `_forces.csv` / `_energy.csv` / `_summary.json` with the *host* Python and writes every figure, including the cross-deck `compare_all.png`. It exists because `postprocess_odb.py` plots only inside `try: import matplotlib`, and Abaqus' bundled Python has none — so no run of this project ever produced a PNG. Run `python REPOST/plots.py [dir]`; `--demo` is a self-check. Figures land in `REPOST/figures/`.
- `REPOST/` — corrected post-processing kit, added after discovering **three bugs** in the six already-completed job results:
  1. stale post-processing scripts were being pulled from `D:\temp` instead of the repo
  2. `find_report()` report auto-resolution picked the wrong material's report for 5 of the 6 jobs
  3. SDV sampling only read the last simulation frame instead of walking every frame (fixed by `hotspot.py`)
  - `REPOST/reports/` and `obd results/` archive the six completed jobs' outputs: force/energy CSVs, JSON summaries, viewport screenshots.

## 6. Verification philosophy

Verification is treated as central to the project, not an afterthought. ~11 root-level `verify_*.py` scripts plus `semgrit/verify.py` and `_check_presets.py`, deliberately implemented as **independent, non-shared-code** cross-checks (so a bug in the main pipeline can't also be baked into its own verifier):
- `verify_all.py` — 26+ unit/integration/export-roundtrip checks on the SEM measurement pipeline.
- `verify_vumat_grind.py` / `verify_vumat_grind2.py` — compiles the Fortran with gfortran and checks results against closed-form algebra, JH94 published benchmarks, and bit-identity with `vumat_jh2.for`.
- `verify_hybrid_deck.py` — confirms the VUMAT's actual ductile/brittle branch choice agrees with what the deck's geometry predicts it should choose.
- `verify_rigid_deck.py` (re-derives geometry from raw node coordinates) and `verify_rigid_deck2.py` (independent keyword-grammar state machine + numerical mass/inertia integration) — two separately-implemented checks that must both pass.
- `verify_envelope.py`, `verify_inp_geometry.py`, `verify_pipeline_A.py` / `verify_pipeline_B.py`, `verify_colab.py` — additional pipeline-stage checks.
- Deliberate **negative controls** exist (e.g. an inward-wound STEP solid must be rejected by the geometry checker) so the suite cannot pass vacuously by never actually triggering failure paths.

## 7. Project history / status timeline (most recent 3 commits, Aug 2026)

1. **`ae2ca76`** (Aug 14, 2026) — "SEM micrographs to a hybrid ductile/brittle Abaqus grinding model." Built the entire hybrid pipeline end-to-end **in code**: `vumat_grind.for`/`vumat_grind2.for`, `semgrit`/`semgrit_multi` packages, all `RUN_ME*`/`PROBE_*` deliverables. At this point, **Abaqus had not yet successfully run any of it** — one prior attempt died at preprocessing due to the 4-per-line `*User Material` card error (now fixed and gated by verification). Every number as of this commit is a statement about code correctness, not simulation results.
2. **`21f5ff1`** (Aug 16, 2026) — "Add REPOST kit; archive the six completed job results." The six planned jobs (sandstone + SiC × single/multi/energy variants) **did successfully run** (61/61 frames, 54,080 elements each). However, the post-processing of those runs had 3 bugs (see REPOST section above), all now fixed. Results archived under `obd results/`.
   - **Open physics question flagged here**: `single_abrasive2_sic` (silicon carbide run) hit a peak Mises stress of ~40 GPa — 2.8× SiC's Hugoniot Elastic Limit (HEL). Not yet determined whether this is legitimate uncapped JH-2 intact-surface behavior, or a strain-gradient-enhancement (SGE) artifact caused by the length scale `h` clamping to ~0 in the rubbing zone. The diagnostic that would resolve this is `SDV19` (FSGE, the strain-gradient enhancement factor) — `hotspot.py` (which reads every frame, not just the last) was written specifically to check this but had not yet been run against the archived `.odb` results as of this commit.
3. **`1a126de`** (Aug 16, 2026, current HEAD) — "Track CAE_SCRIPTS: it was ignored as generated output, but is not." Housekeeping fix: the three hand-written `CAE_SCRIPTS/*.py` fixup/verification scripts were wrongly excluded by `.gitignore` (mistaken for pipeline-generated output); now tracked and documented.

## 8. Known open items / next steps (as of this document)

**Data integrity — read before quoting any archived number:**
- **Four of the six archived "results" are two datasets filed twice.** The energy-criterion and multi-abrasive force CSVs are md5-identical in *both* materials, though their summary JSONs name different `.odb` files. Only **four** distinct runs exist. Either the wrong `.odb` was post-processed twice, or `SWMODE` did not take effect. `REPOST/plots.py` detects this by hashing and hatches the affected bars.
- **Artificial (hourglass) energy is 31–39% of internal energy** on every archived run (bar: <5%), and kinetic energy is 320–56,000× internal (mass scaling dominating). These should be resolved before any force or specific-energy figure is quoted.


- Run `hotspot.py` against the archived six `.odb` results (if still available) to resolve the 40 GPa SiC Mises-stress question via SDV19 (FSGE).
- Johnson-Cook ductile constants remain **placeholders** for both sandstone and silicon carbide (except the derived `A` constant) — pending calibration against real nanoindentation/scratch test data.
- The three `dc` (critical depth of cut) formulas in the codebase are not reconciled into one authoritative form; different scripts may use different ones — check which is active in a given deck before trusting `dc`-dependent results.
- `CAE_SCRIPTS/` existing as hand-written (non-pipeline-generated) fixup scripts suggests CAE-side model setup (contact facet degeneracies, Explicit/element-deletion settings) is still a manual friction point when importing/preparing decks by hand in Abaqus/CAE, rather than something `build_deck.py` fully automates yet. (Note: `fix_coincident_facets.py` would delete **0** elements from the current decks — it targets an older per-grit architecture. The shipped meshes are clean: no zero-area facets, no duplicate nodes, no coincident facets.)
- **The refined decks have not been run.** The 16 rebuilt decks pass every build gate, but the wall clock is now ~9.3 h per sandstone deck and ~153 h per SiC deck on 8 cores (SiC's dilatational wave speed is 6.7× sandstone's, so its stable increment is 6.7× smaller on the same mesh — this is material physics, not a mesh defect). Plan SiC runs accordingly, or raise the core count.
- **The `hs`/`force_vs_h` and `branch_map` figures still have no data**, because the archived runs predate the SDV-reading postprocessor. They will appear on the next run.

## 9. Key file paths for quick reference

- `README.md` — top-level `semgrit` package documentation
- `memory/memory/*.md` — project decision history (e.g. `sem-grinding-wheel-abaqus-goal.md`, `sem-wheel-scale-bug.md`, `wheel-all-discrete-rigid.md`)
- `semgrit/hybrid.py`, `semgrit/materials.py`, `semgrit/build_deck.py`, `semgrit/rigid_wheel.py`
- `vumat_grind.for`, `vumat_grind2.for`, `vumat_jh2.for`
- `_derive_grind2.py`
- `RUN_ME/README.md`, `RUN_ME_SIC/README.md`, `RUN_ME/MANIFEST.json`
- `CAE_SCRIPTS/fix_coincident_facets.py`, `set_explicit_no_deletion.py`, `verify_explicit_settings.py`
- `REPOST/hotspot.py`, `REPOST/postprocess_odb.py`
- `PROBE_sandstone/EXPECTED.md`, `PROBE_silicon_carbide/EXPECTED.md`
