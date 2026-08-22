# Improvements — single abrasive, meshing, launchers, figures

What was asked for, what was found, what changed, and how each change was
checked. Every number below was measured from the files in this repo, before
and after.

---

## 1. "Single abrasive: make sure only workpiece and abrasive exist"

**It did not.** The deck named `1_single_abrasive` contained a full bond rim:

| | before | after |
|---|---|---|
| bond rim facets (R3D4) | **2,812** | **0** |
| grit facets (R3D3) | 116 | 116 |
| rigid elements in the wheel | 2,928 | **116** |

The rim was 96% of the rigid mesh, and it was not inert:

* it was listed in the contact surface —
  `*Surface, name=A_WHEEL_SURF` named `ES_BOND_OUTER` as well as `ES_GRITS`;
* the step used `*Contact Inclusions, ALL EXTERIOR`, so all 2,812 facets were
  in general contact regardless;
* `rim_mass_properties()` derived the rigid body's entire mass and inertia from
  the rim volume.

It also could not be switched off: `build_rim_shell()` was called
unconditionally, and `rim_depth_mm` was validated to be strictly positive.

**Change.** A new `include_bond` flag on `DeckParams`, threaded to
`write_rigid_wheel_inp()`, off for the single-abrasive deck. The header now
reads `ABRASIVE : ONE discrete rigid body -- grit facets (R3D3) ONLY, no bond
rim.`

The rim's **mass and inertia are deliberately kept**. They stand for the
spindle the grit is mounted on; every reference-node DOF is velocity-driven, so
the tensor only has to stay positive definite, and dropping it would change the
rigid-body dynamics for no benefit. `verify_rigid_deck2.py` still integrates
that tensor numerically and matches it to **0.000%**.

**On custom trajectory.** `semgrit_multi/trajectory.py` can read a measured
groove and does so correctly, but it feeds the *chip-thickness field*, not the
boundary condition — the simulated grit still rides the ideal arc. That is
stated three times in its own docstring and is a deliberate design decision, so
it was left alone. `deck_amplitudes()` writes the `*Amplitude` tables that would
drive a reference node along a measured path; wiring them in replaces the
rotation-plus-infeed BC and changes the mechanics, so it remains a decision
rather than a default.

---

## 2 & 3. "Check the meshing of both is correct / accurate results"

Three real defects, all measured from the shipped `.inp`.

### a. Element aspect ratio was 51:1

The surface brick was 0.30 × 1.5385 × 0.030 µm. `build_deck.py` already warned
above 10:1 in a comment that named this exact consequence — a groove "only a
couple of elements wide" — and then shipped 51:1 anyway.

### b. The critical depth of cut was spanned by fewer than two elements on SiC

`dc` is the whole point of the model: below it removal is ductile, above it
brittle. The transition cannot be located more precisely than the element that
straddles it.

| | sandstone | SiC |
|---|---|---|
| `dc` | 87.75 nm | 52.92 nm |
| depth element (before) | 0.030 µm | 0.030 µm |
| **elements across `dc`** | **2.93** | **1.76** |

Both decks shipped the *same* mesh, so the SiC deck had been sized for
sandstone's `dc` and never revisited.

### c. The cutting-direction element was larger than the cut was deep

0.30 µm elements against a 0.20 µm depth of cut — chip curl and segmentation
cannot form at that ratio.

### What changed

The depth element is now **derived from each material's own `dc`**
(`ELEMENTS_PER_DC = 5.0`), not hardcoded, so the next material is sized
correctly without editing anything. The cutting and axial sizes came down to
0.15 µm.

Making the axial direction fine everywhere would have been unaffordable, so
`WorkpieceBlock` gained **axial grading** (`width_band_mm`, `width_growth`) to
match the depth grading it already had: fine columns down the groove lane,
coarsening to the side faces. The groove is a few microns wide in a 20 µm block,
so the fine columns are only paid for where something happens.

| | before | sandstone after | SiC after |
|---|---|---|---|
| element (µm) | 0.30 × 1.5385 × 0.030 | 0.15 × 0.15 × 0.0175 | 0.15 × 0.15 × 0.0106 |
| **aspect ratio** | **51.3 : 1** | **8.5 : 1** | **14.2 : 1** |
| **elements across `dc`** | 2.93 / 1.76 | **5.00** | **5.00** |
| workpiece elements | 54,080 | 612,480 | 872,320 |
| est. wall clock, 8 cores | 0.48 h | **9.3 h** | **152.8 h** |

**The cost is real and is not a mesh problem.** SiC's dilatational wave speed is
1.19e7 mm/s against sandstone's 1.76e6 — 6.7× — so the same mesh gives SiC a
6.7× smaller stable increment. Refining further will not help and coarsening is
what produced the 1.76 elements/`dc` in the first place. Raising the mass
scaling would buy time, but these runs already report kinetic energy 320–56,000×
internal, so more inertia is the last thing they need.

The ductile/brittle split moved accordingly, and now behaves as the physics
predicts — SiC, the harder and tougher material, comes out far more brittle:

| deck | cut elements | ductile | brittle |
|---|---|---|---|
| sandstone multi | 4,322 | 2,554 (59%) | 1,768 |
| SiC multi | 6,796 | 1,432 (21%) | 5,364 |

**Not changed, deliberately:** the workpiece dimensions. Shrinking the block
would have brought SiC to ~20 h, but fewer grits would engage and the
multi-abrasive case would stop being comparable to the archived runs.

**Mesh quality was already clean** and stayed that way — no zero-area facets,
no duplicate nodes, no coincident facets. `CAE_SCRIPTS/fix_coincident_facets.py`
would delete 0 elements from these decks; it targets an older per-grit
architecture.

---

## 4. "run.bat should be correct for all"

25 launcher scripts audited. `double=both` and `interactive` were already
correct everywhere, and `3_energy_criterion` correctly used `vumat_grind2.for`.
Five real problems:

**a. Blocker — the PROBE decks could not run at all.** Both `PROBE_sandstone/`
and `PROBE_silicon_carbide/` passed `user=vumat_grind.for`, and neither folder
contained that file. Abaqus resolves `user=` against the working directory, so
both aborted at compile. This is the deck that proves the field variable reaches
the material points; `EXPECTED.md` notes that if it does not, packages 2 and 3
silently run 100% ductile for hours per material. Root cause was in the
generator — `_make_probe.py` wrote the launcher but never copied the subroutine,
unlike `_make_run_packages.py`, which did. Fixed at the generator.

**b. No script cd'd to its own directory.** Double-clicking worked; invoking
`RUN_ME\1_single_abrasive\run.bat` from the repo root ran Abaqus in the root,
where neither the deck nor the `.for` resolves. Added `cd /d "%~dp0"` (and
`cd "$(dirname "$0")"` in the `.sh`).

**c. A failed datacheck did not stop the solve.** The scripts ran datacheck and
then the solve unconditionally, so a deck that died in preprocessing still fired
the solve and reserved the licence tokens. Added `if errorlevel 1 exit /b 1`,
and `set -e` in the shell versions.

**d. The ablation arms shipped no post-processing script.** The README documented
`abaqus python <name>_postprocess_odb.py`, which for the six ablation folders did
not exist — the comparison the model's claim rests on had no scripted extraction
path. Now written into each arm.

**e. Wasted licence tokens.** The ablation arms ran datacheck at `cpus=8`
(12 tokens) where the other decks correctly used `cpus=1` (5). Datacheck is
serial. Aligned.

Also fixed: the README's copy-paste block said `cd forced_brittle` where the
real folders are numbered (`3_forced_brittle`), so pasting it failed.

---

## 5. Figures

### The core problem: no plots were ever produced

`postprocess_odb.py` has drawn figures all along — inside
`try: import matplotlib`. Abaqus' bundled Python has no matplotlib, so the whole
block has been silently skipped on every run. That is why the six completed jobs
are documented in `obd results/` by **photographs of a screen**.

**Fix: `REPOST/plots.py`**, which reads the CSVs those runs already wrote using
the host Python (matplotlib 3.10.9 is installed), so the figures exist without
Abaqus and without re-running anything. 13 figures were generated from the
existing archive immediately.

### What was wrong with the figures themselves

* **The energy plot conveyed nothing.** ALLKE/ALLIE reaches 56,330 on these
  runs; on a shared linear axis that is one flat line and three traces pinned at
  zero. Now log-scaled, five channels, spanning eight decades legibly.
* **Units were left in SI**, so plots carried a `1e-6` corner offset and
  nanometre data sat in the fifth decimal. Now scaled at the axis: µs, mN, nm.
* **`force_vs_h` contradicted itself** — x-axis in mm, `dc` legend in nm.
* **The branch map used `coolwarm`**, a continuous diverging map, for four
  discrete categories: deleted and unset rendered as near-identical pale blues,
  ductile and brittle as near-identical pale reds. That is the one distinction
  the figure exists to show. Now a discrete `ListedColormap` with a labelled
  colourbar, and axes in µm rather than element indices.
* **Nothing was annotated.** First contact, peak time and the engaged window
  were all computed, written to JSON, and left off the plot. Now drawn.
* **Red/green** was used for ductile/brittle — invisible to deuteranopes. Now
  blue/vermillion (Okabe–Ito), which also survives greyscale.
* **`compare_all.png` did not exist.** There was no cross-deck figure anywhere;
  the closest thing was an ASCII table printed to a terminal. Now four panels —
  specific energy, peak force, material removed, artificial-energy fraction —
  with all six jobs side by side.

### Two problems the figures now report rather than hide

Both were found while building the comparison plot, and both were previously
unrecorded:

1. **Four of the six archived "results" are two datasets filed twice.** The
   energy-criterion and multi-abrasive force CSVs are **md5-identical** in both
   materials, although their summary JSONs name different `.odb` files. Only
   **four** of the six runs are distinct. Either the wrong `.odb` was
   post-processed twice, or `SWMODE` did not take effect. `plots.py` hashes the
   CSVs, hatches the duplicated bars and says so in the figure title.

2. **Artificial (hourglass) energy is 31–39% of internal energy** on every run,
   where under 5% is the usual bar — hourglassing is carrying about a third of
   the load. Kinetic energy runs 320–56,000× internal, so mass scaling is
   dominating. Both thresholds are now drawn on the figures. **These numbers
   should be resolved before any force or specific energy from those runs is
   quoted.**

---

## Verification

Everything was re-run, not assumed.

| gate | result |
|---|---|
| `verify_all.py` | **36 passed, 0 failed** |
| `verify_rigid_deck.py` (bondless deck) | **0 failures, 0 warnings** |
| `verify_rigid_deck2.py` (bondless deck) | **0 failures, 0 warnings** |
| `_check_presets.py` | both decks reproduce |
| `REPOST/plots.py --demo` | self-check passes |
| all 16 rebuilt decks | every build gate passed before shipping |

Per-deck gates on the rebuild: 35, 82, 38, 81, 39, 34, 34, 37 checks, all PASS,
for each material. `_make_run_packages.py` refuses to ship a package that does
not verify, and did not have to.

**Two verifier bugs were fixed** — found because the bondless deck exposed them:

* `verify_rigid_deck.py` and `verify_rigid_deck2.py` both crashed on a deck with
  no R3D4 shell. The bond-derived checks are now skipped **and say so on the
  console**, rather than passing silently — a gate that reports PASS on a file it
  never read is worse than one that admits it did not look.
* `verify_rigid_deck2.py` rejected `*Initial Conditions` as an unrecognised
  keyword, failing the multi-abrasive and energy decks — the two carrying the
  novel physics. That keyword is how the swept chip thickness reaches the
  material points. The verifier was wrong, not the decks.

**Regression checked:** the ungraded axial mesh is byte-identical to the previous
`linspace`, so decks that do not opt in are unaffected.

**A third bug was fixed on the way through.** `3_energy_criterion` and all six
ablation arms failed `verify_rigid_deck2.py`'s "report `size_bytes` matches the
file on disk" check. Those decks are the multi-abrasive deck with `PROPS(56)`
(or `(58)`) rewritten, so they are a few hundred bytes different — but they
shipped a verbatim copy of the multi deck's report, which therefore described a
different file. The report is now re-stamped with the arm's own size and path.
Both decks now report **0 failures**.

### The `.inp` files are no longer tracked in git

Refining the mesh took the decks from ~10 MB to 88–155 MB, and **eleven of them
are over GitHub's 100 MB hard file limit** — a push with them in is rejected
outright. They are also the most regenerable artefact in the tree, so
`.gitignore` now excludes `RUN_ME*/**/*.inp` (plus the `.npy` field arrays and
the copied `grain_library.pkl`). Rebuild them with:

```
python _make_run_packages.py --all
```

The build is deterministic and seeded, and refuses to ship a package that does
not pass its gates. **Everything around the decks is still tracked** —
`run.bat`/`run.sh`, the subroutines, every report JSON, the placements CSVs, the
post-processing scripts, `MANIFEST.json` and `README.md` — so the deliverable's
shape, provenance and exact commands remain in git and reviewable.

### Known-outstanding

* Johnson–Cook constants remain placeholders for both materials except `A`.
* The 40 GPa SiC peak Mises question (2.8× HEL) is still open; `hotspot.py` is
  the tool for it and has not been run against the archived `.odb` files.
