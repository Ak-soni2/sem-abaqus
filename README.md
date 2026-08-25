# semgrit — SEM abrasive-grain measurement → Abaqus grinding-wheel model

Measures individual abrasive grains in SEM micrographs, reconstructs each one as a
3D solid, and assembles them onto a grinding wheel (or an angular sector of one)
exported as an Abaqus `.inp` deck.

Replaces `SEM_WHEEL (1).ipynb` / `sem_wheel (1).py`.

> **New to the project? Read [`THEORY.md`](THEORY.md) first.** It is a
> self-contained account of the physics and mathematics — every equation the
> model uses, derived or cited, with what may and may not be quoted from the
> results. This file is how to *run* it; `context.md` is how the code is laid
> out; `THEORY.md` is what it *means*.

## Install

```bash
pip install -r requirements.txt
```

## Use

```bash
# 1. measure grains in every SEM image
python -m semgrit analyze "*.tif" -o results --verify

# 2. build a 30 deg wheel sector and export for Abaqus
python -m semgrit wheel results -o wheel \
    --diameter 100 --width 10 --sector 30 --rim-depth 2 \
    --areal-density 40 --grain-element R3D3 --verify

# 3. same wheel, but also as CAD for SOLIDWORKS (small sector -- see the note below)
python -m semgrit wheel results -o cad \
    --diameter 100 --width 10 --sector 3 --rim-depth 1 \
    --areal-density 40 --step --step-max-grains 150 --stl --verify

# individual grain solids as CAD, laid out on a grid for inspection
python -m semgrit analyze "*.tif" -o results --step --step-max-grains 200

# full verification suite (unit + integration + export round-trips)
python verify_all.py

# the constitutive law: compiles vumat_grind.for and exercises it on a single
# material point against closed-form algebra, the published JH-2 benchmarks,
# and vumat_jh2.for itself (needs gfortran)
python verify_vumat_grind.py

# a hybrid deck: does the subroutine agree with the card about which elements
# are ductile?  With no argument it builds its own reference deck.
python verify_hybrid_deck.py [<deck>.inp]
```

`--sector` takes any angle in `(0, 360]`: `360` full wheel, `180` half, `30`, `25`, …
`--rim-depth` models only the outer annulus, which is normally what you want —
grinding is confined to a shallow surface layer, and this is what makes a
high-grain-count model tractable.

## Outputs

| File | Contents |
|---|---|
| `<image>_grains.csv` | 25 shape descriptors per grain, in microns |
| `<image>_report.json` | calibration, segmentation settings, validation results |
| `<image>_segmentation.png` | overlay; green = interior grain, red = border-truncated |
| `summary.json` | pooled size distribution across all images |
| `grain_library.pkl` | reusable measured 3D grain solids |
| `<name>.inp` | Abaqus deck |
| `<name>_placements.csv` | per-grain position, orientation, protrusion |
| `<name>.step` | CAD solid bodies — **open in SOLIDWORKS** (`--step`) |
| `<name>.stl` | mesh body, opens anywhere (`--stl`) |
| `grains.step` | individual grain solids on a grid (`analyze --step`) |

## CAD output for SOLIDWORKS

`--step` writes a STEP (ISO 10303-21) **faceted B-rep**: `FACETED_BREP` solids with
planar `FACE_SURFACE`s bounded by `POLY_LOOP`s. SOLIDWORKS imports this as genuine
solid bodies — one body per grain plus the bond rim, in a multibody part — so you can
section, measure and apply features. An STL by contrast arrives as a mesh body you
cannot really work with, which is why STEP is the default CAD route.

Nothing is lost by faceting: the grains genuinely *are* polyhedra, built by lofting a
measured outline between planar rings, so every face is exactly planar. The rim is
faceted at the same angular resolution as its hex mesh, so the CAD and the FE model
describe the same body.

**Use a small sector for CAD.** A CAD kernel is far slower per body than an FE solver:
10,000 grains is routine in Abaqus and unusable in SOLIDWORKS. `--step-max-grains`
defaults to 200 and samples grains evenly across the sector (not the first N, which
would clump them at θ=0). The `.inp` always keeps every grain regardless; the cap only
affects the CAD file, and the exact counts are reported.

Rough guide: `--sector 3 --step-max-grains 150` gives a 9 MB STEP with 151 solid
bodies. Expect roughly 60 kB and ~1,100 STEP entities per grain.

Verified end to end — the STEP is re-parsed from disk and compared against the FE
model: bond-rim volume agrees to **7e-13** relative, grain volumes to **7e-8**, and each
grain's radial protrusion to **5e-11 mm**.

## How it works

**1. Calibration.** Pixel size is read from the instrument's own metadata
(Zeiss SmartSEM TIFF tag 34118, `AP_IMAGE_PIXEL_SIZE`). The burnt-in scale bar is
measured independently — tick-centre to tick-centre — and its label recovered by
snapping to the nearest 1‑2‑5 value, giving an OCR-free cross-check. If the two
disagree by more than 5%, the run **fails** rather than proceeding.

**2. Segmentation.** Multi-Otsu foreground, watershed on a continuous elevation map
built from the distance transform plus the intensity gradient, seeded by h-maxima.
The watershed deliberately over-segments; splits are then **kept only where the
image supports them** (a real intensity edge along the boundary, or a narrow neck
between two particles). All thresholds are in microns, so 5 kX and 10 kX images are
treated identically.

**3. Measurement.** Descriptors are computed on the *actual* outline, never a convex
hull: area, Crofton perimeter, equivalent diameter, max/min Feret, axes, aspect
ratio, circularity, solidity, convexity, corner count and sharpest corner angle.
Border-truncated grains are measured but excluded from distributions.

**4. Grain solids.** Each grain becomes a lofted polyhedron: the measured silhouette
at mid-height, tapering to smaller faces above and below. Two properties are
verified numerically for every grain:
- its maximum projected cross-section reproduces the measured outline (to ~1e‑6),
- its mesh volume matches the closed-form frustum volume (to ~1e‑15).

Tetrahedral (C3D4) meshing needs no extra dependencies: the outline is ear-clipped,
each triangle swept into a prism, each prism split into 3 tets using a globally
sorted vertex rule that keeps shared diagonals — and therefore the mesh — conforming.

**5. Wheel.** Wheel axis on Z. The rim is a structured C3D8 hex mesh of the annular
sector. Grain count comes from areal density or from abrasive concentration
(C100 = 25 vol%) combined with the measured mean grain volume. Positions come from
jittered-grid sampling whose cell size *guarantees* the minimum separation, so
grains cannot interpenetrate. Grains are seated tip-outward with a bounded random
tilt, at a protrusion drawn from a seeded truncated normal.

**6. Abaqus export.** Each distinct grain shape is one `*Part`; each grain on the
wheel is an `*Instance` of it. A 30° sector with 10,472 grains from 548 shapes is
549 parts — 5.9 MB instead of gigabytes. Node sets are emitted for the outer
surface, bore, axial faces and, for sectors, the two cut faces (meshed identically
so entry *k* pairs with entry *k*, ready for `*Equation` cyclic symmetry).

### Instance transform convention

`*Instance` gets one translation line and one rotation line. Abaqus applies the
**translation first**, then the rotation about the axis given in assembly
coordinates; the axis is written so it passes through the grain's own final centre,
leaving the centre fixed. This is the pattern Abaqus/CAE itself emits.

Because a wrong assumption would silently misplace every grain, `verify_inp_roundtrip`
re-reads the written file and checks physical invariants — each grain's protrusion
above the bond surface and its angular position. `--grain-parts baked` is an escape
hatch that writes pre-rotated coordinates and no instance transform, so placement
cannot depend on the convention at all.

## Cutting edge radius

`--edge-radius-um R` blunts every grain to a specified cutting edge radius. **Off by
default** (`0`), which leaves the geometry mathematically sharp.

Why it matters for a grinding simulation:

- A sharp edge has **zero radius**, which is a **stress singularity** in FEA. The
  computed stress grows without bound as the mesh is refined, so results never
  converge and depend on element size.
- The edge radius physically sets the **minimum chip thickness** and the
  ploughing-to-cutting transition, so it is a first-order parameter for grinding,
  not a cosmetic detail.

It is a **user parameter, not measured from the image**, deliberately: edge radius is
a dressing and wear property that changes during wheel life, and a top-down SEM view
of loose grit cannot resolve it reliably. A reasonable starting point is ~10% of the
measured `d50`.

Applied in two places, both verified against analytic ground truth:
- the outline's convex corners, by morphological opening (a square of side *s*
  blunted by *r* loses exactly (4−π)r² — reproduced to 0.04%);
- the meridional profile, by tangent circular fillets, so the ridges where the side
  surface meets the top and bottom caps are rounded too (arc points land on the
  requested radius to 7e-16).

Measured effect on the hardest image, per grain:

| | sharp | `--edge-radius-um 0.35 --arc-segments 1 --max-vertices 12` |
|---|---|---|
| faces | 128 | 210 |
| convex edges >60° | **44.1** | **5.9** |
| convex edges >30° | 78.7 | 72.4 |
| worst edge | 172° | 158° |
| silhouette area removed | — | 13.9% |

`--arc-segments 1` is a chamfer: cheapest, and enough to remove the singularity.
`2`–`3` give a progressively smoother radius at roughly 2× the faces. Lower
`--max-vertices` and raise `--simplify-um` to claw the element count back — that is
the simplification lever, and it is why the blunted grain above costs only 1.6× the
faces rather than 6×.

### Limits, stated plainly

- A handful of edges (~2–3 per grain) stay above 90°. They sit at the flat cap
  boundaries and are artifacts of the loft, not outline corners. For **R3D3 rigid**
  grains this is immaterial — a rigid body carries no stress.
- The achieved circumferential radius is smaller than requested (0.278 vs 0.350 µm
  above) because a short cap run clamps the fillet. Both achieved values are
  reported per run and stored per grain.
- On strongly concave outlines the taper is limited by the largest cross-section
  scale that keeps the polygon simple, so `--top-scale` can be ignored. This is
  reported as a warning.
- **Tet quality is poor for deformable grains.** The prism-based decomposition
  yields ~20% sliver elements even unblunted, and ~29% blunted. Prefer
  `--grain-element R3D3`, which uses the surface only. If you need C3D4, the
  per-run report gives min element volume and quality so you can judge the
  consequences for the stable time increment.

## What is measured vs modelled

Everything in-plane is measured. **Grain height cannot be obtained from a single
top-down SEM image** and is the one modelled quantity: thickness = ratio × minimum
Feret width, with the ratio drawn from a seeded truncated normal
(default mean 0.70, sd 0.12, clipped to [0.45, 0.95]) representing blocky crushed
grit. Override with `--thickness-ratio` / `--thickness-std`. Everything is seeded,
so a run is reproducible; the same seed produces a byte-identical deck.

Material properties in the deck are representative literature values, clearly
marked. Confirm them against your own data before drawing conclusions.

## Ductile / brittle single-abrasive model

`vumat_grind.for` carries both constitutive laws and picks between them **per material
point**, from the undeformed chip thickness the grit takes at that point's station along
the scratch:

```
h <  dc  ->  Johnson-Cook + strain-gradient enhancement   (ductile: plastic flow)
h >= dc  ->  Johnson-Holmquist II                         (brittle: damage, chipping)
```

Below a critical depth of cut a brittle solid is removed by plastic flow, not fracture.
A pure JH-2 deck cannot show that, so it reports brittle chipping at every depth.

**The critical depth.** Two published forms, both provided, selected by `dc_form`:

| | |
|---|---|
| `dc = λc (H/E)^½ (Kc/H)²` | form 1 |
| `dc = λc (E/H) (Kc/H)²` | form 2, Bifano, Dow & Scattergood (1991), λc = 0.15 |

They differ by `(E/H)^1.5` — about 17× on this sandstone, giving 5.3 nm against 88 nm —
so **λc belongs to one form and is not transferable to the other**. Say which you used.
Toughness goes in as MPa·√mm; `hybrid.kic_from_mpa_sqrt_m` does the ×31.6.

**How the subroutine learns `h`.** A VUMAT sees no kinematics. With one grit the
trajectory is exact, so `h` is a closed-form function of the tangential station `u`:

```
h(u) = H0 + HG·u - u²/(2·RTIP)        HG = -v_r / (ω·RTIP)
```

The linear term is the rubbing → ploughing → shearing wedge, produced here by the radial
infeed rather than a table feed; the quadratic term is the sagitta of the grit's circular
path, 15 nm across a 48 µm block, which matters only because `dc` is of the same order.
`semgrit/hybrid.py` computes `H0, HG, RTIP` and writes them into the material card, so
the Fortran carries no process knowledge and stays verifiable on a single material point.
`H0` is pinned to the deck's own tangency rather than re-derived, so no disagreement
about which tip is tallest can reach a nanometre-scale threshold.

**Strain-gradient enhancement**, from Yadav et al. (IJMS 2022 eq. 8, IJMTM 2024 eq. 25,
IJMS 2026 eq. 7 — the same Taylor/GND hardening with different characteristic lengths):

```
σe = σ_JC · sqrt(1 + (r'·η·b·(M α G)² / σ_JC²)^Λ)      η = 4 εp / h
```

Thin cuts raise the strain gradient, which raises the flow stress. That size effect is
*why* thin cuts are ductile, so the JC branch would understate the ductile regime
without it. The gradient length is floored at the Burgers vector, or `h → 0` in the
rubbing zone returns an infinite flow stress.

**Damage degrades the surface, not the stress.** Scaling the stress tensor by `(1-D)`
each increment — which `vumat_jc_damage (1).for` does — compounds to `(1-D)^k`, and the
per-increment loss soon cancels the elastic increment exactly: the point parks under
yield and never fails. Measured on the material-point driver, plastic strain froze at
3.1e-4 and `D` at 0.0017 for the remaining 39,000 increments. Both branches therefore
move their strength surface with damage, as JH-2 always has.

**State variables** (20, `*Depvar, delete=12`). 1–12 are JH-2's, unchanged, so existing
post-processing still works. Then **13** branch (1 ductile, 2 brittle), **14** `h`,
**15** `dc`, **16** temperature, **17** σ_JC, **18** σe, **19** the SGE amplification,
**20** the init latch. Plot SDV13 to see where the transition sits.

```bash
abaqus job=grind input=single_abrasive_hybrid.inp user=vumat_grind.for double=both cpus=8
```

`double=both` is required: `h` is compared against a few nanometres on a 25 mm radius.

Build one from Python, or use notebook section **B**:

```python
from semgrit.build_deck import build_deck, hybrid_single_grit
from semgrit.hybrid import HybridParams, kic_from_mpa_sqrt_m
hp = HybridParams(enabled=True, kic=kic_from_mpa_sqrt_m(0.3))
build_deck(hybrid_single_grit(hybrid=hp), solids, "out")
```

`h_source` = 2 or 3 forces ductile or brittle everywhere on an otherwise identical deck,
which is how you find out how much of a result the switch caused.

**Limits, stated plainly.** The Johnson-Cook constants shipped for sandstone are
placeholders and the deck says so in its own header — `A` is tied to that material's JH-2
quasi-static compressive strength so the two branches meet at the transition, but `B, n,
C, m` and `D1..D5` are order-of-magnitude values. `h` is prescribed from one grit's
trajectory, not measured, so it is not valid for many interacting grits or for grits that
wear during the run. The branch is latched at the first increment and does not migrate.

## Fixes relative to the original notebook

| Issue | Original | Now |
|---|---|---|
| **Scale calibration** | Bright-object search locked onto the white databar strip (1019 px instead of the 68 px bar) → **every measurement 14.9× too small at 10 kX, 29.9× at 5 kX** | Read from TIFF metadata; scale bar cross-checks to +0.4% |
| Scale bar label | Hardcoded `2.0 µm`; wrong for `B4C_16` and `DIAMOND_14`, which show 1 µm | Recovered by 1‑2‑5 snapping |
| Databar removal | Canny row-density heuristic | Row grey-level-count test (databar rows have ≤14 levels, micrograph rows ≥91) |
| Data destruction | Blanked a fixed rectangle, deleting ~5% of real grain area | Removed |
| Border grains | Truncated grains measured as whole (38–67% of regions) | Flagged and excluded from distributions |
| Watershed | Run on a *binary* Canny image — all plateaus, so boundaries fell arbitrarily | Continuous distance + gradient elevation; boundary/interior gradient ratio improved from ~1.7 to ~7 |
| Thresholds | `min_distance=30 px`, `area<800 px` — magnification-dependent | All in microns |
| Dead code | Cells 15–22 computed then discarded; cells 31–33 duplicated | Removed |
| Grain shape | Convex hull simplified to 8–10 vertices, inflating size and erasing the sharp concave cutting features | Real outline, tolerance in microns |
| Grain height | `uniform(0.35, 0.55) × width`, unseeded | Documented seeded model on measured Feret width |
| Mesh validity | Random vertex jitter → self-intersecting, non-watertight, no checks | 548/548 solids verified watertight, positive volume, exact analytic volume |
| Rotations | About the **world origin**, flinging grains far from their measured positions | About the grain centroid |
| Grain overlap | None — grains freely interpenetrated | Guaranteed by construction; audited by KD-tree |
| Wheel output | `concatenate` into one STL (not a valid solid; grains fused into the body) | Abaqus parts/instances, grains as separate bodies for contact/tie |
| CAD | STL only, non-watertight, unusable in CAD | STEP faceted B-rep, verified closed and outward-oriented, imports to SOLIDWORKS as solid bodies |
| Sector | Not possible | Any angle, with cut-face node sets for cyclic symmetry |
| Reproducibility | No seeds anywhere | Fully seeded; identical seed → identical bytes |
| Batch | Colab upload + `input()` prompts, one image at a time | CLI over all 14 images in ~29 s |

## Verification

`python verify_all.py` runs 31 checks: unit tests against analytic ground truth
(circle/ellipse/square/L-shape measurements, prism-to-tet partition by Monte Carlo,
mesh conformity, axis-angle round trip, sampler spacing, hex volumes vs analytic for
six sector angles, STEP write/read volume round-trip, hex-boundary closure),
integration across all images, wheel assembly at five sector/element/part
combinations, STEP and STL exports checked against the FE model, plus Abaqus
round-trips and determinism.

Two checks are deliberately negative controls, so the suite cannot pass vacuously:
an inward-wound STEP solid must be rejected by the audit, and a saturated packing must
reject colliding grains rather than let them interpenetrate.

Perimeter uses the Crofton estimator (4 directions), benchmarked against squares
rotated 0–45°: mean −0.74%, RMS 2.72%, worst 5.14% — versus 5.15% RMS for a traced
contour and 10.68% for Crofton with 2 directions. Circularity therefore carries
roughly ±10% uncertainty on a single grain.
