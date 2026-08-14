# wheel_single_grit.inp — one big grit, everything else identical to FINAL_RIGID

Same D50 wheel, same 2.000 mm arc, same 12 µm rim, same 30 µm width, same
48 × 15 × 6 µm workpiece at 0.3 µm, same 30 m/s, same 54 µm pass, still one discrete
rigid body with the workpiece as the only deformable part. **Only the grit population
changed: 1 instead of 712.**

## Load it

Abaqus/CAE → **File → Run Script…** → `wheel_single_grit_import_into_cae.py`

## The grit

Library grain **0**, the largest of the 96:

| | |
|---|---|
| volume | 193.9 µm³ |
| height | 6.54 µm |
| across | 10.35 µm (aspect 1.05 — chunky and equiaxed, like your screenshot) |
| facets | 116 R3D3 |
| protrusion above the bond | **3.769 µm** |
| seated | 11.15° off the outward radial, so tip-outward |
| position | θ = 2.326209°, i.e. **b = +15 µm** from the block centre |

Seating went through the same code path as the 712-grit deck — tip outward, random spin,
tilt within 35°, protrusion from the same truncated normal, then the radius solved so
the furthest vertex sits exactly that far above the *curved* rim.

## ⚠ Rotation direction now matters

With 712 grits it didn't — grits covered the whole arc. With one it decides whether the
grit cuts or walks away. The grit starts at the trailing end of the block:

```
VR3 = -1200 rad/s   ->  grit is dragged 39 µm across the 48 µm block and exits  ✓
VR3 = +1200 rad/s   ->  grit leaves the block after 9 µm
```

Full BC on `A_WHEEL_REF` (ref node **2875**, origin, on the axis):

```
VR3 = -1200 rad/s              (= 11,459 rpm, 30 m/s at r = 25)
V1 = V2 = V3 = VR1 = VR2 = 0
```

To cut at depth `ae` (mm), add the radial infeed. The workpiece sits at θ = 2.291831°,
so radially inward is

```
V1 = -0.999200 * ae / t_step
V2 = -0.039989 * ae / t_step
```

**Keep `ae` below 3.76 µm** — that is the bond-rim-to-workpiece clearance, and it is
twice what the 712-grit deck allowed because this grain protrudes further. 1.0–1.5 µm is
a sensible single-grit depth of cut.

Fix the workpiece with `A_WP_BACK_FACE` plus whichever of `A_WP_SIDE_A` / `A_WP_SIDE_B` /
`A_WP_END_A` / `A_WP_END_B` you want.

## Contents

| | |
|---|---|
| `WHEEL` | one rigid body — 2,812 R3D4 bond shell + **116 R3D3** grit facets, ref node 2875 |
| `WORKPIECE` | 160,000 C3D8R, 0.3 µm, material `STONE` (placeholder `*Elastic` — swap in your JH-2 VUMAT) |

Ground face at r = 25.003764 mm, exactly tangent to the grit apex, penetration 0.

Surfaces: `A_GRITS_SURF` (the grit), `A_WHEEL_SURF` (grit + bond outer face),
`A_WP_GROUND_SURF`. There is no `A_GRITS_ENGAGE_SURF` in this deck — with one grit it
would be identical to `A_GRITS_SURF`.

## Run time

```
dt = 6.305e-11 s     t_step = 1.8e-6 s     28,549 increments     4.57e9 element-increments
```

The same conservative formula as the 712-grit deck gives **2.4 h on 4 cores, 1.2 h on 8**.
It should land well under that here: contact has 116 rigid facets to track instead of
82,656, and the formula charges a flat 60 % for contact regardless. The estimate is kept
identical so the two decks are directly comparable.

If you want it shorter, `t_step` is linear and 45 µm of travel (`t_step = 1.5e-6 s`) is
enough for the grit to clear the block.

## Verification

```
python verify_rigid_deck.py  SINGLE_GRIT/wheel_single_grit.inp    # 54 checks, 0 failures
python verify_rigid_deck2.py SINGLE_GRIT/wheel_single_grit.inp    # 29 checks, 0 failures
```

Both verifiers, unchanged from the 712-grit deck, pass on this one too. Notable results:
all 2,812 shell quad normals outward; grit closed with outward normals; tangency gap
0.001 nm with zero penetration; all 60 grit vertices inside the block footprint;
mass and inertia matching numerical integration to 0.000 %.
