# wheel_rigid_2mm.inp — geometry only, all rigid except the workpiece

Two parts, two instances, nothing else. No step, no interaction, no BC, no output —
those are yours.

## Load it

Abaqus/CAE → **File → Run Script…** → `wheel_rigid_2mm_import_into_cae.py`

Not File → Import → Part (that reads the `*Part` blocks and skips the `*Assembly`, so
the wheel arrives bare) and not File → Import → Model (it does not accept `.inp`).

## What is in it

| | |
|---|---|
| `WHEEL` | **one** discrete rigid body — bond rim shell (2,812 R3D4) **plus all 712 grits** (79,844 R3D3), single reference node **44161** at the origin on the axis |
| `WORKPIECE` | the only deformable part — 160,000 C3D8R, 0.3 µm elements, material `STONE` |

Wheel: R = 25 mm (D50), **2.000 mm arc** (4.5837°), rim depth 12 µm, width 30 µm.
Sagitta 20 µm = **167 % of the rim depth**, so it reads as an arc, not a rectangle.

Workpiece: 48 × 15 × 6 µm, ground face at r = 25.001860 mm — **exactly tangent** to the
tallest grit inside its footprint, zero initial penetration.

## Driving it

One BC on `A_WHEEL_REF` moves the whole wheel — that is the point of the single rigid
body. For 30 m/s surface speed:

```
VR3 = -1200 rad/s        (= 11,459 rpm; sign sets the cutting direction)
V1 = V2 = V3 = VR1 = VR2 = 0
```

The grits currently *graze* the surface. To cut at depth `ae` (mm), add the radial
infeed to the same BC — the workpiece sits at θ = 2.291831°, so radially inward is

```
V1 = -0.999200 * ae / t_step
V2 = -0.039989 * ae / t_step
```

**Keep `ae` below 1.86 µm.** That is the clearance between the bond rim and the ground
face; past it the bond itself starts hitting the workpiece instead of the grits.
1.0 µm is a sensible single-grit depth of cut here.

Fix the workpiece with `A_WP_BACK_FACE` (encastre) and `A_WP_SIDE_A` / `A_WP_SIDE_B` /
`A_WP_END_A` / `A_WP_END_B` as you prefer.

## Surfaces provided

| name | use |
|---|---|
| `A_GRITS_SURF` | all 712 grits only |
| `A_GRITS_ENGAGE_SURF` | the **56** grits that can reach the block during one pass — use this if contact cost bites |
| `A_WHEEL_SURF` | grits + bond outer face |
| `A_WP_GROUND_SURF` | the workpiece face being ground |

General contact is fine; contact inside a rigid body is ignored, so the grit bases
sitting inside the rim shell cost nothing.

## Step and run time

```
dt (no mass scaling) = 6.305e-11 s      travel 0.054 mm at 30 m/s -> t_step = 1.8e-6 s
increments           = 28,549           element-increments = 4.57e9
```

Estimate **2.4 h on 4 cores, 1.2 h on 8** — assuming 3e5 element-increments/s/core for
C3D8R under a scalar VUMAT, 70 % parallel efficiency, and contact costing a further
60 %. That is a deliberately pessimistic rate; anything faster only shortens the run.
I have no timing data from your machine, so treat it as an estimate, not a measurement.

Wall clock is **linear in step time**. If it finishes early and you want more cut, just
raise `t_step` — up to about 6.6e-5 s before the wheel runs off the 2 mm arc. Going
finer than 0.3 µm elements costs `(0.3/h)^4`, so 0.25 µm would roughly double it.

Mass scaling: not recommended here. Brittle fracture under JH-2 is rate-sensitive, and
the run is short enough not to need it.

## The workpiece material

`STONE` is a placeholder `*Elastic` (50 GPa, ν = 0.25, 2650 kg/m³) so the
`*Solid Section` resolves. Replace it with your JH-2 block:

```
*Material, name=STONE
*Density
 2.65e-09,
*User Material, constants=17
 ... your 17 JH-2 constants ...
*Depvar, delete=12
 12
```

`delete=12` matters — see the note on `vumat_jh2_3 (2).for` vs `vumat_jh2.for`: the
version that shows brittle fracture drives SDV12 to zero at D = 1, and that only
removes elements if the deletion flag is declared here *and* element deletion is on in
Section Controls.

## Verification

```
python verify_rigid_deck.py  FINAL_RIGID/wheel_rigid_2mm.inp    # 54 checks, 0 failures
python verify_rigid_deck2.py FINAL_RIGID/wheel_rigid_2mm.inp    # 29 checks, 0 failures
```

Two independent verifiers, no shared code. The first re-derives the geometry from the
node coordinates (arc, outward normals on all 2,812 shell quads, tangency by facet
clipping, grit-grit interpenetration by point-in-polyhedron). The second walks the file
as a keyword state machine (grammar, data-line arity, context legality, ordering) and
re-derives the mass and inertia by numerical integration — both matched the closed form
to 0.000 %.
