# RUN_ME - three submittable Abaqus packages

**Workpiece material: Sandstone (JH-2, the project's original card)**

| | |
|---|---|
| critical depth of cut `dc` | 87.7500 nm (form 2) |
| depth of cut `ae` | 0.200 um = 2.3 x dc |
| JH-2 `K1, G, HEL, PHEL, T` | 3735.6, 2686, 1982, 1374, 8 MPa |
| density | 2350 kg/m3 |
| hardness `H`, toughness `Kc` | 1000 MPa, 0.3 MPa*sqrt(m) |
| estimated wall clock, 8 cores | 9.27 h per package |

> E = 6500 MPa and nu = 0.21 are exactly this card's K1 and G, so both
> branches share one elasticity and the stable increment is
> unambiguous. A = 90 MPa is the card's own quasi-static uniaxial
> compressive strength. B, n, C, m and D1..D5 are PLACEHOLDERS.

Same wheel, same grain library and same mesh throughout. Packages 2, 3
and the three arms in 4 also share the same **seating** and the same
swept field, so a difference between those is the criterion or the law
and nothing else. That is the comparison to quote.

**Package 1 is not directly comparable to them.** It is one grit, not
twelve, so its seating, its ground radius and its peak chip thickness
all differ -- `MANIFEST.json` in this folder carries the numbers. Treat
it as a smoke test, and as the closed-form check the swept field was
validated against, rather than as an arm of the comparison.

| folder | what | subroutine |
|---|---|---|
| `1_single_abrasive` | chip thickness from the four wedge constants in the card | `vumat_grind.for` |
| `2_multi_abrasive` | chip thickness swept per element, carried as field variable 1 | `vumat_grind.for` |
| `3_energy_criterion` | geometric split plus the local W_p L_c >= PSI Kc^2/E rule | `vumat_grind2.for` |
| `4_ablation` | The same deck as 2_multi_abrasive with PROPS(56) changed and nothing else. forced_brittle is bit-identical to vumat_jh2.for, so the hybrid arm sits bracketed between two known references instead of standing alone. | (three arms) |

## Running them

Each folder has `run.bat` (Windows) and `run.sh`. Or by hand:

```
cd 1_single_abrasive
abaqus job=single_abrasive input=single_abrasive.inp user=vumat_grind.for double=both cpus=8 interactive

cd 2_multi_abrasive
abaqus job=multi_abrasive input=multi_abrasive_field.inp user=vumat_grind.for double=both cpus=8 interactive

cd 3_energy_criterion
abaqus job=energy_criterion input=energy_criterion.inp user=vumat_grind2.for double=both cpus=8 interactive

cd 4_ablation
cd 3_forced_brittle && abaqus job=forced_brittle input=forced_brittle.inp user=vumat_grind.for double=both cpus=8 interactive && cd ..
cd 2_forced_ductile && abaqus job=forced_ductile input=forced_ductile.inp user=vumat_grind.for double=both cpus=8 interactive && cd ..
cd 1_geometric_hybrid && abaqus job=geometric_hybrid input=geometric_hybrid.inp user=vumat_grind.for double=both cpus=8 interactive && cd ..

```

`double=both` is **required**, not a preference: the chip thickness is
compared against a threshold of a few nanometres on a 25 mm radius, a
ratio of 1e-7, and single precision does not have the digits.

The `.for` file must keep its name. A filename with a space or a
bracket makes Abaqus read part of it as a separate argument and abort.

## What to plot afterwards

| SDV | |
|---|---|
| **13** | the branch: 1 ductile, 2 brittle. This is the picture. |
| **14** | the chip thickness that point was given |
| **15** | dc |
| **19** | the strain-gradient amplification, 1 = no size effect |
| **21, 22** | plastic work and the energy ratio (package 3 only); SDV22 reaching 1 is what flips a point |
| 1, 2, 12 | damage, equivalent plastic strain, STATUS |

`RF` and `RM` at `A_WHEEL_REF` are the grinding force. The
`*_postprocess_odb.py` in each folder reads exactly those:

```
abaqus python <name>_postprocess_odb.py <job>.odb
```

That writes the CSVs and the summary JSON. It also tries to draw the
figures, but **Abaqus' bundled Python usually has no matplotlib**, so in
practice it writes data and no pictures -- which is why the first six runs
of this project ended up documented by photographs of a screen. Draw them
with the host Python instead, which needs no Abaqus and no re-run:

```
python REPOST/plots.py <folder holding the CSVs>
```

It also writes `compare_all.png` across every job it finds, and marks on
the figures anything that should stop a number being quoted -- artificial
energy over 5% of internal, and byte-identical duplicate datasets.

## What the three should show

Package 1 puts the transition at one station along one scratch, from a
wedge written into the card. Package 2 computes the chip thickness for
every element from the real grit trajectories, so several grits and
their overlapping grooves are handled properly. Package 3 reaches the
transition from the material point's own plastic work instead of from
geometry at all. **If 2 and 3 disagree, that disagreement is a result**
-- it is the difference between a geometric and an energetic reading of
the same transition.

## Read this before quoting a number

* The **Johnson-Cook constants are placeholders**. `A` is tied to the
  JH-2 card's own quasi-static compressive strength so the two branches
  meet at the transition, but `B, n, C, m` and `D1..D5` are
  order-of-magnitude values. The decks say so in their own headers.
* `lambda_c` belongs to whichever `dc` form was used. The two published
  forms differ by `(E/H)^1.5`, about 17x on this rock, and the energy
  criterion is a third member of the family again.
* `PSI` in package 3 is **mesh-dependent** by construction. Quote the
  element size with it.
* Check `SDV13` against the split the build reported. If they differ,
  the field did not reach the material points.
