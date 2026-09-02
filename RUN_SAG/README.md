# RUN_SAG — decks that mimic the reference paper

Reproduces the experiment in

> Ghosh, G., Sidpara, A. & Bandyopadhyay, P.P. (2021). *Brittle-ductile
> transition in compliant finishing of HVOF sprayed hard WC-Co coating.*
> Int. J. Refractory Metals and Hard Materials **99**, 105610.

Three decks, one per abrasive pad, at the paper's best operating point
(T = 0.4 mm, N = 1050 rpm) on HVOF-sprayed WC-12Co.

| pad | elements | passes | 8-core estimate | paper observed |
|---|---|---|---|---|
| `30um/micro_30um.inp` | 107,604 | 11 | **~96 h** | plastic flow **+ brittle fracture** |
| `15um/micro_15um.inp` | 191,296 | 16 | ~162 h | plastic flow, **fewer** fractures |
| `6um/micro_6um.inp` | 483,181 | 24 | ~392 h | **pure ductile, no fracture** |

**Run the 30 µm deck first.** It is the cheapest, and it is the one that should
show brittle fracture — so it is the fastest route to knowing whether the model
reproduces the paper at all.

## Running

```
cd 30um
abaqus job=micro_30um input=micro_30um.inp user=../../vumat_grind2.for double=both cpus=8 interactive
```

`double=both` is mandatory. `h` and `dc` are compared at 80 nm against a
millimetre-scale geometry; single precision has ~7 digits and does not have
them. The failure is silent — the branch flag comes out wrong, the job does not
crash.

Note the subroutine is **`vumat_grind2.for`**, not `vumat_grind.for`. These
decks carry 58 constants and use the local energy criterion; `vumat_grind.for`
reads 56 and would misread the card.

## What to plot

**SDV13 is the result.** 1 = ductile, 2 = brittle.

Plot it **after every pass**, not only at the end. The energy criterion
accumulates, so the interesting output is *the pass at which a point flips* —
that is the transition, and it is what distinguishes the three pads.

| SDV | meaning |
|---|---|
| **13** | **branch: 1 ductile, 2 brittle** |
| 15 | dc actually used, m |
| 2 | equivalent plastic strain |
| 12 | deletion flag |
| 21, 22 | the energy criterion's own accumulators (grind2 only) |

## What would confirm the model, and what would falsify it

**Confirmation.** The 30 µm deck shows brittle regions, the 15 µm deck shows
fewer, and the 6 µm deck stays entirely ductile. That ordering is the paper's
central observation, and it is what these decks are built to test.

**Falsification.** Any brittle element in the 6 µm deck is a real disagreement:
the paper reports pure ductile removal for that pad, and the model predicts
h/dc = 0.002 there. If the 30 µm deck stays entirely ductile through all 11
passes, the energy criterion accumulates too slowly — which is a finding about
the criterion, not a bug in the deck.

## Where the numbers come from

Everything below is from the paper except one quantity, which is calibrated
against the paper's own measurements.

**From the paper**: wheel 125 mm × 10 mm; T and N from Table 3; WC-Co
E = 200 GPa, H = 11.02 GPa, Kc = 7.78 MPa·√m, BHN = 581 (Table 2);
ν_w = 0.25, ν_t = 0.24; mean carbide 1.36 µm; pad densities 1750 / 720 /
312 mm⁻² with active = 0.8 × areal; **dc = 60–100 nm, measured** (section 4.2).

**Calibrated — the one free parameter.** The paper gives eq. (4) for the
backing pad's modulus from its shore hardness but never prints the shore
hardness, so E_t has to come from somewhere. Two independent sources:

- the user's hand-built CAE deck carries neo-Hookean C10 = 0.0575 MPa,
  i.e. E = 6·C10 = **0.345 MPa**;
- inverting the contact chain for the modulus that reproduces the paper's
  stated per-grain forces (1e-4…1e-5 N, section 4.1) gives **0.43 MPa**.

Those agree to 25 % — a hand-built card and a published measurement arrived at
separately. Pinning it tighter uses the paper's strongest statement, its
headline result: the 6 µm pad is pure ductile with 60–100 nm chips. Requiring
the groove width in that band *and* the 30 µm per-grain force ≤ 1e-4 N leaves
shore 25.6–27.9, so **C10 = 0.16606 MPa (E_t = 0.996 MPa)**. Only that value
satisfies both constraints; the 0.345 MPa card satisfies one. Run
`python _make_sag_paper.py --compare` to see the comparison, or
`--user-c10` to build with the CAE card instead.

## dc is MEASURED, not Bifano's

Bifano's expression on this material's own E, H and Kc gives **1.357 µm** — 17×
the measured 60–100 nm. That is the paper's own conclusion, and this project
reproduces it to 0.96 % through an independent implementation. The decks use
the **measured 80 nm**; `verify_sag_deck.py` recomputes Bifano from the card
and fails a deck that quietly fell back on it.

The energy criterion inherits the measurement: with PSI = 0 the subroutine
derives `dc·E·H/Kc²`, making the threshold exactly

    W_p · L_c  >=  H · dc  =  0.8816 MPa·mm  =  881.6 J/m²

## Repeated passes, and why

The criterion accumulates plastic work **per point**. A grain sliding once
along fresh material leaves every point with exactly one pass, so it can never
trip the threshold however far it goes — this was a real flaw in the first
version of these decks and it is why they now make repeated passes over one
track, reversing each time.

The pass count is **derived, not chosen**: from each grain's own tangential
work per pass against `H·dc`, with a 1.5× margin. That gives 11 / 16 / 24
passes for the 30 / 15 / 6 µm pads — coarser pads need fewer, which is the
physics.

The paper's 10 s spot test puts roughly **20,000** grain crossings over each
point. So these decks test whether work accumulates in the right *direction*
and at the right *rate*, not whether a single pass fractures.

## Caveats to state before quoting a number

- **The Johnson-Cook constants for WC-Co are placeholders** except `A`, which
  is derived from the JH-2 card's own quasi-static strength. `B, n, C, m` and
  `D1..D5` are order-of-magnitude values. **SDV13 is quotable; force
  magnitudes are not** until those are calibrated.
- **The JH-2 card is derived, not published.** K1 and G are exactly the
  paper's E = 200 GPa at ν = 0.22; HEL = 0.6H is the usual cemented-carbide
  ratio. A, B, C, N, M, D1, D2 are placeholders.
- **The energy criterion is mesh-dependent by construction.** It is
  regularised by the element length, so PSI is calibrated *for a mesh*. These
  decks use 16 nm through the depth (dc/5) and say so in their own headers.
  `verify_sag_deck.py` reports the sensitivity: the required work density
  scales 4× over a 4× refinement while the triggering *energy* stays fixed to
  1e-16.
- **The elliptical contact does not fit on the wheel.** Combining the paper's
  eqs. (6) and (7) forces a contact width of 10.2–11.1 mm on a 10 mm face, so
  the patch is clipped at every operating point. The clipped *area* is only
  0.2–3.7 %, so this is a caveat on the width rather than a refutation.

## Rebuilding

The `.inp` files are gitignored (61 MB for the 6 µm deck).

```
python _make_sag_paper.py --all            # all three pads
python _make_sag_paper.py                  # just the 6 um headline case
python _make_sag_paper.py --all --macro    # also the MACRO contact decks
python verify_sag_deck.py RUN_SAG/30um/micro_30um.inp
```

`SUMMARY.json` carries every derived number for each deck.
