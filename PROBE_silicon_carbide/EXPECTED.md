# What to accept

Material: **silicon_carbide**, dc = **52.9200 nm**, card copied from `multi_abrasive_field.inp`.

Every number below was produced by running **this card with these
field values** through the compiled VUMAT outside Abaqus, so the
table is a verified prediction rather than a hope. Any disagreement
is therefore Abaqus' field-variable plumbing and nothing else.

| element | h/dc | h [mm] | SDV14 must be | SDV13 must be | SDV19 |
|---|---|---|---|---|---|
| 1 | 0.10 | 5.292000e-06 | 5.292000e-06 | 1 (ductile) | 5.3960 |
| 2 | 0.50 | 2.646000e-05 | 2.646000e-05 | 1 (ductile) | 2.7303 |
| 3 | 0.90 | 4.762800e-05 | 4.762800e-05 | 1 (ductile) | 2.1445 |
| 4 | 0.99 | 5.239080e-05 | 5.239080e-05 | 1 (ductile) | 2.0668 |
| 5 | 1.01 | 5.344920e-05 | 5.344920e-05 | 2 (brittle) | 0.0000 |
| 6 | 1.10 | 5.821200e-05 | 5.821200e-05 | 2 (brittle) | 0.0000 |
| 7 | 2.00 | 1.058400e-04 | 1.058400e-04 | 2 (brittle) | 0.0000 |
| 8 | 5.00 | 2.646000e-04 | 2.646000e-04 | 2 (brittle) | 0.0000 |

SDV19 is the strain-gradient amplification. It is above 1 in the ductile
branch and rises as h falls (5.3960 at h = 0.10 dc, 2.0668 at h = 0.99 dc);
it is exactly 0 in the brittle branch, which never evaluates it.

## The three outcomes

**All eight SDV14 equal the injected h, and SDV13 flips between
element 4 and element 5.** The field route works. Every
field-carrying deck in the project is a genuine hybrid deck, and
the roadmap's Phase 2 onwards is aimed at code that executes.

**Every SDV14 is zero.** Field variable 1 is not arriving. Then
`hloc` is 0, `0 < dc`, and `RUN_ME*/2_multi_abrasive` and
`3_energy_criterion` are running **100% ductile** -- with a
plausible chip, a clean `.sta` and nothing in any output to say
so. Fix that before spending 5.3 h per material on anything.

**The job dies at preprocessing.** Read the `.dat`. Compare
against `error/single_abrasive.dat`, which is the same failure
this project already had once (`*User Material` four values to a
line instead of eight).

## Also worth reading off the same job

* `STATUS` -- if nothing ever deletes, the three-way agreement
  between `*Depvar, delete=12`, `ELEMENT DELETION=YES` and the
  VUMAT zeroing SDV12 is broken, and no chip will ever separate
  in the real decks either.
* `SDV13` on elements 4 and 5 (h = 0.99 dc and 1.01 dc) -- an
  off-by-one or a `<=` where a `<` belongs shows up only there.
* `SDV19`, the strain-gradient amplification. At h = 0.1 dc it
  should be well above 1; at h = 5 dc close to 1.
* `ALLAE/ALLIE` -- hourglass energy on a single-element test
  should be negligible. If it is not, the enhanced hourglass
  control is not doing its job at this aspect ratio.
