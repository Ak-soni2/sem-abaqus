"""Turn the RUN_SAG decks into self-contained, copyable run folders.

    python _make_sag_packages.py            # all three pads
    python _make_sag_packages.py --only 30  # just one

Each folder ends up with everything a run needs and nothing that has to be
fetched from elsewhere:

    micro_<pad>.inp          the deck
    vumat_grind2.for         the subroutine, COPIED IN (not referenced)
    run.bat / run.sh         submit, with a datacheck gate first
    postprocess_odb.py       reads the .odb with Abaqus' own Python
    README.md                what to expect, and what would falsify it
    EXPECTED.md              the prediction, written BEFORE the run

The subroutine is copied rather than referenced by a relative path because the
folder is meant to be moved to a work directory, and ``user=../../x.for``
breaks the moment it is. Copying costs 100 kB and removes the failure.

EXPECTED.md exists so the prediction is on record before any result is seen.
A model that is only interpreted after the fact can accommodate anything.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "RUN_SAG")
SUB = "vumat_grind2.for"

RUN_BAT = r"""@echo off
cd /d "%~dp0"
rem  {pad} um pad -- SAG on HVOF WC-Co, after Ghosh et al. 2021.
rem
rem  Three stages, in order, and the first two are cheap:
rem
rem    0. abaqus verify -user_exp
rem       Can Abaqus build a user subroutine on THIS machine at all? If the
rem       Fortran toolchain is not wired up, everything after fails with a
rem       message that does not say so.
rem    1. datacheck -- seconds. Reads every keyword and the material card.
rem       The one real submission this project ever made died here, on
rem       *User Material written four values to a line instead of eight.
rem    2. the solve, ~{hours} h on 8 cores.
rem
rem  The datacheck submits under its OWN job name, {job}_check. A datacheck
rem  is a real submission: it writes <job>.lck and LEAVES it, because the
rem  job never completed. A solve reusing that name then aborts with
rem  "Detected lock file" on a job that has never been run. Renaming means
rem  the conflicting lock is never created -- better than deleting locks,
rem  since a lock is sometimes real and blindly removing one would clobber
rem  a live job.
rem
rem  double=both is REQUIRED. h and dc are compared at 80 nm against a
rem  millimetre geometry; single precision has ~7 digits and does not have
rem  them. The failure is SILENT -- the branch flag comes out wrong and the
rem  job does not crash.
rem
rem  The subroutine is vumat_grind2.for, NOT vumat_grind.for: this deck
rem  carries 58 constants and the local energy criterion. vumat_grind.for
rem  reads 56 and would misinterpret the card.
rem  Every abaqus line is `call`ed. On Windows `abaqus` is abaqus.bat,
rem  and running one batch file from another WITHOUT `call` transfers
rem  control permanently -- the caller never resumes. Without it this
rem  script ran the verify and then silently stopped, submitting
rem  nothing, with no error to say so.
call abaqus verify -user_exp
if errorlevel 1 (
  echo.
  echo  Abaqus cannot build a user subroutine on this machine.
  echo  Check that the Fortran compiler is on PATH and linked to Abaqus.
  exit /b 1
)
call abaqus job={job}_check input={job}.inp user={sub} double=both cpus=1 datacheck
if errorlevel 1 (
  echo.
  echo  DATACHECK FAILED -- read {job}_check.dat, the error is a keyword or the
  echo  material card, not the physics. Nothing has been solved yet.
  exit /b 1
)
echo.
echo  datacheck passed. Starting the solve -- about {hours} h on 8 cores.
echo.
call abaqus job={job} input={job}.inp user={sub} double=both cpus=8 interactive
if errorlevel 1 exit /b 1
call abaqus python postprocess_odb.py
"""

RUN_SH = r"""#!/bin/sh
# {pad} um pad -- SAG on HVOF WC-Co, after Ghosh et al. 2021.
# See run.bat for why each stage is here. double=both is required.
cd "$(dirname "$0")" || exit 1
abaqus verify -user_exp || {{
  echo "Abaqus cannot build a user subroutine on this machine." >&2
  exit 1
}}
abaqus job={job}_check input={job}.inp user={sub} double=both cpus=1 datacheck || {{
  echo "DATACHECK FAILED -- read {job}_check.dat. Nothing has been solved yet." >&2
  exit 1
}}
echo "datacheck passed. Solving -- about {hours} h on 8 cores."
abaqus job={job} input={job}.inp user={sub} double=both cpus=8 interactive || exit 1
abaqus python postprocess_odb.py
"""

POST = r'''"""Read the .odb with Abaqus' own Python and write the numbers out.

    abaqus python postprocess_odb.py

Runs INSIDE Abaqus, which has no matplotlib -- so this writes CSV and JSON
only, and plotting is left to the host Python. That is not a limitation being
worked around; it is the reason nothing in this project plots from inside
Abaqus.

The result is SDV13, the branch each material point took: 1 ductile, 2
brittle. What matters is not the final frame but HOW IT EVOLVES -- the pass at
which points start flipping. So the brittle fraction is written for every
frame of every step.
"""
import json
import os

from odbAccess import openOdb

JOB = "%(job)s"
DUCTILE, BRITTLE = 1, 2


def main():
    odb = openOdb(JOB + ".odb", readOnly=True)
    inst = odb.rootAssembly.instances["WORK-1"]
    n_el = len(inst.elements)

    rows = [("step", "frame", "time", "n_active", "n_ductile", "n_brittle",
             "brittle_fraction", "peak_mises", "peak_peeq")]
    per_step = {}
    for sname in odb.steps.keys():
        step = odb.steps[sname]
        for f, frame in enumerate(step.frames):
            fo = frame.fieldOutputs
            if "SDV13" not in fo:
                continue
            sdv = fo["SDV13"].getSubset(region=inst).values
            status = (fo["SDV12"].getSubset(region=inst).values
                      if "SDV12" in fo else None)
            nd = nb = na = 0
            for i, v in enumerate(sdv):
                if status is not None and status[i].data < 0.5:
                    continue
                na += 1
                if v.data >= 1.5:
                    nb += 1
                elif v.data >= 0.5:
                    nd += 1
            mis = 0.0
            if "S" in fo:
                for v in fo["S"].getSubset(region=inst).values:
                    if v.mises > mis:
                        mis = v.mises
            pk = 0.0
            if "PEEQ" in fo:
                for v in fo["PEEQ"].getSubset(region=inst).values:
                    if v.data > pk:
                        pk = v.data
            frac = (float(nb) / na) if na else 0.0
            rows.append((sname, f, frame.frameValue, na, nd, nb, frac,
                         mis, pk))
            per_step[sname] = frac
    odb.close()

    with open(JOB + "_sdv13.csv", "w") as fh:
        for r in rows:
            fh.write(",".join(str(x) for x in r) + "\n")

    passes = [k for k in per_step if k.startswith("PASS")]
    passes.sort(key=lambda s: int(s[4:]))
    first_brittle = None
    for k in passes:
        if per_step[k] > 0.0:
            first_brittle = int(k[4:])
            break
    out = dict(job=JOB, elements=n_el,
               brittle_fraction_by_step=per_step,
               passes=len(passes),
               first_pass_with_brittle=first_brittle,
               final_brittle_fraction=(per_step[passes[-1]] if passes
                                       else None))
    with open(JOB + "_summary.json", "w") as fh:
        json.dump(out, fh, indent=1)

    print("")
    print("  %%s: %%d elements" %% (JOB, n_el))
    print("  brittle fraction by pass:")
    for k in passes:
        print("    %%-8s %%.4f" %% (k, per_step[k]))
    if first_brittle is None:
        print("")
        print("  NO BRITTLE ELEMENTS in any pass.")
        print("  For the 6 um pad that AGREES with the paper (pure ductile).")
        print("  For the 30 um pad it does NOT: the paper sees fracture there,")
        print("  so the energy criterion is accumulating too slowly.")
    else:
        print("")
        print("  first brittle at PASS%%d of %%d" %% (first_brittle,
                                                    len(passes)))
    print("")
    print("  wrote %%s_sdv13.csv and %%s_summary.json" %% (JOB, JOB))


if __name__ == "__main__":
    main()
'''

EXPECTED = """# EXPECTED — written before the run

This is the prediction for the **{pad} µm** pad, recorded before any result is
seen. A model interpreted only after the fact can accommodate any outcome.

## The paper's observation for this pad

{observation}

## What this deck predicts

| quantity | predicted |
|---|---|
| per-grain normal force | {fn:.3e} N |
| groove width | {groove:.1f} nm |
| paper's measured chip | {chip_lo:.0f}–{chip_hi:.0f} nm |
| indentation depth | {indent:.4f} nm |
| h / dc | {hdc:.5f} |
| passes to reach H·dc | {passes} |
| **SDV13 outcome** | **{outcome}** |

## What would falsify this

{falsify}

## What CANNOT be checked against the paper

- **Force magnitudes.** The WC-Co Johnson-Cook constants `B, n, C, m` and
  `D1..D5` are placeholders — order-of-magnitude values, not measurements.
  Only `A` is derived, from the JH-2 card's own quasi-static strength. Until
  those are calibrated against nanoindentation or scratch data, no force this
  deck prints is quotable.
- **Surface roughness.** The paper's S_a = 21 nm is a statistical surface from
  roughly 20,000 grain crossings. This deck runs {passes} passes of one grain.
  There is no S_a to compare.
- **MRR.** The paper's 0.58 mm³/min comes from its analytical eq. (16) over a
  128 mm² spot. This deck removes material from one patch over microseconds.

**SDV13 — the branch map — is the result.** Everything else is diagnostic.
"""

PADS = {
    6.0: dict(
        observation=(
            "Section 4.2: *no brittle fracture is found in this case and the "
            "material is removed through pure ductile mode.* The ratio "
            "k = d_WC/d_g = 4.41 is below the paper's threshold of 5."),
        outcome="ENTIRELY DUCTILE — SDV13 = 1 everywhere, all passes",
        falsify=(
            "**Any brittle element.** The paper reports pure ductile removal "
            "for this pad and the model predicts h/dc = 0.002, so a single "
            "SDV13 = 2 is a real disagreement, not noise.")),
    15.0: dict(
        observation=(
            "Section 4.2: *a few traces of brittle fracture are identified... "
            "the amounts of brittle fracture are significantly lower than "
            "that of 30 um pad.* k = 11."),
        outcome="MOSTLY DUCTILE, with some brittle elements late in the passes",
        falsify=(
            "Either extreme falsifies it: entirely ductile through all passes, "
            "or a brittle fraction at or above the 30 um deck's. This pad "
            "should sit strictly between the other two.")),
    30.0: dict(
        observation=(
            "Section 4.2: *although the material is mostly removed through "
            "plastic flow, several traces of brittle fracture are observed.* "
            "k = 22, far above the threshold of 5."),
        outcome="BRITTLE REGIONS PRESENT — SDV13 = 2 somewhere",
        falsify=(
            "**Entirely ductile through all passes.** That would mean the "
            "energy criterion accumulates too slowly to reach H·dc under a "
            "realistic per-grain load — a finding about the criterion's rate, "
            "not a bug in the deck. It is the single most informative negative "
            "result available here, which is why this pad is the one to run "
            "first.")),
}

README = """# SAG {pad} µm pad — run folder

Self-contained. Copy the whole folder to your work directory and run
`run.bat` (Windows) or `run.sh` (Linux). Nothing is referenced outside it.

```
run.bat
```

That runs three stages: `abaqus verify -user_exp` (can this machine build
a subroutine at all), a **datacheck** (seconds — catches keyword and material
card errors before the queue), then the solve, then the postprocessor.

**Estimated {hours} h on 8 cores.**

## Read EXPECTED.md first

It records what this deck predicts, written before the run. The point of the
exercise is a falsifiable test, and that only works if the prediction is on
record beforehand.

## What to look at

`postprocess_odb.py` runs automatically and writes:

- `{job}_sdv13.csv` — brittle fraction for **every frame of every pass**
- `{job}_summary.json` — the pass at which brittle elements first appear

**SDV13 is the result**: 1 = ductile, 2 = brittle. Plot it in CAE after each
pass, not only at the end — the energy criterion accumulates, so *when* a point
flips is the physics.

## Requirements

- Abaqus with a working Fortran toolchain (`abaqus verify -user_exp`)
- `double=both` — non-negotiable, and the failure is silent if omitted
- `vumat_grind2.for` (58 constants, energy criterion) — **already in this
  folder**, not `vumat_grind.for`

## Honest scope

This deck tests **one** thing: whether SDV13 flips in the right order across
the three pads. It is not a numerical reproduction of the paper — force
magnitudes rest on placeholder Johnson-Cook constants, and surface roughness
needs orders of magnitude more grain passes than are simulated here.
See `EXPECTED.md` for the full list of what cannot be compared.
"""


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", type=float, default=0.0,
                    help="just one pad, e.g. --only 30")
    a = ap.parse_args(argv)

    summ = json.load(open(os.path.join(SRC, "SUMMARY.json"), encoding="utf-8"))
    sub_src = os.path.join(HERE, SUB)
    if not os.path.exists(sub_src):
        print("missing %s" % sub_src)
        return 1

    hours = {6.0: 392, 15.0: 162, 30.0: 96}
    made = []
    for rec in summ["decks"]:
        dg = rec["grain_um"]
        if a.only and abs(dg - a.only) > 1e-9:
            continue
        tag = "%gum" % dg
        d = os.path.join(SRC, tag)
        job = "micro_%s" % tag
        inp = os.path.join(d, job + ".inp")
        if not os.path.exists(inp):
            print("missing %s -- run _make_sag_paper.py --all first" % inp)
            return 1

        shutil.copy2(sub_src, os.path.join(d, SUB))
        fmt = dict(pad=("%g" % dg), job=job, sub=SUB, hours=hours[dg])
        with open(os.path.join(d, "run.bat"), "w", newline="\r\n") as fh:
            fh.write(RUN_BAT.format(**fmt))
        with open(os.path.join(d, "run.sh"), "w", newline="\n") as fh:
            fh.write(RUN_SH.format(**fmt))
        with open(os.path.join(d, "postprocess_odb.py"), "w",
                  newline="\n") as fh:
            fh.write(POST % dict(job=job))
        with open(os.path.join(d, "README.md"), "w", encoding="utf-8",
                  newline="\n") as fh:
            fh.write(README.format(**fmt))

        info = PADS[dg]
        with open(os.path.join(d, "EXPECTED.md"), "w", encoding="utf-8",
                  newline="\n") as fh:
            fh.write(EXPECTED.format(
                pad=("%g" % dg), observation=info["observation"],
                outcome=info["outcome"], falsify=info["falsify"],
                fn=rec["load_per_grain_n"], groove=rec["groove_width_nm"],
                chip_lo=rec["chip_measured_nm"][0],
                chip_hi=rec["chip_measured_nm"][1],
                indent=rec["indentation_nm"], hdc=rec["h_over_dc"],
                passes=rec["micro"]["slide_passes"]))

        mb = os.path.getsize(inp) / 1e6
        made.append((tag, mb, hours[dg]))
        print("  %-6s %-22s %7.1f MB   ~%d h" % (tag, job + ".inp", mb,
                                                 hours[dg]))

    print()
    print("each folder now holds: the .inp, %s, run.bat, run.sh," % SUB)
    print("postprocess_odb.py, README.md and EXPECTED.md -- copy it whole.")
    print()
    print("RUN THE 30 um FOLDER FIRST: cheapest (~96 h) and it is the one")
    print("that SHOULD show brittle fracture, so it is the fastest route to")
    print("knowing whether the mechanism reproduces the paper at all.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
