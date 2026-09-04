"""Build SAG_V2 -- the corrected decks, under names nothing else in this
project has ever used.

    python _make_sag_v2.py            # all three pads
    python _make_sag_v2.py --only 30

WHY A NEW NAME
--------------
Five rounds of fixes have shipped decks called ``micro_<pad>``, and a stale
copy of one is indistinguishable from a fresh one at the console: both print
``Abaqus JOB micro_30um``. That has now cost two full runs -- a job was started
against an old file twice, and the only way to tell was to compute the expected
total mass by hand and compare it against the packager's report.

So every name here is versioned. The job is ``sagv2_<pad>um``. If the console
says that, it is this deck. If it says ``micro_<pad>um``, it is an old one.
There is nothing to check by hand.

WHAT IS DIFFERENT FROM THE DECKS THAT FAILED
--------------------------------------------
* The cut is the MEASURED chip thickness -- 295 nm on the 30 um pad -- not the
  static Brinell indentation, which was 0.18 nm: three orders below what the
  paper measured for the same pad, one ninetieth of a single element, and
  smaller than a WC unit cell.
* The grain is seated 5% of that depth clear of the surface, so it closes the
  gap in the first 5% of the ramp. It used to be seated 2% of the BLOCK DEPTH
  clear, which was 200x the entire travel.
* The rigid grain carries its real diamond mass.
* The grain is displacement-driven, not force-driven: a massless rigid body
  cannot take a *Cload, and a 4e-16 tonne grain cannot be positioned by one.
* The datacheck runs under its own job name, so it cannot lock out the solve.
* Every ``abaqus`` line in run.bat is ``call``ed, and the verify flag is
  ``user_exp``.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "SAG_V2")
SUB = "vumat_grind2.for"
TAG = "sagv2"

RUN_BAT = """@echo off
cd /d "%~dp0"
rem  SAG_V2 -- {pad} um pad, WC-Co, after Ghosh et al. 2021.
rem
rem  JOB NAME: {job}
rem  If the console prints "Abaqus JOB {job}" you are running THIS deck.
rem  Anything else -- micro_{pad}um in particular -- is an older copy.
rem
rem  Stages, in order:
rem    0. abaqus verify -user_exp    can this machine build a subroutine
rem    1. datacheck, seconds, under its OWN job name so it cannot leave a
rem       lock file that blocks the solve
rem    2. the solve, ~{hours} h on 8 cores
rem    3. the postprocessor
rem
rem  Every abaqus line is `call`ed: on Windows abaqus is abaqus.bat, and one
rem  batch file running another without `call` never returns.
rem
rem  double=both is REQUIRED. The cut is compared against dc at tens of
rem  nanometres on a millimetre geometry; single precision has ~7 digits and
rem  the failure is SILENT.
rem
rem  The subroutine is vumat_grind2.for -- 58 constants, energy criterion.
rem  vumat_grind.for reads 56 and would misread the card.
call abaqus verify -user_exp
if errorlevel 1 (
  echo.
  echo  Abaqus cannot build a user subroutine on this machine.
  exit /b 1
)
call abaqus job={job}_check input={job}.inp user={sub} double=both cpus=1 datacheck
if errorlevel 1 (
  echo.
  echo  DATACHECK FAILED -- read {job}_check.dat. Nothing has been solved.
  exit /b 1
)
echo.
echo  datacheck passed. Solving {job} -- about {hours} h on 8 cores.
echo.
call abaqus job={job} input={job}.inp user={sub} double=both cpus=8 interactive
if errorlevel 1 exit /b 1
call abaqus python postprocess_odb.py
"""

RUN_SH = """#!/bin/sh
# SAG_V2 -- {pad} um pad. Job name: {job}. See run.bat for why each stage.
cd "$(dirname "$0")" || exit 1
abaqus verify -user_exp || {{ echo "no user-subroutine toolchain" >&2; exit 1; }}
abaqus job={job}_check input={job}.inp user={sub} double=both cpus=1 datacheck \\
  || {{ echo "DATACHECK FAILED -- read {job}_check.dat" >&2; exit 1; }}
echo "datacheck passed. Solving {job} -- about {hours} h on 8 cores."
abaqus job={job} input={job}.inp user={sub} double=both cpus=8 interactive || exit 1
abaqus python postprocess_odb.py
"""

READY = """# SAG_V2 · {pad} µm pad

**Job name: `{job}`.** If the Abaqus console prints `Abaqus JOB {job}`, this is
the deck running. Anything else — `micro_{pad}um` especially — is an older copy
from a previous round.

Self-contained. Copy the folder, run `run.bat` (Windows) or `run.sh`.

## What this deck cuts

| | |
|---|---|
| cut depth | **{chip:.0f} nm** — the paper's *measured* chip thickness for this pad |
| dc | 80 nm (measured, sections 4.2) |
| **chip / dc** | **{ratio:.2f}** |
| elements through the cut | {nel:.1f} at {el:.1f} nm |
| passes over one track | {passes} |
| estimate | ~{hours} h on 8 cores |

The paper measures chips directly and compares them against dc — that
comparison *is* its result. 6 µm gives 60–100 nm (at dc, pure ductile); 15 µm
gives 160–230; 30 µm gives 240–350 (both above dc, fracture observed).

The static Brinell indentation from eqs. 11–12 is a *different quantity* —
0.18 nm on the 30 µm pad, three orders below the measurement, 1/90th of one
element and smaller than a WC unit cell. A deck built on it runs for days with
a flat energy history. That is what the earlier rounds did.

## First 30 seconds tell you it is working

The packager prints the model mass. It should read **{mass:.5e}** — workpiece
plus the grain. If it reads {wp_mass:.5e}, the grain has no mass and the deck
is an old one.

Then watch KINETIC ENERGY. It should leave zero within the first few output
frames, once the grain closes its {standoff:.1f} nm standoff — about {pct:.0f}%
of the way through step LOAD. Energy that stays at exactly zero past the first
frame means nothing is touching.

## What to plot

**SDV13 is the result**: 1 = ductile, 2 = brittle. Plot it after every pass,
not only at the end — the criterion accumulates, so *when* a point flips is the
physics.

`postprocess_odb.py` runs automatically and writes `{job}_sdv13.csv` (brittle
fraction per frame per pass) and `{job}_summary.json`.

## Expected outcome

{expect}

## Caveats

- Johnson-Cook constants for WC-Co are **placeholders** except `A`. SDV13 is
  quotable; force magnitudes are not.
- The energy criterion is regularised by element length, so PSI is calibrated
  **for this mesh** ({el:.1f} nm through the depth).
"""

EXPECT = {
    6.0: "**Entirely ductile.** SDV13 = 1 everywhere, every pass. The paper "
         "reports pure ductile removal for this pad, and chip/dc = 1.00 sits "
         "exactly at the threshold. Any brittle element is a real "
         "disagreement.",
    15.0: "**Mostly ductile, some brittle.** chip/dc = 2.44. The paper sees "
          "*a few traces of brittle fracture ... significantly lower than the "
          "30 µm pad*. This deck should sit strictly between the other two.",
    30.0: "**Brittle regions present.** chip/dc = 3.69, the deepest cut of the "
          "three, and the paper sees *several traces of brittle fracture* "
          "alongside plastic flow. If this one stays entirely ductile through "
          "all its passes, the energy criterion accumulates too slowly — a "
          "finding about the criterion, not a bug in the deck.",
}


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", type=float, default=0.0)
    a = ap.parse_args(argv)

    from semgrit import sagdeck, sagemit
    from semgrit.quick import measure_images
    import _make_sag_paper as paper

    imgs = sorted(glob.glob(os.path.join(HERE, "B4C_1*.tif")))
    if not imgs:
        print("no B4C_1*.tif for the grain library")
        return 1
    print("measuring %d SEM image(s) ..." % len(imgs))
    got = measure_images(imgs, os.path.join(OUT, "_meas"), log=lambda *x: None)
    solids = got["solids"]
    print("  %d grain solids" % len(solids))
    print()

    hours = {6.0: 392, 15.0: 162, 30.0: 96}
    sub_src = os.path.join(HERE, SUB)
    recs = []
    for dg in (6.0, 15.0, 30.0):
        if a.only and abs(dg - a.only) > 1e-9:
            continue
        p = paper.params(dg, compression_mm=paper.PAPER["best_T_mm"],
                         speed_rpm=paper.PAPER["best_N_rpm"],
                         c10_mpa=paper.C10_CALIBRATED_MPA, passes=0)
        p.name = "%s_%gum" % (TAG, dg)
        pl = sagdeck.plan(p)
        job = "%s_%gum" % (TAG, dg)
        d = os.path.join(OUT, "%gum" % dg)
        os.makedirs(d, exist_ok=True)
        mi = sagemit.write_micro(os.path.join(d, job + ".inp"), pl, solids)

        shutil.copy2(sub_src, os.path.join(d, SUB))
        fmt = dict(pad=("%g" % dg), job=job, sub=SUB, hours=hours[dg])
        with open(os.path.join(d, "run.bat"), "w", newline="\r\n") as fh:
            fh.write(RUN_BAT.format(**fmt))
        with open(os.path.join(d, "run.sh"), "w", newline="\n") as fh:
            fh.write(RUN_SH.format(**fmt))
        with open(os.path.join(d, "postprocess_odb.py"), "w",
                  newline="\n") as fh:
            fh.write(__import__("_make_sag_packages").POST % dict(job=job))

        import math
        wp = mi["volume_mm3"] * 14500e-12
        with open(os.path.join(d, "README.md"), "w", encoding="utf-8",
                  newline="\n") as fh:
            fh.write(READY.format(
                pad=("%g" % dg), job=job, chip=mi["chip_depth_nm"],
                ratio=mi["chip_over_dc"], nel=mi["elements_through_chip"],
                el=mi["element_depth_mm"] * 1e6, passes=mi["n_passes"],
                hours=hours[dg], mass=wp + mi["grain_mass_tonne"],
                wp_mass=wp, standoff=mi["standoff_mm"] * 1e6,
                pct=100.0 * mi["standoff_over_indent"], expect=EXPECT[dg]))

        recs.append(dict(pad=dg, job=job,
                         path=os.path.relpath(mi["path"], HERE),
                         chip_nm=mi["chip_depth_nm"],
                         chip_over_dc=mi["chip_over_dc"],
                         elements=mi["elements"], passes=mi["n_passes"],
                         mb=round(mi["bytes"] / 1e6, 2),
                         model_mass_tonne=wp + mi["grain_mass_tonne"],
                         hours_8core=hours[dg]))
        print("%-6s %-16s %6.0f nm cut  chip/dc %.2f  %8s el  %6.2f MB"
              % ("%gum" % dg, job + ".inp", mi["chip_depth_nm"],
                 mi["chip_over_dc"], format(mi["elements"], ","),
                 mi["bytes"] / 1e6))

    with open(os.path.join(OUT, "SUMMARY.json"), "w", encoding="utf-8") as fh:
        json.dump(dict(tag=TAG, decks=recs), fh, indent=1)
    print()
    print("wrote %s" % OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
