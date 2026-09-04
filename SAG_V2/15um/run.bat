@echo off
cd /d "%~dp0"
rem  SAG_V2 -- 15 um pad, WC-Co, after Ghosh et al. 2021.
rem
rem  JOB NAME: sagv2_15um
rem  If the console prints "Abaqus JOB sagv2_15um" you are running THIS deck.
rem  Anything else -- micro_15um in particular -- is an older copy.
rem
rem  Stages, in order:
rem    0. abaqus verify -user_exp    can this machine build a subroutine
rem    1. datacheck, seconds, under its OWN job name so it cannot leave a
rem       lock file that blocks the solve
rem    2. the solve, ~162 h on 8 cores
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
call abaqus job=sagv2_15um_check input=sagv2_15um.inp user=vumat_grind2.for double=both cpus=1 datacheck
if errorlevel 1 (
  echo.
  echo  DATACHECK FAILED -- read sagv2_15um_check.dat. Nothing has been solved.
  exit /b 1
)
echo.
echo  datacheck passed. Solving sagv2_15um -- about 162 h on 8 cores.
echo.
call abaqus job=sagv2_15um input=sagv2_15um.inp user=vumat_grind2.for double=both cpus=8 interactive
if errorlevel 1 exit /b 1
call abaqus python postprocess_odb.py
