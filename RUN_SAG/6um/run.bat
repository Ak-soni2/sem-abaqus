@echo off
cd /d "%~dp0"
rem  6 um pad -- SAG on HVOF WC-Co, after Ghosh et al. 2021.
rem
rem  Three stages, in order, and the first two are cheap:
rem
rem    0. abaqus verify -user_explicit
rem       Can Abaqus build a user subroutine on THIS machine at all? If the
rem       Fortran toolchain is not wired up, everything after fails with a
rem       message that does not say so.
rem    1. datacheck -- seconds. Reads every keyword and the material card.
rem       The one real submission this project ever made died here, on
rem       *User Material written four values to a line instead of eight.
rem    2. the solve, ~392 h on 8 cores.
rem
rem  double=both is REQUIRED. h and dc are compared at 80 nm against a
rem  millimetre geometry; single precision has ~7 digits and does not have
rem  them. The failure is SILENT -- the branch flag comes out wrong and the
rem  job does not crash.
rem
rem  The subroutine is vumat_grind2.for, NOT vumat_grind.for: this deck
rem  carries 58 constants and the local energy criterion. vumat_grind.for
rem  reads 56 and would misinterpret the card.
abaqus verify -user_explicit
if errorlevel 1 (
  echo.
  echo  Abaqus cannot build a user subroutine on this machine.
  echo  Check that the Fortran compiler is on PATH and linked to Abaqus.
  exit /b 1
)
abaqus job=micro_6um input=micro_6um.inp user=vumat_grind2.for double=both cpus=1 datacheck
if errorlevel 1 (
  echo.
  echo  DATACHECK FAILED -- read micro_6um.dat, the error is a keyword or the
  echo  material card, not the physics. Nothing has been solved yet.
  exit /b 1
)
echo.
echo  datacheck passed. Starting the solve -- about 392 h on 8 cores.
echo.
abaqus job=micro_6um input=micro_6um.inp user=vumat_grind2.for double=both cpus=8 interactive
if errorlevel 1 exit /b 1
abaqus python postprocess_odb.py
