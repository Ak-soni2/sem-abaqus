@echo off
cd /d "%~dp0"
rem  Submit this deck. THREE steps, in order.
rem
rem    0. abaqus verify -user_exp
rem       Does the Fortran toolchain work at all? Every gate in this
rem       project compiled with gfortran, so this is the ONLY check that
rem       Abaqus can build a user subroutine on this machine.
rem    1. datacheck -- seconds. Reads the card and every keyword. The only
rem       real submission this project has made died here, on *User
rem       Material written four values to a line instead of eight.
rem    2. the solve.
rem
rem  double=both is REQUIRED, not a preference: the chip thickness is
rem  compared against a few nanometres on a 25 mm radius, a ratio of 1e-7,
rem  and single precision does not have the digits.
rem
rem  LICENCE: cpus=8 needs int(5*8^0.422) = 12 Abaqus tokens, against 5
rem  at cpus=1. Every wall clock in the README is the 8-core figure.
abaqus verify -user_exp
abaqus job=multi_abrasive input=multi_abrasive_field.inp user=vumat_grind.for double=both cpus=1 datacheck
if errorlevel 1 exit /b 1
abaqus job=multi_abrasive input=multi_abrasive_field.inp user=vumat_grind.for double=both cpus=8 interactive
