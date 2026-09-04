@echo off
cd /d "%~dp0"
rem Step 0: does the Fortran toolchain work at all?
call abaqus verify -user_exp
rem Step 1: preprocessing only. Seconds. Reads the card.
call abaqus job=probe_check input=probe.inp user=vumat_grind.for double=both cpus=1 datacheck
if errorlevel 1 exit /b 1
rem Step 2: solve it.
call abaqus job=probe input=probe.inp user=vumat_grind.for double=both cpus=1 interactive
