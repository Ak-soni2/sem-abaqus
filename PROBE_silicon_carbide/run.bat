@echo off
rem Step 0: does the Fortran toolchain work at all?
abaqus verify -user_explicit
rem Step 1: preprocessing only. Seconds. Reads the card.
abaqus job=probe input=probe.inp user=vumat_grind.for double=both cpus=1 datacheck
rem Step 2: solve it.
abaqus job=probe input=probe.inp user=vumat_grind.for double=both cpus=1 interactive
