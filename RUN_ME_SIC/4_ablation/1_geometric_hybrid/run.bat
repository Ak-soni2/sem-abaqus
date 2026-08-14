@echo off
abaqus verify -user_explicit
abaqus job=geometric_hybrid input=geometric_hybrid.inp user=vumat_grind.for double=both cpus=8 datacheck
abaqus job=geometric_hybrid input=geometric_hybrid.inp user=vumat_grind.for double=both cpus=8 interactive
