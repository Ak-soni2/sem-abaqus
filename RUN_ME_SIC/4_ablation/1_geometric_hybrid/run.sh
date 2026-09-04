#!/bin/sh
set -e
cd "$(dirname "$0")"
abaqus verify -user_exp
abaqus job=geometric_hybrid input=geometric_hybrid.inp user=vumat_grind.for double=both cpus=1 datacheck
abaqus job=geometric_hybrid input=geometric_hybrid.inp user=vumat_grind.for double=both cpus=8 interactive
