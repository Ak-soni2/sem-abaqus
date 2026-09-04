#!/bin/sh
set -e
cd "$(dirname "$0")"
abaqus verify -user_exp
abaqus job=forced_ductile_check input=forced_ductile.inp user=vumat_grind.for double=both cpus=1 datacheck
abaqus job=forced_ductile input=forced_ductile.inp user=vumat_grind.for double=both cpus=8 interactive
