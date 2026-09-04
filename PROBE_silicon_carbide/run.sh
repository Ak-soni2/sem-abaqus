#!/bin/sh
set -e
cd "$(dirname "$0")"
abaqus verify -user_exp
abaqus job=probe input=probe.inp user=vumat_grind.for double=both cpus=1 datacheck
abaqus job=probe input=probe.inp user=vumat_grind.for double=both cpus=1 interactive
