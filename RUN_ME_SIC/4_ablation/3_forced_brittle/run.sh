#!/bin/sh
set -e
cd "$(dirname "$0")"
abaqus verify -user_explicit
abaqus job=forced_brittle input=forced_brittle.inp user=vumat_grind.for double=both cpus=1 datacheck
abaqus job=forced_brittle input=forced_brittle.inp user=vumat_grind.for double=both cpus=8 interactive
