#!/bin/sh
# SAG_V2 -- 30 um pad. Job name: sagv2_30um. See run.bat for why each stage.
cd "$(dirname "$0")" || exit 1
abaqus verify -user_exp || { echo "no user-subroutine toolchain" >&2; exit 1; }
abaqus job=sagv2_30um_check input=sagv2_30um.inp user=vumat_grind2.for double=both cpus=1 datacheck \
  || { echo "DATACHECK FAILED -- read sagv2_30um_check.dat" >&2; exit 1; }
echo "datacheck passed. Solving sagv2_30um -- about 96 h on 8 cores."
abaqus job=sagv2_30um input=sagv2_30um.inp user=vumat_grind2.for double=both cpus=8 interactive || exit 1
abaqus python postprocess_odb.py
