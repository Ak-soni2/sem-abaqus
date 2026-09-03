#!/bin/sh
# 6 um pad -- SAG on HVOF WC-Co, after Ghosh et al. 2021.
# See run.bat for why each stage is here. double=both is required.
cd "$(dirname "$0")" || exit 1
abaqus verify -user_explicit || {
  echo "Abaqus cannot build a user subroutine on this machine." >&2
  exit 1
}
abaqus job=micro_6um input=micro_6um.inp user=vumat_grind2.for double=both cpus=1 datacheck || {
  echo "DATACHECK FAILED -- read micro_6um.dat. Nothing has been solved yet." >&2
  exit 1
}
echo "datacheck passed. Solving -- about 392 h on 8 cores."
abaqus job=micro_6um input=micro_6um.inp user=vumat_grind2.for double=both cpus=8 interactive || exit 1
abaqus python postprocess_odb.py
