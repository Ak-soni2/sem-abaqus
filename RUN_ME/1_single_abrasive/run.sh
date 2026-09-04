#!/bin/sh
set -e
cd "$(dirname "$0")"
#  Submit this deck. THREE steps, in order.
#
#    0. abaqus verify -user_exp
#       Does the Fortran toolchain work at all? Every gate in this
#       project compiled with gfortran, so this is the ONLY check that
#       Abaqus can build a user subroutine on this machine.
#    1. datacheck -- seconds. Reads the card and every keyword. The only
#       real submission this project has made died here, on *User
#       Material written four values to a line instead of eight.
#    2. the solve.
#
#  double=both is REQUIRED, not a preference: the chip thickness is
#  compared against a few nanometres on a 25 mm radius, a ratio of 1e-7,
#  and single precision does not have the digits.
#
#  LICENCE: cpus=8 needs int(5*8^0.422) = 12 Abaqus tokens, against 5
#  at cpus=1. Every wall clock in the README is the 8-core figure.
abaqus verify -user_exp
abaqus job=single_abrasive_check input=single_abrasive.inp user=vumat_grind.for double=both cpus=1 datacheck
abaqus job=single_abrasive input=single_abrasive.inp user=vumat_grind.for double=both cpus=8 interactive
