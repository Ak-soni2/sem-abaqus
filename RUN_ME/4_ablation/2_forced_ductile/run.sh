#!/bin/sh
abaqus verify -user_explicit
abaqus job=forced_ductile input=forced_ductile.inp user=vumat_grind.for double=both cpus=8 datacheck
abaqus job=forced_ductile input=forced_ductile.inp user=vumat_grind.for double=both cpus=8 interactive
