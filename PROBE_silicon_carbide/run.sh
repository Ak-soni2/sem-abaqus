#!/bin/sh
abaqus verify -user_explicit
abaqus job=probe input=probe.inp user=vumat_grind.for double=both cpus=1 datacheck
abaqus job=probe input=probe.inp user=vumat_grind.for double=both cpus=1 interactive
