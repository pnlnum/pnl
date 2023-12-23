#!/bin/bash

CWD="$(pwd)/$0"
CWD=$(dirname "$CWD")
PNLDIR="$CWD"

SRC_DIRS=". examples src/interpol src/librand src/linalg src/math src/mpi src/objects src/optim src/roots src/sort"

# We ignore libamos, libcephes, sobol files, lp_solve
ALL_SRC=$(find $SRC_DIRS -depth 1 \( -name '*.c' -o -name '*.h' \) -and ! -path '*lp_solve*' -and ! -path '*sobol*')
n_lines=$(cat $ALL_SRC | wc -l)
size_src=$(du -ch $ALL_SRC | grep total | awk '{print $1}')
echo "pour PNL : $n_lines lignes de code, soit $size_src"
