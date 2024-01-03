#!/bin/bash

BUILD=build-Debug
TMPDIR=/tmp
SED=gsed

# MPI functions to be removed
FUNCTIONS_TO_BE_REMOVED="pnl_object_load pnl_object_load_into_list pnl_object_mpi_bcast pnl_object_mpi_irecv pnl_object_mpi_isend pnl_object_mpi_pack pnl_object_mpi_pack_size pnl_object_mpi_recv pnl_object_mpi_reduce pnl_object_mpi_send pnl_object_mpi_ssend pnl_object_mpi_unpack pnl_object_save pnl_rng_create_from_file pnl_rng_save_to_file pnl_rng_state_mpi_pack pnl_rng_state_mpi_pack_size pnl_rng_state_mpi_unpack" 

FUNCTIONS_TO_BE_REMOVED=$(echo $FUNCTIONS_TO_BE_REMOVED | sed 's! !\\|!g')
FUNCTIONS_TO_BE_REMOVED="\($FUNCTIONS_TO_BE_REMOVED\)"
echo $FUNCTIONS_TO_BE_REMOVED

# Get the list of all _pnl function symbols and discard any MPI related stuff
nm $BUILD/src/libpnl.dylib | grep ' T ' | awk '{print $3;}' | grep -v "$FUNCTIONS_TO_BE_REMOVED" | grep _pnl | grep -v '\.' | sed 's/^_/    /' > $TMPDIR/new.def

# Compute the difference with the current pnl.def
diff -u -w src/pnl.def $TMPDIR/new.def > $TMPDIR/def.patch

# Delete the first hook
"$SED" -i~ '0,/@@/{/@@/d}' $TMPDIR/def.patch
"$SED" -i~ '3,/@@/{/@@/b;d;}' $TMPDIR/def.patch

# Apply the patch
patch -p0 src/pnl.def < $TMPDIR/def.patch
