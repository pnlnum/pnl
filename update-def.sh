#!/bin/bash

BUILD=build-Debug
TMPDIR=/tmp
SED=gsed

# Get the list of all _pnl function symbols
nm $BUILD/lib/libpnl.dylib | grep ' T ' | awk '{print $3;}' | grep _pnl | grep -v '\.' | sed 's/^_/    /' > $TMPDIR/new.def

# Compute the difference with the current pnl.def
diff -u -w src/pnl.def $TMPDIR/new.def > $TMPDIR/def.patch

# Delete the first hook
"$SED" -i~ '0,/@@/{/@@/d}' $TMPDIR/def.patch
"$SED" -i~ '3,/@@/{/@@/b;d;}' $TMPDIR/def.patch

# Apply the patch
patch -p0 src/pnl.def < $TMPDIR/def.patch
