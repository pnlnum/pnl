#!/bin/bash

# This script creates a release of PNL.
# Make sure, the script is located at the top level of a git repository of PNL

CWD="$(pwd)/$0"
PNLGIT=$(dirname $"CWD")
MINGW_PREFIX="/usr/local/opt/mingw-w64/toolchain-x86_64"

# Make sure lib{blas,lapack}.{lib,dll} can be found in LIBLAPACK_DIR
LIBLAPACK_DIR="/usr/local/src/libs/lapack-3.7.1/build-win64/lib/"

DATE=$(date +%Y%m%d)
# Local temporary directory
LOCAL_TMPDIR="$HOME/tmp/pnl-$DATE"
PNL_DIR=""
PNL_WINDIR=""
SOURCE_ONLY=0
VERSION=0

help() {
    echo "release-pnl.sh [--source] version
    --source    Only create an archive of the source files and manual.pdf.
    version     Number of the release to create.
    "
}

process_options() {
    # Handle options
    while test -n "$1"; do
        case "$1" in
        --source)
            SOURCE_ONLY=1
            echo "Creating a source only archive"
            shift
            ;;
        --help | -h)
            help
            exit 0
            ;;
        --)
            shift
            break
            ;;
        *)
            break
            ;;
        esac
    done
    # Handle arguments
    if [ -z "$1" ]; then
        echo "Version number is missing"
        help
        exit 1
    else
        VERSION=$1
        echo "Creating Version $VERSION"
    fi
}

# Set directories
# Arg 1 : version number
set_dirs() {
    PNL_DIR="$LOCAL_TMPDIR/pnl-$VERSION"
    PNL_WINDIR="$LOCAL_TMPDIR/pnl-win64-$VERSION"
    PNL_WINDIR="$LOCAL_TMPDIR/pnl-win64-$VERSION"
}

# Perform git archive
git_archive() {
    cwd=$(pwd)
    [[ -d "$LOCAL_TMPDIR" ]] && rm -rf "$LOCAL_TMPDIR"
    mkdir -p "$PNL_DIR"
    cd "$PNLGIT"
    git archive --format=tar master | tar -x -C "$PNL_DIR"
    mv $PNL_DIR/docs/manual-html $PNL_DIR
    mv $PNL_DIR/docs/pnl-manual.pdf $PNL_DIR
    rm -rf "$PNL_DIR/docs"
    cd "$cwd"
}

create_win_version() {
    cwd=$(pwd)
    mkdir -p "$PNL_WINDIR"
    mkdir -p "$LOCAL_TMPDIR/build-win"
    cd "$LOCAL_TMPDIR/build-win"
    PATH="$MINGW_PREFIX/bin:$MINGW_PREFIX/x86_64-w64-mingw32/bin:$MINGW_PREFIX/x86_64-w64-mingw32/lib:$LIBLAPACK_DIR:$PATH"
    cmake -DCROSS_COMPILE=ON -DCMAKE_BUILD_TYPE=Release -DPNL_INSTALL_PREFIX="$PNL_WINDIR" "$PNL_DIR"
    make || { echo "make failed" && exit 1; }
    make install || { echo "make install failed" && exit 1; }

    cd "$PNL_WINDIR"
    mv $PNL_WINDIR/lib/libpnl.dll.a $PNL_WINDIR/lib/libpnl.lib
    # Copy the manual
    cp -r "$PNL_DIR/manual-html" "$PNL_WINDIR"
    cp -r "$PNL_DIR/pnl-manual.pdf" "$PNL_WINDIR"
    # Copy all required dll's
    LIBS="libgcc_s_seh-1.dll libgfortran-5.dll libquadmath-0.dll"
    for lib in $LIBS; do
        cp $MINGW_PREFIX/x86_64-w64-mingw32/lib/$lib $PNL_WINDIR/lib
        cp -r $MINGW_PREFIX/x86_64-w64-mingw32/lib/$lib $LOCAL_TMPDIR/build-win/examples
    done
    BINS="libwinpthread-1.dll"
    for bin in $BINS; do
        cp $MINGW_PREFIX/x86_64-w64-mingw32/bin/$bin $PNL_WINDIR/lib
        cp -r $MINGW_PREFIX/x86_64-w64-mingw32/bin/$bin $LOCAL_TMPDIR/build-win/examples
    done
    LIBS="libblas.dll liblapack.dll"
    for lib in $LIBS; do
        cp $LIBLAPACK_DIR/$lib $PNL_WINDIR/lib
        cp -r $LIBLAPACK_DIR/$lib $LOCAL_TMPDIR/build-win/examples
    done
    cp lib/libpnl.dll $LOCAL_TMPDIR/build-win/examples

    cd "$LOCAL_TMPDIR"
    [[ -f "pnl-win64-$VERSION.zip" ]] && rm "pnl-win64-$VERSION.zip"
    zip -q -r "pnl-win64-$VERSION.zip" "pnl-win64-$VERSION"
    # Restore working directory
    cd "$cwd"
}

create_archive() {
    cd "$LOCAL_TMPDIR"
    tar -czf "pnl-$VERSION.tar.gz" "pnl-$VERSION"
}

process_options $@
set_dirs
git_archive
# compile_manual "$PNL_DIR"
create_archive

if [ $SOURCE_ONLY == 1 ]; then
    exit 0
else
    create_win_version
fi
