#!/bin/bash

# This script creates a release of PNL.
# Make sure, the script is located at the top level of a git repository of PNL

CWD="$(pwd)/$0"
PNLGIT=$(dirname $"CWD")

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
}

# Perform git archive
git_archive() {
    cwd=$(pwd)
    [[ -d "$LOCAL_TMPDIR" ]] && rm -rf "$LOCAL_TMPDIR"
    mkdir -p "$PNL_DIR"
    cd "$PNLGIT"
    git archive --format=tar master | tar -x -C "$PNL_DIR"
    rm -rf "$PNL_DIR/docs"
    cd "$cwd"
}

# Compile the manual
# Arg 1: top level directory
compile_manual() {
    cd $1
    make -C man
    make -C man clean
}

create_win_version() {
    cwd=$(pwd)
    mkdir -p "$PNL_WINDIR"
    cd "$PNL_WINDIR/.."
    mkdir -p "$PNL_WINDIR/../build-win"
    cd "$PNL_WINDIR/../build-win"
    source /usr/local/src/mingw64-build/withmingw
    cmake -DCROSS_COMPILE=ON -DCMAKE_BUILD_TYPE=Release -DPNL_INSTALL_PREFIX="$PNL_WINDIR" "$PNL_DIR"
    make
    make install
    cd "$PNL_WINDIR"
    mv lib/libpnl.dll.a lib/libpnl.lib
    # Copy the manual
    cp -r "$PNL_DIR/man" "$PNL_WINDIR"
    # Copy all required dll's
    LIBS="libblas.dll libgcc_s_seh-1.dll libgfortran-3.dll liblapack.dll libquadmath-0.dll"
    for lib in $LIBS; do
        cp /usr/local/mingw64/x86_64-w64-mingw32/bin/$lib lib
    done
    cd "$LOCAL_TMPDIR"
    zip -r "pnl-win64-$VERSION.zip" "pnl-win64-$VERSION"
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
