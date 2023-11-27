# PNL

[![Build Status](https://travis-ci.org/pnlnum/pnl.svg?branch=master)](https://travis-ci.org/pnlnum/pnl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/pnlnum/pnl?branch=master&svg=true)](https://ci.appveyor.com/project/pnlnum/pnl/branch/master)

This is PNL, a library for scientific computing. PNL is free software:
you can redistribute it and/or modify it under the terms of the GNU Lesser
General Public License.

PNL is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

## Installation

### Binary releases for Windows

Binary releases are compatible with Visual C++.
If you downloaded a binary release, just unzip the package and go to
Section [Using the library with Visual Studio](#under-visual-studio).

### Compiling the library

To compile PNL, you need CMake version >= 3.0. [Get CMake](http://cmake.org/cmake/resources/software.html).

Before compiling the library, users should bear in mind that

- Linear Algebra routines mainly rely on Blas & Lapack. If these two
  libraries are not found on the machine, the versions shipped with PNL are
  used. For a better performance, one should consider using `Atlas` under
  Linux.
- When an MPI library is detected on the computer, some MPI bindings
  are compiled within the library to enable the direct manipulation of
  `PnlObjects`.

To actually compile the library, just use CMake.

#### Under Unix

```shell
mkdir build
cd  build
cmake /relative/path/to/pnl
make
make install
```

The `make install` command installs

- the header files to `<prefix>/include`
- the library to `<prefix>/lib`
- the CMake config file `PnlConfig.cmake` to `<prefix>/lib/cmake/pnl`
- the CMakeuser.incl file to `<prefix>/share/pnl` (see [below](#with-a-makefile))

**The default value for `prefix` is the current build directory, but it can be changed by calling cmake as**

```shell
cmake -DPNL_INSTALL_PREFIX=some/new/prefix /relative/path/to/pnl
```

Some useful variables to modify the behavior of cmake.

- `-DCMAKE_BUILD_TYPE=Release/Debug`. Default is Debug. Choose Debug for building a development release without optimization and with debugging symbols. Choose Release for building an optimized version.

- `-DPNL_INSTALL_PREFIX=<path>`. Directory where to install the library. Default is to use the building directory as the installation prefix.

- `-DBUILD_SHARED_LIBS=ON/OFF`. Default is ON. If ON, PNL is built as a shared library and if OFF as a static library.

- `-DLAPACK_LIBRARIES=<path>`. Full path to a Lapack library (not just its directory). Lapack is detected automatically but the user can specify a particular library.

- `-DBLAS_LIBRARIES=<path>`. Full path to a Blas library (not just its directory). Blas is detected automatically but the user can specify a particular library. **Note that you must specify both `BLAS_LIBRARIES` and `LAPACK_LIBRARIES` or none of them.**

- `-DWITH_MPI=ON/OFF`. Default is ON. If ON, build the MPI bindings.

- `-DPNL_ENABLE_TESTS=ON/OFF`. Default is ON. If OFF, no test is compiled nor registered for running with `ctest`. This is intended to be used when PNL is included as a sub project instead of being compiled as an external library. In such a case, the typical usage is

    ```
    set(PNL_ENABLE_TESTS OFF CACHE BOOL "Disable PNL test.")
    add_subdirectory(some/path/to/pnl)
    include_directories(some/path/to/pnl/src/include)
    
    # Define some targets
    
    # For every target, add
    target_link_libraries(my_target pnl)
    # On Windows, to ensure all the dll's are copied next to the executable 
    # (on other platforms, pnl_add_postbuild does nothing).
    pnl_add_postbuild(my_target)
    ```

#### Under Windows

Use CMake to create the Visual C++ solution for the library. Once done,
open the solution in Visual C++ and be sure to call `generate` for both the
`ALL_BUILD` and `INSTALL` projects.

See the Unix section for the description of useful variables to modify
CMake's behaviour.

The generation of the `INSTALL` project takes care of installing the
headers and the library files (`.lib` and `.dll`) in the `build-dir` you
have specified in CMake (ie. the directory containing the Visual C++
solution).

### Getting the documentation

If you have cloned the git repository, you need to compile the
documentation yourself by going to the directory `man` and

- for the pdf version, run `make` (you need a `LaTeX` compiler).
- for the html version, run `make html` (you need `tex4ht`).

The directory `examples`, which is actually used for non regression tests,
contains usage examples for all the functionalities provided by the
library.

## Using the library

### With CMake

In your regular `CMakeLists.txt`, add

```shell
find_package(Pnl REQUIRED)
set(LIBS ${LIBS} ${PNL_LIBRARIES})
include_directories(${PNL_INCLUDE_DIRS})
# Deactivate PNL debugging stuff on Release builds
if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    add_definitions(-DPNL_RANGE_CHECK_OFF)
endif()

# For Windows only. To make sure all the required dll's are copied
# next to every executable, add the following instruction where
# <my_exec> is a target defined by add_executable.
pnl_add_postbuild(my_exec)
```

Note the instruction `pnl_add_postbuild`, which takes care of post build instructions when creating Visual Studio solutions from CMake. It basically copies all the required `.dll` to the directory containing the executable.

To build your project, call CMake with the following extra flag

```shell
-DCMAKE_PREFIX_PATH=path/to/pnl/install or path/to/pnl/build
```

A complete though basic `CMakeLists.txt` is available [there](perso/CMakeLists-example.txt).

### With a Makefile

You can use PNL in your own codes by creating `Makefile` containing

```shell
## Extra linker flags. Can be empty
LDFLAGS=

## Extra compiler flags. Can be empty.
CFLAGS=

## list of executables to create
BINS=pipo

## For each executable, create the variables
pipo_SRC= list_of_source_files

## This line must be the last one
include <prefix>/share/pnl/CMakeuser.incl
```

The `CMakeuser.incl` file can also be found in the root of the `build` directory.

See Section 1.3.1 of the [manual](https://pnlnum.github.io/pnl/manual-html/pnl-manual.html)  for more details on the syntax of this Makefile.

### With Visual Studio

If you want to use the previously compiled library in a new Visual C++
project without using `CMake`, you have to go through the followings steps

1. Set the configuration of the solution to __64 bits__.

        Solution properties -> Configuration

1. Add `pnl/include` as an additional include directory

        Project properties -> C/C++ -> General -> Additional Include Directories

1. Add `pnl/lib` as an additional library directory

        Project properties -> Linker -> General -> Additional Library Directories

1. Add `pnl.lib` as a dependency library

        Project properties -> Linker -> Input -> Additional Dependencies

1. To run your executable, copy all the `.dll` files from  `pnl/lib/` to the folder containing the executable.
