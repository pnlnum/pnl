This is PNL, a library for scientific computing. PNL is free software:
you can redistribute it and/or modify it under the terms of the GNU Lesser
General Public License. 

PNL is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.


# Installation

## Binary releases for Windows

Binary releases are compatible with Visual C++ up to VS2013. 
If you downloaded a binary release, just unzip the package and go to
Section [Using the library](#under-windows-with-visual).


## Compiling the library


To compile PNL, you need CMake version >= 2.8. [Get CMake](http://cmake.org/cmake/resources/software.html).

Before compiling the library, users should bear in mind that 
- Linear Algebra routines mainly rely on Blas & Lapack. If these two
  libraries are not found on the machine, the versions shipped with PNL are
  used. For a better performance, one should consider using `Atlas` under
  Linux.
- When an MPI library is detected on the computer, some MPI bindings
  are compiled within the library to enable the direct manipulation of
  `PnlObjects`. 


To actually compile the library, just use CMake. 

### Under Unix

```
mkdir pnl-build
cd  pnl-build
cmake /relative/path/to/pnl
make
make install
```

Note the command `make install`, which installs
- the header files to `pnl-build/include`
- the library and the CMake config file to `pnl-build/lib`


Some useful variables to modify the behaviour of cmake.

```
-DCMAKE_BUILD_TYPE        Release OR Debug. Default is Debug.
        Choose Debug for building a development release without
        optimization and with debugging symbols.
        Choose Release for building an optimized version.

-DPNL_INSTALL_PREFIX     The path where to install the library.
        Default is to use the building directory as the installation
        prefix.

-DLAPACK_LIBRARIES        Full path of a Lapack library (not just its
        directory). Lapack is detected automatically but the user
        can specify a particular library.

-DBLAS_LIBRARIES          Full path of a Blas library (not just its
        directory). Blas is detected automatically but the user
        can specify a particular library.

    Note that you must specify both BLAS_LIBRARIES and
    LAPACK_LIBRARIES or none of them.

-DUSE_MPI                 ON/OFF. Default is ON. If ON, build the MPI
        bindings.

-DUSE_INTERNAL_BLAS       ON/OFF. Default is OFF. If ON, use the
        internal Blas & Lapack libraries shipped with the PNL
        source code.

-DLINK_TYPE               SHARED or STATIC. Default is SHARED.
                          Determine which type of library to build.

```

### Under Windows

Use CMake to create the Visual C++ solution for the library. Once done,
open the solution in Visual C++ and be sure to call `generate` for both the
`ALL_BUILD` and `INSTALL` projects.


See the Unix section for the description of useful variables to modify
CMake's behaviour.

The generation of the `INSTALL` project takes care of installing the
headers and the library files (`.lib` and `.dll`) in the `build-dir` you
have specified in CMake (ie. the directory containing the Visual C++
solution).



## Getting the documentation

If you have cloned the git repository, you need to compile the
documentation yourself by going to the directory `man` and 
- for the pdf version, run `make` (you need a `LaTeX` compiler). 
- for the html version, run `make html` (you need `tex4ht`). 

The directory `examples`, which is actually used for non regression tests,
contains usage examples for all the functionalities provided by the
library.

# Using the library

## Under Unix

### Using CMake

In your regular `CMakeLists.txt`, add
```
find_package(Pnl REQUIRED)
set(LIBS ${LIBS} ${PNL_LIBRARIES})
include_directories(${PNL_INCLUDE_DIRS})
# Deactivate PNL debugging stuff on Release builds
if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    add_definitions(-DPNL_RANGE_CHECK_OFF)
endif()
```

To build your project, call CMake with the following extra flag
```
-DCMAKE_PREFIX_PATH=path/to/build-dir
```

A complete though basic CMakeLists.txt is avaible [there](perso/CMakeLists-example.txt).

### Using a Makefile


In whatever project you want to use PNL, create a `Makefile` containing
```
## Extra linker flags. Can be empty
LDFLAGS=

## Extra compiler flags. Can be empty.
CFLAGS=

## list of executables to create
BINS=pipo

## For each executable, create the variables
pipo_SRC= list_of_source_files

## This line must be the last one
include /path/to/build-dir/CMakeuser.incl
```
See the [manual](https://pnlnum.github.io/pnl/manual-html/pnl-manual.html) section 1.3.1 for more details on the syntax of this Makefile.



## Under Windows with Visual

If you want to use the previously compiled library in a new Visual C++
project, you have to go through the followings steps

1. Set the configuration of the solution to __64 bits__.
```
Solution properties -> Configuration
```
2. Add `build-dir/include` as an additional include directory
```
Project properties -> C/C++ -> General -> Additional Include Directories
```
3. Add `build-dir/lib` as an additional library directory
```
Project properties -> Linker -> General -> Additional Library Directories
```
4. Add `pnl.lib` as a dependency library
```
Project properties -> Linker -> Input -> Additional Dependencies
```
5. To run your executable, copy the file `build-dir/lib/pnl.dll` to the folder containing the executable.
