# Changelog

All notable changes to PNL will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [1.9.2] - 2017-12-24

### Fixed

- Fix access to global variables under Windows. This was causing a segfault in `random_test`.

### Added

- Add `BUILD_SHARED_LIBS` as a CMake option to set the type of library to build. Default value in `ON`. When set to `OFF`, PNL is built as a static library.

## [1.9.1] - 2017-12-13

### Added

- Add option `PNL_ENABLE_TESTS=OFF` to CMake to disable tests. Default value is s`ON`.

### Changed

- Update the `CMakeuser.incl` mechanism.
- Update some installation path to be more consistent with Linux standards:
  - `CMakeuser.incl` is installed to `<prefix>/share/pnl`,
  - `PnlConfig.cmake` is installed to `<prefix>/lib/cmake/pnl/`.

### Fixed

- Fix minor bugs in `pnl_sp_mat_isequal` and `mtherr`.
- Fix some compilation warnings in LP_Solve.

## [1.9.0] - 2017-11-30

### Added

- Read sparse matrices from file.
- Update read functions to deal with comments.
- Add floating point comparison functions for vectors and matrices.

## [1.8.0] - 2017-10-04

### Added

- Add a linear programming routine based on the Simplex implementation of LP_Solve.

### Fixed

- Fix a few bugs.

## [1.7.5] - 2017-01-18

### Changed

- Improve CMake's integration for Visual Studio.
- The binary version pnl-win64-1.7.5 is compatible with Visual Studio up to VS2013 but not with VS2015. If you use VS2015, download the source files and compile them yourself.

### Fixed

- Fix minor bugs.

## [1.7.4] - 2016-03-21

### Added

- Provide a `PnlConfig.cmake` file to use CMake `find_package` command in config mode.
- Enable in-place `MPI` reductions of `PNLObject`s.

### Changed

- Improve the efficiency of creating and evaluating `PnlBasis` objects.

### Fixed

- Fix compiling errors in the examples under MSVC.

## [1.7.3] - 2016-02-15

### Added

- Add non tensor functions to `PnlBasis` objects.
- Add Complex error functions based on the `Faddeeva` package.

## [1.7.2] - 2015-09-10

### Added

- A CMake module to include the library in other projects.
- The non central chi squared distribution to the random number generation toolbox.

### Changed

- Improves and fixes a bug in the pnl_basis module.

## [1.7.1] -  2014-07-10

### Fixed

- Minor fixes in the examples
- Improvement of the inline facility detection in CMake. This fixes a bug of the 1.7.0 under Windows.

## [1.7.0] -  2014-04-28

### Added

- Random number generation. All methods know how to handle QMC generators. Three new distributions have been added: log-normal, inverse gaussian and asymmetric double exponential distributions.
- Linear Algebra. Computations of the logarithm and exponential of complex matrices. Sparse matrix object.
- Optimization. Newton's algorithm with Armijo line search.

### Changed

- Three new distributions are added: log-normal, inverse gaussian and asymmetric double exponential distributions.
- The evaluation of a multivariate polynomial is greatly improved thanks to the use of sparse storage of these objects. This new implementation relies on sparse matrices.
- The internal structure of the top level object is modified to keep track of the number references on to a given object, which has enabled us to improve by a great deal the memory management of lists and arrays.

## [1.6.0] - 2013-03-22

### Added

- Add a few functions for rounding (some are part of C99).

### Changed

- Improve Blas & Lapack detection on Mac OS X.

### Fixed

- Fix a minor bug in two dimensional FFT.

## [1.5.2] -  2013-01-23

### Changed

- The compilation chain has moved from the `Autotools` to CMake.

## [1.5.1] -  2012-10-16

### Added

- `MPI_Reduce` binding for `PnlVect` and `PnlMat` objects.
- 2D FFT.
- The user can define new function bases and register them for further usage as native bases.
- The PnlObject structure has three new function members: clone, copy, new. It is mainly useful for inheritance purposes.

## [1.5.0] -  2012-03-19

### Added

- Linear solver using a QR decomposition with column pivoting
- Inverse of permutations and its application in place
- Approximations for inverse hyperbolic functions (part of C99 but missing under Visual)

### Changed

- Most macros are renamed with a pnl_ prefix to avoid name clashes

## [1.4.1] -  2011-10-11

### Added

- New functions to extract data from Hmatrices as vectors or matrices
- New PnlArray type: array of PnlObjects. It is often a good alternative to using Hmatrices.

### Changed

- Redesign of how Hmatrices are stored in memory

### Fixed

- A bug in `pnl_mat_{lower,upper}_syslin` has been fixed

## [1.4.0] -  2011-09-09

### Added

- Runge Kutta Fehler 45 integrator for n dimensional ODEs
- LU decomposition of tridiagonal matrices
- Cholesky factorization with complete pivoting for positive semidefinite matrices
- New threadsafe Sobol generators (from John Burkardt)

### Changed

- Design of a new unit test framework
- New organisation of the manual
- Update Lapack to version 3.2.1
- Integration of the Mersenne Twister Dynamic Creator version 0.6.1

## [1.3.3] -  2011-03-10

### Added

- Generation of Bessel random variables
- A generator for parallel computing (MPI) is added : Dynamic created Mersenne Twister
- Sparse polynomial bases based on hyperbolic sets of indices. Bases can be centered and normalized.
- Save/Load interface based on MPI Pack/Unpack

### Fixed

- Random number generators become thread-safe.

## [1.2.0] - 2010-07-15

### Added

- Bindings for MPI

### Changed

- Update of random number generators to work on 64-bit machines