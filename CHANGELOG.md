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

- Add `PNL_ENABLE_TESTS=OFF` as an option to CMake to disable tests. Default value is s`ON`.

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