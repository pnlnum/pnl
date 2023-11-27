cmake_minimum_required(VERSION 3.15)
cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules;${CMAKE_MODULE_PATH}")
set(PROJECT_LIBRARY_NAME ${PROJECT_NAME})
file(STRINGS "VERSION" PNL_VERSION)


# Prefix
if(PNL_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX ${PNL_INSTALL_PREFIX} CACHE STRING "Installation prefix." FORCE)
else(PNL_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR} CACHE STRING "Installation prefix." FORCE)
endif(PNL_INSTALL_PREFIX)

if(CROSS_COMPILE)
    include(cross-compile)
    set(MINGW ON)
endif()

project(pnl C CXX)

set(CMAKE_MACOSX_RPATH ON)
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(pnl_INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib" )
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_FILES_DIRECTORY "/CMakeFiles")


option(WITH_MPI "Compiles the MPI bindings for PNL. Default = on" ON)
option(WITHOUT_DEPRECATED "Do not build deprecated functions. Default = off" OFF)
option(PNL_ENABLE_TESTS "Enable PNL tests. Default = on" ON)
option(BUILD_SHARED_LIBS "Build PNL as a shared library. Default = on. If off, PNL is built as a static library" ON)


# Release or Debug
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES Clang OR "${CMAKE_CXX_COMPILER_ID}" MATCHES GNU)
    set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wall")
endif()
if ("${CMAKE_C_COMPILER_ID}" MATCHES Clang OR "${CMAKE_C_COMPILER_ID}" MATCHES GNU)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
endif()

if (CMAKE_BUILD_TYPE MATCHES "Debug")
    set(PNL_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
    set(PNL_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
elseif (CMAKE_BUILD_TYPE MATCHES "Release")
    add_definitions(-DPNL_RANGE_CHECK_OFF)
    set(PNL_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
    set(PNL_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif ()


# Do we have inline
include(CMakeMacroExternInline)

set(PNL_DLL OFF)
if (WIN32)
    add_definitions(-DMSDOS)
    if (MSVC)
        add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
    endif (MSVC)
    if (BUILD_SHARED_LIBS)
        set(PNL_DLL ON)
    endif (BUILD_SHARED_LIBS)
endif (WIN32)


# Detect Blas & Lapack
# Has the user supplied a Blas Library?
if ( (BLAS_LIBRARIES AND (NOT LAPACK_LIBRARIES)) OR ((NOT BLAS_LIBRARIES) AND LAPACK_LIBRARIES) )
    message(FATAL_ERROR "--> You must specify both BLAS_LIBRARIES and LAPACK_LIBRARIES or none of them.")
endif ()

if (BLAS_LIBRARIES AND LAPACK_LIBRARIES)
    set(BLAS_FOUND true)
    set(LAPACK_FOUND true)
else ()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endif (BLAS_LIBRARIES AND LAPACK_LIBRARIES)

if (BLAS_FOUND AND LAPACK_FOUND)
    set(LIBS ${LIBS} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
else()
    message(FATAL_ERROR "No Blas & Lapack installation found.")
endif (BLAS_FOUND AND LAPACK_FOUND)

# Check if Fortran is available
# enable_language(Fortran OPTIONAL)
if (CMAKE_Fortran_COMPILER_WORKS)
    message("Fortran compiler found ${CMAKE_Fortran_COMPILER}")
    set(PNL_HAVE_FORTRAN_COMPILER ON)
    include(FortranCInterface)
    FortranCInterface_HEADER(src/include/pnl/FC.h MACRO_NAMESPACE "FC_")
else ()
    set(PNL_HAVE_FORTRAN_COMPILER OFF)
endif (CMAKE_Fortran_COMPILER_WORKS)
# Test if dpsrtrf is present
include(CheckFunctionExists)
CHECK_FUNCTION_EXISTS(dpstrf PNL_HAVE_DPSTRF)
if (NOT PNL_HAVE_DPSTRF)
    CHECK_FUNCTION_EXISTS(dpstrf_ PNL_HAVE_DPSTRF)
endif (NOT PNL_HAVE_DPSTRF)


if (WITH_MPI)
    set(MPI_FIND_QUIETLY true)
    find_package(MPI)
    if (MPI_C_FOUND)
        set(LIBS ${LIBS} ${MPI_C_LIBRARIES})
        include_directories(${MPI_C_INCLUDE_PATH})
        set(PNL_CPP_FLAGS "${PNL_CPP_FLAGS} -I${MPI_C_INCLUDE_PATH}")
    endif (MPI_C_FOUND)
endif (WITH_MPI)

# Add -lm and -ldl if needed
if (UNIX)
    set(LIBS ${LIBS} m dl)
endif (UNIX)

# Stuff for CMakeuser.incl

# Check math function from the standard library
set(CMAKE_REQUIRED_LIBRARIES ${LIBS})

# Check some math functions
include(CheckFunctionExists)
include(CheckSymbolExists)
check_symbol_exists(tgamma math.h PNL_HAVE_TGAMMA)
check_symbol_exists(lgamma math.h PNL_HAVE_LGAMMA)
check_symbol_exists(exp10 math.h PNL_HAVE_EXP10)
check_symbol_exists(trunc math.h PNL_HAVE_TRUNC)
check_symbol_exists(round math.h PNL_HAVE_ROUND)
check_symbol_exists(isfinite math.h PNL_HAVE_ISFINITE)
check_symbol_exists(finite math.h PNL_HAVE_FINITE)
check_symbol_exists(isnan math.h PNL_HAVE_ISNAN)
check_symbol_exists(isinf math.h PNL_HAVE_ISINF)

add_subdirectory(src)

if (PNL_ENABLE_TESTS)
    enable_testing()
    add_subdirectory(examples)
    if (UNIX)
        STRING(REGEX REPLACE ";" " -I" PNL_CPP_FLAGS "${PNL_CPP_FLAGS}" )
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeuser.incl.in ${CMAKE_CURRENT_BINARY_DIR}/CMakeuser.incl @ONLY)
        add_custom_target(CMakeuser ALL
            COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/split_linker_command.sh ${CMAKE_INSTALL_PREFIX} ${CMAKE_CURRENT_BINARY_DIR})
        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CMakeuser.incl DESTINATION share/pnl)
    endif (UNIX)
endif (PNL_ENABLE_TESTS)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/src/include/pnl/pnl_config.h)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/src/include/pnl/pnl_config.h DESTINATION include/pnl)

## Output configuration stuff ##
message("C Compiler: ${CMAKE_C_COMPILER}")
message("C FLAGS: ${PNL_C_FLAGS}")
if (CMAKE_CXX_COMPILER_WORKS)
    message("C++ Compiler: ${CMAKE_CXX_COMPILER}")
    message("C++ FLAGS: ${PNL_CXX_FLAGS}")
endif (CMAKE_CXX_COMPILER_WORKS)

message("Installation directory : ${CMAKE_INSTALL_PREFIX}")
if (BLAS_FOUND AND LAPACK_FOUND)
    message("Blas: ${BLAS_LIBRARIES}")
    message("Lapack: ${LAPACK_LIBRARIES}")
else ()
    message("Using internal Blas & Lapack")
endif (BLAS_FOUND AND LAPACK_FOUND)

if (MPI_C_FOUND)
    message("MPI: ${MPI_C_LIBRARIES}")
else ()
    message ("No MPI installation found.")
endif (MPI_C_FOUND)
