cmake_minimum_required(VERSION 2.6)
project(pnl C)
set(PROJECT_LIBRARY_NAME ${PROJECT_NAME})
set(LT_VERSION_CURRENT 1)
set(LT_VERSION_AGE 5)
set(LT_VERSION_REVISION 2)
set(BUILD_SHARED_LIBS ON)
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR})
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(pnl_INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib" )
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(pnl_USE_RPATH ON)
set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)


option(USE_MPI "Compiles the MPI bindings for PNL. Default = on" ON)
option(USE_INTERNAL_BLAS "Use internal Blas & Lapack instead of the system's ones. Default = off" OFF)

set(CMAKE_MODULE_PATH "${pnl_SOURCE_DIR}/CMakeModules;${CMAKE_MODULE_PATH}")
include(CMakeMacroLibtoolFile)

# Release or Debug
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif ()

enable_language(CXX OPTIONAL)
# Compiler Flags
if (UNIX)
    set(CMAKE_C_FLAGS_DEBUG "-Wall -g")
    set(CMAKE_C_FLAGS_RELEASE "-Wall -O2")
    set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g -fno-inline")
    set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O2")
endif (UNIX)

if (CMAKE_BUILD_TYPE MATCHES "Debug")
    set(MY_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
    set(MY_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
elseif (CMAKE_BUILD_TYPE MATCHES "Release")
    add_definitions(-DPNL_RANGE_CHECK_OFF)
    set(MY_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
    set(MY_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif ()

# Do we have inline
if (UNIX AND NOT APPLE)
    set(HAVE_INLINE ON)
else ()
    set(HAVE_INLINE OFF)
endif (UNIX AND NOT APPLE)

# Detect Blas & Lapack
if (NOT USE_INTERNAL_BLAS)
    find_package(BLAS)
    find_package(LAPACK)
    if (BLAS_FOUND AND LAPACK_FOUND)
        set(LIBS ${LIBS} ${BLAS_LIBRARIES})
        set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
    else ()
        set(USE_INTERNAL_BLAS ON)
    endif (BLAS_FOUND AND LAPACK_FOUND)
endif (NOT USE_INTERNAL_BLAS)

if (USE_MPI)
    find_package(MPI)
    set(LIBS ${LIBS} ${MPI_C_LIBRARIES})
    include_directories(${MPI_C_INCLUDE_PATH})
    set(MY_CPP_FLAGS "${MY_CPP_FLAGS} -I${MPI_C_INCLUDE_PATH}")
endif (USE_MPI)

# Add -lm if needed
if (UNIX)
    set(LIBS ${LIBS} m)
endif (UNIX)

LIST( REMOVE_DUPLICATES LIBS )

# Prepare for creating libtool .la
if (UNIX)
    foreach( _lib ${LIBS} )
        # Remove path to the library and suffixes. Transformation example: libXXX-2.0.(so|dylib) => XXX-2.0
        STRING(REGEX REPLACE ".*/lib(.*)\\.(dylib|so)" "\\1" _lib_no_path "${_lib}" )
        STRING(REGEX REPLACE "(.*)/lib[^/]*\\.(dylib|so)" "\\1" _lib_path "${_lib}" )
        if (NOT ${_lib_no_path} STREQUAL ${_lib})
            set(PNL_LDFLAGS "${PNL_LDFLAGS} ${CMAKE_LIBRARY_PATH_FLAG}${_lib_path}" )
        endif ()
        if (NOT ${_lib_no_path} STREQUAL ${_lib})
            set(PNL_LDFLAGS "${PNL_LDFLAGS} ${CMAKE_LINK_LIBRARY_FLAG}${_lib_no_path}" )
        endif ()
    endforeach( _lib )
    if (APPLE)
        set(PNL_LDFLAGS "${PNL_LDFLAGS} -lblas -llapack" )
    endif (APPLE)
    set(PNL_LDFLAGS "${PNL_LDFLAGS} -lm" )
endif (UNIX)

# Stuff for CMakeuser.incl
STRING(REGEX REPLACE ";" " -I" MY_CPP_FLAGS "${MY_CPP_FLAGS}" )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeuser.incl.in ${CMAKE_CURRENT_BINARY_DIR}/CMakeuser.incl @ONLY)

# Check math function from the standard library
set(CMAKE_REQUIRED_LIBRARIES ${LIBS})
CHECK_FUNCTION_EXISTS(tgamma HAVE_TGAMMA)
CHECK_FUNCTION_EXISTS(lgamma HAVE_LGAMMA)
CHECK_FUNCTION_EXISTS(exp10 HAVE_EXP10)
CHECK_FUNCTION_EXISTS(trunc HAVE_TRUNC)
CHECK_FUNCTION_EXISTS(isfinite HAVE_ISFINITE)
CHECK_FUNCTION_EXISTS(finite HAVE_FINITE)
CHECK_FUNCTION_EXISTS(isnan HAVE_ISNAN)
CHECK_FUNCTION_EXISTS(isinf HAVE_ISINF)

# Fortran stuff
enable_language(Fortran OPTIONAL)
if (CMAKE_Fortran_COMPILER_WORKS)
    message ("Fortran compiler found ${CMAKE_Fortran_COMPILER}")
    set(HAVE_FORTRAN_COMPILER ON)
    include(FortranCInterface)
    FortranCInterface_HEADER(src/include/pnl/FC.h MACRO_NAMESPACE "FC_")
    # test for dpstrf (Block Cholesky for non definite matrix)
    CHECK_FORTRAN_FUNCTION_EXISTS(dpstrf HAVE_DPSTRF)
else ()
    set(HAVE_FORTRAN_COMPILER OFF)
    if (APPLE) 
        set (HAVE_DPSTRF ON)
    endif (APPLE)
endif (CMAKE_Fortran_COMPILER_WORKS)


enable_testing()
add_subdirectory(src)
add_subdirectory(examples)

message ("Generate configure file ${CMAKE_CURRENT_BINARY_DIR}/config.h")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/src/include/pnl/config.h)



## Output configuration stuff ##
message("C Compiler :             ${CMAKE_C_COMPILER}")
message("     FLAGS :             ${MY_C_FLAGS}")
if (CMAKE_CXX_COMPILER_WORKS)
    message("C++ Compiler :           ${CMAKE_CXX_COMPILER}")
    message("     FLAGS   :           ${MY_CXX_FLAGS}")
endif (CMAKE_CXX_COMPILER_WORKS)

message("Installation directory : ${CMAKE_INSTALL_PREFIX}") 
if (BLAS_FOUND AND LAPACK_FOUND)
    message("Blas:        ${BLAS_LIBRARIES}")
    message("Lapack:      ${LAPACK_LIBRARIES}")
else ()
    message("Using internal Blas & Lapack")
endif (BLAS_FOUND AND LAPACK_FOUND)

if (MPI_C_FOUND)
    message("MPI:        ${MPI_C_LIBRARIES}")
else ()
    message ("No MPI installation found.")
endif (MPI_C_FOUND)
