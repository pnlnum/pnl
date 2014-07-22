cmake_minimum_required(VERSION 2.8.9)
set(PROJECT_LIBRARY_NAME ${PROJECT_NAME})
set(VERSION_CURRENT 1)
set(VERSION_AGE 6)
set(VERSION_REVISION 0)
set(BUILD_SHARED_LIBS ON)

# Prefix 
if (NOT PREFIX)
    set (PREFIX ${CMAKE_CURRENT_BINARY_DIR})
endif (NOT PREFIX)
set(CMAKE_INSTALL_PREFIX ${PREFIX} CACHE STRING "Installation prefix." FORCE)

# Compilers
if (CC)
    set (CMAKE_C_COMPILER ${CC})
endif (CC)
if (CXX)
    set (CMAKE_CXX_COMPILER ${CXX})
endif (CXX)

project(pnl C)

if (NOT LINK_TYPE) 
    set(LINK_TYPE SHARED)
endif (NOT LINK_TYPE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(pnl_INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib" )
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_MACOSX_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(pnl_USE_RPATH ON)
set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
set(CMAKE_FILES_DIRECTORY "/CMakeFiles")


option(USE_MPI "Compiles the MPI bindings for PNL. Default = on" ON)
option(WITHOUT_DEPRECATED "Build deprecated functions. Default = off" OFF)
option(USE_INTERNAL_BLAS "Use internal Blas & Lapack instead of the system's ones. Default = off" OFF)

set(CMAKE_MODULE_PATH "${pnl_SOURCE_DIR}/CMakeModules;${CMAKE_MODULE_PATH}")

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
    set(PNL_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
    set(PNL_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
elseif (CMAKE_BUILD_TYPE MATCHES "Release")
    add_definitions(-DPNL_RANGE_CHECK_OFF)
    set(PNL_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
    set(PNL_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif ()

# Do we have inline
include(CMakeMacroExternInline)

if (WIN32)
    add_definitions (-DMSDOS)
    if (MSVC)
        set (USE_INTERNAL_BLAS ON)
    endif (MSVC)
endif (WIN32)


# Detect Blas & Lapack
if (NOT USE_INTERNAL_BLAS)
    # Has the user supplied a Blas Library?
    if ( (BLAS_LIBRARIES AND (NOT LAPACK_LIBRARIES)) OR ((NOT BLAS_LIBRARIES) AND LAPACK_LIBRARIES) )
        message (FATAL_ERROR "--> You must specify both BLAS_LIBRARIES and LAPACK_LIBRARIES or none of them.")
    endif ()

    if (BLAS_LIBRARIES AND LAPACK_LIBRARIES) 
        set (BLAS_FOUND true)
        set (LAPACK_FOUND true)
    else ()
        find_package(BLAS)
        find_package(LAPACK)
    endif (BLAS_LIBRARIES AND LAPACK_LIBRARIES) 
    
    if (BLAS_FOUND AND LAPACK_FOUND)
        set(LIBS ${LIBS} ${BLAS_LIBRARIES})
        set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
    else ()
        set(USE_INTERNAL_BLAS ON)
    endif (BLAS_FOUND AND LAPACK_FOUND)
endif (NOT USE_INTERNAL_BLAS)

if (USE_INTERNAL_BLAS)
    set(PNL_HAVE_FORTRAN_COMPILER OFF)
    set (PNL_HAVE_DPSTRF ON)
else (USE_INTERNAL_BLAS)
    # Check if Fortran is available
    enable_language(Fortran OPTIONAL)
    if (CMAKE_Fortran_COMPILER_WORKS)
        message ("Fortran compiler found ${CMAKE_Fortran_COMPILER}")
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
endif (USE_INTERNAL_BLAS)


if (USE_INTERNAL_BLAS)
endif (USE_INTERNAL_BLAS)

if (USE_MPI)
    set (MPI_FIND_QUIETLY true)
    find_package(MPI)
    if (MPI_C_FOUND)
        set(LIBS ${LIBS} ${MPI_C_LIBRARIES})
        include_directories(${MPI_C_INCLUDE_PATH})
        set(PNL_CPP_FLAGS "${PNL_CPP_FLAGS} -I${MPI_C_INCLUDE_PATH}")
    endif (MPI_C_FOUND)
endif (USE_MPI)

# Add -lm if needed
if (UNIX)
    set(LIBS ${LIBS} m)
endif (UNIX)

# Prepare for creating libtool .la
# include(CMakeMacroLibtoolFile)
# if (UNIX)
#     LIST( REMOVE_DUPLICATES LIBS )
#     foreach( _lib ${LIBS} )
#         # Remove path to the library and suffixes. Transformation example: libXXX-2.0.(so|dylib) => XXX-2.0
#         STRING(REGEX REPLACE ".*/lib(.*)\\.(dylib|so)" "\\1" _lib_no_path "${_lib}" )
#         STRING(REGEX REPLACE "(.*)/lib[^/]*\\.(dylib|so)" "\\1" _lib_path "${_lib}" )
#         if (NOT ${_lib_no_path} STREQUAL ${_lib})
#             set(PNL_LDFLAGS "${PNL_LDFLAGS} ${CMAKE_LIBRARY_PATH_FLAG}${_lib_path}" )
#         endif ()
#         if (NOT ${_lib_no_path} STREQUAL ${_lib})
#             set(PNL_LDFLAGS "${PNL_LDFLAGS} ${CMAKE_LINK_LIBRARY_FLAG}${_lib_no_path}" )
#         endif ()
#     endforeach( _lib )
#     if (APPLE)
#         set(PNL_LDFLAGS "${PNL_LDFLAGS} -lblas -llapack" )
#     endif (APPLE)
#     set(PNL_LDFLAGS "${PNL_LDFLAGS} -lm" )
# endif (UNIX)

# Stuff for CMakeuser.incl
STRING(REGEX REPLACE ";" " -I" PNL_CPP_FLAGS "${PNL_CPP_FLAGS}" )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeuser.incl.in ${CMAKE_CURRENT_BINARY_DIR}/CMakeuser.incl @ONLY)

# Check math function from the standard library
set(CMAKE_REQUIRED_LIBRARIES ${LIBS})

# Check some math functions
include(CheckFunctionExists)
include(CheckSymbolExists)
CHECK_FUNCTION_EXISTS(tgamma PNL_HAVE_TGAMMA)
CHECK_FUNCTION_EXISTS(lgamma PNL_HAVE_LGAMMA)
CHECK_FUNCTION_EXISTS(exp10 PNL_HAVE_EXP10)
CHECK_FUNCTION_EXISTS(trunc PNL_HAVE_TRUNC)
CHECK_FUNCTION_EXISTS(round PNL_HAVE_ROUND)
CHECK_SYMBOL_EXISTS(isfinite "math.h" PNL_HAVE_ISFINITE)
CHECK_FUNCTION_EXISTS(finite PNL_HAVE_FINITE)
CHECK_FUNCTION_EXISTS(isnan PNL_HAVE_ISNAN)
CHECK_FUNCTION_EXISTS(isinf PNL_HAVE_ISINF)

add_subdirectory(src)

if (NOT SUB_CMAKE_CALL)
    enable_testing()
    add_subdirectory(examples)
endif (NOT SUB_CMAKE_CALL)


message ("Generate configure file ${CMAKE_CURRENT_BINARY_DIR}/src/include/pnl/pnl_config.h")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake 
    ${CMAKE_CURRENT_BINARY_DIR}/src/include/pnl/pnl_config.h)
install(FILES 
    ${CMAKE_CURRENT_BINARY_DIR}/src/include/pnl/pnl_config.h
    DESTINATION include/pnl
)




# Add the uninstall target
if (NOT SUB_CMAKE_CALL)
    include(CMakeUninstallTarget)
    WRITE_UNINSTALL_TARGET_SCRIPT()
    add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif (NOT SUB_CMAKE_CALL)

# distclean 
if (UNIX AND (NOT SUB_CMAKE_CALL))
add_custom_target (distclean
    COMMAND echo "Cleaning all generated files"
    COMMAND make clean
    COMMAND rm -f "${CMAKE_CURRENT_BINARY_DIR}/src/include/pnl/pnl_config.h"
    COMMAND rm -f "${CMAKE_CURRENT_BINARY_DIR}/src/include/pnl/FC.h"
    COMMAND find "${CMAKE_CURRENT_BINARY_DIR}/" -name 'CMakeFiles' -type d | xargs rm -rf 
    COMMAND find "${CMAKE_CURRENT_BINARY_DIR}/" -name 'CTestTestfile.cmake' | xargs rm -f 
    COMMAND find "${CMAKE_CURRENT_BINARY_DIR}/" -name 'cmake_install.cmake' | xargs rm -f 
    COMMAND find "${CMAKE_CURRENT_BINARY_DIR}/" -name 'cmake_uninstall.cmake' | xargs rm -f 
    COMMAND rm -f "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt" 
    COMMAND rm -f "${CMAKE_CURRENT_BINARY_DIR}/CMakeuser.incl" 
    COMMAND rm -f "${CMAKE_CURRENT_BINARY_DIR}/install_manifest.txt"
    COMMAND rm -f "${CMAKE_CURRENT_BINARY_DIR}/libpnl.la"
    COMMAND find "${CMAKE_CURRENT_BINARY_DIR}/src" -name 'Makefile' | xargs rm -f 
    )
endif (UNIX AND (NOT SUB_CMAKE_CALL))

## Output configuration stuff ##
message("C Compiler :             ${CMAKE_C_COMPILER}")
message("     FLAGS :             ${PNL_C_FLAGS}")
if (CMAKE_CXX_COMPILER_WORKS)
    message("C++ Compiler :           ${CMAKE_CXX_COMPILER}")
    message("     FLAGS   :           ${PNL_CXX_FLAGS}")
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
