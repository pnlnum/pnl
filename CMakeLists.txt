cmake_minimum_required(VERSION 2.8.9)
cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules;${CMAKE_MODULE_PATH}")
set(PROJECT_LIBRARY_NAME ${PROJECT_NAME})
file(STRINGS "VERSION" PNL_VERSION)
set(BUILD_SHARED_LIBS ON)


# Prefix 
if(PNL_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX ${PNL_INSTALL_PREFIX} CACHE STRING "Installation prefix." FORCE)
else(PNL_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR} CACHE STRING "Installation prefix." FORCE)
endif(PNL_INSTALL_PREFIX)

# Compilers
if (CC)
    set (CMAKE_C_COMPILER ${CC})
endif (CC)
if (CXX)
    set (CMAKE_CXX_COMPILER ${CXX})
endif (CXX)


if(CROSS_COMPILE)
    include(cross-compile)
    set(MINGW ON)
endif()

project(pnl C CXX)

if (NOT LINK_TYPE) 
    set(LINK_TYPE SHARED)
endif (NOT LINK_TYPE)

set(CMAKE_MACOSX_RPATH ON)
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(pnl_INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib" )
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_FILES_DIRECTORY "/CMakeFiles")


option(USE_MPI "Compiles the MPI bindings for PNL. Default = on" ON)
option(WITHOUT_DEPRECATED "Do not build deprecated functions. Default = off" OFF)
option(USE_INTERNAL_BLAS "Use internal Blas & Lapack instead of the system's ones. Default = off" OFF)


# Release or Debug
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif ()



enable_language(CXX OPTIONAL)

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
    # enable_language(Fortran OPTIONAL)
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


if (USE_MPI)
    set (MPI_FIND_QUIETLY true)
    find_package(MPI)
    if (MPI_C_FOUND)
        set(LIBS ${LIBS} ${MPI_C_LIBRARIES})
        include_directories(${MPI_C_INCLUDE_PATH})
        set(PNL_CPP_FLAGS "${PNL_CPP_FLAGS} -I${MPI_C_INCLUDE_PATH}")
    endif (MPI_C_FOUND)
endif (USE_MPI)

# Add -lm and -ldl if needed
if (UNIX)
    set(LIBS ${LIBS} m dl)
endif (UNIX)

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
