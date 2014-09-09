# - Find the Pnl Library
#
# Usage:
#   find_package(Pnl [REQUIRED] [QUIET])
#
# It sets the following variables:
#   PNL_FOUND           ... true if Pnl is found on the sytem
#   PNL_LIBRARIES       ... full path to the Pnl library
#   PNL_INCLUDE_DIRS    ... Pnl include directory
#
# The following variables will be checked by the function
#   PNL_ROOT            ... if set, headers and libraries are exclusively
#                           searched under this path
#   PNL_LIBRARY         ... pnl library path to use
#   PNL_INCLUDE_DIR     ... pnl include directory

if (PNL_ROOT)
    find_library(
        PNL_LIBRARY
        NAMES "pnl"
        PATHS ${PNL_ROOT}
        PATH_SUFFIXES "lib" "lib64"
        NO_DEFAULT_PATH
        )

    find_path(
        PNL_INCLUDE_DIR
        NAMES "pnl/pnl_random.h"
        PATHS ${PNL_ROOT}
        PATH_SUFFIXES "include"
        NO_DEFAULT_PATH
        )
else (PNL_ROOT)
    find_library(
        PNL_LIBRARY
        NAMES "pnl"
        PATH_SUFFIXES "lib" "lib64"
        NO_DEFAULT_PATH
        )

    find_path(
        PNL_INCLUDE_DIR
        NAMES "pnl/pnl_random.h"
        PATH_SUFFIXES "include"
        NO_DEFAULT_PATH
        )
endif (PNL_ROOT)

# Handle the QUIETLY and REQUIRED arguments and set PNL_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PNL 
    REQUIRED_VARS PNL_LIBRARY PNL_INCLUDE_DIR)
if (PNL_FOUND)
    set(PNL_INCLUDE_DIRS ${PNL_INCLUDE_DIR})
    set(PNL_LIBRARIES ${PNL_LIBRARY})
endif (PNL_FOUND)
mark_as_advanced(PNL_INCLUDE_DIRS PNL_LIBRARIES)
