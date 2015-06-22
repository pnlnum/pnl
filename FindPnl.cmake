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
        NAMES "pnl" "libpnl"
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
        NAMES "pnl" "libpnl"
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

# if (PNL_LIBRARY AND PNL_INCLUDE_DIR)
#     set(PNL_FOUND TRUE)
#     set(PNL_INCLUDE_DIRS ${PNL_INCLUDE_DIR})
#     set(PNL_LIBRARIES ${PNL_LIBRARY})
# else (PNL_LIBRARY AND PNL_INCLUDE_DIR)
#     set(PNL_FOUND FALSE)
#     message(STATUS
#         "PNL not found. Consider specifying either (PNL_LIBRARY, PNL_INCLUDE_DIR) or PNL_ROOT."
#         )
# endif (PNL_LIBRARY AND PNL_INCLUDE_DIR)


# Handle the QUIETLY and REQUIRED arguments and set PNL_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Pnl
    REQUIRED_VARS 
        PNL_LIBRARY PNL_INCLUDE_DIR
    FAIL_MESSAGE 
        "Pnl not found. Specify either (PNL_LIBRARY, PNL_INCLUDE_DIR) or PNL_ROOT"
)
if (PNL_FOUND)
    set(PNL_INCLUDE_DIRS ${PNL_INCLUDE_DIR})
    set(PNL_LIBRARIES ${PNL_LIBRARY})
    message(STATUS "PNL Include: ${PNL_INCLUDE_DIRS}")
    message(STATUS "PNL Libraries: ${PNL_LIBRARIES}")
endif (PNL_FOUND)
mark_as_advanced(PNL_INCLUDE_DIRS PNL_LIBRARIES)
