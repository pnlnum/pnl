include_directories(${PROJECT_BINARY_DIR}/src/include ${PROJECT_SOURCE_DIR}/src/include)

set(SUBDIRS 
    objects
    libamos
    librand
    fft
    math
    linalg
    sort
    specfun
    libcephes
    libinteg
    interpol
    roots
    libminpack
    optim)

if (NOT BLAS_FOUND OR NOT LAPACK_FOUND)
    set(SUBDIRS ${SUBDIRS} libf2c libblas liblapack)
endif (NOT BLAS_FOUND OR NOT LAPACK_FOUND)

if (MPI_C_FOUND)
    set(SUBDIRS ${SUBDIRS} mpi)
endif (MPI_C_FOUND)

set(CMAKE_POSITION_INDEPENDENT_CODE ON) 

SET(pnl_LIB_SRCS
   cumulfunc.c
   finance_function.c
   pnl.def
)
install(FILES pnl.def DESTINATION lib)

add_subdirectory(include/pnl)
foreach(subdir ${SUBDIRS})
    add_subdirectory (${subdir})
    set (pnl_LIB_SRCS ${pnl_LIB_SRCS} $<TARGET_OBJECTS:${subdir}>)
endforeach ()

# Add build inline stuff
set (pnl_LIB_SRCS ${pnl_LIB_SRCS} include/pnl/build-inline.c)


add_library(pnl ${LINK_TYPE} ${pnl_LIB_SRCS})
target_link_libraries(pnl ${LIBS})

if (MINGW)
    set_target_properties(pnl PROPERTIES 
        LINK_FLAGS "-Wl,--output-def,${CMAKE_CURRENT_BINARY_DIR}/libpnl.def")
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libpnl.def DESTINATION lib)
endif (MINGW)

set_target_properties(pnl PROPERTIES 
    SOVERSION ${PNL_VERSION}
    VERSION ${PNL_VERSION})
install(TARGETS pnl DESTINATION lib)

set(PNL_INCLUDE_DIRS)
if (MPI_C_FOUND)
    set(PNL_INCLUDE_DIRS ${MPI_C_INCLUDE_PATH})
endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/PnlConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/PnlConfig.cmake @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/PnlConfig.cmake DESTINATION lib/pnl)
