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
)

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

if (UNIX)
    set_target_properties(pnl PROPERTIES 
        SOVERSION ${VERSION_CURRENT}
        VERSION ${VERSION_CURRENT}.${VERSION_AGE}.${VERSION_REVISION}
    )
endif(UNIX)
install(TARGETS pnl DESTINATION lib)

