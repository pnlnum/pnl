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

add_subdirectory(include/pnl)

if (NOT BLAS_FOUND OR NOT LAPACK_FOUND)
    set(SUBDIRS ${SUBDIRS} libf2f libblas liblapack)
endif (NOT BLAS_FOUND OR NOT LAPACK_FOUND)

if (MPI_C_FOUND)
    set(SUBDIRS ${SUBDIRS} mpi)
endif (MPI_C_FOUND)



SET(pnl_LIB_SRCS
   cumulfunc.c
   finance_function.c
)

foreach(subdir ${SUBDIRS})
    include(${subdir}/CMakeLists.txt)
    foreach (srcfile ${${subdir}_SRCS})
        set (pnl_LIB_SRCS ${pnl_LIB_SRCS} ${subdir}/${srcfile})
    endforeach ()
endforeach ()

# Add build inline stuff
set (pnl_LIB_SRCS ${pnl_LIB_SRCS} include/pnl/build-inline.c)


add_library(pnl SHARED ${pnl_LIB_SRCS})
target_link_libraries(pnl ${LIBS})

set_target_properties(pnl PROPERTIES 
    LT_VERSION_CURRENT ${LT_VERSION_CURRENT}
    LT_VERSION_AGE ${LT_VERSION_AGE}
    LT_VERSION_REVISION ${LT_VERSION_REVISION}
    LT_DEPENDENCY_LIBS ${PNL_LDFLAGS}
    SOVERSION ${LT_VERSION_CURRENT}
    VERSION ${LT_VERSION_CURRENT}.${LT_VERSION_AGE}.${LT_VERSION_REVISION}

)
CREATE_LIBTOOL_FILE(pnl /lib)
install(TARGETS pnl DESTINATION lib)

