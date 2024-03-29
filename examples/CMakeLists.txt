include_directories(${PROJECT_BINARY_DIR}/src/include ${PROJECT_SOURCE_DIR}/src/include)
enable_testing()
add_custom_target(link_data
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/Data ${CMAKE_CURRENT_BINARY_DIR}/Data
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/Data_specfun ${CMAKE_CURRENT_BINARY_DIR}/Data_specfun
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/Data_basis ${CMAKE_CURRENT_BINARY_DIR}/Data_basis
)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure)

SET(band_matrix_test_SRCS
   band_matrix_test.c
   tests_utils.c
)
add_executable(band_matrix_test ${band_matrix_test_SRCS})
target_link_libraries(band_matrix_test pnl)

SET(sp_matrix_test_SRCS
   spmatrix_test.c
   tests_utils.c
)
add_executable(sp_matrix_test ${sp_matrix_test_SRCS})
target_link_libraries(sp_matrix_test pnl)


SET(basis_test_SRCS
   basis_test.c
   tests_utils.c
)
add_executable(basis_test ${basis_test_SRCS})
target_link_libraries(basis_test pnl)


SET(cumulfunc_test_SRCS
   cumulfunc_test.c
   tests_utils.c
)
add_executable(cumulfunc_test ${cumulfunc_test_SRCS})
target_link_libraries(cumulfunc_test pnl)


SET(complex_test_SRCS
   complex_test.c
   tests_utils.c
)
add_executable(complex_test ${complex_test_SRCS})
target_link_libraries(complex_test pnl)


SET(fft_test_SRCS
   fft_test.c
   tests_utils.c
)
add_executable(fft_test ${fft_test_SRCS})
target_link_libraries(fft_test pnl)


SET(finance_function_test_SRCS
   finance_function_test.c
   tests_utils.c
)
add_executable(finance_function_test ${finance_function_test_SRCS})
target_link_libraries(finance_function_test pnl)


SET(integration_test_SRCS
   integration_test.c
   tests_utils.c
)
add_executable(integration_test ${integration_test_SRCS})
target_link_libraries(integration_test pnl)


SET(lapack_test_SRCS
   lapack_test.c
   tests_utils.c
)
add_executable(lapack_test ${lapack_test_SRCS})
target_link_libraries(lapack_test pnl)
add_dependencies(lapack_test link_data)


SET(laplace_test_SRCS
   laplace_test.c
   tests_utils.c
)
add_executable(laplace_test ${laplace_test_SRCS})
target_link_libraries(laplace_test pnl)


SET(list_test_SRCS
   list_test.c
   tests_utils.c
)
add_executable(list_test ${list_test_SRCS})
target_link_libraries(list_test pnl)


SET(mathtools_test_SRCS
   mathtools_test.c
   tests_utils.c
)
add_executable(mathtools_test ${mathtools_test_SRCS})
target_link_libraries(mathtools_test pnl)


SET(ode_test_SRCS
   ode_test.c
   tests_utils.c
)
add_executable(ode_test ${ode_test_SRCS})
target_link_libraries(ode_test pnl)

SET(optim_test_SRCS
   optim_test.c
   tests_utils.c
)
add_executable(optim_test ${optim_test_SRCS})
target_link_libraries(optim_test pnl)


SET(random_test_SRCS
   random_test.c
   tests_utils.c
)
add_executable(random_test ${random_test_SRCS})
target_link_libraries(random_test pnl)


SET(root_test_SRCS
   root_test.c
   tests_utils.c
)
add_executable(root_test ${root_test_SRCS})
target_link_libraries(root_test pnl)


SET(matrix_test_SRCS
   matrix_test.c
   tests_utils.c
)
add_executable(matrix_test ${matrix_test_SRCS})
target_link_libraries(matrix_test pnl)
add_dependencies(matrix_test link_data)


SET(solver_iterativ_test_SRCS
   solver_iterativ_test.c
   tests_utils.c
)
add_executable(solver_iterativ_test ${solver_iterativ_test_SRCS})
target_link_libraries(solver_iterativ_test pnl)
add_dependencies(solver_iterativ_test link_data)


SET(special_func_test_SRCS
   special_func_test.c
   tests_utils.c
)
add_executable(special_func_test ${special_func_test_SRCS})
target_link_libraries(special_func_test pnl)
add_dependencies(special_func_test link_data)


SET(tridiag_test_SRCS
   tridiag_test.c
   tests_utils.c
)
add_executable(tridiag_test ${tridiag_test_SRCS})
target_link_libraries(tridiag_test pnl)


SET(vector_test_SRCS
   vector_test.c
   tests_utils.c
)
add_executable(vector_test ${vector_test_SRCS})
target_link_libraries(vector_test pnl)


SET(hmat_test_SRCS
   hmat_test.c
   tests_utils.c
)
add_executable(hmat_test ${hmat_test_SRCS})
target_link_libraries(hmat_test pnl)


SET(array_test_SRCS
   array_test.c
   tests_utils.c
)
add_executable(array_test ${array_test_SRCS})
target_link_libraries(array_test pnl)

if (PNL_ENABLE_TESTS)
    add_test(band_matrix_test ${CMAKE_CURRENT_BINARY_DIR}/band_matrix_test)
    add_test(sp_matrix_test ${CMAKE_CURRENT_BINARY_DIR}/sp_matrix_test)
    add_test(basis_test ${CMAKE_CURRENT_BINARY_DIR}/basis_test)
    add_test(cumulfunc_test ${CMAKE_CURRENT_BINARY_DIR}/cumulfunc_test)
    add_test(complex_test ${CMAKE_CURRENT_BINARY_DIR}/complex_test)
    add_test(fft_test ${CMAKE_CURRENT_BINARY_DIR}/fft_test)
    add_test(finance_function_test ${CMAKE_CURRENT_BINARY_DIR}/finance_function_test)
    add_test(integration_test ${CMAKE_CURRENT_BINARY_DIR}/integration_test)
    add_test(lapack_test ${CMAKE_CURRENT_BINARY_DIR}/lapack_test)
    add_test(laplace_test ${CMAKE_CURRENT_BINARY_DIR}/laplace_test)
    add_test(list_test ${CMAKE_CURRENT_BINARY_DIR}/list_test)
    add_test(mathtools_test ${CMAKE_CURRENT_BINARY_DIR}/mathtools_test)
    add_test(ode_test ${CMAKE_CURRENT_BINARY_DIR}/ode_test)
    add_test(optim_test ${CMAKE_CURRENT_BINARY_DIR}/optim_test)
    add_test(random_test ${CMAKE_CURRENT_BINARY_DIR}/random_test)
    add_test(root_test ${CMAKE_CURRENT_BINARY_DIR}/root_test)
    add_test(matrix_test ${CMAKE_CURRENT_BINARY_DIR}/matrix_test)
    add_test(solver_iterativ_test ${CMAKE_CURRENT_BINARY_DIR}/solver_iterativ_test)
    add_test(special_func_test ${CMAKE_CURRENT_BINARY_DIR}/special_func_test)
    add_test(tridiag_test ${CMAKE_CURRENT_BINARY_DIR}/tridiag_test)
    add_test(vector_test ${CMAKE_CURRENT_BINARY_DIR}/vector_test)
    add_test(hmat_test ${CMAKE_CURRENT_BINARY_DIR}/hmat_test)
    add_test(array_test ${CMAKE_CURRENT_BINARY_DIR}/array_test)
endif (PNL_ENABLE_TESTS)

if (NOT USE_INTERNAL_BLAS)
    SET(speed_test_SRCS
        speed_test.c
        speed_matrix.c
        speed_basis.c
        tests_utils.c
        )
    add_executable(speed_test ${speed_test_SRCS})
    target_link_libraries(speed_test pnl)
endif ()

if (MPI_C_FOUND)
    SET(mpibindings_test_SRCS
        mpibindings_test.c
        tests_utils.c
    )
    add_executable(mpibindings_test ${mpibindings_test_SRCS})
    target_link_libraries(mpibindings_test pnl)

    SET(rng_mpi_test_SRCS
        rng-mpi.c
        tests_utils.c
    )
    add_executable(rng_mpi_test ${rng_mpi_test_SRCS})
    target_link_libraries(rng_mpi_test pnl)

    SET(saveload_test_SRCS
        saveload_test.c
        tests_utils.c
    )
    add_executable(saveload_test ${saveload_test_SRCS})
    target_link_libraries(saveload_test pnl)
endif (MPI_C_FOUND)


