include $(top_srcdir)/Make.incl

AM_CPPFLAGS = $(PNL_INCLUDES)
some_extra_libs=@LAPACK_STATIC_LIBS@ @BLAS_STATIC_LIBS@

check_PROGRAMS=         band_matrix_test    basis_test     \
complex_test            cumulfunc_test      dft_test       \
finance_function_test   integration_test    lapack_test    \
laplace_test            list_test           mathtools_test \
ode_test                optim_test          random_test    \
root_test               saveload_test       matrix_test    \
solver_iterativ_test    special_func_test   speed_test     \
tridiag_test            vector_test

TESTS=$(check_PROGRAMS) 
build-tests: $(TESTS)

band_matrix_test_SOURCES=band_matrix_test.c
band_matrix_test_LDADD=../src/libpnl.la ${some_extra_libs}

basis_test_SOURCES=basis_test.c tests_utils.c
basis_test_LDADD=../src/libpnl.la ${some_extra_libs}

complex_test_SOURCES=complex_test.c tests_utils.c
complex_test_LDADD=../src/libpnl.la ${some_extra_libs}

cumulfunc_test_SOURCES=cumulfunc_test.c
cumulfunc_test_LDADD=../src/libpnl.la ${some_extra_libs}

dft_test_SOURCES=dft_test.c tests_utils.c
dft_test_LDADD=../src/libpnl.la ${some_extra_libs}

finance_function_test_SOURCES=finance_function_test.c tests_utils.c
finance_function_test_LDADD=../src/libpnl.la ${some_extra_libs}

integration_test_SOURCES=integration_test.c tests_utils.c
integration_test_LDADD=../src/libpnl.la ${some_extra_libs}

lapack_test_SOURCES=lapack_test.c
lapack_test_LDADD=../src/libpnl.la ${some_extra_libs}

laplace_test_SOURCES=laplace_test.c tests_utils.c
laplace_test_LDADD=../src/libpnl.la ${some_extra_libs}

list_test_SOURCES=list_test.c
list_test_LDADD=../src/libpnl.la ${some_extra_libs}

mathtools_test_SOURCES=mathtools_test.c tests_utils.c
mathtools_test_LDADD=../src/libpnl.la ${some_extra_libs}

matrix_test_SOURCES=matrix_test.c tests_utils.c
matrix_test_LDADD=../src/libpnl.la ${some_extra_libs}

ode_test_SOURCES=ode_test.c
ode_test_LDADD=../src/libpnl.la ${some_extra_libs}

optim_test_SOURCES=optim_test.c
optim_test_LDADD=../src/libpnl.la ${some_extra_libs}

random_test_SOURCES=random_test.c
random_test_LDADD=../src/libpnl.la ${some_extra_libs}

root_test_SOURCES=root_test.c
root_test_LDADD=../src/libpnl.la ${some_extra_libs}

solver_iterativ_test_SOURCES=solver_iterativ_test.c
solver_iterativ_test_LDADD=../src/libpnl.la ${some_extra_libs}

special_func_test_SOURCES=special_func_test.c
special_func_test_LDADD=../src/libpnl.la ${some_extra_libs}

speed_test_SOURCES=speed_test.c speed_matrix.c speed_basis.c tests_utils.c
speed_test_LDADD=../src/libpnl.la ${some_extra_libs}

tridiag_test_SOURCES=tridiag_test.c
tridiag_test_LDADD=../src/libpnl.la ${some_extra_libs}

vector_test_SOURCES=vector_test.c tests_utils.c
vector_test_LDADD=../src/libpnl.la ${some_extra_libs}

if HAVE_MPI
check_PROGRAMS+=pnl-mpi_test rng-mpi_test saveload_test

pnl_mpi_test_SOURCES=mpibindings_test.c
pnl_mpi_test_CFLAGS=@MPI_CFLAGS@ $(AM_CFLAGS)
pnl_mpi_test_LDADD=../src/libpnl.la ${some_extra_libs}

rng_mpi_test_SOURCES=rng-mpi.c
rng_mpi_test_CFLAGS=@MPI_CFLAGS@ $(AM_CFLAGS)
rng_mpi_test_LDADD=../src/libpnl.la ${some_extra_libs}

saveload_test_SOURCES=saveload_test.c
saveload_test_CFLAGS=@MPI_CFLAGS@ $(AM_CFLAGS)
saveload_test_LDADD=../src/libpnl.la ${some_extra_libs}

endif


