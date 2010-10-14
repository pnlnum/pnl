include $(top_srcdir)/Make.incl

AM_CPPFLAGS = $(PNL_INCLUDES)
some_extra_libs=@LAPACK_STATIC_LIBS@ @BLAS_STATIC_LIBS@

bin_PROGRAMS=testpnl

testpnl_SOURCES=cumulfunc_test.c  main.c  vector_test.c matrix_test.c \
random_test.c  speed_matrix.c  \
solver_iterativ_test.c basis_test.c dft_test.c finance_function_test.c \
complex_test.c  mathtools_test.c bessel_test.c integration_test.c \
root_test.c band_matrix_test.c special_func_test.c \
laplace_test.c lapack_test.c tridiag_test.c optim_intpoints_bfgs_test.c \
list_test.c

testpnl_LDADD=../src/libpnl.la ${some_extra_libs}
testpnl_LDFLAGS=@GSL_LIBS@ @BLAS_LIBS@ @LAPACK_LIBS@
testpnl_CPPFLAGS=@HAVE_GSL@ $(AM_CPPFLAGS)

if HAVE_MPI
bin_PROGRAMS+=pnl-mpi rng-mpi

pnl_mpi_SOURCES=mpibindings_test.c
pnl_mpi_CFLAGS=@MPI_CFLAGS@ $(AM_CFLAGS)
pnl_mpi_LDADD=../src/libpnl.la ${some_extra_libs}
pnl_mpi_LDFLAGS=@BLAS_LIBS@ @LAPACK_LIBS@

rng_mpi_SOURCES=rng-mpi.c
rng_mpi_CFLAGS=@MPI_CFLAGS@ $(AM_CFLAGS)
rng_mpi_LDADD=../src/libpnl.la ${some_extra_libs}
rng_mpi_LDFLAGS=@BLAS_LIBS@ @LAPACK_LIBS@
endif

