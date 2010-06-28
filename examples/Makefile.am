include $(top_srcdir)/Make.incl

AM_CPPFLAGS = $(PNL_INCLUDES)

bin_PROGRAMS=testpnl pnl-mpi

testpnl_SOURCES=cumulfunc_test.c  main.c  vector_test.c matrix_test.c \
random_test.c  speed_matrix.c  \
solver_iterativ_test.c basis_test.c dft_test.c finance_function_test.c \
complex_test.c  mathtools_test.c bessel_test.c integration_test.c \
root_test.c band_matrix_test.c special_func_test.c \
laplace_test.c lapack_test.c tridiag_test.c optim_intpoints_bfgs_test.c

testpnl_LDADD=../src/libpnl.la
testpnl_LDFLAGS=@GSL_LIBS@ 

pnl_mpi_SOURCES=mpibindings_test.c
pnl_mpi_LDADD=../src/libpnl.la
