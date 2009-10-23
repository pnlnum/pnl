include $(top_srcdir)/Make.incl

AM_CPPFLAGS = $(PNL_INCLUDES)

bin_PROGRAMS=testpnl

testpnl_SOURCES=cumulfunc_test.c  main.c  vector_test.c matrix_test.c \
random_test.c  speed_matrix.c special_matrix_test.c perm_test.c \
solver_iterativ_test.c basis_test.c dft_test.c finance_function_test.c \
complex_test.c  mathtools_test.c bessel_test.c integration_test.c \
list_test.c root_test.c band_matrix_test.c special_func_test.c \
laplace_test.c lapack_test.c

testpnl_LDADD=../src/libpnl.la

testpnl_LDFLAGS=@BLAS_LIBS@ @LAPACK_LIBS@ @GSL_LIBS@ 
