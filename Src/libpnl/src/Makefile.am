SUBDIRS = libamos librand libcsparse libfftpack \
linalg sort include specfun libcephes

include $(top_srcdir)/Make.incl
AM_CPPFLAGS = $(PNL_INCLUDES)

lib_LTLIBRARIES = libpnl.la
libpnl_la_SOURCES = basis.c pnl_complex.c cumulfunc.c mathtools.c \
finance_function.c integration_GK.c \
brent.c roots.c integ.c laplace.c

libpnl_la_LDFLAGS=-avoid-version @BLAS_LIBS@ @LAPACK_LIBS@
libpnl_la_LIBADD=libamos/libamos.la librand/librand.la\
libcsparse/libcsparse.la linalg/liblinalg.la sort/libsort.la \
libfftpack/libfftpack.la specfun/libspecfun.la libcephes/libcephes.la

