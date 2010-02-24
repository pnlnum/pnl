SUBDIRS = libamos librand libcsparse libfftpack \
linalg sort include specfun libcephes \
@INTERNAL_BLAS_DIR@ @INTERNAL_LAPACK_DIR@ @INTERNAL_F2C_DIR@

include $(top_srcdir)/Make.incl
AM_CPPFLAGS = $(PNL_INCLUDES)

lib_LTLIBRARIES = libpnl.la
libpnl_la_SOURCES = infnan.c pnl_complex.c cumulfunc.c mathtools.c \
basis.c finance_function.c integration_GK.c \
brent.c roots.c integ.c laplace.c interpolation.c

libpnl_la_LDFLAGS=-avoid-version @BLAS_LIBS@ @LAPACK_LIBS@
libpnl_la_LIBADD=libamos/libamos.la librand/librand.la\
libcsparse/libcsparse.la linalg/liblinalg.la sort/libsort.la \
libfftpack/libfftpack.la specfun/libspecfun.la libcephes/libcephes.la \
@INTERNAL_BLAS_LA@ @INTERNAL_LAPACK_LA@ @INTERNAL_F2C_LA@


