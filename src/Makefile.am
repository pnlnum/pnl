SUBDIRS = libamos librand fft math \
linalg sort include specfun libcephes \
libinteg interpol roots \
@INTERNAL_BLAS_DIR@ @INTERNAL_LAPACK_DIR@ @INTERNAL_F2C_DIR@

include $(top_srcdir)/Make.incl
AM_CPPFLAGS = $(PNL_INCLUDES)

lib_LTLIBRARIES = libpnl.la
libpnl_la_SOURCES =  cumulfunc.c \
finance_function.c object.c

libpnl_la_LDFLAGS=-avoid-version @BLAS_LIBS@ @LAPACK_LIBS@
libpnl_la_LIBADD=libamos/libamos.la librand/librand.la\
math/libmath.la linalg/liblinalg.la sort/libsort.la \
fft/libfft.la specfun/libspecfun.la libcephes/libcephes.la \
libinteg/libinteg.la interpol/libinterpol.la\
roots/libroots.la \
@INTERNAL_BLAS_LA@ @INTERNAL_LAPACK_LA@ @INTERNAL_F2C_LA@


