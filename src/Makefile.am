SUBDIRS = objects libamos librand fft math \
linalg sort include/pnl specfun libcephes \
libinteg interpol roots libminpack optim


include $(top_srcdir)/Make.incl
AM_CPPFLAGS = $(PNL_INCLUDES)

lib_LTLIBRARIES = libpnl.la
libpnl_la_SOURCES =  cumulfunc.c finance_function.c

libpnl_la_LDFLAGS=-avoid-version @BLAS_LIBS@ @LAPACK_LIBS@ @STATIC_FLAG@
libpnl_la_LIBADD=libamos/libamos.la librand/librand.la\
math/libmath.la linalg/liblinalg.la sort/libsort.la \
fft/libfft.la specfun/libspecfun.la libcephes/libcephes.la \
libinteg/libinteg.la interpol/libinterpol.la optim/liboptim.la\
roots/libroots.la libminpack/libminpack.la \
objects/libobjects.la include/pnl/libinline.la


if PNL_HAVE_MPI
SUBDIRS+=mpi
libpnl_la_LIBADD+=mpi/libpnlmpi.la
endif

if !PNL_HAVE_BLASLAPACK
SUBDIRS+=libf2c libblas liblapack
libpnl_la_LIBADD+=libf2c/libf2c.la libblas/libblas.la liblapack/liblapack.la
endif

emacs-tags:
	find . \( -name "*.[chCH]" -o -name "*.cpp" \) -print | etags -
