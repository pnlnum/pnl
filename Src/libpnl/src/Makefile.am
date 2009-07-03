SUBDIRS = libamos libcephes librand libcsparse libfftpack \
linalg sort include

include $(top_srcdir)/Make.incl
AM_CPPFLAGS = $(PNL_INCLUDES)

lib_LTLIBRARIES = libpnl.la
libpnl_la_SOURCES = basis.c pnl_complex.c cumulfunc.c mathtools.c \
finance_function.c pde_tools.c integration_GK.c \
sparse_grid_constructor.c gridsparse_constructor.c gridsparse_functions.c \
brent.c roots.c special_func_integral.c integ.c laplace.c

libpnl_la_LDFLAGS=-avoid-version	
libpnl_la_LIBADD=libamos/libamos.la libcephes/libcephes.la librand/librand.la\
libcsparse/libcsparse.la linalg/liblinalg.la sort/libsort.la \
libfftpack/libfftpack.la


