SUBDIRS = 

include $(top_srcdir)/Make.incl
AM_CPPFLAGS = $(PNL_INCLUDES)

noinst_LTLIBRARIES = liblinalg.la
liblinalg_la_SOURCES = \
vector_double.c matrix.c vector_int.c vector_complex.c \
tridiag_matrix.c band_matrix.c pnl_TVector.c pnl_TMatrix.c perm.c \
solver_iterativ.c lapack-wrap.c dpstrf.c
liblinalg_la_LDFLAGS=-avoid-version

noinst_HEADERS = pnl_TVector_source.c pnl_TMatrix_source.c 


