SUBDIRS = 

include $(top_srcdir)/Make.incl
AM_CPPFLAGS = $(PNL_INCLUDES)

noinst_LTLIBRARIES = liblinalg.la
liblinalg_la_SOURCES = matrix.c vector_uint.c vector_complex.c tridiag_matrix.c \
band_matrix.c special_matrix.c pnl_TVector.c pnl_TMatrix.c pnl_Tarray.c perm.c \
solver_iterativ.c lapack-wrap.c
liblinalg_la_LDFLAGS=-avoid-version	

noinst_HEADERS = pnl_TVector_source.c pnl_TMatrix_source.c


