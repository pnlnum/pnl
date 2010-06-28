SUBDIRS = 

include $(top_srcdir)/Make.incl
AM_CPPFLAGS = $(PNL_INCLUDES)

noinst_LTLIBRARIES = libmpi.la
libmpi_la_SOURCES = mpi_ops.c
libmpi_la_CFLAGS = @MPI_CFLAGS@
libmpi_la_LDFLAGS=-avoid-version @MPI_LIBS@

