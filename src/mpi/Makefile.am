SUBDIRS = 

include $(top_srcdir)/Make.incl
AM_CPPFLAGS = $(PNL_INCLUDES)

noinst_LTLIBRARIES = libpnlmpi.la
libpnlmpi_la_SOURCES = mpi_ops.c
libpnlmpi_la_CFLAGS = @MPI_CFLAGS@ $(AM_CFLAGS)
libpnlmpi_la_LDFLAGS=-avoid-version @MPI_LIBS@

