SUBDIRS = 

include $(top_srcdir)/Make.incl
AM_CPPFLAGS = $(PNL_INCLUDES)

noinst_LTLIBRARIES = libfft.la
libfft_la_SOURCES = fftpack.c fftpack_wrap.c laplace.c
libfft_la_LDFLAGS=-avoid-version	

