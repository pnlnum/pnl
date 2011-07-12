SUBDIRS = 

include $(top_srcdir)/Make.incl
AM_CPPFLAGS = $(PNL_INCLUDES)

noinst_LTLIBRARIES = librand.la
librand_la_SOURCES = random.c  distributions.c  mt.c dcmt.c \
					 sobol_i4.c sobol_i8.c
librand_la_LDFLAGS=-avoid-version	

