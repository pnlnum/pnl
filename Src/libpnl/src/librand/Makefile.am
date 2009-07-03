SUBDIRS = 

include $(top_srcdir)/Make.incl
AM_CPPFLAGS = $(PNL_INCLUDES)

noinst_LTLIBRARIES = librand.la
librand_la_SOURCES = random.c     rk_mt.c  rk_primitive.c  rk_sobol.c
librand_la_LDFLAGS=-avoid-version	

