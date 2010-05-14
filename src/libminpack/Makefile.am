SUBDIRS = 

include $(top_srcdir)/Make.incl
AM_CPPFLAGS = $(PNL_INCLUDES)

noinst_LTLIBRARIES = libminpack.la
libminpack_la_SOURCES =chkder.c  enorm.c   hybrj.c   qrfac.c   r1updt.c \
dogleg.c  fdjac1.c  hybrd.c   lmdif.c   qrsolv.c rwupdt.c dpmpar.c \
fdjac2.c  lmder.c   lmpar.c   qform.c   r1mpyq.c covar.c  minpack-wrap.c

libminpack_la_LDFLAGS=-avoid-version	

