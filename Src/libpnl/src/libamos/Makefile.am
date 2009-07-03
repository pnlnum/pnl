SUBDIRS =

include $(top_srcdir)/Make.incl
AM_CPPFLAGS = $(PNL_INCLUDES)

noinst_LTLIBRARIES = libamos.la libdlamch.la
libamos_la_SOURCES = \
amos_wrappers.c  xerror.c  zairy.c  zbesj.c  zbiry.c  zdiv.c   zmlri.c  zseri.c  zunhj.c  zunk1.c\
d1mach.c         zabs.c    zasyi.c  zbesk.c  zbknu.c  zexp.c   zmlt.c   zshch.c  zuni1.c  zunk2.c\
dgamln.c         zacai.c   zbesh.c  zbesy.c  zbuni.c  zkscl.c  zrati.c  zsqrt.c  zuni2.c  zuoik.c\
fdump.c          zacon.c   zbesi.c  zbinu.c  zbunk.c  zlog.c   zs1s2.c  zuchk.c  zunik.c  zwrsk.c

libdlamch_la_SOURCES = dlamch.c
libdlamch_la_CFLAGS = -O0
libamos_la_LDFLAGS=-avoid-version	
libamos_la_LIBADD = libdlamch.la
