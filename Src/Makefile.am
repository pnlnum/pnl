SUBDIRS=common mod opt libpnl include

include $(top_srcdir)/Make.incl
AM_CPPFLAGS = $(PNL_INCLUDES) $(PREMIA_MINIMUM_INCLUDES)

lib_LTLIBRARIES = libpremia.la
libpremia_la_SOURCES=premia_obj.c
libpremia_la_LDFLAGS=-avoid-version
libpremia_la_LIBADD = common/libcommon.la mod/libmod.la opt/libopt.la libpnl/src/libpnl.la -lstdc++ 


bin_PROGRAMS=PremiaTerm
PremiaTerm_SOURCES=premia.c
PremiaTerm_LDADD= libpremia.la
nodist_EXTRA_PremiaTerm_SOURCES=dummy.cpp

include_HEADERS=premia_obj.h config.h configunix.h

all-local: installpnl libpremia.la

installpnl:
	( cd libpnl; make install )

## this file must be compiled without optimization
premia.o: premia.c
	$(CC) $(DEFS) $(DEFAULT_INCLUDES) $(PREMIA_MINIMUM_INCLUDES) $(PNL_INCLUDES) -c -Wall -g -o $@ $<

config.h:
	$(LN_S) $(top_srcdir)/Src/config.h .

