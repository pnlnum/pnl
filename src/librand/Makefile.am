SUBDIRS = 

include $(top_srcdir)/Make.incl
AM_CPPFLAGS = $(PNL_INCLUDES)

noinst_LTLIBRARIES = librand.la
librand_la_SOURCES = random.c  distributions.c  mt.c dcmt.c \
					 sobol_i4.c sobol_i8.c
librand_la_LDFLAGS=-avoid-version	

# Do not compile sobol_i4.c with optimization flags
sobol_i4.lo: sobol_i4.c
	$(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(CC) \
		$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) -Wall \
		-O0 -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ $< 
	mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Plo

# Do not compile sobol_i8.c with optimization flags
sobol_i8.lo: sobol_i8.c
	$(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(CC) \
		$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) -Wall \
		-O0 -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ $< 
	mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Plo
