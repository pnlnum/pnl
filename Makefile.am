SUBDIRS= . Src

datarootdir     = @datarootdir@
pkgconfigdir = @pkgconfigdir@

all-local: copy-configh

copy-configh:
	(cd Src; \
	test -r config.h || $(LN_S) $(top_srcdir)/../Src/config.h config.h)

install-pkgconfig:
	test -z "$(DESTDIR)$(pkgconfigdir)" || $(mkdir_p) "$(DESTDIR)$(pkgconfigdir)"	
	$(INSTALL_DATA) premia.pc "$(DESTDIR)$(pkgconfigdir)"
