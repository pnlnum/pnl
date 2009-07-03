all-local: copy-libpnlheaders copy-premiaheaders

copy-libpnlheaders:
	LIBPNLHEADERS="$(top_srcdir)/Src/libpnl/src/include/*.h"; \
	for h in $$LIBPNLHEADERS; do\
		basefile=`basename $$h`; \
		test -r $$basefile || $(LN_S) $$h $$basefile; \
	done

copy-premiaheaders:
	for h in var.h tools.h premia_types.h optype.h error_msg.h; do \
		test -r $$h || $(LN_S) $(top_srcdir)/Src/common/$$h $$basefile; \
	done
	test -r premia_obj.h || $(LN_S) $(top_srcdir)/Src/premia_obj.h premia_obj.h

clean:
	$(RM) -f *.h

cleanall: clean
