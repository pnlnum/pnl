
do: doxy

# We use a workaround because doxygen only expands macros in the first
# template file! Seems to be a bug but not yet solved...
doxy:
	(cd ../src; \
	cat linalg/pnl_TMatrix.c linalg/pnl_TVector.c > linalg/pnl_Talltemplates.c; \
	doxygen ../doc/libpnldoxy.rc ; \
	rm -f linalg/pnl_Talltemplates.c; )
