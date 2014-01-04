
do: doxy

# We use a workaround because doxygen only expands macros in the first
# template file! Seems to be a bug but not yet solved...
doxy:
	(cd ../src; \
	cat linalg/matrix.c linalg/vector.c linalg/sp_matrix.c  > linalg/pnl_Talltemplates.c; \
	doxygen ../doc/libpnldoxy.rc ; \
	rm -f linalg/pnl_Talltemplates.c; )
