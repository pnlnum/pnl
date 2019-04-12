#!/bin/bash

# Example
# Transform
#   int pnl_vect_complex_isequal_rel(const PnlVectComplex *x, const PnlVectComplex *y, double relerr)
# into
#   item \describefun{int}{pnl_vect_complex_isequal_rel}{const \PnlVectComplex \ptr x, const \PnlVectComplex \ptr y, double relerr} 

sed -E 's/ *extern *//g' | \
sed -E 's/([a-zA-Z0-9_ \*]+) +([a-zA-Z0-9_]+) *\((.*)\);?/\\item \\describefun{\1}{\2}{\3}/g' | \
sed -E 's/\*/\\ptr /g' | \
sed -E 's/(Pnl[a-zA-Z]+)/\\\1/g' | \
sed -E 's/dcomplex/\\dcomplex/g'
