add_library(linalg OBJECT
   vector_double.c
   matrix.c
   vector_int.c
   vector_complex.c
   tridiag_matrix.c
   band_matrix.c
   pnl_TVector.c
   pnl_TMatrix.c
   perm.c
   solver_iterativ.c
   lapack-wrap.c
   dpstrf.c
   zpstrf.c
)
