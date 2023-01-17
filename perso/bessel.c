#include "pnl/pnl_specfun.h"
#include "pnl/pnl_matrix.h"

int main()
{
  PnlMatComplex *A = pnl_mat_complex_create_from_scalar(2, 3, Complex(2, 3));
  PnlMatComplex *B = pnl_mat_complex_create(2, 3);
  MLET_COMPLEX(B, 1, 2) = MGET_COMPLEX(A, 1, 1);
  /* double lambda = 2.3; */
  /* double x = 1560; */
  /* pnl_deactivate_mtherr(); */
  /* double result=pnl_bessel_i(lambda, x); */
  return 0;
}
