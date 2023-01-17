#include "pnl/pnl_basis.h"

int main()
{
  int dim = 5;
  int deg = 4;
  double q = 0.7;

  PnlBasis *B;

  B = pnl_basis_create_from_degree(PNL_BASIS_CANONICAL, deg, dim);
  printf("rows %d, cols %d\n", B->T->m, B->T->n);
  pnl_mat_int_print(B->T);
  pnl_basis_free(&B);
  printf("\n\n");
  B = pnl_basis_create_from_prod_degree(PNL_BASIS_CANONICAL, deg, dim);
  printf("rows %d, cols %d\n", B->T->m, B->T->n);
  pnl_mat_int_print(B->T);
  pnl_basis_free(&B);
  printf("\n\n");
  B = pnl_basis_create_from_hyperbolic_degree(PNL_BASIS_CANONICAL, deg, q, dim);
  printf("rows %d, cols %d\n", B->T->m, B->T->n);
  pnl_mat_int_print(B->T);
  pnl_basis_free(&B);
  printf("\n");

  return 0;
}
