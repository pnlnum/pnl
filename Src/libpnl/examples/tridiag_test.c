
/************************************************************************/
/* Copyright Jérôme Lelong <jerome.lelong@gmail.com>                    */
/*                                                                      */
/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as       */
/* published by the Free Software Foundation, either version 3 of the   */
/* License, or (at your option) any later version.                      */
/*                                                                      */
/* This program is distributed in the hope that it will be useful, but  */
/* WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    */
/* Lesser General Public License for more details.                      */
/*                                                                      */
/* You should have received a copy of the GNU Lesser General Public     */
/* License along with this program.  If not, see                        */
/* <http://www.gnu.org/licenses/>.                                      */
/************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "config.h"
#include "pnl_tridiag_matrix.h"
#include "pnl_random.h"
#include "tests.h"

static PnlTriDiagMat* create_random_tridiag (n, gen)
{
  PnlVect *dl, *du, *d;
  PnlTriDiagMat *M;
  d = pnl_vect_create (n);
  du = pnl_vect_create (n);
  dl = pnl_vect_create (n);
  pnl_vect_rand_uni (d, n, 0., 1., gen);
  pnl_vect_rand_uni (du, n-1, 0., 1., gen);
  pnl_vect_rand_uni (dl, n-1, 0., 1., gen);
  M = pnl_tridiagmat_create_from_ptr (n, dl->array, d->array, du->array);
  pnl_vect_free (&d);
  pnl_vect_free (&dl);
  pnl_vect_free (&du);
  return M;
}

static void triadiag_add_test ()
{
  PnlTriDiagMat *A, *B;
  PnlMat *fullA, *fullB;
  int n = 10;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;

  printf ("mat_plus tridiag test\n");

  pnl_rand_init (gen, 1, 1);
  A = create_random_tridiag (n, gen);
  B = create_random_tridiag (n, gen);
  fullA = pnl_tridiagmat_to_matrix (A);
  fullB = pnl_tridiagmat_to_matrix (B);

  printf ("A = ");
  pnl_mat_print_nsp (fullA); printf ("\n");
  printf ("B = ");
  pnl_mat_print_nsp (fullB); printf ("\n");

  pnl_tridiagmat_plus_tridiagmat (A, B);
  pnl_mat_free (&fullA);
  fullA = pnl_tridiagmat_to_matrix (A);
  printf ("C = ");
  pnl_mat_print_nsp (fullA); printf ("\n");

  pnl_tridiagmat_free (&A);
  pnl_tridiagmat_free (&B);
  pnl_mat_free (&fullA);
  pnl_mat_free (&fullB);
}

static void triadiag_mv_test ()
{
  PnlVect *y, *x;
  PnlTriDiagMat *M;
  PnlMat *full;
  int n = 10;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;

  printf ("mat_mult_vect tridiag test\n");

  pnl_rand_init (gen, 1, 1);
  x = pnl_vect_create (n);
  y = pnl_vect_create (n);
  pnl_vect_rand_normal (x, n, gen);
  M = create_random_tridiag (n, gen);
  full = pnl_tridiagmat_to_matrix (M);

  printf ("M = ");
  pnl_mat_print_nsp (full); printf ("\n");

  printf ("x = "); pnl_vect_print_nsp (x); printf ("\n");

  pnl_tridiagmat_mult_vect_inplace (y, M, x);
  printf ("y = "); pnl_vect_print_nsp (y); printf ("\n");

  pnl_vect_free (&x);
  pnl_vect_free (&y);
  pnl_tridiagmat_free (&M);
  pnl_mat_free (&full);
}

static void triadiag_lAxpby_test ()
{
  PnlVect *y, *x;
  PnlTriDiagMat *M;
  PnlMat *full;
  int n = 10;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;

  printf ("lAxpby tridiag test\n");

  pnl_rand_init (gen, 1, 1);
  x = pnl_vect_create (n);
  y = pnl_vect_create (n);
  pnl_vect_rand_normal (x, n, gen);
  pnl_vect_rand_normal (y, n, gen);
  M = create_random_tridiag (n, gen);
  full = pnl_tridiagmat_to_matrix (M);

  printf ("M = ");
  pnl_mat_print_nsp (full); printf ("\n");

  printf ("y = "); pnl_vect_print_nsp (y); printf ("\n");
  printf ("x = "); pnl_vect_print_nsp (x); printf ("\n");

  pnl_tridiagmat_lAxpby (1.5, M, x, 3., y);
  printf ("y = "); pnl_vect_print_nsp (y); printf ("\n");

  pnl_vect_free (&x);
  pnl_vect_free (&y);
  pnl_tridiagmat_free (&M);
  pnl_mat_free (&full);
}

static void triadiag_syslin_test ()
{
  PnlVect *b, *x;
  PnlTriDiagMat *M;
  PnlMat *full;
  int n = 10;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;

  printf ("Syslin tridiag test\n");

  pnl_rand_init (gen, 1, 1);
  x = pnl_vect_create (n);
  b = pnl_vect_create (n);
  pnl_vect_rand_normal (b, n, gen);
  M = create_random_tridiag (n, gen);
  full = pnl_tridiagmat_to_matrix (M);

  printf ("M = ");
  pnl_mat_print_nsp (full); printf ("\n");

  printf ("b = "); pnl_vect_print_nsp (b); printf ("\n");
  pnl_tridiagmat_lu_syslin (x, M, b);
  printf ("x = "); pnl_vect_print_nsp (x); printf ("\n");

  pnl_vect_free (&x);
  pnl_vect_free (&b);
  pnl_tridiagmat_free (&M);
  pnl_mat_free (&full);
}


void tridiag_matrix_test()
{
  triadiag_add_test ();
  triadiag_mv_test ();
  triadiag_lAxpby_test ();
  triadiag_syslin_test ();
}
