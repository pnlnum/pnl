
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

#include "pnl/pnl_tridiag_matrix.h"
#include "pnl/pnl_random.h"
#include "tests.h"

static void pnl_tridiag_mat_print_as_full (const PnlTridiagMat *T)
{
  PnlMat *Tfull = pnl_tridiag_mat_to_mat (T);
  pnl_mat_print_nsp (Tfull);
  pnl_mat_free (&Tfull);
}

static PnlTridiagMat* create_random_tridiag (n, gen)
{
  PnlVect *dl, *du, *d;
  PnlTridiagMat *M;
  d = pnl_vect_create (n);
  du = pnl_vect_create (n);
  dl = pnl_vect_create (n);
  pnl_vect_rand_uni (d, n, 0., 1., gen);
  pnl_vect_rand_uni (du, n-1, 0., 1., gen);
  pnl_vect_rand_uni (dl, n-1, 0., 1., gen);
  M = pnl_tridiag_mat_create_from_ptr (n, dl->array, d->array, du->array);
  pnl_vect_free (&d);
  pnl_vect_free (&dl);
  pnl_vect_free (&du);
  return M;
}

static void tridiag_create_test ()
{
  int n;
  PnlTridiagMat *A;
  printf ("Create tridiag test\n");
  n = 5;
  A = pnl_tridiag_mat_create_from_two_double (n, -2.5, 4);
  pnl_tridiag_mat_print_as_full (A); printf("\n");
  pnl_tridiag_mat_free (&A);
}

static void tridiag_get_test ()
{
  PnlTridiagMat *T;
  int n, gen;
  printf ("Get tridiag test\n");
  n = 5;
  gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init (gen, n, 1);
  T = create_random_tridiag (n, gen);
  pnl_tridiag_mat_print_as_full (T);
  printf ("T(%i,%i+1) = %f\n", 3, 3, pnl_tridiag_mat_get (T, 3, 1));
  printf ("T(%i,%i) = %f\n", 3, 3, pnl_tridiag_mat_get (T, 3, 0));
  printf ("T(%i,%i-1) = %f\n", 3, 3, pnl_tridiag_mat_get (T, 3, -1));
  printf ("\n");
  pnl_tridiag_mat_free (&T);
}

static void tridiag_add_test ()
{
  PnlTridiagMat *A, *Aclone, *B;
  int n = 5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;

  printf ("mat_plus tridiag test\n");

  pnl_rand_init (gen, 1, 1);
  Aclone = pnl_tridiag_mat_create (0);
  A = create_random_tridiag (n, gen);
  B = create_random_tridiag (n, gen);
  pnl_tridiag_mat_clone (Aclone, A);

  printf ("A = ");
  pnl_tridiag_mat_print_as_full (A); printf ("\n");
  printf ("B = ");
  pnl_tridiag_mat_print_as_full (B); printf ("\n");

  pnl_tridiag_mat_clone (A, Aclone);
  pnl_tridiag_mat_plus_tridiag_mat (A, B);
  printf ("A+B = ");
  pnl_tridiag_mat_print_as_full (A); printf ("\n");

  pnl_tridiag_mat_clone (A, Aclone);
  pnl_tridiag_mat_minus_tridiag_mat (A, B);
  printf ("A-B = ");
  pnl_tridiag_mat_print_as_full (A); printf ("\n");

  pnl_tridiag_mat_clone (A, Aclone);
  pnl_tridiag_mat_mult_tridiag_mat_term (A, B);
  printf ("A.*B = ");
  pnl_tridiag_mat_print_as_full (A); printf ("\n");
  pnl_tridiag_mat_free (&A);
  pnl_tridiag_mat_free (&Aclone);
  pnl_tridiag_mat_free (&B);
}

static void tridiag_mv_test ()
{
  PnlVect *y, *x;
  PnlTridiagMat *M;
  int n = 5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;

  printf ("mat_mult_vect tridiag test\n");

  pnl_rand_init (gen, 1, 1);
  x = pnl_vect_create (n);
  y = pnl_vect_create (n);
  pnl_vect_rand_normal (x, n, gen);
  M = create_random_tridiag (n, gen);

  printf ("M = ");
  pnl_tridiag_mat_print_as_full (M); printf ("\n");

  printf ("x = "); pnl_vect_print_nsp (x); printf ("\n");

  pnl_tridiag_mat_mult_vect_inplace (y, M, x);
  printf ("y = "); pnl_vect_print_nsp (y); printf ("\n");

  pnl_vect_free (&x);
  pnl_vect_free (&y);
  pnl_tridiag_mat_free (&M);
}

static void tridiag_lAxpby_test ()
{
  PnlVect *y, *x;
  PnlTridiagMat *M;
  int n = 5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;

  printf ("lAxpby tridiag test\n");

  pnl_rand_init (gen, 1, 1);
  x = pnl_vect_create (n);
  y = pnl_vect_create (n);
  pnl_vect_rand_normal (x, n, gen);
  pnl_vect_rand_normal (y, n, gen);
  M = create_random_tridiag (n, gen);

  printf ("M = ");
  pnl_tridiag_mat_print_as_full (M); printf ("\n");

  printf ("y = "); pnl_vect_print_nsp (y); printf ("\n");
  printf ("x = "); pnl_vect_print_nsp (x); printf ("\n");

  pnl_tridiag_mat_lAxpby (1.5, M, x, 3., y);
  printf ("y = "); pnl_vect_print_nsp (y); printf ("\n");

  pnl_vect_free (&x);
  pnl_vect_free (&y);
  pnl_tridiag_mat_free (&M);
}

static void tridiag_syslin_test ()
{
  PnlVect *b, *x;
  PnlTridiagMat *M;
  int n = 5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;

  printf ("Syslin tridiag test\n");

  pnl_rand_init (gen, 1, 1);
  x = pnl_vect_create (n);
  b = pnl_vect_create (n);
  pnl_vect_rand_normal (b, n, gen);
  M = create_random_tridiag (n, gen);

  printf ("M = ");
  pnl_tridiag_mat_print_as_full (M); printf ("\n");

  printf ("b = "); pnl_vect_print_nsp (b); printf ("\n");
  pnl_tridiag_mat_syslin (x, M, b);
  printf ("x = "); pnl_vect_print_nsp (x); printf ("\n");

  pnl_vect_free (&x);
  pnl_vect_free (&b);
  pnl_tridiag_mat_free (&M);
}

static void triadiag_scalar_prod_test ()
{
  PnlVect *y, *x;
  PnlTridiagMat *M;
  PnlMat *full;
  int n = 5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;

  printf ("x' * A * y tridiag test\n");

  pnl_rand_init (gen, 1, 1);
  x = pnl_vect_create (n);
  y = pnl_vect_create (n);
  pnl_vect_rand_normal (x, n, gen);
  pnl_vect_rand_normal (y, n, gen);
  M = create_random_tridiag (n, gen);
  full = pnl_tridiag_mat_to_mat (M);

  printf ("M = ");
  pnl_mat_print_nsp (full); printf ("\n");

  printf ("x = "); pnl_vect_print_nsp (x); printf ("\n");
  printf ("y = "); pnl_vect_print_nsp (y); printf ("\n");
  printf ("x' * M * y = %f\n", pnl_tridiag_mat_scalar_prod (x, M, y));
  pnl_vect_free (&x);
  pnl_vect_free (&y);
  pnl_tridiag_mat_free (&M);
  pnl_mat_free (&full); 
}

static void tridiag_mat_op_tridiag_mat_test ()
{
  PnlTridiagMat *T1, *T2;
  int gen, n;

  gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  n= 5;
  pnl_rand_init (gen, n, 1);

  T1 = create_random_tridiag (n, gen);
  T2 = create_random_tridiag (n, gen);

  printf ("T1 = "); pnl_tridiag_mat_print_as_full (T1); printf ("\n");
  printf ("T2 = "); pnl_tridiag_mat_print_as_full (T2); printf ("\n");
  pnl_tridiag_mat_plus_tridiag_mat (T1, T2);
  printf ("T1 + T2 = "); pnl_tridiag_mat_print_as_full (T1); printf ("\n");

  pnl_tridiag_mat_free (&T1);
  pnl_tridiag_mat_free (&T2);
}



int main ()
{
  tridiag_create_test ();
  tridiag_get_test ();
  tridiag_add_test ();
  tridiag_mv_test ();
  tridiag_lAxpby_test ();
  tridiag_syslin_test ();
  triadiag_scalar_prod_test ();
  tridiag_mat_op_tridiag_mat_test ();
  return OK;
}
