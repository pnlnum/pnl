
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
#include "tests_utils.h"

/* static void pnl_tridiag_mat_print_as_full (const PnlTridiagMat *T) */
/* { */
/*   PnlMat *Tfull = pnl_tridiag_mat_to_mat (T); */
/*   pnl_mat_print_nsp (Tfull); */
/*   pnl_mat_free (&Tfull); */
/* } */

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

static void tridiag_get_test ()
{
  PnlTridiagMat *T;
  int n, gen;
  double abserr;
  abserr = 1E-18;
  n = 5;
  gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init (gen, n, 1);
  T = create_random_tridiag (n, gen);
  pnl_test_eq_abs (pnl_tridiag_mat_get (T, 3, 1), T->DU[3], abserr, "tridiag_mat_get (+1)", "");
  pnl_test_eq_abs (pnl_tridiag_mat_get (T, 3, 0), T->D[3], abserr, "tridiag_mat_get (0)", "");
  pnl_test_eq_abs (pnl_tridiag_mat_get (T, 3, -1), T->DL[3-1], abserr, "tridiag_mat_get (-1)", "");
  pnl_tridiag_mat_free (&T);
}

static void tridiag_add_test ()
{
  PnlTridiagMat *A, *Aclone, *B;
  int i, n = 5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  double abserr = 1E-12;


  pnl_rand_init (gen, 1, 1);
  Aclone = pnl_tridiag_mat_create (0);
  A = create_random_tridiag (n, gen);
  B = create_random_tridiag (n, gen);
  pnl_tridiag_mat_clone (Aclone, A);

  pnl_tridiag_mat_clone (A, Aclone);
  pnl_tridiag_mat_plus_tridiag_mat (A, B);
  for ( i=0 ; i<A->size-1 ; i++ )
    {
      if ( pnl_cmp_eq_abs(A->DL[i], Aclone->DL[i] + B->DL[i], abserr) ) { pnl_test_set_fail ("tridiag_mat_plus_tridiag_mat (-1)", A->DL[i], Aclone->DL[i] + B->DL[i] ); goto J1;}
      if ( pnl_cmp_eq_abs(A->DU[i], Aclone->DU[i] + B->DU[i], abserr) ) { pnl_test_set_fail ("tridiag_mat_plus_tridiag_mat (+1)", A->DU[i], Aclone->DU[i] + B->DU[i] ); goto J1;}
      if ( pnl_cmp_eq_abs(A->D[i], Aclone->D[i] + B->D[i], abserr) ) { pnl_test_set_fail ("tridiag_mat_plus_tridiag_mat (0)", A->D[i], Aclone->D[i] + B->D[i] ); goto J1; }
    }
  if ( pnl_cmp_eq_abs(A->D[A->size-1], Aclone->D[A->size-1] + B->D[A->size-1], abserr) ) { pnl_test_set_fail ("tridiag_mat_plus_tridiag_mat (0)", A->D[A->size-1], Aclone->D[A->size-1] + B->D[A->size-1] ); goto J1;}

  pnl_test_set_ok ("tridiag_mat_plus_tridiag_mat");
J1:

  pnl_tridiag_mat_clone (A, Aclone);
  pnl_tridiag_mat_minus_tridiag_mat (A, B);
  for ( i=0 ; i<A->size-1 ; i++ )
    {
      if ( pnl_cmp_eq_abs(A->DL[i], Aclone->DL[i] - B->DL[i], abserr) ) { pnl_test_set_fail ("tridiag_mat_minus_tridiag_mat (-1)", A->DL[i], Aclone->DL[i] - B->DL[i] ); goto J2;}
      if ( pnl_cmp_eq_abs(A->DU[i], Aclone->DU[i] - B->DU[i], abserr) ) { pnl_test_set_fail ("tridiag_mat_minus_tridiag_mat (+1)", A->DU[i], Aclone->DU[i] - B->DU[i] ); goto J2;}
      if ( pnl_cmp_eq_abs(A->D[i], Aclone->D[i] - B->D[i], abserr) ) { pnl_test_set_fail ("tridiag_mat_minus_tridiag_mat (0)", A->D[i], Aclone->D[i] - B->D[i] ); goto J2; }
    }
  if ( pnl_cmp_eq_abs(A->D[A->size-1], Aclone->D[A->size-1] - B->D[A->size-1], abserr) ) { pnl_test_set_fail ("tridiag_mat_minus_tridiag_mat (0)", A->D[A->size-1], Aclone->D[A->size-1] - B->D[A->size-1] ); goto J2;}

  pnl_test_set_ok ("tridiag_mat_minus_tridiag_mat");
J2:

  pnl_tridiag_mat_clone (A, Aclone);
  pnl_tridiag_mat_mult_tridiag_mat_term (A, B);
  for ( i=0 ; i<A->size-1 ; i++ )
    {
      if ( pnl_cmp_eq_abs(A->DL[i], Aclone->DL[i] * B->DL[i], abserr) ) { pnl_test_set_fail ("tridiag_mat_mult_tridiag_mat (-1)", A->DL[i], Aclone->DL[i] * B->DL[i] ); goto J3;}
      if ( pnl_cmp_eq_abs(A->DU[i], Aclone->DU[i] * B->DU[i], abserr) ) { pnl_test_set_fail ("tridiag_mat_mult_tridiag_mat (+1)", A->DU[i], Aclone->DU[i] * B->DU[i] ); goto J3;}
      if ( pnl_cmp_eq_abs(A->D[i], Aclone->D[i] * B->D[i], abserr) ) { pnl_test_set_fail ("tridiag_mat_mult_tridiag_mat (0)", A->D[i], Aclone->D[i] * B->D[i] ); goto J3; }
    }
  if ( pnl_cmp_eq_abs(A->D[A->size-1], Aclone->D[A->size-1] * B->D[A->size-1], abserr) ) { pnl_test_set_fail ("tridiag_mat_mult_tridiag_mat (0)", A->D[A->size-1], Aclone->D[A->size-1] * B->D[A->size-1] ); goto J3;}

  pnl_test_set_ok ("tridiag_mat_mult_tridiag_mat_term");
J3:
  pnl_tridiag_mat_free (&A);
  pnl_tridiag_mat_free (&Aclone);
  pnl_tridiag_mat_free (&B);
}

static void tridiag_mv_test ()
{
  PnlVect *y, *x, *Fy;
  PnlMat *FM;
  PnlTridiagMat *M;
  int n = 5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;

  pnl_rand_init (gen, 1, 1);
  x = pnl_vect_create (n);
  y = pnl_vect_create (n);
  pnl_vect_rand_normal (x, n, gen);
  M = create_random_tridiag (n, gen);
  FM = pnl_tridiag_mat_to_mat (M);

  pnl_tridiag_mat_mult_vect_inplace (y, M, x);
  Fy = pnl_mat_mult_vect (FM, x);
  pnl_test_vect_eq_abs (y, Fy, 1E-18, "tridiag_mat_mutl_vect", "");

  pnl_vect_free (&x);
  pnl_vect_free (&y);
  pnl_vect_free (&Fy);
  pnl_mat_free (&FM);
  pnl_tridiag_mat_free (&M);
}

static void tridiag_lAxpby_test ()
{
  PnlVect *y, *x, *Fy;
  PnlTridiagMat *M;
  PnlMat *FM;
  int n = 5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;

  pnl_rand_init (gen, 1, 1);
  x = pnl_vect_create (n);
  y = pnl_vect_create (n);
  pnl_vect_rand_normal (x, n, gen);
  pnl_vect_rand_normal (y, n, gen);
  Fy = pnl_vect_copy (y);
  M = create_random_tridiag (n, gen);
  FM = pnl_tridiag_mat_to_mat (M);

  pnl_tridiag_mat_lAxpby (1.5, M, x, 3., y);
  pnl_mat_lAxpby (1.5, FM, x, 3., Fy);
  pnl_test_vect_eq_abs (y, Fy, 1E-12, "tridiag_mat_lAxpby", "");

  pnl_vect_free (&x);
  pnl_vect_free (&y);
  pnl_vect_free (&Fy);
  pnl_mat_free (&FM);
  pnl_tridiag_mat_free (&M);
}

static void tridiag_syslin_test ()
{
  PnlVect *b, *x, *Mx;
  PnlTridiagMat *M, *Mcopy;
  int n = 5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;

  pnl_rand_init (gen, 1, 1);
  x = pnl_vect_create (n);
  b = pnl_vect_create (n);
  pnl_vect_rand_normal (b, n, gen);
  M = create_random_tridiag (n, gen);

  Mcopy = pnl_tridiag_mat_copy (M);
  pnl_tridiag_mat_syslin (x, M, b);
  Mx = pnl_tridiag_mat_mult_vect (Mcopy, x);
  pnl_test_vect_eq_abs (Mx, b, 1E-12, "tridiag_mat_syslin", "");

  pnl_vect_free (&x);
  pnl_vect_free (&Mx);
  pnl_vect_free (&b);
  pnl_tridiag_mat_free (&M);
  pnl_tridiag_mat_free (&Mcopy);
}

static void tridiag_lu_syslin_test ()
{
  PnlVect *b, *x, *Mx;
  PnlTridiagMat *M;
  PnlTridiagMatLU *LU;
  int n = 5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;

  pnl_rand_init (gen, 1, 1);
  x = pnl_vect_create (n);
  b = pnl_vect_create (n);
  pnl_vect_rand_normal (b, n, gen);
  M = create_random_tridiag (n, gen);

  LU = pnl_tridiag_mat_lu_new ();
  pnl_tridiag_mat_lu_compute (LU, M);
  pnl_tridiag_mat_lu_syslin (x, LU, b);
  Mx = pnl_tridiag_mat_mult_vect (M, x);
  pnl_test_vect_eq_abs (Mx, b, 1E-12, "tridiag_mat_syslin", "");

  pnl_vect_free (&x);
  pnl_vect_free (&Mx);
  pnl_vect_free (&b);
  pnl_tridiag_mat_free (&M);
  pnl_tridiag_mat_lu_free (&LU);
}


static void triadiag_scalar_prod_test ()
{
  PnlVect *y, *x;
  PnlTridiagMat *M;
  PnlMat *full;
  int n = 5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;


  pnl_rand_init (gen, 1, 1);
  x = pnl_vect_create (n);
  y = pnl_vect_create (n);
  pnl_vect_rand_normal (x, n, gen);
  pnl_vect_rand_normal (y, n, gen);
  M = create_random_tridiag (n, gen);
  full = pnl_tridiag_mat_to_mat (M);
  pnl_test_eq_abs ( pnl_tridiag_mat_scalar_prod (M, x, y), 
                    pnl_mat_scalar_prod (full, x, y),
                    1E-12, "tridiag_mat_scalar_prod", "");

  
  pnl_vect_free (&x);
  pnl_vect_free (&y);
  pnl_tridiag_mat_free (&M);
  pnl_mat_free (&full); 
}

int main (int argc, char *argv[])
{
  pnl_test_init (argc, argv);
  tridiag_get_test ();
  tridiag_add_test ();
  tridiag_mv_test ();
  tridiag_lAxpby_test ();
  tridiag_syslin_test ();
  tridiag_lu_syslin_test ();
  triadiag_scalar_prod_test ();
  exit (pnl_test_finalize("Tridiag matrices"));
}
