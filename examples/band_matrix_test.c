
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
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_band_matrix.h"
#include "pnl/pnl_random.h"
#include "tests_utils.h"

/** 
 * Creates an invertible matrix of size nxn using the exponential function
 * 
 * @param A Output parameter, contains a n x n intervible matrix.
 * @param n size of the matrix
 * @param gen index of the generator
 */
static void create_invertible_matrix (PnlMat *A, int n, int gen)
{
  PnlMat *B;

  B = pnl_mat_create (n, n);
  pnl_mat_rand_uni2 (B, n, n, 0., 1., gen);

  pnl_mat_exp (A, B);
  pnl_mat_free (&B);
}

/** 
 * Creates a symmetric positive definite matrix of size nxn
 * 
 * @param S Output parameter, contains a n x n symmetric positive definite matrix.
 * @param n size of the matrix
 * @param gen index of the generator
 */
static void create_sym_pos_matrix (PnlMat *S, int n, int gen)
{
  PnlVect *b;
  double g;
  int i;

  b = pnl_vect_create (n);
  pnl_vect_rand_uni (b, n, 0., 1., gen);

  /* S is set to a diagonal matrix with positive eigenvalues */
  pnl_mat_set_all (S, 0.);
  for ( i=0 ; i<S->n ; i++ )
    {
      do
        {
          g = fabs (pnl_rand_normal (gen));
        }
      while (g < 1E-4);
      pnl_mat_set (S, i, i, g);
    }

  pnl_mat_dger (1., b, b, S);
  pnl_vect_free (&b);
}


/** 
 * Creates a random band matrix
 * 
 * @param m nb rows
 * @param n nb cols
 * @param nl nb of lower diagonals
 * @param nu nb of upper diagonals
 * @param gen index of the generator
 * 
 * @return a band matrix
 */
PnlBandMat* create_band_mat (int m, int n, int nl, int nu, int gen)
{
  PnlMat *M;
  PnlBandMat *BM;

  M = pnl_mat_create (m, n);
  pnl_mat_rand_normal (M, m, n, gen);
  BM = pnl_band_mat_create_from_mat (M, nl, nu);
  pnl_mat_free (&M);
  return (BM);   
}

/** 
 * Creates a random invertible band matrix
 * 
 * @param n size of the matrix
 * @param nl nb of lower diagonals
 * @param nu nb of upper diagonals
 * @param gen index of the generator
 * 
 * @return a band matrix
 */
PnlBandMat* create_invertible_band_mat (int n, int nl, int nu, int gen)
{
  PnlMat *M;
  PnlBandMat *BM;

  M = pnl_mat_create (n, n);
  create_invertible_matrix (M, n, gen);
  BM = pnl_band_mat_create_from_mat (M, nl, nu);
  pnl_mat_free (&M);
  return (BM);   
}
/** 
 * Creates a symmetric positive definite band matrix
 * 
 * @param n size of the matrix
 * @param nu nb of lower or upper diagonals
 * @param gen index of the generator
 * 
 * @return a band matrix
 */
PnlBandMat* create_sys_pos_band_mat (int n, int nu, int gen)
{
  PnlMat *M;
  PnlBandMat *BM;

  M = pnl_mat_create (n, n);
  create_sym_pos_matrix (M, n, gen);
  BM = pnl_band_mat_create_from_mat (M, nu, nu);
  pnl_mat_free (&M);
  return (BM);   
}

static void basic_band_mat_test ()
{
  PnlBandMat *BM, *BMclone;
  int m, n, nl, nu, gen, i, j;
  double abserr = 1E-12;
  gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init (gen, 1, 1);
  m = n = 5;
  nl = 2; nu = 1;

  BMclone = pnl_band_mat_create (m, n, nl, nu);
  BM = create_band_mat (m, n, nl, nu, gen);
  pnl_band_mat_clone (BMclone, BM);


  pnl_band_mat_clone (BM, BMclone);
  pnl_band_mat_mult_scalar (BM, 2.0);
  for ( i=0 ; i<BM->m ; i++ )
    {
      for ( j=MAX(0, i-BM->nl) ; j<MIN(BM->n, i+BM->nu) ; j++ )
        {
          if (pnl_cmp_eq_abs(pnl_band_mat_get(BM, i, j), 
                             2. * pnl_band_mat_get(BMclone, i, j), abserr))
            {
              pnl_test_set_fail ("band_mat_mult_double", pnl_band_mat_get(BM, i, j),
                                 2. * pnl_band_mat_get(BMclone, i, j));
              goto J1;
            }
        }
    }
  pnl_test_set_ok ("band_mat_mult_double");

J1:
  pnl_band_mat_clone (BM, BMclone);
  pnl_band_mat_plus_scalar (BM, 2.0);
  for ( i=0 ; i<BM->m ; i++ )
    {
      for ( j=MAX(0, i-BM->nl) ; j<MIN(BM->n, i+BM->nu) ; j++ )
        {
          if (pnl_cmp_eq_abs(pnl_band_mat_get(BM, i, j), 
                             2. + pnl_band_mat_get(BMclone, i, j), abserr))
            {
              pnl_test_set_fail ("band_mat_plus_double", pnl_band_mat_get(BM, i, j),
                                 2. + pnl_band_mat_get(BMclone, i, j));
              goto J2;
            }
        }
    }
  pnl_test_set_ok ("band_mat_plus_double");

J2:
  pnl_band_mat_free (&BM);
  pnl_band_mat_free (&BMclone);
}

static void band_mat_ops_test ()
{
  PnlBandMat *BA, *BB, *BAclone;
  int m, n, nl, nu, gen, i, j;
  double abserr = 1E-12;
  gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init (gen, 1, 1);
  m = n = 6;
  nl = 2; nu = 4;


  BAclone = pnl_band_mat_create (m, n, nl, nu);
  BA = create_band_mat (m, n, nl, nu, gen);
  BB = create_band_mat (m, n, nl, nu, gen);
  pnl_band_mat_clone (BAclone, BA);

  pnl_band_mat_clone (BA, BAclone);
  pnl_band_mat_plus_band_mat (BA, BB);
  for ( i=0 ; i<BA->m ; i++ )
    {
      for ( j=MAX(0, i-BA->nl) ; j<MIN(BA->n, i+BA->nu) ; j++ )
      for ( j=0 ; j<BA->n ; j++ )
        {
          double Aij, Bij, Rij;
          Aij = pnl_band_mat_get (BAclone, i, j);
          Bij = pnl_band_mat_get (BB, i, j);
          Rij = pnl_band_mat_get (BA, i, j);
          if (pnl_cmp_eq_abs(Aij + Bij, Rij, abserr))
            {
              pnl_test_set_fail ("band_mat_plus_band_mat", Rij, Aij + Bij);
              goto J1;
            }
        }
    }
  pnl_test_set_ok ("band_mat_plus_band_mat");
J1:
  pnl_band_mat_clone (BA, BAclone);
  pnl_band_mat_mult_band_mat_term (BA, BB);
  for ( i=0 ; i<BA->m ; i++ )
    {
      for ( j=MAX(0, i-BA->nl) ; j<MIN(BA->n, i+BA->nu) ; j++ )
        {
          double Aij, Bij, Rij;
          Aij = pnl_band_mat_get (BAclone, i, j);
          Bij = pnl_band_mat_get (BB, i, j);
          Rij = pnl_band_mat_get (BA, i, j);
          if (pnl_cmp_eq_abs(Aij * Bij, Rij, abserr))
            {
              pnl_test_set_fail ("band_mat_mult_band_mat", Rij, Aij * Bij);
              goto J2;
            }
        }
    }
  pnl_test_set_ok ("band_mat_mult_band_mat");
J2:
  pnl_band_mat_free (&BA);
  pnl_band_mat_free (&BB);
  pnl_band_mat_free (&BAclone);
}

static void band_mat_syslin_test ()
{
  PnlBandMat *S, *Scopy, *B, *Bcopy;
  PnlVect *b, *x, *Sx;
  PnlVectInt *p;
  double abserr = 1E-12;
  int n = 5, nu = 2, nl = 3;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init (gen, n, n);
  b = pnl_vect_create (n);
  Sx = pnl_vect_create (n);
  S = create_sys_pos_band_mat (n, nu, gen);
  pnl_vect_rand_normal (b, n, gen);

  Scopy = pnl_band_mat_copy (S);
  x = pnl_vect_create(0);

  pnl_band_mat_clone (S, Scopy);
  pnl_band_mat_syslin (x, S, b);

  pnl_band_mat_mult_vect_inplace (Sx, Scopy, x); 
  pnl_test_vect_eq_abs (Sx, b, abserr, "band_mat_syslin (symmetric)", "");

  pnl_band_mat_clone (S, Scopy);
  p = pnl_vect_int_create (n);
  pnl_band_mat_lu (S, p);
  pnl_band_mat_lu_syslin (x, S, p, b);
  pnl_band_mat_mult_vect_inplace (Sx, Scopy, x); 
  pnl_test_vect_eq_abs (Sx, b, abserr, "band_mat_lu_syslin (symmetric)", "");
  pnl_vect_int_free (&p);

  B = create_invertible_band_mat (n, nl, nu, gen);
  Bcopy = pnl_band_mat_copy (B);

  pnl_band_mat_clone (B, Bcopy);
  pnl_band_mat_syslin (x, B, b);
  pnl_band_mat_mult_vect_inplace (Sx, Bcopy, x); 
  pnl_test_vect_eq_abs (Sx, b, abserr, "band_mat_syslin", "");


  pnl_band_mat_clone (B, Bcopy);
  p = pnl_vect_int_create (n);
  pnl_band_mat_lu (B, p);
  pnl_band_mat_lu_syslin (x, B, p, b);
  pnl_band_mat_mult_vect_inplace (Sx, Bcopy, x); 
  pnl_test_vect_eq_abs (Sx, b, abserr, "band_mat_lu_syslin", "");

  pnl_vect_int_free (&p);
  pnl_band_mat_free (&B);
  pnl_band_mat_free (&Bcopy);
  pnl_band_mat_free (&S);
  pnl_band_mat_free (&Scopy);
  pnl_vect_free (&b);
  pnl_vect_free (&Sx);
  pnl_vect_free (&x);
}

int main (int argc, char *argv[])
{
  pnl_test_init (argc, argv);
  basic_band_mat_test ();
  band_mat_ops_test ();
  band_mat_syslin_test ();
  exit(pnl_test_finalize ("Band matrices"));
}
