
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
#include "tests.h"

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
  pnl_mat_set_double (S, 0.);
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
  int m, n, nl, nu, gen;
  gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init (gen, 1, 1);
  m = n = 5;
  nl = 2; nu = 1;

  printf ("BandMat basic operations.\n");

  BMclone = pnl_band_mat_create (m, n, nl, nu);
  BM = create_band_mat (m, n, nl, nu, gen);
  printf ("BM = "); pnl_band_mat_print_as_full (BM); printf ("\n");
  pnl_band_mat_clone (BMclone, BM);

  printf ("M(%i, %i) = %f\n", 2, 2, pnl_band_mat_get (BM, 2, 2));
  printf ("M(%i, %i) = %f\n", 2, 3, pnl_band_mat_get (BM, 2, 3));
  printf ("M(%i, %i) = %f\n", 3, 1, pnl_band_mat_get (BM, 3, 1));
  printf("\n");

  pnl_band_mat_clone (BM, BMclone);
  pnl_band_mat_mult_double (BM, 2.0);
  printf ("2.0 * BM = "); pnl_band_mat_print_as_full (BM); printf ("\n");

  pnl_band_mat_clone (BM, BMclone);
  pnl_band_mat_plus_double (BM, 2.0);
  printf ("2.0 + BM = "); pnl_band_mat_print_as_full (BM); printf ("\n");

  pnl_band_mat_free (&BM);
  pnl_band_mat_free (&BMclone);
}

static void band_mat_ops_test ()
{
  PnlBandMat *BA, *BB, *BAclone;
  int m, n, nl, nu, gen;
  gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init (gen, 1, 1);
  m = n = 6;
  nl = 2; nu = 4;

  printf ("BandMat term by term operations.\n");

  BAclone = pnl_band_mat_create (m, n, nl, nu);
  BA = create_band_mat (m, n, nl, nu, gen);
  BB = create_band_mat (m, n, nl, nu, gen);
  printf ("BA = "); pnl_band_mat_print_as_full (BA); printf ("\n");
  printf ("BB = "); pnl_band_mat_print_as_full (BB); printf ("\n");
  pnl_band_mat_clone (BAclone, BA);

  pnl_band_mat_clone (BA, BAclone);
  pnl_band_mat_plus_band_mat (BA, BB);
  printf ("BA + BB  = "); pnl_band_mat_print_as_full (BA); printf ("\n");

  pnl_band_mat_clone (BA, BAclone);
  pnl_band_mat_mult_band_mat_term (BA, BB);
  printf ("BA .* BB  = "); pnl_band_mat_print_as_full (BA); printf ("\n");

  pnl_band_mat_free (&BA);
  pnl_band_mat_free (&BB);
  pnl_band_mat_free (&BAclone);
}

static void band_mat_syslin_test ()
{
  PnlBandMat *S, *Scopy, *B, *Bcopy;
  PnlVect *b, *x;
  PnlVectInt *p;
  int n = 5, nu = 2, nl = 3;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init (gen, n, n);
  b = pnl_vect_create (n);
  S = create_sys_pos_band_mat (n, nu, gen);
  pnl_vect_rand_normal (b, n, gen);

  Scopy = pnl_band_mat_copy (S);
  printf ("S = "); pnl_band_mat_print_as_full (S);
  printf ("\nb = "); pnl_vect_print_nsp (b);
  printf ("\n");
  x = pnl_vect_create(0);

  printf("test of pnl_band_mat_syslin (symmetric band matrix) : \n");
  pnl_band_mat_clone (S, Scopy);
  pnl_band_mat_syslin (x, S, b);
  pnl_vect_print(x);
  
  printf("test of pnl_band_mat_lu_syslin (symmetric band matrix) : \n");
  pnl_band_mat_clone (S, Scopy);
  p = pnl_vect_int_create (n);
  pnl_band_mat_lu (S, p);
  pnl_band_mat_lu_syslin (x, S, p, b);
  pnl_vect_print(x);
  pnl_vect_int_free (&p);


  B = create_invertible_band_mat (n, nl, nu, gen);
  Bcopy = pnl_band_mat_copy (B);

  printf("test of pnl_mat_syslin : \n");
  printf ("A = "); pnl_band_mat_print_as_full (B);
  pnl_band_mat_clone (B, Bcopy);
  pnl_band_mat_syslin (x, B, b);
  printf("x = "); pnl_vect_print_nsp(x); printf("\n");

  printf("test of pnl_band_mat_lu_syslin : \n");
  pnl_band_mat_clone (B, Bcopy);
  p = pnl_vect_int_create (n);
  pnl_band_mat_lu (B, p);
  pnl_band_mat_lu_syslin (x, B, p, b);
  printf("x = "); pnl_vect_print_nsp(x); printf("\n");

  pnl_vect_int_free (&p);
  pnl_band_mat_free (&B);
  pnl_band_mat_free (&Bcopy);
  pnl_band_mat_free (&S);
  pnl_band_mat_free (&Scopy);
  pnl_vect_free (&b);
  pnl_vect_free (&x);
}

int main ()
{
  basic_band_mat_test ();
  band_mat_ops_test ();
  band_mat_syslin_test ();
  return OK;
}
