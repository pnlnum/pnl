
/************************************************************************/
/* Copyright J�r�me Lelong <jerome.lelong@gmail.com>                    */
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
#include "pnl_matrix.h"
#include "pnl_band_matrix.h"
#include "pnl_random.h"
#include "tests.h"

extern void create_sym_pos_matrix (PnlMat *S, int n, int gen);
extern void create_invertible_matrix (PnlMat *A, int n, int gen);


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
PnlBandMat* create_bandmat (int m, int n, int nl, int nu, int gen)
{
  PnlMat *M;
  PnlBandMat *BM;

  M = pnl_mat_create (m, n);
  pnl_mat_rand_normal (M, m, n, gen);
  BM = pnl_bandmat_create_from_mat (M, nl, nu);
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
PnlBandMat* create_invertible_bandmat (int n, int nl, int nu, int gen)
{
  PnlMat *M;
  PnlBandMat *BM;

  M = pnl_mat_create (n, n);
  create_invertible_matrix (M, n, gen);
  BM = pnl_bandmat_create_from_mat (M, nl, nu);
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
PnlBandMat* create_sys_pos_bandmat (int n, int nu, int gen)
{
  PnlMat *M;
  PnlBandMat *BM;

  M = pnl_mat_create (n, n);
  create_sym_pos_matrix (M, n, gen);
  BM = pnl_bandmat_create_from_mat (M, nu, nu);
  pnl_mat_free (&M);
  return (BM);   
}

static void basic_bandmat_test ()
{
  PnlBandMat *BM, *BMclone;
  int m, n, nl, nu, gen;
  gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init (gen, 1, 1);
  m = n = 5;
  nl = 2; nu = 1;

  printf ("BandMat basic operations.\n");

  BMclone = pnl_bandmat_create (m, n, nl, nu);
  BM = create_bandmat (m, n, nl, nu, gen);
  printf ("BM = "); pnl_bandmat_print_as_full (BM); printf ("\n");
  pnl_bandmat_clone (BMclone, BM);

  printf ("M(%i, %i) = %f\n", 2, 2, pnl_bandmat_get (BM, 2, 2));
  printf ("M(%i, %i) = %f\n", 2, 3, pnl_bandmat_get (BM, 2, 3));
  printf ("M(%i, %i) = %f\n", 3, 1, pnl_bandmat_get (BM, 3, 1));
  printf("\n");

  pnl_bandmat_clone (BM, BMclone);
  pnl_bandmat_mult_double (BM, 2.0);
  printf ("2.0 * BM = "); pnl_bandmat_print_as_full (BM); printf ("\n");

  pnl_bandmat_clone (BM, BMclone);
  pnl_bandmat_plus_double (BM, 2.0);
  printf ("2.0 + BM = "); pnl_bandmat_print_as_full (BM); printf ("\n");

  pnl_bandmat_free (&BM);
  pnl_bandmat_free (&BMclone);
}

static void bandmat_ops_test ()
{
  PnlBandMat *BA, *BB, *BAclone;
  int m, n, nl, nu, gen;
  gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init (gen, 1, 1);
  m = n = 6;
  nl = 2; nu = 4;

  printf ("BandMat term by term operations.\n");

  BAclone = pnl_bandmat_create (m, n, nl, nu);
  BA = create_bandmat (m, n, nl, nu, gen);
  BB = create_bandmat (m, n, nl, nu, gen);
  printf ("BA = "); pnl_bandmat_print_as_full (BA); printf ("\n");
  printf ("BB = "); pnl_bandmat_print_as_full (BB); printf ("\n");
  pnl_bandmat_clone (BAclone, BA);

  pnl_bandmat_clone (BA, BAclone);
  pnl_bandmat_plus_bandmat (BA, BB);
  printf ("BA + BB  = "); pnl_bandmat_print_as_full (BA); printf ("\n");

  pnl_bandmat_clone (BA, BAclone);
  pnl_bandmat_mult_bandmat_term (BA, BB);
  printf ("BA .* BB  = "); pnl_bandmat_print_as_full (BA); printf ("\n");

  pnl_bandmat_free (&BA);
  pnl_bandmat_free (&BB);
  pnl_bandmat_free (&BAclone);
}

static void bandmat_syslin_test ()
{
  PnlBandMat *S, *Scopy, *B, *Bcopy;
  PnlVect *b, *x;
  PnlVectInt *p;
  int n = 5, nu = 2, nl = 3;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init (gen, n, n);
  b = pnl_vect_create (n);
  S = create_sys_pos_bandmat (n, nu, gen);
  pnl_vect_rand_normal (b, n, gen);

  Scopy = pnl_bandmat_copy (S);
  printf ("S = "); pnl_bandmat_print_as_full (S);
  printf ("\nb = "); pnl_vect_print_nsp (b);
  printf ("\n");
  x = pnl_vect_create(0);

  printf("test of pnl_bandmat_syslin (symmetric band matrix) : \n");
  pnl_bandmat_clone (S, Scopy);
  pnl_bandmat_syslin (x, S, b);
  pnl_vect_print(x);
  
  printf("test of pnl_bandmat_lu_syslin (symmetric band matrix) : \n");
  pnl_bandmat_clone (S, Scopy);
  p = pnl_vect_int_create (n);
  pnl_bandmat_lu (S, p);
  pnl_bandmat_lu_syslin (x, S, p, b);
  pnl_vect_print(x);
  pnl_vect_int_free (&p);


  B = create_invertible_bandmat (n, nl, nu, gen);
  Bcopy = pnl_bandmat_copy (B);

  printf("test of pnl_mat_syslin : \n");
  printf ("A = "); pnl_bandmat_print_as_full (B);
  pnl_bandmat_clone (B, Bcopy);
  pnl_bandmat_syslin (x, B, b);
  printf("x = "); pnl_vect_print_nsp(x); printf("\n");

  printf("test of pnl_bandmat_lu_syslin : \n");
  pnl_bandmat_clone (B, Bcopy);
  p = pnl_vect_int_create (n);
  pnl_bandmat_lu (B, p);
  pnl_bandmat_lu_syslin (x, B, p, b);
  printf("x = "); pnl_vect_print_nsp(x); printf("\n");

  pnl_vect_int_free (&p);
  pnl_bandmat_free (&B);
  pnl_bandmat_free (&Bcopy);
  pnl_bandmat_free (&S);
  pnl_bandmat_free (&Scopy);
  pnl_vect_free (&b);
  pnl_vect_free (&x);
}

void band_matrix_test()
{
  basic_bandmat_test ();
  bandmat_ops_test ();
  bandmat_syslin_test ();
}
