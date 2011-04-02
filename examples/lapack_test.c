
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

#include "pnl/pnl_matrix.h"
#include "pnl/pnl_matrix_complex.h"
#include "pnl/pnl_random.h"
#include "tests_utils.h"


/**
 * Creates a diagonalizable matrix of size n
 * @param n an integer
 * @return a matrix
 */
static PnlMat* pnl_mat_create_diagonalizable (int n)
{
  PnlMat *A;
  PnlMat *P, *invP;
  int i, j;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;

  pnl_rand_init (gen, n, n);
  A = pnl_mat_create (n, n);
  P = pnl_mat_create (n, n);
  invP = pnl_mat_create (n, n);
  pnl_mat_rand_uni2 (A, n, n, 0., 1., gen);

  /* P = exp (rand) is always invertible */
  pnl_mat_exp(P, A);
  pnl_mat_inverse(invP, P);

  /* P = P * D */
  for ( j=0 ; j<n ; j++ )
    {
      double dj = fabs(pnl_rand_uni (gen));
      for ( i=0 ; i<n ; i++ )
        {
          MLET (P, i, j) *= dj;
        }
    }
  pnl_mat_mult_mat_inplace (A, P, invP);
  pnl_mat_free (&P);
  pnl_mat_free (&invP);
  return A;
}

static void pnl_mat_eigen_test ()
{
  PnlMat *A, *B, *P, *Pinv, *V;
  PnlVect *v;
  int n = 4;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init (gen, n, n);
  A = pnl_mat_create (n, n);
  B = pnl_mat_create (n, n);
  P = pnl_mat_create (n, n);
  Pinv = pnl_mat_create (n, n);
  pnl_mat_rand_uni2 (A, n, n, 0., 1., gen);
  pnl_mat_dgemm ('N', 'T', 1., A, A, 0., B);
  v = pnl_vect_create (0);

  pnl_mat_eigen (v, P, B, TRUE);
  V = pnl_mat_create_diag (v);
  pnl_mat_inverse(Pinv, P);
  pnl_mat_mult_mat_inplace (A, P, V);  /* P V */
  pnl_mat_mult_mat_inplace (V, A, Pinv); /* P V P^-1 */
  pnl_test_mat_eq_abs(V, B, 1E-8, "eigen", "");
  

  pnl_mat_free (&V);
  pnl_mat_free(&A);

  A = pnl_mat_create_diagonalizable (n);
  pnl_mat_eigen (v, P, A, TRUE);
  V = pnl_mat_create_diag (v);
  pnl_mat_inverse(Pinv, P);
  pnl_mat_mult_mat_inplace (B, P, V);  /* P V */
  pnl_mat_mult_mat_inplace (V, B, Pinv); /* P V P^-1 */
  pnl_test_mat_eq_abs(V, A, 1E-8, "eigen", "");

  pnl_vect_free (&v);
  pnl_mat_free (&V);
  pnl_mat_free (&A);
  pnl_mat_free (&B);
  pnl_mat_free (&P);
  pnl_mat_free (&Pinv);
}

static void pnl_mat_log_test ()
{
  PnlMat *A, *B, *expB;
  int n = 5;
  A = pnl_mat_create_diagonalizable (n);
  B = pnl_mat_create (n,n);
  expB = pnl_mat_create (n,n);
  pnl_mat_log (B, A);
  pnl_mat_exp (expB, B);
  pnl_test_mat_eq_abs( expB, A, 1E-8, "logm / expm", "");
  pnl_mat_free (&B);
  pnl_mat_free (&expB);
  pnl_mat_free (&A);
}

static void pnl_mat_ls_test ()
{
  PnlMat *A, *B, *Bcopy, *AB;
  PnlVect *b, *Ab, *bcopy;
  int m, n, gen, nrhs;
  m = 5; n = 6; nrhs = 3;
  gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  A = pnl_mat_create (m, n);
  B = pnl_mat_create (m, nrhs);
  b = pnl_vect_create (m);
  pnl_rand_init (gen, m, n);
  pnl_mat_rand_uni2 (A, m, n, 0, 1, gen);
  pnl_mat_rand_uni2 (B, m, nrhs, 0, 1, gen);
  pnl_vect_rand_uni (b, m, 0, 1, gen);

  Bcopy = pnl_mat_copy (B);
  if (pnl_mat_ls_mat (A, B) != OK) printf("error in pnl_mat_ls\n");
  AB = pnl_mat_mult_mat (A, B);
  pnl_test_mat_eq_abs( AB, Bcopy, 1E-8, "mat_ls_mat", "");

  bcopy = pnl_vect_copy (b);
  if (pnl_mat_ls (A, b) != OK) printf("error in pnl_mat_ls\n");
  Ab = pnl_mat_mult_vect (A, b);
  pnl_test_vect_eq_abs( Ab, bcopy, 1E-8, "mat_ls", "");

  pnl_mat_free (&A);
  pnl_mat_free (&AB);
  pnl_mat_free (&B);
  pnl_mat_free (&Bcopy);
  pnl_vect_free (&b);
  pnl_vect_free (&bcopy);
}

static void pnl_mat_qr_test ()
{
  PnlMat *A, *Ap, *Q, *R, *QR;
  PnlPermutation *p;
  int m, gen;
  m = 5;
  gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  p = pnl_permutation_create (m);
  A = pnl_mat_create (m, m);
  Ap = pnl_mat_create (m, m);
  Q = pnl_mat_create (m, m);
  R = pnl_mat_create (m, m);
  QR = pnl_mat_create (m, m);
  pnl_rand_init (gen, m, m);
  pnl_mat_rand_uni2 (A, m, m, 0, 1, gen);

  pnl_mat_qr (Q, R, NULL, A);
  pnl_mat_mult_mat_inplace (QR, Q, R);
  pnl_test_mat_eq_abs(QR, A, 1E-8, "QR", "");

  pnl_mat_qr (Q, R, p, A);
  pnl_mat_col_permute (Ap, A, p);
  pnl_mat_mult_mat_inplace (QR, Q, R);
  pnl_test_mat_eq_abs(QR, Ap, 1E-8, "QR with pivoting", "");

  pnl_mat_free (&A);
  pnl_mat_free (&Ap);
  pnl_mat_free (&Q);
  pnl_mat_free (&R);
  pnl_mat_free (&QR);
  pnl_permutation_free (&p);
}


int main (int argc, char *argv[])
{
  pnl_test_init(argc, argv);
  pnl_mat_eigen_test ();
  pnl_mat_log_test ();
  pnl_mat_ls_test ();
  pnl_mat_qr_test ();
  exit (pnl_test_finalize("Lapack wrappers"));
}

