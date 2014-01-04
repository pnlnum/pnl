
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
  pnl_test_mat_eq_abs(V, B, 1E-8, "eigen sym", "");
  pnl_mat_free (&V);
  
  pnl_mat_eigen (v, NULL, B, FALSE);
  V = pnl_mat_create_diag (v);
  pnl_mat_mult_mat_inplace (A, P, V);  /* P V */
  pnl_mat_mult_mat_inplace (V, A, Pinv); /* P V P^-1 */
  pnl_test_mat_eq_abs(V, B, 1E-8, "eigen sym without P", "");
  pnl_mat_free (&V);
  pnl_mat_free(&A);

  A = pnl_mat_create_diagonalizable (n);
  pnl_mat_eigen (v, P, A, TRUE);
  V = pnl_mat_create_diag (v);
  pnl_mat_inverse(Pinv, P);
  pnl_mat_mult_mat_inplace (B, P, V);  /* P V */
  pnl_mat_mult_mat_inplace (V, B, Pinv); /* P V P^-1 */
  pnl_test_mat_eq_abs(V, A, 1E-8, "eigen", "");
  pnl_mat_free (&V);

  pnl_mat_eigen (v, NULL, A, FALSE);
  V = pnl_mat_create_diag (v);
  pnl_mat_mult_mat_inplace (B, P, V);  /* P V */
  pnl_mat_mult_mat_inplace (V, B, Pinv); /* P V P^-1 */
  pnl_test_mat_eq_abs(V, A, 1E-8, "eigen witouht P", "");

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

static void pnl_mat_complex_exp_test ()
{
  PnlMatComplex *A, *expA, *sol;
  A = pnl_mat_complex_create_from_file ("Data/Acomplex.txt");
  sol = pnl_mat_complex_create_from_file ("Data/expAcomplex.txt");
  expA = pnl_mat_complex_new ();
  pnl_mat_complex_exp (expA, A);

  pnl_test_mat_complex_eq_abs(expA, sol, 1E-8, "complex expm", "");
  pnl_mat_complex_free (&A);
  pnl_mat_complex_free (&expA);
  pnl_mat_complex_free (&sol);
}

static void pnl_mat_complex_log_test ()
{
  PnlMatComplex *A, *logA, *explogA;
  A = pnl_mat_complex_create_from_file ("Data/expAcomplex.txt");
  logA = pnl_mat_complex_new ();
  explogA = pnl_mat_complex_new ();
  pnl_mat_complex_log (logA, A);
  pnl_mat_complex_exp (explogA, logA);

  pnl_test_mat_complex_eq_abs(explogA, A, 1E-8, "complex logm", "");
  pnl_mat_complex_free (&A);
  pnl_mat_complex_free (&logA);
  pnl_mat_complex_free (&explogA);
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
  pnl_vect_free (&Ab);
  pnl_vect_free (&b);
  pnl_vect_free (&bcopy);
}

static void pnl_mat_pchol_test ()
{
  PnlMat *S, *L, *PS, *PSP, *LL_t;
  PnlVectInt *p;
  int rank;
  double tol = 1E-9;
  int m = 5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init (gen, m, m);
  p = pnl_vect_int_create (m);
  S = pnl_mat_create (m, m);
  PS = pnl_mat_create (m, m);
  PSP = pnl_mat_create (m, m);
  pnl_mat_rand_uni2 (S, m, m, 0, 1, gen);
  L = pnl_mat_create (m, m);
  pnl_mat_dgemm ('T', 'N', 1, S, S, 0., L);
  pnl_mat_clone (S, L);
  pnl_mat_pchol (L, tol, &rank, p);

  LL_t = pnl_mat_create (m,m);
  pnl_mat_dgemm ('N', 'T', 1, L, L, 0., LL_t);

  pnl_mat_row_permute (PS, S, p);
  pnl_mat_col_permute (PSP, PS, p);

  pnl_test_mat_eq_abs(PSP, LL_t, 1E-8, "Cholesky with pivoting", "");

  pnl_mat_free (&S);
  pnl_mat_free (&PS);
  pnl_mat_free (&PSP);
  pnl_mat_free (&LL_t);
  pnl_mat_free (&L);
  pnl_vect_int_free (&p);
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
  pnl_mat_complex_exp_test ();
  pnl_mat_complex_log_test ();
  pnl_mat_ls_test ();
  pnl_mat_pchol_test ();
  pnl_mat_qr_test ();
  exit (pnl_test_finalize("Lapack wrappers"));
}

