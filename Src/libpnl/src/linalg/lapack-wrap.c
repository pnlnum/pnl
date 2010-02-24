
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "config.h"
#include "pnl_internals.h"
#include "pnl_mathtools.h"
#include "pnl_matrix.h"
#include "pnl_vector.h"
#include "pnl_machine.h"

extern double pnl_dlamch (char *cmach);

static int pnl_mat_is_sym (const PnlMat *A);

/* lapack headers */
extern int C2F(dgetrf) (int *m, int *n, double *A, int *lda, int *ipvi, int *info);
extern int C2F(dgeev) (char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr, double *wi, double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);
extern int C2F(dsyev) (char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);
extern int C2F(dgelsy)(int *m, int *n, int *nrhs, double *A, int *lda, double *B, int *ldb, int *jpvt, double *rcond, int *rank, double *work, int *lwork, int *info);
extern int C2F(dgeqrf) (int *m, int *n, double *A, int *lda, double *tau, double *work, int *lwork, int *info);
extern int C2F(dorgqr) (int *m, int *n, int *k, double *Q, int *lda, double *tau, double *work, int *lwork, int *info);
extern int C2F(dgeqp3) (int *m, int *n, double *A, int *lda, int *jpvt, double *tau, double *work, int *lwork, int *info);
                      
/* Pnl wrappers */
static int pnl_dgeev (PnlVect *v, PnlMat *P, const PnlMat *A, int with_eigenvectors);
static int pnl_dsyev (PnlVect *v, PnlMat *P, const PnlMat *A, int with_eigenvectors);


/**
 * Checks if a (real) matrix is symmetric
 * @param A a real matrix
 * @return TRUE or FALSE
 */
static int pnl_mat_is_sym (const PnlMat *A)
{
  int i, j;

  if (A->m != A->n) return FALSE;
  for ( i=0 ; i<A->m ; i++ )
    for ( j=0 ; j<i ; j++ )
      {
        if (MGET(A, i, j) != MGET(A, j, i)) return FALSE;
      }
  return TRUE;
}

/**
 * Puts 0 in the lower triangular part of a square matrix
 * @param A a real matrix
 */
static void pnl_mat_make_upper (PnlMat *A)
{
  int i, j;
  for ( i=1 ; i<A->m ; i++ )
    {
      for ( j=0 ; j<i ; j++ )
        {
          pnl_mat_set (A, i, j, 0.);
        }
    }
}

/**
 * computes a P A = LU factoristion. On exit A contains the L and U
 * matrices. Note that the diagonal elemets of L are all 1.
 *
 * @param A the matrix to decompose.
 * @param p a PnlPermutation.
 */
void pnl_mat_lu (PnlMat *A, PnlPermutation *p)
{
  int info, i, N = A->n;
  int *invpiv;
  
  CheckIsSquare(A);
  invpiv = MALLOC_INT(N);
  pnl_permutation_init (p);

  pnl_mat_sq_transpose (A);
  C2F(dgetrf) (&N, &N, A->array, &N, invpiv, &info);
  if ( info != 0 )
    {
      PNL_ERROR ("LU decomposition cannot be computed", "pnl_mat_lu");
    }
  pnl_mat_sq_transpose (A);
  /* Fortran indices start at 1 */
  for ( i=0 ; i<N ; i++ ) invpiv[i]--;
  /* the permutation is computed the other round */
  for ( i=0 ; i<N ; i++ )
    {
      int ipiv = invpiv[i];
      if (ipiv != i)
        {
          int tmp = p->array[i] ;
          p->array[i] = p->array[ipiv];
          p->array[ipiv] = tmp;
        }
    }
  free (invpiv);
}

/**
 * computes a P A = LU factoristion. On exit A contains the L and U
 * matrices. Note that the diagonal elemets of L are all 1.
 *
 * @param Q an orthogonal matrix on exit
 * @param R an upper triangular matrix on exit
 * @param p a PnlPermutation. If p is NULL no permutation is computed.
 * @param A the matrix to decompose. PA = QR
 *
 */
void pnl_mat_qr (PnlMat *Q, PnlMat *R, PnlPermutation *p, const PnlMat *A)
{
  double *tau, *work, qlwork;
  int     lwork, info, m, i;
  CheckIsSquare (A);

  tau = NULL;
  work = NULL;
  m = A->m;
  pnl_mat_clone (R, A);
  pnl_mat_sq_transpose (R);
  tau=MALLOC_DOUBLE(m);
  
  if ( p == NULL )
    {
      lwork = -1;
      C2F(dgeqrf) (&m, &m, R->array, &m, tau, &qlwork, &lwork, &info);
      if ( info != 0 )
        {
          PNL_ERROR ("QR decomposition cannot be computed", "pnl_mat_qr");
        }
      lwork = (int) qlwork;
      work=MALLOC_DOUBLE(lwork);
      C2F(dgeqrf) (&m, &m, R->array, &m, tau, work, &lwork, &info);
      if ( info != 0 )
        {
          PNL_ERROR ("QR decomposition cannot be computed", "pnl_mat_qr");
        }
    }
  else
    {
      for (i = 0; i < m ; i++) p->array[i]= 0;
      lwork = -1;
      C2F(dgeqp3)(&m, &m, R->array, &m, p->array, tau, &qlwork, &lwork, &info);
      if ( info != 0 )
        {
          PNL_ERROR ("QR decomposition cannot be computed", "pnl_mat_qr");
        }
      lwork = (int) qlwork;
      work=MALLOC_DOUBLE(lwork);
      C2F(dgeqp3)(&m, &m, R->array, &m, p->array, tau, work, &lwork, &info);
      if ( info != 0 )
        {
          PNL_ERROR ("QR decomposition cannot be computed", "pnl_mat_qr");
        }
    }

  /* extract Q */
  pnl_mat_clone (Q, R);
  C2F(dorgqr)(&m, &m, &m, Q->array, &m, tau, work, &lwork, &info);
  if ( info != 0 )
    {
      PNL_ERROR ("QR decomposition cannot be computed", "pnl_mat_qr");
    }


  /* reverse the matrix storage */
  pnl_mat_sq_transpose (R);
  pnl_mat_sq_transpose (Q);

  /* extract R */
  pnl_mat_make_upper (R);

  
  FREE (work);
  FREE (tau);
}


/**
 * Computes the eigenvalues and eigenvectors of a real matrix
 *
 * @param v a vector containing the eigenvalues on exit
 * @param P a matrix containing the eigenvectors on exit
 * @param A a matrix
 * @param with_eigenvectors can be TRUE to compute the eigenvectors or FALSE
 * if they are nt required, in this latter case P can be NULL
 */
void pnl_mat_eigen (PnlVect *v, PnlMat *P, const PnlMat *A, int with_eigenvectors)
{
  int is_sym;
  int info;

  is_sym = pnl_mat_is_sym (A);
  if (is_sym == TRUE) info = pnl_dsyev(v, P, A, with_eigenvectors);
  else info = pnl_dgeev(v, P, A, with_eigenvectors);

  if (info == FAIL)
    {
      PNL_ERROR ("Error", "pnl_mat_eigen");
    }
}

/**
 * Wrapper to dgeev
 * Computes the eigenvalues and eigenvectors of a real symmetric matrix
 *
 * @param v a vector containing the eigenvalues on exit
 * @param P a matrix containing the eigenvectors on exit (P is orthonormal)
 * @param A a real symmetric matrix
 * @param with_eigenvectors can be TRUE to compute the eigenvectors or FALSE
 * if they are nt required, in this latter case P can be NULL
 */
static int pnl_dsyev (PnlVect *v, PnlMat *P, const PnlMat *A, int with_eigenvectors)
{
  int n=A->n;
  int info, lwork;
  double *work=NULL;
  double qlwork;

  /* Clone A, because dsyev modifies its input argument */
  pnl_mat_clone (P, A);
  pnl_vect_resize (v, n);
  
  lwork = -1;
  C2F(dsyev)((with_eigenvectors==TRUE)?"V":"N", "L", &n, P->array, &n, v->array, 
             &qlwork, &lwork, &info);
  lwork = (int) qlwork;
  if ( (work=MALLOC_DOUBLE(lwork)) == NULL ) goto err;
  C2F(dsyev)((with_eigenvectors==TRUE)?"V":"N", "L", &n, P->array, &n, v->array, 
             work, &lwork, &info);

  if (info != 0) 
    {
      printf("Error: convergence problem in dsyev\n");
      goto err;
    }

  /* extract results */ 
  pnl_mat_sq_transpose (P);

  free(work);
  return OK;

 err:
  free(work);  
  return FAIL;
}

/**
 * Wrapper to dgeev
 * Computes the eigenvalues and eigenvectors of a real non-symmetric matrix
 *
 * @param v a vector containing the eigenvalues on exit
 * @param P a matrix containing the eigenvectors on exit
 * @param A a real non-symmetric matrix
 * @param with_eigenvectors can be TRUE to compute the eigenvectors or FALSE
 * if they are nt required, in this latter case P can be NULL
 */
static int pnl_dgeev (PnlVect *v, PnlMat *P, const PnlMat *A, int with_eigenvectors)
{
  int n=A->n;
  int info, lwork,i;
  double *work=NULL, *wi=NULL, *input=NULL;
  double qlwork;

  wi = MALLOC_DOUBLE(n);
  input = MALLOC_DOUBLE(A->mn);
  if ( wi == NULL || input == NULL ) goto err;
  if ( with_eigenvectors == TRUE ) { pnl_mat_resize (P, n, n); }
  pnl_vect_resize (v, n);
  /* Copy A, because dgeev modifies its input argument */
  memcpy (input, A->array, A->mn * sizeof(double));
  
  
  lwork = -1;
  C2F(dgeev)((with_eigenvectors==TRUE)?"V":"N", "N", &n, input, &n, v->array, wi,
             (with_eigenvectors==FALSE)?NULL:P->array, &n, NULL, &n,
             &qlwork, &lwork, &info);
  lwork = (int) qlwork;
  if ( (work=MALLOC_DOUBLE(lwork)) == NULL ) goto err;
  C2F(dgeev)((with_eigenvectors==FALSE)?"N":"V", "N", &n, input, &n, v->array, wi,
             (with_eigenvectors==FALSE)?NULL:P->array, &n, NULL, &n,
             work, &lwork, &info);

  if (info != 0) 
    {
      printf("Error: convergence problem in dgeev\n");
      goto err;
    }

  /* Check that result is real */ 
  for (i = 0 ; i < n ; i++)
    {
      if (wi[i] != 0.) { printf ("Some eigenvalues are complex\n"); goto err;}
    }

  pnl_mat_sq_transpose (P);

  free(wi); free(work); free(input);
  return OK;

 err:
  free(wi); free(work); free(input);
  return FAIL;
}

/**
 * Matrix logarithm B = log( A).
 *
 * Note : A must be diagonalizable.
 *
 * @param A a real diagonalizable matrix
 * @param B contains log(A) on return
 */
void pnl_mat_log (PnlMat *B, const PnlMat *A)
{
  int n, i, j;
  PnlMat *P, *invP;
  PnlVect *D;
  CheckIsSquare (A);
  n = A->n;

  P = pnl_mat_create(n,n);
  invP = pnl_mat_create(n,n);
  D = pnl_vect_create(n);

  pnl_mat_eigen(D, P, A, TRUE);

  for ( i=0 ; i<n ; i++ )
    {
      if (pnl_vect_get(D, i) <= 0)
        {
          PNL_ERROR ("Negative eigenvalues", "pnl_mat_log");
        }
    }

  /* Compute inv(P). If P is not invertible, it means that the matrix A is
     not diagonalizable.
  */
  pnl_mat_inverse (invP, P);
  
  /* compute P = P * diag(log(D)) */
  for ( i=0 ; i<n ; i++ )
    {
      for ( j=0 ; j<n ; j++ )
        {
          pnl_mat_set ( P, i, j, pnl_mat_get (P, i, j) * log (pnl_vect_get (D, j)) );
        }
    }

  /* compute B = P * inv(P) */
  pnl_mat_mult_mat_inplace (B, P, invP);

  pnl_mat_free (&P);
  pnl_mat_free (&invP);
  pnl_vect_free (&D);
}

/**
 * Solves A * X = B in the least square sense
 *
 * A an m x n matrix with m >= n
 * X on exit is an n x nrhs matrix
 * B an m x nrhs matrix
 *
 * @return FAIL or OK
 */
int pnl_mat_ls_mat (const PnlMat *A, PnlMat *B)
{
  int m ,n ,nrhs, lda, ldb, rank, info, lwork, i, j;
  double *work, qwork, rcond, *X;
  int *jpvt;
  PnlMat *tA;

  tA = pnl_mat_transpose (A);
  
  m = A->m; n = A->n;
  nrhs = B->n;
  lda = A->m;
  ldb = MAX(A->m, A->n);
  jpvt = NULL;
  work = NULL;
  rcond = MAX(A->m,A->n) * pnl_dlamch("eps");

  /* En large B to contain the solution and  transpose it because
     of Blas storage.
     X is matrix of size nrhs x ldb 
  */
  if ( (X = MALLOC_DOUBLE(ldb * nrhs)) == NULL) return FAIL;
  for ( i=0 ; i<B->m ; i++ )
    for ( j=0 ; j<B->n ; j++ )
      {
        X[j*ldb+i] = MGET (B, i, j);
      }
  

  lwork = -1;
  C2F(dgelsy) (&m, &n, &nrhs, tA->array, &lda, X, &ldb, jpvt, &rcond, &rank, &qwork, &lwork, &info);

  if ( info != 0 ) goto err;
  lwork = (int) qwork;
  if ((work = MALLOC_DOUBLE(lwork))==NULL || (jpvt = MALLOC_INT(n))==NULL ) goto err;
  for ( i=0 ; i<n ; i++ ) jpvt[i] = 0;

  C2F(dgelsy) (&m, &n, &nrhs, tA->array, &lda, X, &ldb, jpvt, &rcond, &rank, work, &lwork, &info);

  pnl_mat_resize (B, A->n, nrhs);
  for ( i=0 ; i<A->n ; i++ )
    for ( j=0 ; j<nrhs ; j++ )
      {
        MLET(B, i, j) = X[j*ldb+i];
      }
  
  FREE (work); FREE (jpvt); FREE(X);
  pnl_mat_free (&tA);
  return OK;

 err:
  pnl_mat_free (&tA);
  FREE (work); FREE(jpvt);
  return FAIL;
}

/**
 * Solves A * x = b in the least square sense
 *
 * A an m x n matrix with m >= n
 * X on exit is an vector of size n
 * B is a vector of size m
 *
 * @return FAIL or OK
 */
int pnl_mat_ls (const PnlMat *A, PnlVect *b)
{
  PnlMat *B;
  int status;
  B = pnl_mat_create_from_ptr (b->size, 1, b->array);

  status = pnl_mat_ls_mat (A, B);
  pnl_vect_resize (b, A->n);
  memcpy (b->array, B->array, b->size*sizeof(double));

  pnl_mat_free (&B);
  return status;
}


