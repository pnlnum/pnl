
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

#include "pnl/pnl_config.h"
#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_machine.h"
#include "pnl/clapack.h"

extern double pnl_dlamch (char *cmach);

static int pnl_mat_is_sym (const PnlMat *A);

/* Pnl wrappers */
static int pnl_dgeev (PnlVect *v, PnlMat *P, const PnlMat *A, int with_eigenvectors);
static int pnl_dsyev (PnlVect *v, PnlMat *P, const PnlMat *A, int with_eigenvectors);

#define BASE_DOUBLE
#include "pnl/pnl_templates_on.h"
#include "lapack-wrap_source.c"
#include "pnl/pnl_templates_off.h"
#undef  BASE_DOUBLE


#define BASE_PNL_COMPLEX
#include "pnl/pnl_templates_on.h"
#include "lapack-wrap_source.c"
#include "pnl/pnl_templates_off.h"
#undef  BASE_PNL_COMPLEX


/**
 * Check if a (real) matrix is symmetric
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
        if (PNL_MGET(A, i, j) != PNL_MGET(A, j, i)) return FALSE;
      }
  return TRUE;
}

/**
 * Check if a (complex) matrix is Hermitian
 * @param A a real matrix
 * @return TRUE or FALSE
 */
static int pnl_mat_complex_is_he (const PnlMatComplex *A)
{
  int i, j;

  if (A->m != A->n) return FALSE;
  for ( i=0 ; i<A->m ; i++ )
    {
      const dcomplex Aii = PNL_MGET(A, i, i);
      if ( Aii.i != 0. ) return FALSE;
      for ( j=0 ; j<i ; j++ )
        {
          const dcomplex Aij = PNL_MGET(A,i,j);
          const dcomplex conj_Aji = Conj(PNL_MGET(A, j, i));
          if ( (Aij.r != conj_Aji.r) || (Aij.i != conj_Aji.i) ) return FALSE;
        }
    }
  return TRUE;
}

/**
 * Put 0 in the lower triangular part of a square matrix
 * @param A a real matrix
 */
static void pnl_mat_make_upper (PnlMat *A)
{
  int i, j;
  for ( i=1 ; i<A->m ; i++ )
    {
      for ( j=0 ; j<i ; j++ )
        {
          PNL_MLET (A, i, j) = 0.;
        }
    }
}

/**
 * Compute a A P = QR factorisation. 
 *
 * @param Q an orthonormal matrix on exit
 * @param R an upper triangular matrix on exit
 * @param p a PnlPermutation. If p is NULL no permutation is computed.
 * @param A the matrix to decompose. PA = QR
 * @return OK or FAIL
 *
 */
int pnl_mat_qr (PnlMat *Q, PnlMat *R, PnlPermutation *p, const PnlMat *A)
{
  double *tau, *work, qlwork;
  int     lwork, info, m, i;
  CheckIsSquare (A);

  tau = NULL;
  work = NULL;
  m = A->m;
  pnl_mat_clone (R, A);
  /* Convert to column wise storage */
  pnl_mat_sq_transpose (R);
  tau=MALLOC_DOUBLE(m);
  
  if ( p == NULL )
    {
      lwork = -1;
      C2F(dgeqrf) (&m, &m, R->array, &m, tau, &qlwork, &lwork, &info);
      if ( info != 0 )
        {
          PNL_MESSAGE_ERROR ("QR decomposition cannot be computed", "pnl_mat_qr");
          return FAIL;
        }
      lwork = (int) qlwork;
      work=MALLOC_DOUBLE(lwork);
      C2F(dgeqrf) (&m, &m, R->array, &m, tau, work, &lwork, &info);
      if ( info != 0 )
        {
          PNL_MESSAGE_ERROR ("QR decomposition cannot be computed", "pnl_mat_qr");
          return FAIL;
        }
    }
  else
    {
      pnl_vect_int_resize (p, m);
      for (i = 0; i < m ; i++) p->array[i]= 0;
      lwork = -1;
      C2F(dgeqp3)(&m, &m, R->array, &m, p->array, tau, &qlwork, &lwork, &info);
      if ( info != 0 )
        {
          PNL_MESSAGE_ERROR ("QR decomposition cannot be computed", "pnl_mat_qr");
          return FAIL;
        }
      lwork = (int) qlwork;
      work=MALLOC_DOUBLE(lwork);
      C2F(dgeqp3)(&m, &m, R->array, &m, p->array, tau, work, &lwork, &info);
      if ( info != 0 )
        {
          PNL_MESSAGE_ERROR ("QR decomposition cannot be computed", "pnl_mat_qr");
          return FAIL;
        }
      for (i = 0; i < m ; i++) p->array[i]--;
    }

  /* extract Q */
  pnl_mat_clone (Q, R);
  C2F(dorgqr)(&m, &m, &m, Q->array, &m, tau, work, &lwork, &info);
  if ( info != 0 )
    {
      PNL_MESSAGE_ERROR ("QR decomposition cannot be computed", "pnl_mat_qr");
      return FAIL;
    }


  /* Revert to row wise storage */
  pnl_mat_sq_transpose (R);
  pnl_mat_sq_transpose (Q);

  /* extract R */
  pnl_mat_make_upper (R);

  
  FREE (work);
  FREE (tau);
  return OK;
}

/**
 * Solve the linear system A x = b with A P = QR.  *
 * @param x a PnlVect containing the solution on exit
 * @param Q a orthogonal PnlMat 
 * @param R an upper triagular PnlMat 
 * @param p a PnlVectInt (permutation vector)
 * @param b right hand side member
 * @return OK or FAIL
 */
int pnl_mat_qr_syslin (PnlVect *x, const PnlMat *Q, const PnlMat *R,
                       const PnlVectInt *p, const PnlVect *b)
{
  PnlVect *y;
  CheckIsSquare(Q);
  CheckIsSquare(R);
  CheckMatVectIsCompatible (R, b);
  CheckVectMatch (p, b); 
  y = pnl_vect_new ();

  pnl_mat_dgemv ('T', 1., Q, b, 0., y); /* y = Q' b */
  pnl_mat_upper_syslin (x, R, y); /* x = R^-1 Q' b */
  pnl_vect_permute_inverse_inplace (x, p);
  pnl_vect_free (&y);
  return OK;
}

/**
 * Wrapper to dgeev
 * Compute the eigenvalues and eigenvectors of a real symmetric matrix
 *
 * @param v a vector containing the eigenvalues on exit
 * @param P a matrix containing the eigenvectors on exit (P is orthonormal)
 * @param A a real symmetric matrix
 * @param with_eigenvectors can be TRUE to compute the eigenvectors or FALSE
 * if they are not required, in this latter case P can be NULL
 * @return OK or FAIL
 */
static int pnl_dsyev (PnlVect *v, PnlMat *P, const PnlMat *A, int with_eigenvectors)
{
  int n=A->n;
  int info, lwork;
  double *work=NULL;
  double qlwork;
  PnlMat *Acopy;

  /* Clone A, because dsyev modifies its input argument */
  if ( with_eigenvectors == FALSE)
    {
      Acopy = pnl_mat_copy (A);
    }
  else
    {
      pnl_mat_clone (P, A);
      Acopy = P;
    }
  pnl_vect_resize (v, n);
  
  lwork = -1;
  C2F(dsyev)((with_eigenvectors==TRUE)?"V":"N", "L", &n, Acopy->array, &n, v->array, 
             &qlwork, &lwork, &info);
  lwork = (int) qlwork;
  if ( (work=MALLOC_DOUBLE(lwork)) == NULL ) goto err;
  C2F(dsyev)((with_eigenvectors==TRUE)?"V":"N", "L", &n, Acopy->array, &n, v->array, 
             work, &lwork, &info);

  if (info != 0) 
    {
      printf("Error: convergence problem in dsyev\n");
      goto err;
    }

  /* Revert to row wise storage */
  if ( with_eigenvectors ) pnl_mat_sq_transpose (P);


  free(work);
  if ( with_eigenvectors == FALSE) pnl_mat_free (&Acopy);
  return OK;

 err:
  free(work);  
  if ( with_eigenvectors == FALSE) pnl_mat_free (&Acopy);
  return FAIL;
}

/**
 * Wrapper to dgeev
 * Compute the eigenvalues and eigenvectors of a real non symmetric matrix
 *
 * @param v a vector containing the eigenvalues on exit
 * @param P a matrix containing the eigenvectors on exit
 * @param A a real non symmetric matrix
 * @param with_eigenvectors can be TRUE to compute the eigenvectors or FALSE
 * if they are not required, in this latter case P can be NULL
 * @return OK or FAIL
 */
static int pnl_dgeev (PnlVect *v, PnlMat *P, const PnlMat *A, int with_eigenvectors)
{
  int n=A->n;
  int info, lwork,i;
  double *work=NULL, *wi=NULL;
  double qlwork;
  PnlMat *tA;

  wi = MALLOC_DOUBLE(n);
  if ( wi == NULL ) goto err;
  if ( with_eigenvectors == TRUE ) { pnl_mat_resize (P, n, n); }
  pnl_vect_resize (v, n);
  /* 
   * Copy A in a transpose form, because zgeev modifies its input argument
   * and Lapack uses column wise storage.
   */
  tA = pnl_mat_transpose (A);
  
  lwork = -1;
  C2F(dgeev)("N", (with_eigenvectors==TRUE)?"V":"N",  &n, tA->array, &n, v->array, wi,
               NULL, &n, (with_eigenvectors==FALSE)?NULL:P->array, &n,
             &qlwork, &lwork, &info);
  lwork = (int) qlwork;
  if ( (work=MALLOC_DOUBLE(lwork)) == NULL ) goto err;
  C2F(dgeev)("N", (with_eigenvectors==FALSE)?"N":"V",  &n, tA->array, &n, v->array, wi,
              NULL, &n, (with_eigenvectors==FALSE)?NULL:P->array, &n,
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

  /* Convert to row wise storage */
  if ( with_eigenvectors ) pnl_mat_sq_transpose (P);

  free(wi); free(work);
  pnl_mat_free (&tA);
  return OK;

 err:
  free(wi); free(work);
  pnl_mat_free (&tA);
  return FAIL;
}

/**
 * Compute the eigenvalues and eigenvectors of a real matrix
 *
 * @param v a vector containing the eigenvalues on exit
 * @param P a matrix containing the eigenvectors on exit
 * @param A a matrix
 * @param with_eigenvectors can be TRUE to compute the eigenvectors or FALSE
 * if they are not required, in this latter case P can be NULL
 * @return OK or FAIL
 */
int pnl_mat_eigen (PnlVect *v, PnlMat *P, const PnlMat *A, int with_eigenvectors)
{
  int is_sym;
  int info;

  is_sym = pnl_mat_is_sym (A);
  if (is_sym == TRUE) info = pnl_dsyev(v, P, A, with_eigenvectors);
  else info = pnl_dgeev(v, P, A, with_eigenvectors);

  if (info == FAIL)
    {
      PNL_MESSAGE_ERROR ("Error", "pnl_mat_eigen");
      return FAIL;
    }
  return OK;
}

/**
 * Matrix logarithm B = log( A).
 *
 * Note : A must be diagonalizable.
 *
 * @param A a real diagonalizable matrix
 * @param B contains log(A) on return
 * @return OK or FAIL
 */
int pnl_mat_log (PnlMat *B, const PnlMat *A)
{
  int n, i, j;
  PnlMat *P, *invP;
  PnlVect *D;
  CheckIsSquare (A);
  n = A->n;

  P = pnl_mat_create(n,n);
  D = pnl_vect_create(n);

  if ( pnl_mat_eigen(D, P, A, TRUE) != OK )
    {
      pnl_mat_free (&P);
      pnl_vect_free (&D);
      return FAIL;
    }

  for ( i=0 ; i<n ; i++ )
    {
      if (PNL_GET(D, i) <= 0)
        {
          PNL_MESSAGE_ERROR ("Negative eigenvalues", "pnl_mat_log");
          return FAIL;
        }
    }

  /* Compute inv(P). If P is not invertible, it means that the matrix A is
     not diagonalizable.
  */
  invP = pnl_mat_create(n,n);
  if ( pnl_mat_inverse (invP, P) != OK )
    {
      PNL_MESSAGE_ERROR ("matrix is not diagonalizable", "pnl_mat_log");
      return FAIL;
    }
  
  /* compute P = P * diag(log(D)) */
  for ( i=0 ; i<n ; i++ )
    {
      for ( j=0 ; j<n ; j++ )
        {
          PNL_MLET ( P, i, j) = PNL_MGET (P, i, j) * log (PNL_GET (D, j));
        }
    }

  /* compute B = P * inv(P) */
  pnl_mat_mult_mat_inplace (B, P, invP);

  pnl_mat_free (&P);
  pnl_mat_free (&invP);
  pnl_vect_free (&D);
  return OK;
}

/**
 * Wrapper to zhesv
 * Compute the eigenvalues and eigenvectors of an Hermitian matrix
 *
 * @param v a vector containing the eigenvalues on exit
 * @param P a matrix containing the eigenvectors on exit (P is orthonormal)
 * @param A a complex Hermitian matrix
 * @param with_eigenvectors can be TRUE to compute the eigenvectors or FALSE
 * if they are not required, in this latter case P can be NULL
 * @return OK or FAIL
 */
static int pnl_zheev (PnlVectComplex *v, PnlMatComplex *P, const PnlMatComplex *A, int with_eigenvectors)
{
  int n=A->n;
  int i, info, lwork;
  double *rwork = NULL, *w=NULL;
  dcomplex *work=NULL;
  dcomplex qlwork;

  rwork = MALLOC_DOUBLE(MAX(1, 2 * n - 1));
  w = MALLOC_DOUBLE(n);
  /* Clone A, because dsyev modifies its input argument */
  pnl_mat_complex_clone (P, A);
  pnl_vect_complex_resize (v, n);
  
  lwork = -1;
  C2F(zheev)((with_eigenvectors==TRUE)?"V":"N", "L", &n, P->array, &n, w,
             &qlwork, &lwork, rwork, &info);
  lwork = (int) qlwork.r;
  if ( (work=MALLOC_COMPLEX(lwork)) == NULL ) goto err;
  C2F(zheev)((with_eigenvectors==TRUE)?"V":"N", "L", &n, P->array, &n, w,
             work, &lwork, rwork, &info);

  if (info != 0) 
    {
      printf("Error: convergence problem in zhesv\n");
      goto err;
    }

  for ( i=0 ; i<n ; i++ )
    {
      PNL_LET (v, i) = Complex (w[i], 0.);
    }
  /* Revert to row wise storage */
  pnl_mat_complex_sq_transpose (P);

  free(work); free (rwork); free(w);
  return OK;

 err:
  free(work); free (rwork); free(w);
  return FAIL;
}

/**
 * Wrapper to dgeev
 * Compute the eigenvalues and eigenvectors of a complex non Hermitian
 * matrix
 *
 * @param v a vector containing the eigenvalues on exit
 * @param P a matrix containing the eigenvectors on exit
 * @param A a complex matrix
 * @param with_eigenvectors can be TRUE to compute the eigenvectors or FALSE
 * if they are not required, in this latter case P can be NULL
 * @return OK or FAIL
 */
static int pnl_zgeev (PnlVectComplex *v, PnlMatComplex *P, const PnlMatComplex *A, int with_eigenvectors)
{
  int n=A->n;
  int info, lwork;
  double *rwork;
  dcomplex *work=NULL;
  dcomplex qlwork;
  PnlMatComplex *tA;

  rwork = MALLOC_DOUBLE(2*n);
  if ( with_eigenvectors == TRUE ) { pnl_mat_complex_resize (P, n, n); }
  pnl_vect_complex_resize (v, n);
  /* 
   * Copy A in a transpose form, because zgeev modifies its input argument
   * and Lapack uses column wise storage.
   */
  tA = pnl_mat_complex_transpose (A);
  
  
  lwork = -1;
  C2F(zgeev)("N", (with_eigenvectors==TRUE)?"V":"N", &n, tA->array, &n, v->array, 
              NULL, &n, (with_eigenvectors==FALSE)?NULL:P->array, &n,
             &qlwork, &lwork, rwork, &info);
  lwork = (int) qlwork.r;
  if ( (work=MALLOC_COMPLEX(lwork)) == NULL ) goto err;
  C2F(zgeev)("N", (with_eigenvectors==FALSE)?"N":"V", &n, tA->array, &n, v->array, 
             NULL, &n, (with_eigenvectors==FALSE)?NULL:P->array, &n, 
             work, &lwork, rwork, &info);

  if (info != 0) 
    {
      printf("Error: convergence problem in dgeev\n");
      goto err;
    }

  /* Convert to row wise storage */
  pnl_mat_complex_sq_transpose (P);

  free(work); free(rwork); 
  pnl_mat_complex_free (&tA);
  return OK;

 err:
  free(work); free(rwork);
  pnl_mat_complex_free (&tA);
  return FAIL;
}

/**
 * Compute the eigenvalues and eigenvectors of a real matrix
 *
 * @param v a vector containing the eigenvalues on exit
 * @param P a matrix containing the eigenvectors on exit
 * @param A a matrix
 * @param with_eigenvectors can be TRUE to compute the eigenvectors or FALSE
 * if they are not required, in this latter case P can be NULL
 * @return OK or FAIL
 */
int pnl_mat_complex_eigen (PnlVectComplex *v, PnlMatComplex *P, const PnlMatComplex *A, int with_eigenvectors)
{
  int is_sym;
  int info;

  is_sym = pnl_mat_complex_is_he (A);
  if (is_sym == TRUE) info = pnl_zheev(v, P, A, with_eigenvectors);
  else info = pnl_zgeev(v, P, A, with_eigenvectors);

  if (info == FAIL)
    {
      PNL_MESSAGE_ERROR ("Error", "pnl_mat_eigen");
      return FAIL;
    }
  return OK;
}

/**
 * Matrix logarithm B = log( A).
 *
 * Note : A must be diagonalizable.
 *
 * @param A a real diagonalizable matrix
 * @param B contains log(A) on return
 * @return OK or FAIL
 */
int pnl_mat_complex_log (PnlMatComplex *B, const PnlMatComplex *A)
{
  int n, i, j;
  PnlMatComplex *P, *invP;
  PnlVectComplex *D;
  CheckIsSquare (A);
  n = A->n;

  P = pnl_mat_complex_create(n,n);
  D = pnl_vect_complex_create(n);

  if ( pnl_mat_complex_eigen(D, P, A, TRUE) != OK )
    {
      pnl_mat_complex_free (&P);
      pnl_vect_complex_free (&D);
      return FAIL;
    }

  for ( i=0 ; i<n ; i++ )
    {
      dcomplex di = PNL_GET(D, i);
      if ( di.r == 0. && di.i == 0. )
        {
          PNL_MESSAGE_ERROR ("Zero eigenvalue detected", "pnl_mat_complex_log");
          return FAIL;
        }
    }

  /* Compute inv(P). If P is not invertible, it means that the matrix A is
     not diagonalizable.
  */
  invP = pnl_mat_complex_create(n,n);
  if ( pnl_mat_complex_inverse (invP, P) != OK )
    {
      PNL_MESSAGE_ERROR ("matrix is not diagonalizable", "pnl_mat_log");
      return FAIL;
    }
  
  /* compute P = P * diag(log(D)) */
  for ( i=0 ; i<n ; i++ )
    {
      for ( j=0 ; j<n ; j++ )
        {
          PNL_MLET ( P, i, j) = Cmul (PNL_MGET (P, i, j), Clog (PNL_GET (D, j)));
        }
    }

  /* compute B = P * inv(P) */
  pnl_mat_complex_mult_mat_inplace (B, P, invP);

  pnl_mat_complex_free (&P);
  pnl_mat_complex_free (&invP);
  pnl_vect_complex_free (&D);
  return OK;
}

/**
 * Solve A * X = B in the least square sense
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

  /* Convert to column wise storage */
  tA = pnl_mat_transpose (A);
  
  m = A->m; n = A->n;
  nrhs = B->n;
  lda = A->m;
  ldb = MAX(A->m, A->n);
  jpvt = NULL;
  work = NULL;
  rcond = MAX(A->m,A->n) * pnl_dlamch("eps");

  /* En large B to contain the solution and  transpose it because
     of column wise storage.
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
 * Solve A * x = b in the least square sense
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


