
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
 * Cholesky decomposition. Postivity is checked during the
 * transformation, but no test of symmetry.
 *
 * Decomposition done in place. The lower part of the matrix
 * contains the cholesky decomposition. The upper part is
 * set to 0.
 *
 * @param M : a PnlMat pointer.
 * @return OK or FAIL
 */
int pnl_mat_chol (PnlMat *M)
{
  int n, lda, info, i, j;
  CheckIsSquare(M);
  n = M->m;
  lda = M->m;
  /* Because of Fortran column wise storage, we ask for an upper triangular
   * matrix to actually have a lower one */
  C2F(dpotrf)("U", &n, M->array, &lda, &info);
  if (info != 0)
    {
      PNL_MESSAGE_ERROR ("matrix is singular", "pnl_mat_chol");
      return FAIL;
    }
  /* Sets the upper part to 0 */
  for ( i=0 ; i<M->m ; i++ )
    {
      for ( j=i+1 ; j<M->n ; j++ )
        {
          PNL_MLET (M, i, j) = 0.;
        }
    }
  return OK;
}

/** 
 * Computes the Cholesky decompisition with complete pivoting. 
 *
 *       P' * A * P = L  * L'
 *
 * The input matrix must be symmetric and positive semidefinite 
 * 
 * @param M input matrix. On output contains L
 * @param tol tolerance : if a pivot is smaller than this value, it is
 * considered as zero
 * @param rank (output) rank of the matrix
 * @param p An integer vector representing a permutation
 * 
 * @return OK or FAIL
 */
int pnl_mat_pchol (PnlMat *M, double tol, int *rank, PnlVectInt *p)
{
  int n, lda, info, i, j;
  double *work;
  CheckIsSquare(M);
  n = M->m;
  lda = M->m;
  pnl_vect_int_resize (p, n);
  work = malloc (sizeof (double) * 2*n);


  /* Because of Fortran column wise storage, we ask for an upper triangular
   * matrix to actually have a lower one */
  C2F(dpstrf)("U", &n, M->array, &lda, p->array, rank, &tol, work, &info);
  free (work);
  if (info < 0)
    {
      PNL_MESSAGE_ERROR ("matrix has illegal entries", "pnl_mat_pchol");
      return FAIL;
    }
  /* Sets the upper part to 0 */
  for ( i=0 ; i<M->m ; i++ )
    {
      for ( j=i+1 ; j<M->n ; j++ )
        {
          PNL_MLET (M, i, j) = 0.;
        }
    }
  /* C indices start at 0 */
  pnl_vect_int_minus_int (p, 1);
  return OK;
}

/**
 * computes a P A = LU factorisation. On exit A contains the L and U
 * matrices. Note that the diagonal elements of L are all 1.
 *
 * @param A the matrix to decompose.
 * @param p a PnlPermutation.
 * @return OK or FAIL
 */
int pnl_mat_lu (PnlMat *A, PnlPermutation *p)
{
  int info, N = A->n;
  
  CheckIsSquare(A);

  pnl_mat_sq_transpose (A);
  C2F(dgetrf) (&N, &N, A->array, &N, p->array, &info);
  if ( info != 0 )
    {
      PNL_MESSAGE_ERROR ("LU decomposition cannot be computed", "pnl_mat_lu");
      return FAIL;
    }
  pnl_mat_sq_transpose (A);
  /* C indices start at 0 */
  pnl_vect_int_minus_int (p, 1);
  return OK;
}

/**
 * computes a A P = QR factorisation. 
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
 * solves an upper triangular linear system
 *
 * @param x already existing PnlVect that contains the solution on exit
 * @param A an upper triangular matric
 * @param b right hand side member
 * @return OK or FAIL
 */
int pnl_mat_upper_syslin (PnlVect *x, const PnlMat *A, const  PnlVect *b)
{
  int n, nrhs, lda, ldb, info;
  
  CheckIsSquare(A);
  n = A->n;
  nrhs = 1;
  lda = A->m;
  ldb = A->m;
  pnl_vect_clone (x, b);
  /* Beware that Fortran uses a column wise store, we actually consider A^T */
  C2F(dtrtrs)("L","T","N",&n,&nrhs,A->array,&lda,x->array,&ldb,&info);
  if (info != 0)
    {
      PNL_MESSAGE_ERROR ("Matrix is singular", "pnl_mat_upper_syslin");
      return FAIL;
    }
  return OK;
}

/**
 * solves a lower triangular linear system
 *
 * @param x already existing PnlVect that contains the solution on exit
 * @param A a lower triangular matrix
 * @param b right hand side member
 * @return OK or FAIL
 */
int pnl_mat_lower_syslin (PnlVect *x, const PnlMat *A, const  PnlVect *b)
{
  int n, nrhs, lda, ldb, info;
  
  CheckIsSquare(A);
  n = A->n;
  nrhs = 1;
  lda = A->m;
  ldb = A->m;
  pnl_vect_clone (x, b);
  /* Beware that Fortran uses a column wise store, we actually consider A^T */
  C2F(dtrtrs)("U","T","N",&n,&nrhs,A->array,&lda,x->array,&ldb,&info);
  if (info != 0)
    {
      PNL_MESSAGE_ERROR ("Matrix is singular", "pnl_mat_lower_syslin");
      return FAIL;
    }
  return OK;
}

/**
 * solves a symmetric DEFINITE POSITIVE linear system using the Cholesky
 * decomposition of the system A x = b
 *
 * @param chol the Cholesky decomposition of the system as computed by pnl_mat_chol
 * @param b right hand side member. On exit, b contains the solution of the system.
 * @return OK or FAIL
 */
int pnl_mat_chol_syslin_inplace (const PnlMat *chol, PnlVect *b)
{
  int n, nrhs, lda, ldb, info;
  CheckIsSquare(chol);
  CheckMatVectIsCompatible (chol, b);
  n = chol->m;
  lda = chol->m;
  ldb = b->size;
  nrhs = 1;

  C2F(dpotrs)("U", &n, &nrhs, chol->array, &lda, b->array, &ldb, &info);
  if (info != 0)
    {
      PNL_MESSAGE_ERROR ("illegal value", "pnl_mat_chol_syslin");
      return FAIL;
    }
  return OK;
}

/**
 * solves a symmetric DEFINITE POSITIVE linear system using the Cholesky
 * decomposition of the system A x = b
 *
 * @param x already existing PnlVect that contains the solution on exit
 * @param chol the Cholesky decomposition of the system as computed by pnl_mat_chol
 * @param b right hand side member
 * @return OK or FAIL
 */
int pnl_mat_chol_syslin (PnlVect *x, const PnlMat *chol, const  PnlVect *b)
{
  pnl_vect_clone (x,b);
  return pnl_mat_chol_syslin_inplace (chol, x);
}

/**
 * solves a linear system A X = B using a Cholesky factorization of A where B
 * is a matrix. Note that A must be symmetrix positive definite
 *
 * @param A contains the Cholesky decomposition of the matrix A as computed by
 * pnl_mat_chol
 * @param B the r.h.s. matrix of the system of size n x m. On exit B contains
 * the solution X
 * @return OK or FAIL
 */
int pnl_mat_chol_syslin_mat (const PnlMat *A,  PnlMat *B)
{
  int n, nrhs, lda, ldb, info;
  PnlMat *tB;
  CheckIsSquare(A);
  PNL_CHECK (A->m != B->m, "size mismatch", "pnl_mat_chol_syslin_mat");
  n = A->m;
  lda = A->m;
  ldb = B->m;
  nrhs = B->n;

  /* Some tweaks are needed because of Fortran storage in Blas */
  tB = pnl_mat_create (0,0);
  /* Convert to column wise storage */
  pnl_mat_tr (tB, B);
  C2F(dpotrs)("U", &n, &nrhs, A->array, &lda, tB->array, &ldb, &info);
  if (info != 0)
    {
      PNL_MESSAGE_ERROR ("illegal value", "pnl_mat_chol_syslin");
      return FAIL;
    }
  /* Revert to row wise storage */
  pnl_mat_tr (B, tB);
  pnl_mat_free (&tB);
  return OK;
}

/**
 * solves the linear system A x = b with P A = LU. For a symmetric definite
 * positive system, prefer pnl_mat_chol_syslin
 *
 * @param A a PnlMat containing the LU decomposition of A
 * @param p a PnlVectInt.
 * @param b right hand side member. Contains the solution x on exit
 * @return OK or FAIL
 */
int pnl_mat_lu_syslin_inplace (PnlMat *A, const PnlVectInt *p, PnlVect *b)
{
  int i, n, nrhs, lda, ldb, info;
  CheckIsSquare(A);
  CheckMatVectIsCompatible (A, b);
  CheckVectMatch (p, b); 

  n = A->n;
  nrhs = 1;
  lda = A->m;
  ldb = A->m;
  /* Fortran indices start at 1 */
  for ( i=0 ; i<p->size ; i++ ) (p->array[i])++;
  /* Convert to a column wise storage */
  pnl_mat_sq_transpose (A);
  C2F(dgetrs)("N",&n,&nrhs,A->array,&lda,p->array,b->array,&ldb,&info);
  /* Revert to a row wise storage */
  pnl_mat_sq_transpose (A);
  /* Fortran indices start at 1 */
  for ( i=0 ; i<p->size ; i++ ) (p->array[i])--;
  if (info != 0)
    {
      PNL_MESSAGE_ERROR ("Matrix is singular", "pnl_lu_syslin");
      return FAIL;
    }
  return OK;
}

/**
 * solves the linear system A x = b with P A = LU. For a symmetric definite
 * positive system, prefer pnl_mat_chol_syslin
 *
 * @param x existing vector that contains the solution on exit
 * @param LU a PnlMat containing the LU decomposition of A
 * @param b right hand side member
 * @param p a PnlVectInt.
 * @return OK or FAIL
 */
int pnl_mat_lu_syslin (PnlVect *x, PnlMat *LU, const PnlVectInt *p, const PnlVect *b)
{
  pnl_vect_clone (x, b);
  return pnl_mat_lu_syslin_inplace (LU, p, x);
}

/**
 * solves the linear system A x = b with A P = QR.  *
 * @param x a PnlVect containing the solution on exit
 * @param Q a orthogonal PnlMat 
 * @param U an upper triagular PnlMat 
 * @param p a PnlVectInt (permutation vector)
 * @param b right hand side member
 * @return OK or FAIL
 */
int pnl_mat_qr_syslin (PnlVect *x, PnlMat *Q, PnlMat *R, const PnlVectInt *p, PnlVect *b)
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
  return OK;
}

/**
 * solves a linear system A x = b using a LU factorization
 * @param x a PnlVect containing the solution on exit (must have already
 * been created )
 * @param A the matrix of the system
 * @param b the r.h.s. member
 * @return OK or FAIL
 */
int pnl_mat_syslin (PnlVect *x, const PnlMat *A, const PnlVect *b)
{
  PnlMat *LU;
  int status;
  LU = pnl_mat_copy (A);
  pnl_vect_clone (x, b);
  status = pnl_mat_syslin_inplace (LU, x);
  pnl_mat_free (&LU);
  return status;
}

/**
 * solves a linear system A x = b using a LU factorization
 * @param A the matrix of the system. On exit contains the LU decomposition of A.
 * @param b the r.h.s. member
 * @return OK or FAIL
 */
int pnl_mat_syslin_inplace (PnlMat *A, PnlVect *b)
{
  PnlVectInt *p;
  int status;
  CheckIsSquare(A);
  p = pnl_vect_int_create (A->m);
  status = pnl_mat_lu (A, p);
  if ( status != OK ) return FAIL;
  status = pnl_mat_lu_syslin_inplace (A, p, b);
  pnl_vect_int_free (&p);
  return status;
}

/**
 * solves a linear system A X = B using a LU factorization where B is a matrix
 * @param A contains the L and U factors of the PA = LU factoratisation
 * previously computed by pnl_mat_lu
 * @param p the permutation associated to the PA = LU factotisation
 * @param B the r.h.s. matrix of the system of size n x m. On exit B contains
 * the solution X
 * @return OK or FAIL
 */
int pnl_mat_lu_syslin_mat (const PnlMat *A,  const PnlPermutation *p, PnlMat *B)
{
  int i, n, nrhs, lda, ldb, info;
  PnlMat *tB;
  CheckIsSquare(A);
  PNL_CHECK (A->m != B->m, "size mismatch", "pnl_mat_lu_syslin_mat");

  n = A->n;
  nrhs = B->n;
  lda = A->m;
  ldb = A->m;

  /* Some tweaks are needed because of Fortran storage in Blas */
  tB = pnl_mat_create (0,0);
  /* Convert to column wise storage */
  pnl_mat_sq_transpose ( (PnlMat *) A); /* drop the const because A is reverted at the end */
  pnl_mat_tr (tB, B);
  /* shift indices of pivoting */
  for ( i=0 ; i<p->size ; i++ ) (p->array[i])++;

  C2F(dgetrs)("N",&n,&nrhs,A->array,&lda,p->array,tB->array,&ldb,&info);
  /* Revert to row wise storage */
  pnl_mat_sq_transpose ( (PnlMat *) A); /* drop the const because A is reverted at the end */

  pnl_mat_tr (B, tB);
  /* unshift indices of pivoting */
  for ( i=0 ; i<p->size ; i++ ) (p->array[i])--;
  if (info != 0)
    {
      PNL_MESSAGE_ERROR ("Matrix is singular", "pnl_lu_syslin");
      return FAIL;
    }
  pnl_mat_free (&tB);
  return OK;
}

/**
 * solves a linear system A X = B using a LU factorization where B is a matrix
 * @param A the matrix of the system of size n x n. On exit contains the LU decomposition
 * @param B the r.h.s. matrix of the system of size n x m. On exit B contains
 * the solution X
 * @return OK or FAIL
 */
int pnl_mat_syslin_mat (PnlMat *A,  PnlMat *B)
{
  PnlVectInt *p;
  int status;
  CheckIsSquare(A);
  p = pnl_vect_int_create (A->m);
  pnl_mat_lu (A, p);

  status = pnl_mat_lu_syslin_mat (A, p , B);
  pnl_vect_int_free (&p);
  return status;
}

/**
 * inversion of an upper triangular matrix
 *
 * @param A on exit, contains the inverse of B. A must be an already allocated PnlMat
 * @param B an upper triangular matrix
 * @return OK or FAIL
 */
int pnl_mat_upper_inverse(PnlMat *A, const PnlMat *B)
{
  int n, lda, info;
  
  CheckIsSquare(B);
  pnl_mat_clone (A, B);
  n = A->n;
  lda = A->m;
  /* Beware that Fortran uses a column wise storage, we actually consider A^T */
  C2F(dtrtri)("L","N",&n,A->array,&lda,&info);
  if (info != 0)
    {
      PNL_MESSAGE_ERROR ("Matrix is singular", "pnl_mat_upper_inverse");
      return FAIL;
    }
  return OK;
}

/**
 * inversion of a lower triangular matrix
 *
 * @param A on exit, contains the inverse of B. A must be an already allocated PnlMat
 * @param B a lower triangular matrix
 * @return OK or FAIL
 */
int pnl_mat_lower_inverse (PnlMat *A, const PnlMat *B)
{
  int n, lda, info;
  
  CheckIsSquare(B);
  pnl_mat_clone (A, B);
  n = A->n;
  lda = A->m;
  /* Beware that Fortran uses a column wise storage, we actually consider A^T */
  C2F(dtrtri)("U","N",&n,A->array,&lda,&info);
  if (info != 0)
    {
      PNL_MESSAGE_ERROR ("Matrix is singular", "pnl_mat_upper_inverse");
      return FAIL;
    }
  return OK;
}

/**
 * Computes the inverse of a symmetric positive defnite matrix using a Cholesky
 * decomposition
 *
 * @param A a matrix.
 * @param inv a PnlMat (already allocated). contains A^-1 on exit.
 * @return OK or FAIL
 */
int pnl_mat_inverse_with_chol (PnlMat *inv, const PnlMat *A)
{
  int i, j, n, lda, info;
  pnl_mat_clone (inv, A);
  pnl_mat_chol (inv);
  n = A->m;
  lda = A->m;

  C2F(dpotri)("U", &n, inv->array, &lda, &info);
  if (info != 0)
    {
      PNL_MESSAGE_ERROR ("illegal values", "pnl_mat_inverse_with_chol");
      return FAIL;
    }
  /* Now we need to symmetrise inv because the upper part is 0 */
  for ( i=0 ; i<inv->m ; i++ )
    {
      for ( j=0 ; j<i ; j++ )
        {
          PNL_MLET (inv, j, i) = PNL_MGET (inv, i, j);
        }
    }
  return OK;
}

/**
 * Computes the inverse of a matrix using a LU decomposition
 *
 * @param A a matrix.
 * @param inv a PnlMat (already allocated). contains
 * \verbatim A^-1 \endverbatim on exit.
 * @return OK or FAIL
 */
int pnl_mat_inverse (PnlMat *inv, const PnlMat *A)
{
  int n, lda, lwork, info, *ipiv;
  double *work, qwork;
  n = A->m;
  lda = A->m;
  pnl_mat_tr (inv, A);
  ipiv = MALLOC_INT (A->m);
  C2F(dgetrf)(&n, &n, inv->array, &lda, ipiv, &info);
  if (info != 0)
    {
      free (ipiv); 
      PNL_MESSAGE_ERROR ("matrix is singular", "pnl_mat_inverse");
      return FAIL;
    }
  lwork = -1;
  C2F(dgetri)(&n,inv->array,&lda,ipiv,&qwork,&lwork,&info);
  if (info != 0)
    {
      PNL_MESSAGE_ERROR ("Cannot query workspace", "pnl_mat_inverse");
      return FAIL;
    }
  lwork = (int) qwork;
  work = MALLOC_DOUBLE (lwork);
  C2F(dgetri)(&n,inv->array,&lda,ipiv,work,&lwork,&info);
  if (info != 0)
    {
      free (ipiv); free (work);
      PNL_MESSAGE_ERROR ("matrix is singular", "pnl_mat_inverse");
      return FAIL;
    }
  pnl_mat_sq_transpose (inv);
  free (ipiv); free (work);
  return OK;
}

/**
 * Computes the eigenvalues and eigenvectors of a real matrix
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
 * Wrapper to dgeev
 * Computes the eigenvalues and eigenvectors of a real symmetric matrix
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

  /* Revert to row wise storage */
  pnl_mat_sq_transpose (P);

  free(work);
  return OK;

 err:
  free(work);  
  return FAIL;
}

/**
 * Wrapper to dgeev
 * Computes the eigenvalues and eigenvectors of a real non symmetric matrix
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

  /* Convert to row wise storage */
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
      if (pnl_vect_get(D, i) <= 0)
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
          pnl_mat_set ( P, i, j, pnl_mat_get (P, i, j) * log (pnl_vect_get (D, j)) );
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


