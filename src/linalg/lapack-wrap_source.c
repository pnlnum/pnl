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

/* To enable clang completion */
#ifdef PNL_CLANG_COMPLETE
#include <stdlib.h>
#include <string.h>

#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_machine.h"
#include "pnl/clapack.h"
#define BASE_DOUBLE
#include "pnl/pnl_templates_on.h"
#endif



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
int FUNCTION(pnl_mat,chol) (TYPE(PnlMat) *M)
{
  int n, lda, info, i, j;
  CheckIsSquare(M);
  n = M->m;
  lda = M->m;
  /* Because of Fortran column wise storage, we ask for an upper triangular
   * matrix to actually have a lower one */
  PNL_C2F(potrf)("U", &n, M->array, &lda, &info);
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
          PNL_MLET (M, i, j) = ZERO;
        }
    }
  return OK;
}

/** 
 * Compute the Cholesky decompisition with complete pivoting. 
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
int FUNCTION(pnl_mat,pchol) (TYPE(PnlMat) *M, double tol, int *rank, PnlVectInt *p)
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
  PNL_C2F(pstrf)("U", &n, M->array, &lda, p->array, rank, &tol, work, &info);
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
          PNL_MLET (M, i, j) = ZERO;
        }
    }
  /* C indices start at 0 */
  pnl_vect_int_minus_scalar (p, 1);
  return OK;
}

/**
 * Compute a P A = LU factorisation. On exit A contains the L and U
 * matrices. Note that the diagonal elements of L are all 1.
 *
 * @param A the matrix to decompose.
 * @param p a PnlPermutation.
 * @return OK or FAIL
 */
int FUNCTION(pnl_mat,lu) (TYPE(PnlMat) *A, PnlPermutation *p)
{
  int info, N = A->n;
  
  CheckIsSquare(A);

  FUNCTION(pnl_mat,sq_transpose) (A);
  PNL_C2F(getrf) (&N, &N, A->array, &N, p->array, &info);
  if ( info != 0 )
    {
      PNL_MESSAGE_ERROR ("LU decomposition cannot be computed", "pnl_mat_lu");
      return FAIL;
    }
  FUNCTION(pnl_mat,sq_transpose) (A);
  /* C indices start at 0 */
  pnl_vect_int_minus_scalar (p, 1);
  return OK;
}

/**
 * Solve an upper triangular linear system
 *
 * @param x already existing PnlVect that contains the solution on exit
 * @param A an upper triangular matric
 * @param b right hand side member
 * @return OK or FAIL
 */
int FUNCTION(pnl_mat,upper_syslin) (TYPE(PnlVect) *x, const TYPE(PnlMat) *A, const  TYPE(PnlVect) *b)
{
  int n, nrhs, lda, ldb, info;
  
  CheckIsSquare(A);
  n = A->n;
  nrhs = 1;
  lda = A->m;
  ldb = A->m;
  FUNCTION(pnl_vect,clone) (x, b);
  /* Beware that Fortran uses a column wise store, we actually consider A^T */
  PNL_C2F(trtrs)("L","T","N",&n,&nrhs,A->array,&lda,x->array,&ldb,&info);
  if (info != 0)
    {
      PNL_MESSAGE_ERROR ("Matrix is singular", "pnl_mat_upper_syslin");
      return FAIL;
    }
  return OK;
}

/**
 * Solve a lower triangular linear system
 *
 * @param x already existing PnlVect that contains the solution on exit
 * @param A a lower triangular matrix
 * @param b right hand side member
 * @return OK or FAIL
 */
int FUNCTION(pnl_mat,lower_syslin) (TYPE(PnlVect) *x, const TYPE(PnlMat) *A, const  TYPE(PnlVect) *b)
{
  int n, nrhs, lda, ldb, info;
  
  CheckIsSquare(A);
  n = A->n;
  nrhs = 1;
  lda = A->m;
  ldb = A->m;
  FUNCTION(pnl_vect,clone) (x, b);
  /* Beware that Fortran uses a column wise store, we actually consider A^T */
  PNL_C2F(trtrs)("U","T","N",&n,&nrhs,A->array,&lda,x->array,&ldb,&info);
  if (info != 0)
    {
      PNL_MESSAGE_ERROR ("Matrix is singular", "pnl_mat_lower_syslin");
      return FAIL;
    }
  return OK;
}

/**
 * Solve a symmetric DEFINITE POSITIVE linear system using the Cholesky
 * decomposition of the system A x = b
 *
 * @param chol the Cholesky decomposition of the system as computed by pnl_mat_chol
 * @param b right hand side member. On exit, b contains the solution of the system.
 * @return OK or FAIL
 */
int FUNCTION(pnl_mat,chol_syslin_inplace) (const TYPE(PnlMat) *chol, TYPE(PnlVect) *b)
{
  int n, nrhs, lda, ldb, info;
  CheckIsSquare(chol);
  CheckMatVectIsCompatible (chol, b);
  n = chol->m;
  lda = chol->m;
  ldb = b->size;
  nrhs = 1;

  PNL_C2F(potrs)("U", &n, &nrhs, chol->array, &lda, b->array, &ldb, &info);
  if (info != 0)
    {
      PNL_MESSAGE_ERROR ("illegal value", "pnl_mat_chol_syslin");
      return FAIL;
    }
  return OK;
}

/**
 * Solve a symmetric DEFINITE POSITIVE linear system using the Cholesky
 * decomposition of the system A x = b
 *
 * @param x already existing TYPE(PnlVect) that contains the solution on exit
 * @param chol the Cholesky decomposition of the system as computed by pnl_mat_chol
 * @param b right hand side member
 * @return OK or FAIL
 */
int FUNCTION(pnl_mat,chol_syslin) (TYPE(PnlVect) *x, const TYPE(PnlMat) *chol, const  TYPE(PnlVect) *b)
{
  FUNCTION(pnl_vect,clone) (x,b);
  return FUNCTION(pnl_mat,chol_syslin_inplace) (chol, x);
}

/**
 * Solve a linear system A X = B using a Cholesky factorization of A where B
 * is a matrix. Note that A must be symmetrix positive definite
 *
 * @param A contains the Cholesky decomposition of the matrix A as computed by
 * pnl_mat_chol
 * @param B the r.h.s. matrix of the system of size n x m. On exit B contains
 * the solution X
 * @return OK or FAIL
 */
int FUNCTION(pnl_mat,chol_syslin_mat) (const TYPE(PnlMat) *A,  TYPE(PnlMat) *B)
{
  int n, nrhs, lda, ldb, info;
  TYPE(PnlMat) *tB;
  CheckIsSquare(A);
  PNL_CHECK (A->m != B->m, "size mismatch", "pnl_mat_chol_syslin_mat");
  n = A->m;
  lda = A->m;
  ldb = B->m;
  nrhs = B->n;

  /* Some tweaks are needed because of Fortran storage in Blas */
  tB = FUNCTION(pnl_mat,create) (0,0);
  /* Convert to column wise storage */
  FUNCTION(pnl_mat,tr) (tB, B);
  PNL_C2F(potrs)("U", &n, &nrhs, A->array, &lda, tB->array, &ldb, &info);
  if (info != 0)
    {
      PNL_MESSAGE_ERROR ("illegal value", "pnl_mat_chol_syslin");
      return FAIL;
    }
  /* Revert to row wise storage */
  FUNCTION(pnl_mat,tr) (B, tB);
  FUNCTION(pnl_mat,free) (&tB);
  return OK;
}

/**
 * Solve the linear system A x = b with P A = LU. For a symmetric definite
 * positive system, prefer pnl_mat_chol_syslin
 *
 * @param A a TYPE(PnlMat) containing the LU decomposition of A
 * @param p a TYPE(PnlVect)Int.
 * @param b right hand side member. Contains the solution x on exit
 * @return OK or FAIL
 */
int FUNCTION(pnl_mat,lu_syslin_inplace) (TYPE(PnlMat) *A, const PnlVectInt *p, TYPE(PnlVect) *b)
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
  FUNCTION(pnl_mat,sq_transpose) (A);
  PNL_C2F(getrs)("N",&n,&nrhs,A->array,&lda,p->array,b->array,&ldb,&info);
  /* Revert to a row wise storage */
  FUNCTION(pnl_mat,sq_transpose) (A);
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
 * Solve the linear system A x = b with P A = LU. For a symmetric definite
 * positive system, prefer pnl_mat_chol_syslin
 *
 * @param x existing vector that contains the solution on exit
 * @param LU a TYPE(PnlMat) containing the LU decomposition of A
 * @param b right hand side member
 * @param p a TYPE(PnlVect)Int.
 * @return OK or FAIL
 */
int FUNCTION(pnl_mat,lu_syslin) (TYPE(PnlVect) *x, TYPE(PnlMat) *LU, const PnlVectInt *p, const TYPE(PnlVect) *b)
{
  FUNCTION(pnl_vect,clone) (x, b);
  return FUNCTION(pnl_mat,lu_syslin_inplace) (LU, p, x);
}


/**
 * Solve a linear system A x = b using a LU factorization
 * @param A the matrix of the system. On exit contains the LU decomposition of A.
 * @param b the r.h.s. member
 * @return OK or FAIL
 */
int FUNCTION(pnl_mat,syslin_inplace) (TYPE(PnlMat) *A, TYPE(PnlVect) *b)
{
  PnlVectInt *p;
  int status;
  CheckIsSquare(A);
  p = pnl_vect_int_create (A->m);
  status = FUNCTION(pnl_mat,lu) (A, p);
  if ( status != OK ) return FAIL;
  status = FUNCTION(pnl_mat,lu_syslin_inplace) (A, p, b);
  pnl_vect_int_free (&p);
  return status;
}

/**
 * Solve a linear system A x = b using a LU factorization
 * @param x a TYPE(PnlVect) containing the solution on exit (must have already
 * been created )
 * @param A the matrix of the system
 * @param b the r.h.s. member
 * @return OK or FAIL
 */
int FUNCTION(pnl_mat,syslin) (TYPE(PnlVect) *x, const TYPE(PnlMat) *A, const TYPE(PnlVect) *b)
{
  TYPE(PnlMat) *LU;
  int status;
  LU = FUNCTION(pnl_mat,copy) (A);
  FUNCTION(pnl_vect,clone) (x, b);
  status = FUNCTION(pnl_mat,syslin_inplace) (LU, x);
  FUNCTION(pnl_mat,free) (&LU);
  return status;
}

/**
 * Solve a linear system A X = B using a LU factorization where B is a matrix
 * @param A contains the L and U factors of the PA = LU factoratisation
 * previously computed by pnl_mat_lu
 * @param p the permutation associated to the PA = LU factotisation
 * @param B the r.h.s. matrix of the system of size n x m. On exit B contains
 * the solution X
 * @return OK or FAIL
 */
int FUNCTION(pnl_mat,lu_syslin_mat) (const TYPE(PnlMat) *A,  const PnlPermutation *p, TYPE(PnlMat) *B)
{
  int i, n, nrhs, lda, ldb, info;
  TYPE(PnlMat) *tB;
  CheckIsSquare(A);
  PNL_CHECK (A->m != B->m, "size mismatch", "pnl_mat_lu_syslin_mat");

  n = A->n;
  nrhs = B->n;
  lda = A->m;
  ldb = A->m;

  /* Some tweaks are needed because of Fortran storage in Blas */
  tB = FUNCTION(pnl_mat,create) (0,0);
  /* Convert to column wise storage */
  FUNCTION(pnl_mat,sq_transpose) ( (TYPE(PnlMat) *) A); /* drop the const because A is reverted at the end */
  FUNCTION(pnl_mat,tr) (tB, B);
  /* shift indices of pivoting */
  for ( i=0 ; i<p->size ; i++ ) (p->array[i])++;

  PNL_C2F(getrs)("N",&n,&nrhs,A->array,&lda,p->array,tB->array,&ldb,&info);
  /* Revert to row wise storage */
  FUNCTION(pnl_mat,sq_transpose) ( (TYPE(PnlMat) *) A); /* drop the const because A is reverted at the end */

  FUNCTION(pnl_mat,tr) (B, tB);
  /* unshift indices of pivoting */
  for ( i=0 ; i<p->size ; i++ ) (p->array[i])--;
  if (info != 0)
    {
      PNL_MESSAGE_ERROR ("Matrix is singular", "pnl_lu_syslin");
      return FAIL;
    }
  FUNCTION(pnl_mat,free) (&tB);
  return OK;
}

/**
 * Solve a linear system A X = B using a LU factorization where B is a matrix
 * @param A the matrix of the system of size n x n. On exit contains the LU decomposition
 * @param B the r.h.s. matrix of the system of size n x m. On exit B contains
 * the solution X
 * @return OK or FAIL
 */
int FUNCTION(pnl_mat,syslin_mat) (TYPE(PnlMat) *A,  TYPE(PnlMat) *B)
{
  PnlVectInt *p;
  int status;
  CheckIsSquare(A);
  p = pnl_vect_int_create (A->m);
  FUNCTION(pnl_mat,lu) (A, p);

  status = FUNCTION(pnl_mat,lu_syslin_mat) (A, p , B);
  pnl_vect_int_free (&p);
  return status;
}

/**
 * inversion of an upper triangular matrix
 *
 * @param A on exit, contains the inverse of B. A must be an already allocated TYPE(PnlMat)
 * @param B an upper triangular matrix
 * @return OK or FAIL
 */
int FUNCTION(pnl_mat,upper_inverse) (TYPE(PnlMat) *A, const TYPE(PnlMat) *B)
{
  int n, lda, info;
  
  CheckIsSquare(B);
  FUNCTION(pnl_mat,clone) (A, B);
  n = A->n;
  lda = A->m;
  /* Beware that Fortran uses a column wise storage, we actually consider A^T */
  PNL_C2F(trtri)("L","N",&n,A->array,&lda,&info);
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
 * @param A on exit, contains the inverse of B. A must be an already allocated TYPE(PnlMat)
 * @param B a lower triangular matrix
 * @return OK or FAIL
 */
int FUNCTION(pnl_mat,lower_inverse) (TYPE(PnlMat) *A, const TYPE(PnlMat) *B)
{
  int n, lda, info;
  
  CheckIsSquare(B);
  FUNCTION(pnl_mat,clone) (A, B);
  n = A->n;
  lda = A->m;
  /* Beware that Fortran uses a column wise storage, we actually consider A^T */
  PNL_C2F(trtri)("U","N",&n,A->array,&lda,&info);
  if (info != 0)
    {
      PNL_MESSAGE_ERROR ("Matrix is singular", "pnl_mat_upper_inverse");
      return FAIL;
    }
  return OK;
}

/**
 * Compute the inverse of a symmetric positive defnite matrix using a Cholesky
 * decomposition
 *
 * @param A a matrix.
 * @param inv a TYPE(PnlMat) (already allocated). contains A^-1 on exit.
 * @return OK or FAIL
 */
int FUNCTION(pnl_mat,inverse_with_chol) (TYPE(PnlMat) *inv, const TYPE(PnlMat) *A)
{
  int i, j, n, lda, info;
  FUNCTION(pnl_mat,clone) (inv, A);
  FUNCTION(pnl_mat,chol) (inv);
  n = A->m;
  lda = A->m;

  PNL_C2F(potri)("U", &n, inv->array, &lda, &info);
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
 * Compute the inverse of a matrix using a LU decomposition
 *
 * @param A a matrix.
 * @param inv a TYPE(PnlMat) (already allocated). contains
 * \verbatim A^-1 \endverbatim on exit.
 * @return OK or FAIL
 */
int FUNCTION(pnl_mat,inverse) (TYPE(PnlMat) *inv, const TYPE(PnlMat) *A)
{
  int n, lda, lwork, info, *ipiv;
  BASE *work, qwork;
  n = A->m;
  lda = A->m;
  FUNCTION(pnl_mat,tr) (inv, A);
  ipiv = MALLOC_INT (A->m);
  PNL_C2F(getrf)(&n, &n, inv->array, &lda, ipiv, &info);
  if (info != 0)
    {
      free (ipiv); 
      PNL_MESSAGE_ERROR ("matrix is singular", "pnl_mat_inverse");
      return FAIL;
    }
  lwork = -1;
  PNL_C2F(getri)(&n,inv->array,&lda,ipiv,&qwork,&lwork,&info);
  if (info != 0)
    {
      PNL_MESSAGE_ERROR ("Cannot query workspace", "pnl_mat_inverse");
      return FAIL;
    }
#if MULTIPLICITY == 2
  lwork = (int) Creal(qwork);
#else
  lwork = (int) qwork;
#endif
  work = MALLOC_BASE (lwork);
  PNL_C2F(getri)(&n,inv->array,&lda,ipiv,work,&lwork,&info);
  if (info != 0)
    {
      free (ipiv); free (work);
      PNL_MESSAGE_ERROR ("matrix is singular", "pnl_mat_inverse");
      return FAIL;
    }
  FUNCTION(pnl_mat,sq_transpose) (inv);
  free (ipiv); free (work);
  return OK;
}

#ifdef PNL_CLANG_COMPLETE
#include "pnl/pnl_templates_off.h"
#undef  BASE_DOUBLE
#endif

