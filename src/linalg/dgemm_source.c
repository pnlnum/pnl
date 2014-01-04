
/************************************************************************/
/* Copyright Jérôme Lelong <jerome.lelong@gmail.com>                    */
/* Copyright David Pommier <pommier.david@gmail.com                     */
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

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_internals.h"
#define BASE_DOUBLE
#include "pnl/pnl_templates_on.h"
#endif


/**
 * Compute C := alpha * A * B + C with alpha != 0
 *
 * @param alpha a real coefficient
 * @param A a matrix
 * @param B a matrix
 * @param C a matrix. On exit contains the result of alpha * A * B + C
 */
static void FUNCTION(pnl_mat,dgemmNN) (BASE alpha, const TYPE(PnlMat) *A, const TYPE(PnlMat) *B, TYPE(PnlMat) *C)
{
  int i, i1, i2, j, k, p, block;
  BASE Cij, Cij1, Cij2, Bkj, aibk, ai1bk, ai2bk;
  
  block = 3;

#ifndef PNL_RANGE_CHECK_OFF
  if (A->n != B->m || (A->m != C->m || B->n != C->n) )
    {
      PNL_ERROR ("size mismatch", "pnl_mat_dgemm");
    }
#endif

  p = C->m % block;
  for (i=0; i<p; i++)
    {
      for (k=0; k<A->n; k++)
        {
          BASE temp = MULT(alpha, PNL_MGET (A, i, k));
          for (j=0; j<C->n; j++)
            {
              Cij = PNL_MGET(C, i, j);
              PLUSEQ (Cij, MULT(temp, PNL_MGET (B, k, j)));
              PNL_MLET (C, i, j) = Cij;
            }
        }
    }
  
  for (i=p, i1=p+1, i2=p+2; i<C->m; i+=block, i1+=block, i2+=block)
    {
      for (k=0; k<A->n; k++) 
        {
          BASE temp = MULT(alpha, PNL_MGET (A, i, k));
          BASE temp1 = MULT(alpha, PNL_MGET (A, i1, k));
          BASE temp2 = MULT(alpha, PNL_MGET (A, i2, k));
          for (j=0; j<C->n; j++)
            {
              Bkj = PNL_MGET (B, k, j);

              Cij = PNL_MGET(C, i, j);
              aibk = MULT(temp, Bkj);
              PLUSEQ (Cij, aibk);
              PNL_MLET (C, i, j) = Cij;

              Cij1 = PNL_MGET(C, i1, j);
              ai1bk = MULT(temp1, Bkj);
              PLUSEQ (Cij1, ai1bk);
              PNL_MLET (C, i1, j) = Cij1;
              
              Cij2 = PNL_MGET(C, i2, j);
              ai2bk = MULT(temp2, Bkj);
              PLUSEQ (Cij2, ai2bk);
              PNL_MLET (C, i2, j) = Cij2;
            }
        }
    }
}

/**
 * Compute C := alpha * A' * B + C with alpha != 0
 *
 * @param alpha a real coefficient
 * @param A a matrix
 * @param B a matrix
 * @param C a matrix. On exit contains the result of alpha * A' * B + C
 */
static void FUNCTION(pnl_mat,dgemmTN) (BASE alpha, const TYPE(PnlMat) *A, const TYPE(PnlMat) *B, TYPE(PnlMat) *C)
{
  TYPE(PnlMat) * tA;

  /* It seems much faster to create a new matrix which is the transposition of
     A and call dgemmNN, rather than trying to rewrite loops in order to make
     the most of data alignment.

     Atcually, to avoid using all the memory, we should break A, B and C into
     smaller blocks and only copy these smaller blocks of A.
  */
  tA = FUNCTION(pnl_mat,transpose) (A);
  FUNCTION(pnl_mat,dgemmNN) (alpha, tA, B, C);
  FUNCTION(pnl_mat,free) (&tA);
  return;
}

/**
 * Compute C := alpha * A * B' + C with alpha != 0
 *
 * @param alpha a real coefficient
 * @param A a matrix
 * @param B a matrix
 * @param C a matrix. On exit contains the result of alpha * A * B' + C
 */
static void FUNCTION(pnl_mat,dgemmNT) (BASE alpha, const TYPE(PnlMat) *A, const TYPE(PnlMat) *B, TYPE(PnlMat) *C)
{
  int i, i1, i2, j, k, p, block;
  BASE Cij, Cij1, Cij2, Bkj, sum, sum1, sum2;
  

  block = 3;

#ifndef PNL_RANGE_CHECK_OFF
  if (A->n != B->n || (A->m != C->m || B->m != C->n) )
    {
      PNL_ERROR ("size mismatch", "pnl_mat_dgemm");
    }
#endif

  p = C->m % block;
  for (i=0; i<p; i++)
    {
      for (j=0; j<C->n; j++)
        {
          sum = ZERO;
          for (k=0; k<A->n; k++)
            {
              PLUSEQ (sum, MULT(PNL_MGET (A, i, k),
                                PNL_MGET (B, j, k) ) );
            }
          Cij = PNL_MGET(C, i, j);
          PLUSEQ (Cij, MULT(alpha, sum));
          PNL_MLET(C, i, j) =  Cij;
        }
    }
  
 
  for (i=p, i1=p+1, i2=p+2; i<C->m; i+=block, i1+=block, i2+=block)
    {
      for (j=0; j<C->n; j++)
        {
          sum = ZERO;
          sum1 = ZERO;
          sum2 = ZERO;
          for (k=0; k<A->n; k++) 
            {
              Bkj = PNL_MGET (B, j, k);
   
              PLUSEQ (sum, MULT( PNL_MGET(A, i, k), Bkj));
              PLUSEQ (sum1, MULT( PNL_MGET(A, i1, k), Bkj));
              PLUSEQ (sum2, MULT( PNL_MGET(A, i2, k), Bkj));
   
            }
          Cij = PNL_MGET(C, i, j);
          MULTEQ (sum, alpha);
          PLUSEQ (Cij, sum);
          PNL_MLET (C, i, j) = Cij;
          
          Cij1 = PNL_MGET(C, i1, j);
          MULTEQ (sum1, alpha);
          PLUSEQ (Cij1, sum1);
          PNL_MLET (C, i1, j) = Cij1;
           
          Cij2 = PNL_MGET(C, i2, j);
          MULTEQ (sum2, alpha);
          PLUSEQ (Cij2, sum2);
          PNL_MLET (C, i2, j) = Cij2;
        }
    }
}

/**
 * Compute C := alpha * A' * B' + C with alpha != 0
 *
 * @param alpha a real coefficient
 * @param A a matrix
 * @param B a matrix
 * @param C a matrix. On exit contains the result of alpha * A' * B' + C
 */
static void FUNCTION(pnl_mat,dgemmTT) (BASE alpha, const TYPE(PnlMat) *A, const TYPE(PnlMat) *B, TYPE(PnlMat) *C)
{
  TYPE(PnlMat) *tA;
  /* It seems much faster to create a new matrix which is the transposition of
     A and call dgemmNT, rather than trying to rewrite loops in order to make
     the most of data alignment.
     
     Atcually, to avoid using all the memory, we should break A, B and C into
     smaller blocks and only copy this smaller block of A.
  */
  tA = FUNCTION(pnl_mat,transpose) (A);
  FUNCTION(pnl_mat,dgemmNT) (alpha, tA, B, C);
  FUNCTION(pnl_mat,free) (&tA);
}


/**
 * Compute C := alpha * op (A) * op (B) + beta C
 *
 * @param transA a char :
 * 'N' or 'n' for op (A) = A  and 'T' or 't' for op (A) = A'
 * @param transB a char :
 * 'N' or 'n' for op (B) = B  and 'T' or 't' for op (B) = B'
 * @param alpha a real coefficient
 * @param A a matrix
 * @param B a matrix
 * @param beta a real coefficient
 * @param C a matrix. On exit contains the result of alpha * op (A) * op (B) + beta
 * C. When beta equals 0, the size of C is not checked and its content not
 * used; instead C is resized to store the result
 */
void FUNCTION(pnl_mat,dgemm) (char transA, char transB, BASE alpha, const TYPE(PnlMat) *A,
                              const TYPE(PnlMat) *B, BASE beta, TYPE(PnlMat) *C)
{
  int m, n, opA, opB;

  m = A->m; n = B->n;
  opA = 0; opB = 0;
  if (transA == 'T' || transA == 't') { m = A->n; opA = 1; }
  if (transB == 'T' || transB == 't') { n = B->m; opB = 1; }
  
  if ( EQ(beta,ZERO) )
    {
      FUNCTION(pnl_mat,resize) (C, m, n);
      FUNCTION(pnl_mat,set_zero)(C);
    }
  else if ( NEQ(beta, ONE) )
    {
      FUNCTION(pnl_mat,mult_scalar) (C, beta); 
    }
  if ( EQ(alpha,ZERO) ) return;

#ifndef PNL_RANGE_CHECK_OFF
  if (((opA==0 && opB == 0) && (A->n != B->m || (A->m != C->m || B->n != C->n)))
      || ((opA==1 && opB == 0) && (A->m != B->m || (A->n != C->m || B->n != C->n)))
      || ((opA==1 && opB == 1) && (A->m != B->n || (A->n != C->m || B->m != C->n)))
      || ((opA==0 && opB == 1) && (A->n != B->n || (A->m != C->m || B->m != C->n))))    
    {
      PNL_ERROR ("size mismatch", "pnl_mat_dgemm");
    }
#endif

  /* Form alpha * A * B + beta * C */
  if (opA == 0 && opB == 0)
    {
      FUNCTION(pnl_mat,dgemmNN) (alpha, A, B, C);
    }
  else if (opA == 1 && opB == 0)
    /* Form alpha * A' * B + beta * C */
    {
      FUNCTION(pnl_mat,dgemmTN) (alpha, A, B, C);
    }
  else if (opA == 1 && opB == 1)
    /* Form alpha * A' * B' + beta * C */
    {
      FUNCTION(pnl_mat,dgemmTT) (alpha, A, B, C);
    }
  else
    /* Form alpha * A * B' + beta * C */
    {
      FUNCTION(pnl_mat,dgemmNT) (alpha, A, B, C);
    }
}

#ifdef PNL_CLANG_COMPLETE
#include "pnl/pnl_templates_off.h"
#undef  BASE_DOUBLE
#endif


