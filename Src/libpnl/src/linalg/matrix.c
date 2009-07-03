
/*************************************************************************/
/* Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>            */  
/* and Céline Labart, 2007                                               */
/*                                                                       */
/* This program is free software; you can redistribute it and/or modify  */
/* it under the terms of the GNU General Public License as published by  */
/* the Free Software Foundation; either version 3 of the License, or     */
/* (at your option) any later version.                                   */
/*                                                                       */
/* This program is distributed in the hope that it will be useful,       */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/* GNU General Public License for more details.                          */
/*                                                                       */
/* You should have received a copy of the GNU General Public License     */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/*************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define PNL_MATRIX_PRIVATE
#include "pnl_matrix.h"
#include "pnl_mathtools.h"

/**
 * allocates a PnlVectCompact. By default type='d'.
 * @param n size
 * @param x value to fill the vector
 * @return a pointeur to PnlVectCompact
 */
PnlVectCompact* pnl_vect_compact_create (int n, double x)
{
  PnlVectCompact *v;
  if ((v=malloc (sizeof(PnlVectCompact)))==NULL) return NULL;
  v->size = n;
  v->convert = 'd';
  v->val = x;
  return v;
}

/**
 * resize a PnlVectCompact.
 * @param v the Premiavectcompact to be resized
 * @param size new size
 * @param x new value to set
 * @return OK or WRONG
 */
int pnl_vect_compact_resize (PnlVectCompact *v, int size, double x)
{
  if (v->convert == 'a')
    {
      free (v->array); v->array = NULL;
    }
  v->size = size;
  v->convert = 'd';
  v->val = x;
  return OK;
}

/**
 * copies a PnlVectCompact
 *
 * @param v : a constant PnlVectCompact pointer
 * @return  a PnlVectCompact  pointer initialised with v
 */
PnlVectCompact* pnl_vect_compact_copy(const PnlVectCompact *v)
{
  PnlVectCompact *ret = NULL;

  if (v->convert == 'd')
    ret = pnl_vect_compact_create (v->size, v->val);
  else
    {
      if ((ret=malloc(sizeof(PnlVectCompact)))==NULL) return NULL;
      ret->size = v->size;
      ret->convert = 'a';
      if ((ret->array=malloc(v->size*sizeof(double)))==NULL) return NULL;
      memcpy(ret->array, v->array, sizeof(double)*ret->size);
    }
  return ret;
}


/**
 * free a PnlVectCompact
 * @param v address of a PnlVectCompact
 */
void pnl_vect_compact_free (PnlVectCompact **v)
{
  if ((*v)->convert == 'a')
    {
      free ((*v)->array); (*v)->array = NULL;
    }
  free (*v); *v=NULL;
}



/**
 * Convert a PnlVectCompact pointeur to a PnlVect pointeur
 * @param C the PnlVectCompact to be expanded
 * @return a PnlVect pointeur
 */
PnlVect* pnl_vect_compact_to_pnl_vect (const PnlVectCompact *C)
{
  PnlVect *v;
  if (C->convert == 'd')
    {
      v = pnl_vect_create_from_double (C->size, C->val);
    }
  else
    {
      v = pnl_vect_create (C->size);
      memcpy (v->array, C->array, C->size*sizeof(double));
    }
  return v;
}

/**
 * access function
 * @param C a PnlVectCompact
 * @param i index
 * @return C[i]
 */
double pnl_vect_compact_get (const PnlVectCompact *C, int i)
{
  if (C->convert == 'd') return C->val;
  CheckIndexVect (C, i);
  return C->array[i];
}


/**
 * in-place vector scalar double multiplication
 *
 * @param lhs : left hand side vector
 * @param x : scalar double
 * @return  lhs = lhs*x
 */
void pnl_vect_complex_mult_double(PnlVectComplex *lhs , double x)
{
  fcomplex *lptr;
  int i=0;
  while(i<lhs->size)
    {
      lptr = pnl_vect_complex_lget(lhs, i);
      lptr->r *= x;
      lptr->i *= x;
      i++;
    }
}


/**************************************************
 * PnlMat functions specific to base type double  *
 **************************************************/


/**
 *  matrix multiplication
 *
 * @param rhs1 : first right hand side matrix
 * @param rhs2 : second right hand side matrix
 * @return  rhs1*rhs2
 */

PnlMat* pnl_mat_mult_mat(const PnlMat *rhs1, const PnlMat *rhs2)
{
  PnlMat *lhs;
  lhs = pnl_mat_create (0,0);
  pnl_mat_dgemm ('N', 'N', 1., rhs1, rhs2, 0., lhs);
  return lhs;
}

/**
 *  in-place matrix multiplication
 *
 * @param lhs : left hand side matrix
 * @param rhs1 : first right hand side matrix
 * @param rhs2 : second right hand side matrix
 * @return  lhs=rhs1*rhs2
 */

void pnl_mat_mult_mat_inplace(PnlMat *lhs, const  PnlMat *rhs1, const PnlMat *rhs2)
{
  pnl_mat_dgemm ('N', 'N', 1, rhs1, rhs2, 0., lhs);
}


/**
 * Computes C := alpha * A * ' + C with alpha != 0
 *
 * @param alpha a real coefficient
 * @param A a matrix
 * @param B a matrix
 * @param C a matrix. On exit contains the result of alpha * A * B + C
 */
static void pnl_mat_dgemmNN (double alpha, const PnlMat *A, const PnlMat *B, PnlMat *C)
{
  int i, i1, i2, j, k, p, block;
  double Cij, Cij1, Cij2, Bkj, aibk, ai1bk, ai2bk;
  
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
          const register double temp = alpha * pnl_mat_get (A, i, k);
          for (j=0; j<C->n; j++)
            {
              Cij = pnl_mat_get(C, i, j);
              Cij += temp * pnl_mat_get (B, k, j);
              pnl_mat_set (C, i, j, Cij);
            }
        }
    }
  
  for (i=p, i1=p+1, i2=p+2; i<C->m; i+=block, i1+=block, i2+=block)
    {
      for (k=0; k<A->n; k++) 
        {
          const register double temp = alpha * pnl_mat_get (A, i, k);
          const register double temp1 = alpha * pnl_mat_get (A, i1, k);
          const register double temp2 = alpha * pnl_mat_get (A, i2, k);
          for (j=0; j<C->n; j++)
            {
              Bkj = pnl_mat_get (B, k, j);

              Cij = pnl_mat_get(C, i, j);
              aibk = temp * Bkj;
              Cij += aibk;
              pnl_mat_set (C, i, j, Cij);
              
              Cij1 = pnl_mat_get(C, i1, j);
              ai1bk = temp1 * Bkj;;
              Cij1 += ai1bk;
              pnl_mat_set (C, i1, j, Cij1);
              
              Cij2 = pnl_mat_get(C, i2, j);
              ai2bk = temp2 * Bkj;
              Cij2 += ai2bk;
              pnl_mat_set (C, i2, j, Cij2);
            }
        }
    }
}

/**
 * Computes C := alpha * A' * B + C with alpha != 0
 *
 * @param alpha a real coefficient
 * @param A a matrix
 * @param B a matrix
 * @param C a matrix. On exit contains the result of alpha * A' * B + C
 */
static void pnl_mat_dgemmTN (double alpha, const PnlMat *A, const PnlMat *B, PnlMat *C)
{
  PnlMat * tA;

  /* It seems much faster to create a new matrix which is the transposition of
     A and call dgemmNN, rather than trying to rewrite loops in order to make
     the most of data alignment.

     Atcually, to avoid using all the memory, we should break A, B and C into
     smaller blocks and only copy this smaller block of A.
  */
  tA = pnl_mat_transpose (A);
  pnl_mat_dgemmNN (alpha, tA, B, C);
  pnl_mat_free (&tA);
  return;
}

/**
 * Computes C := alpha * A * B' + C with alpha != 0
 *
 * @param alpha a real coefficient
 * @param A a matrix
 * @param B a matrix
 * @param C a matrix. On exit contains the result of alpha * A * B' + C
 */
static void pnl_mat_dgemmNT (double alpha, const PnlMat *A, const PnlMat *B, PnlMat *C)
{
  int i, i1, i2, j, k, p, block;
  register double Cij, Cij1, Cij2, Bkj, sum, sum1, sum2;
  

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
          sum = 0.;
          for (k=0; k<A->n; k++)
            {
              sum += pnl_mat_get (A, i, k) * pnl_mat_get (B, j, k);
            }
          Cij = pnl_mat_get(C, i, j);
          Cij += alpha * sum;
          pnl_mat_set (C, i, j, Cij);
        }
    }
  
 
  for (i=p, i1=p+1, i2=p+2; i<C->m; i+=block, i1+=block, i2+=block)
    {
      for (j=0; j<C->n; j++)
        {
          sum = 0.;
          sum1 = 0.;
          sum2 = 0.;
          for (k=0; k<A->n; k++) 
            {
              Bkj = pnl_mat_get (B, j, k);
   
              sum += pnl_mat_get (A, i, k) * Bkj;
              sum1 += pnl_mat_get (A, i1, k) * Bkj;
              sum2 += pnl_mat_get (A, i2, k) * Bkj;
   
            }
          Cij = pnl_mat_get(C, i, j);
          sum *= alpha;
          Cij += sum;
          pnl_mat_set (C, i, j, Cij);
          
          Cij1 = pnl_mat_get(C, i1, j);
          sum1 *= alpha;
          Cij1 += sum1;
          pnl_mat_set (C, i1, j, Cij1);
           
          Cij2 = pnl_mat_get(C, i2, j);
          sum2 *= alpha;
          Cij2 += sum2;
          pnl_mat_set (C, i2, j, Cij2);
        }
    }
}

/**
 * Computes C := alpha * A' * B' + C with alpha != 0
 *
 * @param alpha a real coefficient
 * @param A a matrix
 * @param B a matrix
 * @param C a matrix. On exit contains the result of alpha * A' * B' + C
 */
static void pnl_mat_dgemmTT (double alpha, const PnlMat *A, const PnlMat *B, PnlMat *C)
{
  PnlMat *tA;
  /* It seems much faster to create a new matrix which is the transposition of
     A and call dgemmNT, rather than trying to rewrite loops in order to make
     the most of data alignment.
     
     Atcually, to avoid using all the memory, we should break A, B and C into
     smaller blocks and only copy this smaller block of A.
  */
  tA = pnl_mat_transpose (A);
  pnl_mat_dgemmNT (alpha, tA, B, C);
  pnl_mat_free (&tA);
}


/**
 * Computes C := alpha * op (A) * op (B) + beta C
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
void pnl_mat_dgemm (char transA, char transB, double alpha, const PnlMat *A,
                    const PnlMat *B, double beta, PnlMat *C)
{
  int m, n, opA, opB;

  m = A->m; n = B->n;
  opA = 0; opB = 0;
  if (transA == 'T' || transA == 't') { m = A->n; opA = 1; }
  if (transB == 'T' || transB == 't') { n = B->m; opB = 1; }
  
  if (beta == 0.)
    {
      pnl_mat_resize (C, m, n);
      pnl_mat_set_double (C, 0.);
    }
  else if (beta != 1)
    {
      pnl_mat_mult_double (C, beta); 
    }
  if (alpha == 0.) return;

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
      pnl_mat_dgemmNN (alpha, A, B, C);
    }
  else if (opA == 1 && opB == 0)
    /* Form alpha * A' * B + beta * C */
    {
      pnl_mat_dgemmTN (alpha, A, B, C);
    }
  else if (opA == 1 && opB == 1)
    /* Form alpha * A' * B' + beta * C */
    {
      pnl_mat_dgemmTT (alpha, A, B, C);
    }
  else
    /* Form alpha * A * B' + beta * C */
    {
      pnl_mat_dgemmNT (alpha, A, B, C);
    }
}


/**
 * Computes Y := a * X + Y
 *
 * @param a real coefficient
 * @param X a matrix
 * @param Y a matrix whose size must be the same as the one of X. Contains the
 * result on exit
 */
void pnl_mat_axpy (double a, const PnlMat *X, PnlMat *Y)
{
  int i;
  CheckMatMatch(X, Y);
  if (a == 0.) return;
  if (a == 1.)
    {
      for (i=0; i<Y->mn; i++)
        {
          Y->array[i] +=  X->array[i];
        }
      return;
    }
  for (i=0; i<Y->mn; i++)
    {
      Y->array[i] += a * X->array[i];
    }
}

/**
 * Computes A := alpha x' * y + A
 *
 * @param alpha a real number
 * @param x a vector
 * @param y a vector
 * @param A a matrix
 */
void pnl_mat_dger (double alpha, const PnlVect *x, const PnlVect *y, PnlMat *A)
{
  int i, j;
  double ai, alpha_xi;
#ifdef PNL_RANGE_CHECK_OFF
  if ((x->size != y->size) || (A->m != A->n) || (A->m != x->size))
    {
      PNL_ERROR ("size mismatch", "pnl_mat_dger");
    }
#endif
  if (alpha == 0.) return;

  for (i=0; i<x->size; i++)
    {
      alpha_xi = alpha * pnl_vect_get (x, i);
      for (j=0; j<x->size; j++)
        {
          ai = pnl_mat_get (A, i, j);
          ai += alpha_xi * pnl_vect_get (x, j);
          pnl_mat_set (A, i, j, ai);
        }
    }
}

/**
 *  Cmputes A * x
 *
 * @param A : matrix
 * @param x : vector
 * @return  y = A * x
 */
PnlVect* pnl_mat_mult_vect(const PnlMat *A, const PnlVect *x)
{
  PnlVect *y;
  y = pnl_vect_create (0);
  pnl_mat_dgemv ('N', 1., A, x, 0., y);
  return y;
}

/**
 *  Computes A' * x
 *
 * @param A : matrix
 * @param x : vector
 * @return  y = A' * x
 */
PnlVect* pnl_mat_mult_vect_transpose(const PnlMat *A, const PnlVect *x)
{
  PnlVect *y;
  y = pnl_vect_create (A->n);
  pnl_mat_dgemv ('T', 1., A, x, 0., y);
  return y;
}

/**
 *  Computes y = A * x
 *
 * @param A : matrix
 * @param x : vector
 * @param y : vector
 * @return  lhs=mat*rhs
 */
void pnl_mat_mult_vect_inplace(PnlVect *y, const PnlMat *A, const PnlVect *x)
{
  pnl_mat_dgemv ('N', 1., A, x, 0., y);
}


/**
 *  Computes y = A' * x
 *
 * @param A : matrix
 * @param x : vector
 * @param y : vector
 * @return  y = A' * x
 */
void pnl_mat_mult_vect_transpose_inplace(PnlVect *y, const PnlMat *A, const PnlVect *x)
{
  pnl_mat_dgemv ('T', 1., A, x, 0., y);
}

/**
 * Computes y := alpha * A * x + beta * y
 *
 * @param A a  matrix
 * @param x a vector such that A * x is legal
 * @param alpha  a double
 * @param beta  a double. When b==0, the content of y is not used, instead y
 * is resized to store A * x
 * @param y a vector containing alpha * A * x + beta * y on exit.
 */
void pnl_mat_lAxpby(double alpha, const PnlMat *A, const PnlVect *x , double beta, PnlVect * y)
{
  pnl_mat_dgemv ('N', alpha, A, x, beta, y);
}


/**
 * Computes y := alpha * op(A) * x + beta * y
 *     where op(X) = X or op(X) = X'
 *
 * @param A a  matrix
 * @param trans a char. If trans=='T', op(X) = X'. If trans=='N', op(X) = X.
 * @param x a vector such that op(A) * x is legal
 * @param alpha  a double
 * @param beta  a double. When b==0, the content of y is not used, instead y
 * is resized to store op(A) * x
 * @param y a vector containing alpha * op(A) * x + beta * y on exit.
 */
void pnl_mat_dgemv (char trans, double alpha, const PnlMat *A,
                    const PnlVect *x , double beta, PnlVect * y)
{
  double yi, temp;
  int i,j;
  
  if (beta == 0.)
    {
      if (trans=='N' || trans=='n') pnl_vect_resize (y, A->m);
      else pnl_vect_resize (y, A->n);
      pnl_vect_set_double (y, 0.);
    }
  else if (beta != 1) pnl_vect_mult_double (y, beta);
  if (alpha == 0.) return;
  

  /* Form alpha * A * x + beta * y */
  if (trans=='N' || trans=='n')
    {
#ifndef PNL_RANGE_CHECK_OFF
      if (A->n != x->size || A->m != y->size)
        {
          PNL_ERROR ("size mismatch", "pnl_mat_dgemv");
        }
#endif
      for (i=0; i<y->size; i++)
        {
          temp = 0.;
          for (j=0; j<A->n; j++) temp += pnl_mat_get (A, i, j) * pnl_vect_get (x, j);
          yi = pnl_vect_get (y, i);
          yi += alpha * temp; 
          pnl_vect_set (y, i, yi);
        }
    }
  /* Form alpha * A' * x + beta * y */
  else   if (trans=='T' || trans=='t')
    {
#ifndef PNL_RANGE_CHECK_OFF
      if (A->m != x->size || A->n != y->size)
        {
          PNL_ERROR ("size mismatch", "pnl_mat_dgmev");
        }
#endif
      for (j=0; j<A->m; j++)
        {
          temp = alpha * pnl_vect_get (x, j);
          for (i=0; i<y->size; i++)
            {
              yi =  pnl_vect_get (y, i);
              yi += pnl_mat_get (A, j, i) * temp;
              pnl_vect_set (y, i, yi);
            }
        }
    }
  else
    {
      PNL_ERROR ("invalid value for trans", "pnl_mat_dgemv");
    }
}

/**
 * Computes y' A x 

 * @param A : matrix
 * @param x : vector
 * @param y : vector
 * @return y' A x 
 */
double pnl_mat_scalar_prod_A(const PnlMat *A, const PnlVect *x , const PnlVect * y)
{
  double *yarray,*Aarray,*xarray;
  double temp,sum;
  int i,j;
  
#ifndef PNL_RANGE_CHECK_OFF
  if (A->n != x->size || A->m != y->size)
    {
      PNL_ERROR ("size mismatch", "pnl_mat_scalar_prod_A");
    }
#endif

  yarray = y->array;
  Aarray = A->array;
  sum=0;
  for (i=0; i<y->size; i++)
    {
      temp = 0.;
      xarray=x->array;
      for (j=0; j<A->n; j++) 
        {
          temp += (*Aarray)*(*xarray);
          Aarray++; xarray++;
        }
      sum += (*yarray)*temp;
      yarray ++;
    }
  return sum;
}


/**
 * Prints a matrix to a file in a format which can be read by Nsp
 *
 * @param fic a file descriptor.
 * @param M a PnlMat pointer.
 */
void pnl_mat_fprint_nsp (FILE *fic, const PnlMat *M)
{
  int i, j;
  fprintf(fic, "[ ");
  for (i=0;i<M->m;i++)
    {
      for (j=0;j<M->n-1;j++)
        {
          fprintf (fic, "%f, " , pnl_mat_get (M, i, j));
        }
      fprintf (fic, "%f " , pnl_mat_get (M, i, j));
      if (i < M->m - 1) fprintf (fic, "; \n");
    }
  fprintf(fic, " ] ");
}



/**
 * inversion of an upper triangular matrix
 *
 * @param A : at the end of the function contains the
 * inverse of <tt>B</tt>. A must be an already allocated PnlMat
 * @param B : the matrix to be inverted
 */
void pnl_mat_upper_inverse(PnlMat *A, const PnlMat *B)
{
  int i, j, k;
  double y, sum;
  CheckIsSquare(B);
  pnl_mat_resize (A, B->m, B->n);
  pnl_mat_set_double (A, 0.0);
  for (k=0; k<B->m; k++)
    {
      for (i=k; i>=0; i--)
        {
          if (k==i) y=1; else y=0;
          sum = 0.0;
          for (j=i+1; j <B->n; j++)
            sum+= pnl_mat_get (B, i, j) * pnl_mat_get (A, j, k);
          pnl_mat_set (A, i, k, (y - sum) / pnl_mat_get (B, i, i));
        }
    }
}

/**
 * inversion of a lower triangular matrix
 *
 * @param A : at the end of the function contains the
 * inverse of <tt>B</tt>. A must be an already allocated PnlMat
 * @param B : the matrix to be inverted
 */
void pnl_mat_lower_inverse (PnlMat *A, const PnlMat *B)
{
  int k;
  PnlVect *b = pnl_vect_create_from_double (B->m, 0.0);
  PnlVect *x = pnl_vect_create (0);
  CheckIsSquare(B);
  pnl_mat_resize (A, B->m, B->n);
  pnl_mat_set_double (A, 0.0);
  for (k=0; k<B->m; k++)
    {
      pnl_vect_set (b, k, 1.0);
      pnl_mat_lower_syslin (x, B, b);
      pnl_mat_set_col (A, x, k);
      pnl_vect_set (b, k, 0.0);
    }
  pnl_vect_free (&b);
  pnl_vect_free (&x);
}

/**
 * solves an upper triangular linear system
 *
 * @param x already existing PnlVect that contains the solution on exit
 * @param U the matrix of the system
 * @param b right hand side member
 */
void pnl_mat_upper_syslin (PnlVect *x, const PnlMat *U, const  PnlVect *b)
{
  int i, j;
  double sum;
  CheckIsSquare(U);
  pnl_vect_resize (x, U->n);
  for (i=U->m-1; i>=0; i--)
    {
      sum = 0.0;
      for (j=i+1; j <U->n; j++)
        sum+= pnl_mat_get (U, i, j) * pnl_vect_get (x,j);
      *(pnl_vect_lget (x, i)) = (pnl_vect_get (b, i) - sum) / pnl_mat_get (U, i, i);
    }
}


/**
 * solves a lower triangular linear system
 *
 * @param x already existing PnlVect that contains the solution on exit
 * @param L the matrix of the system
 * @param b right hand side member
 */
void pnl_mat_lower_syslin (PnlVect *x, const PnlMat *L, const  PnlVect *b)
{
  int i, j;
  double sum;
  CheckIsSquare(L);
  pnl_vect_resize (x, L->n);
  for (i=0; i<L->n; i++)
    {
      sum = 0.0;
      for (j=0; j < i; j++)
        sum+= pnl_mat_get (L,i,j) * pnl_vect_get (x,j);
      pnl_vect_set(x, i, (pnl_vect_get (b, i) - sum) / pnl_mat_get (L, i, i));
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
 * @param M a PnlMat pointer.
 * @param round_off eigenvalue smaller that -round_off are set to
 * +round_off. Useful for badly contionned matrices.
 */
static void pnl_mat_chol_aux(PnlMat *M, double round_off)
{
  int i=0,j,k;
  double somme;
  int dim=M->m;
  int i_dim=0, j_dim=0;

  CheckIsSquare(M);
  while (i<dim)
    {
      /* i_dim = i*dim;*/
      for (j=i, j_dim=i_dim;j<dim;j++, j_dim+=dim)
        {
          /* j_dim = j*dim; */
          somme=M->array[i_dim+j];
          for (k=i-1;k>=0;k--)
            somme -= M->array[i_dim+k] * M->array[j_dim+k];
          if (i==j)
            {
              if (somme<=-round_off)
                {
                  PNL_ERROR("matrix not positive", "pnl_mat_chol");
                }
              else
                {
                  somme = MAX( somme, round_off );
                  M->array[j_dim+i]=sqrt(somme);
                }
            }
          else M->array[j_dim+i]=somme/M->array[i_dim+i];
          
        }
      i++; i_dim+=dim;
    }
  i=0;j=0;
  for(i=0;i<dim;i++)
    {
      for(j=i+1;j<dim;j++)
        pnl_mat_set (M, i, j, 0.0);
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
 */
void pnl_mat_chol(PnlMat *M)
{
  pnl_mat_chol_aux (M, 0.0);
}

/**
 * Robust Cholesky decomposition. Postivity (up to -PRECISION) is checked
 * during the transformation, but no test of symmetry.
 *
 * Decomposition done in place. The lower part of the matrix
 * contains the cholesky decomposition. The upper part is
 * set to 0.
 *
 * @param M : a PnlMat pointer.
 */
void pnl_mat_chol_robust(PnlMat *M)
{
  pnl_mat_chol_aux (M, PRECISION);
}


/**
 * solves a symmetric DEFINITE POSITIVE linear system using the Cholesky
 * decomposition of the system A x = b
 *
 * @param x already existing PnlVect that contains the solution on exit
 * @param chol the Cholesky decomposition of the system as computed by pnl_mat_chol
 * @param b right hand side member
 */
void pnl_mat_chol_syslin (PnlVect *x, const PnlMat *chol, const  PnlVect *b)
{
  pnl_vect_clone (x,b);
  pnl_mat_chol_syslin_inplace (chol, x);
}

/**
 * solves a symmetric DEFINITE POSITIVE linear system using the Cholesky
 * decomposition of the system A x = b
 *
 * @param chol the Cholesky decomposition of the system as computed by pnl_mat_chol
 * @param b right hand side member. On exit, b contains the solution of the system.
 */
void pnl_mat_chol_syslin_inplace (const PnlMat *chol, PnlVect *b)
{
  int i, j;
  double sum, Bj;
  CheckIsSquare(chol);
  CheckMatVectIsCompatible (chol, b);

  /* solve L y = b, store x into b. */
  for (i=0; i<chol->n; i++)
    {
      sum = 0.0;
      for (j=0; j < i; j++)
        sum+= pnl_mat_get (chol,i,j) * pnl_vect_get (b,j);
      pnl_vect_set (b, i, (pnl_vect_get (b, i) - sum) /  pnl_mat_get (chol, i, i));
    }

  /* solve L' x = y, where y is the solution computed above stored into b */
  for (j=chol->n-1; j>=0; j--)
    {
      Bj = pnl_vect_get(b, j) / pnl_mat_get (chol, j, j);
      pnl_vect_set (b, j, Bj);
      for (i=j-1; i>=0; i--)
        {
          pnl_vect_set (b, i, pnl_vect_get (b, i) - Bj *  pnl_mat_get (chol, j, i));
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
  int i, j, k, i_pivot;
  double aij, ajj, aik, ajk, max;
  int N = A->n;
  CheckIsSquare(A);

  pnl_permutation_init (p);
  for (j=0; j<N; j++)
    {
      /* Find maximum in the j-th column */
      max = fabs (pnl_mat_get (A, j, j));
      i_pivot = j;
        
      for (i=j+1; i<N; i++)
        {
          aij = fabs (pnl_mat_get (A, i, j));
          if (aij > max)
            {
              max = aij;
              i_pivot = i;
            }
        }
        
      if (i_pivot != j)
        {
          pnl_mat_swap_rows (A, j, i_pivot);
          pnl_permutation_swap (p, j, i_pivot);
        }
        
      ajj = pnl_mat_get (A, j, j);
      if (ajj != 0.0)
        {
          for (i=j+1; i<N; i++)
            {
              aij = pnl_mat_get (A, i, j) / ajj;
              pnl_mat_set (A, i, j, aij);
                
              for (k=j+1; k<N; k++)
                {
                  aik = pnl_mat_get (A, i, k);
                  ajk = pnl_mat_get (A, j, k);
                  pnl_mat_set (A, i, k, aik - aij * ajk);
                }
            }
        }
      else
        {
          PNL_ERROR("can't find a LU decomposition", "pnl_mat_lu");
        }
    }
}

/**
 * solves the linear system A x = b with P A = LU. For a symmetric definite
 * positive system, prefer pnl_mat_chol_syslin
 *
 * @param LU a PnlMat containing the LU decomposition of A
 * @param p a PnlPermutation.
 * @param b right hand side member. Contains the solution x on exit
 */
void pnl_mat_lu_syslin_inplace (const PnlMat *LU, const PnlPermutation *p, PnlVect *b)
{
  int i, j;
  double Bi, temp, Aii;
  CheckIsSquare(LU);
  CheckMatVectIsCompatible (LU, b);
  CheckVectMatch (p, b); 

  /* apply the permutation */
  pnl_vect_permute_inplace (b, p);
  
  /* solve L y = b, store x into b. Remember that the diagonal of L is filled
     with 1 */
  for (i=0; i<LU->n; i++)
    {
      temp = 0.;
      for (j=0; j < i; j++)
        {
          temp += pnl_mat_get (LU, i, j) * pnl_vect_get (b, j);
        }
      Bi = pnl_vect_get(b, i);
      pnl_vect_set(b, i, Bi - temp);
    }

  /* solve U x = y, where y is the solution computed above stored into b */
  for (i=LU->m-1; i>=0; i--)
    {
      Bi = pnl_vect_get(b, i);
      /* using two descreasing loops improves performance */
      for (j=LU->m-1; j>i; j--)
        {
          Bi -= pnl_mat_get (LU, i, j) * pnl_vect_get (b, j);
        }
      Aii = pnl_mat_get (LU, i, i);
      Bi /= Aii;
      pnl_vect_set (b, i, Bi);
    }
}

/**
 * solves the linear system A x = b with P A = LU. For a symmetric definite
 * positive system, prefer pnl_mat_chol_syslin
 *
 * @param x existing vector that contains the solution on exit
 * @param LU a PnlMat containing the LU decomposition of A
 * @param b right hand side member
 * @param p a PnlPermutation.
 */
void pnl_mat_lu_syslin (PnlVect *x, const PnlMat *LU,
                        const PnlPermutation *p, const PnlVect *b)
{
  pnl_vect_clone (x, b);
  pnl_mat_lu_syslin_inplace (LU, p, x);
}



/**
 * solves a linear system A x = b using a LU factorization
 * @param x a PnlVect containing the solution on exit (must have already
 * been created )
 * @param A the matrix of the system
 * @param b the r.h.s. member
 */
void pnl_mat_syslin (PnlVect *x, const PnlMat *A, const PnlVect *b)
{
  PnlMat *LU;
  LU = pnl_mat_copy (A);
  pnl_vect_clone (x, b);
  pnl_mat_syslin_inplace (LU, x);
  pnl_mat_free (&LU);
}

/**
 * solves a linear system A x = b using a LU factorization
 * @param A the matrix of the system. On exit contains the LU decomposition of A.
 * @param b the r.h.s. member
 */
void pnl_mat_syslin_inplace (PnlMat *A, PnlVect *b)
{
  PnlPermutation *p;
  CheckIsSquare(A);
  p = pnl_permutation_create (A->m);
  pnl_mat_lu (A, p);
  pnl_mat_lu_syslin_inplace (A, p, b);
  pnl_permutation_free (&p);
}

/**
 * solves a linear system A X = B using a LU factorization where B is a matrix
 * @param A the matrix of the system of size n x n. On exit contains the LU decomposition
 * @param B the r.h.s. matrix of the system of size n x m. On exit B contains
 * the solution X
 */
void pnl_mat_syslin_mat (PnlMat *A,  PnlMat *B)
{
  PnlPermutation *p;
  PnlVect *b;
  int i;

  CheckIsSquare(A);
  b = pnl_vect_create (0);
  p = pnl_permutation_create (A->m);
  pnl_mat_lu (A, p);

  for (i=0; i<B->n; i++)
    {
      pnl_mat_get_col (b, B, i);
      pnl_mat_lu_syslin_inplace (A, p, b);
      pnl_mat_set_col (B, b, i);
    }
  pnl_vect_free (&b);
  pnl_permutation_free (&p);
}

/*
 * This function solves a linear system A * X = B in place using a LU
 * factorization. This implementation is strongly inspired from CBlas and does
 * not require any working space unlike pnl_mat_syslin_mat which is faster but
 * does require a working space of size equal to the number of lines of B
 *
 * This function is currently not exported in any .h file 
 */
void pnl_mat_syslin_mat2 (PnlMat *A,  PnlMat *B)
{
  PnlPermutation *p;
  int i, j, k;

  CheckIsSquare(A);
  p = pnl_permutation_create (A->m);
  pnl_mat_lu (A, p);
  
  /* solve L Y = B, store the result in B */
  for (i=0; i<B->m; i++)
    {
      for (k=i+1; k < B->m; k++)
        {
          double aki = pnl_mat_get (A, k, i);
          for (j=0; j<B->n; j++)
            {
              *(pnl_mat_lget (B, k, j)) -= aki * pnl_mat_get (B, i, j);
            }
        }
    }
  
  /* solve L' X = Y, where y is the solution computed above stored into B */
  for (i=B->m-1; i>=0; i--)
    {
      double aii = pnl_mat_get (A, i, i);
      for (j=0; j<B->n; j++)
        {
          *(pnl_mat_lget (B, i, j)) /= aii;
        }
      for (k=0; k < i; k++)
        {
          double aki = pnl_mat_get (A, k, i);
          for (j=0; j<B->n; j++)
            {
              *(pnl_mat_lget (B, k, j)) -= aki * pnl_mat_get (B, i, j);
            }
        }
    }
}

/**
 * solves a linear system A X = B using a Cholesky factorization of A where B
 * is a matrix. Note that A must be symmetrix positive definite
 *
 * @param A the matrix of the system of size n x n. On exit contains the
 * Cholesky decomposition 
 * @param B the r.h.s. matrix of the system of size n x m. On exit B contains
 * the solution X
 */
void pnl_mat_chol_syslin_mat (PnlMat *A,  PnlMat *B)
{
  PnlVect *b;
  int i;

  CheckIsSquare(A);
  b = pnl_vect_create (0);
  pnl_mat_chol (A);

  for (i=0; i<B->n; i++)
    {
      pnl_mat_get_col (b, B, i);
      pnl_mat_chol_syslin_inplace (A, b);
      pnl_mat_set_col (B, b, i);
    }
  pnl_vect_free (&b);
}

/*
 * This function solves a linear system A * X = B in place when A>0 using a
 * Cholecky factorization. This implementation is strongly inspired from CBlas
 * and does not require any working space unlike pnl_mat_chol_syslin_mat which
 * is faster but does require a working space of size equal to the number of
 * lines of B
 *
 * This function is currently not exported in any .h file 
 */
void pnl_mat_chol_syslin_mat2 (PnlMat *A,  PnlMat *B)
{
  int i, j, k;
  CheckMatMatch(A, B);
  pnl_mat_chol (A);

  /* solve L Y = B, store the result in B */
  for (i=0; i<B->m; i++)
    {
      double aii = pnl_mat_get (A, i, i);
      for (j=0; j<B->n; j++)
        {
          *(pnl_mat_lget (B, i, j)) /= aii;
        }
      for (k=i+1; k < B->m; k++)
        {
          double aki = pnl_mat_get (A, k, i);
          for (j=0; j<B->n; j++)
            {
              *(pnl_mat_lget (B, k, j)) -= aki * pnl_mat_get (B, i, j);
            }
        }
    }
  
  /* solve L' X = Y, where y is the solution computed above stored into B */
  for (i=B->m-1; i>=0; i--)
    {
       double aii = pnl_mat_get (A, i, i);
      for (j=0; j<B->n; j++)
        {
          *(pnl_mat_lget (B, i, j)) /= aii;
        }
      for (k=0; k < i; k++)
        {
           double aik = pnl_mat_get (A, i, k); /* the upper of the matrix is 0, need
                                                  to take A(i,k) instead of A(k,i) */
          for (j=0; j<B->n; j++)
            {
              *(pnl_mat_lget (B, k, j)) -= aik * pnl_mat_get (B, i, j);
            }
        }
    }
}


/**
 * computes the inverse of a matrix A given its LU decomposition and the
 * associated permutation
 *
 * @param LU LU decomposition of P A computed by pnl_mat_lu.
 * @param p a PnlPermutation.
 * @param inverse a PnlMat (already allocated). contains
 * \verbatim A^-1 \endverbatim on exit.
 */
void pnl_mat_lu_inverse (PnlMat *inverse, const PnlMat *LU, const PnlPermutation*p)
{
  int k;
  PnlVect *b = pnl_vect_create_from_double (LU->m, 0.0);
  PnlVect *x = pnl_vect_create (0);
  pnl_mat_resize (inverse, LU->m, LU->n);
  for (k=0; k<LU->m; k++)
    {
      pnl_vect_set (b, k, 1.0);
      pnl_mat_lu_syslin (x, LU, p, b);
      pnl_mat_set_col (inverse, x, k);
      pnl_vect_set (b, k, 0.0);
    }
  pnl_vect_free (&b);
  pnl_vect_free (&x);
}

/**
 * Matrix exponential B = exp( A)
 *
 * @param A a matrix
 * @param B contains exp(A) on return
 *
 * This algorithms has been written in C by Jérôme Lelong
 * following the Fortran implementation given in EXPOKIT
 * Roger B. Sidje Department of Mathematics, University of Queensland 
 * Brisbane, QLD-4072, Australia, (c) 1996-1999 All Rights Reserved
 * The work resulting from EXPOKIT has been published in ACM-Transactions 
 * on Mathematical Software, 24(1):130-156, 1998.
 *
 *  The bibtex record of the citation:
 *
 * \verbatim
 * ARTICLE{EXPOKIT,
 *        AUTHOR  = {Sidje, R. B.},
 *        TITLE   = {{\sc Expokit.} {S}oftware {P}ackage for
 *                 {C}omputing {M}atrix {E}xponentials},
 *        JOURNAL = {ACM Trans. Math. Softw.},
 *        VOLUME  = {24},
 *        NUMBER  = {1},
 *        PAGES   = {130-156}
 *        YEAR    = {1998}
 * }
 * \endverbatim
 */
void pnl_mat_exp (PnlMat *B, const PnlMat *A)
{
  int i,j,k,ih2,ip,iq,iused,ifree,iodd,icoef,ideg,iput,iget,lwork, ns;
  double *work;
  PnlMat Mwork, Mwork1, Mwork2;
  double hnorm,scale,scale2,cp,cq;
  
  CheckIsSquare (A);
  
  /* Mwork is used as a container to pnl_mat_xxx routines */
  
  pnl_mat_resize (B, A->m, A->n);
  
  icoef = 0;
  ideg = 6;
  ih2 = icoef + (ideg+1);
  ip  = ih2 + A->mn;
  iq  = ip + A->mn;
  ifree = iq + A->mn;
  lwork = 4*A->m*A->m + ideg + 1;
  if ( (work = malloc (lwork * sizeof(double))) == NULL ) abort();

  /*  scaling: seek ns such that ||t*H/2^ns|| < 1/2;  */
  /*  and set scale = t/2^ns ... */
  for (i=0; i<A->m; i++) work[i] = 0.;
  for (j=0; j<A->m; j++)
    {
      for (i=0; i<A->m; i++)
        work[i] += work[i] + fabs( pnl_mat_get(A,i,j) );      
    }

  hnorm = 0.;
  for (i=0; i<A->m; i++)
    {
      hnorm = MAX( hnorm, work[i] );
    }

  if (hnorm == 0.)
    /* matrix is full of zeros */
    {
      pnl_mat_set_id (B);
      return;
    }
  ns = MAX( 0, (int)(log(hnorm)/log(2.)) + 2 );
  scale = 1. / pow(2., ns);
  scale2 = scale*scale;

  /*  compute Pade coefficients ... */
  i = ideg+1;
  j = 2*ideg+1;
  work[icoef] = 1.0;
  for (k=1; k<ideg+1; k++)
    {
      work[icoef+k] = (work[icoef+k-1]*( i-k )) / ( k*(j-k) );
    }

  Mwork = pnl_mat_create_wrap_array (&(work[ih2]), A->m, A->n);
  pnl_mat_dgemm ('N', 'N', scale2, A, A, 0., &Mwork);  /* H2 = scale2*H*H */
  
  /* initialize p (numerator) and q (denominator) */
  cp = work[icoef+ideg-1];
  cq = work[icoef+ideg];
  for (j=0; j<A->m; j++)
    {
      for (i=0; i<A->m; i++)
        {
          work[ip + j*A->m + i] = 0.;
          work[iq + j*A->m + i] = 0.;
        }
      work[ip + j*(A->m+1)] = cp; /* sets the diagonal */
      work[iq + j*(A->m+1)] = cq; /* sets the diagonal */
    }

  /* Apply Horner rule */
  iodd = 1;
  k = ideg - 1;
  while (k > 0)
    {
      iused = iodd*iq + (1-iodd)*ip;
      Mwork = pnl_mat_create_wrap_array (&(work[ifree]), A->m, A->m);
      Mwork1 = pnl_mat_create_wrap_array (&(work[iused]), A->m, A->m);
      Mwork2 = pnl_mat_create_wrap_array (&(work[ih2]), A->m, A->m);
      pnl_mat_dgemm('N', 'N', 1., &Mwork1, &Mwork2, 0., &Mwork);
      
      for (j=0; j<A->m; j++)
        /* add work[icoef+k-1]; to the diagonal */
        work[ifree+j*(A->m+1)] += work[icoef+k-1];
      ip = (1-iodd)*ifree + iodd*ip;
      iq = iodd*ifree + (1-iodd)*iq;
      ifree = iused;
      iodd = 1-iodd;
      k--;
    }
 
  /* Obtain (+-)(I + 2 * (p/q))  */
  Mwork = pnl_mat_create_wrap_array (&(work[ifree]), A->m, A->m);
  Mwork1 = pnl_mat_create_wrap_array (&(work[ip]), A->m, A->m);
  pnl_mat_dgemm ('N', 'N', scale, &Mwork1, A, 0., &Mwork);
  ip = ifree;

  Mwork1 = pnl_mat_create_wrap_array (&(work[ip]), A->m, A->m);
  Mwork2 = pnl_mat_create_wrap_array (&(work[iq]), A->m, A->m);

  pnl_mat_axpy (-1., &Mwork1, &Mwork2 );
  pnl_mat_syslin_mat (&Mwork2, &Mwork1);

  pnl_mat_mult_double (&Mwork1, 2.);
  for (j=0; j<A->m; j++) work[ip+j*(A->m+1)] += 1.;
  iput = ip;

  if (ns == 0 && iodd == 1)
    {
      pnl_mat_mult_double (&Mwork1, -1.);
    }
  else
    {
      /* squaring : exp(t*H) = (exp(t*H))^(2^ns) */
      iodd = 1;
      for (k=0; k<ns; k++)
        {
          iget = iodd*ip + (1-iodd)*iq;
          iput = (1-iodd)*ip + iodd*iq;
          Mwork = pnl_mat_create_wrap_array (&(work[iget]), A->m, A->m);
          Mwork1 = pnl_mat_create_wrap_array (&(work[iput]), A->m, A->m);
          pnl_mat_dgemm ('N', 'N', 1., &Mwork, &Mwork, 0., &Mwork1);
          iodd = 1-iodd;
        }
    }

  /* the solution is located at work[iput] */
  memcpy (B->array, &(work[iput]), A->mn * sizeof(double));
  
  free (work); work = NULL;
}


