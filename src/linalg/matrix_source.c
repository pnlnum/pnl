/*
 * This code was originally written for matrices of doubles and has been
 * translated using a template system by
 *    David Pommier <pommier.david@gmail.com>
 *    Jérôme Lelong <jerome.lelong@gmail.com>
 */


/************************************************************************/
/* Copyright Jérôme Lelong <jerome.lelong@gmail.com>                    */
/* Copyright Céline Labart                                              */
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
 * Initialize a PnlMat
 */
void FUNCTION(pnl_mat,init)(TYPE(PnlMat) *o)
{
  o->object.parent_type = PNL_TYPE_MATRIX;
  o->object.type = CONCAT2(PNL_TYPE_MATRIX_, BASE_TYPE);
  o->object.label = FUNCTION(pnl_mat,label);
  o->object.destroy = (DestroyFunc *) pnl_mat_object_free;
  o->object.constructor = (NewFunc *) FUNCTION(pnl_mat,new);
  o->object.clone = (CloneFunc *) FUNCTION(pnl_mat,clone);
  o->object.copy = (CopyFunc *) FUNCTION(pnl_mat,copy);
  o->m = 0;
  o->n = 0;
  o->mn = 0;
  o->mem_size = 0;
  o->owner = 1;
  o->array = NULL;
}

/**
 * Create an empty PnlMat
 * @return a TYPE(PnlMat) pointer
 */
TYPE(PnlMat) * FUNCTION(pnl_mat,new)()
{
  TYPE(PnlMat) *o;
  if ( (o = (TYPE(PnlMat) *) pnl_mat_object_new ()) == NULL) return NULL;
  FUNCTION(pnl_mat,init)(o);
  return o;
}

/**
 * Test if 2 matrices are equal
 *
 * @param M1 a matrix
 * @param M2 a matrix
 * @return  TRUE or FALSE
 */
int FUNCTION(pnl_mat,eq)(const TYPE(PnlMat) *M1, const TYPE(PnlMat) *M2)
{
  int i;
  if ( (M1->m != M2->m) || (M1->n != M2->n) ) return FALSE;
  for ( i=0 ; i<M1->mn ; i++ )
    {
      if ( NEQ(M1->array[i],M2->array[i]) ) return FALSE;
    }
  return TRUE;
}

/**
 * Test if all entries of a  matrix are equal to a given BASE element
 *
 * @param M a matrix
 * @param x a BASE type element
 * @return  TRUE or FALSE
 */
int FUNCTION(pnl_mat,eq_all)(const TYPE(PnlMat) *M, BASE x)
{
  int i;
  for ( i=0 ; i<M->mn ; i++ )
    {
      if ( NEQ(M->array[i],x) ) return FALSE;
    }
  return TRUE;
}

/**
 * Create a PnlMat
 * @param m number of rows
 * @param n number of columns
 * @return a TYPE(PnlMat) pointer. The content of the array is not set.
 */
TYPE(PnlMat) * FUNCTION(pnl_mat,create)(int m, int n)
{
  TYPE(PnlMat) *M;
  if ((M = FUNCTION(pnl_mat,new)()) == NULL) return NULL;
  if (m>0 && n>0)
    {
      M->m = m;
      M->n = n;
      M->mn = m * n;
      if ((M->array = malloc (M->mn * sizeof(BASE))) == NULL) return NULL;
      M->mem_size = M->mn;
      M->owner = 1;
    }
  return M;
}

/**
 * Create a TYPE(PnlMat) filled with a constant value
 * @param m number of rows
 * @param n number of columns
 * @param x used to fill the matrix with
 * @return a TYPE(PnlMat) pointer
 */
TYPE(PnlMat)* FUNCTION(pnl_mat,create_from_scalar)(int m, int n, BASE x)
{
  TYPE(PnlMat) *M;
  int i;

  if ((M=FUNCTION(pnl_mat,create)(m, n))==NULL) return NULL;

  for (i=0; i<M->mn; i++) { M->array[i] = x; }
  return M;
}

/**
 * Create a new TYPE(PnlMat) pointer.
 *
 * @param m number of rows
 * @param n number of columns
 * @return  a TYPE(PnlMat) pointer with all entries set to 0
 */
TYPE(PnlMat)* FUNCTION(pnl_mat,create_from_zero)(int m, int n)
{
  TYPE(PnlMat) *M = FUNCTION(pnl_mat,create)(m, n);
  FUNCTION(pnl_mat,set_zero) (M);
  return M;
}

/**
 * Create a TYPE(PnlMat) from a C array
 * @param m number of rows
 * @param n number of columns
 * @param x an array of BASE used to fill the TYPE(PnlMat). should be of length
 * <tt> mxn </tt>. No test can be done about it.
 * @return a TYPE(PnlMat) pointer
 */
TYPE(PnlMat)* FUNCTION(pnl_mat,create_from_ptr)(int m, int n, const BASE* x)
{
  TYPE(PnlMat) *M;

  if ((M=FUNCTION(pnl_mat,create)(m, n))==NULL) return NULL;

  memcpy(M->array, x, M->mn*sizeof(BASE));
  return M;
}

/**
 * Create a new TYPE(PnlMat) using a list of values.
 *
 * @param m number of rows
 * @param n number of columns
 * @param ... is a list of values o type BASE used to fill the array. The list
 * must be of length size. The data are stored row-wise.
 * @return  a TYPE(PnlVect)pointer
 */
TYPE(PnlMat) * FUNCTION(pnl_mat,create_from_list)(int m, int n, ...)
{
  TYPE(PnlMat) * M;
  va_list ap;
  int i;

  if ((M=FUNCTION(pnl_mat,create)(m, n))==NULL) return NULL;
  va_start (ap, n);

  for (i=0; i<M->mn; i++)
    {
      BASE val ;
      val = va_arg (ap, BASE);
      M->array[i] = val;
    }
  va_end(ap);
  return M;
}

/**
 * Wrap a C array into a PnlMat
 * @param m number of rows
 * @param n number of columns
 * @param x an array of BASE used to fill the TYPE(PnlMat). should be of length
 * <tt> mxn </tt>. No test can be done about it.
 * @return a TYPE(PnlMat) pointer
 */
TYPE(PnlMat) FUNCTION(pnl_mat,wrap_array)(const BASE* x, int m, int n)
{
  TYPE(PnlMat) M;
  FUNCTION(pnl_mat,init)(&M);
  M.m = m;
  M.n = n;
  M.mn = m*n;
  M.mem_size = 0;
  M.owner = 0;
  M.array = (BASE *) x;
  return M;
}

/**
 * Read a matrix from a file and creates the corresponding TYPE(PnlMat)
 * @param file  the file to be read
 * @ return a TYPE(PnlMat)
 */
TYPE(PnlMat)* FUNCTION(pnl_mat,create_from_file )(const char * file)
{
  char car, prev = '\0', empty=1;
  TYPE(PnlMat) *M;
  int m, n, count, mn;
  BASE *data;
  FILE *FIC = NULL;

  if ((FIC = fopen( file, "r")) == NULL )
    {
      PNL_ERROR("Cannot open file", "FUNCTION(pnl_mat,create_from_file)");
    }

  /* first pass to determine dimensions */
  m = 0; n = 1;
  while((car=fgetc(FIC))!='\n' && car != EOF)
    {
      if (isdigit(car) || car == '-' || car == '.')
        {
          empty = 0;
          if (prev == ' ' || prev == '\t' ) ++n;
        }
      prev = car;
   }
  n = n / MULTIPLICITY;
  /*if (!empty && car =='\n' && isdigit(prev)) ++n; */
  if (!empty) ++m;
  empty = 1;
  while((car=fgetc(FIC))!= EOF)
    {
      if( car=='\n' )
        {
          if (!empty) { ++m; empty = 1;}
          else break;
        }
      else if (empty && isdigit(car)) empty=0;
    }
  /* If no data has been read */
  if (m==0 || n==0)
    {
      M = FUNCTION(pnl_mat,create)(0,0);
      fclose (FIC);
      return M;
    }

  mn = m*n;
  if ((M = FUNCTION(pnl_mat,create )(m,n))==NULL)
    {
      PNL_ERROR("Allocation error", "FUNCTION(pnl_mat,create_from_file)");
    }

  data = FUNCTION(pnl_mat,lget )(M, 0, 0);

  /* second pass to read data */
  rewind (FIC);
  count = 0;
  while ( fscanf (FIC,IN_FORMAT, IN_PUT_FORMAT(data)) > 0 && count < mn ) { data++; count++;}
  fclose (FIC);
  return M;
}

/**
 * Free a TYPE(PnlMat)
 *
 * @param v adress of a TYPE(PnlMat)*. v is set to NULL at exit.
 */
void FUNCTION(pnl_mat,free)(TYPE(PnlMat) **v)
{
  PnlMatObject *o;
  o = PNL_MAT_OBJECT(*v);
  pnl_mat_object_free (&o);
  *v = NULL;
}

/**
 * Copy a TYPE(PnlMat)
 *
 * @param v : a constant TYPE(PnlMat) pointer
 * @return a TYPE(PnlMat) pointer initialised with v
 */
TYPE(PnlMat)* FUNCTION(pnl_mat,copy)(const TYPE(PnlMat) *v)
{
  TYPE(PnlMat) *ret;

  if ((ret=FUNCTION(pnl_mat,create)(v->m, v->n))==NULL) return NULL;
  memcpy(ret->array, v->array, sizeof(BASE)*ret->mn);
  return ret;
}

/**
 * Clone a TYPE(PnlMat)
 *
 * @param M : a constant TYPE(PnlMat) pointer
 * @param clone : needs to be already allocated. NOT TESTED
 */
void FUNCTION(pnl_mat,clone)(TYPE(PnlMat) *clone, const TYPE(PnlMat) *M)
{
  PNL_CHECK(clone->owner == 0 && clone->mn != M->mn, "owner == 0 and size mismatch", "pnl_mat_clone");
  FUNCTION(pnl_mat,resize)(clone, M->m,M->n);
  memcpy(clone->array, M->array, sizeof(BASE)*M->mn);
}

/**
 * Resize a TYPE(PnlMat).  If the new size is smaller than the current one, no
 * memory is freed. If the new size is larger than the current mem_size, a new
 * pointer is allocated. The old data are kept.
 *
 * @param M : a pointer to an already existing TYPE(PnlMat)
 * @param m : new nb rows
 * @param n : new nb columns
 *
 * @return OK or FAIL. When returns OK, the matrix is changed.
 */
int FUNCTION(pnl_mat,resize)(TYPE(PnlMat) *M, int m, int n)
{
  return pnl_mat_object_resize (PNL_MAT_OBJECT(M), m, n);
}

/**
 * in-place set matrix constant
 *
 * @param lhs : left hand side matrix
 * @param x : scalar
 * @return  lhs = x
 */
void FUNCTION(pnl_mat,set_all) (TYPE(PnlMat) *lhs, BASE x)
{
  int i;
  for (i=0; i<lhs->mn; i++) lhs->array[i] = x;
}

/**
 * in-place set matrix to zero
 *
 * @param lhs : left hand side matrix
 * @return  lhs = x
 */
void FUNCTION(pnl_mat, set_zero) (TYPE(PnlMat) *lhs)
{
  memset (lhs->array, 0, sizeof(BASE) * lhs->mn);
}

/**
 * in-place set matrix to identity
 *
 * @param lhs : left hand side matrix
 */
void FUNCTION(pnl_mat,set_id)(TYPE(PnlMat) *lhs)
{
  int i;
  BASE zero = ZERO, one = ONE;
  CheckIsSquare (lhs);
  for (i=0; i<lhs->mn; i++) lhs->array[i] = zero;
  for (i=0; i<lhs->m; i++)
    {
      PNL_MLET(lhs, i, i) = one;
    }
}

/**
 * in-place set matrix to identity
 *
 * @param M a matrix
 * @param x the value used to fill the diagonal
 * @param d the index of the diagonal (if d>0, we consider the d-th upper
 * diagonal, if d<0, we consider the (-d)-th lower diagonal)
 *
 */
void FUNCTION(pnl_mat,set_diag)(TYPE(PnlMat) *M, BASE x, int d)
{
  int i;
  CheckIsSquare (M);
  for ( i=MAX(0, -d) ; i<M->m + MIN(0, -d) ; i++ )
    {
      PNL_MLET (M, i, i+d) = x;
    }
}

/**
 * Fill in a TYPE(PnlMat) from a C array
 * @param M an already existing TYPE(PnlMat)
 * @param x an array of BASE used to fill the TYPE(PnlMat). should be of length
 * <tt> mxn </tt>. No test can be done on this.
 * @return a TYPE(PnlMat) pointer
 */
void FUNCTION(pnl_mat,set_from_ptr)(TYPE(PnlMat) *M, const BASE* x)
{
  memcpy(M->array, x, M->mn*sizeof(BASE));
}

/**
 *  Non square matrix transposition
 *
 * @param M a TYPE(PnlMat) pointer.
 * @param tM on exit contains the transpose of M
 */
void FUNCTION(pnl_mat,tr)(TYPE(PnlMat)*tM, const TYPE(PnlMat) *M)
{
  int i,j;
  FUNCTION(pnl_mat,resize)(tM, M->n,M->m);
  for ( i=0 ; i<tM->m ; i++)
    {
      for ( j=0 ; j<tM->n ; j++)
        {
          PNL_MLET(tM, i, j) = PNL_MGET(M, j, i);
        }
    }
}

/**
 *  Non square matrix transposition
 *
 * @param M a TYPE(PnlMat) pointer.
 * @return the transpose of M
 */
TYPE(PnlMat)* FUNCTION(pnl_mat,transpose)(const TYPE(PnlMat) *M)
{
  TYPE(PnlMat)*tM;
  tM = FUNCTION(pnl_mat,create)(0,0);
  FUNCTION(pnl_mat,tr)(tM, M);
  return tM;
}

/**
 * in-place transposition of square matrices
 *
 * @param M : a PnlMat pointer.
 */
void FUNCTION(pnl_mat,sq_transpose) (TYPE(PnlMat) *M)
{
  BASE a, b;
  int i,j;
  CheckIsSquare(M);
  for (i=0;i<M->m;i++)
    {
      for (j=i+1;j<M->n;j++)
        {
          b = PNL_MGET (M, j, i);
          a = PNL_MGET (M, i, j);
          PNL_MLET (M, i, j) = b;
          PNL_MLET(M, j, i) = a;
        }
    }
}

/**
 * Compute the trace of a square matrix
 *
 * @param M : a PnlMat pointer.
 */
BASE FUNCTION(pnl_mat,trace) (const TYPE(PnlMat) *M)
{

  int i;
  BASE trace = ZERO;
  CheckIsSquare(M);
  for ( i=0 ; i<M->m ; i++ )
    {
      trace = PLUS (trace, PNL_MGET(M, i, i));
    }
  return trace;
}

/**
 * Print a matrix to a file
 *
 * @param fic a file descriptor.
 * @param M a TYPE(PnlMat) pointer.
 */
void FUNCTION(pnl_mat,fprint )(FILE *fic, const TYPE(PnlMat) *M)
{
  int i, j;
  for (i=0;i<M->m;i++)
    {
      for (j=0;j<M->n;j++)
        {
          fprintf (fic,OUT_FORMAT " ",OUT_PUT_FORMAT(PNL_MGET(M, i, j)));
        }
      fprintf (fic, "\n");
    }
}

/**
 * Print a matrix to the standard output
 *
 * @param M a TYPE(PnlMat) pointer.
 */
void FUNCTION(pnl_mat,print )(const TYPE(PnlMat) *M)
{ FUNCTION(pnl_mat,fprint)(stdout, M);}

/**
 * Print a TYPE(PnlMat) to a file in a format compatible with Nsp
 *
 * @param fic a file descriptor
 * @param M a(constant) TYPE(PnlMat)ptr.
 */
void FUNCTION(pnl_mat, fprint_nsp)(FILE *fic, const TYPE(PnlMat) * M)
{
  int i, j;
  fprintf(fic,"[ ");
  for (i=0; i<M->m; i++)
    {
      for (j=0;j<M->n;j++)
        {
          fprintf(fic,OUT_FORMAT,OUT_PUT_FORMAT(PNL_MGET(M, i, j)));
          if (j != M->n-1) fprintf (fic,", ");
        }
      if (i != M->m-1) fprintf(fic,"; \n  ");
    }
  fprintf(fic," ]; \n");
}

/**
 * Print a TYPE(PnlMat) in a format compatible with Nsp
 *
 * @param M a(constant) TYPE(PnlMat)ptr.
 */
void FUNCTION(pnl_mat, print_nsp)(const TYPE(PnlMat) * M)
{ FUNCTION(pnl_mat,fprint_nsp)(stdout, M);}

/**
 * Create a square matrix by specifying its diagonal terms.
 *
 *@param V : a TYPE(PnlVect)
 *@return a squared matrix whose diagonal is the vector V
 *and null elsewhere
 */
TYPE(PnlMat)* FUNCTION(pnl_mat,create_diag)(const TYPE(PnlVect) *V)
{
  return FUNCTION(pnl_mat,create_diag_from_ptr )(V->array, V->size);
}

/**
 * @param x a double array
 * @param d size of the square matrix to create
 * @return a squared matrix whose diagonal is the array x
 * and 0 anywhereelse
 */
TYPE(PnlMat)* FUNCTION(pnl_mat,create_diag_from_ptr)(const BASE x[], int d)
{
  TYPE(PnlMat)* lhs;
  int i;
  if((lhs=FUNCTION(pnl_mat,create_from_scalar)(d,d,ZERO))==NULL)
    return NULL;
  for(i=0;i<d;i++) { PNL_MLET(lhs, i, i) =  x[i]; }
  return lhs;
}

/**
 * replace the ith row of a matrix M by a vector V
 *
 * @param M : left hand side matrix (needs to be already allocated)
 * @param V : right hand side vector
 * @param i : integer corresponding to the row of M which
 * becomes V
 * @return matrix M which contains whose ith row is becoming V
 */
void FUNCTION(pnl_mat,set_row)(TYPE(PnlMat) *M, const TYPE(PnlVect) *V, int i)
{
  PNL_CHECK(i>M->m,"index out of range", "pnl_mat_set_row");
  PNL_CHECK(M->n != V->size,"incompatible data", "pnl_mat_set_row");
  memcpy (&(M->array[i*M->n]), V->array, V->size*sizeof(BASE));
}

/**
 * Swap two rows of a matrix
 *
 * @param M TYPE(PnlMat)
 * @param i first index
 * @param j second index
 */
void FUNCTION(pnl_mat,swap_rows )(TYPE(PnlMat) *M, int i, int j)
{
  int k, p;
  register BASE Mik, Mjk;
  PNL_CHECK (i >= M->m,"first index is out of range", "FUNCTION(pnl_mat,swap_rows");
  PNL_CHECK (j >= M->m,"second index is out of range", "FUNCTION(pnl_mat,swap_rows");
  p = M->n % 3;
  for (k=0; k<p; k++)
    {
      Mik = PNL_MGET(M, i, k);
      Mjk = PNL_MGET(M,j, k);
      PNL_MLET(M, i, k) = Mjk;
      PNL_MLET(M, j, k) = Mik;
    }
  for (k=p; k<M->n; k+=3)
    {
      Mik = PNL_MGET(M, i, k);
      Mjk = PNL_MGET(M,j, k);
      PNL_MLET(M, i, k) = Mjk;
      PNL_MLET(M, j, k) = Mik;

      Mik = PNL_MGET(M, i, k+1);
      Mjk = PNL_MGET(M,j, k+1);
      PNL_MLET(M, i, k+1) = Mjk;
      PNL_MLET(M, j, k+1) = Mik;

      Mik = PNL_MGET(M, i, k+2);
      Mjk = PNL_MGET(M,j, k+2);
      PNL_MLET(M, i, k+2) = Mjk;
      PNL_MLET(M, j, k+2) = Mik;
    }
}

/**
 * replace the ith column of a matrix M by a vector V
 *
 * @param M : left hand side matrix (needs to be already allocated)
 * @param V : right hand side vector
 * @param j : integer corresponding to the column of M which
 * becomes V
 * @return matrix M which contains whose ith row is becoming V
 */
void FUNCTION(pnl_mat,set_col)(TYPE(PnlMat) *M, const TYPE(PnlVect) *V, int j)
{
  int i;
  PNL_CHECK (j>=M->n, "index of range", "mat_set_col");
  PNL_CHECK (M->m!=V->size, "incompatible size", "mat_set_col");
  for (i=0; i<M->m; i++)
    {
      PNL_MLET(M, i, j) = PNL_GET(V, i);
    }
}

/**
 * Extract the ith row of M. Makes a copy.
 *
 * @param V a PnlVect which contains the row on exit
 * @param M a PnlMat
 * @param i the index of the row to be extracted
 */
void FUNCTION(pnl_mat,get_row)(TYPE(PnlVect) *V, const TYPE(PnlMat) *M, int i)
{
  PNL_CHECK(i>=M->m,"index out of range", "pnl_mat_get_row");
  FUNCTION(pnl_vect,resize)(V, M->n);
  memcpy (V->array,&(M->array[i*M->n]),M->n*sizeof(BASE));
}

/**
 * Extract a row of a matrix and wrap it into a vector.
 * @param M a matrix
 * @param i the index of the row to be extracted
 * @return a vector (not a pointer) whose array pointer is the address of the
 * first element of the ith row of M. No copying is done.
 */
TYPE(PnlVect) FUNCTION(pnl_vect,wrap_mat_row)(const TYPE(PnlMat) *M, int i)
{
  TYPE(PnlVect) V;
  PNL_CHECK (i>=M->m,"index out of range", "pnl_vect_wrap_mat_row");
  FUNCTION(pnl_vect,init)(&V);
  V.size = M->n;
  V.mem_size = 0;
  V.owner = 0;
  V.array = &(M->array[i*M->n]);
  return V;
}

/**
 * Wrap a vector into a PnlMat
 * @param V a vector
 * @return a matrix (not a pointer) whose array pointer is the address of the
 * first element of the vector V. No copying is done.
 */
TYPE(PnlMat) FUNCTION(pnl_mat,wrap_vect)(const TYPE(PnlVect) *V)
{
  TYPE(PnlMat) M;
  FUNCTION(pnl_mat,init)(&M);
  M.m = V->size;
  M.n = 1;
  M.mn = V->size;
  M.mem_size = 0;
  M.owner = 0;
  M.array = V->array;
  return M;
}

/**
 * Copy the ith column of M into V
 *
 * @param V a PnlVect resized within the function
 * @param M the matrix to be extracted
 * @param j the index of the column
 */
void FUNCTION(pnl_mat, get_col)(TYPE(PnlVect) *V, const TYPE(PnlMat) *M, int j)
{
  int i;
  PNL_CHECK (j >= M->n, "index out of range", "pnl_mat_get_col");
  FUNCTION(pnl_vect, resize)(V,M->m);
  for (i=0; i<M->m; i++)
    {
      PNL_LET (V, i) = PNL_MGET (M, i, j);
    }
}

/**
 * in-place map function
 *
 * @param x a matrix
 * @param f : the function to be applied term by term
 * @return  x = f(x)
 */
void FUNCTION(pnl_mat, map_inplace)(TYPE(PnlMat) *x, BASE(*f)(BASE))
{
  BASE  *xi;
  int i;
  xi = x->array;
  for (i=0; i<x->mn; i++)
    {
      *xi = f(*xi);
      xi++;
    }
}

/**
 * map matrix componentwise
 *
 * @param lhs : each component lhs(i) contains f(rhs(i))
 * @param rhs : right hand side vector
 * @param f   : BASE function
 */
void FUNCTION(pnl_mat,map)(TYPE(PnlMat) *lhs, const TYPE(PnlMat) *rhs, BASE(*f)(BASE))
{
  FUNCTION(pnl_mat,clone)(lhs, rhs);
  FUNCTION(pnl_mat,map_inplace)(lhs, f);
}

/**
 * map matrix componentwise
 *
 * @param lhs  each component lhs(i) contains f(lhs(i),rhs(i))
 * @param rhs  a matrix
 * @param f a function taking two real numbers and returning a real number
 */
void FUNCTION(pnl_mat,map_mat_inplace)(TYPE(PnlMat) *lhs, const TYPE(PnlMat) *rhs, BASE(*f)(BASE,BASE))
{
  int i;
  CheckMatMatch(lhs, rhs);
  for (i=0; i<lhs->mn; i++)
    {
      lhs->array[i] = (*f) (lhs->array[i], rhs->array[i]);
    }
}

/**
 * map matrix componentwise
 *
 * @param lhs : each component lhs(i) contains f(rhs1(i),rhs2(i))
 * @param rhs1 a matrix
 * @param rhs2 a matrix
 * @param f a function taking two real numbers and returning a real number
 */
void FUNCTION(pnl_mat,map_mat)(TYPE(PnlMat) *lhs, const TYPE(PnlMat) *rhs1, const TYPE(PnlMat) *rhs2, BASE(*f)(BASE,BASE))
{
  FUNCTION(pnl_mat,clone) (lhs, rhs1);
  FUNCTION(pnl_mat,map_mat_inplace) (lhs, rhs2, f);
}

typedef struct {
  union {BASE x; TYPE(PnlMat) *V ;};
  char type;
} TYPE(cell);
/**
 * Find the indices (i, j) for which f(M(i,j)) == 1
 *
 * @param indi (output) a vector of integers
 * @param indj (output) a vector of integers
 * @param type is a string composed of the letters 'r' (real) and 'm' (matrix)
 * to specify the type of the extra arguments
 * @param f a function returning an integer (typically a test function)
 * @return OK or FAIL if something went wrong
 */
int FUNCTION(pnl_mat,find) (PnlVectInt *indi, PnlVectInt *indj, char* type, int(*f)(BASE *), ...)
{
  va_list ap;
  TYPE(cell) *args;
  int i, j, k, m, n, count, nvar;
  BASE val, *t;
  m = n = -1;

  nvar = strlen (type);
  if ((args = malloc (sizeof(cell) * nvar)) == NULL) return FAIL;
  if ((t = malloc (sizeof(BASE) * nvar)) == NULL) return FAIL;

  va_start (ap, f);

  for ( k=0; k<nvar ; k++ )
    {
      switch (type[k])
        {
          case 'r' :
            val = va_arg (ap, BASE);
            args[k].x = val; args[k].type = 'r';
            break;
          case 'm' :
            args[k].V= va_arg (ap, TYPE(PnlMat) *);
            args[k].type = 'm';
            if ( m == -1 ) { m = args[k].V->m; n = args[k].V->n; }
            else { PNL_CHECK ( (m != args[k].V->m) || (n != args[k].V->n) ,
                               "incompatible size", "pnl_mat_find"); }
            break;
          default:
            return FAIL;
        }
    }
  va_end(ap);

  /*
   * 2 passes are needed.
   * The first one to determine the size of index
   */
  for ( i=0, count=0 ; i<m ; i++ )
    {
      for ( j=0 ; j<n ; j++ )
        {
          for ( k=0 ; k<nvar ; k++ )
            {
              if ( args[k].type == 'r' ) t[k] = args[k].x;
              else t[k] = PNL_MGET(args[k].V, i, j);
            }
          if ( f(t) == 1 ) { count++; }
        }
    }
  /*
   * Second pass to extract the indices for which f == 1
   */
  if ( indj == NULL )
    {
      pnl_vect_int_resize (indi, count);
      for ( i=0, count=0 ; i<m*n ; i++ )
        {
          for ( k=0 ; k<nvar ; k++ )
            {
              if ( args[k].type == 'r' ) t[k] = args[k].x;
              else t[k] = (args[k].V)->array[i];
            }
          if ( f(t) == 1 ) { PNL_LET(indi, count) = i; count++; }
        }
    }
  else
    {
      pnl_vect_int_resize (indi, count);
      pnl_vect_int_resize (indj, count);

      for ( i=0, count=0 ; i<m ; i++ )
        {
          for ( j=0 ; j<n ; j++ )
            {
              for ( k=0 ; k<nvar ; k++ )
                {
                  if ( args[k].type == 'r' ) t[k] = args[k].x;
                  else t[k] = PNL_MGET(args[k].V, i, j);
                }
              if ( f(t) == 1 )
                {
                  PNL_LET(indi, count) = i; PNL_LET(indj, count) = j;
                  count++;
                }
            }
        }
    }

  free (args);
  free (t);
  return OK;
}

/**
 * Extract a subblock from a matrix
 *
 * @param M_sub on exit M_sub = M(i:i+len_i-1, j:j+len_j-1)
 * @param M a matrix
 * @param i line index
 * @param len_i number of elements to be extracted on the line
 * @param j column index
 * @param len_j number of elements to be extracted on the column
 *
 */
void FUNCTION(pnl_mat,extract_subblock)(TYPE(PnlMat) *M_sub, const TYPE(PnlMat) *M, int i, int len_i, int j, int len_j)
{
  int k, l;
  PNL_CHECK (i+len_i > M->m, "size exceeded", "pnl_mat_extract_subblock");
  PNL_CHECK (j+len_j > M->n, "size exceeded", "pnl_mat_extract_subblock");
  FUNCTION(pnl_mat,resize)(M_sub, len_i, len_j);

  for ( k=0 ; k<len_i ; k++ )
    {
      for ( l=0 ; l<len_j ; l++ )
        {
          PNL_MLET(M_sub, k, l) = PNL_MGET (M, i + k, j + l);
        }
    }
}

/**
 * Set a subblock in a matrix
 *
 * @param M a matrix
 * @param block a matrix
 * @param i line index
 * @param j column index. (i,j) represents the top left hand side anchor
 * inside the matrix M
 *
 */
void FUNCTION(pnl_mat,set_subblock)(TYPE(PnlMat) *M, const TYPE(PnlMat) *block, int i,  int j)
{
  int k, l;
  PNL_CHECK (i+block->m > M->m, "size exceeded", "pnl_mat_set_subblock");
  PNL_CHECK (j+block->n > M->n, "size exceeded", "pnl_mat_set_subblock");

  for ( k=0 ; k<block->m ; k++ )
    {
      for ( l=0 ; l<block->n ; l++ )
        {
          PNL_MLET(M, i+k, j+l) = PNL_MGET (block, k, l);
        }
    }
}

/** 
 * 
 * Add a row before row i and fill it with the content of r
 * 
 * @param M a matrix, modified on output
 * @param i the position at which the row is added
 * @param r a vector to fill the new row with. If NULL is passed, the matrix
 * is enlarged but row i is left uninitialized
 */
void FUNCTION(pnl_mat,add_row)(TYPE(PnlMat) *M, int i, const TYPE(PnlVect) *r)
{
  int new_size;
  PNL_CHECK (i > M->m, "size exceeded", "pnl_mat_add_row");
  PNL_CHECK (i < 0, "size exceeded", "pnl_mat_add_row");
  PNL_CHECK (r->size != M->n, "incompatible size", "pnl_mat_add_row");

  new_size = M->mn + M->n;
  if ( M->mem_size < new_size )
    {
      M->array = realloc (M->array, new_size * sizeof(BASE));
      M->mem_size = new_size;
    }
  /* Move the second block, starting at row i */
  if ( i < M->m ) memmove (M->array + (i+1) * M->n, M->array + i * M->n, (M->m-i) * M->n * sizeof(BASE));
  /* Copy r into row i if r != NULL */
  if ( r != NULL ) memcpy (M->array + i * M->n, r->array, M->n * sizeof(BASE));
  M->m ++;
  M->mn += M->n;
}


/** 
 * 
 * Delete row with index i
 * 
 * @param M a matrix
 * @param i the index of the row to be deleted
 */
void FUNCTION(pnl_mat,del_row)(TYPE(PnlMat) *M, int i)
{
  PNL_CHECK (i >= M->m, "size exceeded", "pnl_mat_del_row");
  PNL_CHECK (i < 0, "size exceeded", "pnl_mat_del_row");
  if ( i < M->m - 1 )
    {
      int rem_len = (M->m - 1 - i) * M->n;
      memmove (M->array + (i * M->n), M->array + ((i+1) * M->n), rem_len * sizeof(BASE));
    }
  M->m --;
  M->mn -= M->n;
}


/**
 * in-place matrix scalar addition
 *
 * @param lhs : left hand side matrix
 * @param  x : scalar
 * @return  lhs = lhs+x
 */
void FUNCTION(pnl_mat,plus_scalar)(TYPE(PnlMat) *lhs, BASE x)
{
  int i;
  for ( i=0 ; i<lhs->mn ; i++ ) { PLUSEQ (lhs->array[i], x); }
}

/**
 * in-place matrix scalar substraction
 *
 * @param lhs : left hand side matrix
 * @param  x : scalar
 * @return  lhs = lhs-x
 */
void FUNCTION(pnl_mat,minus_scalar)(TYPE(PnlMat) *lhs, BASE x)
{
  int i;
  for ( i=0 ; i<lhs->mn ; i++ ) { MINUSEQ (lhs->array[i], x); }
}

/**
 * in-place matrix scalar multiplication
 *
 * @param lhs : left hand side matrix
 * @param x : scalar
 * @return  lhs = lhs*x
 */
void FUNCTION(pnl_mat,mult_scalar)(TYPE(PnlMat) *lhs, BASE x)
{
  int i;
  for ( i=0 ; i<lhs->mn ; i++ ) { MULTEQ (lhs->array[i], x); }
}

/**
 * in-place matrix scalar division
 *
 * @param lhs : left hand side matrix
 * @param x : scalar
 * @return  lhs = lhs*x
 */
void FUNCTION(pnl_mat,div_scalar)(TYPE(PnlMat) *lhs, BASE x)
{
  int i;
  for ( i=0 ; i<lhs->mn ; i++ ) { DIVEQ (lhs->array[i], x); }
}

/**
 * in-place matrix matrix addition
 *
 * @param lhs : left hand side matrix
 * @param rhs : rigth hand side matrix
 * @return  lhs = lhs+rhs
 */
void FUNCTION(pnl_mat,plus_mat)(TYPE(PnlMat) *lhs, const TYPE(PnlMat) *rhs)
{
  int i;
  CheckMatMatch(lhs, rhs);
  for (i=0; i<lhs->mn; i++)
    {
      lhs->array[i] = PLUS (lhs->array[i], rhs->array[i]);
    }
}

/**
 * in-place matrix matrix substraction
 *
 * @param lhs : left hand side matrix
 * @param rhs : rigth hand side matrix
 * @return  lhs = lhs+rhs
 */
void FUNCTION(pnl_mat,minus_mat)(TYPE(PnlMat) *lhs, const TYPE(PnlMat) *rhs)
{
  int i;
  CheckMatMatch(lhs, rhs);
  for (i=0; i<lhs->mn; i++)
    {
      lhs->array[i] = MINUS (lhs->array[i], rhs->array[i]);
    }
}

/**
 * in-place matrix matrix term by term product
 *
 * @param lhs : left hand side matrix
 * @param rhs : rigth hand side matrix
 * @return  lhs = lhs.*rhs
 */
void FUNCTION(pnl_mat,mult_mat_term)(TYPE(PnlMat) *lhs, const TYPE(PnlMat) *rhs)
{
  int i;
  CheckMatMatch(lhs, rhs);
  for (i=0; i<lhs->mn; i++)
    {
      lhs->array[i] = MULT (lhs->array[i], rhs->array[i]);
    }
}

/**
 * in-place matrix matrix term by term division
 *
 * @param lhs : left hand side matrix
 * @param rhs : rigth hand side matrix
 * @return  lhs = lhs.*rhs
 */
void FUNCTION(pnl_mat,div_mat_term)(TYPE(PnlMat) *lhs, const TYPE(PnlMat) *rhs)
{
  int i;
  CheckMatMatch(lhs, rhs);
  for (i=0; i<lhs->mn; i++)
    {
      lhs->array[i] = DIV (lhs->array[i], rhs->array[i]);
    }
}

/**
 * Compute y := alpha * op(A) * x + beta * y
 *     where op(X) = X or op(X) = X'
 *
 * @param A a  matrix
 * @param trans a char. If trans=='T', op(X) = X'. If trans=='N', op(X) = X.
 * @param x a vector such that op(A) * x is legal
 * @param alpha  a BASE constant
 * @param beta  a BASE constant. When b==0, the content of y is not used, instead y
 * is resized to store op(A) * x
 * @param y a vector containing alpha * op(A) * x + beta * y on exit.
 */
void FUNCTION(pnl_mat,dgemv) (char trans, BASE alpha, const TYPE(PnlMat) *A,
                    const TYPE(PnlVect) *x , BASE beta, TYPE(PnlVect) * y)
{
  BASE yi, temp;
  int i,j;

  if ( EQ(beta,ZERO) )
    {
      if (trans=='N' || trans=='n') FUNCTION(pnl_vect,resize) (y, A->m);
      else FUNCTION(pnl_vect,resize)(y, A->n);
      FUNCTION(pnl_vect,set_all) (y, ZERO);
    }
  else if ( NEQ(beta, ONE) ) FUNCTION(pnl_vect,mult_scalar) (y, beta);
  if ( EQ(alpha, ZERO) ) return;


  /* Form alpha * A * x + beta * y */
  if (trans=='N' || trans=='n')
    {
      PNL_CHECK (A->n != x->size || A->m != y->size, "size mismatch", "pnl_mat_dgemv");
      for (i=0; i<y->size; i++)
        {
          temp = ZERO;
          for ( j=0 ; j<A->n ; j++)
            {
              temp = PLUS (temp, MULT(PNL_MGET (A, i, j),
                                      PNL_GET (x, j) ) );
            }
          yi = PNL_GET (y, i);
          yi = PLUS( yi, MULT(alpha, temp));
          PNL_LET (y, i) = yi;
        }
    }
  /* Form alpha * A' * x + beta * y */
  else   if (trans=='T' || trans=='t')
    {
      PNL_CHECK (A->m != x->size || A->n != y->size, "size mismatch", "pnl_mat_dgmev");
      for (j=0; j<A->m; j++)
        {
          temp = MULT(alpha, PNL_GET (x, j));
          for (i=0; i<y->size; i++)
            {
              yi =  PNL_GET (y, i);
              yi = PLUS(yi, MULT(PNL_MGET (A, j, i), temp));
              PNL_LET (y, i) = yi;
            }
        }
    }
  else
    {
      PNL_ERROR ("invalid value for trans", "pnl_mat_dgemv");
    }
}

/**
 *  Computes A * x
 *
 * @param A : matrix
 * @param x : vector
 * @return  y = A * x
 */
TYPE(PnlVect)* FUNCTION(pnl_mat,mult_vect)(const TYPE(PnlMat) *A, const TYPE(PnlVect) *x)
{
  TYPE(PnlVect) *y;
  y = FUNCTION(pnl_vect,create) (0);
  FUNCTION(pnl_mat,dgemv) ('N', ONE, A, x, ZERO, y);
  return y;
}

/**
 *  Computes A' * x
 *
 * @param A : matrix
 * @param x : vector
 * @return  y = A' * x
 */
TYPE(PnlVect)* FUNCTION(pnl_mat,mult_vect_transpose)(const TYPE(PnlMat) *A, const TYPE(PnlVect) *x)
{
  TYPE(PnlVect) *y;
  y = FUNCTION(pnl_vect,create) (A->n);
  FUNCTION(pnl_mat,dgemv) ('T', ONE, A, x, ZERO, y);
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
void FUNCTION(pnl_mat,mult_vect_inplace)(TYPE(PnlVect) *y, const TYPE(PnlMat) *A, const TYPE(PnlVect) *x)
{
  FUNCTION(pnl_mat,dgemv) ('N', ONE, A, x, ZERO, y);
}


/**
 *  Computes y = A' * x
 *
 * @param A : matrix
 * @param x : vector
 * @param y : vector
 * @return  y = A' * x
 */
void FUNCTION(pnl_mat,mult_vect_transpose_inplace)(TYPE(PnlVect) *y, const TYPE(PnlMat) *A, const TYPE(PnlVect) *x)
{
  FUNCTION(pnl_mat,dgemv) ('T', ONE, A, x, ZERO, y);
}

/**
 * Compute y := alpha * A * x + beta * y
 *
 * @param A a  matrix
 * @param x a vector such that A * x is legal
 * @param alpha  a BASE constant
 * @param beta  a BASE constant. When b==0, the content of y is not used, instead y
 * is resized to store A * x
 * @param y a vector containing alpha * A * x + beta * y on exit.
 */
void FUNCTION(pnl_mat,lAxpby)(BASE alpha, const TYPE(PnlMat) *A, const TYPE(PnlVect) *x, BASE beta, TYPE(PnlVect) * y)
{
  FUNCTION(pnl_mat,dgemv) ('N', alpha, A, x, beta, y);
}

/**
 * Compute Y := a * X + Y
 *
 * @param a BASE coefficient
 * @param X a matrix
 * @param Y a matrix whose size must be the same as the one of X. Contains the
 * result on exit
 */
void FUNCTION(pnl_mat,axpy) (BASE a, const TYPE(PnlMat) *X, TYPE(PnlMat) *Y)
{
  int i;
  CheckMatMatch(X, Y);
  if ( EQ(a, ZERO) ) return;
  if ( EQ(a, ONE) )
    {
      for ( i=0 ; i<Y->mn ; i++ )
        {
          Y->array[i] = PLUS(Y->array[i], X->array[i]);
        }
      return;
    }
  for ( i=0 ; i<Y->mn ; i++)
    {
        Y->array[i] = PLUS(Y->array[i], MULT(a, X->array[i]));
    }
}

/**
 * Compute A := alpha x' * y + A
 *
 * @param alpha a BASE number
 * @param x a vector
 * @param y a vector
 * @param A a matrix
 */
void FUNCTION(pnl_mat,dger) (BASE alpha, const TYPE(PnlVect) *x, const TYPE(PnlVect) *y, TYPE(PnlMat) *A)
{
  int i, j;
  PNL_CHECK ((x->size != y->size) || (A->m != A->n) || (A->m != x->size), "size mismatch", "pnl_mat_dger");
  if ( EQ(alpha,ZERO) ) return;

  for ( i=0 ; i<x->size ; i++ )
    {
      const BASE alpha_xi = MULT(alpha, PNL_GET (x, i));
      for (j=0; j<x->size; j++)
        {
          BASE ai = PNL_MGET (A, i, j);
          ai = PLUS(ai, MULT(alpha_xi, PNL_GET (x, j)));
          PNL_MLET(A, i, j) = ai;
        }
    }
}

/**
 * Compute x' A y 

 * @param A : matrix
 * @param x : vector
 * @param y : vector
 * @return x' A y 
 */
BASE FUNCTION(pnl_mat,scalar_prod)(const TYPE(PnlMat) *A, const TYPE(PnlVect) *x , const TYPE(PnlVect) * y)
{
  BASE sum;
  int i,j;
  
#ifndef PNL_RANGE_CHECK_OFF
  if (A->n != y->size || A->m != x->size)
    {
      PNL_ERROR ("size mismatch", "pnl_mat_scalar_prod_A");
    }
#endif

  sum=ZERO;
  for (i=0; i<x->size; i++)
    {
      BASE temp = ZERO;
      for (j=0; j<A->n; j++) 
        {
          temp = PLUS(temp, MULT(PNL_MGET(A,i,j), PNL_GET(y,j)));
        }
      sum = PLUS(sum, MULT(PNL_GET(x,i), temp));
    }
  return sum;
}


/**
 *  matrix multiplication
 *
 * @param rhs1 : first right hand side matrix
 * @param rhs2 : second right hand side matrix
 * @return  rhs1*rhs2
 */

TYPE(PnlMat)* FUNCTION(pnl_mat,mult_mat)(const TYPE(PnlMat) *rhs1, const TYPE(PnlMat) *rhs2)
{
  TYPE(PnlMat) *lhs;
  lhs = FUNCTION(pnl_mat,create) (0,0);
  FUNCTION(pnl_mat,dgemm) ('N', 'N', ONE, rhs1, rhs2, ZERO, lhs);
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

void FUNCTION(pnl_mat,mult_mat_inplace)(TYPE(PnlMat) *lhs, const  TYPE(PnlMat) *rhs1, const TYPE(PnlMat) *rhs2)
{
  FUNCTION(pnl_mat,dgemm) ('N', 'N', ONE, rhs1, rhs2, ZERO, lhs);
}

/**
 * sum matrix componentwise
 *
 * @param lhs : left hand side matrix
 * @return  sum=sum(lhs(i,j))
 */
BASE FUNCTION(pnl_mat,sum)(const TYPE(PnlMat) *lhs)
{
  BASE sum;
  int i;
  sum=ZERO;
  for(i=0;i<lhs->mn;i++) { sum = PLUS(sum, lhs->array[i]); }
  return sum;
}

/**
 * prod matrix componentwise
 *
 * @param lhs : left hand side matrix
 * @return  sum=sum(lhs(i,j))
 */
BASE FUNCTION(pnl_mat,prod)(const TYPE(PnlMat) *lhs)
{
  BASE p;
  int i;
  p=ONE;
  for(i=0;i<lhs->mn;i++) { p = MULT(p, lhs->array[i]); }
  return p;
}

/**
 * sum rows or columns matrix componentwise
 *
 * @param A : left hand side matrix
 * @param y a vector containing on exit y(j) = sum(A(i,j), i=1..m) if c='r'
 * and y(i) = sum(A(i,j), j=1..m) if c='c'
 * @param a :char equals 'r' (if we want to sum rows) we get
 * a row vector  or 'c' (if we want to sum columns) we get a
 * column vector
 */
void FUNCTION(pnl_mat,sum_vect)(TYPE(PnlVect) *y, const TYPE(PnlMat) *A, char a)
{

  BASE sum, yj;
  int i,j;
  if (a=='r')
    {
      /* the loops are done in a non natural order to make the most of
         data alignment while using matrix with row major storage */
      FUNCTION(pnl_vect,resize) (y, A->n);
      FUNCTION(pnl_vect,set_all)(y, ZERO);
      for (i=0; i<A->m; i++)
        {
          for (j=0; j<A->n; j++)
            {
              yj = PNL_GET (y, j);
              yj = PLUS (yj, PNL_MGET (A, i, j));
              PNL_LET (y, j) = yj;
            }
        }
    }
  else if (a=='c')
    {
      FUNCTION(pnl_vect,resize) (y, A->m);
      for (i=0; i<A->m; i++)
        {
          sum = ZERO;
          for (j=0; j<A->n; j++)
            {
              sum = PLUS (sum, PNL_MGET (A, i, j));
            }
          PNL_LET (y, i) = sum;
        }
    }
  else
    {
      PNL_ERROR ("wrong char arg ('c' or 'r')", "pnl_mat_sum_vect");
    }
}

/**
 * prod rows or columns matrix componentwise
 *
 * @param A : left hand side matrix
 * @param y a vector containing on exit y(j) = prod(A(i,j), i=1..m) if c='r'
 * and y(i) = prod(A(i,j), j=1..m) if c='c'
 * @param a :char equals 'r' (if we want to prod rows) we get
 * a row vector  or 'c' (if we want to prod columns) we get a
 * column vector
 */
void FUNCTION(pnl_mat,prod_vect)(TYPE(PnlVect) *y, const TYPE(PnlMat) *A, char a)
{

  BASE prod, yj;
  int i,j;
  if (a=='r')
    {
      /* the loops are done in a non natural order to make the most of
         data alignment while using matrix with row major storage */
      FUNCTION(pnl_vect,resize) (y, A->n);
      FUNCTION(pnl_vect,set_all) (y, ONE);
      for (i=0; i<A->m; i++)
        {
          for (j=0; j<A->n; j++)
            {
              yj = PNL_GET (y, j);
              yj = MULT (yj, PNL_MGET (A, i, j));
              PNL_LET (y, j) = yj;
            }
        }
    }
  else if (a=='c')
    {
      FUNCTION(pnl_vect,resize) (y, A->m);
      for (i=0; i<A->m; i++)
        {
          prod = ONE;
          for (j=0; j<A->n; j++)
            {
              prod = MULT (prod, PNL_MGET (A, i, j));
            }
          PNL_LET (y, i) = prod;
        }
    }
  else
    {
      PNL_ERROR ("wrong char arg ('c' or 'r')", "pnl_mat_prod_vect");
    }
}

/**
 * cumulative sum rows or columns matrix componentwise
 *
 * @param A : left hand side matrix
 * On exit A(i,j) = sum(A(k,j), k=1..i) if c='r'
 * and A(i,j) = sum(A(i,k), k=1..j) if c='c'
 * @param a :char equals 'r' (if we want to sum rows) we get
 * a row vector  or 'c' (if we want to sum columns) we get a
 * column vector
 */
void FUNCTION(pnl_mat,cumsum)(TYPE(PnlMat) *A, char a)
{

  BASE sum, Aij, Ai1j;
  int i, i1, j;
  if (a=='r')
    {
      /* the loops are done in a non natural order to make the most of
         data alignment while using matrix with row major storage */
      for (i=1, i1=0; i<A->m; i++, i1++)
        {
          for (j=0; j<A->n; j++)
            {
              Aij = PNL_MGET (A, i, j);
              Ai1j = PNL_MGET(A, i1, j);
              PNL_MLET(A, i, j) = PLUS(Ai1j, Aij);
            }
        }
    }
  else if (a=='c')
    {
      for (i=0; i<A->m; i++)
        {
          sum = PNL_MGET (A, i, 0);
          for (j=1; j<A->n; j++)
            {
              sum = PLUS(sum, PNL_MGET(A, i, j));
              PNL_MLET (A, i, j) =  sum;
            }
        }
    }
  else
    {
      PNL_ERROR ("wrong char arg ('c' or 'r')", "pnl_mat_cumsum");
    }
}

/**
 * cumulative prod rows or columns matrix componentwise
 *
 * @param A : left hand side matrix
 * On exit A(i,j) = prod(A(k,j), k=1..i) if c='r'
 * and A(i,j) = prod(A(i,k), k=1..j) if c='c'
 * @param a :char equals 'r' (if we want to prod rows) we get
 * a row vector  or 'c' (if we want to prod columns) we get a
 * column vector
 */
void FUNCTION(pnl_mat,cumprod)(TYPE(PnlMat) *A, char a)
{

  BASE prod, Aij, Ai1j;
  int i, i1, j;
  if (a=='r')
    {
      /* the loops are done in a non natural order to make the most of
         data alignment while using matrix with row major storage */
      for (i=1, i1=0; i<A->m; i++, i1++)
        {
          for (j=0; j<A->n; j++)
            {
              Aij = PNL_MGET (A, i, j);
              Ai1j = PNL_MGET(A, i1, j);
              PNL_MLET(A, i, j) = MULT(Ai1j, Aij);
            }
        }
    }
  else if (a=='c')
    {
      for (i=0; i<A->m; i++)
        {
          prod = PNL_MGET (A, i, 0);
          for (j=1; j<A->n; j++)
            {
              prod = MULT(prod, PNL_MGET(A, i, j));
              PNL_MLET (A, i, j) = prod;
            }
        }
    }
  else
    {
      PNL_ERROR ("wrong char arg ('c' or 'r')", "pnl_mat_cumprod");
    }
}

#if defined ORDERED
extern void FUNCTION(pnl_array, min_index)(const BASE *a, int n, int incr,
                                           BASE *min_out, int *imin_out);
extern void FUNCTION(pnl_array, max_index)(const BASE *a, int n, int incr,
                                           BASE *max_out, int *imax_out);
extern void FUNCTION(pnl_array, minmax_index)(const BASE *a, int n, int incr,
                                              BASE *min_out, BASE *max_out,
                                              int *imin_out, int *imax_out);

/**
 * Return the maxima of the components of a matrix row-wise or column-wise
 *
 * @param A a matrix
 * @param d can be 'c' (out(i) = max(A(i,:)) or 'r' (out(i) = max(A(:,i))
 * @param out a vector containing on exit the maxima
 * @param index a vector of integers containing on exit the indices of the
 * maxima. if NULL, no index is computed.
 */
void FUNCTION(pnl_mat, max_index)(TYPE(PnlVect) *out, PnlVectInt *index, const TYPE(PnlMat) *A, char d)
{
  int length, n, incr, i, lda, fake;
  switch (d)
    {
    case 'c' : n = A->m; length = A->n; lda = A->n;  incr = 1; break;
    case 'r' : length = A->m; n = A->n, lda = 1; incr = A->n; break;
    case '*' : length = A->mn; n = 1, lda = 1; incr = 1; break;
    default : PNL_ERROR("unknow direction", "pnl_mat_max_index"); break;
    }
  if (index!=NULL) pnl_vect_int_resize (index, n);
  FUNCTION(pnl_vect, resize)(out, n);

  for (i=0; i<n; i++)
    {
      FUNCTION(pnl_array, max_index)(&(A->array[i*lda]), length, incr, &(out->array[i]), &fake);
      if (index!=NULL) index->array[i] = fake;
    }
}

/**
 * Return the minima of the components of a matrix row-wise or column-wise
 *
 * @param A a matrix
 * @param d can be 'c' (out(i) = min(A(i,:)) or 'r' (out(i) = min(A(:,i))
 * @param out a vector containing on exit the minima
 * @param index a vector of integers containing on exit the indices of the
 * minima. if NULL, no index is computed.
 */
void FUNCTION(pnl_mat, min_index)(TYPE(PnlVect) *out, PnlVectInt *index, const TYPE(PnlMat) *A,  char d)
{
  int length, n, incr, i, fake, lda;
  switch (d)
    {
    case 'c' : n = A->m; length = A->n; lda = A->n; incr = 1; break;
    case 'r' : length = A->m; n = A->n; lda = 1; incr = A->n; break;
    case '*' : length = A->mn; n = 1, lda = 1; incr = 1; break;
    default : PNL_ERROR("unknow direction", "pnl_mat_min_index"); break;
    }
  if (index!=NULL) pnl_vect_int_resize (index, n);
  FUNCTION(pnl_vect, resize)(out, n);

  for (i=0; i<n; i++)
    {
      FUNCTION(pnl_array, min_index)(&(A->array[i*lda]), length, incr, &(out->array[i]), &fake);
      if (index!=NULL) index->array[i] = fake;
    }
}

/**
 * Return the minima of the components of a matrix row-wise or column-wise
 *
 * @param A a matrix
 * @param d can be 'c' (out(i) = min(A(i,:)) or 'r' (out(i) = min(A(:,i))
 * @param out_min a vector containing on exit the minima
 * @param out_max a vector containing on exit the maxima
 * @param index_min a vector of integers containing on exit the indices of the
 * minima. if NULL, no index is computed.
 * @param index_max a vector of integers containing on exit the indices of the
 * maxima. if NULL, no index is computed.
 */
void FUNCTION(pnl_mat, minmax_index)(TYPE(PnlVect) *out_min, TYPE(PnlVect) *out_max, PnlVectInt *index_min, PnlVectInt *index_max, const TYPE(PnlMat) *A, char d)
{
  int length, n, incr, i, lda, fake_min, fake_max;
  switch (d)
    {
    case 'c' : n = A->m; length = A->n; lda = A->n; incr = 1; break;
    case 'r' : length = A->m; n = A->n; lda = 1; incr = A->n; break;
    case '*' : length = A->mn; n = 1, lda = 1; incr = 1; break;
    default : PNL_ERROR("unknow direction", "pnl_mat_minmax_index"); break;
    }
  if (index_min != NULL) pnl_vect_int_resize (index_min, n);
  if (index_max != NULL) pnl_vect_int_resize (index_max, n);
  FUNCTION(pnl_vect, resize)(out_min, n);
  FUNCTION(pnl_vect, resize)(out_max, n);

  for (i=0; i<n; i++)
    {
      FUNCTION(pnl_array, minmax_index)(&(A->array[i*lda]), length, incr, &(out_min->array[i]), &(out_max->array[i]), &fake_min, &fake_max);
      if (index_min != NULL) index_min->array[i] = fake_min;
      if (index_max != NULL) index_max->array[i] = fake_max;
    }
}

/**
 * Return the maxima of the components of a matrix row-wise or column-wise
 *
 * @param A a matrix
 * @param d can be 'c' (out(i) = max(A(i,:)) or 'r' (out(i) = max(A(:,i))
 * @param out a vector containing on exit the maxima
 */
void FUNCTION(pnl_mat,max)(TYPE(PnlVect) *out, const TYPE(PnlMat) *A,  char d)
{
  FUNCTION (pnl_mat, max_index)(out, NULL, A, d);
}

/**
 * Return the minima of the components of a matrix row-wise or column-wise
 *
 * @param A a matrix
 * @param d can be 'c' (out(i) = min(A(i,:)) or 'r' (out(i) = min(A(:,i))
 * @param out a vector containing on exit the minima
 */
void FUNCTION(pnl_mat,min)(TYPE(PnlVect) *out, const TYPE(PnlMat) *A, char d)
{
  FUNCTION (pnl_mat, min_index)(out, NULL, A, d);
}

/**
 * Return the minima and maxima of the components of a matrix row-wise or column-wise
 *
 * @param A a matrix
 * @param d can be 'c' (min_out(i) = min(A(i,:)) or 'r' (min_out(i) = min(A(:,i))
 * @param min_out a vector containing on exit the minima
 * @param max_out a vector containing on exit the maxima
 */
void FUNCTION(pnl_mat,minmax)( TYPE(PnlVect) *min_out,
                              TYPE(PnlVect) *max_out, const TYPE(PnlMat) *A, char d)
{
  FUNCTION (pnl_mat, minmax_index)(min_out, max_out, NULL, NULL, A, d);
}

static int FUNCTION(__pnl, cmp_i) ( const void *a, const void *b)
{
  if ( * (BASE *) a < * (BASE *) b) return -1;
  if ( * (BASE *) a == * (BASE *) b) return 0;
  return 1;
}

static int FUNCTION(__pnl, cmp_d) ( const void *a, const void *b)
{
  if ( * (BASE *) a > * (BASE *) b) return -1;
  if ( * (BASE *) a == * (BASE *) b) return 0;
  return 1;
}

static void FUNCTION(pnl_mat, qsort_aux)(TYPE(PnlMat) * A, PnlMatInt *t, int use_index, char dir, char order)
{
  int i, j, lda, loops, n, incr, *index = NULL;
  cmp_func f;

  lda = incr = loops = n = 0;   /* avoid warnings */

  if (use_index == TRUE) { pnl_mat_int_resize (t, A->m, A->n); index=t->array; }

  switch (order)
    {
    case 'i' : f = FUNCTION(__pnl, cmp_i); break;
    case 'd' : f = FUNCTION(__pnl, cmp_d); break;
    default : PNL_ERROR ("unknow order", "pnl_mat_qsort"); break;
    }

  switch (dir)
    {
    case 'r' : lda = A->n; incr = 1; n = A->m; loops = A->n;
      if (use_index == TRUE){
        for (i=0; i<t->m; i++)
          for (j=0; j<t->n; j++) { pnl_mat_int_set (t, i, j, i); }
      }
      break;
    case 'c' : lda = 1; incr = A->n; n = A->n; loops = A->m;
      if (use_index == TRUE){
        for (i=0; i<t->m; i++)
          for (j=0; j<t->n; j++) { pnl_mat_int_set (t, i, j, j); }
      }
      break;
    default: PNL_ERROR("unknown direction", "pnl_mat_qsort_aux"); break;
    }
  for (i=0; i< loops; i++)
    pnl_qsort ( A->array + i * incr, n, sizeof(BASE), lda, index + i * incr, lda, use_index, f);
}

/**
 * Quick sort function for matrices
 *
 * @param A a TYPE(PnlMat), contains the sorted matrix on exit
 * @param order can be 'i' or 'd' for increasing or decreasing order.
 * @param dir can be 'r' or 'c'
 */
void FUNCTION(pnl_mat, qsort)(TYPE(PnlMat) * A, char dir, char order)
{
  FUNCTION(pnl_mat, qsort_aux)(A, NULL, FALSE, dir, order);
}

/**
 * Quick sort function for matrices with permutation index
 *
 * @param A a TYPE(PnlMat), contains the sorted matrix on exit
 * @param t a PnlMatInt, on exit contains the permutation used to sort A
 * @param order can be 'i' or 'd' for increasing or decreasing order.
 * @param dir can be 'r' or 'c'
 */
void FUNCTION(pnl_mat, qsort_index)(TYPE(PnlMat) * A, PnlMatInt *t, char dir, char order)
{
  FUNCTION(pnl_mat, qsort_aux)(A, t, TRUE, dir, order);
}

#endif


/****************************
 *** PnlHmat functions ***
 ****************************/
extern void pnl_hmat_compute_pdims (int *pdims, int ndim, const int *dims);
extern int pnl_hmat_compute_linear_index (PnlHmatObject *H, int *tab);

TYPE(PnlHmat)* FUNCTION(pnl_hmat,init)(TYPE(PnlHmat) *o)
{
  o->object.parent_type = PNL_TYPE_HMATRIX;
  o->object.type = CONCAT2(PNL_TYPE_HMATRIX_, BASE_TYPE);
  o->object.label = FUNCTION(pnl_hmat,label);
  o->object.destroy = (DestroyFunc *) pnl_hmat_object_free;
  o->object.constructor = (NewFunc *) FUNCTION(pnl_hmat,new);
  o->object.clone = (CloneFunc *) FUNCTION(pnl_hmat,clone);
  o->object.copy = (CopyFunc *) FUNCTION(pnl_hmat,copy);
  o->ndim = 0;
  o->mn = 0;
  o->dims = NULL;
  o->pdims = NULL;
  o->array = NULL;
  return o;
}

TYPE(PnlHmat)* FUNCTION(pnl_hmat,new)()
{
  TYPE(PnlHmat) *o;
  if ( (o = (TYPE(PnlHmat) *) pnl_hmat_object_new ()) == NULL) return NULL;
  FUNCTION(pnl_hmat,init)(o);
  return o;
}

TYPE(PnlHmat)* FUNCTION(pnl_hmat,create)(int ndim, const int *dims)
{
  TYPE(PnlHmat) *H;
  if ((H=FUNCTION(pnl_hmat,new)())==NULL) return NULL;
  H->ndim = ndim;
  if (ndim>0)
    {
      if ((H->dims=malloc(sizeof(int)*ndim))==NULL) return NULL;
      if ((H->pdims=malloc(sizeof(int)*ndim))==NULL) return NULL;
      memcpy(H->dims, dims, ndim*sizeof(int));
      pnl_hmat_compute_pdims (H->pdims, ndim, dims);
      H->mn = H->dims[0] * H->pdims[0];
      if((H->array=malloc(H->mn*sizeof(BASE)))==NULL) return NULL;
    }
  else
    {
      H->dims = NULL;
      H->pdims = NULL;
      H->mn = 0;
      H->array = (BASE*)NULL;
    }
  return H;
}

TYPE(PnlHmat)* FUNCTION(pnl_hmat,create_from_scalar)(int ndim, const int *dims, BASE x)
{
  TYPE(PnlHmat) *H;
  int i;

  if ((H=FUNCTION(pnl_hmat,create)(ndim,dims))==NULL) return NULL;

  for ( i=0 ; i<H->mn ; i++ ) { H->array[i] = x; }
  return H;
}

TYPE(PnlHmat)* FUNCTION(pnl_hmat,create_from_ptr)(int ndim, const int *dims, const BASE *x)
{
  TYPE(PnlHmat) *H;

  if ((H=FUNCTION(pnl_hmat,create)(ndim,dims))==NULL) return NULL;

  memcpy(H->array, x, H->mn*sizeof(BASE));
  return H;
}

void FUNCTION(pnl_hmat,free)(TYPE(PnlHmat) **H)
{
  PnlHmatObject *o;
  o = PNL_HMAT_OBJECT (*H);
  pnl_hmat_object_free (&o);
  *H = NULL;
}

/**
 * Copy a TYPE(PnlHmat)
 *
 * @param H : a constant TYPE(PnlHmat) pointer
 * @return a TYPE(PnlHmat) pointer initialised with H
 */
TYPE(PnlHmat)* FUNCTION(pnl_hmat,copy)(const TYPE(PnlHmat) *H)
{
  TYPE(PnlHmat) *ret;
  if ((ret=FUNCTION(pnl_hmat,create)(H->ndim, H->dims))==NULL) return NULL;
  memcpy(ret->array, H->array, sizeof(BASE)*ret->mn);
  return ret;
}

/**
 * Clone a TYPE(PnlHmat)
 *
 * @param H : a constant TYPE(PnlHmat) pointer
 * @param clone : needs to be already allocated. NOT TESTED
 */
void FUNCTION(pnl_hmat,clone)(TYPE(PnlHmat) *clone, const TYPE(PnlHmat) *H)
{
  FUNCTION(pnl_hmat,resize)(clone, H->ndim,H->dims);
  memcpy(clone->array, H->array, sizeof(BASE)*H->mn);
}

/**
 * Resize a TYPE(PnlHmat).
 *
 * If the new size is smaller than the current one, no
 * memory is free. If the new size is larger than the
 * current one, more space is allocated. Note that for the
 * sake of efficiency the old data are not copied.
 *
 * @param H : a pointer to an already existing TYPE(PnlHmat) (mn
 * must be initialised)
 * @param ndim : new nb dimensions
 * @param dims : new pointer to the dimensions array
 *
 * @return OK or FAIL. When returns OK, the hmatrix is changed.
 */
int FUNCTION(pnl_hmat,resize)(TYPE(PnlHmat) *H, int ndim, const int *dims)
{
  return pnl_hmat_object_resize (PNL_HMAT_OBJECT(H), ndim, dims);
}


/**
 * Set the value of self[tab]=x
 *
 * @param self : a TYPE(PnlHmat)
 * @param tab : coordinates array
 * @param x : self[tab]=x
 */
void FUNCTION(pnl_hmat,set)(TYPE(PnlHmat) *self, int *tab, BASE x)
{
  int index;
  CheckIndexHmat(self,tab);

  index = pnl_hmat_compute_linear_index (PNL_HMAT_OBJECT(self), tab);
  self->array[index]=x;
}

/**
 * Get the value of self[tab]
 *
 * @param self : a TYPE(PnlHmat)
 * @param tab : coordinates array;
 * @return  self[tab]
 */
BASE FUNCTION(pnl_hmat,get)(const TYPE(PnlHmat) *self, int *tab)
{
  int index;
  CheckIndexHmat(self,tab);

  index = pnl_hmat_compute_linear_index (PNL_HMAT_OBJECT(self), tab);
  return self->array[index];
}

/**
 * Return the address of self[tab] for use as a lvalue.
 *
 * @param self : a TYPE(PnlHmat)
 * @param tab : coordinates array
 * @return  &(self[i,j])
 */
BASE* FUNCTION(pnl_hmat,lget)(TYPE(PnlHmat) *self, int *tab)
{
  int index;
  CheckIndexHmat(self,tab);

  index = pnl_hmat_compute_linear_index (PNL_HMAT_OBJECT(self), tab);
  return &(self->array[index]);
}

/**
 * Print a hmatrix
 *
 * @param H : a TYPE(PnlHmat) pointer.
 */
void FUNCTION(pnl_hmat,print)(const TYPE(PnlHmat) *H)
{
  BASE *lptr=H->array;
  int *ptr=NULL;
  int i,j,s;
  int nd;
  int nd_1;
  int *index=NULL;
  int d=H->ndim;
  int l=d-3;
  s=0;
  ptr=H->dims+(d-1);/*last cell of H->dims */
  nd=*(ptr);
  nd_1=*(ptr-1);
  if ((index=malloc((d-2)*sizeof(int)))==NULL) abort();
  for (i=0; i<d-2; i++) index[i]=0;
  while (s<H->mn)
    {
      printf("cell=[ ");
      for ( i=0 ; i<d-2 ; i++) printf("%d ",index[i]);
      printf("]\n");
      for ( i=0 ; i<nd_1 ; i++)
        {
          for (j=0 ; j<nd ;j++)
            {
              printf(OUT_FORMAT ,OUT_PUT_FORMAT(*lptr)); printf (" ");
              lptr++;s++;
            }
          printf("\n");
        }
      if(l>=0)
        {
          if(index[l]<H->dims[l]-1)
            (index[l])++;
          else
            {
              while(l>=0 && index[l]==H->dims[l]-1)
                {index[l]=0; l--;}
              if (l==-1 && s != H->mn) abort();
              if (l>=0)
                {
                  (index[l])++; l=d-3;
                }
            }
        }
    }
  free(index);
}

/**
 * in-place hypermatrix hypermatrix addition
 *
 * @param lhs : left hand side hypermatrix
 * @param rhs : rigth hand side hypermatrix
 * @return  lhs = lhs+rhs
 */
void FUNCTION(pnl_hmat,plus_hmat)(TYPE(PnlHmat) *lhs, const TYPE(PnlHmat) *rhs)
{

  int i;
  CheckHmatMatch(lhs, rhs);
  for ( i=0 ; i<lhs->mn ; i++ )
    {
      lhs->array[i] = PLUS(lhs->array[i], rhs->array[i]);
    }
}

/**
 * in-place hyper matrix scalar multiplication
 *
 * @param lhs : left hand side hyper matrix
 * @param x : scalar
 * @return  lhs = lhs*x
 */
void FUNCTION(pnl_hmat,mult_scalar)(TYPE(PnlHmat) *lhs, BASE x)
{
  int i;
  for ( i=0 ; i<lhs->mn ; i++ )
    {
      lhs->array[i] = MULT(lhs->array[i], x);
    }
}

/** 
 * 
 * Create a PnlMat wrapper to the block pointed by t. t must be of size
 * M->ndim - 2
 * 
 * @param M
 * @param t an array of size M->ndim - 2, this is not checked in the
 * function
 */
TYPE(PnlMat) FUNCTION(pnl_mat,wrap_hmat)(TYPE(PnlHmat) *M, int *t)
{
  int i, leading_index;
  TYPE(PnlMat) mat;
  leading_index = 0;
  for ( i=0 ; i<M->ndim-2 ; i++ ) { leading_index += t[i] * M->pdims[i]; }
  mat = FUNCTION(pnl_mat,wrap_array) (M->array + leading_index, M->dims[M->ndim-2], M->dims[M->ndim-1]);
  return mat;
}

/** 
 * 
 * Create a PnlVect wrapper to the block pointed by t. t must be of size
 * M->ndim - 1
 * 
 * @param M
 * @param t an array of size M->ndim - 1, this is not checked in the
 * function
 */
TYPE(PnlVect) FUNCTION(pnl_vect,wrap_hmat)(TYPE(PnlHmat) *M, int *t)
{
  int i, leading_index;
  TYPE(PnlVect) v;
  leading_index = 0;
  for ( i=0 ; i<M->ndim-1 ; i++ ) { leading_index += t[i] * M->pdims[i]; }
  v = FUNCTION(pnl_vect,wrap_array) (M->array + leading_index, M->dims[M->ndim-1]);
  return v;
}


#ifdef PNL_CLANG_COMPLETE
#include "pnl/pnl_templates_off.h"
#undef  BASE_DOUBLE
#endif
