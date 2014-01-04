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
#include "pnl/pnl_sp_matrix.h"
#include "pnl/pnl_internals.h"
#define BASE_DOUBLE
#include "pnl/pnl_templates_on.h"
#endif

/**
 * Initialize a PnlSpMat
 */
void FUNCTION(pnl_sp_mat,init)(TYPE(PnlSpMat) *o)
{
  o->object.parent_type = PNL_TYPE_SP_MATRIX;
  o->object.type = CONCAT2(PNL_TYPE_SP_MATRIX_, BASE_TYPE);
  o->object.label = FUNCTION(pnl_sp_mat,label);
  o->object.destroy = (DestroyFunc *) pnl_sp_mat_object_free;
  o->object.constructor = (NewFunc *) FUNCTION(pnl_sp_mat,new);
  o->object.clone = (CloneFunc *) FUNCTION(pnl_sp_mat,clone);
  o->object.copy = (CopyFunc *) FUNCTION(pnl_sp_mat,copy);
  o->m = 0;
  o->n = 0;
  o->nz = 0;
  o->nzmax = 0;
  o->array = NULL;
  o->I = NULL;
  o->J = NULL;
}

/**
 * Create an empty PnlSpMat
 * @return a TYPE(PnlSpMat) pointer
 */
TYPE(PnlSpMat) * FUNCTION(pnl_sp_mat,new)()
{
  TYPE(PnlSpMat) *o;
  if ( (o = (TYPE(PnlSpMat) *) pnl_sp_mat_object_new ()) == NULL) return NULL;
  FUNCTION(pnl_sp_mat,init)(o);
  return o;
}

/**
 * Free a TYPE(PnlSpMat)
 *
 * @param v adress of a TYPE(PnlSpMat)*. v is set to NULL at exit.
 */
void FUNCTION(pnl_sp_mat,free)(TYPE(PnlSpMat) **v)
{
  PnlSpMatObject *o;
  o = PNL_SP_MAT_OBJECT(*v);
  pnl_sp_mat_object_free (&o);
  *v = NULL;
}

/**
 * Copy a TYPE(PnlSpMat)
 *
 * @param v : a constant TYPE(PnlSpMat) pointer
 * @return a TYPE(PnlSpMat) pointer initialised with v
 */
TYPE(PnlSpMat)* FUNCTION(pnl_sp_mat,copy)(const TYPE(PnlSpMat) *v)
{
  TYPE(PnlSpMat) *ret;

  if ((ret=FUNCTION(pnl_sp_mat,new)())==NULL) return NULL;
  FUNCTION(pnl_sp_mat,clone)(ret, v);
  return ret;
}

/**
 * Clone a TYPE(PnlSpMat)
 *
 * @param M : a constant TYPE(PnlSpMat) pointer
 * @param clone : needs to be already allocated. NOT TESTED
 */
void FUNCTION(pnl_sp_mat,clone)(TYPE(PnlSpMat) *clone, const TYPE(PnlSpMat) *M)
{
  FUNCTION(pnl_sp_mat,resize)(clone, M->m,M->n, M->nz);
  memcpy(clone->I, M->I, sizeof(int) * (M->m + 1));
  memcpy(clone->J, M->J, sizeof(int) * M->nz);
  memcpy(clone->array, M->array, sizeof(BASE) * M->nz);
}

/**
 * Resize a TYPE(PnlSpMat).  If the new size is smaller than the current one, no
 * memory is freed. If the new size is larger than the current nzmax, a new
 * pointer is allocated. The old data are kept only when (m,n) are unchanged
 * and we increase nzmax.
 *
 * @param M a pointer to an already existing TYPE(PnlSpMat)
 * @param m new nb of rows
 * @param n new nb of columns
 * @param nzmax new maximum number of non-zero elements.
 *
 * @return OK or FAIL. When returns OK, the matrix M is changed.
 */
int FUNCTION(pnl_sp_mat,resize)(TYPE(PnlSpMat) *M, int m, int n, int nzmax)
{
  return pnl_sp_mat_object_resize (PNL_SP_MAT_OBJECT(M), m, n, nzmax);
}

/** 
 * Create a sparse matrix of size m x n with all entries set to 0.
 * Note that the created matrix can hold as many as nz entries.
 * 
 * @param m number of rows 
 * @param n number of cols
 * @param nzmax maximum number of non zeros elements
 */
TYPE(PnlSpMat)* FUNCTION(pnl_sp_mat,create)(int m, int n, int nzmax)
{
  int i;
  TYPE(PnlSpMat) *M = FUNCTION(pnl_sp_mat,new)();  
  M->I = malloc (sizeof(int) * (m + 1));
  for ( i=0 ; i<m+1 ; i++ ) M->I[i] = 0;
  M->J = malloc (sizeof(int) * nzmax);
  for ( i=0 ; i<nzmax ; i++ ) M->J[i] = 0;
  M->array = malloc (sizeof(BASE) * nzmax);
  M->m = m; M->n = n; M->nz = 0; 
  M->nzmax = nzmax;
  return M;
}

/** 
* 
* Set M(,i, j) = x. May create this entry or overwrite it.
* 
* @param M sparse matrix
* @param i row index
* @param j col index
* @param x value
*/
void FUNCTION(pnl_sp_mat,set)(TYPE(PnlSpMat) *M, int i, int j, BASE x)
{
  int k, k0, k1, k2;
  CheckIndexSpMat(M, i, j);

  k1 = M->I[i];
  k2 = M->I[i+1];
  for ( k=k1 ; k<k2 ; k++ ) 
    { 
      if ( M->J[k] == j ) { M->array[k] = x; return; }
      if ( M->J[k] > j ) { break; }
    }

  /* If we arrive here, it means that the index (i,j) is not already in M, we
   * need to add it */
  if ( M->nz == M->nzmax ) FUNCTION(pnl_sp_mat,resize)(M, M->m, M->n, M->nz + 1);
  M->nz++;
  k0 = k;
  for ( k=M->nz ; k>k0 ; k-- )
    {
      M->array[k] = M->array[k-1]; M->J[k] = M->J[k-1];
    }
  for ( k=i+1 ; k<M->m+1 ; k++ ) { M->I[k]++; }
  M->array[k0] = x; M->J[k0] = j; 
}

/** 
* 
* Get M(,i, j) in x. Return TRUE if M has an entry (i,j) and FALSE otherwise
* 
* @param M sparse matrix
* @param i row index
* @param j col index
*/
BASE FUNCTION(pnl_sp_mat,get)(const TYPE(PnlSpMat) *M, int i, int j)
{
  int k, k1, k2;
  CheckIndexSpMat(M, i, j);

  k1 = M->I[i];
  k2 = M->I[i+1];
  for ( k=k1 ; k<k2 ; k++ ) 
    { 
      if ( M->J[k] == j ) { return M->array[k]; }
      if ( M->J[k] > j ) { break; }
    }
  return ZERO;
}

/** 
 * Create a dense matrix from a sparse one
 * 
 * @param Sp a sparse matrix
 */
TYPE(PnlMat)* FUNCTION(pnl_mat,create_from_sp_mat)(const TYPE(PnlSpMat) *Sp)
{
  int i, j;
  TYPE(PnlMat) *M = FUNCTION(pnl_mat,create) (Sp->m, Sp->n);
  FUNCTION(pnl_mat,set_zero)(M);

  for ( i=0 ; i<M->m ; i++ )
    {
      int k1 = Sp->I[i];
      int k2 = Sp->I[i+1];
      for ( j=k1 ; j<k2 ; j++ )
        {
          PNL_MLET(M, i, Sp->J[j]) = Sp->array[j];
        }
    }
  return M;
}

/** 
 * Create a sparse matrix from a dense one
 * 
 * @param M a dense matrix
 */
TYPE(PnlSpMat)* FUNCTION(pnl_sp_mat,create_from_mat)(const TYPE(PnlMat) *M)
{
  int i, j, nz;
  TYPE(PnlSpMat) *Sp;

  for ( i=0, nz=0 ; i<M->mn ; i++ ) { if ( NEQ(M->array[i], ZERO) ) nz++; }
  Sp = FUNCTION(pnl_sp_mat,create)(M->m, M->n, nz);

  for ( i=0 ; i<M->m ; i++ ) 
    {
      Sp->I[i] = Sp->nz;
      for ( j=0 ; j<M->n ; j++ )
        {
          BASE Mij = PNL_MGET(M, i, j);
          if ( NEQ(Mij, ZERO) ) 
            {
              Sp->array[Sp->nz] = Mij;
              Sp->J[Sp->nz] = j;
              Sp->nz++;
            }
        }
    }
  Sp->I[M->m] = Sp->nz;
  return Sp;
}

/** 
 * Compare two sparse matrices
 * 
 * @param Sp1  a sparse matrix
 * @param Sp2  a sparse matrix
 * 
 * @return  TRUE or FALSE
 */
int FUNCTION(pnl_sp_mat, eq)(const TYPE(PnlSpMat) *Sp1, const TYPE(PnlSpMat) *Sp2)
{
  int k;
  if ( (Sp1->m != Sp2->m) || (Sp1->n != Sp1->n) || (Sp1->nz != Sp2->nz) ) return FALSE;
  for ( k=0 ; k<Sp1->m ; k++ ) { if ( Sp1->I[k] != Sp2->I[k] ) return FALSE; }
  for ( k=0 ; k<Sp1->nz ; k++ ) 
    {
      if ( (Sp1->J[k] != Sp2->J[k]) || NEQ(Sp1->array[k], Sp2->array[k]) ) return FALSE; 
    }
  return TRUE;
}

#ifdef PNL_CLANG_COMPLETE
#include "pnl/pnl_templates_off.h"
#undef  BASE_DOUBLE
#endif


