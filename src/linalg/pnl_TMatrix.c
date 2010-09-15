
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
#include <math.h>
#include <ctype.h>
#include <stdarg.h>

#include "config.h"
#include "pnl_mathtools.h"
#include "pnl_matrix.h"
#include "pnl_internals.h"

static char pnl_matrix_label[] = "PnlMatObject";
static char pnl_hmatrix_label[] = "PnlHmatObject";

/**
 * Creates a PnlMatObject which is the parent type of all matrices
 */
PnlMatObject* pnl_mat_object_new ()
{
  PnlMatObject *o;
  if ( (o = malloc (sizeof (PnlMatObject))) == NULL) return NULL;
  o->m = 0;
  o->n = 0;
  o->mn = 0;
  o->mem_size = 0;
  o->owner = 1;
  o->array = NULL;
  o->object.type = PNL_TYPE_MATRIX;
  o->object.parent_type = PNL_TYPE_MATRIX;
  o->object.label = pnl_matrix_label;
  o->object.destroy = (destroy_func *) pnl_mat_object_free;
  return o;
}

/**
 * Frees a PnlMatObject
 */
void pnl_mat_object_free (PnlMatObject **o)
{
    if (*o != NULL)
    {
      if ((*o)->array != NULL && (*o)->owner == 1)
        { free((*o)->array); (*o)->array = NULL; }
      free(*o);
    }
}

/**
 * Resizes a PnlMatObject.  If the new size is smaller than the current one, no
 * memory is freed. If the new size is larger than the current mem_size, a new
 * pointer is allocated. The old data are kept.
 *
 * @param M : a pointer to an already existing PnlMatObject 
 * @param m : new nb rows
 * @param n : new nb columns
 *
 * @return OK or FAIL. When returns OK, the matrix is changed. 
 */
int pnl_mat_object_resize(PnlMatObject *M, int m, int n)
{
  int mn;
  size_t sizeof_base = 0;
  mn = m*n;
  if (M->owner == 0) return OK;
  if (mn < 0) return FAIL;
  if (mn == 0) /* free array */
    {
      M->m = M->n = M->mn = M->mem_size = 0;
      if (M->array != NULL) { free(M->array); M->array = NULL; }
      return OK;
    }

  if (M->mem_size >= mn) 
    {
      /* If the new size is smaller, we do not reduce the size of the
         allocated block. It may change, but it may allow to grow the matrix
         quicker */
      M->m=m; M->n=n; M->mn = mn;
      return OK;
    }

  /* Now, M->mem_size < mn */
  M->m = m; M->n = n;
  M->mn = M->mem_size = mn;
  switch (PNL_GET_TYPE(M))
    {
    case PNL_TYPE_MATRIX_DOUBLE : sizeof_base = sizeof(double);
      break;
    case PNL_TYPE_MATRIX_COMPLEX : sizeof_base = sizeof(dcomplex);
      break;
    case PNL_TYPE_MATRIX_INT : sizeof_base = sizeof(int);
      break;
    default : PNL_MESSAGE (1, "Unknown type in pnl_mat_object_resize.\n");
      return FAIL;
    }
  if ((M->array=realloc(M->array,mn*sizeof_base)) == NULL) return FAIL;
  return OK;
}


/**
 * Creates a PnlHmatObject which is the parent type of all hyper-matrices
 */
PnlHmatObject* pnl_hmat_object_new ()
{
  PnlHmatObject *o;
  if ( (o = malloc (sizeof (PnlHmatObject))) == NULL) return NULL;
  o->ndim = 0;
  o->dims = NULL;
  o->mn = 0;
  o->array = 0;
  o->object.type = PNL_TYPE_HMATRIX;
  o->object.parent_type = PNL_TYPE_HMATRIX;
  o->object.label = pnl_hmatrix_label;
  return o;
}

/**
 * Frees a Hmat object
 */
void pnl_hmat_object_free(PnlHmatObject **H)
{
  if (*H != NULL)
    {
      free((*H)->array);
      free((*H)->dims);
      free(*H);
      *H=NULL;
    }
}


/**
 * resizes a PnlHmatObject.
 *
 * If the new size is smaller than the current one, no
 * memory is free. If the new size is larger than the
 * current one, more space is allocated. Note that for the
 * sake of efficiency the old data are not copied.
 *
 * @param H : a pointer to an already existing PnlHmatObject
 * @param ndim : new nb dimensions
 * @param dims : new pointer to the dimensions array
 *
 * @return OK or FAIL. When returns OK, the hmatrix is changed. 
 */
int pnl_hmat_object_resize(PnlHmatObject *H, int ndim, const int *dims)
{
  int i, s=1;
  size_t sizeof_base = 0;
  const int *ptr;
  ptr=dims;
  for(i=0;i<ndim;i++) { s*=(*ptr); ptr++; }

  if (H->mn == s) /*nothing to do, just adjust ndim and dims*/
    {
      H->ndim=ndim;
      if(H->ndim> ndim) if ((H->dims=realloc(H->dims,sizeof(int)*ndim))==NULL) return FAIL;
      memcpy(H->dims, dims, ndim*sizeof(int));
      return OK;
    }
  if (s< 0) return FAIL;
  if (s==0) /* free array */
    {
      H->ndim =  H->mn = 0;
      H->dims=NULL;
      free(H->array); H->array = NULL;
      return OK;
    }
  H->ndim=ndim; H->mn=s;
  if ((H->dims=realloc(H->dims,sizeof(int)*ndim))==NULL) return FAIL;
  memcpy(H->dims, dims, ndim*sizeof(int));

  switch (PNL_GET_TYPE (H))
    {
    case PNL_TYPE_HMATRIX_DOUBLE : sizeof_base = sizeof(double);
      break;
    case PNL_TYPE_HMATRIX_COMPLEX : sizeof_base = sizeof(dcomplex);
      break;
    case PNL_TYPE_HMATRIX_INT : sizeof_base = sizeof(int);
      break;
    default : PNL_MESSAGE (1, "Unknown type in pnl_hmat_object_resize.\n");
      return FAIL;
    }
  if (H->array==NULL)
    {
      if ((H->array = malloc(H->mn*sizeof_base))==NULL)
        return FAIL;
    }else
    {
      if ((H->array = realloc(H->array,H->mn*sizeof_base))==NULL)
        return FAIL;
    }
  return OK;
}


typedef int(*cmp_func)(const void *, const void *);

#define BASE_DOUBLE
#include "pnl_templates_on.h"
static char pnl_mat_label[] = "PnlMatrixDouble";
static char pnl_hmat_label[] = "PnlHmatrixDouble";
#include "pnl_TMatrix_source.c"
#include "dgemm_source.c"
#include "pnl_templates_off.h"
#undef  BASE_DOUBLE


#define BASE_PNL_COMPLEX
#include "pnl_templates_on.h"
static char pnl_mat_complex_label[] = "PnlMatrixComplex";
static char pnl_hmat_complex_label[] = "PnlHmatrixComplex";
#include "pnl_TMatrix_source.c"
#include "dgemm_source.c"
#include "pnl_templates_off.h"
#undef  BASE_PNL_COMPLEX

#define BASE_INT
#include "pnl_templates_on.h"
static char pnl_mat_int_label[] = "PnlMatrixInt";
static char pnl_hmat_int_label[] = "PnlHmatrixInt";
#include "pnl_TMatrix_source.c"
#include "dgemm_source.c"
#include "pnl_templates_off.h"
#undef  BASE_INT
