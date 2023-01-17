
/************************************************************************/
/* Copyright J�r�me Lelong <jerome.lelong@gmail.com>                    */
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
#include <string.h>

#include "pnl/pnl_config.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_sp_matrix.h"
#include "pnl/pnl_internals.h"


static char pnl_sp_matrix_label[] = "PnlSpMatObject";

/**
 * Create a new sparse matrix object
 */
PnlSpMatObject *pnl_sp_mat_object_new()
{
  PnlSpMatObject *o;
  if ((o = malloc(sizeof(PnlSpMatObject))) == NULL) return NULL;
  o->object.type = PNL_TYPE_SP_MATRIX;
  o->object.parent_type = PNL_TYPE_SP_MATRIX;
  o->object.nref = 0;
  o->object.label = pnl_sp_matrix_label;
  o->object.destroy = (DestroyFunc *) pnl_sp_mat_object_free;
  o->object.constructor = (NewFunc *) pnl_sp_mat_object_new;
  o->object.clone = NULL;
  o->object.copy = NULL;
  o->m = o->n = o->nz = 0;
  o->nzmax = 0;
  o->array = NULL;
  o->I = NULL;
  o->J = NULL;
  return o;
}

/**
 * Free a sparse matrix object
 */
void pnl_sp_mat_object_free(PnlSpMatObject **o)
{
  if ((*o != NULL))
    {
      if ((*o)->array != NULL)
        {
          free((*o)->array);
          (*o)->array = NULL;
        }
      if ((*o)->I != NULL)
        {
          free((*o)->I);
          (*o)->I = NULL;
        }
      if ((*o)->J != NULL)
        {
          free((*o)->J);
          (*o)->J = NULL;
        }
      (*o)->m = (*o)->n = (*o)->nz = 0;
      free(*o);
    }
}

/**
 * Resize a PnlSpMatObject.  If the new size is smaller than the current one, no
 * memory is freed. If the new size is larger than the current nzmax, a new
 * pointer is allocated. The old data are kept only when m is left unchanged
 * and we increase nzmax.
 *
 * @param o a pointer to an already existing PnlSpMatObject
 * @param m new nb of rows
 * @param n new nb of columns
 * @param nz new maximum number of non-zero elements.
 *
 * @return PNL_OK or PNL_FAIL. When returns PNL_OK, the matrix M is changed.
 */
int pnl_sp_mat_object_resize(PnlSpMatObject *o, int m, int n, int nz)
{
  if (m < 0 || n < 0 || nz < 0) return PNL_FAIL;
  if (m == 0 || n == 0)
    {
      o->m = o->n = 0;
      o->nz = o->nzmax = 0;
      if (o->array)
        {
          free(o->array);
          o->array = NULL;
        }
      return PNL_OK;
    }

  /* Increase the number of rows if needed */
  if (m > o->m)
    {
      int i;
      o->I = realloc(o->I, (m + 1) * sizeof(int));
      for (i = o->m + 1 ; i < m + 1 ; i++)
        {
          o->I[i] = o->I[o->m];
        }
    }
  o->m = m;

  /* Increase the maximum number of non-zero elements */
  if (nz > o->nzmax)
    {
      size_t sizeof_base = 0;
      switch (PNL_GET_TYPE(o))
        {
        case PNL_TYPE_SP_MATRIX_DOUBLE:
          sizeof_base = sizeof(double);
          break;
        case PNL_TYPE_SP_MATRIX_COMPLEX :
          sizeof_base = sizeof(dcomplex);
          break;
        case PNL_TYPE_SP_MATRIX_INT :
          sizeof_base = sizeof(int);
          break;
        default :
          PNL_MESSAGE(1, "Unknown type in pnl_sp_mat_object_resize.\n");
          return PNL_FAIL;
        }
      o->array = realloc(o->array, nz * sizeof_base);
      o->J = realloc(o->J, nz * sizeof(int));
      o->nzmax = nz;
    }
  if (o->n > n)
    {
      int i;
      /* Some values in J must be out of range (greater than n-1) so we empty
       * the matrix by setting all the elements to zero */
      o->nz = 0;
      for (i = 0 ; i < m + 1 ; i++)
        {
          o->I[i] = 0;
        }
      for (i = 0 ; i < o->nzmax ; i++)
        {
          o->J[i] = 0;
        }
    }
  o->n = n;
  return PNL_OK;
}

#define BASE_DOUBLE
#include "pnl/pnl_templates_on.h"
static char pnl_sp_mat_label[] = "PnlSpMatrixDouble";
#include "sp_matrix_source.c"
#include "pnl/pnl_templates_off.h"
#undef  BASE_DOUBLE


#define BASE_PNL_COMPLEX
#include "pnl/pnl_templates_on.h"
static char pnl_sp_mat_complex_label[] = "PnlSpMatrixComplex";
#include "sp_matrix_source.c"
#include "pnl/pnl_templates_off.h"
#undef  BASE_PNL_COMPLEX

#define BASE_INT
#include "pnl/pnl_templates_on.h"
static char pnl_sp_mat_int_label[] = "PnlSpMatrixInt";
#include "sp_matrix_source.c"
#include "pnl/pnl_templates_off.h"
#undef  BASE_INT

/**
 * Test if two integer sparse matrices are equal
 *
 * @param Sp1  a sparse matrix
 * @param Sp2  a sparse matrix
 * @return PNL_TRUE or PNL_FALSE
 */
int pnl_sp_mat_int_isequal(const PnlSpMatInt *Sp1, const PnlSpMatInt *Sp2)
{
  int k;
  if ((Sp1->m != Sp2->m) || (Sp1->n != Sp2->n) || (Sp1->nz != Sp2->nz)) return PNL_FALSE;
  for (k = 0; k < Sp1->m; k++)
    {
      if (Sp1->I[k] != Sp2->I[k]) return PNL_FALSE;
    }
  for (k = 0; k < Sp1->nz; k++)
    {
      if ((Sp1->J[k] != Sp2->J[k]) || (Sp1->array[k] != Sp2->array[k])) return PNL_FALSE;
    }
  return PNL_TRUE;
}