
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

#include "config.h"
#include "pnl_matrix.h"

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
  return o;
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
