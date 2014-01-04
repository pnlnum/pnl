
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


#include "pnl/pnl_config.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_internals.h"

static char pnl_vector_label[] = "PnlVectObject";

/**
 * Create a PnlVectObject which is the parent type of all vectors
 *
 * @return a PnlVectObject
 */
PnlVectObject* pnl_vect_object_new ()
{
  PnlVectObject *o;
  if ( (o = malloc (sizeof (PnlVectObject))) == NULL) return NULL;
  o->size = 0;
  o->mem_size = 0;
  o->owner = 1;
  o->array = NULL;
  o->object.type = PNL_TYPE_VECTOR;
  o->object.parent_type = PNL_TYPE_VECTOR;
  o->object.nref = 0;
  o->object.label = pnl_vector_label;
  o->object.destroy = (DestroyFunc *) pnl_vect_object_free;
  o->object.constructor = (NewFunc *) pnl_vect_object_new;
  o->object.clone = NULL;
  o->object.copy = NULL;
  return o;
}

/**
 * free a TYPE(PnlVect)pointer and set the data pointer to
 * NULL
 *
 * @param v address of the pointer to free
 */
void pnl_vect_object_free(PnlVectObject **v)
{
  if (*v != NULL)
    {
      if ((*v)->array != NULL && (*v)->owner == 1)
        { free((*v)->array); (*v)->array = NULL; }
      (*v)->size = 0;
      free(*v);
    }
}

/**
 * Resize a PnlVectObject. If the new size is smaller than the
 * current one, no memory is freed and the datas are
 * kept. If the new size is larger than the current one, a
 * new pointer is allocated. The old datas are kept. 
 *
 * @param v a pointer to an already existing PnlVectObject. If v->owner=0,
 * nothing is done
 * @param size the new size of the array
 */
int pnl_vect_object_resize(PnlVectObject * v, int size)
{

  size_t sizeof_base = 0;
  if (v->owner == 0) return OK;
  if (size < 0) return FAIL;
  if (size == 0)
    {
      if (v->mem_size > 0) free (v->array);
      v->size = 0;
      v->mem_size = 0;
      v->array=NULL;
      return OK;
    }
  
  if (v->mem_size >= size)
    {
      /* If the new size is smaller, we do not reduce the size of the
         allocated block. It may change, but it allows to grow the vector
         quicker */
      v->size=size; return OK;
    }

  /* Now, v->mem_size < size */
  switch (PNL_GET_TYPE(v))
    {
    case PNL_TYPE_VECTOR_DOUBLE : sizeof_base = sizeof(double);
      break;
    case PNL_TYPE_VECTOR_COMPLEX : sizeof_base = sizeof(dcomplex);
      break;
    case PNL_TYPE_VECTOR_INT : sizeof_base = sizeof(int);
      break;
    default : PNL_MESSAGE (1, "Unknown type in pnl_vect_object_resize.\n");
      return FAIL;
    }
  if ((v->array=realloc(v->array,size*sizeof_base)) == NULL) return FAIL;
  v->size = size;
  v->mem_size = size;
  return OK;
}


#define BASE_DOUBLE
#include "pnl/pnl_templates_on.h"
static char pnl_vect_label[] = "PnlVectorDouble";
#include "vector_source.c"
#include "pnl/pnl_templates_off.h"
#undef BASE_DOUBLE
 
#define BASE_PNL_COMPLEX
#include "pnl/pnl_templates_on.h"
static char pnl_vect_complex_label[] = "PnlVectorComplex";
#include "vector_source.c"
#include "pnl/pnl_templates_off.h"
#undef BASE_PNL_COMPLEX


#define BASE_INT
#include "pnl/pnl_templates_on.h"
static char pnl_vect_int_label[] = "PnlVectorInt";
#include "vector_source.c"
#include "pnl/pnl_templates_off.h"
#undef  BASE_INT

