
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
#include "pnl_vector.h"

static char pnl_vector_label[] = "PnlVectObject";

/**
 * Creates a PnlVectObject which is the parent type of all vectors
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
  o->object.label = pnl_vector_label;
  return o;
}


#define BASE_DOUBLE
#include "pnl_templates_on.h"
static char pnl_vect_label[] = "PnlVectorDouble";
#include "pnl_TVector_source.c"
#include "pnl_templates_off.h"
#undef BASE_DOUBLE
 
#define BASE_PNL_COMPLEX
#include "pnl_templates_on.h"
static char pnl_vect_complex_label[] = "PnlVectorComplex";
#include "pnl_TVector_source.c"
#include "pnl_templates_off.h"
#undef BASE_PNL_COMPLEX


#define BASE_INT
#include "pnl_templates_on.h"
static char pnl_vect_int_label[] = "PnlVectorInt";
#include "pnl_TVector_source.c"
#include "pnl_templates_off.h"
#undef  BASE_INT

