#include <stdlib.h>

#include "config.h"
#include "pnl_vector.h"
#include "pnl_matrix.h"

static char pnl_vector_label[] = "PnlVectObject";

/**
 * Creates a PnlVectObject which is the parent type of all vectors
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

