
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
#include "pnl/pnl_array.h"
#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"

static char pnl_array_label[] = "PnlArray";

/**
 * Creates an empty array
 */
PnlArray* pnl_array_new ()
{
  PnlArray *o;
  if ( (o = malloc (sizeof(PnlArray))) == NULL ) return NULL;
  o->array = NULL;
  o->size = 0;
  o->mem_size = 0;
  o->object.type = PNL_TYPE_ARRAY;
  o->object.parent_type = PNL_TYPE_ARRAY;
  o->object.label = pnl_array_label;
  o->object.destroy = (destroy_func *) pnl_array_free;
  return o;
}

/** 
 * Creates an array of size n
 * 
 * @param n size of the array
 * 
 * @return  an array
 */
PnlArray* pnl_array_create (int n)
{
  PnlArray *T;
  PnlObject **A;
  if ( (A=malloc (n*sizeof (PnlObject *))) == NULL) return NULL;
  T = pnl_array_new ();
  T->size = T->mem_size = n;
  T->array = A;
  return T;
}

/**
 * Resizes a PnlArray. If the new size is smaller than the
 * current one, no memory is freed and the datas are
 * kept. If the new size is larger than the current one, a
 * new pointer is allocated. The old datas are kept. 
 *
 * @param v a pointer to an already existing PnlArray. 
 * @param size the new size of the array
 */
int pnl_object_resize(PnlArray * v, int size)
{

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
  if ((v->array=realloc(v->array,size*sizeof(PnlObject *))) == NULL) return FAIL;
  v->size = size;
  v->mem_size = size;
  return OK;
}


/**
 * Returns the adress of the i-th element of an array. No copy is made
 *
 * @param T a PnlArray
 * @param i an interger index
 *
 * @return a PnlObject*
 */
PnlObject* pnl_array_get (PnlArray *T, int i)
{
  PNL_CHECK (i >= T->size, "index exceeded", "pnl_array_get");
  return T->array[i];
}

/**
 * Sets the i-th element of an array to the object O. No copy is made
 *
 * @param T a PnlArray
 * @param i an interger index
 * @param O a PnlObject
 *
 * @return a PnlObject*
 */
void pnl_array_set (PnlArray *T, int i, PnlObject *O)
{
  PNL_CHECK (i >= T->size, "index exceeded", "pnl_array_set");
  T->array[i] = O;
}

/**
 * Frees a PnlArray
 *
 * @param T the address of a PnlArray
 */
void pnl_array_free (PnlArray **T)
{
  int i=0;
  if (*T == NULL) return;

  for ( i=0 ; i< (*T)->size ; i++ )
    {
      PnlObject *O = (*T)->array[i];
      if ( O != NULL ) O->destroy ((void **) &O);
    }
  if ( (*T)->array != NULL ) free ((*T)->array);
  free (*T);
  *T = NULL;
}

/**
 * Prints the typename of everything element stored in the array.
 * When the PnlObject structure has a copy_func field, we will use it to
 * actually print the content of each element
 *
 * @param T a array
 */
void pnl_array_print (PnlArray *T)
{
  int i;
  for ( i=0 ; i<T->size ; i++ )
    {
      const PnlObject *C;
      const char *name;
      C = pnl_array_get (T, i);
      name =  PNL_GET_TYPENAME(C);
      printf ("T[%d] : %s\n", i, name);
    }
}

