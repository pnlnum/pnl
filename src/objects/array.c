
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
 * Create an empty array
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
  o->object.nref = 0;
  o->object.label = pnl_array_label;
  o->object.destroy = (DestroyFunc *) pnl_array_free;
  o->object.constructor = (NewFunc *) pnl_array_new;
  o->object.clone = (CloneFunc *) NULL;
  o->object.copy = (CopyFunc *) NULL;
  return o;
}

/** 
 * Create an array of size n. Each element is set to NULL.
 * 
 * @param n size of the array
 * 
 * @return  an array
 */
PnlArray* pnl_array_create (int n)
{
  int i;
  PnlArray *T;
  PnlObject **A;
  if ( (A=malloc (n*sizeof (PnlObject *))) == NULL) return NULL;
  for ( i=0 ; i<n ; i++ ) A[i] = NULL;
  T = pnl_array_new ();
  T->size = T->mem_size = n;
  T->array = A;
  return T;
}

/** 
 * Copy an Array
 * 
 * @param A
 * 
 * @return 
 */
PnlArray* pnl_array_copy (const PnlArray *A)
{
  int i;
  PnlArray *C;
  C = pnl_array_create (A->size);
  for ( i=0 ; i<C->size; i++ )
    {
      PnlObject *Ai = pnl_array_get (A, i);
      PnlObject *Ci = Ai->copy(Ai);
      pnl_array_set (C, i, Ci);
    }
  return C;
}

/** 
 * Clone an Array
 *
 * For each element, it tries to call its clone member, but A[i] and C[i]
 * are not of the same type, C[i] is destroyed and A[i]->copy is called to
 * re--create C[i].
 * 
 * @param C destination
 * @param A source
 */
void pnl_array_clone (PnlArray *C, const PnlArray *A)
{
  int i;
  pnl_array_resize (C, A->size);
  for ( i=0 ; i<C->size; i++ )
    {
      PnlObject *Ai = pnl_array_get (A, i);
      PnlObject *Ci = pnl_array_get (C, i);
      /* By construction, empty cells are set to NULL */
      if ( Ci != NULL  && PNL_GET_TYPE(Ci) != PNL_GET_TYPE(Ai) )
        {
          void *Ci_void = (void *)Ci;
          if ( Ci->nref > 1 )
            Ci->nref --;
          else
            Ci->destroy(&Ci_void); 
          Ci = NULL;
        }
      else if ( Ci == NULL )
        {
          Ci = Ai->copy(Ai);
          pnl_array_set (C, i, Ci);
        }
      else
        {
          Ai->clone(Ci, Ai);
        }
    }
}

/**
 * Resize a PnlArray. 
 *
 * If the new size is smaller than the current one, no memory is freed and
 * the datas are kept.
 *
 * If the new size is larger than the current one, a new pointer is
 * allocated. The old datas are kept and the new cells are set to NULL.
 *
 * @param v a pointer to an already existing PnlArray. 
 * @param size the new size of the array
 */
int pnl_array_resize(PnlArray * v, int size)
{

  int i, old_size = v->size;
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
  for ( i=old_size ; i<size ; i++ ) v->array[i] = NULL;
  v->size = size;
  v->mem_size = size;
  return OK;
}

/**
 * Return the adress of the i-th element of an array. No copy is made
 *
 * @param T a PnlArray
 * @param i an interger index
 *
 * @return a PnlObject*
 */
PnlObject* pnl_array_get (const PnlArray *T, int i)
{
  PNL_CHECK (i >= T->size, "index exceeded", "pnl_array_get");
  return T->array[i];
}

/**
 * Set the i-th element of an array to the object O. No copy is made
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
  /* update nref for O */
  if ( O != NULL ) O->nref ++;

}

/**
 * Free a PnlArray
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
      if ( O != NULL ) 
        {
          if ( O->nref > 1 )
            O->nref --;
          else
            {
              void *O_void = (void *) O;
              O->destroy (&O_void);
              O = NULL;
            }
        }
    }
  if ( (*T)->array != NULL ) free ((*T)->array);
  free (*T);
  *T = NULL;
}

/**
 * Print the typename of everything element stored in the array.
 * When the PnlObject structure has a CopyFunc field, we will use it to
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

