
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
#include "pnl/pnl_list.h"
#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"

static char pnl_list_label[] = "PnlList";

/**
 * Creates an empty list
 */
PnlList* pnl_list_new ()
{
  PnlList *o;
  if ( (o = malloc (sizeof(PnlList))) == NULL ) return NULL;
  o->first = NULL;
  o->last = NULL;
  o->curcell = NULL;
  o->icurcell = NULLINT;
  o->len = 0;
  o->object.type = PNL_TYPE_LIST;
  o->object.parent_type = PNL_TYPE_LIST;
  o->object.label = pnl_list_label;
  o->object.destroy = (destroy_func *) pnl_list_free;
  o->object.new = (new_func *) pnl_list_new;
  o->object.clone = (clone_func *) NULL;
  o->object.copy = (copy_func *) NULL;
  return o;
}

/**
 * Creates an empty Cell
 */
PnlCell* pnl_cell_new ()
{
  PnlCell *o;
  if ( (o = malloc (sizeof(PnlCell))) == NULL ) return NULL;
  o->prev = NULL;
  o->next = NULL;
  o->self = NULL;
  return o;
}

/**
 * Returns the adress of the i-th element of a list. No copy is made
 *
 * @param L a PnlList
 * @param i an interger index
 *
 * @return a PnlObject*
 */
PnlObject* pnl_list_get (PnlList *L, int i)
{
  int j;
  PnlCell *C;
  PNL_CHECK (i >= L->len, "index exceeded", "pnl_list_get");

  /*
   * First, we chek if we are linearly going through the list. If yes, rather
   * than starting at the beginning of the list, we use the curcell member to
   * speed up the access
   */
  if ( L->icurcell != NULLINT  && i == L->icurcell + 1 )
    {
      L->curcell = L->curcell->next;
      L->icurcell ++;
    }
  else
    {
      C = L->first;
      for ( j=0 ; j<i ; j++ )
        {
          C = C->next;
        }
      L->curcell = C;
      L->icurcell = i;
    }
  return L->curcell->self;
}

/**
 * Inserts a new object in the first position of a List
 *
 * @param L an already existing PnlList
 * @param o a PnlObject
 */
void pnl_list_insert_first (PnlList *L, PnlObject *o)
{
  PnlCell *C;

  C = pnl_cell_new ();

  /* Downward linkage */
  C->next = L->first;
  C->prev = NULL;
  C->self = o;

  /* Upward linkage */
  if (L->len >0)
    {
      PnlCell *second;
      second = C->next;
      second->prev = C;
    }

  /* Put C on top of L */
  if (L->len == 0) L->last = C;
  L->len++;
  L->first = C;

}


/**
 * Inserts a new object in the last position of a List
 *
 * @param L an already existing PnlList
 * @param o a PnlObject
 */
void pnl_list_insert_last (PnlList *L, PnlObject *o)
{
  PnlCell *C;
  C = pnl_cell_new ();

  /* Upward linkage */
  C->prev = L->last;
  C->next = NULL;
  C->self = o;

  /* Downward linkage */
  if (L->len > 0)
    {
      PnlCell *lastbut;
      lastbut = L->last;
      lastbut->next = C;
    }

  /* Put C at the bottom of L */
  if (L->len == 0) L->first = C;
  L->len++;
  L->last = C;
}

/**
 * Frees a PnlList
 *
 * @param L the address of a PnlList
 */
void pnl_list_free (PnlList **L)
{
  int i=0;
  PnlCell *node, *next;

  if (*L == NULL) return;
  node = (*L)->first;

  for ( i=0 ; i< (*L)->len ; i++ )
    {
      next = node->next;
      pnl_cell_free (&node);
      node = next;
    }
  free (*L);
  *L = NULL;
}

/**
 * Frees a PnlCell
 *
 * @param c the address of a PnlCell
 */
void pnl_cell_free (PnlCell **c)
{
  PnlObject *O;
  if (*c == NULL) return;

  /*
   * destroy the content of the cell
   */
  O = (*c)->self;
  if ( O != NULL )
    {
      (*(O->destroy)) ((void **) &O);
      O = NULL;
    }
  /*
   * Update downward and upward linkage
   */
  if ( (*c)->prev != NULL ) { (*c)->prev->next = (*c)->next; }
  if ( (*c)->next != NULL ) { (*c)->next->prev = (*c)->prev; }
  free (*c); *c = NULL;
}

/**
 * Removes the last cell of a List
 *
 * @param L an already existing PnlList
 */
void pnl_list_remove_last (PnlList *L)
{
  PnlCell *last_but, *last;

  if (L->len == 0) return;
  if (L->len == 1)
    {
      pnl_cell_free (&(L->first));
      L->first = NULL;
      L->last = NULL;
      L->len = 0;
      return;
    }

  last = L->last;
  last_but = last->prev;

  last_but->next = NULL;

  L->len--;
  L->last = last_but;

  pnl_cell_free (&last);
}

/**
 * Removes the first cell of a List
 *
 * @param L an already existing PnlList
 */
void pnl_list_remove_first (PnlList *L)
{
  PnlCell *C;

  if (L->len == 0) return;
  if (L->len == 1)
    {
      pnl_cell_free (&(L->first));
      L->first = NULL;
      L->last = NULL;
      L->len = 0;
      return;
    }

  C = L->first;
  L->first = C->next;
  L->first->prev = NULL;
  L->len--;

  pnl_cell_free (&C);
}

/**
 * Concatenate two lists.
 *
 * @param L1 (input/output) On input contains a PnlList. On output, it
 * contains the result of the concatenation
 * @param L2 (input) Contains a PnlList which will be appended to L1. L2 is
 * destroyed on output, which means that you do not want to free it yourself
 */
void pnl_list_concat (PnlList *L1, PnlList *L2)
{
  PnlCell *last1, *first2;

  last1 = L1->last;
  first2 = L2->first;

  L1->last = L2->last;
  last1->next = first2;
  first2->prev = last1;
  L1->len += L2->len;

  free (L2);
}

/**
 * Prints the typename of everything element stored in the list.
 * When the PnlObject structure has a copy_func field, we will use it to
 * actually print the content of each element
 *
 * @param L a list
 */
void pnl_list_print (PnlList *L)
{
  int i;
  for ( i=0 ; i<L->len ; i++ )
    {
      const PnlObject *C;
      const char *name;
      C = pnl_list_get (L, i);
      name =  PNL_GET_TYPENAME(C);
      printf ("L[%d] : %s\n", i, name);
    }
}

/**
 * Removes and frees cell i of a list
 *
 * @param L a PnlList
 * @param i the index of the cell to be removed
 */
void pnl_list_remove_i (PnlList *L, int i)
{
  PnlCell *C;
  int k;
  PNL_CHECK(L->len <= i,"index of out range", "pnl_list_remove_i");

  if ( i==0 ) { pnl_list_remove_first (L); return; }
  if ( i==L->len ) { pnl_list_remove_last (L); return; }

  C = L->first;
  for ( k=0 ; k<i ; k++ ) { C = C->next; }
  pnl_cell_free (&C);
  L->len--;
}
