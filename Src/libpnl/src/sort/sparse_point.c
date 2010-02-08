
/************************************************************************/
/* Copyright David Pommier <pommier.david@gmail.com>                    */
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

#include "pnl_list.h"


/**
 * allocates a contains.
 * @param ind key 
 * @param val value
 * @return a pointeur to PnlSparsePoint
 */
PnlSparsePoint * pnl_sparse_point_create(const PnlVectUint * ind,int val)
{
  PnlSparsePoint *C;
  if((C=malloc(sizeof(PnlSparsePoint)))==NULL) return NULL;
  C->index=pnl_vect_uint_copy(ind);
  C->value=val;
  return C;
}

/**
 * allocates a contains.
 * @param ind key 
 * @param val value
 * @return a pointeur to PnlSparsePoint
 */
PnlSparsePoint * pnl_sparse_point_clone(PnlVectUint * ind,int val)
{
  PnlSparsePoint *C;
  if((C=malloc(sizeof(PnlSparsePoint)))==NULL) return NULL;
  C->index=malloc(sizeof(PnlVectUint));
  C->index->owner=0;
  C->index->size=ind->size;
  C->index->array=&(ind->array[0]);
  //C->index=ind;
  C->value=val;
  return C;
}

/**
 * allocates a contains - copy constructor.
 * @param C2 contains pointer 
 * @return a pointeur to PnlSparsePoint
 */
PnlSparsePoint * pnl_sparse_point_copy(const PnlSparsePoint *C2)
{
  PnlSparsePoint *C;
  if((C=malloc(sizeof(PnlSparsePoint)))==NULL) return NULL;
  C->index=pnl_vect_uint_copy(C2->index);
  C->value=C2->value;
  return C;
}

/**
 * free a contains
 * @param C address of a contains
 */
void pnl_sparse_point_free(PnlSparsePoint **C)
{
  if (*C != NULL)
    {
      pnl_vect_uint_free(&((*C)->index));
      free(*C);
      *C=NULL;
    }
}

/**
 * Prints a contains to a file 
 *
 * @param fic a file descriptor.
 * @param C a Contains pointer.
 */
void pnl_sparse_point_fprint(FILE *fic,PnlSparsePoint *C)
{
  pnl_vect_uint_print(C->index);
  fprintf(fic," Index  %d \n",C->value);
} 


/**
 * Add - do nothing in this case 
 *
 * @param C a PnlSparsePoint pointer, C.Value Value.
 * @param C2 a Contains pointer.
 */
void  pnl_sparse_point_add(PnlSparsePoint *C,const PnlSparsePoint *C2)
{}

/**
 * Less compute relation C1<C2
 *
 * @param C1 a PnlSparsePoint pointer.
 * @param C2 a Contains pointer.
 * @return a int C1<C2
 */
int  pnl_sparse_point_less(const PnlSparsePoint *C1,const PnlSparsePoint *C2)
{return pnl_vect_uint_less(C1->index,C2->index);}

/**
 * Equal compute relation C1==C2
 *
 * @param C1 a PnlSparsePoint pointer.
 * @param C2 a Contains pointer.
 * @return a int C1==C2
 */
int  pnl_sparse_point_equal(const PnlSparsePoint *C1,const PnlSparsePoint *C2)
{return pnl_vect_uint_equal(C1->index,C2->index);}

