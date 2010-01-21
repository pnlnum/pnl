
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

#include "pnl_vector.h"
#include "pnl_perm.h"
#include "config.h"
#include "pnl_mathtools.h"

/**
 * Allocation of a PnlPermutation
 *
 * @param n length of the permutation
 */
PnlPermutation* pnl_permutation_create (int n)
{
  PnlPermutation *p;
  if ((p=malloc(sizeof(PnlPermutation)))==NULL) return NULL;
  p->size = n;
  if (p->size > 0)
    {
      if ((p->array=malloc(sizeof(int)*n))==NULL) return NULL;
    }
  else p->array = NULL;
  return p;
}

PnlPermutation* pnl_permutation_create_from_ptr (int n, const int *t)
{
  int i;
  PnlPermutation *p;
  p = pnl_permutation_create (n);
  for (i=0; i<n; i++)
    p->array[i] = t[i];
  return p;
}

/**
 * Initializes an existing permutation to the identity permutation
 *
 * @param p a PnlPermutation
 */
void pnl_permutation_init (PnlPermutation *p)
{
  int i;
  for (i=0; i<p->size; i++)
    {
      p->array[i] = i;
    }
}

/**
 * Swaps two elements of a permutation
 *
 * @param p a PnlPermutation
 * @param i first index
 * @param j second index
 */
void pnl_permutation_swap (PnlPermutation *p, int i, int j)
{
#ifndef PNL_RANGE_CHECK_OFF
  if (i >= p->size)
    {
      PNL_ERROR("first index is out of range", "pnl_permutation_swap");
    }
    
  if (j >= p->size)
    {
      PNL_ERROR("second index is out of range", "pnl_permutation_swap");
    }
#endif
  
  if (i != j)
    {
      int tmp = p->array[i];
      p->array[i] = p->array[j];
      p->array[j] = tmp;
    }
}

/**
 * Frees a PnlPermutation
 *
 * @param p address of a permutation
 */
void pnl_permutation_free (PnlPermutation **p)
{
  if ( *p == NULL ) return ;
  if ( (*p)->size >0 ) free ((*p)->array);
  free (*p); *p=NULL;
}

/**
 * Applies a PnlPermutation to a PnlVect
 * px[i] = x[p[i]]
 *
 * @param px at exit contains the permutated vector
 * @param x the vector to permute
 * @param p a permutation
 */
void pnl_vect_permute (PnlVect *px, const PnlVect *x, const PnlPermutation *p)
{
  int i, k;
#ifndef  PNL_RANGE_CHECK_OFF
  if (x->size != p->size)
    {
      PNL_ERROR("incompatible permutation size", "pnl_vect_permute");
    }
#endif
  pnl_vect_resize (px, x->size);
  for (i=0; i<x->size; i++)
    {
      k = p->array[i];
      pnl_vect_set (px, i, pnl_vect_get (x, k));
    }
}

/**
 * Applies a Permutation to an array in place
 * x[i] = x[p[i]]
 *
 * @param x the vector to permute. Contains the permuted vector on exit
 * @param p a permutation
 *
 * This algorithm comes
 * From Knuth "Sorting and Searching", Volume 3 (3rd ed), Section 5.2
 *  Exercise 10 (answers), p 617
 * It is based on the decomposition of any permutation into disjoined
 * cycles. First we search for the cycles and then apply each of them inplace,
 * which is much easier.
 */
static void pnl_permute_inplace (double *x, const int *p, int n)
{
  int i, k, pk;
  double t;
  for (i = 0; i < n; i++)
    {
      k = p[i];

      /* we are looking for the last element of the current cycle */
      while (k > i) k = p[k];

      pk = p[k];
      if (k == i && pk != i)
        {
          /* we have found the last of element of the current cycle.
             if k <i, then i was not a cycle leader
             if pk == 1, this is not a true cycle, but instead  i is a fixed
             point
             now we also have pk = p[i] 
          */
          
          /* shuffle the elements of the cycle starting from i */
          t = x[i];
          while (pk != i)
            {
              x[k] = x[pk];
              k = pk;
              pk = p[k];
            }
          x[k] = t;
        }
    }
}

/**
 * Applies a Permutation to a PnlVect in place
 *
 * @param x the vector to permute. Contains the permuted vector on exit
 * @param p a permutation
 */
void pnl_vect_permute_inplace (PnlVect *x, const PnlPermutation *p)
{
#ifndef PNL_RANGE_CHECK_OFF
  if (x->size != p->size)
    {
      PNL_ERROR("incompatible permutation size", "pnl_vect_permute");
    }
#endif
  pnl_permute_inplace (x->array, p->array, x->size);
}
  
/**
 * Prints a permutation to a file
 *
 * @param fic a File decriptor
 * @param p the permutation to print
 */
void pnl_permutation_fprint (FILE *fic, const PnlPermutation *p)
{
  int i;
  for (i=0; i<p->size; i++)
    {
      fprintf (fic, "%d ", p->array[i]);
    }
  fprintf (fic, "\n");
}

/**
 * Prints a permutation to the standard output
 *
 * @param p the permutation to print
 */
void pnl_permutation_print (const PnlPermutation *p)
{ pnl_permutation_fprint (stdout, p); }

