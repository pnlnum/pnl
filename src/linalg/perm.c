
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

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_perm.h"
#include "pnl/pnl_config.h"
#include "pnl/pnl_mathtools.h"

/**
 * Create an empty PnlPermutation
 */
PnlPermutation* pnl_permutation_new ()
{
  PnlPermutation *p;
  p = pnl_vect_int_create (0);
  return p;
}

/**
 * Create a PnlPermutation
 *
 * @param n length of the permutation
 */
PnlPermutation* pnl_permutation_create (int n)
{
  PnlPermutation *p;
  p = pnl_vect_int_create (n);
  return p;
}

/** 
 * Compute the inverse of a permutation vector
 * 
 * @param inv contains the inverse on output
 * @param p a permutation
 * 
 */
void pnl_permutation_inverse (PnlPermutation *inv, const PnlPermutation *p)
{
  int i;
  pnl_vect_int_resize (inv, p->size);
  
  for ( i=0 ; i<p->size ; i++ )
    {
      inv->array[p->array[i]] = i;
    }
}

/**
 * Free a PnlPermutation
 *
 * @param p address of a permutation
 */
void pnl_permutation_free (PnlPermutation **p)
{
  pnl_vect_int_free (p);
}

/**
 * Apply a PnlPermutation to a PnlVect
 * px[i] = x[p[i]]
 *
 * @param px on exit contains the permuted vector
 * @param x the vector to permute
 * @param p a permutation
 */
void pnl_vect_permute (PnlVect *px, const PnlVect *x, const PnlPermutation *p)
{
  int i, k;
  PNL_CHECK (x->size != p->size, "incompatible permutation size", "pnl_vect_permute");
  pnl_vect_resize (px, x->size);
  for (i=0; i<x->size; i++)
    {
      k = p->array[i];
      PNL_LET (px, i) = PNL_GET (x, k);
    }
}

/**
 * Apply the inverse of a PnlPermutation to a PnlVect
 * px[p[i]] = x[i]
 *
 * @param px on exit contains the permuted vector
 * @param x the vector to permute
 * @param p a permutation
 */
void pnl_vect_permute_inverse (PnlVect *px, const PnlVect *x, const PnlPermutation *p)
{
  int i;
  PNL_CHECK (x->size != p->size, "incompatible permutation size", "pnl_vect_permute_inverse");
  pnl_vect_resize (px, x->size);
  for (i=0; i<x->size; i++)
    {
      const int k = p->array[i];
      LET (px, k) = GET (x, i);
    }
}

/**
 * Apply a Permutation to an array in place
 * x[i] = x[p[i]]
 *
 * @param x a C array of real values to permute. On exit, contains the permuted data.
 * @param p a C arary of integers representing a permutation
 * @param n the size of the array x (it is also the size of p)
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
 * Apply the inverse of a Permutation to an array in place
 * x[p[i]] = x[i]
 *
 * @param x a C array of real values to permute. On exit, contains the permuted data.
 * @param p a C arary of integers representing a permutation
 * @param n the size of the array x (it is also the size of p)
 *
 * This algorithm comes
 * From Knuth "Sorting and Searching", Volume 3 (3rd ed), Section 5.2
 *  Exercise 10 (answers), p 617
 * It is based on the decomposition of any permutation into disjoined
 * cycles. First we search for the cycles and then apply each of them inplace,
 * which is much easier.
 */
static void pnl_permute_inverse_inplace (double *x, const int *p, int n)
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
          
          /* shuffle the elements of the cycle starting from k */
          t = x[k];
          while (pk != i)
            {
              const double tmp = x[pk];
              x[pk] = t;
              t = tmp;
              k = pk;
              pk = p[k];
            }
          x[pk] = t;
        }
    }
}

/**
 * Apply a Permutation to a PnlVect in place
 *
 * @param x the vector to permute. Contains the permuted vector on exit
 * @param p a permutation
 */
void pnl_vect_permute_inplace (PnlVect *x, const PnlPermutation *p)
{
  PNL_CHECK (x->size != p->size, "incompatible permutation size", "pnl_vect_permute");
  pnl_permute_inplace (x->array, p->array, x->size);
}
  
/**
 * Apply the inverse of a Permutation to a PnlVect in place
 *
 * @param x the vector to permute. Contains the permuted vector on exit
 * @param p a permutation
 */
void pnl_vect_permute_inverse_inplace (PnlVect *x, const PnlPermutation *p)
{
  PNL_CHECK (x->size != p->size, "incompatible permutation size", "pnl_vect_permute_inverse");
  pnl_permute_inverse_inplace (x->array, p->array, x->size);
}

/**
 * Apply a PnlPermutation to the columns of a matrix
 * pX[] = x[:,p]
 *
 * @param pX on exit contains the permuted matrix
 * @param X the matrix to permute
 * @param p a permutation
 */
void pnl_mat_col_permute (PnlMat *pX, const PnlMat *X, const PnlPermutation *p)
{
  int i, j, k;
  PNL_CHECK (X->n != p->size, "incompatible permutation size", "pnl_mat_col_permute");
  pnl_mat_resize (pX, X->m, X->n);
  for ( j=0 ; j<X->n ; j++ )
    {
      k = p->array[j];
      for ( i=0 ; i<X->m ; i++ )
        {
          MLET (pX, i, j) = MGET(X, i, k);
        }
    }
}

/**
 * Apply a PnlPermutation to the rows of a matrix
 * pX[] = x[p,:]
 *
 * @param pX on exit contains the permuted matrix
 * @param X the matrix to permute
 * @param p a permutation
 */
void pnl_mat_row_permute (PnlMat *pX, const PnlMat *X, const PnlPermutation *p)
{
  int i, j, k;
  PNL_CHECK (X->m != p->size, "incompatible permutation size", "pnl_mat_col_permute");
  pnl_mat_resize (pX, X->m, X->n);
  for ( i=0 ; i<X->m ; i++ )
    {
      k = p->array[i];
      for ( j=0 ; j<X->n ; j++ )
        {
          MLET (pX, i, j) = MGET(X, k, j);
        }
    }
}

/**
 * Print a permutation to a file
 *
 * @param fic a File decriptor
 * @param p the permutation to print
 */
void pnl_permutation_fprint (FILE *fic, const PnlPermutation *p)
{
  pnl_vect_int_fprint (fic, p);
}

/**
 * Print a permutation to the standard output
 *
 * @param p the permutation to print
 */
void pnl_permutation_print (const PnlPermutation *p)
{ 
  pnl_vect_int_print (p);
}

