
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
PnlPermutation* pnl_permutation_new ()
{
  PnlPermutation *p;
  p = pnl_vect_int_new ();
  return p;
}

/**
 * Creates a PnlPermutation
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
 * Frees a PnlPermutation
 *
 * @param p address of a permutation
 */
void pnl_permutation_free (PnlPermutation **p)
{
  pnl_vect_int_free (p);
}

/**
 * Prints a permutation to a file
 *
 * @param fic a File decriptor
 * @param p the permutation to print
 */
void pnl_permutation_fprint (FILE *fic, const PnlPermutation *p)
{
  pnl_vect_int_fprint (fic, p);
}

/**
 * Prints a permutation to the standard output
 *
 * @param p the permutation to print
 */
void pnl_permutation_print (const PnlPermutation *p)
{ 
  pnl_vect_int_print (p);
}

