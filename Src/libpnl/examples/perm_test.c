
/*************************************************************************/
/* Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>            */
/*                                                                       */
/* This program is free software; you can redistribute it and/or modify  */
/* it under the terms of the GNU General Public License as published by  */
/* the Free Software Foundation; either version 3 of the License, or     */
/* (at your option) any later version.                                   */
/*                                                                       */
/* This program is distributed in the hope that it will be useful,       */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/* GNU General Public License for more details.                          */
/*                                                                       */
/* You should have received a copy of the GNU General Public License     */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pnl_vector.h"
#include "pnl_perm.h"
#include "tests.h"

static void permute_vect_test ()
{
  PnlVect *x, *px;
  PnlPermutation *p;
  int t[10] = { 2, 4, 7, 9, 1, 0, 3, 6, 5, 8};

  printf("Test of permutation functions\n");
  x = pnl_vect_create_from_double(10, 1.);
  pnl_vect_cumsum (x);
  px = pnl_vect_create (0);
  p = pnl_permutation_create_from_ptr (10, t);

  printf("original vector\n");
  pnl_vect_print (x);
  printf ("permutation\n");
  pnl_permutation_print (p);
  printf ("permuted vector\n");
  pnl_vect_permute (px, x, p);
  pnl_vect_print (px);
  printf ("permuted vector in place\n");
  pnl_vect_permute_inplace (x, p);
  pnl_vect_print (x);

  pnl_vect_free (&x);
  pnl_vect_free (&px);
  pnl_permutation_free (&p);
}


void perm_test()
{
  permute_vect_test ();
}
