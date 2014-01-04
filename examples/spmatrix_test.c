
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_sp_matrix.h"
#include "pnl/pnl_random.h"
#include "tests_utils.h"

static PnlSpMat* create_random_sp (int m, int n, PnlRng *rng) 
{
  int it;
  int nz = floor (pnl_rng_uni (rng) * m * n * 0.5);
  PnlSpMat *M = pnl_sp_mat_create (m, n, nz);

  for ( it=0 ; it<nz ; it++ )
    {
      int i = (int) floor (pnl_rng_uni (rng) * m);
      int j = (int) floor (pnl_rng_uni (rng) * n);
      double x = pnl_rng_normal (rng);
      pnl_sp_mat_set (M, i, j, x);
    }
  return M;
}

static void sp_get_set_test ()
{
  PnlSpMat *Sp;
  PnlMat *M;
  int i, j, nb_fail;
  PnlRng * rng = pnl_rng_create(PNL_RNG_MERSENNE);
  int m = 9, n = 11;
  double abserr = 1E-18;
  pnl_rng_sseed (rng, 0);
  Sp = create_random_sp (m, n, rng);
  M = pnl_mat_create_from_sp_mat (Sp);

  nb_fail = 0;
  for ( i=0 ; i<m ; i++ )
    {
      for ( j=0 ; j<n ; j++ )
        {
          nb_fail += fabs(pnl_sp_mat_get (Sp, i, j) -  MGET(M, i, j)) > abserr;
        }
    }
  if ( nb_fail == 0 )
    {
      pnl_test_set_ok ("sp_mat_get_set");
    }
  else
    {
      pnl_test_set_fail ("sp_mat_get_set", nb_fail, 0);
    }
  pnl_sp_mat_free (&Sp);
  pnl_mat_free (&M);
  pnl_rng_free (&rng);
}

int main (int argc, char *argv[])
{
  pnl_test_init (argc, argv);
  sp_get_set_test ();
  /* sp_add_test (); */
  exit (pnl_test_finalize("Sparse matrices"));
}
