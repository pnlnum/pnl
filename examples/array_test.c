
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

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_array.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_tridiag_matrix.h"
#include "tests_utils.h"


static void array_ops ()
{
  PnlArray *T;
  PnlMat *M;
  PnlTridiagMat *tridiag;

  M = pnl_mat_create_from_double (2, 3, 3.5);
  tridiag = pnl_tridiag_mat_create_from_double (4, 0.5);
  T = pnl_array_create (2);
  pnl_array_set (T, 0, PNL_OBJECT(M));
  pnl_array_set (T, 1, PNL_OBJECT(tridiag));

  if ( pnl_array_get (T, 0) != PNL_OBJECT(M) || 
       pnl_array_get (T, 1) != PNL_OBJECT(tridiag))
    {
      pnl_test_set_fail ("pnl_array_{g,s}et", 0., 0.);
      goto J1;
    }
  pnl_test_set_ok ("pnl_array_{g,s}et");
J1:
  pnl_array_free (&T);
 }



int main (int argc, char **argv)
{
  pnl_test_init (argc, argv);
  array_ops ();
  exit (pnl_test_finalize ("Array operations"));
}
