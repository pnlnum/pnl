
/************************************************************************/
/* Copyright Jérôme Lelong <jerome.lelong@gmail.com>                    */
/* Copyright David Pommier <david.pommier@gmail.com>                    */
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
#include "pnl/pnl_complex.h"
#include "pnl/pnl_random.h"
#include "tests_utils.h"

struct complex_tests
{
  char *label;
  dcomplex (*f) (dcomplex);
  double arg_r;
  double arg_i;
  double res_r;
  double res_i;
};

struct complex_tests list_tst [] =
{
#include "complex_test.dat"
    { NULL, NULL, 0, 0, 0, 0}
};


int main (int argc, char *argv[])
{
  int i;
  double tol  = 1E6 * DBL_EPSILON;
  pnl_test_init (argc, argv);
  if ( argc == 2 && strcmp (argv[1], "-v") == 0 )
    {
      verbose = 1;
    }
  for ( i=0 ; list_tst[i].f != NULL ; i++ )
    {
      struct complex_tests t = list_tst[i];
      dcomplex arg = Complex(t.arg_r, t.arg_i);
      dcomplex res = (t.f)(arg);
      pnl_test_eq (res.r, t.res_r, tol, t.label, 
                  "computed at (%g,%g)", t.arg_r, t.arg_i);
      pnl_test_eq (res.i, t.res_i, tol, t.label, 
                  "computed at (%g,%g)", t.arg_r, t.arg_i);
    }
  exit(pnl_test_finalize("Complex Functions"));
}
