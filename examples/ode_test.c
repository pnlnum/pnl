
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
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_integration.h"

#include "tests_utils.h"


/*
 * Computes  0.25 * Y * ( 1 - Y / 20 )
 */
static void r4_f1 (int neqn, double t, const double *y, double *yp, void *param)
{
  yp[0] = 0.25 * y[0] * ( 1.0 - y[0] / 20.0 );
}

/*
 * Computes the exact solution to the ODE associated to r4_f1
 */
static double r4_y1x (double t)
{
  double value;
  value = 20.0 / ( 1.0 + 19.0 * exp ( -0.25 * t ) );
  return value;
}

void test01 ( )
{
# define NEQN 1
  int flag, i, n_step;
  double abserr, relerr;
  double t, t_out, t_start, t_stop;
  double y[NEQN], yp[NEQN];
  PnlODEFunc f;

  abserr = 0.;
  relerr = 0.;

  t_start = 0.0;
  t_stop = 20.0;

  n_step = 10;

  t = 0.0;
  t_out = 0.0;
  y[0] = 1.0;
  r4_f1 ( 1, t, y, yp, NULL );

  f.F = r4_f1;
  f.neqn = 1;
  f.params = NULL;

  for ( i = 0; i < n_step; i++ )
    {
      t = ((n_step - i) * t_start + i * t_stop) / ( double ) n_step;
      t_out = ((n_step - (i + 1) ) * t_start + (i + 1) * t_stop)  / ( double ) n_step;
      if ( pnl_ode_rkf45 ( &f, y, t, t_out, relerr, abserr, &flag ) == FAIL )
        {
          pnl_test_set_fail ("pnl_ode_rkf45 failed", 0., 0.);
        }
      pnl_test_eq_abs (y[0], r4_y1x (t_out), 1E-5, "pnl_ode_rkf45", "at time %t", t_out);
    }
# undef NEQN
}



int main (int argc, char **argv)
{
  pnl_test_init (argc, argv);
  test01 ();
  exit (pnl_test_finalize ("ODE"));
}
