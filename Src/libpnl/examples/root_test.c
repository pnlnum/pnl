
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
#include "pnl_root.h"

#include "tests.h"

static double epsabs = 0.0001;
static double epsrel = 0.00001;
static int N_max = 1000;

#define FUNC  f_cos
#define FDF_FUNC fdf_cos 

static double x_square(double x, void *p) { return x*x -2.; }

static void x_fdf_square(double x, double *f, double *df, void *p)
{
  *f = x*x -2.;
  *df = 2*x;
}

static double f_cos (double x, void *p) {return cos (x);}
static void fdf_cos (double x, double *f, double *df, void *p) {*f = cos (x); *df = -sin(x);}

static void bisection_test ()
{
  double x1, x2, r;
  PnlFunc func;
  int status;
  
  x1 = 0;
  x2 = 3;
  func.function = FUNC;
  func.params = NULL;

  status = pnl_root_bisection(&func, x1, x2, epsrel, epsabs, N_max, &r);
  printf("bisection : root is %f, status==OK %d\n", r, status==OK);

}

static void newton_test ()
{
  double x0, r;
  PnlFuncDFunc func;
  int status;
  
  x0 = 0.5;
  func.function = FDF_FUNC;
  func.params = NULL;

  status = pnl_root_newton(&func, x0, epsrel, epsabs, N_max, &r);
  printf("newton : root is %f, status==OK %d\n", r, status==OK);
}

static void brent_test()
{
  double x1, x2, tol, r;
  PnlFunc func;

  x1 = -1;
  x2 = 3;
  tol = 0.001;
  func.function = FUNC;
  func.params = NULL;

  r = pnl_root_brent(&func, x1, x2, &tol);
  printf("brent : root is %f (tol = %f) \n", r, tol);
}

static void find_root_test ()
{
  double x1, x2, tol, r;
  PnlFuncDFunc func;
  int status;
  
  x1 = -1;
  x2 = 3;
  tol = 0.001;
  func.function = FDF_FUNC;
  func.params = NULL;

  status = pnl_find_root(&func, x1, x2, tol, N_max, &r);
  printf("find_root : root is %f, status==OK %d\n", r, status==OK);
}

void root_test()
{
  brent_test () ;
  bisection_test ();
  newton_test ();
  find_root_test ();
}
