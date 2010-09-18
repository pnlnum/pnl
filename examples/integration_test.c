
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

#include "tests.h"

static double x_square(double x, void *p)
{ return tan(x) ;}

static double test_xy(double x,double y, void *p)
{ return cos(x+y);}

/*
 * Function used to test integration of a function with singular points 
 * The integral over (0,10) is equal to 1 - cos(2.5) + exp(-2.5) - exp(-5.0)
 */
static double singular (double x, void *p)
{
  if ( x > 0 && x < 2.5 ) return sin(x);
  else if ( x >= 2.5 && x <= 5. ) return exp(-x);
  else return 0.0;
}

/*
 * Function used to test integration over a non finite interval (0, +oo)
 * Its integral over (0 +oo) is the Euler Gamma constant
 */
static double indefinite (double x, void *p)
{
  return -exp(-x) * log(x);
}

static void integration_qag_test()
{
  double result,abserr;
  int neval;
  PnlFunc func;
  func.function = x_square;
  func.params = NULL;
  pnl_integration_qag(&func,0.0,1.0,0.00001,0.000001,0,&result,&abserr,&neval);
  printf(" integration (QAGS finite interval with no singularity )\n\tres %f error %f , iter %d \n",result,abserr,neval); 

  func.function = indefinite;
  func.params = NULL;
  pnl_integration_qag(&func,0,PNL_POSINF,0.00001,0.000001,0,&result,&abserr,&neval);
  printf(" integration (QAGI infinite interval) \n\tres %f true value %f , iter %d \n",result,M_EULER,neval); 

  func.function = singular;
  func.params = NULL;
  pnl_integration_qag(&func,0.0,10,0.00001,0.000001,0,&result,&abserr,&neval);
  printf(" integration (QAGS with singularities) \n\tres %f true value %f , iter %d \n",
         result, 1 - cos(2.5) + exp(-2.5) - exp(-5.0),neval); 
}

/*
 * Test of integration routines for functions with known singular points
 */
static void integration_qagp_test()
{
  double result,abserr;
  PnlVect *pts;
  int neval;
  PnlFunc func;

  pts = pnl_vect_create_from_list (2, 2.5, 5.);
  func.function = singular;
  func.params = NULL;
  pnl_integration_qagp(&func,0.0,10,pts,0.00001,0.000001,0,&result,&abserr,&neval);
  printf(" integration (QAGP known singularities) \n\tres %f true value %f , iter %d \n",
         result, 1 - cos(2.5) + exp(-2.5) - exp(-5.0),neval); 
  pnl_vect_free (&pts);
}
static void integration_qng_test()
{
  double result,abserr;
  int neval;
  PnlFunc func;
  func.function = x_square;
  func.params = NULL;
  pnl_integration_qng(&func,0.0,1.0,0.00001,0.000001,&result,&abserr,&neval);
  printf(" integration (qng) res %f error %f , iter %d \n",result,abserr,neval); 
}

static void integration_qng_2d_test()
{
  double result,abserr;
  int neval;
  PnlFunc2D func;
  func.function = test_xy;
  func.params = NULL;
  pnl_integration_qng_2d(&func,-1.0,1.0,-1.0,1.0,0.00001,0.000001,&result,&abserr,&neval);
  printf(" integration 2D (qng) res %f error %f , iter %d \n",result,abserr,neval); 
}


static void integ_test()
{
  double result;
  PnlFunc func;
  int n = 50;
  func.function = x_square;
  func.params = NULL;
  result = pnl_integration(&func,0.0,1.0,n, "rect");
  printf(" integration (rectangle rule) : %f, iter %d \n", result, n); 
  result = pnl_integration(&func,0.0,1.0,n, "trap");
  printf(" integration (trapezoidal rule) : %f, iter %d \n", result, n); 
  result = pnl_integration(&func,0.0,1.0,n, "simpson");
  printf(" integration (simpson rule) : %f, iter %d \n", result, n); 
}

static void integ_2d_test()
{
  double result;
  int nx, ny;
  PnlFunc2D func;
  nx = 50;
  ny = 50;
  func.function = test_xy;
  func.params = NULL;
  result = pnl_integration_2d(&func,-1.0,1.0,-1.0,1.0,nx, ny, "rect");
  printf(" integration 2D (rectangle rule) : %f, iter %d \n", result, nx * ny); 
  result = pnl_integration_2d(&func,-1.0,1.0,-1.0,1.0,nx, ny, "trap");
  printf(" integration 2D (trapezoidal rule) : %f, iter %d \n", result, nx *ny); 
  result = pnl_integration_2d(&func,-1.0,1.0,-1.0,1.0,nx, ny, "simpson");
  printf(" integration 2D (simpson rule) : %f, iter %d \n", result, nx *ny); 
}

void integration_test()
{
  integration_qag_test();
  integration_qagp_test();
  integration_qng_test();
  integration_qng_2d_test();
  integ_test();
  integ_2d_test();
}
