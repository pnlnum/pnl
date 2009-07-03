
/*************************************************************************/
/* Written and (C) by the PNL core team                                  */
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
#include "pnl_integration.h"

#include "tests.h"

static double x_square(double x, void *p)
{ return tan(x);}

static double test_xy(double x,double y, void *p)
{ return cos(x+y);}


static void integration_GK_test()
{
  double result,abserr;
  int neval;
  PnlFunc func;
  func.function = x_square;
  func.params = NULL;
  pnl_integration_GK(&func,0.0,1.0,0.00001,0.000001,&result,&abserr,&neval);
  printf(" integration res %f error %f , iter %d \n",result,abserr,neval); 
}

static void integration_GK2D_test()
{
  double result,abserr;
  int neval;
  PnlFunc2D func;
  func.function = test_xy;
  func.params = NULL;
  pnl_integration_GK2D(&func,-1.0,1.0,-1.0,1.0,0.00001,0.000001,&result,&abserr,&neval);
  printf(" integration 2D res %f error %f , iter %d \n",result,abserr,neval); 
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

static void integ_2D_test()
{
  double result;
  int nx, ny;
  PnlFunc2D func;
  nx = 50;
  ny = 50;
  func.function = test_xy;
  func.params = NULL;
  result = pnl_integration_2D(&func,-1.0,1.0,-1.0,1.0,nx, ny, "rect");
  printf(" integration 2D (rectangle rule) : %f, iter %d \n", result, nx * ny); 
  result = pnl_integration_2D(&func,-1.0,1.0,-1.0,1.0,nx, ny, "trap");
  printf(" integration 2D (trapezoidal rule) : %f, iter %d \n", result, nx *ny); 
  result = pnl_integration_2D(&func,-1.0,1.0,-1.0,1.0,nx, ny, "simpson");
  printf(" integration 2D (simpson rule) : %f, iter %d \n", result, nx *ny); 
}

void integration_test()
{
  integration_GK_test();
  integration_GK2D_test();
  integ_test();
  integ_2D_test();
}
