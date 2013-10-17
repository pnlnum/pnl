
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

#include <string.h>
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_integration.h"

/**
 * Evaluate an integral over [x0, x1] with a non-adaptive method
 *
 * @param *F a PnlFunc to be integrated 
 * @param x0 left corner of domain
 * @param x1 rightt corner of domain
 * @param n the number of steps to be used
 * @param meth a character string can be "rect" (rectangle rule), "trap"
 * (trapezoidal rule), "simpson" (Simpson's rule)
 * @return the value of the integral
 */
double pnl_integration (PnlFunc *F, double x0, double x1, int n, char *meth)
{
  double h, xh, sum;

  h = (x1 - x0) / (double)n;
  sum = 0.;

  /* sum over the strictly inner points */
  for (xh = x0 + h; xh < x1; xh += h) { sum += PNL_EVAL_FUNC(F, xh); }

  if (strcmp (meth, "rect") == 0) { sum += PNL_EVAL_FUNC(F, x0); }
  else if (strcmp (meth, "trap") == 0) { sum += (PNL_EVAL_FUNC(F, x0) + PNL_EVAL_FUNC(F, x1)) / 2.; }
  else if (strcmp (meth, "simpson") == 0)
    {
      double sum2 = 0.;
      for (xh = x0 + h/2; xh < x1; xh += h) { sum2 += PNL_EVAL_FUNC(F, xh); }
      sum = 2. * sum + PNL_EVAL_FUNC(F, x0) + PNL_EVAL_FUNC(F, x1) + 4. * sum2;
      sum /= 6.;
    }
  else { PNL_ERROR("unknow method", "pnl_integration"); }

  sum *= h;
  return sum;
}


static double xsav, y0sav, y1sav;
static int nysav;
static char *methsav;
static const PnlFunc2D* globalfunc;

static double func1D(double y, void *params) 
{
  return globalfunc->F(xsav, y, params);
}

static double int_1d(double x, void *params) 
{
  PnlFunc wrap_1d;
  wrap_1d.F = func1D;
  wrap_1d.params = params;
  xsav=x;
  return pnl_integration (&wrap_1d, y0sav, y1sav, nysav, methsav);
}
/**
 * integrate a function on [x0,x1]x[y0, y1] 
 * @param *F a PnlFunc of two variables to be integrated 
 * @param x0 bottom left corner of domain
 * @param x1 bottom right corner of domain
 * @param y0 top left corner of domain
 * @param y1 top right corner of domain
 * @param nx the number of steps to discretise [x0, x1]
 * @param ny the number of steps to discretise [y0, y1]
 * @param meth a character string can be "rect" (rectangle rule), "trap"
 * (rapezoidal rule), "simpson" (Simpson's rule)
 * @return the value of the integral
 */
double pnl_integration_2d (PnlFunc2D *F, double x0, double x1,
                          double y0,double y1, int nx, int ny, char *meth)
{
  PnlFunc func_1d;
  func_1d.F = &int_1d;
  func_1d.params   = NULL;
  globalfunc       = F;
  y0sav            = y0;
  y1sav            = y1;
  nysav            = ny;
  methsav          = meth;
  return pnl_integration (&func_1d, x0, x1, nx, meth);
}
