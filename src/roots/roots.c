
/************************************************************************/
/* Copyright Jérôme Lelong <jerome.lelong@gmail.com>                    */
/* Copyright David Pommier <pommier.david@gmail.com>                    */
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


#include <math.h>
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_root.h"


/**
 * Find the root of a function
 *
 * @param x_min a double pointer to an already existing PnlVect
 * @param x_max a double the root is suppose to be in [x_min, x_max] 
 * @param tol a double speed of derivative decrease, if derivative is less
 * than tolerance, then it is converged 
 * @param N_max a int maximal number of iteration
 * @param res the root, if the algorithm has converged
 * @param func a function pointer which computes f(x_anc) df(x_anc) 
 * @return OK or FAIL
 */
int pnl_find_root (PnlFuncDFunc *func, double x_min, double x_max, double tol, int N_max, double * res)
{
  int iter; 
  double func_low, func_high, func_current, diff_func_current;
  double rts, dx_anc, dx,xl, xh;

  PNL_EVAL_FUNC_DFUNC (func, x_min, &func_low, &diff_func_current);
  if (func_low == 0.0)
    {
      *res=x_min;return OK;
    }
  PNL_EVAL_FUNC_DFUNC (func, x_max, &func_high,&diff_func_current);
  if (func_high == 0.0) 
    {
      *res=x_max;return OK;
    }
  if ((func_low > 0.0 && func_high > 0.0) || (func_low < 0.0 && func_high < 0.0)) 
    return FAIL;
  /*Root must be bracketed by x1 and x2  */
  if (func_low < 0.0)
    {
      xl = x_min; xh = x_max;
    }
  else
    {
      xl = x_max; xh = x_min;
    }
  rts = 0.5 * (x_min + x_max);
  dx_anc = fabs(x_max - x_min);
  dx = dx_anc;
  PNL_EVAL_FUNC_DFUNC (func, rts, &func_current, &diff_func_current);
  for (iter=0;iter<N_max;iter++)
    {
      if (((rts - xh) * diff_func_current - func_current) * ((rts - xl) * diff_func_current - func_current) >= 0.0
          || fabs(2.0 * func_current) > fabs(dx_anc * diff_func_current))
        {
          dx_anc = dx;
          dx = 0.5 * (xh - xl);
          rts = xl + dx;
        }
      else
        {
          dx_anc = dx;
          dx = func_current / diff_func_current;
          rts -= dx;
        }
      if (fabs(dx) < tol)
        {
          *res= rts;return OK;
        }
      PNL_EVAL_FUNC_DFUNC (func, rts, &func_current,&diff_func_current);
      if (func_current < 0.0)
        xl = rts;
      else 
        xh= rts;
    }
  return FAIL;
  /* Maximum number of iterations exceeded */
};


/**
 * Find the root of a function using Newton's algorithm
 *
 * @param x0 initial guess
 * @param epsrel if the relative improvement over the root is less than this value,
 * then break;
 * @param epsabs if the absolute improvement over the root is less than this value,
 * then break;
 * @param res contains the root on exit
 * @param N_max maximum number of iterations
 * @param func a F_DF pointer
 * @return OK or FAIL if N_max reached
 */
int pnl_root_newton (PnlFuncDFunc *func, double x0, double epsrel, double epsabs, int N_max, double * res)
{
  int i;
  double f, df, root, root_prev;

  root = x0;

  for (i=0; i<N_max; i++)
    {
      PNL_EVAL_FUNC_DFUNC (func, root, &f, &df);
      if (df == 0.)
        {
          PNL_ERROR ("div by zero", "pnl_root_newton");
        }
      root_prev = root;
      root -= f / df;
      if (fabs (root - root_prev) < epsrel * fabs(root_prev) + epsabs) break;
    }
  *res = root;
  if (i==N_max)
    return FAIL; /* maximum number of iterations reached */
  return OK;
}

/**
 * Find the root of a function using a bisection method
 *
 * @param xmin lower bound
 * @param xmax upper bound
 * @param epsrel if the relative improvement over the root is less than this value,
 * then break;
 * @param epsabs if the absolute improvement over the root is less than this value,
 * then break;
 * @param res contains the root on exit
 * @param N_max maximum number of iterations
 * @param func a function pointer
 * @return OK or FAIL if N_max reached
 */
int pnl_root_bisection (PnlFunc *func, double xmin, double xmax, double epsrel, double epsabs,
                        int N_max, double *res)
{
  int i;
  double a, b, c, fa, fb, fc;

  fa = PNL_EVAL_FUNC (func, xmin);
  fb = PNL_EVAL_FUNC (func, xmax);

  if ((fa < 0. && fb < 0.) || (fa > 0. && fb > 0.))
    {
      PNL_ERROR ("Root must be bracketed", "pnl_root_bisection");
    }

  if (fa > 0.)
    {
      a = xmax; b = xmin;
    }
  else
    {
      a = xmin; b = xmax;
    }
    
  fa = PNL_EVAL_FUNC (func, a);
  fb = PNL_EVAL_FUNC (func, b);

  for (i=0; i<N_max; i++)
    {
      c = (a + b) / 2.;
      fc = PNL_EVAL_FUNC (func, c);

      if (fc < 0.)  { a = c;  fa = fc; }
      else { b = c; fb = fc; }
      if ( fabs (b - a) < epsabs + epsrel * fabs (a)) break;
    }
  *res = (a + b) / 2.;
  if (i==N_max) return FAIL;
  return OK;
}
