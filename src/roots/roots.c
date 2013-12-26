
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
#include "pnl/pnl_matrix.h"


/**
 * Find the root of a function by combining Newton's method with the bisection
 * method
 *
 * @param x_min a double pointer to an already existing PnlVect
 * @param x_max a double the root is suppose to be in [x_min, x_max] 
 * @param tol a double speed of derivative decrease, if derivative is less
 * than tolerance, then it is converged 
 * @param max_iter a int maximal number of iteration
 * @param res the root, if the algorithm has converged
 * @param func a function pointer which computes f(x_anc) df(x_anc) 
 * @return OK or FAIL
 */
int pnl_root_newton_bisection (PnlFuncDFunc *func, double x_min, double x_max, double tol, int max_iter, double *res)
{
  int i; 
  double func_low, func_high, func_current, diff_func_current;
  double rts, dx_anc, dx,xl, xh;

  PNL_EVAL_FUNC_FDF (func, x_min, &func_low, &diff_func_current);
  if (func_low == 0.0)
    {
      *res=x_min; return OK;
    }
  PNL_EVAL_FUNC_FDF (func, x_max, &func_high,&diff_func_current);

  if (func_high == 0.0) 
    {
      *res=x_max; return OK;
    }

  /* Root is not bracketed by x1 and x2  */
  if ( (func_low > 0.0 && func_high > 0.0) 
       || (func_low < 0.0 && func_high < 0.0) ) 
    return FAIL;

  if (func_low < 0.0) { xl = x_min; xh = x_max; }
  else { xl = x_max; xh = x_min; }

  rts = 0.5 * (x_min + x_max);
  dx_anc = fabs(x_max - x_min);
  dx = dx_anc;
  PNL_EVAL_FUNC_FDF (func, rts, &func_current, &diff_func_current);

  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( ( ((rts - xh) * diff_func_current - func_current) 
             * ((rts - xl) * diff_func_current - func_current) >= 0.0 )
          || fabs(2.0 * func_current) > fabs(dx_anc * diff_func_current) )
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
          *res = rts; return OK;
        }
      PNL_EVAL_FUNC_FDF (func, rts, &func_current, &diff_func_current);
      if ( func_current < 0.0 ) xl = rts; else xh= rts;
    }
  /* Maximum number of iterations exceeded */
  return FAIL;
};

/**
 * Find the root of a function using Newton's algorithm with the Armijo line
 * search to ensure the absolute value of the function decreases along the
 * iterations.  The descent direction at step k is given by 
 *
 *     \f[ d_k = f(x_k) / f'(x_k) \f]
 *
 * Determine \f$ \alpha_k = \max\{2^{-j}, j \ge 0\} \f$ s.t. 
 *
 *   \f[  f(x_k + \alpha_k d_k) \le f(x_k) (1 - \omega \alpha_k) \f]
 *
 * where in this implementation \f$ \omega = 10^{-4} \f$.
 *
 * @param x0 initial guess
 * @param x_eps if the relative improvement over the root is less than this value,
 * then stop;
 * @param fx_eps if |f(x)| < fx_eps * then stop;
 * @param res contains the root on exit
 * @param max_iter maximum number of iterations
 * @param func a F_DF pointer
 * @return OK or FAIL if N_max reached
 */
int pnl_root_newton (PnlFuncDFunc *func, double x0, double x_eps, double fx_eps, int max_iter, double *res)
{
  int i;
  double omega, gamma;
  double f, df, root;

  omega = 1E-4;
  gamma = 0.5;
  root = x0;
  PNL_EVAL_FUNC_FDF (func, root, &f, &df);

  for ( i=0; i<max_iter; i++ )
    {
      double dx, t;
      double norm0, norm;
      t = 1.;
      if (df == 0.)
        {
          PNL_ERROR ("div by zero", "pnl_root_newton");
        }
      dx = f /df;
      norm0 = fabs (f);
      /* Armijo line search */
      while ( t > DBL_EPSILON )
        {
          double x_linesearch;
          x_linesearch = root - t * dx;
          PNL_EVAL_FUNC_FDF (func, x_linesearch, &f, &df);
          norm = fabs (f);
          if ( norm < norm0 * (1. - omega * t) )
            {
              root = x_linesearch;
              break;
            }
          t *= gamma;
        }
      if (fabs (dx) < x_eps * fabs(root) || norm < fx_eps) break;
    }
  *res = root;
  /* maximum number of iterations reached */
  if ( i == max_iter ) return FAIL; 
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

/** 
 * Print the state of a Newton iteration
 * 
 * @param iter index of the iteration
 * @param x 
 * @param Fx 
 * @param verbose
 */
static void print_newton_iter (int iter, const PnlVect *x, const PnlVect *Fx, int verbose)
{
  if ( verbose == TRUE )
    {
      printf ("iter = %d\n", iter);
      printf ("\tx = "); pnl_vect_print_asrow (x);
      printf ("\tf(x) = "); pnl_vect_print_asrow (Fx);
      printf ("\t|f(x)| = %f\n", pnl_vect_norm_two (Fx));
    }
}

/** 
 * Run Newton's algorithm to find the zero of a function f using the Armijo line
 * search to make sure the absolute value of the function decreases from one
 * step to the following. The descent direction at step k is given by 
 *
 *     \f[ d_k = (\nabla f(x_k))^{-1} f(x_k) \f]
 *
 * Determine \f$ \alpha_k = \max\{2^{-j}, j \ge 0\} \f$ s.t. 
 *
 *   \f[  f(x_k + \alpha_k d_k) \le f(x_k) (1 - \omega \alpha_k) \f]
 *
 * where in this implementation \f$ \omega = 10^{-4} \f$.
 * 
 * @param func the function of which to find the zero. Note that, the Jacobian
 * function must also be provided
 * @param x0 starting value
 * @param x_eps minimum relative movement of x before the algorithms stops.
 * @param fx_eps the algorithm stops if |f(x)| < fs_eps
 * @param max_iter Maximum number of iterations
 * @param verbose boolean (TRUE or FALSE). If TRUE, the details of each
 * iteration are printed
 * @param res Holds the solution on output
 * 
 * @return 
 */
int pnl_multiroot_newton (PnlRnFuncRnDFunc *func, const PnlVect *x0, double x_eps, double fx_eps, int max_iter, int verbose, PnlVect *res)
{
  int iter, n;
  double gamma, omega;
  PnlMat *DFx;
  PnlVect *dx, *Fx, *x_linesearch;


  n = x0->size;
  omega = 1E-4;
  gamma = 0.5;
  pnl_vect_resize (res, n);
  pnl_vect_clone (res, x0);
  dx = pnl_vect_create (n);
  x_linesearch = pnl_vect_create (n);
  Fx = pnl_vect_create (n);
  DFx = pnl_mat_create (n, n);

  PNL_CHECK ( func->FDF == NULL && (func->F == NULL || func->DF == NULL), 
              "Incomplete function definition", "pnl_multiroot_newton" );
  PNL_CHECK ( x_eps < 0.0 || fx_eps < 0.0, 
              "Negative precision required", "pnl_multiroot_newton" );

  /*
   * Evaluate the function before entering the loop
   */
  PNL_EVAL_RNFUNCRN_F_DF (func, res, Fx, DFx);

  for ( iter=0 ; iter<max_iter ; iter++ )
    {
      int nb_tries_ls = 1;
      double norm, norm0, t;
      t = 1.;

      print_newton_iter (iter, res, Fx, verbose);
      if ( pnl_mat_syslin (dx, DFx, Fx) == FAIL )
        {
          PNL_ERROR ("Non invertible Jacobian", "pnl_multiroot_newton");
        }
      norm0 = pnl_vect_norm_two (Fx);
      /* Armijo line search */
      while ( t > DBL_EPSILON )
        {
          pnl_vect_clone (x_linesearch, res);
          pnl_vect_axpby (-t , dx, 1., x_linesearch); 
          PNL_EVAL_RNFUNCRN_F_DF (func, x_linesearch, Fx, DFx);
          norm = pnl_vect_norm_two (Fx);
          if ( norm < norm0 * (1. - omega * t) )
            {
              pnl_vect_clone (res, x_linesearch);
              break;
            }
          t *= gamma;
          nb_tries_ls++;
        }
      if ( (verbose == TRUE) && (nb_tries_ls > 1) )
        {
          printf ("\tNeed to redude step size %d times.\n", nb_tries_ls);
        }
      if ( (pnl_vect_norm_two (dx) < x_eps * pnl_vect_norm_two (res))
           || (norm < fx_eps) )
        break;
    }
  print_newton_iter (iter, res, Fx, verbose);

  pnl_vect_free (&dx);
  pnl_vect_free (&Fx);
  pnl_mat_free (&DFx);

  /* maximum number of iterations reached */
  if ( iter == max_iter ) return FAIL; 
  return OK;
}
