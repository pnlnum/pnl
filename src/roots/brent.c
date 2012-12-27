#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_root.h"


/************************************************************************/
/* Copyright PNL developpment team                                      */
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

/* from NETLIB c/brent.shar */

/*************************************************************************
 *			    C math library
 * function ZEROIN - obtain a function zero within the given range
 *
 * Input
 *	double zeroin(ax,bx,f,info,Tol,Maxit)
 *	double ax;			Root will be seeked for within
 *	double bx;			a range [ax,bx]
 *	double (*f)(double x, void *info); Name of the function whose zero
 *					will be seeked for
 *	void *info;			Add'l info passed to f
 *	double *Tol;			Acceptable tolerance for the root
 *					value.
 *					May be specified as 0.0 to cause
 *					the program to find the root as
 *					accurate as possible
 *
 *	int *Maxit;			Max. iterations
 *
 *
 * Output
 *	Zeroin returns an estimate for the root with accuracy
 *	4*EPSILON*abs(x) + tol
 *	*Tol returns estimated precision
 *	*Maxit returns actual # of iterations, or -1 if maxit was
 *      reached without convergence.
 *
 * Algorithm
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *	computations. M., Mir, 1980, p.180 of the Russian edition
 *
 *	The function makes use of the bisection procedure combined with
 *	the linear or quadric inverse interpolation.
 *	At every step program operates on three abscissae - a, b, and c.
 *	b - the last and the best approximation to the root
 *	a - the last but one approximation
 *	c - the last but one or even earlier approximation than a that
 *		1) |f(b)| <= |f(c)|
 *		2) f(b) and f(c) have opposite signs, i.e. b and c confine
 *		   the root
 *	At every step Zeroin selects one of the two new approximations, the
 *	former being obtained by the bisection procedure and the latter
 *	resulting in the interpolation (if a,b, and c are all different
 *	the quadric interpolation is utilized, otherwise the linear one).
 *	If the latter (i.e. obtained by the interpolation) point is
 *	reasonable (i.e. lies within the current interval [b,c] not being
 *	too close to the boundaries) it is accepted. The bisection result
 *	is used in the other case. Therefore, the range of uncertainty is
 *	ensured to be reduced at least by the factor 1.6
 *
 ************************************************************************
 */


/**
 * Search for the root of F in the interval [x1, x2] with a
 * given tolerance
 *
 * @param F a PnlFunc wrapper for the true function under investigation
 * @param x1 lower bound for the root
 * @param x2 upper bound for the root
 * @param tol on input gives the tolerance below with the algorithm stops. On
 * exit contains the maximum error made on the root.
 */
double pnl_root_brent( PnlFunc * F, double  x1, double  x2, double *tol) 
{ 
  int iter, ITMAX=100; 
  double a, b, c, prev_step, new_step;
  double fa, fb, fc;
  double p, q, r, s, tol1;

  a = x1;
  b = x2;
  c = a;
  fa=PNL_EVAL_FUNC(F,a);
  fb=PNL_EVAL_FUNC(F,b);
  fc=fa; 
  prev_step = b - a;
  
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    {
      PNL_ERROR("Roots must be bracketed", "pnl_root_brent"); 
    }

  /* Test if one the endpoints is the root */
  if (fa == 0.)
    {
      *tol = 0.;
      return a;
    }
  if (fb == 0.)
    {
      *tol = 0.;
      return b;
    }
  
  for (iter=1;iter<=ITMAX;iter++)
    {
      prev_step =  b - a;
      if (fabs(fc) < fabs(fb))
        { 
          a = b; b = c; c = a; 
          fa = fb; fb = fc; fc = fa; 
        }
      
      tol1 = 2.0 * DBL_EPSILON * fabs(b) + 0.5 * *tol; 
      new_step = 0.5 * (c - b); 
      if (fabs (new_step) <= tol1 || fb == 0.0)
        {
          *tol = fabs (c - b);
          return b; 
        }
      /* decide if the interpolation can be tried. if prev_step was
         large enough and in the right direction */
      if (fabs (prev_step) >= tol1 && fabs(fa) > fabs(fb))
        {
          
          s = fb / fa; 
          if (a == c) 
            {
              /* if we only have two distinct points, only linear
                 interpolation can be applied */
              p = 2.0 * new_step * s; 
              q = 1.0 - s; 
            }
          else
            {
              /* Quadratic inverse interpolation */
              q = fa / fc; 
              r = fb / fc; 
              p = s * (2.0 * new_step * q * (q - r) - (b - a) * (r - 1.0)); 
              q = (q - 1.0) * (r - 1.0) * (s - 1.0); 
            }
          /* p was calculated with the oppposite sign; make p positive and
             assign the possible minus to q */
          if (p > 0.0)
            {
              q = -q; 
            }
          else
            {
              p = -p;
            }

          /* if b+p/q falls in [b,c] and isn't too large, it is accepted. If
          p/q is too large the the bisection procedure can reduce [b,c] more
          significantly */
          if (2.0 * p < 3.0 * new_step * q - fabs (tol1 * q)
              && 2. * p < fabs (prev_step * q))
            { 
              new_step = p / q; 
            }
          else
            {
              new_step = prev_step = 0.5 * (c - b); 
            }
        }

      /* adjust the step to be not less than tolerance */
      if (fabs(new_step) < tol1)
        {
          new_step = (new_step > 0 ? tol1 : -tol1);
        }

      a = b; 
      fa = fb;
      b += new_step;
      fb=PNL_EVAL_FUNC(F,b);

      /* adjust c to have the opposite sign of b */
      if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0))
        {
          c = a;
          fc = fa;          
        }
    }
  *tol = fabs(c-b);
  return b;
    
} 
