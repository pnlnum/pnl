
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
#include "pnl/pnl_specfun.h"


extern double hyp2f0 (double, double, double, int, double*);

/**
 * Incomplete Gamma function
 *
 * @param a the power of the Gamma function
 * @param x the lower bound of the integral
 * @return Gamma (a, x)
 */
double pnl_sf_gamma_inc (double a, double x)
{
  if (x == 0.)
    {
      return pnl_sf_gamma (a);
    }
  else if (a < 0)
    {
      double shift = -pow (x, a) * exp (-x);
      return (shift +  pnl_sf_gamma_inc (a + 1., x)) / a;
    }
  else if (a == 0.)
    {
      return pnl_sf_expint_En (1, x);
    }
  else
    {
      return  pnl_sf_gamma_inc_Q (a, x) * pnl_sf_gamma (a);      
    }
}

/* double pnl_sf_expint_Ei (double x)
 * {
 *   return - pnl_sf_expint_En (1, -x);
 * } */

double pnl_sf_log_erf (double x)
{
  return log (pnl_sf_erf (x));
}

double pnl_sf_log_erfc (double x)
{
  return log (pnl_sf_erfc (x));
}

/*
 * This a wrapper to hyp2f0 from the Cephes library, but it does not return
 * the same result as GSL which uses the relation between 2F0 and U as
 * implemented in the function pnl_sf_hyperg_2F0
 */

/* double pnl_sf_hyperg_2F0 (double a, double b, double x)
 * {
 *   int type;
 *   double err1, res1;
 *   type = 1;
 *   res1 = hyp2f0 (a, b, x, type, &err1);
 *   return res1;
 * } */

double pnl_sf_hyperg_2F0 (double a, double b, double x)
{
  if (x < 0.)
    {
      double mult;
      mult = pow (-1. / x, a);
      return mult * pnl_sf_hyperg_U (a, 1. + a - b, -1 / x);
    }
  else if (x == 0.)
    {
      return 1.;
    }
  else
    {
      PNL_ERROR ("x must be <= 0", "hyperg_2F0");
    }
}

double pnl_sf_hyperg_0F1 (double v, double x)
{

  if ( x==0 ) return 1.0;

  if ( x < 0 )
    {
      double lgam, bes, m;
      int sign;
      pnl_sf_log_gamma_sgn(v, &lgam, &sign);
      bes = pnl_bessel_j (v-1, 2*sqrt(-x));
      m = log(-x) * 0.5 * (1.0 - v);
      return sign * bes * exp (lgam + m); 

    }
  else
    {
      double lgam, bes, m;
      int sign;
      pnl_sf_log_gamma_sgn(v, &lgam, &sign);
      bes = pnl_bessel_i (v-1, 2*sqrt(x));
      m = log(x) * 0.5 * (1.0 - v);
      return sign * bes * exp (lgam + m); 
    }
}

