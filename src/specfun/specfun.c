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


#include <errno.h>
#include <math.h>
#include "pnl/pnl_complex.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_specfun.h"

#define MAX_INT_FACT 170

/* 
 * This code originally comes from specfun.f written and copyrighted by 
 *      Shanjie Zhang and Jianming Jin
 *
 * but permission is granted to use code in programs
 * "Computation of Special Functions", 1996, John Wiley & Sons, Inc.
 *
 * The fortran code has been hand translated into C by Jérôme Lelong
 */
/**
 * Compute the Gamma function for complex argument
 *
 * gr + gi i = Gamma(x + i y) if kf = 1
 * gr + gi i = log(Gamma(x + i y)) if kf = 0
 *
 * @param x real part of the argument
 * @param y imaginary of the argument
 * @param gr holds the real part of the output
 * @param gi holds the imaginary part of the output
 * @param kf an integer flag. If kf==1, Gamma is computed, if kf==0, log(Gamma)
 * is computed
 */
static void sp_cgamma (double x, double y, int kf, double *gr, double *gi)
{
  double x1, _y1, x0, z1, z2, th, gr1, gi1, t;
  double sr, si, th1, th2, g0;
  double tmp;
  int na, k, j;
  double a[10] = 
    { .08333333333333333, -2.777777777777778e-3,
      7.936507936507937e-4, -5.952380952380952e-4,
      8.417508417508418e-4, -1.917526917526918e-3,
      6.410256410256410e-3, -.02955065359477124,
      .1796443723688307, -1.39243221690590 };

  if (y==0. && x == (double) ((int) x) && x <= 0.)
    {
      *gr = 1.0e300;
      *gi = 0.;
      return;
    } 
  if (x < 0.)
    {
      x1 = x;
      _y1 = y;
      x = -x;
      y = -y;
    } 
  else
    {
      _y1 = 0.;
      x1 = x;
    } 
  x0 = x;
  na = 0;
  if (x < 7.)
    {
      na = (int) (7 - x); 
      x0 = x + na;
    } 
  z1 = sqrt (x0 * x0 + y * y);
  th = atan (y / x0);
  *gr = (x0 - .5) * log (z1) - th * y - x0 + .5 * log (M_2PI);
  *gi = th * (x0 - .5) + y * log (z1) - y;

  for ( k=1 ; k<=10 ; k++ )
    {
      t = pow (z1, 1- 2 * k);
      *gr += a[k-1] * t * cos ((2. * k - 1.) * th);
      *gi -= a[k-1] * t * sin ((2. * k - 1.) * th);
    } 
  if (x <= 7.)
    {
      gr1 = 0.;
      gi1 = 0.;
      for ( j=0 ; j <= na - 1 ; j++ )
        {
          tmp = x + j;
          gr1 += .5 * log ( tmp * tmp + y * y); 
          gi1 += atan (y / (x + j));
        } 
      *gr -= gr1;
      *gi -= gi1;
    } 
  if (x1 < 0.)
    {
      z1 = sqrt (x * x + y * y);
      th1 = atan (y / x);
      sr = -sin (M_PI * x) * cosh (M_PI * y);
      si = -cos (M_PI * x) * sinh (M_PI * y);
      z2 = sqrt (sr * sr + si * si);
      th2 = atan (si / sr);
      if (sr < 0.) th2 += M_PI;
      *gr = log (M_PI / (z1 * z2)) - *gr;
      *gi = -th1 - th2 - *gi;
      x = x1;
      y = _y1;
    } 
  if (kf == 1)
    {
      g0 = exp (*gr);
      *gr = g0 * cos (*gi);
      *gi = g0 * sin (*gi);
    } 
} 

/**
 * Gamma(z), the Gamma function
 * @param z  a complex number
 * @return  Gamma(z), the Gamma function of z
 */
dcomplex Ctgamma (dcomplex z)
{
  dcomplex c;
  int kf = 1;
  sp_cgamma ((z.r), (z.i), kf, &(c.r), &(c.i));
  return c;
} 

/**
 * log(Gamma(z)), the logarithm of the Gamma function
 * @param z  a complex number
 * @return  log(Gamma(z))
 */
dcomplex Clgamma (dcomplex z)
{
  dcomplex c;
  int kf = 0;
  sp_cgamma ((z.r), (z.i), kf, &(c.r), &(c.i));
  return c;
} 

/** factorial function
 * @param n an integer
 * @return the factorial of n as a double
 */
double pnl_sf_fact (int n)  
{
  if ( n <= 0 ) return 1.0;
  if ( n >= MAX_INT_FACT ) return PNL_INF;
  else
    { 
      int i;
      double z = 1.;
      for ( i=n ; i>1 ; i-- ) z *= i;
      return z;
    }
}

/** 
 * Compute the binomial coefficient 
 *
 *      / n \        n!      
 *      |   |  = --------- 
 *      \ k /    k! (n-k)! 
 *
 * @param n a non-negative integer
 * @param p a non-negative integer smaller that n
 * 
 * @return 
 */
double pnl_sf_choose (int n, int p)
{
  int k = MAX(p, n-p);

  if ( n-p<0 || n<0 || p<0 ) return 0.;

  if ( k < 64 )
    {
      int i;
      double prod = 1.;
      for ( i=n ; i>k ; i-- )
        {
          prod *= (double) (i) / (double) (i-k);
        }
      return prod;
    }
  else
    {
      double ln_choose = pnl_sf_log_gamma (n+1) - pnl_sf_log_gamma (p+1) 
        - pnl_sf_log_gamma (n-p+1);
      return exp (ln_choose);
    }
}
