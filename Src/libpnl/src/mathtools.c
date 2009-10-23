/*************************************************************************/
/* Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>            */
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

#include <math.h>

#include "config.h"
#include "pnl_mathtools.h"


/** nearest integer round off function
 * @param s the value to be rounded
 * @return the rounded value casted to an integer. 1.5 is rounded to 1 and
 * -1.5 is rounded to -1
 */
int intapprox(double s)
{
  int r = (int) s;
  if (s > r + 0.5) { r++; }
  if (s < r - 0.5) { r--; }
  return r;
}


#ifndef HAVE_TRUNC
double trunc(double x)
{
  return (double) ((int) x);
}
#endif

/** Computes the binomial coefficients with an error smaller than 1 for large
 * values of the parameter n
 * @param n an integer
 * @param p an integer
 * @return C(n, p)
 */
double Cnp(int n, int p) 
{                          
  
  double z, iter;
  int i;
  z=0.0;

  if ((n-p<= -1) || (n<0) || (p<0))
  {
    return z;
  }
  else
  {
    if (p==0)
      z=1.0;
    else
    {
      z=1.0;
      iter=z;
      i=0;
      while(i<=p-1)
      {
        iter=iter*(double)(n-i)/(p-i);
        i=i+1;
      }
      if ((iter-floor(iter))>=0.5)
        z=floor(iter)+1;
      else
        z=floor(iter);
    }
  }

  return z;
}

/** factorial function
 * @param n an integer
 * @return the factorial of n as a double
 */
double pnl_fact(int n)  
{
  int i;
  double z=1.0;
  if (n<=0) { z=1.0; }
  else { for ( i=1 ; i<=n ; i++) z=z*i; }
  return z;
}


#ifndef HAVE_LGAMMA
extern double amos_dgamln (double *z__, int *ierr);
/**
  * Computes the natural logarithm of the Gamma function
  * @param x  a real number > 0
  * @return  log(Gamma(x))
  */ 

double lgamma(double x)
{
  int ierr;
  if ( x <= 0 ) return NAN;
  return amos_dgamln (&x, &ierr);
}
#endif

/* the tgamma function is part of C99, not available under
 * Windows */
#ifndef HAVE_TGAMMA
/**
  * Computes the true Gamma function 
  * @param x  a real number > -1
  * @return  Gamma(x)
  */
double tgamma(double x)
{
  return (x > 0) ? exp(lgamma(x)) : M_PI/tgamma(1-x)/sin(M_PI*(1-x));
} 
#endif

/**
  * Interger power function
  * @param x a real number
  * @param n an integer
  * @ return x^n
  */
double pnl_pow_i (double x, int n)
{
  double y = 1.;

  if ( x == 0 && n > 0 ) return 0.;

  if ( n < 0 ) 
    {
      x = 1. / x;
      n = -n;
    }

  while (n)
    {
       if (PNL_IS_ODD (n)) y *= x; /* y = x * x^(n-1) */
       n >>= 1;
       x *= x;
    } 
  return y;
} 
