
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pnl/pnl_complex.h"
#include "pnl/pnl_mathtools.h"


/**
 * Re(z) the real part
 * @param z  a complex number
 * @return  the real part 
 */
double Creal( dcomplex z ) { return z.r;}

/**
 *  Im(z) imazinary part 
 * @param z  a complex number
 * @return  the imaginary part  Im(z)
 */
double Cimag( dcomplex z ) {return z.i;}

/**
 * a+b 
 * @param a  a complex number
 * @param b  a complex number 
 * @return the sum:  a+b
 */
dcomplex Cadd(dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r=a.r+b.r;
  c.i=a.i+b.i;
  return c;
}


/**
 * a+b 
 * @param z  a complex number
 * @param x  a real number 
 * @return the sum:  z + x
 */
dcomplex CRadd(dcomplex z, double x)
{
  dcomplex c;
  c.r=z.r+x;
  c.i=z.i;
  return c;
}

/**
 * a+b 
 * @param z  a complex number
 * @param x  a real number 
 * @return the sum:  z + x
 */
dcomplex RCadd(double x, dcomplex z)
{
  dcomplex c;
  c.r=z.r+x;
  c.i=z.i;
  return c;
}


/**
 *  a-b 
 * @param a  a complex number
 * @param b  a complex number 
 * @return  a-b , substraction
 */
dcomplex Csub(dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r=a.r-b.r;
  c.i=a.i-b.i;
  return c;
}

/**
 *  a-b 
 * @param a  a complex number
 * @param b  a real number 
 * @return  a-b , substraction
 */
dcomplex CRsub(dcomplex a, double b)
{
  dcomplex c;
  c.r=a.r-b;
  c.i=a.i;
  return c;
}

/**
 *  a-b 
 * @param a  a complex number
 * @param b  a real number 
 * @return  a-b , substraction
 */
dcomplex RCsub(double a, dcomplex b)
{
  dcomplex c;
  c.r=a-b.r;
  c.i=-b.i;
  return c;
}

/**
 *  -z 
 * @param z  a complex number
 * @return  -z, unary minus
 */
dcomplex Cminus(dcomplex z)
{
  dcomplex c;
  c.r = - z.r;
  c.i = - z.i;
  return c;
}


/**
 *  a*b 
 * @param a  a complex number
 * @param b  a complex number 
 * @return  a*b , the product
 */
dcomplex Cmul(dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r=a.r*b.r-a.i*b.i;
  c.i=a.i*b.r+a.r*b.i;
  return c;
}


/**
 *  multiplication by a real number
 * @param a  a complex number
 * @param x  a double
 * @return  x*a
 */
dcomplex RCmul(double x, dcomplex a)
{
  dcomplex c;
  c.r=x*a.r;
  c.i=x*a.i;
  return c;
}

/**
 *  multiplication by a real number
 * @param a  a complex number
 * @param x  a double
 * @return  x*a
 */
dcomplex CRmul(dcomplex a, double x)
{
  dcomplex c;
  c.r=x*a.r;
  c.i=x*a.i;
  return c;
}

/**
 * constructor
 * @param re real part
 * @param im imaginary part
 * @return  re + i im  
 */
dcomplex Complex(double re, double im)
{
  dcomplex c;
  c.r=re;
  c.i=im;
  return c;
}

/**
 * constructor  r exp(i theta) 
 * @param r : module 
 * @param theta : angle  
 * @return  r exp(i theta)
 */
dcomplex Complex_polar(double r, double theta)
{
  dcomplex c;
  if (r<0)
    {
      PNL_ERROR ("rho < 0", "Complex_polar");
    }
  if (r == 0.) return CZERO;
  c.r=r*cos(theta);
  c.i=r*sin(theta);
  return c;
}

/**
 * Print a complex numbers
 *
 * @param z a complex number
 */
void Cprintf( dcomplex z)
{
  printf("%f + %f i", z.r, z.i);
}

/**
 * conjugate operator  conf(z)
 * @param z a complex number
 * @return  the conjugate of z
 */
dcomplex Conj(dcomplex z)
{
  dcomplex c;
  c.r = z.r;
  c.i = -z.i;
  return c;
}

/**
 *  1/a
 * @param a  a complex number
 * @return  1/a, the inverse
 */
dcomplex Cinv(dcomplex a)
{
  return RCdiv (1., a);
}

/**
 *  a/b
 * @param a  a complex number
 * @param b  a complex number
 * @return a/b, the division
 */
dcomplex Cdiv(dcomplex a, dcomplex b)
{
  dcomplex c;
  double r,d;
  if (fabs (b.r) >= fabs (b.i))
    {
      r = b.i / b.r;
      d = b.r * (1. + r * r);
      c.r = (a.r + r * a.i) / d;
      c.i = (a.i- r * a.r) / d;
    }
  else
    {
      r = b.r / b.i;
      d = b.i * (1. + r * r);
      c.r = (a.r * r + a.i) / d;
      c.i = (a.i * r - a.r) / d;
    }
  return c;
}

/**
 *  a/b
 * @param a  a real number
 * @param b  a complex number
 * @return a/b, the division
 */
dcomplex RCdiv(double a, dcomplex b)
{
  dcomplex c;
  double r,d;
  if (fabs(b.r) >= fabs(b.i))
    {
      r = b.i / b.r;
      d = b.r * (1. + r * r);
      c.r = a / d;
      c.i = -r * a / d;
    }
  else
    {
      r = b.r / b.i;
      d = b.i * (1. + r * r);
      c.r = a * r / d;
      c.i = -a / d;
    }
  return c;
}


/**
 *  a/b
 * @param a  a complex number
 * @param b  a real number
 * @return a/b, the division
 */
dcomplex CRdiv(dcomplex a, double b)
{
  dcomplex c;
  c = a;
  c.r /= b;
  c.i /= b;
  return c;
}
/**
 * Re(z)^2 + Im(z)^2 
 * @param z  a complex number
 * @return Re(z)^2 + Im(z)^2 
 */
double Csqr_norm(dcomplex z)
{
  return z.r * z.r + z.i * z.i;
}

/**
 * Modulus function
 * @param z a complex number
 * @return sqrt(x^2 + y^2)
 */
double Cabs(dcomplex z)
{
  double x, y, t;
  x = fabs (z.r);
  y = fabs (z.i);
  if (x == 0.) return y;
  if (y == 0.) return x;
  if (x > y)
    {
      t = y / x;
      return x * sqrt (1. + t * t);
    }
  t = x / y;
  return y * sqrt (1. + t * t);
}

/**
 * sqrt(z) , one square root (with positive real part)
 *
 * c = sqrt(z)
 *
 * If z.r >= 0
 *      c.r = sqrt (2 (|z.r| + sqrt(z.r^2 + z.i^2))) / 2
 *      c.i = z.i / sqrt (2 (|z.r| + sqrt(z.r^2 + z.i^2))) 
 * else
 *      c.r = |z.i| / sqrt (2 (|z.r| + sqrt(z.r^2 + z.i^2))) 
 *      c.i = sign(z.i) sqrt (2 (|z.r| + sqrt(z.r^2 + z.i^2))) / 2
 *
 * @param z  a complex number
 * @return  sqrt(z) , one square root (with positive real part)
 */
dcomplex Csqrt(dcomplex z)
{
  dcomplex c;
  double x,y,t;
  if ((z.r == 0.0) && (z.i == 0.0))
    {
      return CZERO;;
    }

  x=fabs(z.r);
  y=fabs(z.i);
  t = sqrt(2. * (x + Cabs (z)));
  if (z.r >= 0.0)
    {
      c.r = 0.5 * t;
      c.i = z.i / t;
    }
  else
    {
      c.r = y / t;
      c.i = 0.5 * ((z.i >= 0) ? t : -t);
    }
  return c;
}

/**
 *  log(z)  
 * @param z  a complex number
 * @return log(z) , the log function
 */
dcomplex Clog(dcomplex z)
{
  dcomplex c;
  if (z.r == 0. && z.i == 0.)
    {
      PNL_ERROR("z==0", "Clog");
    }
  if (z.i == 0.0)
    {
      c.r = log( fabs(z.r) );
      if (z.r > 0.0)
        {
          c.i = 0.0;
        }
      else
        {
          c.i = M_PI;
        }
      return c;
    }
  if (z.r == 0.0)
    {
      c.r = log( fabs(z.i) );
      c.i = M_PI_2;
      if (z.i < 0.0)
        {
          c.i = -c.i;
        }
      return c;
    }
  
  /* z.r and z.i != 0. */
  c.i = Carg (z);
  c.r = 0.5 * log ( z.r * z.r + z.i * z.i);
  return c;
}

/**
 *  the exponential function exp(z)
 * @param z  a complex number
 * @return  exp(z)
 */
dcomplex Cexp(dcomplex z)
{
  return Complex_polar (exp(z.r), z.i);
}

/**
 * the unitary exponential function
 * @param t   a real number
 * @return exp(i t), the unitary exponential
 */
dcomplex CIexp(double t)
{
  return Complex ( cos ( t ), sin ( t ) );
}

/**
 *  z^y  the power function
 * @param z  a complex number
 * @param y  a complex exponential factor
 * @return  z^y the power function 
 */
dcomplex Cpow(dcomplex z, dcomplex y)
{
  double logr, theta;
  if (y.r == 0. && y.i == 0.)
    {
      return CONE;
    }
  /* now y is not nul */
  if (z.r == 0. && z.i == 0.)
    {
      return CZERO;
    }
  if (y.r == 1. && y.i == 0.)
    {
      return z;
    }
  else if (y.r == -1. && y.i == 0.)
    {
      return Cinv(z);
    }

  theta = Carg(z);
  logr = log (Cabs(z));
  return Complex_polar (exp (y.r * logr - y.i * theta),
                        y.i * logr + y.r * theta);
}

/**
 *  z^y  the power function
 * @param z  a complex number
 * @param y  a real exponential factor
 * @return  z^y the power function
 */
dcomplex Cpow_real (dcomplex z, double y)
{
  double logr, theta;
  if (y == 0. )
    {
      return CONE;
    }
  /* now y is not nul */
  if (z.r == 0. && z.i == 0.)
    {
      return CZERO;
    }
  if (y == 1. )
    {
      return z;
    }
  else if (y == -1. )
    {
      return Cinv(z);
    }
  theta = Carg(z);
  logr = log (Cabs(z));
  return Complex_polar (exp(y * logr), y * theta);
}
 
/**
 *  cos(z)
 * @param z  a complex number
 * @return  cos(z)
 */
dcomplex Ccos (dcomplex z)
{
  dcomplex c;
  c.r = cos (z.r) * cosh (z.i);
  c.i = - sin (z.r) * sinh (z.i);

  return c;
}

/**
 *  sin(z)
 * @param z  a complex number
 * @return  sin(z)
 */
dcomplex Csin (dcomplex z)
{
  dcomplex c;
  c.r = sin (z.r) * cosh (z.i);
  c.i = cos (z.r) * sinh (z.i);

  return c;
}

/**
 * tan(z)
 * @param z  a complex number
 * @return  tan(z)
 */
dcomplex Ctan (dcomplex z)
{
  return Cdiv (Csin (z), Ccos (z));
}

/**
 * cotan(z)
 * @param z  a complex number
 * @return  cotan(z)
 */
dcomplex Ccotan (dcomplex z)
{
  return Cdiv (Ccos (z), Csin (z));
}

/**
 *  cosh(z)
 * @param z  a complex number
 * @return  cosh(z)
 */
dcomplex Ccosh(dcomplex z)
{
  dcomplex c;
    
  c.r = cos(z.i) * cosh (z.r);
  c.i = sin(z.i) * sinh (z.r);
  return c;
}


/**
 * sinh(z)
 * @param z  a complex number
 * @return sinh(z)
 */
dcomplex Csinh(dcomplex z)
{
  dcomplex c;
    
  c.r = cos(z.i) * sinh (z.r);
  c.i = sin(z.i) * cosh (z.r);
  return c;
}

/**
 * tanh
 * based on the formula tanh(z) = -i * tan(i * z)
 * @param z  a complex number
 * @return tanh(z)
 */
dcomplex Ctanh(dcomplex z)
{
  dcomplex tmp, tmp2;
  tmp.r = - z.i;
  tmp.i = z.r;
  tmp2 = Ctan (tmp);

  tmp.r = tmp2.i;
  tmp.i = -tmp2.r;
  return tmp;
}

/**
 * cotanh
 * based on the formula cotanh(z) = i * cotan(i * z)
 * @param z  a complex number
 * @return cotanh(z)
 */
dcomplex Ccotanh(dcomplex z)
{
  dcomplex tmp, tmp2;
  tmp.r = - z.i;
  tmp.i = z.r;
  tmp2 = Ccotan (tmp);

  tmp.r = -tmp2.i;
  tmp.i = tmp2.r;
  return tmp;
}

/**
 *  arg(z) , argument
 * @param z  a complex number
 * @return  arg(z) 
 */
double Carg(dcomplex z)       /* arg(z) in [-pi,pi] */
{
  return atan2 (z.i, z.r);
}

/*
 * All the functions C_op_xxx have written and (C) by David Pommier
 */

/**
 *  a+ i b 
 * @param a  a complex number
 * @param b  a complex number 
 * @return the sum:  a+i b 
 */
dcomplex C_op_apib(dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r=a.r-b.i;
  c.i=a.i+b.r;
  return c;
}

/**
 *  a- i b 
 * @param a  a complex number
 * @param b  a complex number 
 * @return  a- i b 
 */
dcomplex C_op_amib(dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r=a.r+b.i;
  c.i=a.i-b.r;
  return c;
}

/**
 *  a+conj(b) 
 * @param a  a complex number
 * @param b  a complex number 
 * @return  a+conj(b)
 */
dcomplex C_op_apcb(dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r=a.r+b.r;
  c.i=a.i-b.i;
  return c;
}

/**
 *  a-conj(b) 
 * @param a  a complex number
 * @param b  a complex number 
 * @return  a-conj(b)
 */
dcomplex C_op_amcb(dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r=a.r-b.r;
  c.i=a.i+b.i;
  return c;
}

/**
 *  d(a+b) 
 * @param d a real number
 * @param a  a complex number
 * @param b  a complex number 
 * @return  d(a+b)
 */
dcomplex C_op_dapb(double d,dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r=d*(a.r+b.r);
  c.i=d*(a.i+b.i);
  return c;
}

/**
 *  d(a-b) 
 * @param d a real number
 * @param a  a complex number
 * @param b  a complex number 
 * @return  d(a-b)
 */
dcomplex C_op_damb(double d,dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r=d*(a.r-b.r);
  c.i=d*(a.i-b.i);
  return c;
}

/**
 *  d(a+ib) 
 * @param d a real number
 * @param a  a complex number
 * @param b  a complex number 
 * @return  d (a+ i b)
 */
dcomplex C_op_dapib(double d,dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r=d*(a.r-b.i);
  c.i=d*(a.i+b.r);
  return c;
}

/**
 *  d(a-ib) 
 * @param d a real number
 * @param a  a complex number
 * @param b  a complex number 
 * @return  d (a- i b)
 */
dcomplex C_op_damib(double d,dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r=d*(a.r+b.i);
  c.i=d*(a.i-b.r);
  return c;
}

/**
 *  d(a+conj(b)) 
 * @param d a real number
 * @param a  a complex number
 * @param b  a complex number 
 * @return  d(a+conj(b))
 */
dcomplex C_op_dapcb(double d,dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r=d*(a.r+b.r);
  c.i=d*(a.i-b.i);
  return c;
}

/**
 *  d(a-conj(b)) 
 * @param d a real number
 * @param a  a complex number
 * @param b  a complex number 
 * @return  d(a-conj(b))
 */
dcomplex C_op_damcb(double d,dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r=d*(a.r-b.r);
  c.i=d*(a.i+b.i);
  return c;
}
/* */

/**
 * i d(a+b) 
 * @param d a real number
 * @param a  a complex number
 * @param b  a complex number 
 * @return  i d(a+b)
 */
dcomplex C_op_idapb(double d,dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r=-d*(a.i+b.i);
  c.i=d*(a.r+b.r);
  return c;
}

/**
 * i d(a-b) 
 * @param d a real number
 * @param a  a complex number
 * @param b  a complex number 
 * @return  i d(a-b)
 */
dcomplex C_op_idamb(double d,dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r=d*(-a.i+b.i);
  c.i=d*(a.r-b.r);
  return c;
}

/**
 * i d(a+conj(b)) 
 * @param d a real number
 * @param a  a complex number
 * @param b  a complex number 
 * @return  i d(a+conj(b))
 */
dcomplex C_op_idapcb(double d,dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r=d*(-a.i+b.i);
  c.i=d*(a.r+b.r);
  return c;
}

/**
 * i d(a-conj(b)) 
 * @param d a real number
 * @param a  a complex number
 * @param b  a complex number 
 * @return i d(a-conj(b))
 */
dcomplex C_op_idamcb(double d,dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r=-d*(a.i+b.i);
  c.i=d*(a.r-b.r);
  return c;
}

