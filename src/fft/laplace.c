
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

#include "pnl/pnl_complex.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_specfun.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_fft.h"
#include "pnl/pnl_laplace.h"

#define ALTERNATE(n) (((n)&1) ? -1 : 1) /* (-1)^n */


/**
 * Euler acceleration to recover a CDF using the Laplace transform of its density
 *
 * @param f the Laplace transform with complex values of a real valued
 * density function
 * @param t the point at which the orignal is to be evaluated
 * @param h step size used to discretize the integral
 * @param N series truncation
 * @param M Euler averaging length.
 */
double pnl_ilap_cdf_euler(PnlCmplxFunc *f, double t, double h, int N, int M)
{
  int    i, Cnp;
  double sum, run_sum;

  sum = 0.;

  for(i=1;i<N+1;i++)
    {
      sum += sin(i*h*t) * Creal (PNL_EVAL_FUNC(f, Complex (0., h * i))) / (double)i;
    }
  run_sum = sum; /* partial sum of sn */
  sum = 0.0; /* partial exponential sum */
  Cnp = 1; /* binomial coefficients */

  for(i=0;i<M+1;i++)
    {
      sum += run_sum * (double) Cnp ;
      run_sum += sin ((i + N + 1) * h * t) *
        Creal (PNL_EVAL_FUNC(f, Complex (0., h * (i + N + 1)))) / (double) (i + N + 1);
      Cnp = (Cnp * (M - i)) / (i + 1);
    }

  return(2.0 / M_PI * sum / pow(2.0,M) + h * t / M_PI);
}

/**
 * Euler acceleration to invert a Laplace transform
 *
 * @param f the Laplace transform with complex values of a real valued
 * function
 * @param t the point at which the orignal is to be evaluated
 * @param N series truncation
 * @param M Euler averaging length.
 */
double pnl_ilap_euler(PnlCmplxFunc *f, double t, int N, int M)
{
  int i, Cnp;
  double  sum, a, pit, run_sum;
  double A;

  A = 13.; /* MIN (13., 13. / t); */

  a = A/(2.0*t);
  pit = M_PI/t;
  sum = 0.5 * Creal (PNL_EVAL_FUNC(f, Complex (a, 0.)));

  for(i=1;i<N+1;i++)
    sum += ALTERNATE(i) * Creal (PNL_EVAL_FUNC(f, Complex (a, pit * i)));

  run_sum = sum; /* partial sum of sn */
  sum = 0.0; /* partial exponential sum */
  Cnp = 1; /* binomial coefficients */

  for(i=0;i<M+1;i++)
    {
      sum +=  run_sum * (double) Cnp ;
      run_sum += ALTERNATE(i+N+1) * Creal (PNL_EVAL_FUNC(f, Complex (a, pit * (i + N + 1))));
      Cnp = (Cnp * (M - i)) / (i + 1);
    }
  return exp (A / 2. - M * M_LN2) * sum / t ;
}

/**
 * FFT algorithm to invert a Laplace transform
 * @param res a real vector containing the result of the inversion. We know that
 * the imaginary part of the inversion vanishes.
 * @param f : the Laplace transform to be inverted
 * @param T : the time horizon up to which the function is to be recovered
 * @param eps : precision required on [0, T]
 */
void pnl_ilap_fft(PnlVect *res, PnlCmplxFunc *f, double T, double eps)
{
  PnlVectComplex *fft;
  int             i, N, size;
  double          h, time_step, a;
  double          f_a, omega;
  dcomplex        mul, fac;

  h = M_PI / (2 * T); 
  a = h * log (1 + 1. / eps) / (M_2PI); 

  N = MAX (sqrt (exp (a * T) / eps), h / (M_2PI * eps)) ;
  N = pow (2., ceil (log (N) / M_LN2) );
  time_step = M_2PI / (N * h);
  
  fft = pnl_vect_complex_create (N);
  size = floor (T / time_step) + 1;
  pnl_vect_resize (res, size);
  
  fac = CIexp (-M_2PI / N);
  mul = fac;
  f_a = Creal (PNL_EVAL_FUNC (f, Complex (a, 0)));
  omega = h;
  
  for (i=0 ; i<N ; i++)
    {
      pnl_vect_complex_set (fft, i, Cmul(mul,PNL_EVAL_FUNC (f, Complex (a, - omega))));
      omega += h;
      mul = Cmul(mul,fac);
    }
  pnl_fft_inplace (fft);
  mul = CONE;

  for (i=0 ; i<size ; i++)
    {
      double res_i;
      res_i = Creal (Cmul (pnl_vect_complex_get (fft, i), mul));
      mul = Cmul (mul, fac);
      res_i = (h / M_PI) * exp (a * (i + 1) * time_step) * (res_i + .5 * f_a);
      PNL_LET (res, i) = res_i;
    }
  pnl_vect_complex_free (&fft);
}

/**
 * Gave Stehfest inverse Laplace formula
 * Unlike the Euler method, the inversion is here performed on the real line
 *
 * @param fhat the real valued Laplace transform of f
 * @param t the point at which f is to be recovered
 * @param n the number of iterations of the method
 * @return f(t)
 */
double pnl_ilap_gs_basic (PnlFunc *fhat, double t, int n)
{
  int k, Cnk;
  double alpha, f, fact1, fact2;

  alpha = M_LN2 / t;
  f = 0.;
  Cnk = 1;

  for ( k=0 ; k<n+1 ; k++ )
    {
      f += ALTERNATE(k) * Cnk * PNL_EVAL_FUNC(fhat, (n + k) * alpha);
      Cnk = (Cnk * (n - k)) / (k + 1);
    }
  fact1 = n + 1.;
  fact2 = 1.;
  for ( k=1 ; k<n ; k++)
    {
      fact1 *= n + 1 + k;
      fact2 *= k;
    }
  f *= alpha * (double) fact1 / (double) fact2;
  return f;
}


/**
 * Auxiliary function used in pnl_ilap_gs
 * See the recurrence relation on
 * \verbatim \tilde f_m(t, j) \endverbatim in
 * THE FOURIER-SERIES METHOD FOR INVERTING TRANSFORMS
 *           OF PROBABILITY DISTRIBUTIONS
 * by Abate and Whitt
 * @param fhat the real valued Laplace transform of f
 * @param work is a workspace used by the function to avoid recomputing the same
 * quantity twice.
 * @param iwork is a boolean matrix used to know if the corresponding element
 * of the matrix work has already been computed
 * @param t the point at which f is to be recovered
 * @param m an index
 * @param j an index
 * @return \verbatim \tilde f_m(t, j) \endverbatim
 */
static double f_tilde (PnlFunc *fhat, PnlMat *work, PnlMatInt *iwork, int m, double t, int j)
{
  double x;
  if (PNL_MGET(iwork, m-1, j) == 0)
    { /*
       * In this case, the corresponding f_m(,j) has not been computed yet,
       * so we do it and store the result in work(m-1, j)
      */
      if (j==0)
        {
          double alpha = M_LN2 / t;
          x = m * alpha * PNL_EVAL_FUNC(fhat, m * alpha);
          PNL_MSET (work, m-1, j, x);
        }
      else
        {
          x = (1. + (double) m / (double) j) * f_tilde (fhat, work, iwork, m, t, j - 1) -
            ((double) m / (double) j) * f_tilde (fhat, work, iwork, m + 1, t, j - 1);
        }
      PNL_MSET (work, m-1, j, x);
      PNL_MSET (iwork, m-1, j, 1);
      return x;
    }
  else
    {
      return PNL_MGET (work, m-1, j);
    }
}


/**
 * Inverse Laplace transform
 * Optimal linear combination for the Gaver Stehfest's method.
 * Unlike the Euler method, the inversion is here performed on the real line
 *
 * @param fhat the real valued Laplace transform of f
 * @param t the point at which f is to be recovered
 * @param n the number of iterations of the method
 * @return f(t)
 */
double pnl_ilap_gs (PnlFunc *fhat, double t, int n)
{
  int k;
  double f, Cnk;
  PnlMat *work;
  PnlMatInt *iwork;

  f = 0.;
  Cnk = 1. / pnl_sf_fact (n-1);
  work = pnl_mat_create (2 * n, n+1); 
  iwork = pnl_mat_int_create_from_scalar (2 * n, n+1, 0);
  for ( k=1 ; k<n+1 ; k++ )
    {
      f += ALTERNATE(n-k) * pow(k, n) * Cnk * f_tilde (fhat, work, iwork, k, t, k);
      Cnk = (Cnk * (n - k)) / (k + 1);
    }
  pnl_mat_free (&work);
  pnl_mat_int_free (&iwork);
  return f;
}

#undef ALTERNATE
