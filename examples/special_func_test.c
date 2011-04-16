
/************************************************************************/
/* Copyright Jérôme Lelong <jerome.lelong@gmail.com>                    */
/* Copyright David Pommier <david.pommier@gmail.com>                    */
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_specfun.h"
#include "pnl/pnl_random.h"
#include "tests_utils.h"

static void exp_int_test ()
{
#ifdef HAVE_GSL
  int      n;
  double Gamma_gsl,IP_gsl,IQ_gsl,En_gsl,Gamma_neg_gsl;
  double a,x;
  double Gamma,IP,IQ,En,Gamma_neg;

  a=0.2;
  x=5.0;
  
  Gamma = pnl_sf_gamma_inc(a,0);
  IP = pnl_sf_gamma_inc_P(a,x);
  IQ = pnl_sf_gamma_inc_Q(a,x);
  /* Ei=pnl_sf_expint_Ei(x); */
  Gamma_neg=pnl_sf_gamma_inc(-a,x);

  printf( "\nTest of Exponential Integrals\n");

  Gamma_gsl=gsl_sf_gamma_inc(a,0);
  Gamma_neg_gsl=gsl_sf_gamma_inc(-a,x);
  IQ_gsl=gsl_sf_gamma_inc_Q (a,x);
  IP_gsl=gsl_sf_gamma_inc_P (a,x);
  /* Ei_gsl=gsl_sf_expint_Ei (x); */

  for(n=1;n<10;n++)
    {
      En=pnl_sf_expint_En(n,x);
      En_gsl=gsl_sf_expint_En (n,x);
      printf( "  E_%d      = %f - %f = %f \n",n,En,En_gsl,En-En_gsl);
    }
  printf( "Gamma     = %f - %f = %f \n",Gamma,Gamma_gsl,Gamma-Gamma_gsl);
  printf( "Gamma_neg = %f - %f = %f \n",Gamma_neg,Gamma_neg_gsl,Gamma_neg-Gamma_neg_gsl);
  printf( "  IQ      = %f - %f = %f \n",IQ,IQ_gsl,IQ-IQ_gsl);
  printf( "  IP      = %f - %f = %f \n",IP,IP_gsl,IP-IP_gsl);
  /* printf( "  Ei      = %f - %f = %f \n",Ei,Ei_gsl,Ei-Ei_gsl); */
  printf( "  En      = %f - %f = %f \n",En,En_gsl,En-En_gsl);
#else
  printf("Tests for Exponential Integrals only available with GSL\n");
#endif
}


struct d2d_test
{
  char *label;
  double (*f)(double);
  double arg;
  double res;
};

struct dd2d_test
{
  char *label;
  double (*f)(double, double);
  double arg;
  double nu;
  double res;
};

struct dd2c_test
{
  char *label;
  dcomplex (*f)(double, double);
  double arg;
  double nu;
  double res_r, res_i;
};

struct dc2c_test
{
  char *label;
  dcomplex (*f)(double, dcomplex);
  double arg_r, arg_i;
  double nu;
  double res_r, res_i;
};


struct d2d_test list_gamma_tst [] =
{
#include "gamma_test.dat"
    { NULL, NULL, 0, 0}
};

struct dd2d_test list_real_bessel_tst [] =
{
#include "real_bessel_test.dat"
    { NULL, NULL, 0, 0}
};

struct dd2c_test list_real_besselh_tst [] =
{
#include "real_besselh_test.dat"
    { NULL, NULL, 0, 0}
};

struct dc2c_test list_complex_bessel_tst [] =
{
#include "complex_bessel_test.dat"
    { NULL, NULL, 0, 0, 0, 0, 0}
};


void d2d_funcs_test (struct d2d_test *tst)
{
  int i;
  double tol  = 1E-5;
  for ( i=0 ; tst[i].f != NULL ; i++ )
    {
      struct d2d_test t = tst[i];
      double res = (t.f)(t.arg);
      pnl_test_eq (res, t.res, tol, t.label, 
                  "computed at %g", t.arg);
    }
}

void dd2d_funcs_test (struct dd2d_test *tst)
{
  int i;
  double tol  = 1E-5;
  for ( i=0 ; tst[i].f != NULL ; i++ )
    {
      struct dd2d_test t = tst[i];
      double res = (t.f)(t.nu, t.arg);
      pnl_test_eq (res, t.res, tol, t.label, 
                  "computed at (%g, %g)", t.nu, t.arg);
    }
}

void dd2c_funcs_test (struct dd2c_test *tst)
{
  int i;
  double tol  = 1E-5;
  for ( i=0 ; tst[i].f != NULL ; i++ )
    {
      struct dd2c_test t = tst[i];
      dcomplex res = (t.f)(t.nu, t.arg);
      pnl_test_eq (res.r, t.res_r, tol, t.label, 
                  "computed at (%g,%g+i %g)", t.nu, t.arg);
      pnl_test_eq (res.i, t.res_i, tol, t.label, 
                  "computed at (%g,%g+i %g)", t.nu, t.arg);
    }
}


void dc2c_funcs_test (struct dc2c_test *tst)
{
  int i;
  double tol  = 1E-4;
  for ( i=0 ; tst[i].f != NULL ; i++ )
    {
      struct dc2c_test t = tst[i];
      dcomplex arg = Complex(t.arg_r, t.arg_i);
      dcomplex res = (t.f)(t.nu, arg);
      pnl_test_eq (res.r, t.res_r, tol, t.label, 
                  "computed at (%g,%g+i %g)", t.nu, t.arg_r, t.arg_i);
      pnl_test_eq (res.i, t.res_i, tol, t.label, 
                  "computed at (%g,%g+i %g)", t.nu, t.arg_r, t.arg_i);
    }
}



void hyperg_test ()
{
#ifdef HAVE_GSL
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  double x, a, b, c, r, r_gsl;
  pnl_rand_init (1, 1, gen);
  x = fabs (pnl_rand_normal (gen));
  a = fabs (pnl_rand_normal (gen));
  b = fabs (pnl_rand_normal (gen));
  c = fabs (pnl_rand_normal (gen));
  printf( "\nTest error of Hypergeometric functions\n");
  printf ("a = %f; b = %f; c = %f; x = %f\n", a,b,c,x);
  r_gsl = gsl_sf_hyperg_0F1(c,x);
  r = pnl_sf_hyperg_0F1(c,x);
  printf ("error on 0F1 %f - %f = %f\n", r, r_gsl, r-r_gsl);
  r_gsl = gsl_sf_hyperg_1F1(a,b,x);
  r = pnl_sf_hyperg_1F1(a,b,x);
  printf ("error on 1F1 %f - %f = %f\n", r, r_gsl, r-r_gsl);
  r = pnl_sf_hyperg_2F0(a,b,-x);
  r_gsl = gsl_sf_hyperg_2F0(a,b,-x);
  printf ("error on 2F0 %f - %f = %f\n", r, r_gsl, r-r_gsl);

  x = pnl_rand_uni (gen);
  r_gsl = gsl_sf_hyperg_2F1(a,b,c,x);
  r = pnl_sf_hyperg_2F1(a,b,c,x);
  printf ("error on 2F1 %f - %f = %f\n", r, r_gsl, r-r_gsl);
  r_gsl = gsl_sf_hyperg_U(a,b,x);
  r = pnl_sf_hyperg_U(a,b,x);
  printf ("error on U %f - %f = %f\n", r, r_gsl, r-r_gsl);
#else
  printf( "Test error of Hypergoemtric functions only available with GSL\n");  
#endif
}

int main (int argc, char **argv)
{
  pnl_test_init (argc, argv);
  d2d_funcs_test (list_gamma_tst);
  dd2d_funcs_test (list_real_bessel_tst);
  dd2c_funcs_test (list_real_besselh_tst);
  dc2c_funcs_test (list_complex_bessel_tst);

  exp_int_test();
  hyperg_test ();
  exit (pnl_test_finalize("Special functions"));
}
