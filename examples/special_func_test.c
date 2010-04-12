
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

#include "pnl_specfun.h"
#include "pnl_random.h"

#ifdef HAVE_GSL
extern double gsl_sf_expint_Ei (double x);
extern double gsl_sf_expint_En (int n,double x);
extern double gsl_sf_gamma_inc(double a,double x);
extern double gsl_sf_gamma_inc_Q (double a,double x);
extern double gsl_sf_gamma_inc_P (double a,double x);
extern double gsl_sf_gamma (double);
extern double gsl_sf_lngamma (double);
extern double gsl_sf_log_erfc (double x);
extern double gsl_sf_erf (double x);
extern double gsl_sf_erfc (double x);
extern double gsl_sf_hyperg_2F1 (double a, double b, double c, double x);
extern double gsl_sf_hyperg_1F1 (double a, double b, double x);
extern double gsl_sf_hyperg_2F0 (double a, double b, double x);
extern double gsl_sf_hyperg_U (double a, double b, double x);
extern double gsl_sf_hyperg_0F1 (double c, double x);
#endif


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


void gamma_test ()
{
#ifdef HAVE_GSL
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  double x;
  pnl_rand_init (1, 1, gen);
  x = fabs (pnl_rand_normal (gen));
  printf( "\nTest error of Gamma functions\n");
  printf ("error on gamma     : %f\n", gsl_sf_gamma (x) - pnl_sf_gamma (x));
  printf ("error on log_gamma : %f\n", gsl_sf_lngamma (x) - pnl_sf_log_gamma (x));
  printf ("error on erf       : %f\n", gsl_sf_erf (x) - pnl_sf_erf (x));
  printf ("error on erfc      : %f\n", gsl_sf_erfc (x) - pnl_sf_erfc (x));
  printf ("error on log_erfc  : %f\n", gsl_sf_log_erfc (x) - pnl_sf_log_erfc (x));
#else
  printf( "Test error of Gamma functions only available with GSL\n");  
#endif
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

static void complex_bessel_test ()
{
  dcomplex z, c;
  double v = -1.5;
  z = Complex (5., 3.);
  printf("Test des fonctions de bessel complexes \n");
  c = pnl_complex_bessel_i (v, z);
  printf("bessel_i(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_i_scaled (v, z);
  printf("bessel_i_scaled(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_j (v, z);
  printf("bessel_j(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_j_scaled (v, z);
  printf("bessel_j_scaled(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_y (v, z);
  printf("bessel_y(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_y_scaled (v, z);
  printf("bessel_y_scaled(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_k (v, z);
  printf("bessel_k(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_k_scaled (v, z);
  printf("bessel_k_scaled(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_h1 (v, z);
  printf("bessel_h1(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_h1_scaled (v, z);
  printf("bessel_h1_scaled(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_h2 (v, z);
  printf("bessel_h2(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_h2_scaled (v, z);
  printf("bessel_h2_scaled(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);

}

static void real_bessel_test ()
{
  double z, c;
  dcomplex zc;
  double v = -1.5;
  z = 5.;
  printf("Test des fonctions de bessel reelles \n");
  c = pnl_bessel_i (v, z);
  printf("bessel_i(%f, %f) = %f\n", v, z, c);
  c = pnl_bessel_i_scaled (v, z);
  printf("bessel_i_scaled(%f, %f) = %f\n", v, z, c);
  c = pnl_bessel_j (v, z);
  printf("bessel_j(%f, %f) = %f\n", v, z, c);
  c = pnl_bessel_j_scaled (v, z);
  printf("bessel_j_scaled(%f, %f) = %f\n", v, z, c);
  c = pnl_bessel_y (v, z);
  printf("bessel_y(%f, %f) = %f\n", v, z, c);
  c = pnl_bessel_y_scaled (v, z);
  printf("bessel_y_scaled(%f, %f) = %f\n", v, z, c);
  c = pnl_bessel_k (v, z);
  printf("bessel_k(%f, %f) = %f\n", v, z, c);
  c = pnl_bessel_k_scaled (v, z);
  printf("bessel_k_scaled(%f, %f) = %f\n", v, z, c);
  zc = pnl_bessel_h1 (v, z);
  printf("bessel_h1(%f, %f) = %f + %f i\n", v, z, CMPLX(zc));
  zc = pnl_bessel_h1_scaled (v, z);
  printf("bessel_h1_scaled(%f, %f) = %f + %f i\n", v, z, CMPLX(zc));
  zc = pnl_bessel_h2 (v, z);
  printf("bessel_h2(%f, %f) = %f + %f i\n", v, z, CMPLX(zc));
  zc = pnl_bessel_h2_scaled (v, z);
  printf("bessel_h2_scaled(%f, %f) = %f + %f i\n", v, z, CMPLX(zc));
}

void special_func_test ()
{
  printf("\n");
  printf("Special function tests.\n");
  printf("\n");
  complex_bessel_test ();
  real_bessel_test ();
  exp_int_test();
  gamma_test ();
  hyperg_test ();
}
