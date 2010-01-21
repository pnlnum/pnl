
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
#include "pnl_complex.h"
#include "pnl_random.h"

static void Csqrt_test ()
{
  dcomplex z, c;
  z = Complex (3., 5.);
  c = Csqrt (z);
  printf("Csqrt(%f + %f i) = %f + %f i\n", z.r, z.i, c.r, c.i);
  z = Complex (5., 3.);
  c = Csqrt (z);
  printf("Csqrt(%f + %f i) = %f + %f i\n", z.r, z.i, c.r, c.i);
}

static void Clog_test ()
{
  dcomplex z, c;
  z = Complex (3., 5.);
  c = Clog (z);
  printf("Clog(%f + %f i) = %f + %f i\n", z.r, z.i, c.r, c.i);
  z = Complex (-5., 3.);
  c = Clog (z);
  printf("Clog(%f + %f i) = %f + %f i\n", z.r, z.i, c.r, c.i);
}

static void Cpow_test ()
{
  dcomplex a, b, c, d;
  a = Complex (3., 5.);
  b = Complex (2., -1.5);
  c = Cpow (a, b);
  d = Cpow_real (a, Creal(b));
  printf("Cpow(%f + %f i, %f + %f i) = %f + %f i\n", CMPLX(a), CMPLX(b), CMPLX(c));
  printf("Cpow(%f + %f i, %f ) = %f + %f i\n", CMPLX(a), Creal(b), CMPLX(d));
  a = Complex (-3., 5.);
  b = Complex (2., -1.5);
  c = Cpow (a, b);
  c = Cpow (a, b);
  printf("Cpow(%f + %f i, %f + %f i) = %f + %f i\n", CMPLX(a), CMPLX(b), CMPLX(c));
  printf("Cpow(%f + %f i, %f ) = %f + %f i\n", CMPLX(a), Creal(b), CMPLX(d));
}

static void Ctrigo_test ()
{
  dcomplex z, c;
  z = Complex (3., 5.);
  c = Ccos (z);
  printf("Ccos(%f + %f i) = %f + %f i\n", z.r, z.i, c.r, c.i);
  c = Csin (z);
  printf("Csin(%f + %f i) = %f + %f i\n", z.r, z.i, c.r, c.i);
  c = Ctan (z);
  printf("Ctan(%f + %f i) = %f + %f i\n", z.r, z.i, c.r, c.i);
  c = Ccotan (z);
  printf("Ccotan(%f + %f i) = %f + %f i\n", z.r, z.i, c.r, c.i);
  c = Ccosh (z);
  printf("Ccosh(%f + %f i) = %f + %f i\n", z.r, z.i, c.r, c.i);
  c = Csinh (z);
  printf("Csinh(%f + %f i) = %f + %f i\n", z.r, z.i, c.r, c.i);
  c = Ctanh (z);
  printf("Ctanh(%f + %f i) = %f + %f i\n", z.r, z.i, c.r, c.i);
  c = Ccotanh (z);
  printf("Ccotanh(%f + %f i) = %f + %f i\n", z.r, z.i, c.r, c.i);
}

static void Cgamma_test ()
{
  dcomplex z, c1, c2;
  z = Complex (3., 2.);
  c1 = Ctgamma (z);
  c2 = Clgamma (z);
  printf("Cgamma(%f + %f i) = %f + %f i\n", CMPLX(z), CMPLX(c1));
  printf("Clgamma(%f + %f i) = %f + %f i\n", CMPLX(z), CMPLX(c2));
  z = Complex (-3., 2.);
  c1 = Ctgamma (z);
  c2 = Clgamma (z);
  printf("Cgamma(%f + %f i) = %f + %f i\n", CMPLX(z), CMPLX(c1));
  printf("Clgamma(%f + %f i) = %f + %f i\n", CMPLX(z), CMPLX(c2));
} 

static void Cdiv_test ()
{
  dcomplex z1, z2;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init (gen, 1, 1);
  z1 = Complex (pnl_rand_uni (gen), pnl_rand_uni (gen));
  z2 = Complex (pnl_rand_uni (gen), pnl_rand_uni (gen));

  printf ("(%.12f + %.12f * %%i) / (%.12f + %.12f * %%i) = (%.12f + %.12f * %%i)\n", CMPLX(z1), CMPLX(z2), CMPLX(Cdiv(z1,z2)));
  printf ("%.12f / (%.12f + %.12f * %%i) = (%.12f + %.12f * %%i)\n", z1.r, CMPLX(z2), CMPLX(RCdiv(z1.r,z2)));
  printf ("(%.12f + %.12f * %%i) / %.12f = (%.12f + %.12f * %%i)\n", CMPLX(z1), z2.r, CMPLX(CRdiv(z1,z2.r)));
} 


void complex_test()
{
  printf("\n");
  printf("TEST DES FONCTIONS PNL_COMPLEX \n");
  printf("\n");
  Csqrt_test ();
  Clog_test ();
  Cpow_test ();
  Ctrigo_test ();
  Cgamma_test ();
  Cdiv_test ();
}
