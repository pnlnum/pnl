
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

static double wrap_pnl_sf_choose (double n, double p)
{
  return pnl_sf_choose ((int) n, (int) p);
}

/* 
 * struct for double(*f)(double)
 */
struct d2d_test
{
  char *label;
  double (*f)(double);
  double arg;
  double res;
};

/* 
 * struct for double(*f)(double, double)
 */
struct dd2d_test
{
  char *label;
  double (*f)(double, double);
  double arg;
  double nu;
  double res;
};

/* 
 * struct for double(*f)(int, double)
 */
struct id2d_test
{
  char *label;
  double (*f)(int, double);
  int n;
  double arg;
  double res;
};

/* 
 * struct for dcomplex(*f)(double, double)
 */
struct dd2c_test
{
  char *label;
  dcomplex (*f)(double, double);
  double arg;
  double nu;
  double res_r, res_i;
};

/* 
 * struct for double(*f)(double, double, double)
 */
struct ddd2d_test
{
  char *label;
  double (*f)(double, double, double);
  double a, b, x;
  double res;
};

/* 
 * struct for double(*f)(double, double, double, double)
 */
struct dddd2d_test
{
  char *label;
  double (*f)(double, double, double, double);
  double a, b, c, x;
  double res;
};

/* 
 * struct for dcomplex(*f)(doubla, dcomplex)
 */
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
#include "Data_specfun/gamma_test.dat"
    { NULL, NULL, 0, 0}
};

struct id2d_test list_expint_tst [] =
{
#include "Data_specfun/expint_test.dat"
    { NULL, NULL, 0, 0, 0}
};

struct dd2d_test list_gammainc_tst [] =
{
#include "Data_specfun/gammainc_test.dat"
    { NULL, NULL, 0, 0, 0}
};

struct dd2d_test list_real_bessel_tst [] =
{
#include "Data_specfun/real_bessel_test.dat"
    { NULL, NULL, 0, 0, 0}
};

struct dd2c_test list_real_besselh_tst [] =
{
#include "Data_specfun/real_besselh_test.dat"
    { NULL, NULL, 0, 0, 0, 0}
};

struct dc2c_test list_complex_bessel_tst [] =
{
#include "Data_specfun/complex_bessel_test.dat"
    { NULL, NULL, 0, 0, 0, 0, 0}
};

struct dd2d_test list_hyper0F1_tst [] =
{
#include "Data_specfun/hyperg0F1_test.dat"
    { NULL, NULL, 0, 0, 0}
};

struct ddd2d_test list_hyper1F1_tst [] =
{
#include "Data_specfun/hyperg1F1_test.dat"
    { NULL, NULL, 0, 0, 0, 0}
};

struct ddd2d_test list_hyperU_tst [] =
{
#include "Data_specfun/hypergU_test.dat"
    { NULL, NULL, 0, 0, 0, 0}
};

struct ddd2d_test list_hyper2F0_tst [] =
{
#include "Data_specfun/hyperg2F0_test.dat"
    { NULL, NULL, 0, 0, 0, 0}
};

struct dddd2d_test list_hyper2F1_tst [] =
{
#include "Data_specfun/hyperg2F1_test.dat"
    { NULL, NULL, 0, 0, 0, 0}
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

void id2d_funcs_test (struct id2d_test *tst)
{
  int i;
  double tol  = 1E-5;
  for ( i=0 ; tst[i].f != NULL ; i++ )
    {
      struct id2d_test t = tst[i];
      double res = (t.f)(t.n, t.arg);
      pnl_test_eq (res, t.res, tol, t.label, 
                  "computed at (%d, %g)", t.n, t.arg);
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
                  "computed at (%.18g, %.18g)", t.nu, t.arg);
    }
}

void ddd2d_funcs_test (struct ddd2d_test *tst)
{
  int i;
  double tol  = 1E-5;
  for ( i=0 ; tst[i].f != NULL ; i++ )
    {
      struct ddd2d_test t = tst[i];
      double res = (t.f)(t.a, t.b, t.x);
      pnl_test_eq (res, t.res, tol, t.label, 
                   "computed at (%g, %g, %g)", t.a, t.b, t.x);
    }
}

void dddd2d_funcs_test (struct dddd2d_test *tst)
{
  int i;
  double tol  = 1E-5;
  for ( i=0 ; tst[i].f != NULL ; i++ )
    {
      struct dddd2d_test t = tst[i];
      double res = (t.f)(t.a, t.b, t.c, t.x);
      pnl_test_eq (res, t.res, tol, t.label, 
                   "computed at (%g, %g, %g, %g)", t.a, t.b, t.c, t.x);
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

int main (int argc, char **argv)
{
  pnl_test_init (argc, argv);
  d2d_funcs_test (list_gamma_tst);
  id2d_funcs_test (list_expint_tst);
  dd2d_funcs_test (list_gammainc_tst);

  dd2d_funcs_test (list_real_bessel_tst);
  dd2c_funcs_test (list_real_besselh_tst);
  dc2c_funcs_test (list_complex_bessel_tst);


  dd2d_funcs_test (list_hyper0F1_tst);
  ddd2d_funcs_test (list_hyper1F1_tst);
  ddd2d_funcs_test (list_hyperU_tst);
  ddd2d_funcs_test (list_hyper2F0_tst);
  dddd2d_funcs_test (list_hyper2F1_tst);

  exit (pnl_test_finalize("Special functions"));
}
