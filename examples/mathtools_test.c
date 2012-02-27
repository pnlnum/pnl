#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_mathtools.h"
#include "tests_utils.h"

struct pow_i_tests
{
  char *label;
  double (*f) (double, int);
  double arg_x;
  int arg_y;
  double res;
};

struct pow_i_tests pow_i_tst [] =
{
#include "mathtools_test.dat"
    { NULL, NULL, 0, 0, 0}
};

static void test_pow_i ()
{
  int i;
  double tol  = 1E6 * DBL_EPSILON;
  for ( i=0 ; pow_i_tst[i].f != NULL ; i++ )
    {
      struct pow_i_tests t = pow_i_tst[i];
      /* double res = (t.f)(t.arg_x, t.arg_y); */
      double res = pow(t.arg_x, t.arg_y);
      pnl_test_eq (res, t.res, tol, t.label, 
                  "computed at (%g,%d)", t.arg_x, t.arg_y);
    }
}

struct invh_test { double arg; double res; };
#define N_acosh 4
#define N_asinh 11
#define N_atanh 9

static void test_invh ()
{ 
  int i;
  double tol = 1E-15;
  struct invh_test acosh_dat[N_acosh] =
    {
        { 1., 0. },
        { 1.05, 0.314924756603848 },
        { 10., 2.99322284612638 },
        { 1E10, 23.7189981105004 }
    };

  struct invh_test asinh_dat[N_asinh] = 
    {
        { 0., 0. },
        { 1E-10, 1E-10 },
        { -1E-10, -1E-10 },
        { 0.1, 0.0998340788992076 },
        { -0.1, -0.0998340788992076 },
        { 1., 0.881373587019543 },
        { -1., -0.881373587019543 },
        { 10, 2.99822295029797 },
        { -10, -2.99822295029797 },
        { 1E10, 23.7189981105004 },
        { -1E10, -23.7189981105004 }
    };

  struct invh_test atanh_dat[N_atanh] =
    {
        { 0., 0. },
        { 1E-10, 1E-10 },
        { -1E-10, -1E-10 },
        { 0.1, 0.100335347731076 },
        { -0.1, -0.100335347731076 },
        { 0.9, 1.47221948958322 },
        { -0.9, -1.47221948958322 },
        { 1, PNL_POSINF },
        { -1, PNL_NEGINF }
    };

  for ( i=0 ; i<N_acosh ; i++ )
    {
      struct invh_test t = acosh_dat[i];
      double res = pnl_acosh (t.arg);
      pnl_test_eq (res, t.res, tol, "pnl_acosh", "computed at %g", t.arg);
    }
  for ( i=0 ; i<N_asinh ; i++ )
    {
      struct invh_test t = asinh_dat[i];
      double res = pnl_asinh (t.arg);
      pnl_test_eq (res, t.res, tol, "pnl_asinh", "computed at %g", t.arg);
    }
  for ( i=0 ; i<N_atanh ; i++ )
    {
      struct invh_test t = atanh_dat[i];
      double res = pnl_atanh (t.arg);
      pnl_test_eq (res, t.res, tol, "pnl_atanh", "computed at %g", t.arg);
    }
}

int main (int argc, char *argv[])
{
  pnl_test_init (argc, argv);
  test_pow_i ();
  test_invh ();
  exit(pnl_test_finalize("Mathtools Functions"));
}

