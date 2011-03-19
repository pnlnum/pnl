#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_mathtools.h"
#include "tests_utils.h"

struct mathtools_tests
{
  char *label;
  double (*f) (double, int);
  double arg_x;
  int arg_y;
  double res;
};

struct mathtools_tests list_tst [] =
{
#include "mathtools_test.dat"
    { NULL, NULL, 0, 0, 0}
};


int main (int argc, char *argv[])
{
  int i;
  double tol  = 1E6 * DBL_EPSILON;
  pnl_test_init ();
  if ( argc == 2 && strcmp (argv[1], "-v") == 0 )
    {
      verbose = 1;
    }
  for ( i=0 ; list_tst[i].f != NULL ; i++ )
    {
      struct mathtools_tests t = list_tst[i];
      /* double res = (t.f)(t.arg_x, t.arg_y); */
      double res = pow(t.arg_x, t.arg_y);
      pnl_test_eq (res, t.res, tol, t.label, 
                  "computed at (%g,%d)", t.arg_x, t.arg_y);
    }
  exit(pnl_test_finalize("Mathtools Functions"));
}

