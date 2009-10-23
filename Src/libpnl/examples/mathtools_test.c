#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl_mathtools.h"

static void pow_i_test ()
{
  int n;
  double y, x;
  printf("Test of the integer power function.\n");
  n = 0;
  x = 0.; y = pnl_pow_i (x, n);
  printf ("%f = %f^%i\n", y, x , n);
  n = 1;
  x = 0.; y = pnl_pow_i (x, n);
  printf ("%f = %f^%i\n", y, x , n);
  n = 4;
  x = 2.; y = pnl_pow_i (x, n);
  printf ("%f = %f^%i\n", y, x , n);
  n = 5;
  x = 2.; y = pnl_pow_i (x, n);
  printf ("%f = %f^%i\n", y, x , n);
  n = 4;
  x = -2.; y = pnl_pow_i (x, n);
  printf ("%f = %f^%i\n", y, x , n);
  n = 5;
  x = -2.; y = pnl_pow_i (x, n);
  printf ("%f = %f^%i\n", y, x , n);
  n = -4;
  x = 2.; y = pnl_pow_i (x, n);
  printf ("%f = %f^%i\n", y, x , n);
  n = -5;
  x = 2.; y = pnl_pow_i (x, n);
} 

void mathtools_test ()
{
  pow_i_test ();
} 
