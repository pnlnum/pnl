
/*
 * Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as 
 * published by  the Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License  along with this program.  If not, see
 * <http://www.gnu.org/licenses/>. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "pnl/pnl_mathtools.h"
#include "tests_utils.h"

int verbose = 0;
int count_tests;
int count_ok;
int count_fail;

/** 
 * Initializes test counting system
 */
void pnl_test_init (int argc, char **argv)
{
  if ( argc == 2 && strcmp (argv[1], "-v") == 0 )
    {
      verbose = 1;
    }
  count_tests = 0;
  count_ok = 0;
  count_fail = 0;
}


/** 
 * Update test counter according to status
 * 
 * @param status an integer. If 0 last test is considered as past and
 * FAILEd otherwise
 */
static void update_count_tests (int status)
{
  count_tests ++;
  if ( status )
    {
      count_fail ++;
    }
  else
    {
      count_ok ++;
    }
}

/** 
 * Prints a summary of the tests
 * 
 * @return 0 if all tests were succeddful and 1 otherwise
 */
int pnl_test_finalize(const char *str)
{
  printf ("%s : %s (TOTAL : %d, PASSED : %d, FAILED %d)\n", str, (count_fail == 0) ? "OK" : "FAIL", count_tests, count_ok, count_fail);
  return ( count_fail >0 );
}


static int test_eq_rel (double x, double y, double relerr)
{
  int status;
  if ( y != 0 )
    {
      status = (fabs (x -y) / fabs(y) > relerr);
    }
  else
    {
      status = (fabs(x) > relerr);
    }
  return status;
}

static int test_eq_abs (double x, double y, double abserr)
{
  int status;
  status = (fabs (x -y) > abserr);
  return status;
}

static int test_eq (double x, double y, double abserr)
{
  int status;
  status = (fabs (x -y) / MAX(1, fabs (y)) > abserr);
  return status;
}


/** 
 * Checks if |x - y| / (max(1, |y|)) < relerr
 * 
 * @param x computed result
 * @param y expected result
 * @param relerr relative error (note that when |y| < 1, it is an abolute
 * error)
 * @param str the name of the tested function
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return FALSE or TRUE
 */
int pnl_test_eq(double x, double y, double relerr, const char *str, const char *fmt, ...)
{
  va_list ap;
  va_start (ap, fmt);
  va_end (ap);
  if ( y >= 1 )
    {
      return pnl_test_eq_rel (x, y, relerr, str, fmt, ap);
    }
  else
    {
      return pnl_test_eq_abs (x, y, relerr, str, fmt, ap);
    }
}

static int pnl_test_eq_aux (double x, double y, double relerr, int(*cmp)(double, double, double), const char *str, const char *fmt, ...)
{
  int status = 0;
  if ( (pnl_isnan (x) && !pnl_isnan(y)) ||
       (pnl_isinf (x) && !pnl_isinf(y)) )
    {
      status = 1;
    }
  else
    {
      status = (*cmp)(x, y, relerr);
    }
  if ( status || verbose )
    {
      printf ("\t%s : ", str);
      printf ( status ? "FAIL" : "OK");
      if ( status ) 
        {
          va_list ap;
          va_start (ap, fmt);
          va_end (ap);
          printf (" (");
          vprintf (fmt, ap);
          printf (" expected %.18g observed %.18g)", y, x);
        }
      printf ("\n");
    }
  update_count_tests (status);
  return (status ? FALSE : TRUE);

}

/** 
 * Checks if |x -y| / |y| < relerr
 * 
 * @param x computed result
 * @param y exepcted value
 * @param relerr maximum relative error
 * @param str the fonctionnality tested
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return FALSE or TRUE
 */
int pnl_test_eq_rel (double x, double y, double relerr, const char *str, const char *fmt, ...)
{
  va_list ap;
  va_start (ap, fmt);
  va_end (ap);
  return pnl_test_eq_aux (x, y, relerr, test_eq_rel, str, fmt, ap);
}

/** 
 * Checks if |x - y|  < abserr
 * 
 * @param x computed result
 * @param y exepcted value
 * @param abserr maximum absolute error
 * @param str the fonctionnality tested
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return FALSE or TRUE
 */
int pnl_test_eq_abs (double x, double y, double abserr, const char *str, const char *fmt, ...)
{
  va_list ap;
  va_start (ap, fmt);
  va_end (ap);
  return pnl_test_eq_aux (x, y, abserr, test_eq_abs, str, fmt, ap);
}

static int pnl_test_array (const double *X, const double *Y, int n, double relerr, int(*cmp)(double, double, double), const char *str, const char *fmt, ...)
{
  int i, status;
  for ( i=0 ; i<n ; i++ )
    {
      const double x = X[i];
      const double y = Y[i];
      status = (*cmp)(x, y, relerr);
      if ( status ) break;
    }
  if ( status || verbose )
    {
      printf ("\t%s : ", str);
      printf ( status ? "FAIL" : "OK");
      if ( status ) 
        {
          va_list ap;
          va_start (ap, fmt);
          printf (" (");
          vprintf (fmt, ap);
          va_end (ap);
          printf (" expected %.18g observed %.18g)", Y[i], X[i]);
        }
      printf ("\n");
    }
  update_count_tests (status);
  return (status ? FALSE : TRUE);

}


/** 
 * Checks if |x(i,j) - y(i,j)| / |y(i,j)| < relerr
 * 
 * @param X computed result (matrix)
 * @param Y expected result (matrix)
 * @param relerr relative error (note that when |y| < 1, it is an abolute
 * error)
 * @param str the name of the tested function
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return FALSE or TRUE
 */
int pnl_test_mat_eq_rel (const PnlMat *X, const PnlMat *Y, double relerr, const char *str, const char *fmt, ...)
{
  va_list ap;
  va_start (ap, fmt);
  va_end (ap);
  if ( X->m != Y->m || X->n != Y->n )
    {
      printf ("%s : ", str);
      printf ("FAIL (size mismatch");
      printf (fmt, ap); printf (")\n");
      return FALSE;
    }
  return pnl_test_array (X->array, Y->array, X->mn, relerr, test_eq_rel, str, fmt, ap);
}

/** 
 * Checks if |x(i,j) - y(i,j)| < abserr
 * 
 * @param X computed result (matrix)
 * @param Y expected result (matrix)
 * @param relerr relative error (note that when |y| < 1, it is an abolute
 * error)
 * @param str the name of the tested function
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return FALSE or TRUE
 */
int pnl_test_mat_eq_abs (const PnlMat *X, const PnlMat *Y, double abserr, const char *str, const char *fmt, ...)
{
  va_list ap;
  va_start (ap, fmt);
  va_end (ap);
  if ( X->m != Y->m || X->n != Y->n )
    {
      printf ("%s : ", str);
      printf ("FAIL (size mismatch");
      printf (fmt, ap); printf (")\n");
      return FALSE;
    }
  return pnl_test_array (X->array, Y->array, X->mn, abserr, test_eq_abs, str, fmt, ap);
}

/** 
 * Checks if |x(i,j) - y(i,j)| / (max(1, |y(i,j)|)) < relerr for all (i,j)
 * 
 * @param X computed result (matrix)
 * @param Y expected result (matrix)
 * @param relerr relative error (note that when |y| < 1, it is an abolute
 * error)
 * @param str the name of the tested function
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return FALSE or TRUE
 */
int pnl_test_mat_eq(const PnlMat *X, const PnlMat *Y, double relerr, const char *str, const char *fmt, ...)
{
  va_list ap;
  va_start (ap, fmt);
  va_end (ap);
  if ( X->m != Y->m || X->n != Y->n )
    {
      printf ("%s : ", str);
      printf ("FAIL (size mismatch");
      printf (fmt, ap); printf (")\n");
      return FALSE;
    }
  return pnl_test_array (X->array, Y->array, X->mn, relerr, test_eq, str, fmt, ap);
}

/** 
 * Checks if |x(i,j) - y(i,j)| / |y(i,j)| < relerr
 * 
 * @param X computed result (vector)
 * @param Y expected result (vector)
 * @param relerr relative error (note that when |y| < 1, it is an abolute
 * error)
 * @param str the name of the tested function
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return FALSE or TRUE
 */
int pnl_test_vect_eq_rel (const PnlVect *X, const PnlVect *Y, double relerr, const char *str, const char *fmt, ...)
{
  va_list ap;
  va_start (ap, fmt);
  va_end (ap);
  if ( X->size != Y->size )
    {
      printf ("%s : ", str);
      printf ("FAIL (size mismatch");
      printf (fmt, ap); printf (")\n");
      return FALSE;
    }
  return pnl_test_array (X->array, Y->array, X->size, relerr, test_eq_rel, str, fmt, ap);
}

/** 
 * Checks if |x(i,j) - y(i,j)| < abserr
 * 
 * @param X computed result (vector)
 * @param Y expected result (vector)
 * @param relerr relative error (note that when |y| < 1, it is an abolute
 * error)
 * @param str the name of the tested function
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return FALSE or TRUE
 */
int pnl_test_vect_eq_abs (const PnlVect *X, const PnlVect *Y, double abserr, const char *str, const char *fmt, ...)
{
  va_list ap;
  va_start (ap, fmt);
  va_end (ap);
  if ( X->size != Y->size )
    {
      printf ("%s : ", str);
      printf ("FAIL (size mismatch");
      printf (fmt, ap); printf (")\n");
      return FALSE;
    }
  return pnl_test_array (X->array, Y->array, X->size, abserr, test_eq_abs, str, fmt, ap);
}

/** 
 * Checks if |x(i,j) - y(i,j)| / (max(1, |y(i,j)|)) < relerr for all (i,j)
 * 
 * @param X computed result (vector)
 * @param Y expected result (vector)
 * @param relerr relative error (note that when |y| < 1, it is an abolute
 * error)
 * @param str the name of the tested function
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return FALSE or TRUE
 */
int pnl_test_vect_eq(const PnlVect *X, const PnlVect *Y, double relerr, const char *str, const char *fmt, ...)
{
  va_list ap;
  va_start (ap, fmt);
  va_end (ap);
  if ( X->size != Y->size )
    {
      printf ("%s : ", str);
      printf ("FAIL (size mismatch");
      printf (fmt, ap); printf (")\n");
      return FALSE;
    }
  return pnl_test_array (X->array, Y->array, X->size, relerr, test_eq, str, fmt, ap);
}



void run_all_test (tst_list *l)
{
  int len=0;
  while (l[len].func != NULL)
    {
      if (strcmp (l[len].label, "all_test") != 0) (l[len].func)();
      len ++;
    }
}

void menu_test (tst_list *l)
{
  int len=0, choice;

  while (l[len].func != NULL)
    {
      printf("%2d. %s\n",  len+1, l[len].label);
      len ++;
    }
  len--;
  printf("Which test do you want to run?\n");

  while (1)
    {
      scanf("%d", &choice); 
      if (choice <1 || choice > len+1) printf("illegal choice\n");
      else break;
    }
  len = 0;
  (l[choice-1].func)();
}
