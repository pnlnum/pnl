
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

int verbose = FALSE;
int count_tests;
int count_ok;
int count_fail;

/** 
 * Initializes test counting system
 */
void pnl_test_init (int argc, char **argv)
{
  if ( argc == 2 &&  
       ( ( strcmp (argv[1], "-v") == 0 ) || 
         ( strcmp (argv[1], "--verbose") == 0 ) )
     )
    {
      verbose = TRUE;
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
  printf ("%s : %s (TOTAL: %d, PASSED: %d, FAILED: %d)\n", str, (count_fail == 0) ? "OK" : "FAIL", count_tests, count_ok, count_fail);
  return ( count_fail >0 );
}

int pnl_test_is_verbose ()
{
  return verbose;
}

/** 
 * Declares a test as passed
 * 
 * @param str a string
 */
void pnl_test_set_ok (const char *str)
{
  update_count_tests (0);
  if ( verbose == TRUE ) 
    {
      printf ("\t%s : OK\n", str);
    }
}

/** 
 * Declares a test as failed
 * 
 * @param str a string
 * @param res computed result
 * @param expected expected result
 */
void pnl_test_set_fail (const char *str, double res, double expected)
{
  update_count_tests (1);
  printf ("\t%s : FAIL (observed %.18f expected %.18f)\n", str, res, expected);
}

/** 
 * Relative comparison of two real numbers
 * 
 * @param x real number
 * @param y real number
 * @param relerr real number 
 * 
 * @return  0 or 1
 */
int pnl_cmp_eq_rel (double x, double y, double relerr)
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

/** 
 * Absolute comparison of two real numbers
 * 
 * @param x real number
 * @param y real number
 * @param relerr real number 
 * 
 * @return  0 or 1
 */
int pnl_cmp_eq_abs (double x, double y, double abserr)
{
  int status;
  status = (fabs (x -y) > abserr);
  return status;
}

/** 
 * Comparison of two real numbers
 * 
 * @param x real number
 * @param y real number
 * @param relerr real number 
 * 
 * @return  0 or 1
 */
int pnl_cmp_eq (double x, double y, double abserr)
{
  int status;
  status = (fabs (x -y) / MAX(1, fabs (y)) > abserr);
  return status;
}

static int pnl_test_eq_aux (double x, double y, double relerr, int(*cmp)(double, double, double), const char *str, const char *fmt, va_list ap)
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
  if ( status || verbose == TRUE )
    {
      printf ("\t%s : ", str);
      printf ( status ? "FAIL" : "OK");
      if ( status ) 
        {
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
  if ( fabs(y) >= 1 )
    {
      return pnl_test_eq_aux (x, y, relerr, pnl_cmp_eq_rel, str, fmt, ap);
    }
  else
    {
      return pnl_test_eq_aux (x, y, relerr, pnl_cmp_eq_abs, str, fmt, ap);
   }
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
  return pnl_test_eq_aux (x, y, relerr, pnl_cmp_eq_rel, str, fmt, ap);
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
  return pnl_test_eq_aux (x, y, abserr, pnl_cmp_eq_abs, str, fmt, ap);
}

/** 
 * Compares two vectors component-wise using the comparison function
 * specified b cmp
 * 
 * @param X a vector (computed result)
 * @param Y a vector (expected result)
 * @param n size of the expected result
 * @param relerr maximum admissible error in the comparison
 * @param cmp the comparison function
 * @param str name of the tested functionnality
 * @param fmt a format string
 * @param ap extra arguments
 * 
 * @return TRUE or FALSE
 */
static int pnl_test_array (const double *X, const double *Y, int n, double relerr, int(*cmp)(double, double, double), const char *str, const char *fmt, va_list ap)
{
  int i, status;
  status = 0;
  for ( i=0 ; i<n ; i++ )
    {
      const double x = X[i];
      const double y = Y[i];
      status = (*cmp)(x, y, relerr);
      if ( status ) break;
    }
  if ( status || verbose == TRUE )
    {
      printf ("\t%s : ", str);
      printf ( status ? "FAIL" : "OK");
      if ( status ) 
        {
          printf (" (");
          vprintf (fmt, ap);
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
 * @param relerr relative error
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
      update_count_tests (1);
      return FALSE;
    }
  return pnl_test_array (X->array, Y->array, X->mn, relerr, pnl_cmp_eq_rel, str, fmt, ap);
}

/** 
 * Checks if |x(i,j) - y(i,j)| < abserr
 * 
 * @param X computed result (matrix)
 * @param Y expected result (matrix)
 * @param abserr absolute error 
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
      update_count_tests (1);
      return FALSE;
    }
  return pnl_test_array (X->array, Y->array, X->mn, abserr, pnl_cmp_eq_abs, str, fmt, ap);
}

/** 
 * Checks if |x(i,j) - y(i,j)| < abserr
 * 
 * @param X computed result (matrix)
 * @param Y expected result (matrix)
 * @param abserr absolute error 
 * @param str the name of the tested function
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return FALSE or TRUE
 */
int pnl_test_mat_complex_eq_abs (const PnlMatComplex *X, const PnlMatComplex *Y, double abserr, const char *str, const char *fmt, ...)
{
  va_list ap;
  va_start (ap, fmt);
  va_end (ap);
  if ( X->m != Y->m || X->n != Y->n )
    {
      printf ("%s : ", str);
      printf ("FAIL (size mismatch");
      printf (fmt, ap); printf (")\n");
      update_count_tests (1);
      return FALSE;
    }
  return pnl_test_array ((double *)X->array, (double *)Y->array, 2 * X->mn, abserr, pnl_cmp_eq_abs, str, fmt, ap);
}

/** 
 * Checks if x(i,j) = y(i,j) for integer matrices
 * 
 * @param X computed result (matrix)
 * @param Y expected result (matrix)
 * @param str the name of the tested function
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return FALSE or TRUE
 */
int pnl_test_mat_int_eq(const PnlMatInt *X, const PnlMatInt *Y, const char *str, const char *fmt, ...)
{
  int i, status;
  va_list ap;
  va_start (ap, fmt);
  va_end (ap);
  status = 0;
  if ( X->m != Y->m || X->n != Y->n )
    {
      printf ("%s : ", str);
      printf ("FAIL (size mismatch");
      printf (fmt, ap); printf (")\n");
      update_count_tests (1);
      return FALSE;
    }
  for ( i=0 ; i<X->mn ; i++ )
    {
      const int x = X->array[i];
      const int y = Y->array[i];
      status  = (x != y);
      if ( status ) break;
    }
  if ( status || verbose == TRUE )
    {
      printf ("\t%s : ", str);
      printf ( status ? "FAIL" : "OK");
      if ( status ) 
        {
          printf (" (");
          vprintf (fmt, ap);
          va_end (ap);
          printf (" expected %d observed %d)", Y->array[i], X->array[i]);
        }
      printf ("\n");
    }
  update_count_tests (status);
  return (status ? FALSE : TRUE);
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
      update_count_tests (1);
      return FALSE;
    }
  return pnl_test_array (X->array, Y->array, X->mn, relerr, pnl_cmp_eq, str, fmt, ap);
}

/** 
 * Checks if |x(i,j) - y(i,j)| / |y(i,j)| < relerr
 * 
 * @param X computed result (vector)
 * @param Y expected result (vector)
 * @param relerr relative error
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
      update_count_tests (1);
      return FALSE;
    }
  return pnl_test_array (X->array, Y->array, X->size, relerr, pnl_cmp_eq_rel, str, fmt, ap);
}

/** 
 * Checks if |x(i,j) - y(i,j)| < abserr
 * 
 * @param X computed result (vector)
 * @param Y expected result (vector)
 * @param abserr absolute error
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
      update_count_tests (1);
      return FALSE;
    }
  return pnl_test_array (X->array, Y->array, X->size, abserr, pnl_cmp_eq_abs, str, fmt, ap);
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
      update_count_tests (1);
      return FALSE;
    }
  return pnl_test_array (X->array, Y->array, X->size, relerr, pnl_cmp_eq, str, fmt, ap);
}

/** 
 * Checks if |x(i,j) - y(i,j)| < abserr
 * 
 * @param X computed result (vector)
 * @param Y expected result (vector)
 * @param abserr absolute error
 * @param str the name of the tested function
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return FALSE or TRUE
 */
int pnl_test_vect_complex_eq_abs (const PnlVectComplex *X, const PnlVectComplex *Y, double abserr, const char *str, const char *fmt, ...)
{
  va_list ap;
  va_start (ap, fmt);
  va_end (ap);
  if ( X->size != Y->size )
    {
      printf ("%s : ", str);
      printf ("FAIL (size mismatch");
      printf (fmt, ap); printf (")\n");
      update_count_tests (1);
      return FALSE;
    }
  return pnl_test_array ((double *)X->array, (double *)Y->array, 2 * X->size, abserr, pnl_cmp_eq_abs, str, fmt, ap);
}

/**
 * Checks if |x(i,j,...) - y(i,j,...)| / |y(i,j,...)| < relerr
 * 
 * @param X computed result (Hmatrix)
 * @param Y expected result (Hmatrix)
 * @param relerr relative error 
 * @param str the name of the tested function
 * @param fmt a forhmat string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return FALSE or TRUE
 */
int pnl_test_hmat_eq_rel (const PnlHmat *X, const PnlHmat *Y, double relerr, const char *str, const char *fmt, ...)
{
  int i;
  va_list ap;
  va_start (ap, fmt);
  va_end (ap);

  if ( X->ndim != Y->ndim ) goto dim_fail;
  for ( i=0 ; i<X->ndim ; i++ )
    {
      if ( X->dims[i] != Y->dims[i] ) goto dim_fail;
    }
  return pnl_test_array (X->array, Y->array, X->mn, relerr, pnl_cmp_eq_rel, str, fmt, ap);

dim_fail:
  printf ("%s : ", str);
  printf ("FAIL (size mismatch");
  printf (fmt, ap); printf (")\n");
  update_count_tests (1);
  return FALSE;
}

/** 
 * Checks if |x(i,j,...) - y(i,j,...)| < abserr
 * 
 * @param X computed result (Hmatrix)
 * @param Y expected result (Hmatrix)
 * @param abserr absolute error 
 * @param str the name of the tested function
 * @param fmt a forhmat string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return FALSE or TRUE
 */
int pnl_test_hmat_eq_abs (const PnlHmat *X, const PnlHmat *Y, double abserr, const char *str, const char *fmt, ...)
{
  int i;
  va_list ap;
  va_start (ap, fmt);
  va_end (ap);

  if ( X->ndim != Y->ndim ) goto dim_fail;
  for ( i=0 ; i<X->ndim ; i++ )
    {
      if ( X->dims[i] != Y->dims[i] ) goto dim_fail;
    }
  return pnl_test_array (X->array, Y->array, X->mn, abserr, pnl_cmp_eq_abs, str, fmt, ap);

dim_fail:
  printf ("%s : ", str);
  printf ("FAIL (size mismatch");
  printf (fmt, ap); printf (")\n");
  update_count_tests (1);
  return FALSE;
}

/** 
 * Checks if |x(i,j,...) - y(i,j,...)| / max (1, y(i,j,...) < abserr
 * 
 * @param X computed result (Hmatrix)
 * @param Y expected result (Hmatrix)
 * @param abserr  error (absolute)
 * @param str the name of the tested function
 * @param fmt a forhmat string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return FALSE or TRUE
 */
int pnl_test_hmat_eq (const PnlHmat *X, const PnlHmat *Y, double abserr, const char *str, const char *fmt, ...)
{
  int i;
  va_list ap;
  va_start (ap, fmt);
  va_end (ap);

  if ( X->ndim != Y->ndim ) goto dim_fail;
  for ( i=0 ; i<X->ndim ; i++ )
    {
      if ( X->dims[i] != Y->dims[i] ) goto dim_fail;
    }
  return pnl_test_array (X->array, Y->array, X->mn, abserr, pnl_cmp_eq, str, fmt, ap);

dim_fail:
  printf ("%s : ", str);
  printf ("FAIL (size mismatch");
  printf (fmt, ap); printf (")\n");
  update_count_tests (1);
  return FALSE;
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
