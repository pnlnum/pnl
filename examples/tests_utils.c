
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
#include "tests_utils.h"

int verbose = PNL_FALSE;
int count_tests;
int count_ok;
int count_fail;

/** 
 * Initialize test counting system
 */
void pnl_test_init (int argc, char **argv)
{
  if ( argc == 2 &&
       ( ( strcmp (argv[1], "-v") == 0 ) ||
         ( strcmp (argv[1], "--verbose") == 0 ) )
     )
    {
      verbose = PNL_TRUE;
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
  count_tests++;
  if (status == PNL_FALSE)
    {
      count_fail++;
    }
  else
    {
      count_ok++;
    }
}

/** 
 * Prints a summary of the tests
 * 
 * @return 0 if all tests were successful and 1 otherwise
 */
int pnl_test_finalize(const char *str)
{
  printf ("%s: %s (TOTAL: %d, PASSED: %d, FAILED: %d)\n", str, (count_fail == 0) ? "OK" : "FAIL", count_tests, count_ok, count_fail);
  return ( count_fail >0 );
}

int pnl_test_is_verbose ()
{
  return verbose;
}

/** 
 * Declare a test as passed
 * 
 * @param str a string
 */
void pnl_test_set_ok (const char *str)
{
  update_count_tests(PNL_TRUE);
  if (verbose) 
    {
      printf ("\t%s: OK\n", str);
    }
}

/** 
 * Declare a test as failed
 * 
 * @param str a string
 * @param res computed result
 * @param expected expected result
 */
void pnl_test_set_fail (const char *str, double res, double expected)
{
  update_count_tests(PNL_FALSE);
  printf ("\t%s: FAIL (observed %.18f expected %.18f)\n", str, res, expected);
}


/**
 * Declare a test as failed
 *
 * @param str a string
 */
void pnl_test_set_fail0(const char *str)
{
  update_count_tests(PNL_FALSE);
  printf ("\t%s: FAIL\n", str);
}


static int pnl_test_eq_aux (double x, double y, double relerr, int(*cmp)(double, double, double), const char *str, const char *fmt, va_list ap)
{
  int status = (*cmp)(x, y, relerr);
  if ((status == PNL_FALSE) || (verbose == PNL_TRUE))
    {
      printf("\t%s: ", str);
      printf(status == PNL_FALSE ? "FAIL" : "OK");
      if (!status)
        {
          printf (" (");
          vprintf (fmt, ap);
          printf (" expected %.18g observed %.18g)", y, x);
        }
      printf ("\n");
    }
  update_count_tests(status);
  return status;
}

/** 
 * Check if |x - y| / (max(1, |y|)) < relerr
 * 
 * @param x computed result
 * @param y expected result
 * @param relerr relative error (note that when |y| < 1, it is an absolute
 * error)
 * @param str the name of the tested function
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return PNL_FALSE or PNL_TRUE
 */
int pnl_test_eq(double x, double y, double relerr, const char *str, const char *fmt, ...)
{
  int status;
  va_list ap;
  va_start(ap, fmt);
  if (fabs(y) >= 1)
    status = pnl_test_eq_aux (x, y, relerr, pnl_isequal_rel, str, fmt, ap);
  else
    status = pnl_test_eq_aux (x, y, relerr, pnl_isequal_abs, str, fmt, ap);
  va_end(ap);
  return status;
}


/** 
 * Check if |x -y| / |y| < relerr
 * 
 * @param x computed result
 * @param y expected value
 * @param relerr maximum relative error
 * @param str the functionality tested
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return PNL_FALSE or PNL_TRUE
 */
int pnl_test_eq_rel (double x, double y, double relerr, const char *str, const char *fmt, ...)
{
  int status;
  va_list ap;
  va_start (ap, fmt);
  status = pnl_test_eq_aux (x, y, relerr, pnl_isequal_rel, str, fmt, ap);
  va_end (ap);
  return status;
}

/** 
 * Check if |x - y|  < abserr
 * 
 * @param x computed result
 * @param y expected value
 * @param abserr maximum absolute error
 * @param str the functionality tested
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return PNL_FALSE or PNL_TRUE
 */
int pnl_test_eq_abs (double x, double y, double abserr, const char *str, const char *fmt, ...)
{
  int status;
  va_list ap;
  va_start (ap, fmt);
  status = pnl_test_eq_aux (x, y, abserr, pnl_isequal_abs, str, fmt, ap);
  va_end (ap);
  return status;
}

/** 
 * Compare two vectors component-wise using the comparison function
 * specified b cmp
 * 
 * @param X a vector (computed result)
 * @param Y a vector (expected result)
 * @param n size of the expected result
 * @param relerr maximum admissible error in the comparison
 * @param cmp the comparison function
 * @param str name of the tested functionality
 * @param fmt a format string
 * @param ... extra arguments for @p fmt
 * 
 * @return PNL_TRUE or PNL_FALSE
 */
static int pnl_test_array (const double *X, const double *Y, int n, double relerr, int(*cmp)(double, double, double), const char *str, const char *fmt, ...)
{
  int i, status;
  va_list ap;
  status = 0;
  va_start (ap, fmt);
  for (i = 0; i < n; i++)
    {
      const double x = X[i];
      const double y = Y[i];
      status = (*cmp)(x, y, relerr);
      if (status == PNL_FALSE) break;
    }
  if ((status == PNL_FALSE) || (verbose == PNL_TRUE))
    {
      printf("\t%s: ", str);
      printf(status == PNL_TRUE ? "OK" : "FAIL");
      if (status == PNL_FALSE)
        {
          printf(" (");
          vprintf(fmt, ap);
          printf(" expected %.18g observed %.18g)", Y[i], X[i]);
        }
      printf("\n");
    }
  update_count_tests(status);
  return status;
}

/** 
 * Check if |x(i,j) - y(i,j)| < abserr
 * 
 * @param X computed result (matrix)
 * @param Y expected result (matrix)
 * @param abserr absolute error 
 * @param str the name of the tested function
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return PNL_FALSE or PNL_TRUE
 */
int pnl_test_mat_eq_abs (const PnlMat *X, const PnlMat *Y, double abserr, const char *str, const char *fmt, ...)
{
  int status;
  va_list ap;
  va_start (ap, fmt);
  if ( X->m != Y->m || X->n != Y->n )
    {
      printf ("%s: ", str);
      printf ("FAIL (size mismatch");
      printf (fmt, ap); printf (")\n");
      update_count_tests(PNL_FALSE);
      return PNL_FALSE;
    }
  status = pnl_test_array (X->array, Y->array, X->mn, abserr, pnl_isequal_abs, str, fmt, ap);
  va_end (ap);
  return status;
}

/** 
 * Check if |x(i,j) - y(i,j)| < abserr
 * 
 * @param X computed result (matrix)
 * @param Y expected result (matrix)
 * @param abserr absolute error 
 * @param str the name of the tested function
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return PNL_FALSE or PNL_TRUE
 */
int pnl_test_mat_complex_eq_abs (const PnlMatComplex *X, const PnlMatComplex *Y, double abserr, const char *str, const char *fmt, ...)
{
  va_list ap;
  va_start (ap, fmt);
  va_end (ap);
  if ( X->m != Y->m || X->n != Y->n )
    {
      printf ("%s: ", str);
      printf ("FAIL (size mismatch");
      printf (fmt, ap); printf (")\n");
      update_count_tests(PNL_FALSE);
      return PNL_FALSE;
    }
  return pnl_test_array ((double *)X->array, (double *)Y->array, 2 * X->mn, abserr, pnl_isequal_abs, str, fmt, ap);
}

/** 
 * Check if x(i,j) = y(i,j) for integer matrices
 * 
 * @param X computed result (matrix)
 * @param Y expected result (matrix)
 * @param str the name of the tested function
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return PNL_FALSE or PNL_TRUE
 */
int pnl_test_mat_int_eq(const PnlMatInt *X, const PnlMatInt *Y, const char *str, const char *fmt, ...)
{
  int status;
  va_list ap;
  va_start(ap, fmt);
  status = pnl_mat_int_isequal(X, Y);
  update_count_tests(status);
  if (status == PNL_FALSE || verbose == PNL_TRUE)
    {
      printf("\t%s: ", str);
      printf(status ? "OK" : "FAIL");
      if (!status)
        {
          printf(" (");
          vprintf(fmt, ap);
          printf("\n  expected:\n");
          pnl_mat_int_print(Y);
          printf("  observed:\n");
          pnl_mat_int_print(X);
          printf(")");
        }
      printf("\n");
    }
  va_end(ap);
  return status;
}

/** 
 * Check if x(i,j) = y(i,j) for integer sparse matrices
 * 
 * @param X computed result (matrix)
 * @param Y expected result (matrix)
 * @param str the name of the tested function
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return PNL_FALSE or PNL_TRUE
 */
int pnl_test_sp_mat_int_eq(const PnlSpMatInt *X, const PnlSpMatInt *Y, const char *str, const char *fmt, ...)
{
  int status;
  va_list ap;
  va_start(ap, fmt);
  status = pnl_sp_mat_int_isequal(X, Y);
  update_count_tests(status);
  if (status == PNL_FALSE || verbose == PNL_TRUE)
    {
      printf("\t%s: ", str);
      printf(status ? "OK" : "FAIL");
      if (!status)
        {
          printf(" (");
          vprintf(fmt, ap);
          printf("\n  expected:\n");
          pnl_sp_mat_int_print(Y);
          printf("  observed:\n");
          pnl_sp_mat_int_print(X);
          printf(")");
        }
      printf("\n");
    }
  va_end(ap);
  return status;
}

/** 
 * Check if |x(i,j) - y(i,j)| < abserr
 * 
 * @param X computed result (vector)
 * @param Y expected result (vector)
 * @param abserr absolute error
 * @param str the name of the tested function
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return PNL_FALSE or PNL_TRUE
 */
int pnl_test_vect_eq_abs (const PnlVect *X, const PnlVect *Y, double abserr, const char *str, const char *fmt, ...)
{
  int status;
  va_list ap;
  va_start (ap, fmt);
  status = pnl_vect_isequal_abs(X, Y, abserr);
  update_count_tests(status);
  if (!status || verbose)
    {
      printf("\t%s: ", str);
      printf(status ? "OK" : "FAIL");
      if (!status)
        {
          printf(" (");
          vprintf(fmt, ap);
          printf("\n  expected; ");
          pnl_vect_print_asrow(Y);
          printf("  observed; ");
          pnl_vect_print_asrow(X);
          printf(")");
        }
      printf("\n");
    }
  va_end(ap);
  return status;
}

/** 
 * Check if |x(i,j) - y(i,j)| < abserr
 * 
 * @param X computed result (vector)
 * @param Y expected result (vector)
 * @param abserr absolute error
 * @param str the name of the tested function
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return PNL_FALSE or PNL_TRUE
 */
int pnl_test_vect_complex_eq_abs (const PnlVectComplex *X, const PnlVectComplex *Y, double abserr, const char *str, const char *fmt, ...)
{
  int status;
  va_list ap;
  va_start (ap, fmt);
  status = pnl_vect_complex_isequal_abs(X, Y, abserr);
  update_count_tests(status);
  if (!status || verbose)
    {
      printf("\t%s: ", str);
      printf(status ? "OK" : "FAIL");
      if (!status)
        {
          printf(" (");
          vprintf(fmt, ap);
          printf("\n  expected; ");
          pnl_vect_complex_print_asrow(Y);
          printf("  observed; ");
          pnl_vect_complex_print_asrow(X);
          printf(")");
        }
      printf("\n");
    }
  va_end(ap);
  return status;
}

/** 
 * Check if |x(i,j,...) - y(i,j,...)| < abserr
 * 
 * @param X computed result (Hmatrix)
 * @param Y expected result (Hmatrix)
 * @param abserr absolute error 
 * @param str the name of the tested function
 * @param fmt a format string to be passed to printf
 * @param ... extra arguments for printf
 * 
 * @return PNL_FALSE or PNL_TRUE
 */
int pnl_test_hmat_eq_abs (const PnlHmat *X, const PnlHmat *Y, double abserr, const char *str, const char *fmt, ...)
{
  int status;
  int i;
  va_list ap;
  va_start (ap, fmt);

  if ( X->ndim != Y->ndim ) goto dim_fail;
  for ( i=0 ; i<X->ndim ; i++ )
    {
      if ( X->dims[i] != Y->dims[i] ) goto dim_fail;
    }
  status = pnl_test_array (X->array, Y->array, X->mn, abserr, pnl_isequal_abs, str, fmt, ap);
  va_end(ap);
  return status;

dim_fail:
  printf ("%s: ", str);
  printf ("FAIL (size mismatch");
  printf (fmt, ap); printf (")\n");
  update_count_tests(PNL_FALSE);
  va_end (ap);
  return PNL_FALSE;
}

int pnl_test_basis_eq(const PnlBasis *observed, const PnlBasis *expected)
{
  if (
    expected->id != observed->id ||
    expected->nb_func != observed->nb_func ||
    expected->nb_variates != observed->nb_variates ||
    expected->isreduced != observed->isreduced ||
    expected->f_params_size != observed->f_params_size
  )
    {
      pnl_test_set_fail0("Bases differ on first level parameters");
      return PNL_FALSE;
    }
  if (expected->isreduced)
    {
      pnl_test_array(expected->center, observed->center, expected->nb_variates, 1E-12, pnl_isequal_abs, "test_basis center", "");
      pnl_test_array(expected->scale, observed->scale, expected->nb_variates, 1E-12, pnl_isequal_abs, "test_basis scale", "");
      return PNL_FALSE;
    }
  if (memcmp(expected->f_params, observed->f_params, expected->f_params_size) != 0)
    {
      pnl_test_set_fail0("PnlBasis: f_params differ");
      return PNL_FALSE;
    }
  if (expected->SpT != NULL && observed->SpT != NULL)
    {
      return pnl_test_sp_mat_int_eq(expected->SpT, observed->SpT, "basis test tensor", "");
    }
  if (expected->SpT == NULL && observed->SpT == NULL)
    {
      pnl_test_set_ok("basis test local");
      return PNL_TRUE;
    }
  else
    {
      pnl_test_set_fail0("basis test");
      return PNL_FALSE;
    }
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
