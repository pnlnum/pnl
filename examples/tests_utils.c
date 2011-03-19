
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
void pnl_init_tests ()
{
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
int pnl_finalize_tests (const char *str)
{
  printf ("%s : %s (TOTAL : %d, PASSED : %d, FAILED %d)\n", str, (count_fail == 0) ? "OK" : "FAIL", count_tests, count_ok, count_fail);
  return ( count_fail >0 );
}

/** 
 * Checks if |x -y| / |y| < relerr
 * 
 * @param x computed result
 * @param y exepcted value
 * @param relerr maximum relative error
 * @param str the fonctionnality tested
 * 
 * @return FALSE or TRUE
 */
int pnl_eq_rel (double x, double y, double relerr, const char *str, const char *fmt, ...)
{
  int status = 0;
  if ( (pnl_isnan (x) && !pnl_isnan(y)) ||
       (pnl_isinf (x) && !pnl_isinf(y)) )
    {
      status = 1;
    }
  else if ( y != 0 )
    {
      status = (fabs (x -y) / fabs(y) > relerr);
    }
  else
    {
      status = (fabs(y) > relerr);
    }
  if ( status || verbose )
    {
      printf ("%s : ", str);
      printf ( status ? "FAIL" : "OK");
      if ( status ) 
        {
          va_list ap;
          va_start (ap, fmt);
          printf (" (");
          vprintf (fmt, ap);
          va_end (ap);
          printf ("  expected %.18g observed %.18g)", y, x);
        }
      printf ("\n");
    }
  update_count_tests (status);
  return (status ? FALSE : TRUE);
}

/** 
 * Checks if |x - y|  < abserr
 * 
 * @param x computed result
 * @param y exepcted value
 * @param abserr maximum absolute error
 * @param str the fonctionnality tested
 * 
 * @return FALSE or TRUE
 */
int pnl_eq_abs (double x, double y, double abserr, const char *str, const char *fmt, ...)
{
  int status = 0;
  if ( (pnl_isnan (x) && !pnl_isnan(y)) ||
       (pnl_isinf (x) && !pnl_isinf(y)) )
    {
      status = 1;
    }
  else 
    {
      status = (fabs (x -y) > abserr);
    }
  if ( status || verbose )
    {
      printf ("%s : ", str);
      printf ( status ? "FAIL" : "OK");
      if ( status ) 
        {
          va_list ap;
          va_start (ap, fmt);
          printf (" (");
          vprintf (fmt, ap);
          va_end (ap);
          printf ("  expected %.18g observed %.18g)", y, x);
        }
      printf ("\n");
    }
  update_count_tests (status);
  return (status ? FALSE : TRUE);
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
