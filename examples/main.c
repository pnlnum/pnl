
/************************************************************************/
/* Copyright Jérôme Lelong <jerome.lelong@gmail.com>                    */
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
#include <string.h>
#include "tests.h"

static void all_test();
extern void random_test();
extern void cumulfunc_test();
extern void vector_test();
extern void matrix_test();
extern void lapack_test();
extern void speed_test();
extern void solver_iterativ_test();
extern void basis_test();
extern void dft_test();
extern void laplace_test();
extern void finance_function_test();
extern void integration_test();
extern void complex_test();
extern void mathtools_test();
extern void root_test();
extern void special_func_test();
extern void tridiag_matrix_test ();
extern void band_matrix_test ();

static tst_list tests[] =
  {
    MAKE_ENUM(all_test),
    MAKE_ENUM(random_test),
    MAKE_ENUM(cumulfunc_test),
    MAKE_ENUM(vector_test),
    MAKE_ENUM(matrix_test),
    MAKE_ENUM(lapack_test),
    MAKE_ENUM(speed_test),
    MAKE_ENUM(solver_iterativ_test),
    MAKE_ENUM(tridiag_matrix_test),
    MAKE_ENUM(band_matrix_test),
    MAKE_ENUM(basis_test),
    MAKE_ENUM(dft_test),
    MAKE_ENUM(laplace_test),
    MAKE_ENUM(finance_function_test),
    MAKE_ENUM(integration_test),
    MAKE_ENUM(mathtools_test),
    MAKE_ENUM(complex_test),
    MAKE_ENUM(root_test),
    MAKE_ENUM(special_func_test),
    MAKE_ENUM(NULL)
  };

void run_all_test (tst_list *l)
{
  int len=0;
  while (l[len].func != NULL)
    {
      if (strcmp (l[len].label, "all_test") != 0) (l[len].func)();
      len ++;
    }
}


static void all_test()
{
  run_all_test (tests);
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

int main ()
{
  menu_test (tests);
  return 0;
}
