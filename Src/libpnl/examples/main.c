
/*************************************************************************/
/* Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>            */
/*                                                                       */
/* This program is free software; you can redistribute it and/or modify  */
/* it under the terms of the GNU General Public License as published by  */
/* the Free Software Foundation; either version 3 of the License, or     */
/* (at your option) any later version.                                   */
/*                                                                       */
/* This program is distributed in the hope that it will be useful,       */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/* GNU General Public License for more details.                          */
/*                                                                       */
/* You should have received a copy of the GNU General Public License     */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/*************************************************************************/

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
extern void special_matrix_test();
extern void basis_test();
extern void dft_test();
extern void laplace_test();
extern void finance_function_test();
extern void SpGrid_test();
extern void integration_test();
extern void complex_test();
extern void bessel_test();
extern void perm_test();
extern void list_test();
extern void root_test();
extern void band_matrix_test();
extern void special_func_test();

static tst_list tests[] =
  {
    MAKE_ENUM(1, all_test),
    MAKE_ENUM(2, random_test),
    MAKE_ENUM(3, cumulfunc_test),
    MAKE_ENUM(4, vector_test),
    MAKE_ENUM(5, matrix_test),
    MAKE_ENUM(6, lapack_test),
    MAKE_ENUM(7, speed_test),
    MAKE_ENUM(8, solver_iterativ_test),
    MAKE_ENUM(9, special_matrix_test),
    MAKE_ENUM(10, basis_test),
    MAKE_ENUM(11,dft_test),
    MAKE_ENUM(12,laplace_test),
    MAKE_ENUM(13,finance_function_test),
    MAKE_ENUM(14,integration_test),
    MAKE_ENUM(15,complex_test),
    MAKE_ENUM(16,bessel_test),    
    MAKE_ENUM(17,perm_test),    
    MAKE_ENUM(18,SpGrid_test),
    MAKE_ENUM(19,list_test),
    MAKE_ENUM(20,root_test),
    MAKE_ENUM(21,band_matrix_test),
    MAKE_ENUM(22,special_func_test),
    MAKE_ENUM(NULL_INT, NULL)
  };

void run_all_test (tst_list *l)
{
  int len=0;
  while (l[len].id != NULL_INT)
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

  while (l[len].id != NULL_INT)
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
