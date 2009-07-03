
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
#include "tests.h"

static void all_test();
extern void random_test();
extern void cumulfunc_test();
extern void vector_test();
extern void matrix_test();
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

static list tests[] =
{
    MAKE_ENUM(1, all_test),
    MAKE_ENUM(2, random_test),
    MAKE_ENUM(3, cumulfunc_test),
    MAKE_ENUM(4, vector_test),
    MAKE_ENUM(5, matrix_test),
    MAKE_ENUM(6, speed_test),
    MAKE_ENUM(7, solver_iterativ_test),
    MAKE_ENUM(8, special_matrix_test),
    MAKE_ENUM(9, basis_test),
    MAKE_ENUM(10,dft_test),
    MAKE_ENUM(11,laplace_test),
    MAKE_ENUM(12,finance_function_test),
    MAKE_ENUM(13,integration_test),
    MAKE_ENUM(14,complex_test),
    MAKE_ENUM(15,bessel_test),    
    MAKE_ENUM(16,perm_test),    
    MAKE_ENUM(17,SpGrid_test),
    MAKE_ENUM(18,list_test),
    MAKE_ENUM(19,root_test),
    MAKE_ENUM(20,band_matrix_test),
    MAKE_ENUM(21,special_func_test),
    MAKE_ENUM(NULL_INT, NULL)
};
      
static void all_test()
{
    int len=0;
    while (tests[len].id != NULL_INT)
    {
        if (tests[len].func != all_test) (tests[len].func)();
        len ++;
    }
}


int main()
{
  int len=0, choice;
    while (tests[len].id != NULL_INT)
    {
        printf("%2d. %s\n",  tests[len].id, tests[len].label);
        len ++;
    }
    len--;
    printf("Which test do you want to run?\n");

    while(1)
    {
        scanf("%d", &choice); 
        if (choice <1 || choice > tests[len].id) printf("illegal choice\n");
        else break;
    }
    len = 0;
    while (tests[len].id != NULL_INT)
      {
        if (tests[len].id == choice)
          {
            (tests[len].func)();
            return 0;
          }
        len++;
      }
    printf("Can't find test %i\n", choice);
    return 0;
}
