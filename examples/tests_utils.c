
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
#include "tests.h"

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
