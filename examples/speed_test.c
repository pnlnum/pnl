
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
#include <time.h>

#include "pnl/pnl_matrix.h"
#include "pnl/pnl_basis.h"
#include "pnl/pnl_random.h"
#include "tests_utils.h"


extern void speed_matrix_test ();
extern void speed_basis_test ();

static
tst_list speed_tests[] =
{
  MAKE_ENUM(speed_matrix_test),
  MAKE_ENUM(speed_basis_test),
  MAKE_ENUM(NULL)
};

int main ()
{
  menu_test(speed_tests); 
  return OK;
}


