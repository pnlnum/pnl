
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

#include "pnl_list.h"
#include "pnl_matrix.h"
#include "pnl_tridiag_matrix.h"


static void list_ops ()
{
  PnlList *L, *L1, *L2;
  PnlMat *M;
  PnlTridiagMat *T;

  printf ("--> List operations : insert first\n");
  
  L = pnl_list_new ();
  M = pnl_mat_create_from_double (2, 3, 3.5);
  T = pnl_tridiag_mat_create_from_double (4, 0.5);
  pnl_list_insert_first (L, PNL_OBJECT(M));
  pnl_list_insert_first (L, PNL_OBJECT(T));
  pnl_list_free (&L);

  printf ("--> List operations : insert last\n");
  
  L = pnl_list_new ();
  M = pnl_mat_create_from_double (2, 3, 3.5);
  T = pnl_tridiag_mat_create_from_double (4, 0.5);
  pnl_list_insert_last (L, PNL_OBJECT(M));
  pnl_list_insert_last (L, PNL_OBJECT(T));
  pnl_list_free (&L);

  printf ("--> List operations : remove last\n");
  
  L = pnl_list_new ();
  M = pnl_mat_create_from_double (2, 3, 3.5);
  T = pnl_tridiag_mat_create_from_double (4, 0.5);
  pnl_list_insert_last (L, PNL_OBJECT(M));
  pnl_list_insert_last (L, PNL_OBJECT(T));
  pnl_list_remove_last (L);
  pnl_list_free (&L);

  printf ("--> List operations : remove first\n");
  
  L = pnl_list_new ();
  M = pnl_mat_create_from_double (2, 3, 3.5);
  T = pnl_tridiag_mat_create_from_double (4, 0.5);
  pnl_list_insert_last (L, PNL_OBJECT(M));
  pnl_list_insert_last (L, PNL_OBJECT(T));
  pnl_list_print (L);
  pnl_list_remove_first (L);
  pnl_list_remove_last (L);
  pnl_list_free (&L);

  printf ("--> List operations : concat\n");

  L1 = pnl_list_new ();
  L2 = pnl_list_new ();
  M = pnl_mat_create_from_double (2, 3, 3.5);
  T = pnl_tridiag_mat_create_from_double (4, 0.5);
  pnl_list_insert_last (L1, PNL_OBJECT(M));
  pnl_list_insert_last (L1, PNL_OBJECT(T));
  pnl_list_print (L1);
  M = pnl_mat_create_from_double (2, 3, 3.5);
  T = pnl_tridiag_mat_create_from_double (4, 0.5);
  pnl_list_insert_first (L2, PNL_OBJECT(M));
  pnl_list_insert_first (L2, PNL_OBJECT(T));
  pnl_list_print (L2);

  pnl_list_concat(L1, L2);
  pnl_list_print(L1);
  pnl_list_free (&L1);
}


void list_test ()
{
  list_ops ();
}
