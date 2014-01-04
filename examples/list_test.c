
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

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_list.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_tridiag_matrix.h"
#include "tests_utils.h"

static void list_insert_test ()
{
  PnlList *L;
  PnlMat *M;
  PnlTridiagMat *T;

  L = pnl_list_new ();
  M = pnl_mat_create_from_scalar (2, 3, 3.5);
  T = pnl_tridiag_mat_create_from_scalar (4, 0.5);
  pnl_list_insert_first (L, PNL_OBJECT(M));
  pnl_list_insert_first (L, PNL_OBJECT(T));

  if ( pnl_list_get (L, 0) != PNL_OBJECT(T)  || 
       pnl_list_get (L, 1) != PNL_OBJECT(M)  )
    {
      pnl_test_set_fail ("pnl_list_insert_firt", 0., 0.);
    }
  else
    {
      pnl_test_set_ok ("pnl_list_insert_firt");
    }
  pnl_list_free (&L);

  L = pnl_list_new ();
  M = pnl_mat_create_from_scalar (2, 3, 3.5);
  T = pnl_tridiag_mat_create_from_scalar (4, 0.5);
  pnl_list_insert_last (L, PNL_OBJECT(M));
  pnl_list_insert_last (L, PNL_OBJECT(T));
  if ( pnl_list_get (L, 0) != PNL_OBJECT(M)  || 
       pnl_list_get (L, 1) != PNL_OBJECT(T)  )
    {
      pnl_test_set_fail ("pnl_list_insert_last", 0., 0.);
    }
  else
    {
      pnl_test_set_ok ("pnl_list_insert_last");
    }
  pnl_list_free (&L);


  L = pnl_list_new ();
  M = pnl_mat_create_from_scalar (2, 3, 3.5);
  T = pnl_tridiag_mat_create_from_scalar (4, 0.5);
  pnl_list_insert_last (L, PNL_OBJECT(M));
  pnl_list_insert_last (L, PNL_OBJECT(T));
  pnl_list_remove_last (L);
  if ( pnl_list_get (L, 0) != PNL_OBJECT(M)  || 
       L->len != 1)
    {
      pnl_test_set_fail ("pnl_list_remove_last", 0., 0.);
    }
  else
    {
      pnl_test_set_ok ("pnl_list_remove_last");
    }
  pnl_list_free (&L);


  L = pnl_list_new ();
  M = pnl_mat_create_from_scalar (2, 3, 3.5);
  T = pnl_tridiag_mat_create_from_scalar (4, 0.5);
  pnl_list_insert_last (L, PNL_OBJECT(M));
  pnl_list_insert_last (L, PNL_OBJECT(T));
  pnl_list_remove_first (L);
  if ( pnl_list_get (L, 0) != PNL_OBJECT(T)  || 
       L->len != 1)
    {
      pnl_test_set_fail ("pnl_list_remove_first", 0., 0.);
    }
  else
    {
      pnl_test_set_ok ("pnl_list_remove_first");
    }
  pnl_list_free (&L);

}

static void list_remove_test ()
{
  PnlList *L;
  PnlMat *M;
  PnlTridiagMat *T;

  L = pnl_list_new ();
  M = pnl_mat_create_from_scalar (2, 3, 3.5);
  T = pnl_tridiag_mat_create_from_scalar (4, 0.5);
  pnl_list_insert_last (L, PNL_OBJECT(M));
  pnl_list_insert_last (L, PNL_OBJECT(T));
  pnl_list_remove_last (L);
  if ( pnl_list_get (L, 0) != PNL_OBJECT(M)  || 
       L->len != 1)
    {
      pnl_test_set_fail ("pnl_list_remove_last", 0., 0.);
    }
  else
    {
      pnl_test_set_ok ("pnl_list_remove_last");
    }
  pnl_list_free (&L);


  L = pnl_list_new ();
  M = pnl_mat_create_from_scalar (2, 3, 3.5);
  T = pnl_tridiag_mat_create_from_scalar (4, 0.5);
  pnl_list_insert_last (L, PNL_OBJECT(M));
  pnl_list_insert_last (L, PNL_OBJECT(T));
  pnl_list_remove_first (L);
  if ( pnl_list_get (L, 0) != PNL_OBJECT(T)  || 
       L->len != 1)
    {
      pnl_test_set_fail ("pnl_list_remove_first", 0., 0.);
    }
  else
    {
      pnl_test_set_ok ("pnl_list_remove_first");
    }
  pnl_list_free (&L);
}

static void list_concat_test ()
{
  PnlList *L1, *L2;
  PnlMat *M1, *M2;
  PnlTridiagMat *T1, *T2;

  L1 = pnl_list_new ();
  L2 = pnl_list_new ();
  M1 = pnl_mat_create_from_scalar (2, 3, 3.5);
  T1 = pnl_tridiag_mat_create_from_scalar (4, 0.5);
  pnl_list_insert_last (L1, PNL_OBJECT(M1));
  pnl_list_insert_last (L1, PNL_OBJECT(T1));
  M2 = pnl_mat_create_from_scalar (2, 3, 3.5);
  T2 = pnl_tridiag_mat_create_from_scalar (4, 0.5);
  pnl_list_insert_first (L2, PNL_OBJECT(M2));
  pnl_list_insert_first (L2, PNL_OBJECT(T2));

  pnl_list_concat(L1, L2);
  if  (L1->len != 4 ||
       pnl_list_get (L1, 0) != PNL_OBJECT(M1)  || 
       pnl_list_get (L1, 1) != PNL_OBJECT(T1)  || 
       pnl_list_get (L1, 2) != PNL_OBJECT(T2)  || 
       pnl_list_get (L1, 3) != PNL_OBJECT(M2)  )
       
    {
      pnl_test_set_fail ("pnl_list_concat", 0., 0.);
    }
  else
    {
      pnl_test_set_ok ("pnl_list_concat");
    }
  pnl_list_free (&L1);
}

static void list_copy_test ()
{
  PnlList *L1, *L2;
  PnlMat *M1;
  PnlVect *T1;

  L1 = pnl_list_new ();
  M1 = pnl_mat_create_from_scalar (2, 3, 3.5);
  T1 = pnl_vect_create_from_scalar (4, 0.5);
  pnl_list_insert_last (L1, PNL_OBJECT(M1));
  pnl_list_insert_last (L1, PNL_OBJECT(T1));

  L2 = pnl_list_copy(L1);
  if  (L2->len != 2 ||
       pnl_test_mat_eq ((PnlMat *)pnl_list_get (L2, 0),
                        M1, 1E-12, "", "") != TRUE  || 
       pnl_test_vect_eq ((PnlVect *)pnl_list_get (L2, 1),
                         T1, 1E-12, "", "") != TRUE )
    {
      pnl_test_set_fail ("pnl_list_copy", 0., 0.);
    }
  else
    {
      pnl_test_set_ok ("pnl_list_copy");
    }
  pnl_list_free (&L1);
  pnl_list_free (&L2);
}

int main (int argc, char **argv)
{
  pnl_test_init (argc, argv);
  list_insert_test ();
  list_remove_test ();
  list_concat_test ();
  list_copy_test ();
  exit (pnl_test_finalize ("List operations"));
}
