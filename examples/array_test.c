
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
#include "pnl/pnl_array.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_tridiag_matrix.h"
#include "tests_utils.h"



/** 
 * Test set and get for arrays
 */
static void array_ops ()
{
  PnlArray *T;
  PnlMat *M;
  PnlTridiagMat *tridiag;

  M = pnl_mat_create_from_scalar (2, 3, 3.5);
  tridiag = pnl_tridiag_mat_create_from_scalar (4, 0.5);
  T = pnl_array_create (2);
  pnl_array_set (T, 0, PNL_OBJECT(M));
  pnl_array_set (T, 1, PNL_OBJECT(tridiag));

  if ( pnl_array_get (T, 0) != PNL_OBJECT(M) || 
       pnl_array_get (T, 1) != PNL_OBJECT(tridiag))
    {
      pnl_test_set_fail ("pnl_array_{g,s}et", 0., 0.);
      goto J1;
    }
  pnl_test_set_ok ("pnl_array_{g,s}et");
J1:
  pnl_array_free (&T);
 }

/** 
 * Test the copy of an array using the copy member function of PnlObjects
 */
static void array_copy ()
{
  PnlArray *T, *C;
  PnlMat *M;
  PnlVect *v;

  M = pnl_mat_create_from_scalar (2, 3, 3.5);
  v = pnl_vect_create_from_scalar (5, M_SQRT2);
  T = pnl_array_create (2);
  pnl_array_set (T, 0, PNL_OBJECT(M));
  pnl_array_set (T, 1, PNL_OBJECT(v));

  C = pnl_array_copy (T);
  if ( C->size != T->size ||
       pnl_test_mat_eq ((PnlMat *) T->array[0], M, 1E-12, "", "") != TRUE ||
       pnl_test_vect_eq ((PnlVect *) T->array[1], v, 1E-12, "", "") != TRUE )
    {
      pnl_test_set_fail ("pnl_array_copy", 0., 0.);
      goto J1;
    }
  pnl_test_set_ok ("pnl_array_copy");
J1:
  pnl_array_free (&T);
  pnl_array_free (&C);
}

/** 
 * Test cloning an array using the clone member function of PnlObjects
 */
static void array_clone ()
{
  PnlArray *T, *C;
  PnlMat *M;
  PnlVect *v;

  M = pnl_mat_create_from_scalar (2, 3, 3.5);
  v = pnl_vect_create_from_scalar (5, M_SQRT2);
  T = pnl_array_create (2);
  pnl_array_set (T, 0, PNL_OBJECT(M));
  pnl_array_set (T, 1, PNL_OBJECT(v));

  C = pnl_array_new ();
  pnl_array_clone (C, T);
  if ( C->size != T->size ||
       pnl_test_mat_eq ((PnlMat *) T->array[0], M, 1E-12, "", "") != TRUE ||
       pnl_test_vect_eq ((PnlVect *) T->array[1], v, 1E-12, "", "") != TRUE )
    {
      pnl_test_set_fail ("pnl_array_copy", 0., 0.);
      goto J1;
    }
  pnl_test_set_ok ("pnl_array_copy");
J1:
  pnl_array_free (&T);
  pnl_array_free (&C);
}


int main (int argc, char **argv)
{
  pnl_test_init (argc, argv);
  array_ops ();
  array_copy ();
  array_clone ();
  exit (pnl_test_finalize ("Array operations"));
}
