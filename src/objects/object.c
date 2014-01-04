
/*************************************************************************/
/* Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>            */
/*                                                                       */
/* This program is free software; you can redistribute it and/or modify  */
/* it under the terms of the GNU Lesser General Public License as        */ 
/* published by  the Free Software Foundation; either version 3 of the   */
/* License, or (at your option) any later version.                       */
/*                                                                       */
/* This program is distributed in the hope that it will be useful,       */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/* GNU Lesser General Public License for more details.                   */
/*                                                                       */
/* You should have received a copy of the GNU Lesser General Public      */
/* License  along with this program.  If not, see                        */
/* <http://www.gnu.org/licenses/>.                                       */ 
/*************************************************************************/


#include "pnl/pnl_object.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_sp_matrix.h"
#include "pnl/pnl_tridiag_matrix.h"
#include "pnl/pnl_band_matrix.h"
#include "pnl/pnl_basis.h"
#include "pnl/pnl_random.h"
#include "pnl/pnl_list.h"
#include "pnl/pnl_array.h"
#include "pnl/pnl_linalgsolver.h"


/**
 * Create a new PnlObject from a PnlType
 * 
 * @param type is one the PNL_TYPE_xxx values used to create a given specialised
 * PnlObject
 * 
 * @return a PnlObject
 */
PnlObject* pnl_object_create (PnlType type)
{
  PnlObject *o;
  /* Do not change the order of the cases, because types are grouped together */
  switch (type)
    {
    case PNL_TYPE_VECTOR:
      o = PNL_OBJECT(pnl_vect_object_new ()); break;
    case PNL_TYPE_VECTOR_DOUBLE:
      o = PNL_OBJECT(pnl_vect_new ()); break;
    case PNL_TYPE_VECTOR_INT:
      o = PNL_OBJECT(pnl_vect_int_new ()); break;
    case PNL_TYPE_VECTOR_COMPLEX:
      o = PNL_OBJECT(pnl_vect_complex_new ()); break;
    case PNL_TYPE_VECTOR_COMPACT:
      o = PNL_OBJECT(pnl_vect_object_new ()); break;
    case PNL_TYPE_MATRIX:
      o = PNL_OBJECT(pnl_mat_object_new ()); break;
    case PNL_TYPE_MATRIX_DOUBLE:
      o = PNL_OBJECT(pnl_mat_new ()); break;
    case PNL_TYPE_MATRIX_INT:
      o = PNL_OBJECT(pnl_mat_int_new ()); break;
    case PNL_TYPE_MATRIX_COMPLEX:
      o = PNL_OBJECT(pnl_mat_complex_new ()); break;
    case PNL_TYPE_TRIDIAG_MATRIX:
      o = PNL_OBJECT(pnl_tridiag_mat_object_new ()); break;
    case PNL_TYPE_TRIDIAG_MATRIX_DOUBLE:
      o = PNL_OBJECT(pnl_tridiag_mat_new ()); break;
    case PNL_TYPE_TRIDIAG_MATRIX_LU:
      o = PNL_OBJECT(pnl_tridiag_mat_lu_object_new ()); break;
    case PNL_TYPE_TRIDIAG_MATRIX_LU_DOUBLE:
      o = PNL_OBJECT(pnl_tridiag_mat_lu_new ()); break;
    case PNL_TYPE_BAND_MATRIX:
      o = PNL_OBJECT(pnl_band_mat_object_new ()); break;
    case PNL_TYPE_BAND_MATRIX_DOUBLE:
      o = PNL_OBJECT(pnl_band_mat_new ()); break;
    case PNL_TYPE_SP_MATRIX:
      o = PNL_OBJECT(pnl_sp_mat_object_new ()); break;
    case PNL_TYPE_SP_MATRIX_DOUBLE:
      o = PNL_OBJECT(pnl_sp_mat_new ()); break;
    case PNL_TYPE_SP_MATRIX_INT:
      o = PNL_OBJECT(pnl_sp_mat_int_new ()); break;
    case PNL_TYPE_SP_MATRIX_COMPLEX:
      o = PNL_OBJECT(pnl_sp_mat_complex_new ()); break;
    case PNL_TYPE_HMATRIX:
      o = PNL_OBJECT(pnl_hmat_object_new ()); break;
    case PNL_TYPE_HMATRIX_DOUBLE:
      o = PNL_OBJECT(pnl_hmat_new ()); break;
    case PNL_TYPE_HMATRIX_INT:
      o = PNL_OBJECT(pnl_hmat_int_new ()); break;
    case PNL_TYPE_HMATRIX_COMPLEX:
      o = PNL_OBJECT(pnl_hmat_complex_new ()); break;
    case PNL_TYPE_BASIS:
      o = PNL_OBJECT(pnl_basis_new ());
    case PNL_TYPE_RNG:
      o = PNL_OBJECT(pnl_rng_new()); break;
    case PNL_TYPE_LIST:
      o = PNL_OBJECT(pnl_list_new()); break;
    case PNL_TYPE_ARRAY:
      o = PNL_OBJECT(pnl_array_new()); break;
    case PNL_TYPE_CG_SOLVER:
      o = PNL_OBJECT(pnl_cg_solver_new ()); break;
    case PNL_TYPE_BICG_SOLVER:
      o = PNL_OBJECT(pnl_bicg_solver_new ()); break;
    case PNL_TYPE_GMRES_SOLVER:
      o = PNL_OBJECT(pnl_gmres_solver_new ()); break;
    default:
      o = NULL;
    }
  return o;
}

static int MessageIsOn = FALSE;
void pnl_message_on () { MessageIsOn = TRUE; }
void pnl_message_off () { MessageIsOn = FALSE; }
int pnl_message_is_on ()
{
  return ( MessageIsOn == TRUE ) ? TRUE : FALSE;
}
