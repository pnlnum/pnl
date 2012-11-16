/************************************************************************/
/* Copyright David Pommier <pommier.david@gmail.com>                    */
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

#include "pnl/pnl_config.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_mathtools.h"

int pnl_vect_int_less(const PnlVectInt * a,const PnlVectInt * b)
{
  int i;
  for( i=0;i<a->size;i++)
    if(a->array[i]!=b->array[i])
      return a->array[i]<b->array[i];
  return FALSE;
}

int pnl_vect_int_equal(const PnlVectInt * a,const PnlVectInt * b)
{
  int i;
  for(i=0;i<a->size;i++)
    if(a->array[i]!=b->array[i])
      return FALSE;
  return TRUE;
}

