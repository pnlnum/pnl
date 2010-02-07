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

#include "config.h"
#include "pnl_vector_uint.h"
#include "pnl_mathtools.h"

uint log2uint(uint x)
{
  uint level, y;
  if(x==0){printf("error in log2uint ");abort();}
  level=0;
  y=x;
  while(y>1)
    {
      level++;
      y>>=1;
    }
  return level;
}

double pnl_dyadic_cast(uint i)
{
  /* \frac{2 (i-2^(lev))+1 }{2^{lev+1}} */
  int l=log2uint(i);
  return (double)(((i-(1<<l))<<1)+1)/ 
    (double)(2<<l);} 

void pnl_vect_uint_dyadic_cast(const PnlVectUint * v_in,PnlVect * v_out)
{
  int i;
  for(i=0;i<v_in->size;i++)
    LET(v_out,i)=pnl_dyadic_cast(v_in->array[i]);
}


int pnl_vect_uint_less(const PnlVectUint * a,const PnlVectUint * b)
{
  int i;
  for( i=0;i<a->size;i++)
    if(a->array[i]!=b->array[i])
      return a->array[i]<b->array[i];
  return false;
}

int pnl_vect_uint_equal(const PnlVectUint * a,const PnlVectUint * b)
{
  int i;
  for(i=0;i<a->size;i++)
    if(a->array[i]!=b->array[i])
      return false;
  return true;
}

