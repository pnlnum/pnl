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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <stdarg.h>
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_complex.h"

static BASE FUNCTION(,_op_plus)(BASE a, BASE b) { return PLUS(a,b);}
static BASE FUNCTION(,_op_minus)(BASE a, BASE b) { return MINUS(a,b);}
static BASE FUNCTION(,_op_mult)(BASE a, BASE b) { return MULT(a,b); }
static BASE FUNCTION(,_op_div)(BASE a, BASE b) { return DIV(a,b); }
static BASE FUNCTION(,_op_inv)(BASE a) { return INV(a); }

static void FUNCTION(pnl_array,apply_op)(BASE *lhs,
                                         BASE  x, 
                                         BASE (*op)(BASE, BASE ),
                                         int n)
{
  int i=0;
  BASE *lptr= lhs;
  while(i<n)
    {
      (*lptr)= op(*lptr,x);
      i++;
      lptr++;
    }
}

/**
 * in-place map function compute 
 *
 * @param lhs left hand side array lhs(i) contains  f(lhs(i))
 * @param f the function to be applied term by term
 * @param n int size
 */
void FUNCTION(pnl_array,map_inplace)(BASE *lhs, 
                                     BASE(*f)(BASE ),
                                     int n)
{
  int i=0;
  BASE *lptr= lhs;
  while(i<n)
    {
      (*lptr)=(*f)(*lptr);
      i++;
      lptr++;
    }
}

/**
 * in-place array scalar addition
 *
 * @param lhs left hand side array  lhs(i) contains lhs(i) + x
 * @param x scalar
 * @param n int size
  */
void FUNCTION(pnl_array,CONCAT2(plus_,BASE))(BASE *lhs , BASE x, int n)
{
  FUNCTION(pnl_array,apply_op)(lhs, x,FUNCTION(,_op_plus),n);
}

/**
 * in-place array scalar substraction
 *
 * @param lhs left hand side array lhs(i) contains lhs(i) -x 
 * @param x scalar
 * @param n int size
 */
void FUNCTION(pnl_array,CONCAT2(minus_,BASE))(BASE *lhs , BASE x, int n)
{
  FUNCTION(pnl_array,apply_op)(lhs, x, FUNCTION(,_op_minus),n);
}

/**
 * in-place array scalar multiplication
 *
 * @param lhs left hand side array lhs(i) contains lhs(i) * x
 * @param x scalar
 * @param n int size
  */
void FUNCTION(pnl_array,CONCAT2(mult_,BASE))(BASE *lhs , BASE x, int n)
{
  FUNCTION(pnl_array,apply_op)(lhs, x, FUNCTION(,_op_mult),n);
}

/**
 * in-place array scalar division
 *
 * @param lhs left hand side array lhs(i) contains  lhs(i)/x
 * @param x scalar
 * @param n int size
 */
void FUNCTION(pnl_array,CONCAT2(div_,BASE))(BASE *lhs , BASE x, int n)
{
  FUNCTION(pnl_array,apply_op)(lhs, x, FUNCTION(,_op_div),n);
}

/**
 * map array componentwise
 *
 * @param lhs each component lhs(i) contains f(rhs(i))
 * @param rhs right hand side array
 * @param f real function 
 * @param n int size
 */
void FUNCTION(pnl_array,map)(BASE *lhs, const BASE *rhs, BASE(*f)(BASE),int n)
{
  BASE *lptr=lhs;
  const BASE *rptr=rhs;
  int i=0;
  while(i<n)
    {
      *lptr =(*f)(*rptr);
      i++;
      lptr++;
      rptr++;
    }
}

/**
 * map array componentwise
 *
 * @param lhs each component lhs(i) contains f(lhs(i),rhs(i))
 * @param rhs right hand side array
 * @param f real function 
 * @param n int size
 */
static void FUNCTION(pnl_array,map_array)(BASE *lhs, const BASE *rhs, BASE(*f)(BASE,BASE),int n)
{
  BASE *lptr=lhs;
  const BASE *rptr=rhs;
  int i=0;
  while(i<n)
    {
      *lptr =(*f)(*lptr, *rptr);
      i++;
      lptr++;
      rptr++;
    }
}

/**
 * in-place term by term array inverse
 *
 * @param lhs left hand side array lhs(i) contains = 1/ lhs(i)
 * @param n int size
 */
void FUNCTION(pnl_array,inv_term)(BASE *lhs,int n)
{
  FUNCTION(pnl_array,map_inplace)(lhs,FUNCTION(,_op_inv),n);
}

/**
 * in-place array term by term addition
 *
 * @param lhs left hand side array lhs(i) contians lhs(i)+rhs(i)
 * @param rhs right hand side array
 * @param n int size
 */
void FUNCTION(pnl_array,plus_array_term)(BASE *lhs, const BASE *rhs,int n)
{
  FUNCTION(pnl_array,map_array)(lhs, rhs, FUNCTION(,_op_plus),n);
}

/**
 * in-place array term by term substraction
 *
 * @param lhs left hand side array lhs(i) contains lhs(i)-rhs(i)
 * @param rhs right hand side array
 * @param n int size
 */
void FUNCTION(pnl_array,minus_array_term)(BASE *lhs, const BASE *rhs,int n)
{
  FUNCTION(pnl_array,map_array)(lhs, rhs, FUNCTION(,_op_minus),n);
}


/**
 * in-place term by term array inverse
 *
 * @param lhs left hand side array lhs contains lhs(i)=lhs(i)/rhs(i)
 * @param rhs right hand side array
 * @param n int size
 */
void FUNCTION(pnl_array,div_array_term)(BASE *lhs, const BASE *rhs,int n)
{
  FUNCTION(pnl_array,map_array)(lhs, rhs,FUNCTION(,_op_div),n);
}

/**
 * in-place array term by term multiplication
 *
 * @param lhs left hand side array lhs(i) contains lhs(i)*rhs(i)
 * @param rhs right hand side array
 * @param n int size
 */
void FUNCTION(pnl_array,mult_array_term)(BASE *lhs, const BASE *rhs,int n)
{
  FUNCTION(pnl_array,map_array)(lhs, rhs, FUNCTION(,_op_mult),n);
}


