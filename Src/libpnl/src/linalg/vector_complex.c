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

#include "config.h"
#include "pnl_vector.h"


/**
 * Creates a Complex Vector from the arrays of real an imaginary parts
 *
 * @param size  size of the vector
 * @param re array of the real parts
 * @param im array of the imaginary parts
 * @return v[i] = re[i] + I im[i]
 */
PnlVectComplex* pnl_vect_complex_create_from_array(int size, const double *re, const double *im)
{
  int i, i_re, i_im;
  PnlVectComplex *v = pnl_vect_complex_create (size);
  double *ptr = (double*) v->array;
  
  i_re = 0; i_im = 1;
  for (i=0; i<size; i++)
    {
      ptr[i_re] = re[i];
      ptr[i_im] = im[i];
      i_re += 2; i_im+=2;
    }
  return v;
}

/**
 * Stores a Complex Vector in two arrays of real an imaginary parts
 *
 * @param v a complex vector
 * @param re array of the real parts. Must already be allocated.
 * @param im array of the imaginary parts. Must already be allocated.
 * v[i] = re[i] + I im[i]
 */
void pnl_vect_complex_split_in_array(PnlVectComplex* v, double *re, double *im)
{
  int i, i_re, i_im;
  int size =v->size;
  double *ptr = (double*)v->array;
  
  i_re = 0; i_im = 1;
  for (i=0; i<size; i++)
    {
      re[i]=ptr[i_re] ;
      im[i]=ptr[i_im] ;
      i_re += 2; i_im+=2;
    }
}


/**
 * Stores a Complex Vector in two Real Vectors of real an imaginary parts
 *
 * @param v a complex vector
 * @param re PnlVect of the real parts. Must already be allocated.
 * @param im PnlVect of the imaginary parts. Must already be allocated.
 * v[i] = re[i] + I im[i]
 */
void pnl_vect_complex_split_in_vect(PnlVectComplex* v, PnlVect *re, PnlVect *im)
{
  pnl_vect_resize (re, v->size);
  pnl_vect_resize (im, v->size);
  pnl_vect_complex_split_in_array (v, re->array, im->array);
}


double pnl_vect_complex_get_real (const PnlVectComplex *v, int i)
{
  return ((double *)(v->array))[2*i];
}

double pnl_vect_complex_get_imag (const PnlVectComplex *v, int i)
{
  return ((double *)(v->array))[2*i+1];
}

double* pnl_vect_complex_lget_real (const PnlVectComplex *v, int i)
{
  return &(((double *)(v->array))[2*i]);
}

double* pnl_vect_complex_lget_imag (const PnlVectComplex *v, int i)
{
  return &(((double *)(v->array))[2*i+1]);
}

void pnl_vect_complex_set_real (const PnlVectComplex *v, int i, double re)
{
  ((double *)(v->array))[2*i] = re;
}

void pnl_vect_complex_set_imag (const PnlVectComplex *v, int i, double im)
{
  ((double *)(v->array))[2*i+1] = im;
}


/**
 * in-place vector scalar double multiplication
 *
 * @param lhs : left hand side vector
 * @param x : scalar double
 * @return  lhs = lhs*x
 */
void pnl_vect_complex_mult_double(PnlVectComplex *lhs , double x)
{
  dcomplex *lptr;
  int i=0;
  while(i<lhs->size)
    {
      lptr = pnl_vect_complex_lget(lhs, i);
      lptr->r *= x;
      lptr->i *= x;
      i++;
    }
}


