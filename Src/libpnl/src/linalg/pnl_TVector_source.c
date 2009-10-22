/*
 * This file contains two parts, each with its own copyright.
 */

/*************************************************************************/
/* Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com> 2007       */
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


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <stdarg.h>
#include "pnl_mathtools.h"


/**
 * creates a new TYPE(PnlVect)pointer.
 *
 * @param size the size of the array
 * @return  a TYPE(PnlVect)pointer
 */
TYPE(PnlVect) * FUNCTION(pnl_vect,create)(int size)
{
  TYPE(PnlVect)* v;
  if ((v = malloc (sizeof (TYPE(PnlVect)))) == NULL) return NULL;
  v->size = size;
  v->mem_size = size;
  v->owner = 1;
  if (size > 0)
    {
      if ((v->array = malloc (size * sizeof (BASE))) == NULL) return NULL;
    }
  else
    v->array = (BASE*) NULL;
  return v;
}


/**
 * creates a new TYPE(PnlVect)pointer.
 *
 * @param size the size of the array
 * @param x value of all component
 * @return  a TYPE(PnlVect)pointer all value at 0
 */
TYPE(PnlVect)* FUNCTION(pnl_vect,CONCAT2(create_from_,BASE))(const int size,BASE x)
{
  int i;
  TYPE(PnlVect) *v;
  if ((v=FUNCTION(pnl_vect, create)(size))==NULL)
    return NULL;
  for ( i=0 ; i<size ; i++ )
    {
      FUNCTION(pnl_vect, set) (v, i, x);
    }
  return v;
}

/**
 * creates a new TYPE(PnlVect)pointer.
 *
 * @param size the size of the array
 * @return  a TYPE(PnlVect)pointer all value at 0
 */
TYPE(PnlVect)* FUNCTION(pnl_vect,create_from_zero)(const int size)
{
  const BASE zero = ZERO;
  return FUNCTION(pnl_vect,CONCAT2(create_from_,BASE)) (size, zero);
}

/**
 * creates a new TYPE(PnlVect)pointer.
 *
 * @param size the size of the array
 * @param x the values used to fill the array. Obviously
 * \c x has to be at least of length size.
 * @return  a TYPE(PnlVect)pointer
 */
TYPE(PnlVect) * FUNCTION(pnl_vect,create_from_ptr)(const int size,const BASE  * x)
{
  TYPE(PnlVect) * v;
  if ((v = FUNCTION(pnl_vect,create) (size)) == NULL)
    return NULL;
  memcpy (v->array, x, v->size * sizeof(BASE));
  return v;
}

/**
 * creates a new TYPE(PnlVect)pointer.
 *
 * @param size the size of the array
 * @param ... is a list of values o type BASE used to fill the array. The list
 * must be of length size
 * @return  a TYPE(PnlVect)pointer
 */
TYPE(PnlVect) * FUNCTION(pnl_vect,create_from_list)(const int size,...)
{
  TYPE(PnlVect) * v;
  va_list ap;
  int i;

  if ((v = FUNCTION(pnl_vect,create)(size)) == NULL)
    return NULL;
  va_start (ap, size);

  for ( i=0; i<size ; i++ )
    {
      BASE val ;
      val = va_arg (ap, BASE);
      FUNCTION(pnl_vect,set)(v, i, val);
    }
  va_end(ap);
  return v;
}



/**
 * creates a TYPE(PnlVect)
 * @param m number of elements
 * @param x an array of BASE used to fill the TYPE(PnlVect). should be of length
 * m. No test can be done about it.
 * @return a TYPE(PnlVect) pointer
 */
TYPE(PnlVect) FUNCTION(pnl_vect,create_wrap_array)(const BASE* x, int m)
{
  TYPE(PnlVect) v;

  v.size = m;
  v.owner = 0;
  v.mem_size = 0;
  v.array = (BASE *) x;
  return v;
}


/**
 * Reads a vector from a file and creates the corresponding PnlVect
 * @param file  the file to be read
 * @return a PnlVect
 */
TYPE(PnlVect)* FUNCTION(pnl_vect,create_from_file) (const char * file)
{
  char car, empty=1;
  TYPE(PnlVect) *v;
  int m, count;
  ATOMIC *atomic_data;
  FILE *FIC = fopen( file, "r");
  if ( FIC == NULL )
    {
      PNL_ERROR("Cannot open file", "pnl_vect_create_from_file");
    }

  /* first pass to determine dimensions */
  m = 0; empty = 1;
  while((car=fgetc(FIC))!= EOF)
    {
      if ( car=='\n' )
        {
          if (!empty) { ++m; empty = 1;}
          else break;
        }
      else if (empty && isdigit(car)) empty=0;
    }
  if (m==0)
    {
      PNL_ERROR ("No matrix found in input file",  "pnl_vect_create_from_file");
    }

  /* need special care when MULTIPLICITY > 1 */
  if ( MULTIPLICITY * ((double)m/MULTIPLICITY) != m)
    {
      PNL_ERROR ("Incorrect vector format",  "pnl_mat_create_from_file");
    }
  if ((v = FUNCTION(pnl_vect,create) (m/MULTIPLICITY)) == NULL)
    {
      PNL_ERROR("Allocation error", "pnl_vect_create_from_file");
    }  

  /* second pass to read data */
  rewind( FIC );
  count = 0;
  atomic_data = (ATOMIC*) FUNCTION(pnl_vect,lget) (v, 0);
  while( fscanf(FIC,IN_FORMAT,atomic_data) > 0 && count<m)   { atomic_data++; count++;}
  fclose( FIC );
  return v;
}

/**
 * free a TYPE(PnlVect)pointer and set the data pointer to
 * NULL
 *
 * @param v address of the pointer to free
 */
void FUNCTION(pnl_vect,free)(TYPE(PnlVect) **v)
{
  if (*v != NULL)
    {
      if ((*v)->array != NULL && (*v)->owner == 1) free((*v)->array);
      free(*v);
      *v=NULL;
    }
}

/**
 * resizes a TYPE(PnlVect). If the new size is smaller than the
 * current one, no memory is freed and the datas are
 * kept. If the new size is larger than the current one, a
 * new pointer is allocated. The old datas are lost. 
 *
 * @param v a pointer to an already existing TYPE(PnlVect). If v->owner=0,
 * nothing is done
 * @param size the new size of the array
 */
int FUNCTION(pnl_vect,resize)(TYPE(PnlVect) * v, int size)
{

  if (v->owner == 0) return OK;
  if (size < 0) return FAIL;
  if (size == 0)
    {
      if (v->mem_size > 0) free (v->array);
      v->size = 0;
      v->mem_size = 0;
      v->array=NULL;
      return OK;
    }
  
  if (v->mem_size >= size)
    {
      /* If the new size is smaller, we do not reduce the size of the
         allocated block. It may change, but il may allow to grow the vector
         quicker */
      v->size=size; return OK;
    }

  /* Now, v->mem_size < size */
  if (v->array != NULL) free (v->array);
  if ((v->array = malloc (size * sizeof(BASE))) == NULL) return FAIL;
  v->size = size;
  v->mem_size = size;
  return OK;
}


/**
 * resizes a PnlVect. If the new size is smaller than the
 * current one, no memory is freed. If the new size is larger
 * than the current mem_size, a new pointer is allocated, the old
 * data are copied and the old pointer is freed.
 *
 * @param v : a pointer to an already existing PnlVect
 * @param size : the new size of the array
 * @param x : the value used to fill the new cells when the
 * array is enlarged
 */
int FUNCTION(pnl_vect,CONCAT2(resize_from_,BASE))(TYPE(PnlVect) *v, int size, BASE x)
{
  int i=0;
  int old_size = v->size;

  if (v->size >0 && v->owner == 0) return FAIL;

  v->size=size;
  if (size <= old_size)
    return OK;
  if (v->array == NULL)
    {
      if ((v->array=malloc(size*sizeof(BASE)))==NULL)
        return FAIL;
    }
  else
    {
      if ((v->array=realloc(v->array, size*sizeof(BASE)))==NULL)
        return FAIL;
    }
  for(i=old_size; i<size; i++)
    {
      FUNCTION(pnl_vect, set)(v, i, x);
    }
  return OK;
}

/**
 * resizes a PnlVect. If the new size is smaller than the current one, no
 * memory is freed. If the new size is larger than the current mem_size, a new
 * pointer is allocated.
 *
 * @param v a pointer to an already existing PnlVect
 * @param size  the new size of the array
 * @param t array used to fill v, must be of size size
 */
int FUNCTION(pnl_vect,resize_from_ptr)(TYPE(PnlVect) *v, int size, const BASE *t)
{
  if (FUNCTION(pnl_vect,resize) (v, size) != OK) return FAIL;
  memcpy (v->array, t, size*sizeof(BASE));
  return OK;
}


/**
 * Put all the components of v to x, v[i]=x for all i
 *
 * @param v a  TYPE(PnlVect)
 * @param x a  Base element
 */
void FUNCTION(pnl_vect,CONCAT2(set_,BASE))(TYPE(PnlVect) * v, BASE x)
{
  int i;
  for(i = 0; i < v->size; i++)
    {v->array[i] = x; }
}
/**
 * Sets v to zero v[i]=0 for all i
 *
 * @param v a TYPE(PnlVect) 
 */
void FUNCTION(pnl_vect, set_zero)(TYPE(PnlVect) * v)
{
  const BASE zero = ZERO ;
  int i;
  for(i = 0; i < v->size; i++)
    {v->array[i] = zero;}
}
/**
 * Sets v to the i-th basis vector
 * v[j]=0 if i!=j and v[i] = 1
 *
 * @param v a constant TYPE(PnlVect)pointer
 * @param i a index of the vector
 */
void FUNCTION(pnl_vect, set_basis)(TYPE(PnlVect) * v, int i)
{
  const BASE zero = ZERO ;
  const BASE one = ONE;
  int k;
  CheckIndexVect(v,i);
  for(k = 0; k < v->size; k++)
    {v->array[k] = zero; }
  v->array[i] = one;
}

/**
 * copies a TYPE(PnlVect)
 *
 * @param v a constant TYPE(PnlVect) pointer
 * @return a TYPE(PnlVect)
 */
TYPE(PnlVect) * FUNCTION(pnl_vect, copy)(const TYPE(PnlVect) * v)
{
  TYPE(PnlVect) *ret;
  if ((ret=FUNCTION(pnl_vect,create)(v->size))==NULL) return NULL;
  memcpy(ret->array, v->array, sizeof(BASE)*ret->size);
  return ret;
}

/**
 * clones a TYPE(PnlVect)
 *
 * @param clone needs to be already allocated.
 * @param v a constant TYPE(PnlVect)pointer
 */
void FUNCTION(pnl_vect, clone)(TYPE(PnlVect) * clone,
                               const TYPE(PnlVect) * v)
{
  if (clone->owner == 0 && clone->size != v->size)
    {
      PNL_ERROR ("owner == 0 and size mismatch", "pnl_vect_clone");
    }
  FUNCTION(pnl_vect,resize)(clone, v->size);
  memcpy(clone->array, v->array, sizeof(BASE)*v->size);
}


/**
 * prints a TYPE(PnlVect)in file fic.
 *
 * @param V a(constant) TYPE(PnlVect) ptr.
 * @param fic a file descriptor.
 */
void FUNCTION(pnl_vect, fprint)(FILE *fic, const TYPE(PnlVect) * V)
{
  int i=0;
  while(i<V->size) 
    {
      fprintf(fic,OUT_FORMAT,OUT_PUT_FORMAT(V->array[i]));
      fprintf(fic, " ");
      i++;
    }
  fprintf(fic,"\n");
}

/**
 * prints a TYPE(PnlVect).
 *
 * @param V a(constant) TYPE(PnlVect)ptr.
 */
void FUNCTION(pnl_vect, print)(const TYPE(PnlVect) * V)
{
  FUNCTION(pnl_vect, fprint)(stdout, V);
}

/**
 * prints a TYPE(PnlVect) in a format compatible with Nsp
 *
 * @param V a(constant) TYPE(PnlVect)ptr.
 */
void FUNCTION(pnl_vect, print_nsp)(const TYPE(PnlVect) * V)
{
  int i;
  printf("[ ");
  for (i=0; i<V->size-1; i++) 
    {
      printf(OUT_FORMAT,OUT_PUT_FORMAT(V->array[i]));
      printf("; ");
    }
  printf(OUT_FORMAT,OUT_PUT_FORMAT(V->array[i]));
  printf(" ]; \n");
}

/**
 * sets the value of v[i]=x
 *
 * @param v a TYPE(PnlVect)
 * @param i index
 * @param x v[i]=x
 */
void
FUNCTION(pnl_vect, set)(TYPE(PnlVect) * v,
                        const int i,    const BASE x)
{
#ifndef PNL_RANGE_CHECK_OFF
  if (i>=v->size || i<0)
    {
      PNL_ERROR ("index out of range", "pnl_vect_set");
    }
#endif
  v->array[i]=x;
}

/**
 * gets the value of v[i]. Can NOT be used as a lvalue
 *
 * @param v a TYPE(PnlVect)
 * @param i index
 * @return  v[i]
 */
BASE FUNCTION(pnl_vect, get)(const TYPE(PnlVect) * v, 
                             const int i)
{
  CheckIndexVect(v,i);
  return v->array[i];
}


/**
 * returns the address of v[i]. Can be used as a lvalue.
 *
 * @param v a TYPE(PnlVect)
 * @param i index
 * @return  &(v[i])
 */
BASE * FUNCTION(pnl_vect, lget)(TYPE(PnlVect) * v ,const int i)
{
#ifndef PNL_RANGE_CHECK_OFF
  if (i>=v->size || i<0)
    {
      PNL_ERROR ("index out of range", "pnl_vect_lget");
    }
#endif
  return &(v->array[i]);
}

/**
 * in-place map function
 *
 * @param lhs left hand side vector
 * @param f the function to be applied term by term
 * @return  lhs = f(lhs)
 */
void FUNCTION(pnl_vect,map_inplace)(TYPE(PnlVect) *lhs, 
                                     BASE(*f)(BASE ))
{
  BASE *lptr;
  int i=0;
  while(i<lhs->size)
    {
      lptr = FUNCTION(pnl_vect, lget)(lhs, i);
      (*lptr)=(*f)(*lptr); i++;
    }
}

#if (defined BASE_DOUBLE || defined BASE_PNL_COMPLEX || defined BASE_INT)
static BASE FUNCTION(,_unary_minus)(BASE a) { return MINUS(ZERO,a); }

/**
 * in-place unary minus vector operator application
 *
 * @param lhs left hand side vector. On exit lhs = -lhs
 */
void FUNCTION(pnl_vect,minus)(TYPE(PnlVect) *lhs)
{
  FUNCTION(pnl_vect,map_inplace)(lhs, FUNCTION(,_unary_minus));
}
#endif

/**
 * in-place vector operator application
 *
 * @param lhs left hand side vector
 * @param x double arg
 * @param op a binary operator, given as a function
 * @return  lhs = lhs op x
 */
static void FUNCTION(pnl_vect,apply_op)(TYPE(PnlVect) *lhs,
                                        BASE  x, 
                                        BASE (*op)(BASE, BASE ))
{
  BASE *lptr;
  int i=0;
  while(i<lhs->size)
    {
      lptr = FUNCTION(pnl_vect, lget)(lhs, i);
      (*lptr)= op(*lptr,x); i++;
    }
}


static BASE FUNCTION(,_op_plus)(BASE a, BASE b) { return PLUS(a,b);}
static BASE FUNCTION(,_op_minus)(BASE a, BASE b) { return MINUS(a,b);}
static BASE FUNCTION(,_op_mult)(BASE a, BASE b) { return MULT(a,b); }
static BASE FUNCTION(,_op_div)(BASE a, BASE b) { return DIV(a,b); }
static BASE FUNCTION(,_op_inv)(BASE a) { return INV(a); }
static double FUNCTION(,_op_sqr_norm)(BASE a) { return SQUARE_NORM(a); }
static double FUNCTION(,_op_abs)(BASE a) { return NORMONE(a); }



/**
 * in-place vector scalar addition
 *
 * @param lhs left hand side vector
 * @param x scalar
 * @return  lhs = lhs+x
 */
void FUNCTION(pnl_vect,CONCAT2(plus_,BASE))(TYPE(PnlVect) *lhs , BASE x)
{
  FUNCTION(pnl_vect,apply_op)(lhs, x,FUNCTION(,_op_plus));
}

#if (defined BASE_DOUBLE || defined BASE_PNL_COMPLEX || defined BASE_INT)
/**
 * in-place vector scalar substraction
 *
 * @param lhs left hand side vector
 * @param x scalar
 * @return  lhs = lhs-x
 */
void FUNCTION(pnl_vect,CONCAT2(minus_,BASE))(TYPE(PnlVect) *lhs , BASE x)
{
  FUNCTION(pnl_vect,apply_op)(lhs, x, FUNCTION(,_op_minus));
}
#endif

/**
 * in-place vector scalar multiplication
 *
 * @param lhs left hand side vector
 * @param x scalar
 * @return  lhs = lhs*x
 */
void FUNCTION(pnl_vect,CONCAT2(mult_,BASE))(TYPE(PnlVect) *lhs , BASE x)
{
  FUNCTION(pnl_vect,apply_op)(lhs, x, FUNCTION(,_op_mult));
}

/**
 * in-place vector scalar division
 *
 * @param lhs left hand side vector
 * @param x scalar
 * @return  lhs = lhs/x
 */
void FUNCTION(pnl_vect,CONCAT2(div_,BASE))(TYPE(PnlVect) *lhs , BASE x)
{
  FUNCTION(pnl_vect,apply_op)(lhs, x, FUNCTION(,_op_div));
}


/**
 * map vector componentwise
 *
 * @param lhs each component lhs(i) contains f(rhs(i))
 * @param rhs right hand side vector
 * @param f real function 
 */
void FUNCTION(pnl_vect,map)(TYPE(PnlVect) *lhs, const TYPE(PnlVect) *rhs, BASE(*f)(BASE))
{
  FUNCTION(pnl_vect,clone)(lhs, rhs);
  FUNCTION(pnl_vect,map_inplace)(lhs, f);
}

/**
 * map vector componentwise
 *
 * @param lhs each component lhs(i) contains f(lhs(i),rhs(i))
 * @param rhs right hand side vector
 * @param f real function 
 */
static void FUNCTION(pnl_vect,map_vect)(TYPE(PnlVect) *lhs, const TYPE(PnlVect) *rhs, BASE(*f)(BASE,BASE))
{
  BASE *lptr;
  int i;
  CheckVectMatch(lhs, rhs);
  i=0;
  while(i<lhs->size)
    {
      lptr = FUNCTION(pnl_vect,lget)(lhs, i);
      *lptr =(*f)(*lptr, FUNCTION(pnl_vect,get)(rhs,i)); i++;
    }
}


/**
 * in-place vector vector addition
 *
 * @param lhs left hand side vector
 * @param rhs rigth hand side vector
 * @return  lhs = lhs+rhs
 */
void FUNCTION(pnl_vect,plus_vect)(TYPE(PnlVect) *lhs, const TYPE(PnlVect) *rhs)
{
  FUNCTION(pnl_vect,map_vect)(lhs, rhs, FUNCTION(,_op_plus));
}

/**
 * in-place vector vector substraction
 *
 * @param lhs left hand side vector
 * @param rhs rigth hand side vector
 * @return  lhs = lhs+rhs
 */
void FUNCTION(pnl_vect,minus_vect)(TYPE(PnlVect) *lhs, const TYPE(PnlVect) *rhs)
{
  FUNCTION(pnl_vect,map_vect)(lhs, rhs, FUNCTION(,_op_minus));
}


/**
 * in-place term by term vector inverse
 *
 * @param lhs left hand side vector
 * @return  lhs = 1 ./ lhs
 */
void FUNCTION(pnl_vect,inv_term)(TYPE(PnlVect) *lhs)
{
  FUNCTION(pnl_vect,map_inplace)(lhs,FUNCTION(,_op_inv));
}

/**
 * in-place term by term vector inverse
 *
 * @param lhs left hand side vector
 * @param rhs right hand side vector
 * @return  lhs = lhs ./ rhs
 */
void FUNCTION(pnl_vect,div_vect_term)(TYPE(PnlVect) *lhs, const TYPE(PnlVect) *rhs)
{
  FUNCTION(pnl_vect,map_vect)(lhs, rhs,FUNCTION(,_op_div));
}

/**
 * in-place vector term by term multiplication
 *
 * @param lhs left hand side vector
 * @param rhs right hand side vector
 * @return  lhs = lhs.*rhs
 */
void FUNCTION(pnl_vect,mult_vect_term)(TYPE(PnlVect) *lhs, const TYPE(PnlVect) *rhs)
{
  FUNCTION(pnl_vect,map_vect)(lhs, rhs, FUNCTION(,_op_mult));
}


/**
 * sum vector componentwise
 *
 * @param lhs left hand side vector
 * @return  sum=sum(lhs)
 */
BASE FUNCTION(pnl_vect,sum)(const TYPE(PnlVect) *lhs)
{
  BASE sum;
  int i=0;
  sum=ZERO;
  while(i<lhs->size)
    {sum = PLUS(sum, FUNCTION(pnl_vect,get)(lhs, i)); i++;}
  return sum;
}

/**
 * cumsum vector componentwise
 *
 * @param lhs : left hand side vector
 */
void FUNCTION(pnl_vect,cumsum)(TYPE(PnlVect) *lhs)
{
  BASE sum;
  int i;
  sum=ZERO;
  for (i=0; i<lhs->size; i++)
    {
      sum = PLUS(sum, FUNCTION(pnl_vect,get)(lhs, i)); 
      FUNCTION(pnl_vect, set) (lhs, i, sum);
    }
}


/**
 * computes the scalar product x' * y
 *
 * @param x : first right hand side  TYPE(PnlVect)
 * @param y : second right hand side TYPE(PnlVect)
 * @return the scalar product x' * y
 * and rhs2
 */
BASE FUNCTION(pnl_vect,scalar_prod)(const TYPE(PnlVect) *x, 
                                    const TYPE(PnlVect) *y)
{
  BASE xi, yi;
  BASE sum=ZERO;
  int i=0;
  CheckVectMatch(x, y);
  while(i<x->size)
    {
      xi = FUNCTION(pnl_vect,get)(x, i); 
      yi = FUNCTION(pnl_vect,get)(y, i);
#if MULTIPLICITY == 1
      sum += xi * yi;
#else
      sum = PLUS(sum,(MULT(xi, yi)));
#endif
      i++;
    }
  return sum;
}

/**
 * computes the product of the components of a vector
 *
 * @param V a a TYPE(PnlVect) 
 * @return the product of the components of V
 */
BASE FUNCTION(pnl_vect,prod)(const TYPE(PnlVect) *V)
{
  int i;
  BASE p=ONE;
  for (i=0; i<V->size; i++)
    {
      p = MULT(p,FUNCTION(pnl_vect,get)(V, i));
    }
  return p;
}

/**
 * computes the cumulative product of the components of a vector
 *
 * @param V a a TYPE(PnlVect) 
 */
void FUNCTION(pnl_vect,cumprod)(TYPE(PnlVect) *V)
{
  int i;
  BASE p=ONE;
  for (i=0; i<V->size; i++)
    {
      p = MULT(p,FUNCTION(pnl_vect,get)(V, i));
      FUNCTION(pnl_vect, set) (V, i, p);      
    }
}





/* Swap  */

/**
 * Swaps two elements of a vector
 *
 * @param v a TYPE(PnlVect)
 * @param i first index
 * @param j second index
 */
void FUNCTION(pnl_vect, swap_elements)(TYPE(PnlVect) * v, int i, int j)
{
#ifndef PNL_RANGE_CHECK_OFF
  if (i>=v->size || j>=v->size)
    {
      PNL_ERROR ("index out of range", "pnl_vect_swap_elements");
    }
#endif
  if (i != j)
    {
      BASE vj = v->array[j];
      BASE *vi = &(v->array[i]);
      v->array[j] = *vi;
      *vi = vj;
    }
}

/**
 * Mirror a vector
 *
 * @param v a TYPE(PnlVect)
 */
void FUNCTION(pnl_vect, reverse)(TYPE(PnlVect) * v)
{
  int i ;
  int size = v->size ;
  for(i = 0 ; i <(size / 2) ; i++)
    {
      int j = size - i - 1 ;
      BASE tmp = v->array[j];
      v->array[j] = v->array[i];
      v->array[i] = tmp;
    }
}

#if defined ORDERED

/**
 * returns the index of the minimum of the elements of an array
 *
 * @param a a vector
 * @param n the size of the vector
 * @param incr the jump between two consecutive elements (useful to treat non
 * contiguous arrays)
 * @param imin_out is a pointer to BASE, contains the index of the minimum on
 * exit
 * @param min_out contains the the minimum on exit
 */
void FUNCTION(pnl_array, min_index)(const BASE *a, int n, int incr,
                                    BASE *min_out, int *imin_out)
{
  BASE min = a[0];
  int j, i, imin = 0;
  for(i = 1, j=incr; i < n; i++, j += incr)
    {
      BASE x = a[j];
      if (x < min)
        {
          min = x; imin = i;
        }
      if (isnan(x))
        {
          *imin_out = i; *min_out = x;
          return;
        }
    }
  *imin_out = imin;
  *min_out = min;
}

/**
 * returns the index of the maximum of the elements of an array
 *
 * @param a a vector
 * @param n the size of the vector
 * @param incr the jump between two consecutive elements (useful to treat non
 * contiguous arrays)
 * @param imax_out is a pointer to BASE, contains the index of the maximum on
 * exit
 * @param max_out contains the the maximum on exit
 */
void FUNCTION(pnl_array, max_index)(const BASE *a, int n, int incr,
                                    BASE *max_out, int *imax_out)
{
  BASE max = a[0];
  int j, i, imax = 0;
  for(i = 1, j=incr; i < n; i++, j += incr)
    {
      BASE x = a[j];
      if (x > max)
        {
          max = x; imax = i;
        }
      if (isnan(x))
        {
          *imax_out = i; *max_out = x;
          return;
        }
    }
  *imax_out = imax;
  *max_out = max;
}

/**
 * returns the indices of the minimum and maximum of the elements of an array
 *
 * @param a a vector
 * @param n the size of the vector
 * @param incr the jump between two consecutive elements (useful to treat non
 * contiguous arrays)
 * @param imin_out is a pointer to BASE, contains the index of the minimum on exit
 * @param imax_out is a pointer to BASE, contains the index of the maximum on
 * exit
 * @param min_out contains the minimum on exit
 * @param max_out contains the maximum on exit
 */
void FUNCTION(pnl_array, minmax_index)(const BASE *a, int n, int incr,
                                       BASE *min_out, BASE *max_out,
                                       int *imin_out, int *imax_out)
{
  int i, j;
  int imin = 0, imax = 0;
  BASE max = a[0];
  BASE min = max;

  for(i = 1, j=incr; i < n; i++, j += incr)
    {
      BASE x = a[j];
      if (x < min)
        {
          min = x; imin = i;
        }
      if (x > max)
        {
          max = x; imax = i;
        }
      if (isnan(x))
        {
          *imin_out = *imax_out = i;
          *min_out = *max_out = x;
          return;
        }
    }
  *imin_out = imin;
  *imax_out = imax;
  *min_out = min;
  *max_out = max;
}

/**
 * returns the maximum of the components of a vector
 *
 * @param V a vector 
 * @return the maximum of the components of V
 */
BASE FUNCTION(pnl_vect,max)(const TYPE(PnlVect) *V)
{
  BASE m;
  int i;

  FUNCTION (pnl_vect, max_index)(V, &m, &i);
  return m;
}

/**
 * returns the minimum of the components of a vector
 *
 * @param V a vector 
 * @return the minimum of the components of V
 */
BASE FUNCTION(pnl_vect,min)(const TYPE(PnlVect) *V)
{
  BASE m;
  int i;

  FUNCTION (pnl_vect, min_index)(V, &m, &i);
  return m;
}



/**
 * returns the minimum and maximum of the elements of a vector
 *
 * @param V a vector
 * @param min_out is a pointer to BASE, contains the minimum on exit
 * @param max_out is a pointer to BASE, contains the maximum on exit
 */
void FUNCTION(pnl_vect, minmax)(const TYPE(PnlVect) *V, BASE *min_out, BASE *max_out)
{
  int im, iM;
  FUNCTION(pnl_vect, minmax_index)(V, min_out, max_out, &im, &iM);
}


/**
 * returns the index of the maximum of the elements of a vector
 *
 * @param m a vector
 * @param imax_out is a pointer to int, contains the index of the maximum on exit
 * @param max_out contains the the maximum on exit
 */
void FUNCTION(pnl_vect, max_index)(const TYPE(PnlVect) *m, BASE *max_out, int *imax_out)
{
  FUNCTION(pnl_array, max_index)(m->array, m->size, 1, max_out, imax_out);
}
  
/**
 * returns the index of the minimum of the elements of a vector
 *
 * @param m a vector
 * @param imin_out is a pointer to BASE, contains the index of the minimum on
 * exit
 * @param min_out contains the the minimum on exit
 */
void FUNCTION(pnl_vect, min_index)(const TYPE(PnlVect) *m, BASE *min_out, int *imin_out)
{
  FUNCTION(pnl_array, min_index)(m->array, m->size, 1, min_out, imin_out);
}

/**
 * returns the indices of the minimum and maximum of the elements of a vector
 *
 * @param m a vector
 * @param imin_out is a pointer to BASE, contains the index of the minimum on exit
 * @param imax_out is a pointer to BASE, contains the index of the maximum on
 * exit
 * @param min_out contains the minimum on exit
 * @param max_out contains the maximum on exit
 */
void FUNCTION(pnl_vect, minmax_index)(const TYPE(PnlVect) *m,
                                      BASE *min_out, BASE *max_out,
                                      int *imin_out, int *imax_out)
{
  FUNCTION(pnl_array, minmax_index)(m->array, m->size, 1, min_out, max_out,
                                    imin_out, imax_out);
}

static int FUNCTION(__pnl, cmp_i) ( const void *a, const void *b)
{
  if ( * (BASE *) a < * (BASE *) b) return -1;
  if ( * (BASE *) a == * (BASE *) b) return 0;
  return 1;
}

static int FUNCTION(__pnl, cmp_d) ( const void *a, const void *b)
{
  if ( * (BASE *) a > * (BASE *) b) return -1;
  if ( * (BASE *) a == * (BASE *) b) return 0;
  return 1;
}


static void FUNCTION(pnl_vect, qsort_aux)(TYPE(PnlVect) * m, PnlVectInt *t, int use_index, char order)
{
  int i, *index = NULL;
  if (use_index == TRUE)
    {
      pnl_vect_int_resize (t, m->size);
      for (i=0; i<t->size; i++) t->array[i] = i;
      index = t->array;
    }

  switch (order)
    {
    case 'i' :
      pnl_qsort ( m->array, m->size, sizeof(BASE), 1, index, 1, use_index, FUNCTION(__pnl, cmp_i));
      break;
    case 'd' :
      pnl_qsort ( m->array, m->size, sizeof(BASE), 1, index, 1, use_index, FUNCTION(__pnl, cmp_d));
      break;
    default:
      PNL_ERROR ("unknow order", "pnl_vect_qsort"); break;
    }
}

/**
 * Quick sort function for vectors
 *
 * @param m a TYPE(PnlVect), contains the sorted vector on exit
 * @param order can be 'i' or 'd' for increasing or decreasing order.
 */
void FUNCTION(pnl_vect, qsort)(TYPE(PnlVect) * m, char order)
{
  FUNCTION(pnl_vect, qsort_aux)(m, NULL, FALSE, order);
}

/**
 * Quick sort function for vectors with permutation index
 *
 * @param m a TYPE(PnlVect), contains the sorted vector on exit
 * @param order can be 'i' or 'd' for increasing or decreasing order.
 * @param t on exit contains the permutation used to sort  the vector
 */
void FUNCTION(pnl_vect, qsort_index)(TYPE(PnlVect) * m, PnlVectInt *t, char order)
{
  FUNCTION(pnl_vect, qsort_aux)(m, t, TRUE, order);
}

#endif

/*************************************************************************/
/* Written and (C) by David Pommier <pommier.david@gmail.com>            */ 
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

/**
 * Computes a x + b y and stores the result in y
 *
 * @param a BASE
 * @param x a vector
 * @param b BASE 
 * @param y a vector
 */
void FUNCTION(pnl_vect,axpby)(BASE a, const TYPE(PnlVect) *x,  BASE b, TYPE(PnlVect) *y)
{
  BASE zero, one, xi, *yi;
  int i;
  zero = ZERO;
  one = ONE;

  if ( EQ(b,zero) )
    {
      FUNCTION(pnl_vect,resize) (y, x->size);
      FUNCTION(pnl_vect,CONCAT2(set_,BASE))(y, zero);
    }
  else if ( NEQ(b,one) ) FUNCTION(pnl_vect,CONCAT2(mult_,BASE)) (y, b);
  if ( EQ(a,zero) ) return;

#ifndef PNL_RANGE_CHECK_OFF
  if (x->size != y->size)
    {
      PNL_ERROR ("size mismatch", "pnl_vect_axpby");
    }
#endif

  for (i=0; i< x->size; i++)
    {
      xi =  FUNCTION(pnl_vect,get)(x, i);
      xi = MULT(xi, a);
      yi = FUNCTION(pnl_vect,lget)(y, i);
      *yi = PLUS(*yi, xi);      
    }
}

/**
 * returns the X norm of a vector 
 *
 * @param V a TYPE(PnlVect)
 * @param f a pointer's function on a norm function
 * @return the X norm of V
 */
double FUNCTION(pnl_vect,norm_x)(const TYPE(PnlVect) *V,double(*f)(BASE))
{
  int i=0;
  double p=0;
  while(i<V->size)
    { p+=(*f)(FUNCTION(pnl_vect,get)(V, i)); i++;}
  return p;
}
/**
 * returns the two norm of a vector 
 *
 * @param V  TYPE(PnlVect) 
 * @return  the square root of the sum of the square value of components of V
 */
double FUNCTION(pnl_vect,norm_two)(const TYPE(PnlVect) *V)
{
  return sqrt(FUNCTION(pnl_vect,norm_x(V,FUNCTION(,_op_sqr_norm))));
}

/**
 * returns the one norm of a vector 
 *
 * @param V a vector 
 * @return the sum of the absolute value of components of V
 */
double FUNCTION(pnl_vect,norm_one)(const TYPE(PnlVect) *V) /*res=\Vert rhs1 \Vert_{l^1} */
{
  return FUNCTION(pnl_vect,norm_x)(V,FUNCTION(,_op_abs));
}
/**
 * returns the infty norm of a vector 
 *
 * @param V a TYPE(PnlVect) 
 * @return the maximum of the absolute value of components of V
 */
double FUNCTION(pnl_vect,norm_infty)(const TYPE(PnlVect) *V)
{
  int i=0;
  double p=0.0;
  double q=0.0;
  while(i<V->size)
    {
      q=NORMONE(FUNCTION(pnl_vect,get)(V,i));
      p=((p<q)?q:p);i++;
    }
  return p;
}

/**
 * Extract a sub vector and wrap it into a vector.
 * @param V a vector
 * @param i the index of first element to be extracted
 * @param s the size of extracted vector
 * @return a vector (not a pointer) whose array pointer is the address of the
 * ith element of V. No copying is done.
 */
TYPE(PnlVect) FUNCTION(pnl_vect, wrap_subvect)(const TYPE(PnlVect) *V, int i,int s)
{ 
  TYPE(PnlVect) ret;
#ifndef PNL_RANGE_CHECK_OFF
  if ((i >= V->size) ||(s >= V->size) || (s<=0))  {PNL_ERROR ("index out of range", "pnl_vect_extract_with_size");}
#endif
  ret.size = s;
  ret.mem_size = 0;
  ret.owner = 0;
  ret.array = &(V->array[i]); 
  return ret;
} 

/**
 * Extract  a sub vector and wrap it into a vector.
 * @param V a vector
 * @param i the index of first element to be extracted
 * @param j the index of last  element to be extracted
 * @return a vector (not a pointer) whose array pointer is the address of the
 * ith element of V. No copying is done.
 */
TYPE(PnlVect) FUNCTION(pnl_vect, wrap_subvect_with_last)(const TYPE(PnlVect) *V, int i,int j)
{
  return FUNCTION(pnl_vect, wrap_subvect)(V,i,j-i+1);
}


