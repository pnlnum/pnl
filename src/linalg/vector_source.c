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

/* To enable clang completion */
#ifdef PNL_CLANG_COMPLETE
#include <stdlib.h>
#include <string.h>

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_internals.h"
#define BASE_DOUBLE
#include "pnl/pnl_templates_on.h"
#endif



/**
 * Initialize a PnlVect with size 0
 */
void FUNCTION(pnl_vect,init)(TYPE(PnlVect) *o)
{
  o->object.parent_type = PNL_TYPE_VECTOR;
  o->object.type = CONCAT2(PNL_TYPE_VECTOR_, BASE_TYPE);
  o->object.label = FUNCTION(pnl_vect,label);
  o->object.destroy = (DestroyFunc *) pnl_vect_object_free;
  o->object.constructor = (NewFunc *) FUNCTION(pnl_vect,new);
  o->object.clone = (CloneFunc *) FUNCTION(pnl_vect,clone);
  o->object.copy = (CopyFunc *) FUNCTION(pnl_vect,copy);
  o->size = 0;
  o->mem_size = 0;
  o->owner = 1;
  o->array = NULL;
}

/**
 * Create a new PnlVect with size 0 and initializes the object
 *
 */
TYPE(PnlVect)* FUNCTION(pnl_vect,new)()
{
  TYPE(PnlVect) *o;
  if ( (o = (TYPE(PnlVect) *) pnl_vect_object_new ()) == NULL) return NULL;
  FUNCTION(pnl_vect,init)(o);
  return o;
}

/**
 * Create a new TYPE(PnlVect) pointer.
 *
 * @param size the size of the array
 * @return  a TYPE(PnlVect) pointer
 */
TYPE(PnlVect) * FUNCTION(pnl_vect,create)(int size)
{
  TYPE(PnlVect)* v;
  if ((v = FUNCTION(pnl_vect,new)()) == NULL) return NULL;
  if (size > 0)
    {
      if ((v->array = malloc (size * sizeof (BASE))) == NULL) return NULL;
      v->size = size;
      v->mem_size = size;
      v->owner = 1;
    }
  return v;
}

/**
 * Create a new TYPE(PnlVect) pointer.
 *
 * @param size the size of the array
 * @param x value of all component
 * @return  a TYPE(PnlVect)pointer all value at 0
 */
TYPE(PnlVect)* FUNCTION(pnl_vect,create_from_scalar)(const int size,BASE x)
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
 * Test if 2 vectors are equal
 *
 * @param v1 a vector
 * @param v2 a vector
 * @return  TRUE or FALSE
 */
int FUNCTION(pnl_vect,eq)(const TYPE(PnlVect) *v1, const TYPE(PnlVect) *v2)
{
  int i;
  if (v1->size != v2->size) return FAIL;
  for ( i=0 ; i<v1->size ; i++ )
    {
      if ( NEQ(PNL_GET(v1, i), PNL_GET(v2, i)) ) return FALSE;
    }
  return TRUE;
}

/**
 * Test if 2 vectors are equal
 *
 * @param v a vector
 * @param x a BASE element
 * @return  TRUE or FALSE
 */
int FUNCTION(pnl_vect,eq_all)(const TYPE(PnlVect) *v, BASE x)
{
  int i;
  for ( i=0 ; i<v->size ; i++ )
    {
      if ( NEQ(PNL_GET(v, i), x) ) return FALSE;
    }
  return TRUE;
}

/**
 * Create a new TYPE(PnlVect) pointer.
 *
 * @param size the size of the array
 * @return  a TYPE(PnlVect) pointer with all entries set to 0
 */
TYPE(PnlVect)* FUNCTION(pnl_vect,create_from_zero)(const int size)
{
  TYPE(PnlVect)* v = FUNCTION(pnl_vect, create)(size);
  FUNCTION(pnl_vect, set_zero)(v);
  return v;
}

/**
 * Create a new TYPE(PnlVect) pointer.
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
 * Create a new TYPE(PnlVect) pointer.
 *
 * @param M a TYPE(PnlMat)
 * @return  a TYPE(PnlVect)pointer
 */
TYPE(PnlVect) * FUNCTION(pnl_vect,create_from_mat)(const TYPE(PnlMat) *M)
{
  TYPE(PnlVect) * v;
  if ((v = FUNCTION(pnl_vect,create) (M->mn)) == NULL)
    return NULL;
  memcpy (v->array, M->array, v->size * sizeof(BASE));
  return v;
}

/**
 * Create a new TYPE(PnlVect)pointer.
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
      PNL_SET(v, i, val);
    }
  va_end(ap);
  return v;
}

/**
 * Read a vector from a file and creates the corresponding PnlVect. The
 * data can be stored column or row wise
 * @param file  the file to be read
 * @return a PnlVect
 */
TYPE(PnlVect)* FUNCTION(pnl_vect,create_from_file) (const char * file)
{
  char car, prev='\0', empty=1;
  TYPE(PnlVect) *v;
  int m, count;
  BASE *data;
  FILE *FIC = fopen( file, "r");
  if ( FIC == NULL )
    {
      PNL_ERROR("Cannot open file", "pnl_vect_create_from_file");
    }

  /*
   * First pass to determine dimensions : need to find out if the data are
   * stored row-wise or column-wise
   */
  /*
   * Count the number of elements on the first line
   */
  m = 1;
  while((car=fgetc(FIC))!='\n' && car != EOF)
    {
      if (isdigit(car) || car == '-' || car == '.')
        {
          empty = 0;
          if (prev == ' ' || prev == '\t' ) ++m;
        }
      prev = car;
    }
  if ( empty ) { m = 0; }
  m = m / MULTIPLICITY;
  if ( m == 1 )
    { 
      /* We consider it as a column vector */
      empty = 1;
      while((car=fgetc(FIC))!= EOF)
        {
          if ( car=='\n' )
            {
              if (!empty) { ++m; empty = 1;}
              else break;
            }
          else if (empty && isdigit(car)) empty=0;
        }
    }
  if (m==0)
    {
      v = FUNCTION(pnl_vect,create) (0);
      fclose (FIC);
      return v;
    }

  if ((v = FUNCTION(pnl_vect,create) (m)) == NULL)
    {
      PNL_ERROR("Allocation error", "pnl_vect_create_from_file");
    }

  /* second pass to read data */
  rewind (FIC);
  count = 0;
  data = FUNCTION(pnl_vect,lget) (v, 0);
  while ( fscanf(FIC,IN_FORMAT,IN_PUT_FORMAT(data)) == MULTIPLICITY && count<m )   { data++; count++;}
  fclose (FIC);
  return v;
}

/**
 * Extract a sub vector from a vector V. The components to be extracted are
 * sepcified by the indices listed in the vector ind
 *
 * @param V_sub (output) V(ind(:))
 * @param V a vector
 * @param ind a vector of integers representing the list of indices to be extracted
 */
void FUNCTION(pnl_vect,extract_subvect_with_ind)(TYPE(PnlVect) *V_sub, const TYPE(PnlVect) *V, const PnlVectInt *ind)
{
  int i;
  FUNCTION(pnl_vect,resize)(V_sub, ind->size);
  for ( i=0 ; i<ind->size ; i++ )
    {
      PNL_CHECK (PNL_GET(ind, i) >= V->size, "index exceeded", "pnl_vect_extract_subvect_with_ind");
      PNL_LET (V_sub, i) = PNL_GET(V, PNL_GET(ind, i));
    }
}

/**
 * Extract a sub vector from a vector V. The components to be extracted are
 * (i:i+len)
 *
 * @param V_sub (output) V(i:i+len-1)
 * @param V a vector
 * @param i starting point of extraction
 * @param len length of extraction
 */
void FUNCTION(pnl_vect,extract_subvect)(TYPE(PnlVect) *V_sub, const TYPE(PnlVect) *V, int i, int len)
{
  int j;
  FUNCTION(pnl_vect,resize)(V_sub, len);
  PNL_CHECK (V->size < len + i, "index exceeded", "pnl_vect_extract_subvect");
  for ( j=0 ; j<len ; j++ )
    {
      PNL_LET (V_sub, j) = PNL_GET(V, j+i);
    }
}

/**
 * Extract a sub vector from a vector V. The components to be extracted are
 * sepcified by the indices listed in the vector ind
 *
 * @param V a vector
 * @param ind a vector of integers representing the list of indices to be extracted
 * @return V(ind(:))
 */
TYPE(PnlVect)* FUNCTION(pnl_vect,create_subvect_with_ind) (const TYPE(PnlVect) *V, const PnlVectInt *ind)
{
  TYPE(PnlVect) *V_sub;
  V_sub  = FUNCTION(pnl_vect,create)(ind->size);
  FUNCTION(pnl_vect,extract_subvect_with_ind)(V_sub, V, ind);
  return V_sub;
}

/**
 * Extract a sub vector from a vector V. The components to be extracted are
 * (i:i+len)
 *
 * @param V a vector
 * @param i starting point of extraction
 * @param len length of extraction
 * @return V(i:i+len-1)
 */
TYPE(PnlVect)* FUNCTION(pnl_vect,create_subvect) (const TYPE(PnlVect) *V, int i, int len)
{
  TYPE(PnlVect) *V_sub;
  V_sub  = FUNCTION(pnl_vect,new)();
  FUNCTION(pnl_vect,extract_subvect)(V_sub, V, i, len);
  return V_sub;
}

/**
 * Create a TYPE(PnlVect)
 * @param m number of elements
 * @param x an array of BASE used to fill the TYPE(PnlVect). should be of length
 * m. No test can be done about it.
 * @return a TYPE(PnlVect) pointer
 */
TYPE(PnlVect) FUNCTION(pnl_vect,wrap_array)(const BASE* x, int m)
{
  TYPE(PnlVect) v;
  FUNCTION (pnl_vect,init) (&v);

  v.size = m;
  v.owner = 0;
  v.mem_size = 0;
  v.array = (BASE *) x;
  return v;
}

/**
 * Create a sub matrix from a matrix M only containing the indices present in
 * (indi ,indj)
 *
 * @param val (output) a vector containing M( (indi(:),indj(:)) )
 * @param M a matrix
 * @param indi a vector of integers
 * @param indj a vector of integers
 */
void FUNCTION(pnl_vect,extract_submat) (TYPE(PnlVect) *val, const TYPE(PnlMat) *M,
                                        const PnlVectInt *indi, const PnlVectInt *indj)
{
  int k;
  int i, j;
  CheckVectMatch (indi, indj);

  FUNCTION(pnl_vect,resize)(val, indi->size);
  for ( k=0 ; k<indi->size ; k++ )
    {
      i = PNL_GET (indi, k);
      j = PNL_GET (indj, k);
      PNL_CHECK (i<0 || j<0 || i>= M->m || j>= M->n, "index of range", "extract_submat");
      PNL_LET (val, k) = PNL_MGET (M, i, j);
    }
}

/**
 * Create a sub matrix from a matrix M only containing the indices present in
 * (indi ,indj)
 *
 * @param M a vector
 * @param indi a vector of integers
 * @param indj a vector of integers. The pairs (indi, indj) are the indices (i,
 * j) of the elements to be extracted
 * @return M(indi(:),indj(:))
 */
TYPE(PnlVect)* FUNCTION(pnl_vect,create_submat) (const TYPE(PnlMat) *M,
                                        const PnlVectInt *indi, const PnlVectInt *indj)
{
  TYPE(PnlVect) *val;
  val  = FUNCTION(pnl_vect,create)(0);
  FUNCTION(pnl_vect,extract_submat)(val, M, indi, indj);
  return val;
}


/**
 * free a TYPE(PnlVect)pointer and set the data pointer to
 * NULL
 *
 * @param v address of the pointer to free
 */
void FUNCTION(pnl_vect,free)(TYPE(PnlVect) **v)
{
  PnlVectObject *o;
  o = PNL_VECT_OBJECT(*v);
  pnl_vect_object_free (&o);
  *v = NULL;
}

/**
 * Resize a TYPE(PnlVect). If the new size is smaller than the
 * current one, no memory is freed and the datas are
 * kept. If the new size is larger than the current one, a
 * new pointer is allocated. The old datas are kept.
 *
 * @param v a pointer to an already existing TYPE(PnlVect). If v->owner=0,
 * nothing is done
 * @param size the new size of the array
 */
int FUNCTION(pnl_vect,resize)(TYPE(PnlVect) * v, int size)
{
  return pnl_vect_object_resize (PNL_VECT_OBJECT(v), size);
}

/**
 * Resize a PnlVect. If the new size is smaller than the
 * current one, no memory is freed. If the new size is larger
 * than the current mem_size, a new pointer is allocated, the old
 * data are copied and the old pointer is freed.
 *
 * @param v : a pointer to an already existing PnlVect
 * @param size : the new size of the array
 * @param x : the value used to fill the new cells when the
 * array is enlarged
 */
int FUNCTION(pnl_vect,resize_from_scalar)(TYPE(PnlVect) *v, int size, BASE x)
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
 * Resize a PnlVect. If the new size is smaller than the current one, no
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
void FUNCTION(pnl_vect,set_all)(TYPE(PnlVect) * v, BASE x)
{
  int i;
  for(i = 0; i < v->size; i++)
    {v->array[i] = x; }
}
/**
 * Set v to zero v[i]=0 for all i
 *
 * @param v a TYPE(PnlVect)
 */
void FUNCTION(pnl_vect, set_zero)(TYPE(PnlVect) * v)
{
  memset (v->array, 0, v->size * sizeof(BASE));
}

/**
 * Copie a TYPE(PnlVect)
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
 * Clone a TYPE(PnlVect)
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
 * Print a TYPE(PnlVect) to a  file with one entry per line.
 *
 * @param V a(constant) TYPE(PnlVect)
 * @param fic a file descriptor.
 */
void FUNCTION(pnl_vect, fprint)(FILE *fic, const TYPE(PnlVect) * V)
{
  int i;
  for ( i=0 ; i<V->size ; i++ )
    {
      fprintf(fic,OUT_FORMAT,OUT_PUT_FORMAT(V->array[i]));
      fprintf(fic, "\n");
    }
}

/**
 * Print a TYPE(PnlVect) with one entry per line
 *
 * @param V a(constant) TYPE(PnlVect).
 */
void FUNCTION(pnl_vect, print)(const TYPE(PnlVect) * V)
{
  FUNCTION(pnl_vect, fprint)(stdout, V);
}

/**
 * Print a TYPE(PnlVect) to a  file with as a row vector
 *
 * @param V a(constant) TYPE(PnlVect)
 * @param fic a file descriptor.
 */
void FUNCTION(pnl_vect, fprint_asrow)(FILE *fic, const TYPE(PnlVect) * V)
{
  int i;
  for ( i=0 ; i<V->size ; i++ )
    {
      fprintf(fic,OUT_FORMAT,OUT_PUT_FORMAT(V->array[i]));
      fprintf(fic, " ");
    }
  fprintf (fic, "\n");
}

/**
 * Print a TYPE(PnlVect) as a row vector
 *
 * @param V a(constant) TYPE(PnlVect).
 */
void FUNCTION(pnl_vect, print_asrow)(const TYPE(PnlVect) * V)
{
  FUNCTION(pnl_vect, fprint_asrow)(stdout, V);
}

/**
 * Print a TYPE(PnlVect) in a format compatible with Nsp to a file
 *
 * @param V a(constant) TYPE(PnlVect).
 * @param fic a file descriptor.
 */
void FUNCTION(pnl_vect, fprint_nsp)(FILE *fic, const TYPE(PnlVect) * V)
{
  int i;
  fprintf(fic,"[ ");
  for (i=0; i<V->size-1; i++)
    {
      fprintf(fic,OUT_FORMAT,OUT_PUT_FORMAT(V->array[i]));
      fprintf(fic,"; ");
    }
  fprintf(fic,OUT_FORMAT,OUT_PUT_FORMAT(V->array[i]));
  fprintf(fic," ]; \n");
}

/**
 * Print a TYPE(PnlVect) in a format compatible with Nsp
 *
 * @param V a(constant) TYPE(PnlVect).
 */
void FUNCTION(pnl_vect, print_nsp)(const TYPE(PnlVect) * V)
{
  FUNCTION(pnl_vect, fprint_nsp)(stdout, V);
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
  int i;
  for ( i=0 ; i<lhs->size ; i++ )
    {
      const BASE xi = PNL_GET(lhs, i);
      PNL_LET(lhs, i) =(*f)(xi);
    }
}

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

static double FUNCTION(,_op_sqr_norm)(BASE a) { return SQUARE_NORM(a); }
static double FUNCTION(,_op_abs)(BASE a) { return NORMONE(a); }

/**
 * in-place vector scalar addition
 *
 * @param v left hand side vector
 * @param x scalar
 * @return  v = v+x
 */
void FUNCTION(pnl_vect,plus_scalar)(TYPE(PnlVect) *v , BASE x)
{
  int i;
  for ( i=0 ; i<v->size ; i++ )
    {
      const BASE vi = v->array[i];
      v->array[i] = PLUS(vi,x);
    }
}

/**
 * in-place vector scalar substraction
 *
 * @param v left hand side vector
 * @param x scalar
 * @return  v = v-x
 */
void FUNCTION(pnl_vect,minus_scalar)(TYPE(PnlVect) *v , BASE x)
{
  int i;
  for ( i=0 ; i<v->size ; i++ )
    {
      const BASE vi = v->array[i];
      v->array[i] = MINUS(vi,x);
    }
}

/**
 * in-place vector scalar multiplication
 *
 * @param v left hand side vector
 * @param x scalar
 * @return  v = v*x
 */
void FUNCTION(pnl_vect,mult_scalar)(TYPE(PnlVect) *v , BASE x)
{
  int i;
  for ( i=0 ; i<v->size ; i++ )
    {
      const BASE vi = v->array[i];
      v->array[i] = MULT(vi,x);
    }
}

/**
 * in-place vector scalar division
 *
 * @param v left hand side vector
 * @param x scalar
 * @return  v = v/x
 */
void FUNCTION(pnl_vect,div_scalar)(TYPE(PnlVect) *v , BASE x)
{
  int i;
  for ( i=0 ; i<v->size ; i++ )
    {
      const BASE vi = v->array[i];
      v->array[i] = DIV(vi,x);
    }
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
void FUNCTION(pnl_vect,map_vect_inplace)(TYPE(PnlVect) *lhs, const TYPE(PnlVect) *rhs, BASE(*f)(BASE,BASE))
{
  int i;
  CheckVectMatch(lhs, rhs);
  for ( i=0 ; i<lhs->size ; i++ )
    {
      PNL_LET(lhs, i) = (*f)(PNL_GET(lhs, i), PNL_GET(rhs,i)); 
    }
}

typedef struct {
  union {BASE x; TYPE(PnlVect) *V ;};
  char type;
} TYPE(cell);

/**
 * Find the indices i for which f == 1
 *
 * @param index (output) a vector of integers
 * @param type a string composed by the letters 'r' and 'v'.
 * @param f a function returning an integer (typically a test function)
 */
int FUNCTION(pnl_vect,find) (PnlVectInt *index, char* type, int(*f)(BASE *), ...)
{
  va_list ap;
  TYPE(cell) *args;
  int i, j, count, size, nvar;
  BASE val, *t;
  size = -1;

  nvar = strlen (type);
  if ((args = malloc (sizeof(cell) * nvar)) == NULL) return FAIL;
  if ((t = malloc (sizeof(BASE) * nvar)) == NULL) return FAIL;

  va_start (ap, f);

  for ( i=0; i<nvar ; i++ )
    {
      switch (type[i])
        {
          case 'r' :
            val = va_arg (ap, BASE);
            args[i].x = val; args[i].type = 'r';
            break;
          case 'v' :
            args[i].V= va_arg (ap, TYPE(PnlVect) *);
            args[i].type = 'v';
            if ( size == -1 ) size = args[i].V->size;
            else { PNL_CHECK ( size != args[i].V->size , "incompatible size", "pnl_vect_find"); }
            break;
        }
    }
  va_end(ap);

  /*
   * 2 passes are needed.
   * The first one to determine the size of index
   */
  for ( i=0, count=0 ; i<size ; i++ )
    {
       for ( j=0 ; j<nvar ; j++ )
         {
           if ( args[j].type == 'r' ) t[j] = args[j].x;
           else t[j] = PNL_GET(args[j].V, i);
         }
       if ( f(t) == 1 ) count++;
    }
  pnl_vect_int_resize (index, count);
  /*
   * Second pass to extract the indices for which f == 1
   */
  for ( i=0, count=0 ; i<size ; i++ )
    {
       for ( j=0 ; j<nvar ; j++ )
         {
           if ( args[j].type == 'r' ) t[j] = args[j].x;
           else t[j] = PNL_GET(args[j].V, i);
         }
       if ( f(t) == 1 ) { PNL_LET (index, count) = i; count++; }
    }

  free (args);
  free (t);
  return OK;
}

/**
 * map vector componentwise
 *
 * @param lhs each component lhs(i) contains f(rhs1(i),rhs2(i))
 * @param rhs1 a vector
 * @param rhs2 a vector
 * @param f real function
 */
void FUNCTION(pnl_vect,map_vect)(TYPE(PnlVect) *lhs, const TYPE(PnlVect) *rhs1, const TYPE(PnlVect) *rhs2, BASE(*f)(BASE,BASE))
{
  FUNCTION(pnl_vect,clone)(lhs,rhs1);
  FUNCTION(pnl_vect,map_vect_inplace) (lhs, rhs2, f);
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
  int i;
  CheckVectMatch(lhs, rhs);
  for ( i=0 ; i<lhs->size ; i++ )
    {
      PNL_LET(lhs, i) = PLUS(PNL_GET(lhs, i), PNL_GET(rhs,i)); 
    }
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
  int i;
  CheckVectMatch(lhs, rhs);
  for ( i=0 ; i<lhs->size ; i++ )
    {
      PNL_LET(lhs, i) = MINUS(PNL_GET(lhs, i), PNL_GET(rhs,i)); 
    }
}

/**
 * in-place term by term vector inverse
 *
 * @param v left hand side vector
 * @return  v = 1 ./ v
 */
void FUNCTION(pnl_vect,inv_term)(TYPE(PnlVect) *v)
{
  int i;
  for ( i=0 ; i<v->size ; i++ )
    {
      const BASE vi = v->array[i];
      v->array[i] = DIV(ONE, vi);
    }
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
  int i;
  CheckVectMatch(lhs, rhs);
  for ( i=0 ; i<lhs->size ; i++ )
    {
      PNL_LET(lhs, i) = DIV(PNL_GET(lhs, i), PNL_GET(rhs,i)); 
    }
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
  int i;
  CheckVectMatch(lhs, rhs);
  for ( i=0 ; i<lhs->size ; i++ )
    {
      PNL_LET(lhs, i) = MULT(PNL_GET(lhs, i), PNL_GET(rhs,i)); 
    }
}

/**
 * sum vector componentwise
 *
 * @param lhs left hand side vector
 * @return  sum=sum(lhs)
 */
BASE FUNCTION(pnl_vect,sum)(const TYPE(PnlVect) *lhs)
{
  int i;
  BASE sum=ZERO;
  for ( i=0 ; i<lhs->size ; i++ ) { sum = PLUS(sum, PNL_GET(lhs, i)); }
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
      sum = PLUS(sum, PNL_GET(lhs, i));
      FUNCTION(pnl_vect, set) (lhs, i, sum);
    }
}

/**
 * Compute the scalar product x' * y
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
  int i;
  CheckVectMatch(x, y);
  for ( i=0 ; i<x->size ; i++)
    {
      xi = PNL_GET(x, i);
      yi = PNL_GET(y, i);
      sum = PLUS(sum,(MULT(xi, yi)));
    }
  return sum;
}


/**
 * Compute the product of the components of a vector
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
      p = MULT(p,PNL_GET(V, i));
    }
  return p;
}

/**
 * Compute the cumulative product of the components of a vector
 *
 * @param V a a TYPE(PnlVect)
 */
void FUNCTION(pnl_vect,cumprod)(TYPE(PnlVect) *V)
{
  int i;
  BASE p=ONE;
  for (i=0; i<V->size; i++)
    {
      p = MULT(p,PNL_GET(V, i));
      FUNCTION(pnl_vect, set) (V, i, p);
    }
}

/**
 * Swap two elements of a vector
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
 * Return TRUE if a[i] < b[i] and a and b have the same length.
 * 
 * @param a
 * @param b
 * 
 * @return  TRUE or FALSE
 */
int FUNCTION(pnl_vect,less) (const TYPE(PnlVect) *a, const TYPE(PnlVect) *b)
{
  int i;
  if (a->size != b->size) return FAIL;
  for ( i=0 ; i<a->size ; i++ )
    {
      if ( (PNL_GET(a, i) >= PNL_GET(b, i)) ) return FALSE;
    }
  return TRUE;
}

/**
 * Return the index of the minimum of the elements of an array
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
      if (pnl_isnan(x))
        {
          *imin_out = i; *min_out = x;
          return;
        }
    }
  *imin_out = imin;
  *min_out = min;
}

/**
 * Return the index of the maximum of the elements of an array
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
      if (pnl_isnan(x))
        {
          *imax_out = i; *max_out = x;
          return;
        }
    }
  *imax_out = imax;
  *max_out = max;
}

/**
 * Return the indices of the minimum and maximum of the elements of an array
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
      if (pnl_isnan(x))
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
 * Return the maximum of the components of a vector
 *
 * @param V a vector
 * @return the maximum of the components of V
 */
BASE FUNCTION(pnl_vect,max)(const TYPE(PnlVect) *V)
{
  BASE m;
  int i;

  FUNCTION (pnl_vect, max_index)(&m, &i, V);
  return m;
}

/**
 * Return the minimum of the components of a vector
 *
 * @param V a vector
 * @return the minimum of the components of V
 */
BASE FUNCTION(pnl_vect,min)(const TYPE(PnlVect) *V)
{
  BASE m;
  int i;

  FUNCTION (pnl_vect, min_index)(&m, &i, V);
  return m;
}

/**
 * Return the minimum and maximum of the elements of a vector
 *
 * @param V a vector
 * @param min_out is a pointer to BASE, contains the minimum on exit
 * @param max_out is a pointer to BASE, contains the maximum on exit
 */
void FUNCTION(pnl_vect, minmax)(BASE *min_out, BASE *max_out, const TYPE(PnlVect) *V)
{
  int im, iM;
  FUNCTION(pnl_vect, minmax_index)(min_out, max_out, &im, &iM, V);
}

/**
 * Return the index of the maximum of the elements of a vector
 *
 * @param m a vector
 * @param imax_out is a pointer to int, contains the index of the maximum on exit
 * @param max_out contains the the maximum on exit
 */
void FUNCTION(pnl_vect, max_index)(BASE *max_out, int *imax_out, const TYPE(PnlVect) *m)
{
  FUNCTION(pnl_array, max_index)(m->array, m->size, 1, max_out, imax_out);
}

/**
 * Return the index of the minimum of the elements of a vector
 *
 * @param m a vector
 * @param imin_out is a pointer to BASE, contains the index of the minimum on
 * exit
 * @param min_out contains the the minimum on exit
 */
void FUNCTION(pnl_vect, min_index)(BASE *min_out, int *imin_out, const TYPE(PnlVect) *m)
{
  FUNCTION(pnl_array, min_index)(m->array, m->size, 1, min_out, imin_out);
}

/**
 * Return the indices of the minimum and maximum of the elements of a vector
 *
 * @param m a vector
 * @param imin_out is a pointer to BASE, contains the index of the minimum on exit
 * @param imax_out is a pointer to BASE, contains the index of the maximum on
 * exit
 * @param min_out contains the minimum on exit
 * @param max_out contains the maximum on exit
 */
void FUNCTION(pnl_vect, minmax_index)(BASE *min_out, BASE *max_out,
                                      int *imin_out, int *imax_out,
                                      const TYPE(PnlVect) *m)
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

/**
 * Compute a x + b y and stores the result in y
 *
 * @param a BASE
 * @param x a vector
 * @param b BASE
 * @param y a vector
 */
void FUNCTION(pnl_vect,axpby)(BASE a, const TYPE(PnlVect) *x,  BASE b, TYPE(PnlVect) *y)
{
  int i;
  BASE zero, one;
  zero = ZERO;
  one = ONE;

  if ( EQ(b,zero) )
    {
      FUNCTION(pnl_vect,resize) (y, x->size);
      FUNCTION(pnl_vect,set_all)(y, zero);
    }
  else if ( NEQ(b,one) ) FUNCTION(pnl_vect,mult_scalar) (y, b);
  if ( EQ(a,zero) ) return;

#ifndef PNL_RANGE_CHECK_OFF
  if (x->size != y->size)
    {
      PNL_ERROR ("size mismatch", "pnl_vect_axpby");
    }
#endif

  for ( i=0 ; i<x->size ; i++ )
    {
      BASE xi =  PNL_GET(x, i);
      BASE *yi = &(PNL_GET(y, i));
      xi = MULT(xi, a);
      *yi = PLUS(*yi, xi);
    }
}

/**
 * Return the X norm of a vector
 *
 * @param V a TYPE(PnlVect)
 * @param f a pointer's function on a norm function
 * @return the X norm of V
 */
static double FUNCTION(pnl_vect,norm_x)(const TYPE(PnlVect) *V,double(*f)(BASE))
{
  int i;
  double p=0;
  for ( i=0 ; i<V->size ; i++ )
    { 
      p += (*f)(PNL_GET(V, i));
    }
  return p;
}

/**
 * Return the two norm of a vector
 *
 * @param V  TYPE(PnlVect)
 * @return  the square root of the sum of the square value of components of V
 */
double FUNCTION(pnl_vect,norm_two)(const TYPE(PnlVect) *V)
{
  return sqrt(FUNCTION(pnl_vect,norm_x(V,FUNCTION(,_op_sqr_norm))));
}

/**
 * Return the one norm of a vector
 *
 * @param V a vector
 * @return the sum of the absolute value of components of V
 */
double FUNCTION(pnl_vect,norm_one)(const TYPE(PnlVect) *V)
{
  return FUNCTION(pnl_vect,norm_x)(V,FUNCTION(,_op_abs));
}

/**
 * Return the infty norm of a vector
 *
 * @param V a TYPE(PnlVect)
 * @return the maximum of the absolute value of components of V
 */
double FUNCTION(pnl_vect,norm_infty)(const TYPE(PnlVect) *V)
{
  int i;
  double p=0., q=0.;
  for ( i=0 ; i<V->size ; i++ )
    {
      q = NORMONE(PNL_GET(V,i));
      p = MAX(p, q);
    }
  return p;
}

/**
 * Extract a sub vector and wrap it into a vector.
 * @param V a vector
 * @param i the index of first element to be extracted
 * @param s the size of extracted vector
 * @return a vector (not a pointer) whose array pointer is the address of the
 * i-th element of V. No copying is done. This is a container for
 * V(i:i+s-1). The length of the xtracted vector is s.
 */
TYPE(PnlVect) FUNCTION(pnl_vect, wrap_subvect)(const TYPE(PnlVect) *V, int i,int s)
{
  TYPE(PnlVect) ret;
#ifndef PNL_RANGE_CHECK_OFF
  if ( (i+s-1 >= V->size) || (s<=0) || (i<0) )
    {
      PNL_ERROR ("index out of range", "pnl_vect_wrap_subvect");
    }
#endif
  FUNCTION(pnl_vect,init)(&ret);
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
 * i-th element of V. No copying is done. This is a container for V(i:j)
 */
TYPE(PnlVect) FUNCTION(pnl_vect, wrap_subvect_with_last)(const TYPE(PnlVect) *V, int i,int j)
{
  return FUNCTION(pnl_vect, wrap_subvect)(V,i,j-i+1);
}

/**
 * Wrap a matrix into a PnlVect
 * @param M a matrix
 * @return a vector (not a pointer) whose array pointer is the address of the
 * first element of the matrix M. No copying is done.
 */
TYPE(PnlVect) FUNCTION(pnl_vect,wrap_mat)(const TYPE(PnlMat) *M)
{
  TYPE(PnlVect) V;
  FUNCTION(pnl_vect,init)(&V);
  V.size = M->mn;
  V.mem_size = 0;
  V.owner = 0;
  V.array = M->array;
  return V;
}

#ifdef PNL_CLANG_COMPLETE
#include "pnl/pnl_templates_off.h"
#undef  BASE_DOUBLE
#endif
