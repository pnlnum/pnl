#ifndef VECTOR_COMPLEX_H
#define VECTOR_COMPLEX_H

#include "pnl_complex.h"

/**
 * \defgroup PnlVectComplex Complex Vector structure for Premia
 */
/*@{*/
typedef struct PnlVectComplex{
  int size;/*!< size of the vector */ 
  int mem_size; /*!< size of the memory block allocated for array */
  fcomplex *array;/*!< pointer to store the data */
  int owner; /*!< 1 if the structure owns its array pointer */
} PnlVectComplex;

#ifdef HAVE_INLINE 
extern inline
fcomplex pnl_vect_complex_get (const PnlVectComplex *self, int i)
{
  CheckIndexVect(self,i);
  return self->array[i];
}

extern inline
fcomplex* pnl_vect_complex_lget (PnlVectComplex *self, int i)
{
  CheckIndexVect(self,i);
  return &(self->array[i]);
}

extern inline
void pnl_vect_complex_set (PnlVectComplex *self, int i, fcomplex x)
{
  CheckIndexVect(self,i);
  self->array[i]=x;
}

extern inline
double pnl_vect_complex_get_real (const PnlVectComplex *v, int i)
{
  return ((double *)(v->array))[2*i];
}

extern inline
double pnl_vect_complex_get_imag (const PnlVectComplex *v, int i)
{
  return ((double *)(v->array))[2*i+1];
}

extern inline
double* pnl_vect_complex_lget_real (const PnlVectComplex *v, int i)
{
  return &(((double *)(v->array))[2*i]);
}

extern inline
double* pnl_vect_complex_lget_imag (const PnlVectComplex *v, int i)
{
  return &(((double *)(v->array))[2*i+1]);
}

extern inline
void pnl_vect_complex_set_real (const PnlVectComplex *v, int i, double re)
{
  ((double *)(v->array))[2*i] = re;
}

extern inline
void pnl_vect_complex_set_imag (const PnlVectComplex *v, int i, double im)
{
  ((double *)(v->array))[2*i+1] = im;
}
#endif

#ifndef PNL_RANGE_CHECK_OFF
#define GET_REAL(v,i) pnl_vect_complex_get_real(v,i)
#define LET_REAL(v,i) *(pnl_vect_complex_lget_real(v,i))
#define GET_IMAG(v,i) pnl_vect_complex_get_imag(v,i)
#define LET_IMAG(v,i) *(pnl_vect_complex_lget_imag(v,i))
#else
#define GET_REAL(v,i) (v->array[i]).re
#define LET_REAL(v,i) (v->array[i]).re
#define GET_IMAG(v,i) (v->array[i]).im
#define LET_IMAG(v,i) (v->array[i]).im
#endif


extern void pnl_vect_complex_set(PnlVectComplex *v, int i, fcomplex x);
extern fcomplex pnl_vect_complex_get(const PnlVectComplex *v, int i);
extern fcomplex* pnl_vect_complex_lget(PnlVectComplex *v, int i);
extern double pnl_vect_complex_get_real (const PnlVectComplex *v, int i);
extern double pnl_vect_complex_get_imag (const PnlVectComplex *v, int i);
extern double* pnl_vect_complex_lget_real (const PnlVectComplex *v, int i);
extern double* pnl_vect_complex_lget_imag (const PnlVectComplex *v, int i);
extern void pnl_vect_complex_set_real (const PnlVectComplex *v, int i, double re);
extern void pnl_vect_complex_set_imag (const PnlVectComplex *v, int i, double im);

extern void pnl_vect_complex_free(PnlVectComplex **v);
extern PnlVectComplex* pnl_vect_complex_create(int size);
extern PnlVectComplex pnl_vect_complex_create_wrap_array(const fcomplex *x, int size);
extern PnlVectComplex* pnl_vect_complex_create_from_fcomplex(int size, fcomplex x);
extern PnlVectComplex* pnl_vect_complex_create_from_array(int size, const double *re, const double *im);
extern PnlVectComplex* pnl_vect_complex_create_from_ptr(int size, const fcomplex* x);
extern PnlVectComplex* pnl_vect_complex_create_from_list(int size,...);
extern PnlVectComplex* pnl_vect_complex_create_from_file (const char * file);

extern int pnl_vect_complex_resize(PnlVectComplex *v, int size);
extern int pnl_vect_complex_resize_from_complex(PnlVectComplex *v, int size, fcomplex x);
extern int pnl_vect_complex_resize_from_ptr(PnlVectComplex *v, int size, const fcomplex *t);
extern PnlVectComplex* pnl_vect_complex_copy(const PnlVectComplex *v);
extern void pnl_vect_complex_clone(PnlVectComplex *clone, const PnlVectComplex *v);
extern PnlVectComplex pnl_vect_complex_wrap_subvect(const PnlVectComplex *V, int i,int s);
extern PnlVectComplex pnl_vect_complex_wrap_subvect_with_last(const PnlVectComplex *V, int i,int j);

extern void pnl_vect_complex_print(const PnlVectComplex *V);
extern void pnl_vect_complex_print_nsp(const PnlVectComplex *V);
extern void pnl_vect_complex_fprint(FILE *fic, const PnlVectComplex *V);

extern void pnl_vect_complex_plus_vect(PnlVectComplex *lhs, const PnlVectComplex *rhs); 
extern void pnl_vect_complex_minus_vect(PnlVectComplex *lhs, const PnlVectComplex *rhs); 
extern void pnl_vect_complex_map_inplace(PnlVectComplex *lhs, fcomplex(*f)(fcomplex)); /*lhs=f(lhs)*/
extern void pnl_vect_complex_minus(PnlVectComplex *lhs);
extern void pnl_vect_complex_plus_fcomplex(PnlVectComplex *lhs, fcomplex x); /*lhs+=x*/
extern void pnl_vect_complex_minus_fcomplex(PnlVectComplex *lhs, fcomplex x); /*lhs-=x*/
extern void pnl_vect_complex_axpby(fcomplex a, const PnlVectComplex *x, fcomplex b, PnlVectComplex *y); 
extern void pnl_vect_complex_mult_fcomplex(PnlVectComplex *lhs, fcomplex x); /*lhs*=x*/
extern void pnl_vect_complex_div_fcomplex(PnlVectComplex *lhs, fcomplex x); /*lhs*=x*/
extern void pnl_vect_complex_inv_term(PnlVectComplex *lhs); /* lhs = 1 ./ lhs*/
extern void pnl_vect_complex_mult_double(PnlVectComplex *lhs , double x);
extern void
pnl_vect_complex_div_vect_term(PnlVectComplex *lhs, const PnlVectComplex *rhs);/* lhs = lhs ./ rhs*/
extern void
pnl_vect_complex_mult_vect_term(PnlVectComplex *lhs, const PnlVectComplex *rhs); /* lhs= lhs.*rhs */ 
extern void pnl_vect_complex_set_fcomplex(PnlVectComplex *v, fcomplex x);/* v[j]= x */
extern void pnl_vect_set_zero(PnlVect * v); /* v[j]= 0 */
extern fcomplex pnl_vect_complex_sum(const PnlVectComplex *lhs);/* sum(x) */
extern void pnl_vect_complex_cumsum(PnlVectComplex *lhs);
extern void
pnl_vect_complex_map(PnlVectComplex *lhs, const PnlVectComplex *rhs,
                fcomplex(*f)(fcomplex));/* lhs(i)=f(rhs(i)) */
extern fcomplex
pnl_vect_complex_scalar_prod(const PnlVectComplex *rhs1, const PnlVectComplex *rhs2); /*rhs1.rhs2*/
extern fcomplex pnl_vect_complex_prod(const PnlVectComplex *V);
extern void pnl_vect_complex_cumprod(PnlVectComplex *V);

extern double pnl_vect_complex_norm_two(const PnlVectComplex *V); /*res=\Vert V \Vert_{l^2} */
extern double pnl_vect_complex_norm_one(const PnlVectComplex *V); /*res=\Vert V \Vert_{l^1} */
extern double pnl_vect_complex_norm_infty(const PnlVectComplex *V); /*res=\Vert V \Vert_{l^\infty} */
extern double pnl_vect_complex_norm_x(const PnlVectComplex *V,double(*f)(fcomplex)); /*res=\Vert V \Vert_{l^X} */
extern void pnl_vect_complex_split_in_array(PnlVectComplex* v, double *re, double *im);
extern void pnl_vect_complex_split_in_vect(PnlVectComplex* v, PnlVect *re, PnlVect *im);
extern void pnl_vect_complex_swap_elements(PnlVectComplex * v, int i, int j); 
extern void pnl_vect_complex_reverse(PnlVectComplex * v);

/*@}*/



#endif /* VECTOR_FCOMPLEX_H */


