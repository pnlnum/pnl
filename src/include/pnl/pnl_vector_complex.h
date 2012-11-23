#ifndef _PNL_VECTOR_COMPLEX_H
#define _PNL_VECTOR_COMPLEX_H

#include "pnl/pnl_object.h"
#include "pnl/pnl_complex.h"


#ifndef _PNL_VECTOR_H
#error "Do not include this file directly. Include pnl_vector.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**
 * \ingroup PnlVectors
 */
/*@{*/

/**
 * \defgroup PnlVectComplex Complex Vector 
 */
/*@{*/

#ifdef PNL_PNL_HAVE_INLINE 
PNL_INLINE_FUNC dcomplex pnl_vect_complex_get (const PnlVectComplex *self, int i)
{
  CheckIndexVect(self,i);
  return self->array[i];
}

PNL_INLINE_FUNC dcomplex* pnl_vect_complex_lget (PnlVectComplex *self, int i)
{
  CheckIndexVect(self,i);
  return &(self->array[i]);
}

PNL_INLINE_FUNC void pnl_vect_complex_set (PnlVectComplex *self, int i, dcomplex x)
{
  CheckIndexVect(self,i);
  self->array[i]=x;
}

PNL_INLINE_FUNC double pnl_vect_complex_get_real (const PnlVectComplex *v, int i)
{
  return ((double *)(v->array))[2*i];
}

PNL_INLINE_FUNC double pnl_vect_complex_get_imag (const PnlVectComplex *v, int i)
{
  return ((double *)(v->array))[2*i+1];
}

PNL_INLINE_FUNC double* pnl_vect_complex_lget_real (const PnlVectComplex *v, int i)
{
  return &(((double *)(v->array))[2*i]);
}

PNL_INLINE_FUNC double* pnl_vect_complex_lget_imag (const PnlVectComplex *v, int i)
{
  return &(((double *)(v->array))[2*i+1]);
}

PNL_INLINE_FUNC void pnl_vect_complex_set_real (const PnlVectComplex *v, int i, double re)
{
  ((double *)(v->array))[2*i] = re;
}

PNL_INLINE_FUNC void pnl_vect_complex_set_imag (const PnlVectComplex *v, int i, double im)
{
  ((double *)(v->array))[2*i+1] = im;
}
#endif

#ifndef PNL_RANGE_CHECK_OFF
#define GET_REAL(v,index) pnl_vect_complex_get_real(v,index)
#define LET_REAL(v,index) *(pnl_vect_complex_lget_real(v,index))
#define GET_IMAG(v,index) pnl_vect_complex_get_imag(v,index)
#define LET_IMAG(v,index) *(pnl_vect_complex_lget_imag(v,index))
#else
#define GET_REAL(v,index) ((v)->array[index]).r
#define LET_REAL(v,index) ((v)->array[index]).r
#define GET_IMAG(v,index) ((v)->array[index]).i
#define LET_IMAG(v,index) ((v)->array[index]).i
#endif


PNL_INLINE_DECL void pnl_vect_complex_set(PnlVectComplex *v, int i, dcomplex x);
PNL_INLINE_DECL dcomplex pnl_vect_complex_get(const PnlVectComplex *v, int i);
PNL_INLINE_DECL dcomplex* pnl_vect_complex_lget(PnlVectComplex *v, int i);
PNL_INLINE_DECL double pnl_vect_complex_get_real (const PnlVectComplex *v, int i);
PNL_INLINE_DECL double pnl_vect_complex_get_imag (const PnlVectComplex *v, int i);
PNL_INLINE_DECL double* pnl_vect_complex_lget_real (const PnlVectComplex *v, int i);
PNL_INLINE_DECL double* pnl_vect_complex_lget_imag (const PnlVectComplex *v, int i);
PNL_INLINE_DECL void pnl_vect_complex_set_real (const PnlVectComplex *v, int i, double re);
PNL_INLINE_DECL void pnl_vect_complex_set_imag (const PnlVectComplex *v, int i, double im);

extern void pnl_vect_complex_free(PnlVectComplex **v);
extern void pnl_vect_complex_init(PnlVectComplex *v);
extern PnlVectComplex* pnl_vect_complex_new();
extern int pnl_vect_complex_eq (const PnlVectComplex *, const PnlVectComplex *);
extern int pnl_vect_complex_eq_dcomplex (const PnlVectComplex *, dcomplex);
extern PnlVectComplex* pnl_vect_complex_create(int size);
extern PnlVectComplex pnl_vect_complex_wrap_array(const dcomplex *x, int size);
extern PnlVectComplex* pnl_vect_complex_create_from_dcomplex(int size, dcomplex x);
extern PnlVectComplex* pnl_vect_complex_create_from_zero(int size);
extern PnlVectComplex* pnl_vect_complex_create_from_array(int size, const double *re, const double *im);
extern PnlVectComplex* pnl_vect_complex_create_from_ptr(int size, const dcomplex* x);
extern PnlVectComplex* pnl_vect_complex_create_from_list(int size,...);
extern PnlVectComplex* pnl_vect_complex_create_from_file (const char * file);
extern PnlVectComplex* pnl_vect_complex_create_subvect_with_ind (const PnlVectComplex *V, const PnlVectInt *ind);
extern void pnl_vect_complex_extract_subvect_ind (PnlVectComplex *V_sub, const PnlVectComplex *V, const PnlVectInt *ind);
extern PnlVectComplex* pnl_vect_complex_create_subvect (const PnlVectComplex *V,  int i, int len);
extern void pnl_vect_complex_extract (PnlVectComplex *V_sub, const PnlVectComplex *V, int i, int len);

extern int pnl_vect_complex_resize(PnlVectComplex *v, int size);
extern int pnl_vect_complex_resize_from_dcomplex(PnlVectComplex *v, int size, dcomplex x);
extern int pnl_vect_complex_resize_from_ptr(PnlVectComplex *v, int size, const dcomplex *t);
extern PnlVectComplex* pnl_vect_complex_copy(const PnlVectComplex *v);
extern void pnl_vect_complex_clone(PnlVectComplex *clone, const PnlVectComplex *v);
extern PnlVectComplex pnl_vect_complex_wrap_subvect(const PnlVectComplex *V, int i,int s);
extern PnlVectComplex pnl_vect_complex_wrap_subvect_with_last(const PnlVectComplex *V, int i,int j);
extern PnlVectComplex pnl_vect_complex_wrap_mat(const PnlMatComplex *M);

extern void pnl_vect_complex_print(const PnlVectComplex *V);
extern void pnl_vect_complex_print_asrow(const PnlVectComplex *V);
extern void pnl_vect_complex_print_nsp(const PnlVectComplex *V);
extern void pnl_vect_complex_fprint(FILE *fic, const PnlVectComplex *V);
extern void pnl_vect_complex_fprint_asrow(FILE *fic, const PnlVectComplex *V);
extern void pnl_vect_complex_fprint_nsp(FILE *fic, const PnlVectComplex *V);

extern void pnl_vect_complex_plus_vect(PnlVectComplex *lhs, const PnlVectComplex *rhs); 
extern void pnl_vect_complex_minus_vect(PnlVectComplex *lhs, const PnlVectComplex *rhs); 
extern void pnl_vect_complex_map_inplace(PnlVectComplex *lhs, dcomplex(*f)(dcomplex));
extern void pnl_vect_complex_map(PnlVectComplex *lhs, const PnlVectComplex *rhs, dcomplex(*f)(dcomplex));
extern void pnl_vect_complex_map_vect_inplace(PnlVectComplex *lhs, const PnlVectComplex *rhs, dcomplex(*f)(dcomplex,dcomplex));
extern void pnl_vect_complex_map_vect(PnlVectComplex *lhs, const PnlVectComplex *rhs1, const PnlVectComplex *rhs2, dcomplex(*f)(dcomplex,dcomplex));
extern int pnl_vect_complex_find(PnlVectInt *ind, char *type, int(*f)(dcomplex *), ...);  
extern void pnl_vect_complex_minus(PnlVectComplex *lhs);
extern void pnl_vect_complex_plus_dcomplex(PnlVectComplex *lhs, dcomplex x); /*lhs+=x*/
extern void pnl_vect_complex_minus_dcomplex(PnlVectComplex *lhs, dcomplex x); /*lhs-=x*/
extern void pnl_vect_complex_axpby(dcomplex a, const PnlVectComplex *x, dcomplex b, PnlVectComplex *y); 
extern void pnl_vect_complex_mult_dcomplex(PnlVectComplex *lhs, dcomplex x); /*lhs*=x*/
extern void pnl_vect_complex_div_dcomplex(PnlVectComplex *lhs, dcomplex x); /*lhs*=x*/
extern void pnl_vect_complex_inv_term(PnlVectComplex *lhs); /* lhs = 1 ./ lhs*/
extern void pnl_vect_complex_mult_double(PnlVectComplex *lhs , double x);
extern void
pnl_vect_complex_div_vect_term(PnlVectComplex *lhs, const PnlVectComplex *rhs);/* lhs = lhs ./ rhs*/
extern void
pnl_vect_complex_mult_vect_term(PnlVectComplex *lhs, const PnlVectComplex *rhs); /* lhs= lhs.*rhs */ 
extern void pnl_vect_complex_set_dcomplex(PnlVectComplex *v, dcomplex x);
extern void pnl_vect_complex_set_zero(PnlVectComplex *v);
extern void pnl_vect_set_zero(PnlVect * v); /* v[j]= 0 */
extern dcomplex pnl_vect_complex_sum(const PnlVectComplex *lhs);/* sum(x) */
extern void pnl_vect_complex_cumsum(PnlVectComplex *lhs);
extern dcomplex
pnl_vect_complex_scalar_prod(const PnlVectComplex *rhs1, const PnlVectComplex *rhs2); /*rhs1.rhs2*/
extern dcomplex pnl_vect_complex_prod(const PnlVectComplex *V);
extern void pnl_vect_complex_cumprod(PnlVectComplex *V);

extern double pnl_vect_complex_norm_two(const PnlVectComplex *V); /*res=\Vert V \Vert_{l^2} */
extern double pnl_vect_complex_norm_one(const PnlVectComplex *V); /*res=\Vert V \Vert_{l^1} */
extern double pnl_vect_complex_norm_infty(const PnlVectComplex *V); /*res=\Vert V \Vert_{l^\infty} */
extern void pnl_vect_complex_split_in_array(const PnlVectComplex* v, double *re, double *im);
extern void pnl_vect_complex_split_in_vect(const PnlVectComplex* v, PnlVect *re, PnlVect *im);
extern void pnl_vect_complex_swap_elements(PnlVectComplex * v, int i, int j); 
extern void pnl_vect_complex_reverse(PnlVectComplex * v);

/*@}*/
/*@}*/



#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_VECTOR_FCOMPLEX_H */
