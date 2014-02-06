#ifndef _PNL_VECTOR_H
#define _PNL_VECTOR_H

#include <stdio.h>
#include <stdlib.h>

#ifndef PNL_RANGE_CHECK_OFF
#define CheckIndexVect(v,i) {                                           \
    if (i>=v->size || i<0) {perror("index out of range"); abort();}}
#define CheckVectMatch(lhs, rhs) { if ((lhs)->size != (rhs)->size)  \
      {perror("non compatible dimensions"); abort();}}
#else
#define CheckIndexVect(v,i) {}
#define CheckVectMatch(lhs, rhs){}
#endif /* PNL_RANGE_CHECK_OFF */

#include "pnl/pnl_config.h"

/**
 * \defgroup PnlVect  a Vector object
 */

#define PNL_GET(v,i) (v)->array[i]
#define PNL_LET(v,i) (v)->array[i]
#define PNL_SET(v,i,x) (v)->array[i]=(x)

#include "pnl/pnl_matvect.h"
  
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

extern PnlVectObject* pnl_vect_object_new ();
extern void pnl_vect_object_free (PnlVectObject **);
extern int pnl_vect_object_resize(PnlVectObject * v, int size);


/*
 * PnlVect
 */

/**
 * \ingroup PnlVect
 */
/*@{*/


#ifdef PNL_HAVE_INLINE 
PNL_INLINE_FUNC double pnl_vect_get (const PnlVect *self, int i)
{
  CheckIndexVect(self,i);
  return self->array[i];
}

PNL_INLINE_FUNC double* pnl_vect_lget (PnlVect *self, int i)
{
  CheckIndexVect(self,i);
  return &(self->array[i]);
}

PNL_INLINE_FUNC void pnl_vect_set (PnlVect *self, int i, double x)
{
  CheckIndexVect(self,i);
  self->array[i]=x;
}
#endif

#ifndef PNL_RANGE_CHECK_OFF
#define GET(v,i) pnl_vect_get(v,i)
#define LET(v,i) *(pnl_vect_lget(v,i))
#else
#define GET(v,i) (v)->array[i]
#define LET(v,i) (v)->array[i]
#endif

PNL_INLINE_DECL void pnl_vect_set(PnlVect *v, int i, double x);
PNL_INLINE_DECL double pnl_vect_get(const PnlVect *v, int i);
PNL_INLINE_DECL double* pnl_vect_lget(PnlVect *v, int i);

extern void pnl_vect_free(PnlVect **v);
extern void pnl_vect_init(PnlVect *v);
extern PnlVect* pnl_vect_new();
extern int pnl_vect_eq (const PnlVect *, const PnlVect *);
extern int pnl_vect_eq_all (const PnlVect *, double);
extern int pnl_vect_less(const PnlVect * a,const PnlVect * b);
extern PnlVect* pnl_vect_create(int size);
extern PnlVect pnl_vect_wrap_array(const double *x, int size);
extern PnlVect* pnl_vect_create_from_zero(int size);
extern PnlVect* pnl_vect_create_from_scalar(int size, double x);
extern PnlVect* pnl_vect_create_from_ptr(int size, const double* x);
extern PnlVect* pnl_vect_create_from_list(int size,...);
extern PnlVect* pnl_vect_create_from_file (const char * file);
extern PnlVect* pnl_vect_create_subvect_with_ind (const PnlVect *V, const PnlVectInt *ind);
extern void pnl_vect_extract_subvect_with_ind (PnlVect *V_sub, const PnlVect *V, const PnlVectInt *ind);
extern PnlVect* pnl_vect_create_subvect (const PnlVect *V, int i, int len);
extern void pnl_vect_extract_subvect (PnlVect *V_sub, const PnlVect *V, int i, int len);
extern int pnl_vect_resize(PnlVect *v, int size);
extern int pnl_vect_resize_from_scalar(PnlVect *v, int size, double x);
extern int pnl_vect_resize_from_ptr(PnlVect *v, int size, const double *t);
extern PnlVect* pnl_vect_copy(const PnlVect *v);
extern void pnl_vect_clone(PnlVect *clone, const PnlVect *v);
extern PnlVect pnl_vect_wrap_subvect(const PnlVect *V, int i,int s);
extern PnlVect pnl_vect_wrap_subvect_with_last(const PnlVect *V, int i,int j);
extern PnlVect pnl_vect_wrap_mat(const PnlMat *M);
extern void pnl_vect_print(const PnlVect *V);
extern void pnl_vect_print_asrow(const PnlVect *V);
extern void pnl_vect_print_nsp(const PnlVect *V);
extern void pnl_vect_fprint(FILE *fic, const PnlVect *V);
extern void pnl_vect_fprint_asrow(FILE *fic, const PnlVect *V);
extern void pnl_vect_fprint_nsp(FILE *fic, const PnlVect *V);
extern void pnl_vect_plus_vect(PnlVect *lhs, const PnlVect *rhs); 
extern void pnl_vect_minus_vect(PnlVect *lhs, const PnlVect *rhs); 


extern void pnl_vect_map_inplace(PnlVect *lhs, double(*f)(double)); 
extern void pnl_vect_map(PnlVect *lhs, const PnlVect *rhs, double(*f)(double));
extern void pnl_vect_map_vect_inplace(PnlVect *lhs, const PnlVect *rhs, double(*f)(double, double));
extern void pnl_vect_map_vect(PnlVect *lhs, const PnlVect *rhs1, const PnlVect *rhs2, double(*f)(double, double));
extern int pnl_vect_find(PnlVectInt *ind, char *type, int(*f)(double *), ...);
extern void pnl_vect_minus(PnlVect *lhs);
extern void pnl_vect_plus_scalar(PnlVect *lhs, double x); 
extern void pnl_vect_minus_scalar(PnlVect *lhs, double x);
extern void pnl_vect_axpby(double a, const PnlVect *x, double b, PnlVect *y); 
extern void pnl_vect_mult_scalar(PnlVect *lhs, double x);
extern void pnl_vect_div_scalar(PnlVect *lhs, double x); 
extern void pnl_vect_inv_term(PnlVect *lhs); 
extern void pnl_vect_div_vect_term(PnlVect *lhs, const PnlVect *rhs);
extern void pnl_vect_mult_vect_term(PnlVect *lhs, const PnlVect *rhs);
extern void pnl_vect_set_all(PnlVect *v, double x);
extern void pnl_vect_set_zero(PnlVect * v);
extern double pnl_vect_sum(const PnlVect *lhs);
extern void pnl_vect_cumsum(PnlVect *lhs);
extern double pnl_vect_scalar_prod(const PnlVect *rhs1, const PnlVect *rhs2);
extern double pnl_vect_dist (const PnlVect *x, const PnlVect *y);
extern int pnl_vect_cross(PnlVect *lhs, const PnlVect *x, const PnlVect *y);
extern double pnl_vect_prod(const PnlVect *V); 
extern void pnl_vect_cumprod(PnlVect *V); 
extern double pnl_vect_max(const PnlVect *V);
extern double pnl_vect_min(const PnlVect *V);
extern void pnl_vect_minmax (double *, double *, const PnlVect *);
extern void pnl_vect_min_index (double *, int *, const PnlVect *);
extern void pnl_vect_max_index (double *, int *, const PnlVect *);
extern void pnl_vect_minmax_index (double *, double *, int *, int *, const PnlVect *);
extern void pnl_vect_qsort (PnlVect *, char);
extern void pnl_vect_qsort_index (PnlVect *, PnlVectInt *, char);

extern double pnl_vect_norm_two(const PnlVect *V);
extern double pnl_vect_norm_one(const PnlVect *V);
extern double pnl_vect_norm_infty(const PnlVect *V);
extern void pnl_vect_swap_elements(PnlVect * v, int i, int j); 
extern void pnl_vect_reverse(PnlVect * v);

/**
 * Compact PnlVect : used for variables that can either contain a single
  * number or a PnlVect.
  * vectors likes x*ones(n,1) are simply stored as a double x 
  */
typedef struct PnlVectCompact 
{
  PnlObject object;
  int size; /*!< size of the vector */
  double val; /*!< single value */
  double *array; /*!< Pointer to double values */
  char convert; /*!< 'a', 'd' : array, double */
} PnlVectCompact;

extern PnlVectCompact* pnl_vect_compact_new ();
extern PnlVectCompact* pnl_vect_compact_create_from_ptr (int n, double const *x);
extern PnlVectCompact* pnl_vect_compact_create (int n, double x);
extern int pnl_vect_compact_resize (PnlVectCompact *v, int size, double x);
extern PnlVectCompact* pnl_vect_compact_copy(const PnlVectCompact *v);
extern void pnl_vect_compact_free (PnlVectCompact **v);
extern PnlVect* pnl_vect_compact_to_pnl_vect (const PnlVectCompact *C);
extern double pnl_vect_compact_get (const PnlVectCompact *C, int i);
extern void pnl_vect_compact_set_all (PnlVectCompact *C, double x);
extern void pnl_vect_compact_set_ptr (PnlVectCompact *C, double *ptr);


/*
 * PnlVectInt
 */


#ifdef PNL_HAVE_INLINE 
PNL_INLINE_FUNC int pnl_vect_int_get (const PnlVectInt *self, int i)
{
  CheckIndexVect(self,i);
  return self->array[i];
}

PNL_INLINE_FUNC int* pnl_vect_int_lget (PnlVectInt *self, int i)
{
  CheckIndexVect(self,i);
  return &(self->array[i]);
}

PNL_INLINE_FUNC void pnl_vect_int_set (PnlVectInt *self, int i, int x)
{
  CheckIndexVect(self,i);
  self->array[i]=x;
}
#endif

PNL_INLINE_DECL void pnl_vect_int_set(PnlVectInt *v, int i, int x);
PNL_INLINE_DECL int pnl_vect_int_get(const PnlVectInt *v, int i);
PNL_INLINE_DECL int* pnl_vect_int_lget(PnlVectInt *v, int i);

#ifndef PNL_RANGE_CHECK_OFF
#define GET_INT(v,i) pnl_vect_int_get(v,i)
#define LET_INT(v,i) *(pnl_vect_int_lget(v,i))
#else
#define GET_INT(v,i) (v)->array[i]
#define LET_INT(v,i) (v)->array[i]
#endif


extern void pnl_vect_int_free(PnlVectInt **v);
extern void pnl_vect_int_init(PnlVectInt *v);
extern PnlVectInt* pnl_vect_int_new();
extern int pnl_vect_int_eq (const PnlVectInt *, const PnlVectInt *);
extern int pnl_vect_int_eq_all (const PnlVectInt *, int);
extern int pnl_vect_int_less(const PnlVectInt * a,const PnlVectInt * b);
extern PnlVectInt* pnl_vect_int_create(int size);
extern PnlVectInt pnl_vect_int_wrap_array(const int *x, int size);
extern PnlVectInt* pnl_vect_int_create_from_scalar(int size, int x);
extern PnlVectInt* pnl_vect_int_create_from_zero(int size);
extern PnlVectInt* pnl_vect_int_create_from_ptr(int size, const int* x);
extern PnlVectInt* pnl_vect_int_create_from_list(int size,...);
extern PnlVectInt* pnl_vect_int_create_from_file (const char * file);
extern PnlVectInt* pnl_vect_int_create_subvect_with_ind (const PnlVectInt *V, const PnlVectInt *ind);
extern void pnl_vect_int_extract_subvect_with_ind (PnlVectInt *V_sub, const PnlVectInt *V, const PnlVectInt *ind);
extern PnlVectInt* pnl_vect_int_create_subvect (const PnlVectInt *V, int i, int len);
extern void pnl_vect_int_extract_subvect (PnlVectInt *V_sub, const PnlVectInt *V, int i, int len);
extern int pnl_vect_int_resize(PnlVectInt *v, int size);
extern int pnl_vect_int_resize_from_scalar(PnlVectInt *v, int size, int x);
extern int pnl_vect_int_resize_from_ptr(PnlVectInt *v, int size, const int *t);
extern PnlVectInt* pnl_vect_int_copy(const PnlVectInt *v);
extern void pnl_vect_int_clone(PnlVectInt *clone, const PnlVectInt *v);
extern PnlVectInt pnl_vect_int_wrap_subvect(const PnlVectInt *V, int i,int s);
extern PnlVectInt pnl_vect_int_wrap_subvect_with_last(const PnlVectInt *V, int i,int j);
extern PnlVectInt pnl_vect_int_wrap_mat(const PnlMatInt *M);
extern void pnl_vect_int_print(const PnlVectInt *V);
extern void pnl_vect_int_print_asrow(const PnlVectInt *V);
extern void pnl_vect_int_print_nsp(const PnlVectInt *V);
extern void pnl_vect_int_fprint(FILE *fic, const PnlVectInt *V);
extern void pnl_vect_int_fprint_asrow(FILE *fic, const PnlVectInt *V);
extern void pnl_vect_int_fprint_nsp(FILE *fic, const PnlVectInt *V);
extern void pnl_vect_int_plus_vect(PnlVectInt *lhs, const PnlVectInt *rhs); 
extern void pnl_vect_int_minus_vect(PnlVectInt *lhs, const PnlVectInt *rhs); 
extern void pnl_vect_int_minus(PnlVectInt *lhs);
extern void pnl_vect_int_map_inplace(PnlVectInt *lhs, int(*f)(int)); /*lhs=f(lhs)*/
extern void pnl_vect_int_map(PnlVectInt *lhs, const PnlVectInt *rhs, int(*f)(int));
extern void pnl_vect_int_map_vect_inplace(PnlVectInt *lhs, const PnlVectInt *rhs, int(*f)(int, int)); 
extern void pnl_vect_int_map_vect(PnlVectInt *lhs, const PnlVectInt *rhs1, const PnlVectInt *rhs2, int(*f)(int, int));
extern int pnl_vect_int_find(PnlVectInt *ind, char *type, int(*f)(int *), ...);
extern void pnl_vect_int_plus_scalar(PnlVectInt *lhs, int x); /*lhs+=x*/
extern void pnl_vect_int_minus_scalar(PnlVectInt *lhs, int x); /*lhs-=x*/
extern void pnl_vect_int_axpby(int a, const PnlVectInt *x, int b, PnlVectInt *y); /* y:=a x + b y */
extern void pnl_vect_int_mult_scalar(PnlVectInt *lhs, int x); /*lhs*=x*/
extern void pnl_vect_int_div_scalar(PnlVectInt *lhs, int x); /*lhs*=x*/
extern void pnl_vect_int_inv_term(PnlVectInt *lhs); /* lhs = 1 ./ lhs*/
extern void pnl_vect_int_div_vect_term(PnlVectInt *lhs, const PnlVectInt *rhs);/* lhs = lhs ./ rhs*/
extern void pnl_vect_int_mult_vect_term(PnlVectInt *lhs, const PnlVectInt *rhs); /* lhs= lhs.*rhs */ 
extern void pnl_vect_int_set_all(PnlVectInt *v, int x);/* v[j]= x */
extern void pnl_vect_int_set_zero(PnlVectInt *v); /* v[j]= 0 */
extern int pnl_vect_int_sum(const PnlVectInt *lhs);/* sum(x) */
extern void pnl_vect_int_cumsum(PnlVectInt *lhs);
extern int pnl_vect_int_scalar_prod(const PnlVectInt *rhs1, const PnlVectInt *rhs2); /*rhs1.rhs2*/
extern int pnl_vect_int_prod(const PnlVectInt *V); /*res=prod(V(i))*/
extern void pnl_vect_int_cumprod(PnlVectInt *V); /*res=prod(V(i))*/
extern int pnl_vect_int_max(const PnlVectInt *V); /*res=max(V)*/
extern int pnl_vect_int_min(const PnlVectInt *V); /*res=min(V)*/
extern void pnl_vect_int_minmax (int *, int *,const PnlVectInt *);
extern void pnl_vect_int_min_index (int *, int *, const PnlVectInt *);
extern void pnl_vect_int_max_index (int *, int *,const PnlVectInt *);
extern void pnl_vect_int_minmax_index (int *, int *, int *, int *,const PnlVectInt *);
extern void pnl_vect_int_qsort (PnlVectInt *, char);
extern void pnl_vect_int_qsort_index (PnlVectInt *, PnlVectInt *, char);

extern double pnl_vect_int_norm_two(const PnlVectInt *V); /*res=\Vert V \Vert_{l^2} */
extern double pnl_vect_int_norm_one(const PnlVectInt *V); /*res=\Vert V \Vert_{l^1} */
extern double pnl_vect_int_norm_infty(const PnlVectInt *V); /*res=\Vert V \Vert_{l^\infty} */

extern int pnl_vect_int_level(PnlVectInt *v, int i);
extern int pnl_vect_int_level_norm_one(PnlVectInt *v);
extern int pnl_vect_int_level_norm_inf(PnlVectInt *v);
extern void pnl_vect_int_dyadic_cast(const PnlVectInt * v_int,PnlVect * v_out);
extern void pnl_vect_int_swap_elements(PnlVectInt * v, int i, int j); 
extern void pnl_vect_int_reverse(PnlVectInt * v);



/*
 * PnlVectComplex
 */

#include "pnl/pnl_complex.h"


#ifdef PNL_HAVE_INLINE 
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
#define GET_COMPLEX(v,index) pnl_vect_complex_get(v,index)
#define LET_COMPLEX(v,index) *(pnl_vect_complex_lget(v,index))
#else
#define GET_REAL(v,index) ((v)->array[index]).r
#define LET_REAL(v,index) ((v)->array[index]).r
#define GET_IMAG(v,index) ((v)->array[index]).i
#define LET_IMAG(v,index) ((v)->array[index]).i
#define GET_COMPLEX(v,index) (v)->array[index]
#define LET_COMPLEX(v,index) (v)->array[index]
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
extern int pnl_vect_complex_eq_all (const PnlVectComplex *, dcomplex);
extern PnlVectComplex* pnl_vect_complex_create(int size);
extern PnlVectComplex pnl_vect_complex_wrap_array(const dcomplex *x, int size);
extern PnlVectComplex* pnl_vect_complex_create_from_scalar(int size, dcomplex x);
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
extern int pnl_vect_complex_resize_from_scalar(PnlVectComplex *v, int size, dcomplex x);
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
extern void pnl_vect_complex_plus_scalar(PnlVectComplex *lhs, dcomplex x); /*lhs+=x*/
extern void pnl_vect_complex_minus_scalar(PnlVectComplex *lhs, dcomplex x); /*lhs-=x*/
extern void pnl_vect_complex_axpby(dcomplex a, const PnlVectComplex *x, dcomplex b, PnlVectComplex *y); 
extern void pnl_vect_complex_mult_scalar(PnlVectComplex *lhs, dcomplex x); /*lhs*=x*/
extern void pnl_vect_complex_div_scalar(PnlVectComplex *lhs, dcomplex x); /*lhs*=x*/
extern void pnl_vect_complex_inv_term(PnlVectComplex *lhs); /* lhs = 1 ./ lhs*/
extern void pnl_vect_complex_mult_double(PnlVectComplex *lhs , double x);
extern void
pnl_vect_complex_div_vect_term(PnlVectComplex *lhs, const PnlVectComplex *rhs);/* lhs = lhs ./ rhs*/
extern void
pnl_vect_complex_mult_vect_term(PnlVectComplex *lhs, const PnlVectComplex *rhs); /* lhs= lhs.*rhs */ 
extern void pnl_vect_complex_set_all(PnlVectComplex *v, dcomplex x);
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

/*
 * Some deprecated names
 */

#include "pnl/pnl_deprecated.h"


#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_VECTOR_H */
