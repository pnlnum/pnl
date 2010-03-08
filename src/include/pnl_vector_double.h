#ifndef _PNL_VECTOR_DOUBLE_H
#define _PNL_VECTOR_DOUBLE_H

#ifndef _PNL_VECTOR_H
#error "Do not include this file directly. Include pnl_vector.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdlib.h>
#include "pnl_matrix.h"

/**
 * \ingroup PnlVectors
 */
/*@{*/
/**
 * \addtogroup PnlVect 
 */

/*@{*/
struct _PnlVect {
  int size;/*!< size of the vector */ 
  double *array;/*!< pointer to store the data */
  int mem_size; /*!< size of the memory block allocated for array */
  int owner; /*!< 1 if the object owns its array member, 0 otherwise */
};


#ifdef HAVE_INLINE 
extern inline
double pnl_vect_get (const PnlVect *self, int i)
{
  CheckIndexVect(self,i);
  return self->array[i];
}

extern inline
double* pnl_vect_lget (PnlVect *self, int i)
{
  CheckIndexVect(self,i);
  return &(self->array[i]);
}

extern inline
void pnl_vect_set (PnlVect *self, int i, double x)
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

extern void pnl_vect_set(PnlVect *v, int i, double x);
extern double pnl_vect_get(const PnlVect *v, int i);
extern double* pnl_vect_lget(PnlVect *v, int i);
extern void pnl_vect_free(PnlVect **v);
extern PnlVect* pnl_vect_create(int size);
extern PnlVect pnl_vect_wrap_array(const double *x, int size);
extern PnlVect* pnl_vect_create_from_zero(int size);
extern PnlVect* pnl_vect_create_from_double(int size, double x);
extern PnlVect* pnl_vect_create_from_ptr(int size, const double* x);
extern PnlVect* pnl_vect_create_from_list(int size,...);
extern PnlVect* pnl_vect_create_from_file (const char * file);
extern int pnl_vect_resize(PnlVect *v, int size);
extern int pnl_vect_resize_from_double(PnlVect *v, int size, double x);
extern int pnl_vect_resize_from_ptr(PnlVect *v, int size, const double *t);
extern PnlVect* pnl_vect_copy(const PnlVect *v);
extern void pnl_vect_clone(PnlVect *clone, const PnlVect *v);
extern PnlVect pnl_vect_wrap_subvect(const PnlVect *V, int i,int s);
extern PnlVect pnl_vect_wrap_subvect_with_last(const PnlVect *V, int i,int j);
extern PnlVect pnl_vect_wrap_mat(const PnlMat *M);
extern void pnl_vect_print(const PnlVect *V);
extern void pnl_vect_print_nsp(const PnlVect *V);
extern void pnl_vect_fprint(FILE *fic, const PnlVect *V);
extern void pnl_vect_fprint_nsp(FILE *fic, const PnlVect *V);
extern void pnl_vect_plus_vect(PnlVect *lhs, const PnlVect *rhs); 
extern void pnl_vect_minus_vect(PnlVect *lhs, const PnlVect *rhs); 

/*lhs+=rhs*/

extern void pnl_vect_map_inplace(PnlVect *lhs, double(*f)(double)); /*lhs=f(lhs)*/
extern void pnl_vect_minus(PnlVect *lhs);
extern void pnl_vect_plus_double(PnlVect *lhs, double x); /*lhs+=x*/
extern void pnl_vect_minus_double(PnlVect *lhs, double x); /*lhs-=x*/
extern void pnl_vect_axpby(double a, const PnlVect *x, double b, PnlVect *y); /* y := a x + b y */
extern void pnl_vect_mult_double(PnlVect *lhs, double x); /*lhs*=x*/
extern void pnl_vect_div_double(PnlVect *lhs, double x); /*lhs*=x*/
extern void pnl_vect_inv_term(PnlVect *lhs); /* lhs = 1 ./ lhs*/
extern void
pnl_vect_div_vect_term(PnlVect *lhs, const PnlVect *rhs);/* lhs = lhs ./ rhs*/
extern void
pnl_vect_mult_vect_term(PnlVect *lhs, const PnlVect *rhs); /* lhs= lhs.*rhs */ 
extern void pnl_vect_set_double(PnlVect *v, double x);/* v[j]= x */
extern void pnl_vect_set_zero(PnlVect * v); /* v[j]= 0 */
extern double pnl_vect_sum(const PnlVect *lhs);/* sum(x) */
extern void pnl_vect_cumsum(PnlVect *lhs);
extern void
pnl_vect_map(PnlVect *lhs, const PnlVect *rhs,
                double(*f)(double));/* lhs(i)=f(rhs(i)) */
extern double
pnl_vect_scalar_prod(const PnlVect *rhs1, const PnlVect *rhs2); /*rhs1.rhs2*/
extern double pnl_vect_prod(const PnlVect *V); /*res=prod(V(i))*/
extern void pnl_vect_cumprod(PnlVect *V); /*res=prod(V(i))*/
extern double pnl_vect_max(const PnlVect *V); /*res=max(V)*/
extern double pnl_vect_min(const PnlVect *V); /*res=min(V)*/
extern void pnl_vect_minmax (const PnlVect *, double *, double *);
extern void pnl_vect_min_index (const PnlVect *, double *, int *);
extern void pnl_vect_max_index (const PnlVect *, double *, int *);
extern void pnl_vect_minmax_index (const PnlVect *, double *, double *, int *, int *);
extern void pnl_vect_qsort (PnlVect *, char);
extern void pnl_vect_qsort_index (PnlVect *, PnlVectInt *, char);

extern double pnl_vect_norm_two(const PnlVect *V); /*res=\Vert V \Vert_{l^2} */
extern double pnl_vect_norm_one(const PnlVect *V); /*res=\Vert V \Vert_{l^1} */
extern double pnl_vect_norm_infty(const PnlVect *V); /*res=\Vert V \Vert_{l^\infty} */
extern void pnl_vect_swap_elements(PnlVect * v, int i, int j); 
extern void pnl_vect_reverse(PnlVect * v);
/*@}*/
/*@}*/



#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_VECTOR_DOUBLE_H */
