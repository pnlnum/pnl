#ifndef _PNL_VECTOR_DOUBLE_H
#define _PNL_VECTOR_DOUBLE_H

#include <stdlib.h>
#include "pnl/pnl_object.h"

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
 * \defgroup PnlVect Double Vector 
 */
/*@{*/


#ifdef PNL_PNL_HAVE_INLINE 
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
extern PnlVect* pnl_vect_create(int size);
extern PnlVect pnl_vect_wrap_array(const double *x, int size);
extern PnlVect* pnl_vect_create_from_zero(int size);
extern PnlVect* pnl_vect_create_from_double(int size, double x);
extern PnlVect* pnl_vect_create_from_ptr(int size, const double* x);
extern PnlVect* pnl_vect_create_from_list(int size,...);
extern PnlVect* pnl_vect_create_from_file (const char * file);
extern PnlVect* pnl_vect_create_subvect_with_ind (const PnlVect *V, const PnlVectInt *ind);
extern void pnl_vect_extract_subvect_with_ind (PnlVect *V_sub, const PnlVect *V, const PnlVectInt *ind);
extern PnlVect* pnl_vect_create_subvect (const PnlVect *V, int i, int len);
extern void pnl_vect_extract_subvect (PnlVect *V_sub, const PnlVect *V, int i, int len);
extern int pnl_vect_resize(PnlVect *v, int size);
extern int pnl_vect_resize_from_double(PnlVect *v, int size, double x);
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
extern void pnl_vect_plus_double(PnlVect *lhs, double x); 
extern void pnl_vect_minus_double(PnlVect *lhs, double x);
extern void pnl_vect_axpby(double a, const PnlVect *x, double b, PnlVect *y); 
extern void pnl_vect_mult_double(PnlVect *lhs, double x);
extern void pnl_vect_div_double(PnlVect *lhs, double x); 
extern void pnl_vect_inv_term(PnlVect *lhs); 
extern void pnl_vect_div_vect_term(PnlVect *lhs, const PnlVect *rhs);
extern void pnl_vect_mult_vect_term(PnlVect *lhs, const PnlVect *rhs);
extern void pnl_vect_set_double(PnlVect *v, double x);
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
  union 
    {
      double val; /*!< single value */
      double *array; /*!< Pointer to double values */
    };
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
extern void pnl_vect_compact_set_double (PnlVectCompact *C, double x);
extern void pnl_vect_compact_set_ptr (PnlVectCompact *C, double *ptr);

/*@}*/
/*@}*/



#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_VECTOR_DOUBLE_H */
