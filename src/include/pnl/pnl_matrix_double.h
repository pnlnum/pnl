#ifndef _PNL_MATRIX_DOUBLE_H
#define _PNL_MATRIX_DOUBLE_H

#ifndef _PNL_MATRIX_H
#error "Do not include this file directly. Include pnl_matrix.h"
#endif

#include "pnl/pnl_perm.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


#ifndef PNL_RANGE_CHECK_OFF

#define MGET(v,i,j) pnl_mat_get((v), (i), (j))
#define MLET(v,i,j) *(pnl_mat_lget((v), (i), (j)))

#else /* PNL_RANGE_CHECK_OFF */

#define MGET(v,i,j) (v)->array[(i)*(v)->n+(j)]
#define MLET(v,i,j) (v)->array[(i)*(v)->n+(j)]

#endif /* PNL_RANGE_CHECK_OFF */


/**
 * \ingroup PnlMatrices
 */
/*@{*/

/**
 * \defgroup PnlMat Double Matrix 
 */
/*@{*/


extern void pnl_mat_init(PnlMat *); 
extern PnlMat* pnl_mat_new(); 
extern int pnl_mat_eq (const PnlMat *, const PnlMat *);
extern int pnl_mat_eq_double (const PnlMat *, double);
extern PnlMat* pnl_mat_create(int m, int n); 
extern PnlMat* pnl_mat_create_from_double(int m, int n, double x);
extern PnlMat* pnl_mat_create_from_ptr(int m, int n, const double* x);
extern PnlMat* pnl_mat_create_from_list(int m, int n, ...); 
extern PnlMat pnl_mat_wrap_array(const double* x,int m, int n);
extern PnlMat* pnl_mat_create_from_file (const char * file);
extern void pnl_vect_extract_submat (PnlVect *M_sub, const PnlMat *M, const PnlVectInt *indi, const PnlVectInt *indj);
extern PnlVect* pnl_vect_create_submat (const PnlMat *M, const PnlVectInt *indi, const PnlVectInt *indj);
extern void pnl_mat_extract_subblock (PnlMat *M_sub, const PnlMat *M, int i, int len_i, int j, int len_j);
extern int pnl_mat_resize(PnlMat *v, int m, int n);
extern void pnl_mat_free(PnlMat **v);
extern PnlMat* pnl_mat_copy(const PnlMat *v);
extern void pnl_mat_clone(PnlMat *clone, const PnlMat *M);
extern void pnl_mat_map_inplace(PnlMat *lhs, double(*f)(double)); /*lhs=f(lhs)*/
extern void pnl_mat_map(PnlMat *lhs, const PnlMat *rhs, double(*f)(double));/* lhs(i)=f(rhs(i)) */
extern void pnl_mat_map_mat_inplace(PnlMat *lhs, const PnlMat *rhs, double(*f)(double,double));
extern void pnl_mat_map_mat(PnlMat *lhs, const PnlMat *rhs1, const PnlMat *rhs2, double(*f)(double,double));
extern int pnl_mat_find(PnlVectInt *indi, PnlVectInt *indj, char* type, int(*f)(double *), ...);
extern void pnl_mat_plus_double(PnlMat *lhs, double x); /*lhs+=x*/
extern void pnl_mat_minus_double(PnlMat *lhs, double x); /*lhs-=x*/
extern void pnl_mat_plus_mat(PnlMat *lhs, const PnlMat *rhs); /*lhs+=rhs*/
extern void pnl_mat_minus_mat(PnlMat *lhs, const PnlMat *rhs); /*lhs-=rhs*/
extern void pnl_mat_mult_double(PnlMat *lhs, double x); /*lhs*=x*/
extern void pnl_mat_div_double(PnlMat *lhs, double x); /*lhs/=x*/
extern PnlMat* pnl_mat_mult_mat(const PnlMat *rhs1, const PnlMat *rhs2);
extern void pnl_mat_mult_mat_inplace(PnlMat *lhs, const PnlMat *rhs1,
                                     const PnlMat *rhs2);/*lhs=rhs1*rhs2*/
extern void pnl_mat_dgemm (char transA, char transB, double alpha, const PnlMat *A,
                           const PnlMat *B, double beta, PnlMat *C); 
extern void pnl_mat_axpy (double a, const PnlMat *X, PnlMat *Y); 
extern void pnl_mat_dger (double alpha, const PnlVect *x, const PnlVect *y, PnlMat *A);
extern PnlVect* pnl_mat_mult_vect(const PnlMat *mat, const PnlVect *vec);
extern void pnl_mat_mult_vect_inplace(PnlVect *lhs, const PnlMat *mat,
                                      const PnlVect *rhs);

extern PnlVect* pnl_mat_mult_vect_transpose(const PnlMat *mat, const PnlVect *vec);
extern void pnl_mat_mult_vect_transpose_inplace(PnlVect *lhs, const PnlMat *mat, const PnlVect *rhs);
extern  void pnl_mat_lAxpby(double l, const PnlMat *A, const PnlVect *x, double b, PnlVect * y);
extern void pnl_mat_dgemv (char trans, double alpha, const PnlMat *A,
                           const PnlVect *x , double beta, PnlVect * y);
extern double pnl_mat_scalar_prod(const PnlMat *A, const PnlVect *x , const PnlVect * y);
extern int pnl_mat_cross(PnlMat *lhs, const PnlMat *A, const PnlMat *B);

extern void pnl_mat_set_double(PnlMat *lhs, double x);
extern void pnl_mat_set_zero(PnlMat *lhs);
extern void pnl_mat_set_id(PnlMat *lhs);
extern void pnl_mat_set_diag(PnlMat *lhs, double x, int d);
extern void pnl_mat_set_from_ptr(PnlMat *lhs, const double *x);
extern void pnl_mat_mult_mat_term(PnlMat *lhs, const PnlMat *rhs); 
extern void pnl_mat_div_mat_term(PnlMat *lhs, const PnlMat *rhs); 
extern void pnl_mat_sq_transpose(PnlMat *M);
extern PnlMat* pnl_mat_transpose(const PnlMat *M);
extern void pnl_mat_tr(PnlMat*, const PnlMat *M);
extern void pnl_mat_print(const PnlMat *M);
extern void pnl_mat_fprint(FILE *fic, const PnlMat *M);
extern void pnl_mat_fprint_nsp (FILE *fic, const PnlMat *M);
extern void pnl_mat_print_nsp (const PnlMat *M);
extern PnlVect pnl_vect_wrap_mat_row(const PnlMat *M, int i);
extern PnlMat pnl_mat_wrap_vect(const PnlVect *v);
extern void pnl_mat_get_row(PnlVect *V, const PnlMat *M, int i);
extern void pnl_mat_get_col(PnlVect *V, const PnlMat *M, int j);
extern PnlMat* pnl_mat_create_diag_from_ptr(const double x[], int d);
extern PnlMat* pnl_mat_create_diag(const PnlVect *V);
extern double pnl_mat_sum(const PnlMat *lhs);
extern double pnl_mat_prod(const PnlMat *lhs);
extern void pnl_mat_sum_vect(PnlVect *y, const PnlMat *A, char c);
extern void pnl_mat_cumsum(PnlMat *A, char a);
extern void pnl_mat_prod_vect(PnlVect *y, const PnlMat *A, char c);
extern void pnl_mat_cumprod(PnlMat *A, char a);
extern void pnl_mat_min (PnlVect *out, const PnlMat *A, char d);
extern void pnl_mat_max (PnlVect *out, const PnlMat *A, char d);
extern void pnl_mat_minmax (PnlVect *m, PnlVect *M, const PnlMat *A, char d);
extern void pnl_mat_min_index (PnlVect *out, PnlVectInt *i, const PnlMat *A, char d);
extern void pnl_mat_max_index (PnlVect *out, PnlVectInt *i, const PnlMat *A, char d);
extern void pnl_mat_minmax_index (PnlVect *m, PnlVect *M, PnlVectInt *im, PnlVectInt *iM, const PnlMat *A, char d);
extern void pnl_mat_qsort_index (PnlMat *A, PnlMatInt *t, char dir, char order);
extern void pnl_mat_qsort (PnlMat *A, char dir, char order);
extern void pnl_mat_set_row(PnlMat *M, const PnlVect *V, int i);
extern void pnl_mat_set_col(PnlMat *M, const PnlVect *V, int j);
extern void pnl_mat_swap_rows (PnlMat *M, int i, int j);
extern int pnl_mat_upper_inverse(PnlMat *A, const PnlMat *B);
extern int pnl_mat_lower_inverse (PnlMat *A, const PnlMat *B);
extern int pnl_mat_chol(PnlMat *x);
extern int pnl_mat_pchol (PnlMat *M, double tol, int *rank, PnlVectInt *p);
extern int pnl_mat_upper_syslin (PnlVect *x, const PnlMat *U, const  PnlVect *b);
extern int pnl_mat_lower_syslin (PnlVect *x, const PnlMat *L, const  PnlVect *b);
extern int pnl_mat_chol_syslin (PnlVect *x, const PnlMat *chol, const  PnlVect *b);
extern int pnl_mat_chol_syslin_inplace (const PnlMat *chol, PnlVect *b);
extern int pnl_mat_lu (PnlMat *A, PnlVectInt *p);
extern int pnl_mat_lu_syslin (PnlVect *x, PnlMat *LU,
                               const PnlVectInt *p, const PnlVect *b);
extern int pnl_mat_lu_syslin_inplace(PnlMat *LU, const PnlVectInt *p, PnlVect *b);
extern int pnl_mat_qr_syslin (PnlVect *x, const PnlMat *Q, const PnlMat *R, const PnlVectInt *p, const PnlVect *b);
extern int pnl_mat_syslin_mat (PnlMat *A,  PnlMat *B);
extern int pnl_mat_lu_syslin_mat (const PnlMat *A,  const PnlPermutation *p, PnlMat *B);
extern int pnl_mat_chol_syslin_mat (const PnlMat *A,  PnlMat *B);
extern int pnl_mat_syslin (PnlVect *x, const PnlMat *A, const PnlVect *b);
extern int pnl_mat_syslin_inplace (PnlMat *A, PnlVect *b);
extern int pnl_mat_inverse (PnlMat *inverse, const PnlMat *A);
extern int pnl_mat_inverse_with_chol (PnlMat *inverse, const PnlMat *A);
extern int pnl_mat_exp (PnlMat *B, const PnlMat *A);
extern int pnl_mat_log (PnlMat *B, const PnlMat *A);
extern int pnl_mat_eigen (PnlVect *v, PnlMat *P, const PnlMat *A, int with_eigenvectors);
extern int pnl_mat_qr (PnlMat *Q, PnlMat *R, PnlVectInt *p, const PnlMat *A);
extern int pnl_mat_ls_mat (const PnlMat *A, PnlMat *B);
extern int pnl_mat_ls (const PnlMat *A, PnlVect *b);

/* inline functions if you are using GCC */
#ifdef PNL_PNL_HAVE_INLINE 
PNL_INLINE_FUNC double pnl_mat_get (const PnlMat *self, int i, int j)
{
  CheckIndexMat(self,i, j);
  return self->array[i*self->n+j];
}

PNL_INLINE_FUNC double* pnl_mat_lget (PnlMat *self, int i, int j)
{
  CheckIndexMat(self,i, j);
  return &(self->array[i*self->n+j]);
}

PNL_INLINE_FUNC void pnl_mat_set (PnlMat *self, int i, int j, double x)
{
  CheckIndexMat(self,i, j);
  self->array[i*self->n+j] = x;
}
#endif

PNL_INLINE_DECL void pnl_mat_set(PnlMat *v, int i, int j, double x);
PNL_INLINE_DECL double pnl_mat_get(const PnlMat *v, int i, int j);
PNL_INLINE_DECL double* pnl_mat_lget(PnlMat *v, int i, int j);

/*@}*/
/*@}*/

/**
 * \ingroup PnlHmatrices
 */
/*@{*/

/**
 * \defgroup PnlHmat Double HyperMatrix 
 */
/*@{*/


extern PnlHmat* pnl_hmat_new(); 
extern PnlHmat* pnl_hmat_create(int ndim, const int *dims); 
extern PnlHmat* pnl_hmat_create_from_double(int ndim, const int *dims, double x); 
extern PnlHmat* pnl_hmat_create_from_ptr(int ndim, const int *dims, const double *x);
extern int pnl_hmat_resize(PnlHmat *v, int ndim, const int *dims);
extern void pnl_hmat_free(PnlHmat **v);
extern PnlHmat* pnl_hmat_copy(const PnlHmat *H);
extern void pnl_hmat_clone(PnlHmat *clone, const PnlHmat *H);
extern void pnl_hmat_print(const PnlHmat *H);
extern void pnl_hmat_plus_hmat(PnlHmat *lhs, const PnlHmat *rhs);/*lhs+=rhs*/
extern void pnl_hmat_mult_double(PnlHmat *lhs, double x);/* lhs *=x;*/

extern void pnl_hmat_set(PnlHmat *H, int *tab, double x);
extern double pnl_hmat_get(const PnlHmat *H, int *tab);
extern double* pnl_hmat_lget(PnlHmat *H, int *tab);
extern PnlMat pnl_mat_wrap_hmat(PnlHmat *H, int *t);
extern PnlVect pnl_vect_wrap_hmat(PnlHmat *H, int *t);

/*@}*/
/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_MATRIX_DOUBLE_H */
