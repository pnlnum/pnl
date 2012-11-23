#ifndef _PNL_MATRIX_COMPLEX_H
#define _PNL_MATRIX_COMPLEX_H

#ifndef _PNL_MATRIX_H
#error "Do not include this file directly. Include pnl_matrix.h"
#endif

#include "pnl/pnl_complex.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/**
 * \ingroup PnlMatrices
 */
/*@{*/

/**
 * \defgroup PnlMatComplex Complex Matrix
 */
/*@{*/




extern void pnl_mat_complex_init(PnlMatComplex *); 
extern PnlMatComplex* pnl_mat_complex_new(); 
extern int pnl_mat_complex_eq (const PnlMatComplex *, const PnlMatComplex *);
extern int pnl_mat_complex_eq_domplex (const PnlMatComplex *, dcomplex);
extern PnlMatComplex* pnl_mat_complex_create(int m, int n); 
extern PnlMatComplex* pnl_mat_complex_create_from_dcomplex(int m, int n, dcomplex x);
extern PnlMatComplex* pnl_mat_complex_create_from_ptr(int m, int n, const dcomplex* x);
extern PnlMatComplex* pnl_mat_complex_create_from_list(int m, int n, ...); 
extern PnlMatComplex pnl_mat_complex_wrap_array(const dcomplex* x, int m, int n);
extern PnlMatComplex* pnl_mat_complex_create_from_file (const char * file);
extern void pnl_vect_complex_extract_submat (PnlVectComplex *M_sub, const PnlMatComplex *M, const PnlVectInt *indi, const PnlVectInt *indj);
extern PnlVectComplex* pnl_vect_complex_create_submat (const PnlMatComplex *M, const PnlVectInt *indi, const PnlVectInt *indj);
extern void pnl_mat_complex_extract_subblock (PnlMatComplex *M_sub, const PnlMatComplex *M, int i, int len_i, int j, int len_j);
extern int pnl_mat_complex_resize(PnlMatComplex *v, int m, int n);
extern void pnl_mat_complex_free(PnlMatComplex **v);
extern PnlMatComplex* pnl_mat_complex_copy(const PnlMatComplex *v);
extern void pnl_mat_complex_clone(PnlMatComplex *clone, const PnlMatComplex *M);
extern void pnl_mat_complex_set_dcomplex(PnlMatComplex *lhs, dcomplex x);/*lhs=x*/
extern void pnl_mat_complex_set_zero(PnlMatComplex *lhs);
extern void pnl_mat_complex_set_id(PnlMatComplex *lhs);
extern void pnl_mat_complex_set_diag(PnlMatComplex *lhs, dcomplex x, int d);
extern void pnl_mat_complex_set_from_ptr(PnlMatComplex *lhs, const dcomplex *x);
extern PnlMatComplex* pnl_mat_complex_transpose(const PnlMatComplex *M);
extern void pnl_mat_complex_tr(PnlMatComplex*, const PnlMatComplex *M);
extern void pnl_mat_complex_sq_transpose(PnlMatComplex *M);
extern void pnl_mat_complex_print(const PnlMatComplex *M);
extern void pnl_mat_complex_print_nsp(const PnlMatComplex *M);
extern void pnl_mat_complex_fprint(FILE *fic,const PnlMatComplex *M);
extern void pnl_mat_complex_fprint_nsp(FILE *fic,const PnlMatComplex *M);
extern PnlMatComplex* pnl_mat_complex_create_diag_from_ptr(const dcomplex x[], int d);
extern PnlMatComplex* pnl_mat_complex_create_diag(const PnlVectComplex *V);
extern void
pnl_mat_complex_set_row(PnlMatComplex *M, const PnlVectComplex *V, int i);/* M(i,:)=V(:) */
extern void
pnl_mat_complex_set_col(PnlMatComplex *M, const PnlVectComplex *V, int j);/* M(:,j)=V(:) */
extern void
pnl_mat_complex_swap_rows (PnlMatComplex *M, int i, int j);
extern void
pnl_mat_complex_get_row(PnlVectComplex *V, const PnlMatComplex *M, int i);/* V(:)=M(i,:) */
extern void pnl_mat_complex_get_col(PnlVectComplex *V, const PnlMatComplex *M, int j);
extern PnlVectComplex pnl_vect_complex_wrap_mat_row(const PnlMatComplex *M, int i);/* M(i,:)=V(:) */
extern PnlMatComplex pnl_mat_complex_wrap_vect(const PnlVectComplex *V);
extern void pnl_mat_complex_map_inplace(PnlMatComplex *lhs, dcomplex(*f)(dcomplex)); /*lhs=f(lhs)*/
extern void pnl_mat_complex_map(PnlMatComplex *lhs, const PnlMatComplex *rhs, dcomplex(*f)(dcomplex));/* lhs(i)=f(rhs(i)) */
extern void pnl_mat_complex_map_mat_inplace(PnlMatComplex *lhs, const PnlMatComplex *rhs, dcomplex(*f)(dcomplex,dcomplex));
extern void pnl_mat_complex_map_mat(PnlMatComplex *lhs, const PnlMatComplex *rhs1, const PnlMatComplex *rhs2, dcomplex(*f)(dcomplex,dcomplex));
extern int pnl_mat_complex_find(PnlVectInt *indi, PnlVectInt *indj, char* type, int(*f)(dcomplex *), ...);
extern void pnl_mat_complex_plus_dcomplex(PnlMatComplex *lhs, dcomplex x); /*lhs+=x*/
extern void pnl_mat_complex_minus_dcomplex(PnlMatComplex *lhs, dcomplex x); /*lhs-=x*/
extern void pnl_mat_complex_plus_mat(PnlMatComplex *lhs, const PnlMatComplex *rhs); /*lhs+=rhs*/
extern void pnl_mat_complex_minus_mat(PnlMatComplex *lhs, const PnlMatComplex *rhs); /*lhs-=rhs*/
extern void pnl_mat_complex_mult_dcomplex(PnlMatComplex *lhs, dcomplex x); /*lhs*=x*/
extern void pnl_mat_complex_div_dcomplex(PnlMatComplex *lhs, dcomplex x); /*lhs/=x*/
extern void
pnl_mat_complex_mult_mat_term(PnlMatComplex *lhs, const PnlMatComplex *rhs); /*lhs=lhs.*rhs*/
extern void
pnl_mat_complex_div_mat_term(PnlMatComplex *lhs, const PnlMatComplex *rhs); /*lhs=lhs./rhs*/
extern void pnl_mat_complex_axpy (dcomplex a, const PnlMatComplex *X, PnlMatComplex *Y); 
extern void pnl_mat_complex_dger (dcomplex alpha, const PnlVectComplex *x, const PnlVectComplex *y, PnlMatComplex *A);
extern PnlVectComplex* pnl_mat_complex_mult_vect(const PnlMatComplex *mat, const PnlVectComplex *vec);
extern void pnl_mat_complex_mult_vect_inplace(PnlVectComplex *lhs, const PnlMatComplex *mat,
                                              const PnlVectComplex *rhs);
extern PnlVectComplex* pnl_mat_complex_mult_vect_transpose(const PnlMatComplex *mat, const PnlVectComplex *vec);
extern void pnl_mat_complex_mult_vect_transpose_inplace(PnlVectComplex *lhs, const PnlMatComplex *mat, const PnlVectComplex *rhs);
extern  void pnl_mat_complex_lAxpby(dcomplex l, const PnlMatComplex *A, const PnlVectComplex *x, dcomplex b, PnlVectComplex * y);
extern void pnl_mat_complex_dgemv (char trans, dcomplex alpha, const PnlMatComplex *A,
                                   const PnlVectComplex *x , dcomplex beta, PnlVectComplex * y);
extern PnlMatComplex* pnl_mat_complex_mult_mat(const PnlMatComplex *rhs1, const PnlMatComplex *rhs2);
extern void pnl_mat_complex_mult_mat_inplace(PnlMatComplex *lhs, const PnlMatComplex *rhs1,
                                          const PnlMatComplex *rhs2);/*lhs=rhs1*rhs2*/
extern void pnl_mat_complex_dgemm (char transA, char transB, dcomplex alpha, const PnlMatComplex *A,
                                const PnlMatComplex *B, dcomplex beta, PnlMatComplex *C);
extern dcomplex pnl_mat_complex_sum(const PnlMatComplex *lhs);
extern dcomplex pnl_mat_complex_prod(const PnlMatComplex *lhs);
extern void pnl_mat_complex_sum_vect (PnlVectComplex *y, const PnlMatComplex *M, char a);
extern void pnl_mat_complex_cumsum(PnlMatComplex *lhs, char a);
extern void pnl_mat_complex_prod_vect (PnlVectComplex *y, const PnlMatComplex *M, char a);
extern void pnl_mat_complex_cumprod(PnlMatComplex *lhs, char a);

extern int pnl_mat_complex_upper_inverse(PnlMatComplex *A, const PnlMatComplex *B);
extern int pnl_mat_complex_lower_inverse (PnlMatComplex *A, const PnlMatComplex *B);
extern int pnl_mat_complex_chol(PnlMatComplex *x);
extern int pnl_mat_complex_pchol (PnlMatComplex *M, double tol, int *rank, PnlVectInt *p);
extern int pnl_mat_complex_upper_syslin (PnlVectComplex *x, const PnlMatComplex *U, const  PnlVectComplex *b);
extern int pnl_mat_complex_lower_syslin (PnlVectComplex *x, const PnlMatComplex *L, const  PnlVectComplex *b);
extern int pnl_mat_complex_chol_syslin (PnlVectComplex *x, const PnlMatComplex *chol, const  PnlVectComplex *b);
extern int pnl_mat_complex_chol_syslin_inplace (const PnlMatComplex *chol, PnlVectComplex *b);
extern int pnl_mat_complex_lu (PnlMatComplex *A, PnlVectInt *p);
extern int pnl_mat_complex_lu_syslin (PnlVectComplex *x, PnlMatComplex *LU,
                               const PnlVectInt *p, const PnlVectComplex *b);
extern int pnl_mat_complex_lu_syslin_inplace(PnlMatComplex *LU, const PnlVectInt *p, PnlVectComplex *b);
extern int pnl_mat_complex_syslin_mat (PnlMatComplex *A,  PnlMatComplex *B);
extern int pnl_mat_complex_lu_syslin_mat (const PnlMatComplex *A,  const PnlPermutation *p, PnlMatComplex *B);
extern int pnl_mat_complex_chol_syslin_mat (const PnlMatComplex *A,  PnlMatComplex *B);
extern int pnl_mat_complex_syslin (PnlVectComplex *x, const PnlMatComplex *A, const PnlVectComplex *b);
extern int pnl_mat_complex_syslin_inplace (PnlMatComplex *A, PnlVectComplex *b);
extern int pnl_mat_complex_inverse (PnlMatComplex *inverse, const PnlMatComplex *A);
extern int pnl_mat_complex_inverse_with_chol (PnlMatComplex *inverse, const PnlMatComplex *A);
extern int pnl_mat_complex_exp (PnlMatComplex *expA, const PnlMatComplex *A);


/* inline functions if you are using GCC */
#ifdef PNL_PNL_HAVE_INLINE 
PNL_INLINE_FUNC dcomplex pnl_mat_complex_get (const PnlMatComplex *self, int i, int j)
{
  CheckIndexMat(self,i, j);
  return self->array[i*self->n+j];
}

PNL_INLINE_FUNC dcomplex* pnl_mat_complex_lget (PnlMatComplex *self, int i, int j)
{
  CheckIndexMat(self,i, j);
  return &(self->array[i*self->n+j]);
}

PNL_INLINE_FUNC void pnl_mat_complex_set (PnlMatComplex *self, int i, int j, dcomplex x)
{
  CheckIndexMat(self,i, j);
  self->array[i*self->n+j] = x;
}

PNL_INLINE_FUNC double pnl_mat_complex_get_real (const PnlMatComplex *self, int i, int j)
{
  CheckIndexMat(self,i, j);
  return self->array[i*self->n+j].r;
}

PNL_INLINE_FUNC double* pnl_mat_complex_lget_real (PnlMatComplex *self, int i, int j)
{
  CheckIndexMat(self,i, j);
  return &(self->array[i*self->n+j].r);
}

PNL_INLINE_FUNC void pnl_mat_complex_set_real (PnlMatComplex *self, int i, int j, double x)
{
  CheckIndexMat(self,i, j);
  self->array[i*self->n+j].r = x; 
}

PNL_INLINE_FUNC double pnl_mat_complex_get_imag (const PnlMatComplex *self, int i, int j)
{
  CheckIndexMat(self,i, j);
  return (self->array[i*self->n+j]).i;
}

PNL_INLINE_FUNC double* pnl_mat_complex_lget_imag (PnlMatComplex *self, int i, int j)
{
  CheckIndexMat(self,i, j);
  return &((self->array[i*self->n+j]).i);
}

PNL_INLINE_FUNC void pnl_mat_complex_set_imag (PnlMatComplex *self, int i, int j, double x)
{
  CheckIndexMat(self,i, j);
  self->array[i*self->n+j].i = x;
}
#endif

PNL_INLINE_DECL void pnl_mat_complex_set(PnlMatComplex *v, int i, int j, dcomplex x);
PNL_INLINE_DECL dcomplex pnl_mat_complex_get(const PnlMatComplex *v, int i, int j);
PNL_INLINE_DECL dcomplex* pnl_mat_complex_lget(PnlMatComplex *v, int i, int j);
PNL_INLINE_DECL void pnl_mat_complex_set_real(PnlMatComplex *v, int i, int j, double x);
PNL_INLINE_DECL double pnl_mat_complex_get_real(const PnlMatComplex *v, int i, int j);
PNL_INLINE_DECL double* pnl_mat_complex_lget_real(PnlMatComplex *v, int i, int j);
PNL_INLINE_DECL void pnl_mat_complex_set_imag(PnlMatComplex *v, int i, int j, double x);
PNL_INLINE_DECL double pnl_mat_complex_get_imag(const PnlMatComplex *v, int i, int j);
PNL_INLINE_DECL double* pnl_mat_complex_lget_imag(PnlMatComplex *v, int i, int j);

#ifndef PNL_RANGE_CHECK_OFF
#define MGET_REAL(v,i,j) pnl_mat_complex_get_real(v,i,j)
#define MLET_REAL(v,i,j) *(pnl_mat_complex_lget_real(v,i,j))
#define MGET_IMAG(v,i,j) pnl_mat_complex_get_imag(v,i,j)
#define MLET_IMAG(v,i,j) *(pnl_mat_complex_lget_imag(v,i,j))
#else
#define MGET_REAL(v,i,j) ((v)->array[(i)*(v)->n+(j)]).r
#define MLET_REAL(v,i,j) ((v)->array[(i)*(v)->n+(j)]).r
#define MGET_IMAG(v,i,j) ((v)->array[(i)*(v)->n+(j)]).i
#define MLET_IMAG(v,i,j) ((v)->array[(i)*(v)->n+(j)]).i
#endif

/*@}*/
/*@}*/


/**
 * \ingroup PnlHmatrices a Hyper Matrix
 */
/*@{*/
/**
 * \defgroup PnlHmatComplex Complex HyperMatrix 
 */
/*@{*/


extern PnlHmatComplex* pnl_hmat_complex_new(); 
extern PnlHmatComplex* pnl_hmat_complex_create(int ndim, const int *dims); 
extern PnlHmatComplex* pnl_hmat_complex_create_from_dcomplex(int ndim, const int *dims, dcomplex x); 
extern PnlHmatComplex* pnl_hmat_complex_create_from_ptr(int ndim, const int *dims, const dcomplex *x);
extern int pnl_hmat_complex_resize(PnlHmatComplex *v, int ndim, const int *dims);
extern void pnl_hmat_complex_free(PnlHmatComplex **v);
extern void pnl_hmat_complex_print(const PnlHmatComplex *H);
extern PnlHmatComplex* pnl_hmat_complex_copy(const PnlHmatComplex *H);
extern void pnl_hmat_complex_clone(PnlHmatComplex *clone, const PnlHmatComplex *H);
extern void pnl_hmat_complex_plus_hmat(PnlHmatComplex *lhs, const PnlHmatComplex *rhs);/*lhs+=rhs*/
extern void pnl_hmat_complex_mult_dcomplex(PnlHmatComplex *lhs, dcomplex x);/* lhs *=x;*/

extern void pnl_hmat_complex_set(PnlHmatComplex *H, int *tab, dcomplex x);
extern dcomplex pnl_hmat_complex_get(const PnlHmatComplex *H, int *tab);
extern dcomplex* pnl_hmat_complex_lget(PnlHmatComplex *H, int *tab);
extern PnlMatComplex pnl_mat_complex_wrap_hmat(PnlHmatComplex *H, int *t);
extern PnlVectComplex pnl_vect_complex_wrap_hmat(PnlHmatComplex *H, int *t);

/*@}*/
/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_MATRIX_COMPLEX_H */
