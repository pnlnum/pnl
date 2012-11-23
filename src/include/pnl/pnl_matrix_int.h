#ifndef _PNL_MATRIX_INT_H
#define _PNL_MATRIX_INT_H

#ifndef _PNL_MATRIX_H
#error "Do not include this file directly. Include pnl_matrix.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */



/**
 * \ingroup PnlMatrices
 */
/*@{*/

/**
 * \defgroup PnlMatInt Int Matrix 
 */
/*@{*/
  


extern void pnl_mat_int_init(PnlMatInt *); 
extern PnlMatInt* pnl_mat_int_new(); 
extern int pnl_mat_int_eq (const PnlMatInt *, const PnlMatInt *);
extern int pnl_mat_int_eq_int (const PnlMatInt *, int);
extern PnlMatInt* pnl_mat_int_create(int m, int n); 
extern PnlMatInt* pnl_mat_int_create_from_int(int m, int n, int x);
extern PnlMatInt* pnl_mat_int_create_from_ptr(int m, int n, const int* x);
extern PnlMatInt* pnl_mat_int_create_from_list(int m, int n, ...); 
extern PnlMatInt pnl_mat_int_wrap_array(const int* x,int m, int n);
extern PnlMatInt* pnl_mat_int_create_from_file (const char * file);
extern void pnl_vect_int_extract_submat (PnlVectInt *sub, const PnlMatInt *M, const PnlVectInt *indi, const PnlVectInt *indj);
extern PnlVectInt* pnl_vect_int_create_submat (const PnlMatInt *M, const PnlVectInt *indi, const PnlVectInt *indj);
extern void pnl_mat_int_extract_subblock (PnlMatInt *M_sub, const PnlMatInt *M, int i, int len_i, int j, int len_j);
extern int pnl_mat_int_resize(PnlMatInt *v, int m, int n);
extern void pnl_mat_int_free(PnlMatInt **v);
extern PnlMatInt* pnl_mat_int_copy(const PnlMatInt *v);
extern void pnl_mat_int_clone(PnlMatInt *clone, const PnlMatInt *M);
extern void pnl_mat_int_set_int(PnlMatInt *lhs, int x);/*lhs=x*/
extern void pnl_mat_int_set_zero(PnlMatInt *lhs);
extern void pnl_mat_int_set_id(PnlMatInt *lhs);
extern void pnl_mat_int_set_diag(PnlMatInt *lhs, int x, int d);
extern void pnl_mat_int_set_from_ptr(PnlMatInt *lhs, const int *x);
extern PnlMatInt* pnl_mat_int_transpose(const PnlMatInt *M);
extern void pnl_mat_int_tr(PnlMatInt *tM, const PnlMatInt *M);
extern void pnl_mat_int_sq_transpose(PnlMatInt *M);
extern void pnl_mat_int_print(const PnlMatInt *M);
extern void pnl_mat_int_print_nsp(const PnlMatInt *M);
extern void pnl_mat_int_fprint(FILE *fic,const PnlMatInt *M);
extern void pnl_mat_int_fprint_nsp(FILE *fic,const PnlMatInt *M);
extern PnlMatInt* pnl_mat_int_create_diag_from_ptr(const int x[], int d);
extern PnlMatInt* pnl_mat_int_create_diag(const PnlVectInt *V);
extern void pnl_mat_int_set_row(PnlMatInt *M, const PnlVectInt *V, int i);/* M(i,:)=V(:) */
extern void pnl_mat_int_set_col(PnlMatInt *M, const PnlVectInt *V, int j);/* M(:,j)=V(:) */
extern void pnl_mat_int_swap_rows (PnlMatInt *M, int i, int j);
extern void pnl_mat_int_get_row(PnlVectInt *V, const PnlMatInt *M, int i);/* V(:)=M(i,:) */
extern void pnl_mat_int_get_col(PnlVectInt *V, const PnlMatInt *M, int j);
extern PnlVectInt pnl_vect_int_wrap_mat_row(const PnlMatInt *M, int i);/* M(i,:)=V(:) */
extern PnlMatInt pnl_mat_int_wrap_vect(const PnlVectInt *V);
extern void pnl_mat_int_map_inplace(PnlMatInt *lhs, int(*f)(int)); /*lhs=f(lhs)*/
extern void pnl_mat_int_map(PnlMatInt *lhs, const PnlMatInt *rhs, int(*f)(int));/* lhs(i)=f(rhs(i)) */
extern void pnl_mat_int_map_mat_inplace(PnlMatInt *lhs, const PnlMatInt *rhs, int(*f)(int,int));
extern void pnl_mat_int_map_mat(PnlMatInt *lhs, const PnlMatInt *rhs1, const PnlMatInt *rhs2, int(*f)(int,int));
extern int pnl_mat_int_find(PnlVectInt *indi, PnlVectInt *indj, char* type, int(*f)(int *), ...);
extern void pnl_mat_int_plus_int(PnlMatInt *lhs, int x); /*lhs+=x*/
extern void pnl_mat_int_minus_int(PnlMatInt *lhs, int x); /*lhs-=x*/
extern void pnl_mat_int_plus_mat(PnlMatInt *lhs, const PnlMatInt *rhs); /*lhs+=rhs*/
extern void pnl_mat_int_minus_mat(PnlMatInt *lhs, const PnlMatInt *rhs); /*lhs-=rhs*/
extern void pnl_mat_int_mult_int(PnlMatInt *lhs, int x); /*lhs*=x*/
extern void pnl_mat_int_div_int(PnlMatInt *lhs, int x); /*lhs/=x*/
extern void
pnl_mat_int_mult_mat_term(PnlMatInt *lhs, const PnlMatInt *rhs); /*lhs=lhs.*rhs*/
extern void
pnl_mat_int_div_mat_term(PnlMatInt *lhs, const PnlMatInt *rhs); /*lhs=lhs./rhs*/
extern void pnl_mat_int_axpy (int a, const PnlMatInt *X, PnlMatInt *Y); 
extern void pnl_mat_int_dger (int alpha, const PnlVectInt *x, const PnlVectInt *y, PnlMatInt *A);
extern PnlVectInt* pnl_mat_int_mult_vect(const PnlMatInt *mat, const PnlVectInt *vec);
extern void pnl_mat_int_mult_vect_inplace(PnlVectInt *lhs, const PnlMatInt *mat,
                                          const PnlVectInt *rhs);

extern PnlVectInt* pnl_mat_int_mult_vect_transpose(const PnlMatInt *mat, const PnlVectInt *vec);
extern void pnl_mat_int_mult_vect_transpose_inplace(PnlVectInt *lhs, const PnlMatInt *mat, const PnlVectInt *rhs);
extern  void pnl_mat_int_lAxpby(int l, const PnlMatInt *A, const PnlVectInt *x, int b, PnlVectInt * y);
extern void pnl_mat_int_dgemv (char trans, int alpha, const PnlMatInt *A,
                               const PnlVectInt *x , int beta, PnlVectInt * y);
extern PnlMatInt* pnl_mat_int_mult_mat(const PnlMatInt *rhs1, const PnlMatInt *rhs2);
extern void pnl_mat_int_mult_mat_inplace(PnlMatInt *lhs, const PnlMatInt *rhs1,
                                     const PnlMatInt *rhs2);/*lhs=rhs1*rhs2*/
extern void pnl_mat_int_dgemm (char transA, char transB, int alpha, const PnlMatInt *A,
                           const PnlMatInt *B, int beta, PnlMatInt *C); 
extern int pnl_mat_int_sum(const PnlMatInt *lhs);
extern int pnl_mat_int_prod(const PnlMatInt *lhs);
extern void pnl_mat_int_sum_vect (PnlVectInt *y, const PnlMatInt *M, char a);
extern void pnl_mat_int_cumsum(PnlMatInt *lhs, char a);
extern void pnl_mat_int_prod_vect (PnlVectInt *y, const PnlMatInt *M, char a);
extern void pnl_mat_int_cumprod(PnlMatInt *lhs, char a);
extern void pnl_mat_int_min (PnlVectInt *out, const PnlMatInt *A, char d);
extern void pnl_mat_int_max (PnlVectInt *out, const PnlMatInt *A, char d);
extern void pnl_mat_int_minmax (PnlVectInt *m, PnlVectInt *M, const PnlMatInt *A, char d);
extern void pnl_mat_int_min_index (PnlVectInt *out, PnlVectInt *i, const PnlMatInt *A, char d);
extern void pnl_mat_int_max_index (PnlVectInt *out, PnlVectInt *i, const PnlMatInt *A, char d);
extern void pnl_mat_int_minmax_index (PnlVectInt *m, PnlVectInt *M, PnlVectInt *im, PnlVectInt *iM, const PnlMatInt *A, char d);
extern void pnl_mat_int_qsort_index (PnlMatInt *A, PnlMatInt *t, char dir, char order);
extern void pnl_mat_int_qsort (PnlMatInt *A, char dir, char order);


/* inline functions if you are using GCC */
#ifdef PNL_PNL_HAVE_INLINE 
PNL_INLINE_FUNC int pnl_mat_int_get (const PnlMatInt *self, int i, int j)
{
  CheckIndexMat(self,i, j);
  return self->array[i*self->n+j];
}

PNL_INLINE_FUNC int* pnl_mat_int_lget (PnlMatInt *self, int i, int j)
{
  CheckIndexMat(self,i, j);
  return &(self->array[i*self->n+j]);
}

PNL_INLINE_FUNC void pnl_mat_int_set (PnlMatInt *self, int i, int j, int x)
{
  CheckIndexMat(self,i, j);
  self->array[i*self->n+j] = x;
}

#endif

PNL_INLINE_DECL void pnl_mat_int_set(PnlMatInt *v, int i, int j, int x);
PNL_INLINE_DECL int pnl_mat_int_get(const PnlMatInt *v, int i, int j);
PNL_INLINE_DECL int* pnl_mat_int_lget(PnlMatInt *v, int i, int j);

/*@}*/
/*@}*/


/**
 * \ingroup PnlHmatrices a Hyper Matrix
 */
/*@{*/
/**
 * \defgroup PnlHmatInt Int HyperMatrix 
 */
/*@{*/



extern PnlHmatInt* pnl_hmat_int_new(); 
extern PnlHmatInt* pnl_hmat_int_create(int ndim, const int *dims); 
extern PnlHmatInt* pnl_hmat_int_create_from_int(int ndim, const int *dims, int x); 
extern PnlHmatInt* pnl_hmat_int_create_from_ptr(int ndim, const int *dims, const int *x);
extern int pnl_hmat_int_resize(PnlHmatInt *v, int ndim, const int *dims);
extern void pnl_hmat_int_free(PnlHmatInt **v);
extern PnlHmatInt* pnl_hmat_int_copy(const PnlHmatInt *H);
extern void pnl_hmat_int_print(const PnlHmatInt *H);
extern void pnl_hmat_int_clone(PnlHmatInt *clone, const PnlHmatInt *H);
extern void pnl_hmat_int_plus_hmat(PnlHmatInt *lhs, const PnlHmatInt *rhs);/*lhs+=rhs*/
extern void pnl_hmat_int_mult_int(PnlHmatInt *lhs, int x);/* lhs *=x;*/

extern void pnl_hmat_int_set(PnlHmatInt *H, int *tab, int x);
extern int pnl_hmat_int_get(const PnlHmatInt *H, int *tab);
extern int* pnl_hmat_int_lget(PnlHmatInt *H, int *tab);
extern PnlMatInt pnl_mat_int_wrap_hmat(PnlHmatInt *H, int *t);
extern PnlVectInt pnl_vect_int_wrap_hmat(PnlHmatInt *H, int *t);

/*@}*/
/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_MATRIX_INT_H */
