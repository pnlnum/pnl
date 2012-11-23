#ifndef MATRIX_UINT_H
#define MATRIX_UINT_H


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl/pnl_vector_uint.h"
#include "pnl/pnl_matrix.h"


/**
 * \ingroup PnlMatrices
 */
/*@{*/

/**
 * \defgroup PnlMatUint Uint Matrix 
 */
/*@{*/

typedef struct PnlMatUint{
  int m; /*!< nb rows */ 
  int n; /*!< nb columns */ 
  int mn; /*!< product m*n */
  int mem_size; /*!< size of the memory block allocated for array */
  uint *array; /*!< pointer to store the data row-wise */
  int owner; /*!< 1 if the owns its array pointer */
} PnlMatUint;
 
extern PnlMatUint* pnl_mat_uint_create(int m, int n); 
extern PnlMatUint* pnl_mat_uint_create_from_uint(int m, int n, uint x);
extern PnlMatUint* pnl_mat_uint_create_from_ptr(int m, int n, const uint* x);
extern PnlMatUint* pnl_mat_uint_create_from_list(int m, int n, ...); 
extern PnlMatUint pnl_mat_uint_wrap_array(const uint* x,int m, int n);
extern PnlMatUint* pnl_mat_uint_create_from_file (const char * file);
extern int pnl_mat_uint_resize(PnlMatUint *v, int m, int n);
extern void pnl_mat_uint_free(PnlMatUint **v);
extern PnlMatUint* pnl_mat_uint_copy(const PnlMatUint *v);
extern void pnl_mat_uint_clone(PnlMatUint *clone, const PnlMatUint *M);
extern void pnl_mat_uint_set_uint(PnlMatUint *lhs, uint x);/*lhs=x*/
extern void pnl_mat_uint_set_id(PnlMatUint *lhs);
extern void pnl_mat_uint_set_diag(PnlMatUint *lhs, uint x, int d);
extern PnlMatUint* pnl_mat_uint_transpose(const PnlMatUint *M);
extern void pnl_mat_uint_sq_transpose(PnlMatUint *M);
extern void pnl_mat_uint_print(const PnlMatUint *M);
extern void pnl_mat_uint_print_nsp(const PnlMatUint *M);
extern void pnl_mat_uint_fprint(FILE *fic,const PnlMatUint *M);
extern void pnl_mat_uint_fprint_nsp(FILE *fic,const PnlMatUint *M);
extern PnlMatUint* pnl_mat_uint_create_diag_from_ptr(const uint x[], int d);
extern PnlMatUint* pnl_mat_uint_create_diag(const PnlVectUint *V);
extern void
pnl_mat_uint_set_row(PnlMatUint *M, const PnlVectUint *V, int i);/* M(i,:)=V(:) */
extern void
pnl_mat_uint_set_col(PnlMatUint *M, const PnlVectUint *V, int j);/* M(:,j)=V(:) */
extern void
pnl_mat_uint_swap_rows (PnlMatUint *M, int i, int j);
extern void
pnl_mat_uint_get_row(PnlVectUint *V, const PnlMatUint *M, int i);/* V(:)=M(i,:) */
extern void pnl_mat_uint_get_col(PnlVectUint *V, const PnlMatUint *M, int j);
extern PnlVectUint pnl_vect_uint_wrap_mat_row(const PnlMatUint *M, int i);/* M(i,:)=V(:) */
extern PnlVectUint pnl_mat_uint_wrap_vect(const PnlMatUint *M);
extern void pnl_mat_uint_map_inplace(PnlMatUint *lhs, uint(*f)(uint)); /*lhs=f(lhs)*/
extern void pnl_mat_uint_map(PnlMatUint *lhs, const PnlMatUint *rhs, uint(*f)(uint));/* lhs(i)=f(rhs(i)) */
extern void pnl_mat_uint_map_mat(PnlMatUint *lhs, const PnlMatUint *rhs, uint(*f)(uint,uint));
extern void pnl_mat_uint_plus_uint(PnlMatUint *lhs, uint x); /*lhs+=x*/
extern void pnl_mat_uint_minus_uint(PnlMatUint *lhs, uint x); /*lhs-=x*/
extern void pnl_mat_uint_plus_mat(PnlMatUint *lhs, const PnlMatUint *rhs); /*lhs+=rhs*/
extern void pnl_mat_uint_minus_mat(PnlMatUint *lhs, const PnlMatUint *rhs); /*lhs-=rhs*/
extern void pnl_mat_uint_mult_uint(PnlMatUint *lhs, uint x); /*lhs*=x*/
extern void pnl_mat_uint_div_uint(PnlMatUint *lhs, uint x); /*lhs/=x*/
extern void
pnl_mat_uint_mult_mat_term(PnlMatUint *lhs, const PnlMatUint *rhs); /*lhs=lhs.*rhs*/
extern void
pnl_mat_uint_div_mat_term(PnlMatUint *lhs, const PnlMatUint *rhs); /*lhs=lhs./rhs*/
extern void pnl_mat_uint_axpy (uint a, const PnlMatUint *X, PnlMatUint *Y); 
extern void pnl_mat_uint_dger (uint alpha, const PnlVectUint *x, const PnlVectUint *y, PnlMatUint *A);
extern PnlVectUint* pnl_mat_uint_mult_vect(const PnlMatUint *mat, const PnlVectUint *vec);
extern void pnl_mat_uint_mult_vect_inplace(PnlVectUint *lhs, const PnlMatUint *mat,
                                          const PnlVectUint *rhs);

extern PnlVectUint* pnl_mat_uint_mult_vect_transpose(const PnlMatUint *mat, const PnlVectUint *vec);
extern void pnl_mat_uint_mult_vect_transpose_inplace(PnlVectUint *lhs, const PnlMatUint *mat, const PnlVectUint *rhs);
extern  void pnl_mat_uint_lAxpby(uint l, const PnlMatUint *A, const PnlVectUint *x, uint b, PnlVectUint * y);
extern void pnl_mat_uint_dgemv (char trans, uint alpha, const PnlMatUint *A,
                               const PnlVectUint *x , uint beta, PnlVectUint * y);
extern PnlMatUint* pnl_mat_uint_mult_mat(const PnlMatUint *rhs1, const PnlMatUint *rhs2);
extern void pnl_mat_uint_mult_mat_inplace(PnlMatUint *lhs, const PnlMatUint *rhs1,
                                          const PnlMatUint *rhs2);/*lhs=rhs1*rhs2*/
extern void pnl_mat_uint_dgemm (char transA, char transB, uint alpha, const PnlMatUint *A,
                                const PnlMatUint *B, uint beta, PnlMatUint *C); 
extern uint pnl_mat_uint_sum(const PnlMatUint *lhs);
extern uint pnl_mat_uint_prod(const PnlMatUint *lhs);
extern void pnl_mat_uint_sum_vect (PnlVectUint *y, const PnlMatUint *M, char a);
extern void pnl_mat_uint_cumsum(PnlMatUint *lhs, char a);
extern void pnl_mat_uint_prod_vect (PnlVectUint *y, const PnlMatUint *M, char a);
extern void pnl_mat_uint_cumprod(PnlMatUint *lhs, char a);
extern void pnl_mat_uint_min (const PnlMatUint *A, PnlVectUint *out, char d);
extern void pnl_mat_uint_max (const PnlMatUint *A, PnlVectUint *out, char d);
extern void pnl_mat_uint_minmax (const PnlMatUint *A, PnlVectUint *m, PnlVectUint *M, char d);
extern void pnl_mat_uint_min_index (const PnlMatUint *A, PnlVectUint *out, PnlVectInt *i, char d);
extern void pnl_mat_uint_max_index (const PnlMatUint *A, PnlVectUint *out, PnlVectInt *i, char d);
extern void pnl_mat_uint_minmax_index (const PnlMatUint *A, PnlVectUint *m, PnlVectUint *M, PnlVectInt *im, PnlVectInt *iM, char d);
extern void pnl_mat_uint_qsort_index (PnlMatUint *A, PnlMatInt *t, char dir, char order);
extern void pnl_mat_uint_qsort (PnlMatUint *A, char dir, char order);
extern void pnl_vect_uint_reverse(PnlVectUint * v);


/* inline functions if you are using GCC */
#ifdef PNL_HAVE_INLINE 
extern inline
uint pnl_mat_uint_get (const PnlMatUint *self, int i, int j)
{
  CheckIndexMat(self,i, j);
  return self->array[i*self->n+j];
}

extern inline
uint* pnl_mat_uint_lget (PnlMatUint *self, int i, int j)
{
  CheckIndexMat(self,i, j);
  return &(self->array[i*self->n+j]);
}

extern inline
void pnl_mat_uint_set (PnlMatUint *self, int i, int j, uint x)
{
  CheckIndexMat(self,i, j);
  self->array[i*self->n+j] = x;
}

#endif

extern void pnl_mat_uint_set(PnlMatUint *v, int i, int j, uint x);
extern uint pnl_mat_uint_get(const PnlMatUint *v, int i, int j);
extern uint* pnl_mat_uint_lget(PnlMatUint *v, int i, int j);

/*@}*/
/*@}*/


/**
 * \ingroup PnlHmatrices a Hyper Matrix
 */
/*@{*/
/**
 * \defgroup PnlHmatUint Uint HyperMatrix 
 */
/*@{*/

typedef struct PnlHmatUint{
  int ndim; /*!< nb dimensions */ 
  int *dims; /*!< pointer to store the value of the ndim dimensions */ 
  int mn; /*!< product dim_1 *...*dim_ndim */
  uint *array; /*!< pointer to store */
} PnlHmatUint;

extern PnlHmatUint* pnl_hmat_uint_create(int ndim, const int *dims); 
extern PnlHmatUint* pnl_hmat_uint_create_from_uint(int ndim, const int *dims, uint x); 
extern PnlHmatUint* pnl_hmat_uint_create_from_ptr(int ndim, const int *dims, const uint *x);
extern int pnl_hmat_uint_resize(PnlHmatUint *v, int ndim, const int *dims);
extern void pnl_hmat_uint_free(PnlHmatUint **v);
extern void pnl_hmat_uint_print(const PnlHmatUint *H);
extern PnlHmatUint* pnl_hmat_uint_copy(const PnlHmatUint *H);
extern void pnl_hmat_uint_clone(PnlHmatUint *clone, const PnlHmatUint *H);
extern void pnl_hmat_uint_plus_hmat(PnlHmatUint *lhs, const PnlHmatUint *rhs);/*lhs+=rhs*/
extern void pnl_hmat_uint_mult_uint(PnlHmatUint *lhs, uint x);/* lhs *=x;*/


extern void pnl_hmat_uint_set(PnlHmatUint *H, int *tab, uint x);
extern uint pnl_hmat_uint_get(const PnlHmatUint *H, int *tab);
extern uint* pnl_hmat_uint_lget(PnlHmatUint *H, int *tab);

/*@}*/
/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* MATRIX_UINT_H */
