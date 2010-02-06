#ifndef MATRIX_COMPLEX_H
#define MATRIX_COMPLEX_H


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl_complex.h"
#include "pnl_vector.h"



/*@}*/

/**
 * \defgroup PnlMatComplex Int matrix structure 
 */
/*@{*/

typedef struct PnlMatComplex{
  int m; /*!< nb rows */ 
  int n; /*!< nb columns */ 
  int mn; /*!< product m*n */
  int mem_size; /*!< size of the memory block allocated for array */
  dcomplex *array; /*!< pointer to store the data row-wise */
  int owner; /*!< 1 if the structure owns its array pointer */
} PnlMatComplex;

/*@}*/

/**
 * \defgroup PnlHMatComplex Int HyperMatrix structure 
 */
/*@{*/

typedef struct PnlHMatComplex{
  int ndim; /*!< nb dimensions */ 
  int *dims; /*!< pointer to store the value of the ndim dimensions */ 
  int mn; /*!< product dim_1 *...*dim_ndim */
  dcomplex *array; /*!< pointer to store */
} PnlHMatComplex;

/*@}*/

/**
 * \addtogroup PnlMatComplex
 */
/*@{*/
 
extern PnlMatComplex* pnl_mat_complex_create(int m, int n); 
extern PnlMatComplex* pnl_mat_complex_create_from_dcomplex(int m, int n, dcomplex x);
extern PnlMatComplex* pnl_mat_complex_create_from_ptr(int m, int n, const dcomplex* x);
extern PnlMatComplex* pnl_mat_complex_create_from_list(int m, int n, ...); 
extern PnlMatComplex pnl_mat_complex_create_wrap_array(const dcomplex* x, int m, int n);
extern PnlMatComplex* pnl_mat_complex_create_from_file (const char * file);
extern int pnl_mat_complex_resize(PnlMatComplex *v, int m, int n);
extern void pnl_mat_complex_free(PnlMatComplex **v);
extern PnlMatComplex* pnl_mat_complex_copy(const PnlMatComplex *v);
extern void pnl_mat_complex_clone(PnlMatComplex *clone, const PnlMatComplex *M);
extern void pnl_mat_complex_set_dcomplex(PnlMatComplex *lhs, dcomplex x);/*lhs=x*/
extern void pnl_mat_complex_set_id(PnlMatComplex *lhs);
extern void pnl_mat_complex_set_diag(PnlMatComplex *lhs, dcomplex x, int d);
extern PnlMatComplex* pnl_mat_complex_transpose(const PnlMatComplex *M);
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
extern PnlVectComplex pnl_mat_complex_wrap_row(const PnlMatComplex *M, int i);/* M(i,:)=V(:) */
extern PnlVectComplex pnl_mat_complex_wrap_vect(const PnlMatComplex *M);
extern void pnl_mat_complex_row_to_vect_inplace(PnlVectComplex * V,const PnlMatComplex *M, int i);/* M(i,:)=V(:) */
extern void pnl_mat_complex_map_inplace(PnlMatComplex *lhs, dcomplex(*f)(dcomplex)); /*lhs=f(lhs)*/
extern void pnl_mat_complex_map(PnlMatComplex *lhs, const PnlMatComplex *rhs, dcomplex(*f)(dcomplex));/* lhs(i)=f(rhs(i)) */
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
/*@}*/


/**
 * \addtogroup PnlHMatComplex
 */
/*@{*/
extern PnlHMatComplex* pnl_hmat_complex_create(int ndim, const int *dims); 
extern PnlHMatComplex* pnl_hmat_complex_create_from_dcomplex(int ndim, const int *dims, dcomplex x); 
extern PnlHMatComplex* pnl_hmat_complex_create_from_ptr(int ndim, const int *dims, const dcomplex *x);
extern int pnl_hmat_complex_resize(PnlHMatComplex *v, int ndim, const int *dims);
extern void pnl_hmat_complex_free(PnlHMatComplex **v);
extern PnlHMatComplex* pnl_hmat_complex_copy(const PnlHMatComplex *H);
extern void pnl_hmat_complex_clone(PnlHMatComplex *clone, const PnlHMatComplex *H);
extern void pnl_hmat_complex_plus_hmat(PnlHMatComplex *lhs, const PnlHMatComplex *rhs);/*lhs+=rhs*/
extern void pnl_hmat_complex_mult_dcomplex(PnlHMatComplex *lhs, dcomplex x);/* lhs *=x;*/
/*@}*/

/* inline functions if you are using GCC */
#ifdef HAVE_INLINE 
/**
 * \addtogroup PnlMatComplex
 */
/*@{*/
extern inline
dcomplex pnl_mat_complex_get (const PnlMatComplex *self, int i, int j)
{
  CheckIndexMat(self,i, j);
  return self->array[i*self->n+j];
}

extern inline
dcomplex* pnl_mat_complex_lget (PnlMatComplex *self, int i, int j)
{
  CheckIndexMat(self,i, j);
  return &(self->array[i*self->n+j]);
}

extern inline
void pnl_mat_complex_set (PnlMatComplex *self, int i, int j, dcomplex x)
{
  CheckIndexMat(self,i, j);
  self->array[i*self->n+j] = x;
}

/*@}*/
#endif



/**
 * \addtogroup PnlMatComplex
 */
/*@{*/

extern void pnl_mat_complex_set(PnlMatComplex *v, int i, int j, dcomplex x);
extern dcomplex pnl_mat_complex_get(const PnlMatComplex *v, int i, int j);
extern dcomplex* pnl_mat_complex_lget(PnlMatComplex *v, int i, int j);

/*@}*/

/**
 * \addtogroup PnlHMatComplex
 */
/*@{*/

extern void pnl_hmat_complex_set(PnlHMatComplex *H, int *tab, dcomplex x);
extern dcomplex pnl_hmat_complex_get(const PnlHMatComplex *H, int *tab);
extern dcomplex* pnl_hmat_complex_lget(PnlHMatComplex *H, int *tab);

/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* MATRIX_COMPLEX_H */
