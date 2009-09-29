#ifndef MATRIX_UINT_H
#define MATRIX_UINT_H

#include "pnl_vector_uint.h"
#include "pnl_matrix.h"



/*@}*/

/**
 * \defgroup PnlMatUint Uint matrix structure for Premia
 */
/*@{*/

typedef struct PnlMatUint{
  int m; /*!< nb rows */ 
  int n; /*!< nb columns */ 
  int mn; /*!< product m*n */
  int mem_size; /*!< size of the memory block allocated for array */
  uint *array; /*!< pointer to store the data row-wise */
  int owner; /*!< 1 if the structure owns its array pointer */
} PnlMatUint;

/*@}*/

/**
 * \defgroup PnlHMatUint Uint HyperMatrix structure for Premia
 */
/*@{*/

typedef struct PnlHMatUint{
  int ndim; /*!< nb dimensions */ 
  int *dims; /*!< pointer to store the value of the ndim dimensions */ 
  int mn; /*!< product dim_1 *...*dim_ndim */
  uint *array; /*!< pointer to store */
} PnlHMatUint;

/*@}*/

/**
 * \addtogroup PnlMatUint
 */
/*@{*/
 
extern PnlMatUint* pnl_mat_uint_create(int m, int n); 
extern PnlMatUint* pnl_mat_uint_create_from_uint(int m, int n, uint x);
extern PnlMatUint* pnl_mat_uint_create_from_ptr(int m, int n, const uint* x);
extern PnlMatUint* pnl_mat_uint_create_from_list(int m, int n, ...); 
extern PnlMatUint pnl_mat_uint_create_wrap_array(const uint* x,int m, int n);
extern PnlMatUint* pnl_mat_uint_create_from_file (const char * file);
extern int pnl_mat_uint_resize(PnlMatUint *v, int m, int n);
extern void pnl_mat_uint_free(PnlMatUint **v);
extern PnlMatUint* pnl_mat_uint_copy(const PnlMatUint *v);
extern void pnl_mat_uint_clone(PnlMatUint *clone, const PnlMatUint *M);
extern void pnl_mat_uint_set_uint(PnlMatUint *lhs, uint x);/*lhs=x*/
extern void pnl_mat_uint_set_id(PnlMatUint *lhs);
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
extern PnlVectUint pnl_mat_uint_wrap_row(const PnlMatUint *M, int i);/* M(i,:)=V(:) */
extern PnlVectUint pnl_mat_uint_wrap_vect(const PnlMatUint *M);
extern void pnl_mat_uint_row_to_vect_inplace(PnlVectUint * V,const PnlMatUint *M, int i);/* M(i,:)=V(:) */
extern void pnl_mat_uint_map_inplace(PnlMatUint *lhs, uint(*f)(uint)); /*lhs=f(lhs)*/
extern void pnl_mat_uint_map(PnlMatUint *lhs, const PnlMatUint *rhs, uint(*f)(uint));/* lhs(i)=f(rhs(i)) */
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
/*@}*/


/**
 * \addtogroup PnlHMatUint
 */
/*@{*/
extern PnlHMatUint* pnl_hmat_uint_create(int ndim, const int *dims); 
extern PnlHMatUint* pnl_hmat_uint_create_from_uint(int ndim, const int *dims, uint x); 
extern PnlHMatUint* pnl_hmat_uint_create_from_ptr(int ndim, const int *dims, const uint *x);
extern int pnl_hmat_uint_resize(PnlHMatUint *v, int ndim, const int *dims);
extern void pnl_hmat_uint_free(PnlHMatUint **v);
extern PnlHMatUint* pnl_hmat_uint_copy(const PnlHMatUint *H);
extern void pnl_hmat_uint_clone(PnlHMatUint *clone, const PnlHMatUint *H);
extern void pnl_hmat_uint_plus_hmat(PnlHMatUint *lhs, const PnlHMatUint *rhs);/*lhs+=rhs*/
extern void pnl_hmat_uint_mult_uint(PnlHMatUint *lhs, uint x);/* lhs *=x;*/
/*@}*/

/* inline functions if you are using GCC */
#ifdef HAVE_INLINE 
/**
 * \addtogroup PnlMatUint
 */
/*@{*/
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

/*@}*/
#endif



/**
 * \addtogroup PnlMatUint
 */
/*@{*/

extern void pnl_mat_uint_set(PnlMatUint *v, int i, int j, uint x);
extern uint pnl_mat_uint_get(const PnlMatUint *v, int i, int j);
extern uint* pnl_mat_uint_lget(PnlMatUint *v, int i, int j);

/*@}*/

/**
 * \addtogroup PnlHMatUint
 */
/*@{*/

extern void pnl_hmat_uint_set(PnlHMatUint *H, int *tab, uint x);
extern uint pnl_hmat_uint_get(const PnlHMatUint *H, int *tab);
extern uint* pnl_hmat_uint_lget(PnlHMatUint *H, int *tab);

/*@}*/

#endif /* MATRIX_UINT_H */


