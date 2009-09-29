#ifndef MATRIX_INT_H
#define MATRIX_INT_H

#include "pnl_vector_int.h"
#include "pnl_matrix.h"


/**
 * \defgroup PnlMatInt Int matrix structure for Premia
 */
/*@{*/

typedef struct PnlMatInt{
  int m; /*!< nb rows */ 
  int n; /*!< nb columns */ 
  int mn; /*!< product m*n */
  int mem_size; /*!< size of the memory block allocated for array */
  int *array; /*!< pointer to store the data row-wise */
  int owner; /*!< 1 if the structure owns its array pointer */
} PnlMatInt;

/*@}*/

/**
 * \defgroup PnlHMatInt Int HyperMatrix structure for Premia
 */
/*@{*/

typedef struct PnlHMatInt{
  int ndim; /*!< nb dimensions */ 
  int *dims; /*!< pointer to store the value of the ndim dimensions */ 
  int mn; /*!< product dim_1 *...*dim_ndim */
  int *array; /*!< pointer to store */
} PnlHMatInt;

/*@}*/

/**
 * \addtogroup PnlMatInt
 */
/*@{*/
 
extern PnlMatInt* pnl_mat_int_create(int m, int n); 
extern PnlMatInt* pnl_mat_int_create_from_int(int m, int n, int x);
extern PnlMatInt* pnl_mat_int_create_from_ptr(int m, int n, const int* x);
extern PnlMatInt* pnl_mat_int_create_from_list(int m, int n, ...); 
extern PnlMatInt pnl_mat_int_create_wrap_array(const int* x,int m, int n);
extern PnlMatInt* pnl_mat_int_create_from_file (const char * file);
extern int pnl_mat_int_resize(PnlMatInt *v, int m, int n);
extern void pnl_mat_int_free(PnlMatInt **v);
extern PnlMatInt* pnl_mat_int_copy(const PnlMatInt *v);
extern void pnl_mat_int_clone(PnlMatInt *clone, const PnlMatInt *M);
extern void pnl_mat_int_set_int(PnlMatInt *lhs, int x);/*lhs=x*/
extern void pnl_mat_int_set_id(PnlMatInt *lhs);
extern PnlMatInt* pnl_mat_int_transpose(const PnlMatInt *M);
extern void pnl_mat_int_sq_transpose(PnlMatInt *M);
extern void pnl_mat_int_print(const PnlMatInt *M);
extern void pnl_mat_int_print_nsp(const PnlMatInt *M);
extern void pnl_mat_int_fprint(FILE *fic,const PnlMatInt *M);
extern void pnl_mat_int_fprint_nsp(FILE *fic,const PnlMatInt *M);
extern PnlMatInt* pnl_mat_int_create_diag_from_ptr(const int x[], int d);
extern PnlMatInt* pnl_mat_int_create_diag(const PnlVectInt *V);
extern void
pnl_mat_int_set_row(PnlMatInt *M, const PnlVectInt *V, int i);/* M(i,:)=V(:) */
extern void
pnl_mat_int_set_col(PnlMatInt *M, const PnlVectInt *V, int j);/* M(:,j)=V(:) */
extern void
pnl_mat_int_swap_rows (PnlMatInt *M, int i, int j);
extern void
pnl_mat_int_get_row(PnlVectInt *V, const PnlMatInt *M, int i);/* V(:)=M(i,:) */
extern void pnl_mat_int_get_col(PnlVectInt *V, const PnlMatInt *M, int j);
extern PnlVectInt pnl_mat_int_wrap_row(const PnlMatInt *M, int i);/* M(i,:)=V(:) */
extern PnlVectInt pnl_mat_int_wrap_vect(const PnlMatInt *M);
extern void pnl_mat_int_row_to_vect_inplace(PnlVectInt * V,const PnlMatInt *M, int i);/* M(i,:)=V(:) */
extern void pnl_mat_int_map_inplace(PnlMatInt *lhs, int(*f)(int)); /*lhs=f(lhs)*/
extern void pnl_mat_int_map(PnlMatInt *lhs, const PnlMatInt *rhs, int(*f)(int));/* lhs(i)=f(rhs(i)) */
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
extern int pnl_mat_int_sum(const PnlMatInt *lhs);
extern int pnl_mat_int_prod(const PnlMatInt *lhs);
extern void pnl_mat_int_sum_vect (PnlVectInt *y, const PnlMatInt *M, char a);
extern void pnl_mat_int_cumsum(PnlMatInt *lhs, char a);
extern void pnl_mat_int_prod_vect (PnlVectInt *y, const PnlMatInt *M, char a);
extern void pnl_mat_int_cumprod(PnlMatInt *lhs, char a);
extern void pnl_mat_int_min (const PnlMatInt *A, PnlVectInt *out, char d);
extern void pnl_mat_int_max (const PnlMatInt *A, PnlVectInt *out, char d);
extern void pnl_mat_int_minmax (const PnlMatInt *A, PnlVectInt *m, PnlVectInt *M, char d);
extern void pnl_mat_int_min_index (const PnlMatInt *A, PnlVectInt *out, PnlVectInt *i, char d);
extern void pnl_mat_int_max_index (const PnlMatInt *A, PnlVectInt *out, PnlVectInt *i, char d);
extern void pnl_mat_int_minmax_index (const PnlMatInt *A, PnlVectInt *m, PnlVectInt *M, PnlVectInt *im, PnlVectInt *iM, char d);
extern void pnl_mat_int_qsort_index (PnlMatInt *A, PnlMatInt *t, char dir, char order);
/*@}*/


/**
 * \addtogroup PnlHMatInt
 */
/*@{*/
extern PnlHMatInt* pnl_hmat_int_create(int ndim, const int *dims); 
extern PnlHMatInt* pnl_hmat_int_create_from_int(int ndim, const int *dims, int x); 
extern PnlHMatInt* pnl_hmat_int_create_from_ptr(int ndim, const int *dims, const int *x);
extern int pnl_hmat_int_resize(PnlHMatInt *v, int ndim, const int *dims);
extern void pnl_hmat_int_free(PnlHMatInt **v);
extern PnlHMatInt* pnl_hmat_int_copy(const PnlHMatInt *H);
extern void pnl_hmat_int_clone(PnlHMatInt *clone, const PnlHMatInt *H);
extern void pnl_hmat_int_plus_hmat(PnlHMatInt *lhs, const PnlHMatInt *rhs);/*lhs+=rhs*/
extern void pnl_hmat_int_mult_int(PnlHMatInt *lhs, int x);/* lhs *=x;*/
/*@}*/

/* inline functions if you are using GCC */
#ifdef HAVE_INLINE 
/**
 * \addtogroup PnlMatInt
 */
/*@{*/
extern inline
int pnl_mat_int_get (const PnlMatInt *self, int i, int j)
{
  CheckIndexMat(self,i, j);
  return self->array[i*self->n+j];
}

extern inline
int* pnl_mat_int_lget (PnlMatInt *self, int i, int j)
{
  CheckIndexMat(self,i, j);
  return &(self->array[i*self->n+j]);
}

extern inline
void pnl_mat_int_set (PnlMatInt *self, int i, int j, int x)
{
  CheckIndexMat(self,i, j);
  self->array[i*self->n+j] = x;
}

/*@}*/
#endif



/**
 * \addtogroup PnlMatInt
 */
/*@{*/

extern void pnl_mat_int_set(PnlMatInt *v, int i, int j, int x);
extern int pnl_mat_int_get(const PnlMatInt *v, int i, int j);
extern int* pnl_mat_int_lget(PnlMatInt *v, int i, int j);

/*@}*/

/**
 * \addtogroup PnlHMatInt
 */
/*@{*/

extern void pnl_hmat_int_set(PnlHMatInt *H, int *tab, int x);
extern int pnl_hmat_int_get(const PnlHMatInt *H, int *tab);
extern int* pnl_hmat_int_lget(PnlHMatInt *H, int *tab);

/*@}*/

#endif /* MATRIX_INT_H */


