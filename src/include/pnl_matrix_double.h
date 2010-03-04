#ifndef _MATRIX_DOUBLE_H
#define _MATRIX_DOUBLE_H


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl_vector.h"
#include "pnl_perm.h"
#include "pnl_matrix_int.h"


#ifndef PNL_RANGE_CHECK_OFF

#define MGET(v,i,j) pnl_mat_get((v), (i), (j))
#define MLET(v,i,j) *(pnl_mat_lget((v), (i), (j)))

#else /* PNL_RANGE_CHECK_OFF */

#define MGET(v,i,j) v->array[(i)*v->n+(j)]
#define MLET(v,i,j) v->array[(i)*v->n+(j)]

#endif /* PNL_RANGE_CHECK_OFF */

/**
 * \defgroup PnlVect Double Vector 
 */
/*@{*/
/* Compact PnlVect : used for variables that can either contain a single
   number or a PnlVect.
   vectors likes x*ones(n,1) are simply stored as a double x */
typedef struct PnlVectCompact {
  int size; /*!< size of the vector */
  union {
    double val; /*!< single value */
    double *array; /*!< Pointer to double values */
  };
  char convert; /*!< 'a', 'd' : array, double */
} PnlVectCompact;

/*@}*/


/**
 * \addtogroup PnlVect
 */
/*@{*/
 
extern PnlVectCompact* pnl_vect_compact_create (int n, double x);
extern int pnl_vect_compact_resize (PnlVectCompact *v, int size, double x);
extern PnlVectCompact* pnl_vect_compact_copy(const PnlVectCompact *v);
extern void pnl_vect_compact_free (PnlVectCompact **v);
extern PnlVect* pnl_vect_compact_to_pnl_vect (const PnlVectCompact *C);
extern double pnl_vect_compact_get (const PnlVectCompact *C, int i);

/*@}*/


/**
 * \ingroup PnlMatrices
 */
/*@{*/

/**
 * \defgroup PnlMat Double Matrix 
 */
/*@{*/

typedef struct PnlMat{
  int m; /*!< nb rows */ 
  int n; /*!< nb columns */ 
  int mn; /*!< product m*n */
  int mem_size; /*!< size of the memory block allocated for array */
  double *array; /*!< pointer to store the data row-wise */
  int owner; /*!< 1 if the owns its array pointer */
} PnlMat;

extern PnlMat* pnl_mat_create(int m, int n); 
extern PnlMat* pnl_mat_create_from_double(int m, int n, double x);
extern PnlMat* pnl_mat_create_from_ptr(int m, int n, const double* x);
extern PnlMat* pnl_mat_create_from_list(int m, int n, ...); 
extern PnlMat pnl_mat_create_wrap_array(const double* x,int m, int n);
extern PnlMat* pnl_mat_create_from_file (const char * file);
extern int pnl_mat_resize(PnlMat *v, int m, int n);
extern void pnl_mat_free(PnlMat **v);
extern PnlMat* pnl_mat_copy(const PnlMat *v);
extern void pnl_mat_clone(PnlMat *clone, const PnlMat *M);
extern void pnl_mat_map_inplace(PnlMat *lhs, double(*f)(double)); /*lhs=f(lhs)*/
extern void pnl_mat_map(PnlMat *lhs, const PnlMat *rhs, double(*f)(double));/* lhs(i)=f(rhs(i)) */
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
extern double pnl_mat_scalar_prod_A(const PnlMat *A, const PnlVect *x , const PnlVect * y);

extern void pnl_mat_set_double(PnlMat *lhs, double x);
extern void pnl_mat_set_id(PnlMat *lhs);
extern void pnl_mat_set_diag(PnlMat *lhs, double x, int d);
extern void pnl_mat_mult_mat_term(PnlMat *lhs, const PnlMat *rhs); 
extern void pnl_mat_div_mat_term(PnlMat *lhs, const PnlMat *rhs); 
extern void pnl_mat_sq_transpose(PnlMat *M);
extern PnlMat* pnl_mat_transpose(const PnlMat *M);
extern void pnl_mat_print(const PnlMat *M);
extern void pnl_mat_fprint(FILE *fic, const PnlMat *M);
extern void pnl_mat_fprint_nsp (FILE *fic, const PnlMat *M);
extern void pnl_mat_print_nsp (const PnlMat *M);
extern PnlVect pnl_mat_wrap_row(const PnlMat *M, int i);
extern PnlVect pnl_mat_wrap_vect(const PnlMat *M);
extern void pnl_mat_row_to_vect_inplace(PnlVect *V, const PnlMat *M, int i);
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
extern void pnl_mat_min (const PnlMat *A, PnlVect *out, char d);
extern void pnl_mat_max (const PnlMat *A, PnlVect *out, char d);
extern void pnl_mat_minmax (const PnlMat *A, PnlVect *m, PnlVect *M, char d);
extern void pnl_mat_min_index (const PnlMat *A, PnlVect *out, PnlVectInt *i, char d);
extern void pnl_mat_max_index (const PnlMat *A, PnlVect *out, PnlVectInt *i, char d);
extern void pnl_mat_minmax_index (const PnlMat *A, PnlVect *m, PnlVect *M, PnlVectInt *im, PnlVectInt *iM, char d);
extern void pnl_mat_qsort_index (PnlMat *A, PnlMatInt *t, char dir, char order);
extern void pnl_mat_qsort_index (PnlMat *A, char dir, char order);
extern void pnl_mat_set_row(PnlMat *M, const PnlVect *V, int i);
extern void pnl_mat_set_col(PnlMat *M, const PnlVect *V, int j);
extern void pnl_mat_swap_rows (PnlMat *M, int i, int j);
extern void pnl_mat_upper_inverse(PnlMat *A, const PnlMat *B);
extern void pnl_mat_lower_inverse (PnlMat *A, const PnlMat *B);
extern void pnl_mat_chol(PnlMat *x);
extern void pnl_mat_chol_robust(PnlMat *x);
extern void pnl_mat_upper_syslin (PnlVect *x, const PnlMat *U, const  PnlVect *b);
extern void pnl_mat_lower_syslin (PnlVect *x, const PnlMat *L, const  PnlVect *b);
extern void pnl_mat_chol_syslin (PnlVect *x, const PnlMat *chol, const  PnlVect *b);
extern void pnl_mat_chol_syslin_inplace (const PnlMat *chol, PnlVect *b);
extern void pnl_mat_lu (PnlMat *A, PnlPermutation *p);
extern void pnl_mat_lu_syslin (PnlVect *x, const PnlMat *LU,
                               const PnlPermutation *p, const PnlVect *b);
extern void pnl_mat_lu_syslin_inplace(const PnlMat *LU, const PnlPermutation *p, PnlVect *b);
extern void pnl_mat_syslin_mat (PnlMat *A,  PnlMat *B);
extern void pnl_mat_chol_syslin_mat (PnlMat *A,  PnlMat *B);
extern void pnl_mat_syslin (PnlVect *x, const PnlMat *A, const PnlVect *b);
extern void pnl_mat_syslin_inplace (PnlMat *A, PnlVect *b);
extern void pnl_mat_chol_inverse (PnlMat *inverse, const PnlMat *A);
extern void pnl_mat_inverse (PnlMat *inverse, const PnlMat *A);
extern void pnl_mat_read_matrix_from_file (PnlMat *M,const char * filename);
extern void pnl_mat_exp (PnlMat *B, const PnlMat *A);
extern void pnl_mat_log (PnlMat *B, const PnlMat *A);
extern void pnl_mat_eigen (PnlVect *v, PnlMat *P, const PnlMat *A, int with_eigenvectors);
extern void pnl_mat_qr (PnlMat *Q, PnlMat *R, PnlPermutation *p, const PnlMat *A);
extern int pnl_mat_ls_mat (const PnlMat *A, PnlMat *B);
extern int pnl_mat_ls (const PnlMat *A, PnlVect *b);

/* inline functions if you are using GCC */
#ifdef HAVE_INLINE 
extern inline
double pnl_mat_get (const PnlMat *self, int i, int j)
{
  CheckIndexMat(self,i, j);
  return self->array[i*self->n+j];
}

extern inline
double* pnl_mat_lget (PnlMat *self, int i, int j)
{
  CheckIndexMat(self,i, j);
  return &(self->array[i*self->n+j]);
}

extern inline
void pnl_mat_set (PnlMat *self, int i, int j, double x)
{
  CheckIndexMat(self,i, j);
  self->array[i*self->n+j] = x;
}
#endif

extern void pnl_mat_set(PnlMat *v, int i, int j, double x);
extern double pnl_mat_get(const PnlMat *v, int i, int j);
extern double* pnl_mat_lget(PnlMat *v, int i, int j);

/*@}*/
/*@}*/

/**
 * \ingroup PnlHMatrices
 */
/*@{*/

/**
 * \defgroup PnlHMat Double HyperMatrix 
 */
/*@{*/

typedef struct PnlHMat{
  int ndim; /*!< nb dimensions */ 
  int *dims; /*!< pointer to store the value of the ndim dimensions */ 
  int mn; /*!< product dim_1 *...*dim_ndim */
  double *array; /*!< pointer to store */
} PnlHMat;

extern PnlHMat* pnl_hmat_create(int ndim, const int *dims); 
extern PnlHMat* pnl_hmat_create_from_double(int ndim, const int *dims, double x); 
extern PnlHMat* pnl_hmat_create_from_ptr(int ndim, const int *dims, const double *x);
extern int pnl_hmat_resize(PnlHMat *v, int ndim, const int *dims);
extern void pnl_hmat_free(PnlHMat **v);
extern PnlHMat* pnl_hmat_copy(const PnlHMat *H);
extern void pnl_hmat_clone(PnlHMat *clone, const PnlHMat *H);
extern void pnl_hmat_print(const PnlHMat *H);
extern void pnl_hmat_plus_hmat(PnlHMat *lhs, const PnlHMat *rhs);/*lhs+=rhs*/
extern void pnl_hmat_mult_double(PnlHMat *lhs, double x);/* lhs *=x;*/

extern void pnl_hmat_set(PnlHMat *H, int *tab, double x);
extern double pnl_hmat_get(const PnlHMat *H, int *tab);
extern double* pnl_hmat_lget(PnlHMat *H, int *tab);

/*@}*/
/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _MATRIX_DOUBLE_H */
