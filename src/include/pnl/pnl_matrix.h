#ifndef _PNL_MATRIX_H
#define _PNL_MATRIX_H



#ifndef PNL_RANGE_CHECK_OFF
#define CheckIndexHmat(H,index) {{int l;                \
      for(l=0; l<(H)->ndim; l++)                        \
        { if((index)[l]>(H)->dims[l] || (index)[l]<0)   \
            {perror("index out of range"); abort();}    \
        }                                               \
    }}
#define CheckIndexMat(v,i,j) {                                          \
    if (i>=v->m || j>=v->n || i<0 || j<0) {perror("index out of range"); abort();}}
#define CheckIsSquare(v) {                                          \
    if (v->m != v->n) {perror("not a square matrix"); abort();}}
#define CheckMatIsCompatible(lhs, rhs) {                \
    if((lhs)->n != (rhs)->m)                            \
      {perror("non compatible dimensions"); abort();}}
#define CheckMatVectIsCompatible(mat, vect){            \
    if((mat)->n != (vect)->size)                        \
      {perror("non compatible dimensions"); abort();}}
#define CheckMatMatch(lhs, rhs) { if ((lhs)->m != (rhs)->m || (lhs)->n != (rhs)->n) \
      {perror("non compatible dimensions"); abort();}}
#define CheckHmatMatch(lhs,rhs) {                                       \
    if  ((lhs)->ndim != (rhs)->ndim) {perror("index out of range"); abort();} \
    {int l;                                                             \
      for(l=0; l<(lhs)->ndim; l++)                                      \
        { if((lhs)->dims[l] != (rhs)->dims[l])                          \
            {perror("index out of range"); abort();}                    \
        }                                                               \
    }}

#else

#define CheckIndexHmat(H,index) {}
#define CheckIndexMat(v,i,j) {}
#define CheckIsSquare(v) {}
#define CheckMatIsCompatible(lhs, rhs) {}
#define CheckMatVectIsCompatible(mat, vect){}
#define CheckMatMatch(lhs, rhs) {}
#define CheckHmatMatch(lhs,rhs) {}

#endif /* PNL_RANGE_CHECK_OFF */

#include "pnl/pnl_config.h"

/**
 * \defgroup PnlMat a Matrix object
 *
 * Matrix are stored row-wise in a one dimensional array.
 * The element (i,j) of a matrix is stored in array[i*n+j]
 */
/*@{*/
#define PNL_MGET(v,i,j) (v)->array[(i)*(v)->n+(j)]
#define PNL_MSET(v,i,j, x) (v)->array[(i)*(v)->n+(j)] = (x)
#define PNL_MLET(v,i,j) (v)->array[(i)*(v)->n+(j)]

#include "pnl/pnl_matvect.h"
#include "pnl/pnl_perm.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

extern PnlMatObject* pnl_mat_object_new ();
extern void pnl_mat_object_free (PnlMatObject **);
extern int pnl_mat_object_resize(PnlMatObject *M, int m, int n);
/*@}*/

/**
 * \defgroup PnlHmat Hyper Matrix object
 *
 * HMatrices are stored as a contiguous memory block following the same
 * scheme as for matrices and applying it recursively. For instance, for a
 * three dimensional Hmatrix with size n1 x n2 x n3, the element (i,j,k) is
 * located at array[k+n1*j+n2*n3*i]
 */
/*@{*/
extern PnlHmatObject* pnl_hmat_object_new ();
extern void pnl_hmat_object_free (PnlHmatObject **);
extern int pnl_hmat_object_resize(PnlHmatObject *H, int ndim, const int *dims);
/*@}*/


/*
 * PnlMat
 */

#ifndef PNL_RANGE_CHECK_OFF
#define MGET(v,i,j) pnl_mat_get((v), (i), (j))
#define MLET(v,i,j) *(pnl_mat_lget((v), (i), (j)))
#else /* PNL_RANGE_CHECK_OFF */
#define MGET(v,i,j) (v)->array[(i)*(v)->n+(j)]
#define MLET(v,i,j) (v)->array[(i)*(v)->n+(j)]
#endif /* PNL_RANGE_CHECK_OFF */


/**
 * \ingroup PnlMat
 */
/*@{*/

extern void pnl_mat_init(PnlMat *); 
extern PnlMat* pnl_mat_new(); 
extern int pnl_mat_eq (const PnlMat *, const PnlMat *);
extern int pnl_mat_eq_all (const PnlMat *, double);
extern PnlMat* pnl_mat_create(int m, int n); 
extern PnlMat* pnl_mat_create_from_scalar(int m, int n, double x);
extern PnlMat* pnl_mat_create_from_zero(int m, int n);
extern PnlMat* pnl_mat_create_from_ptr(int m, int n, const double* x);
extern PnlMat* pnl_mat_create_from_list(int m, int n, ...); 
extern PnlMat pnl_mat_wrap_array(const double* x,int m, int n);
extern PnlMat* pnl_mat_create_from_file (const char * file);
extern void pnl_vect_extract_submat (PnlVect *M_sub, const PnlMat *M, const PnlVectInt *indi, const PnlVectInt *indj);
extern PnlVect* pnl_vect_create_submat (const PnlMat *M, const PnlVectInt *indi, const PnlVectInt *indj);
extern void pnl_mat_extract_subblock (PnlMat *M_sub, const PnlMat *M, int i, int len_i, int j, int len_j);
extern void pnl_mat_set_subblock (PnlMat *M, const PnlMat *block, int i, int j);
extern int pnl_mat_resize(PnlMat *v, int m, int n);
extern void pnl_mat_free(PnlMat **v);
extern PnlMat* pnl_mat_copy(const PnlMat *v);
extern void pnl_mat_clone(PnlMat *clone, const PnlMat *M);
extern void pnl_mat_map_inplace(PnlMat *lhs, double(*f)(double)); /*lhs=f(lhs)*/
extern void pnl_mat_map(PnlMat *lhs, const PnlMat *rhs, double(*f)(double));/* lhs(i)=f(rhs(i)) */
extern void pnl_mat_map_mat_inplace(PnlMat *lhs, const PnlMat *rhs, double(*f)(double,double));
extern void pnl_mat_map_mat(PnlMat *lhs, const PnlMat *rhs1, const PnlMat *rhs2, double(*f)(double,double));
extern int pnl_mat_find(PnlVectInt *indi, PnlVectInt *indj, char* type, int(*f)(double *), ...);
extern void pnl_mat_plus_scalar(PnlMat *lhs, double x); /*lhs+=x*/
extern void pnl_mat_minus_scalar(PnlMat *lhs, double x); /*lhs-=x*/
extern void pnl_mat_plus_mat(PnlMat *lhs, const PnlMat *rhs); /*lhs+=rhs*/
extern void pnl_mat_minus_mat(PnlMat *lhs, const PnlMat *rhs); /*lhs-=rhs*/
extern void pnl_mat_mult_scalar(PnlMat *lhs, double x); /*lhs*=x*/
extern void pnl_mat_div_scalar(PnlMat *lhs, double x); /*lhs/=x*/
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

extern void pnl_mat_set_all(PnlMat *lhs, double x);
extern void pnl_mat_set_zero(PnlMat *lhs);
extern void pnl_mat_set_id(PnlMat *lhs);
extern void pnl_mat_set_diag(PnlMat *lhs, double x, int d);
extern void pnl_mat_set_from_ptr(PnlMat *lhs, const double *x);
extern void pnl_mat_mult_mat_term(PnlMat *lhs, const PnlMat *rhs); 
extern void pnl_mat_div_mat_term(PnlMat *lhs, const PnlMat *rhs); 
extern void pnl_mat_sq_transpose(PnlMat *M);
extern PnlMat* pnl_mat_transpose(const PnlMat *M);
extern void pnl_mat_tr(PnlMat*, const PnlMat *M);
extern double pnl_mat_trace (const PnlMat *M);
extern void pnl_mat_print(const PnlMat *M);
extern void pnl_mat_fprint(FILE *fic, const PnlMat *M);
extern void pnl_mat_fprint_nsp (FILE *fic, const PnlMat *M);
extern void pnl_mat_print_nsp (const PnlMat *M);
extern PnlVect pnl_vect_wrap_mat_row(const PnlMat *M, int i);
extern PnlMat pnl_mat_wrap_vect(const PnlVect *v);
extern void pnl_mat_get_row(PnlVect *V, const PnlMat *M, int i);
extern void pnl_mat_get_col(PnlVect *V, const PnlMat *M, int j);
extern void pnl_mat_add_row(PnlMat *M, int i, const PnlVect *r);
extern void pnl_mat_del_row(PnlMat *M, int i);
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
#ifdef PNL_HAVE_INLINE 
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
/* @} */

/**
 * \ingroup PnlHmat
 */
/*@{*/


extern PnlHmat* pnl_hmat_new(); 
extern PnlHmat* pnl_hmat_create(int ndim, const int *dims); 
extern PnlHmat* pnl_hmat_create_from_scalar(int ndim, const int *dims, double x); 
extern PnlHmat* pnl_hmat_create_from_ptr(int ndim, const int *dims, const double *x);
extern int pnl_hmat_resize(PnlHmat *v, int ndim, const int *dims);
extern void pnl_hmat_free(PnlHmat **v);
extern PnlHmat* pnl_hmat_copy(const PnlHmat *H);
extern void pnl_hmat_clone(PnlHmat *clone, const PnlHmat *H);
extern void pnl_hmat_print(const PnlHmat *H);
extern void pnl_hmat_plus_hmat(PnlHmat *lhs, const PnlHmat *rhs);/*lhs+=rhs*/
extern void pnl_hmat_mult_scalar(PnlHmat *lhs, double x);/* lhs *=x;*/

extern void pnl_hmat_set(PnlHmat *H, int *tab, double x);
extern double pnl_hmat_get(const PnlHmat *H, int *tab);
extern double* pnl_hmat_lget(PnlHmat *H, int *tab);
extern PnlMat pnl_mat_wrap_hmat(PnlHmat *H, int *t);
extern PnlVect pnl_vect_wrap_hmat(PnlHmat *H, int *t);

/*@}*/


/*
 * PnlMatInt
 */

/**
 * \ingroup PnlMat
 */
/*@{*/

#ifndef PNL_RANGE_CHECK_OFF
#define MGET_INT(v,i,j) pnl_mat_int_get((v), (i), (j))
#define MLET_INT(v,i,j) *(pnl_mat_int_lget((v), (i), (j)))
#else /* PNL_RANGE_CHECK_OFF */
#define MGET_INT(v,i,j) (v)->array[(i)*(v)->n+(j)]
#define MLET_INT(v,i,j) (v)->array[(i)*(v)->n+(j)]
#endif /* PNL_RANGE_CHECK_OFF */


extern void pnl_mat_int_init(PnlMatInt *); 
extern PnlMatInt* pnl_mat_int_new(); 
extern int pnl_mat_int_eq (const PnlMatInt *, const PnlMatInt *);
extern int pnl_mat_int_eq_all (const PnlMatInt *, int);
extern PnlMatInt* pnl_mat_int_create(int m, int n); 
extern PnlMatInt* pnl_mat_int_create_from_scalar(int m, int n, int x);
extern PnlMatInt* pnl_mat_int_create_from_zero(int m, int n);
extern PnlMatInt* pnl_mat_int_create_from_ptr(int m, int n, const int* x);
extern PnlMatInt* pnl_mat_int_create_from_list(int m, int n, ...); 
extern PnlMatInt pnl_mat_int_wrap_array(const int* x,int m, int n);
extern PnlMatInt* pnl_mat_int_create_from_file (const char * file);
extern void pnl_vect_int_extract_submat (PnlVectInt *sub, const PnlMatInt *M, const PnlVectInt *indi, const PnlVectInt *indj);
extern PnlVectInt* pnl_vect_int_create_submat (const PnlMatInt *M, const PnlVectInt *indi, const PnlVectInt *indj);
extern void pnl_mat_int_extract_subblock (PnlMatInt *M_sub, const PnlMatInt *M, int i, int len_i, int j, int len_j);
extern void pnl_mat_int_set_subblock (PnlMatInt *M, const PnlMatInt *block, int i, int j);
extern int pnl_mat_int_resize(PnlMatInt *v, int m, int n);
extern void pnl_mat_int_free(PnlMatInt **v);
extern PnlMatInt* pnl_mat_int_copy(const PnlMatInt *v);
extern void pnl_mat_int_clone(PnlMatInt *clone, const PnlMatInt *M);
extern void pnl_mat_int_set_all(PnlMatInt *lhs, int x);/*lhs=x*/
extern void pnl_mat_int_set_zero(PnlMatInt *lhs);
extern void pnl_mat_int_set_id(PnlMatInt *lhs);
extern void pnl_mat_int_set_diag(PnlMatInt *lhs, int x, int d);
extern void pnl_mat_int_set_from_ptr(PnlMatInt *lhs, const int *x);
extern PnlMatInt* pnl_mat_int_transpose(const PnlMatInt *M);
extern void pnl_mat_int_tr(PnlMatInt *tM, const PnlMatInt *M);
extern void pnl_mat_int_sq_transpose(PnlMatInt *M);
extern int pnl_mat_int_trace (const PnlMatInt *M);
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
extern void pnl_mat_int_add_row(PnlMatInt *M, int i, const PnlVectInt *r);
extern void pnl_mat_int_del_row(PnlMatInt *M, int i);
extern PnlVectInt pnl_vect_int_wrap_mat_row(const PnlMatInt *M, int i);/* M(i,:)=V(:) */
extern PnlMatInt pnl_mat_int_wrap_vect(const PnlVectInt *V);
extern void pnl_mat_int_map_inplace(PnlMatInt *lhs, int(*f)(int)); /*lhs=f(lhs)*/
extern void pnl_mat_int_map(PnlMatInt *lhs, const PnlMatInt *rhs, int(*f)(int));/* lhs(i)=f(rhs(i)) */
extern void pnl_mat_int_map_mat_inplace(PnlMatInt *lhs, const PnlMatInt *rhs, int(*f)(int,int));
extern void pnl_mat_int_map_mat(PnlMatInt *lhs, const PnlMatInt *rhs1, const PnlMatInt *rhs2, int(*f)(int,int));
extern int pnl_mat_int_find(PnlVectInt *indi, PnlVectInt *indj, char* type, int(*f)(int *), ...);
extern void pnl_mat_int_plus_scalar(PnlMatInt *lhs, int x); /*lhs+=x*/
extern void pnl_mat_int_minus_scalar(PnlMatInt *lhs, int x); /*lhs-=x*/
extern void pnl_mat_int_plus_mat(PnlMatInt *lhs, const PnlMatInt *rhs); /*lhs+=rhs*/
extern void pnl_mat_int_minus_mat(PnlMatInt *lhs, const PnlMatInt *rhs); /*lhs-=rhs*/
extern void pnl_mat_int_mult_scalar(PnlMatInt *lhs, int x); /*lhs*=x*/
extern void pnl_mat_int_div_scalar(PnlMatInt *lhs, int x); /*lhs/=x*/
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
#ifdef PNL_HAVE_INLINE 
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


/**
 * \ingroup PnlHmat a Hyper Matrix
 */
/*@{*/



extern PnlHmatInt* pnl_hmat_int_new(); 
extern PnlHmatInt* pnl_hmat_int_create(int ndim, const int *dims); 
extern PnlHmatInt* pnl_hmat_int_create_from_scalar(int ndim, const int *dims, int x); 
extern PnlHmatInt* pnl_hmat_int_create_from_ptr(int ndim, const int *dims, const int *x);
extern int pnl_hmat_int_resize(PnlHmatInt *v, int ndim, const int *dims);
extern void pnl_hmat_int_free(PnlHmatInt **v);
extern PnlHmatInt* pnl_hmat_int_copy(const PnlHmatInt *H);
extern void pnl_hmat_int_print(const PnlHmatInt *H);
extern void pnl_hmat_int_clone(PnlHmatInt *clone, const PnlHmatInt *H);
extern void pnl_hmat_int_plus_hmat(PnlHmatInt *lhs, const PnlHmatInt *rhs);/*lhs+=rhs*/
extern void pnl_hmat_int_mult_scalar(PnlHmatInt *lhs, int x);/* lhs *=x;*/

extern void pnl_hmat_int_set(PnlHmatInt *H, int *tab, int x);
extern int pnl_hmat_int_get(const PnlHmatInt *H, int *tab);
extern int* pnl_hmat_int_lget(PnlHmatInt *H, int *tab);
extern PnlMatInt pnl_mat_int_wrap_hmat(PnlHmatInt *H, int *t);
extern PnlVectInt pnl_vect_int_wrap_hmat(PnlHmatInt *H, int *t);

/*@}*/


/*
 * PnlMatComplex
 */

#include "pnl/pnl_complex.h"

/**
 * \ingroup PnlMat
 */
/*@{*/

extern void pnl_mat_complex_init(PnlMatComplex *); 
extern PnlMatComplex* pnl_mat_complex_new(); 
extern int pnl_mat_complex_eq (const PnlMatComplex *, const PnlMatComplex *);
extern int pnl_mat_complex_eq_all (const PnlMatComplex *, dcomplex);
extern PnlMatComplex* pnl_mat_complex_create(int m, int n); 
extern PnlMatComplex* pnl_mat_complex_create_from_scalar(int m, int n, dcomplex x);
extern PnlMatComplex* pnl_mat_complex_create_from_zero(int m, int n);
extern PnlMatComplex* pnl_mat_complex_create_from_ptr(int m, int n, const dcomplex* x);
extern PnlMatComplex* pnl_mat_complex_create_from_mat (const PnlMat *R);
extern PnlMatComplex* pnl_mat_complex_create_from_list(int m, int n, ...); 
extern PnlMatComplex pnl_mat_complex_wrap_array(const dcomplex* x, int m, int n);
extern PnlMatComplex* pnl_mat_complex_create_from_file (const char * file);
extern void pnl_vect_complex_extract_submat (PnlVectComplex *M_sub, const PnlMatComplex *M, const PnlVectInt *indi, const PnlVectInt *indj);
extern PnlVectComplex* pnl_vect_complex_create_submat (const PnlMatComplex *M, const PnlVectInt *indi, const PnlVectInt *indj);
extern void pnl_mat_complex_extract_subblock (PnlMatComplex *M_sub, const PnlMatComplex *M, int i, int len_i, int j, int len_j);
extern void pnl_mat_complex_set_subblock (PnlMatComplex *M, const PnlMatComplex *block, int i, int j);
extern int pnl_mat_complex_resize(PnlMatComplex *v, int m, int n);
extern void pnl_mat_complex_free(PnlMatComplex **v);
extern PnlMatComplex* pnl_mat_complex_copy(const PnlMatComplex *v);
extern void pnl_mat_complex_clone(PnlMatComplex *clone, const PnlMatComplex *M);
extern void pnl_mat_complex_set_all(PnlMatComplex *lhs, dcomplex x);/*lhs=x*/
extern void pnl_mat_complex_set_zero(PnlMatComplex *lhs);
extern void pnl_mat_complex_set_id(PnlMatComplex *lhs);
extern void pnl_mat_complex_set_diag(PnlMatComplex *lhs, dcomplex x, int d);
extern void pnl_mat_complex_set_from_ptr(PnlMatComplex *lhs, const dcomplex *x);
extern PnlMatComplex* pnl_mat_complex_transpose(const PnlMatComplex *M);
extern void pnl_mat_complex_tr(PnlMatComplex*, const PnlMatComplex *M);
extern void pnl_mat_complex_sq_transpose(PnlMatComplex *M);
extern dcomplex pnl_mat_complex_trace (const PnlMatComplex *M);
extern void pnl_mat_complex_print(const PnlMatComplex *M);
extern void pnl_mat_complex_print_nsp(const PnlMatComplex *M);
extern void pnl_mat_complex_fprint(FILE *fic,const PnlMatComplex *M);
extern void pnl_mat_complex_fprint_nsp(FILE *fic,const PnlMatComplex *M);
extern PnlMatComplex* pnl_mat_complex_create_diag_from_ptr(const dcomplex x[], int d);
extern PnlMatComplex* pnl_mat_complex_create_diag(const PnlVectComplex *V);
extern void pnl_mat_complex_set_row(PnlMatComplex *M, const PnlVectComplex *V, int i);/* M(i,:)=V(:) */
extern void pnl_mat_complex_set_col(PnlMatComplex *M, const PnlVectComplex *V, int j);/* M(:,j)=V(:) */
extern void pnl_mat_complex_swap_rows (PnlMatComplex *M, int i, int j);
extern void pnl_mat_complex_get_row(PnlVectComplex *V, const PnlMatComplex *M, int i);/* V(:)=M(i,:) */
extern void pnl_mat_complex_get_col(PnlVectComplex *V, const PnlMatComplex *M, int j);
extern void pnl_mat_complex_add_row(PnlMatComplex *M, int i, const PnlVectComplex *r);
extern void pnl_mat_complex_del_row(PnlMatComplex *M, int i);
extern PnlVectComplex pnl_vect_complex_wrap_mat_row(const PnlMatComplex *M, int i);/* M(i,:)=V(:) */
extern PnlMatComplex pnl_mat_complex_wrap_vect(const PnlVectComplex *V);
extern void pnl_mat_complex_map_inplace(PnlMatComplex *lhs, dcomplex(*f)(dcomplex)); /*lhs=f(lhs)*/
extern void pnl_mat_complex_map(PnlMatComplex *lhs, const PnlMatComplex *rhs, dcomplex(*f)(dcomplex));/* lhs(i)=f(rhs(i)) */
extern void pnl_mat_complex_map_mat_inplace(PnlMatComplex *lhs, const PnlMatComplex *rhs, dcomplex(*f)(dcomplex,dcomplex));
extern void pnl_mat_complex_map_mat(PnlMatComplex *lhs, const PnlMatComplex *rhs1, const PnlMatComplex *rhs2, dcomplex(*f)(dcomplex,dcomplex));
extern int pnl_mat_complex_find(PnlVectInt *indi, PnlVectInt *indj, char* type, int(*f)(dcomplex *), ...);
extern void pnl_mat_complex_plus_scalar(PnlMatComplex *lhs, dcomplex x); /*lhs+=x*/
extern void pnl_mat_complex_minus_scalar(PnlMatComplex *lhs, dcomplex x); /*lhs-=x*/
extern void pnl_mat_complex_plus_mat(PnlMatComplex *lhs, const PnlMatComplex *rhs); /*lhs+=rhs*/
extern void pnl_mat_complex_minus_mat(PnlMatComplex *lhs, const PnlMatComplex *rhs); /*lhs-=rhs*/
extern void pnl_mat_complex_mult_scalar(PnlMatComplex *lhs, dcomplex x); /*lhs*=x*/
extern void pnl_mat_complex_div_scalar(PnlMatComplex *lhs, dcomplex x); /*lhs/=x*/
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
extern int pnl_mat_complex_eigen (PnlVectComplex *v, PnlMatComplex *P, const PnlMatComplex *A, int with_eigenvectors);
extern int pnl_mat_complex_log (PnlMatComplex *B, const PnlMatComplex *A);



/* inline functions if you are using GCC */
#ifdef PNL_HAVE_INLINE 
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
#define MGET_COMPLEX(v,i,j) pnl_mat_get_complex((v), (i), (j))
#define MLET_COMPLEX(v,i,j) *(pnl_mat_lget_complex((v), (i), (j)))
#else
#define MGET_REAL(v,i,j) ((v)->array[(i)*(v)->n+(j)]).r
#define MLET_REAL(v,i,j) ((v)->array[(i)*(v)->n+(j)]).r
#define MGET_IMAG(v,i,j) ((v)->array[(i)*(v)->n+(j)]).i
#define MLET_IMAG(v,i,j) ((v)->array[(i)*(v)->n+(j)]).i
#define MGET_COMPLEX(v,i,j) (v)->array[(i)*(v)->n+(j)]
#define MLET_COMPLEX(v,i,j) (v)->array[(i)*(v)->n+(j)]
#endif


/*@}*/


/**
 * \ingroup PnlHmat a Hyper Matrix
 */
/*@{*/

extern PnlHmatComplex* pnl_hmat_complex_new(); 
extern PnlHmatComplex* pnl_hmat_complex_create(int ndim, const int *dims); 
extern PnlHmatComplex* pnl_hmat_complex_create_from_scalar(int ndim, const int *dims, dcomplex x); 
extern PnlHmatComplex* pnl_hmat_complex_create_from_ptr(int ndim, const int *dims, const dcomplex *x);
extern int pnl_hmat_complex_resize(PnlHmatComplex *v, int ndim, const int *dims);
extern void pnl_hmat_complex_free(PnlHmatComplex **v);
extern void pnl_hmat_complex_print(const PnlHmatComplex *H);
extern PnlHmatComplex* pnl_hmat_complex_copy(const PnlHmatComplex *H);
extern void pnl_hmat_complex_clone(PnlHmatComplex *clone, const PnlHmatComplex *H);
extern void pnl_hmat_complex_plus_hmat(PnlHmatComplex *lhs, const PnlHmatComplex *rhs);/*lhs+=rhs*/
extern void pnl_hmat_complex_mult_scalar(PnlHmatComplex *lhs, dcomplex x);/* lhs *=x;*/

extern void pnl_hmat_complex_set(PnlHmatComplex *H, int *tab, dcomplex x);
extern dcomplex pnl_hmat_complex_get(const PnlHmatComplex *H, int *tab);
extern dcomplex* pnl_hmat_complex_lget(PnlHmatComplex *H, int *tab);
extern PnlMatComplex pnl_mat_complex_wrap_hmat(PnlHmatComplex *H, int *t);
extern PnlVectComplex pnl_vect_complex_wrap_hmat(PnlHmatComplex *H, int *t);

/*@}*/

/*
 * Some deprecated names
 */

#include "pnl/pnl_deprecated.h"

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_MATRIX_H */
