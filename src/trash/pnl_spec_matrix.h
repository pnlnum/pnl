#ifndef SPECIAL_MATRIX_H
#define SPECIAL_MATRIX_H


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl/pnl_perm.h"
#include "pnl/pnl_matrix.h"
#include "cs.h"

#ifndef PNL_RANGE_CHECK_OFF
#define CheckIndexSparseMat(v,i,j){             \
    if (i>=v->m || j>=v->n || j<0  || i<0 )       \
      {perror("index out of range"); abort();\
      }}
#define CheckSparseMatVectIsCompatible(mat, vect){            \
    if((mat)->m != (vect)->size)                        \
      {perror("non compatible dimensions"); abort();}}
#define CheckSparseMatMatch(lhs,rhs){             \
  if(lhs->nzmax!=rhs->nzmax && lhs->n!=rhs->n &&  \
     lhs->m!=rhs->m && lhs->nz!=rhs->nz)           \
      {perror("non compatible dimensions"); abort();}}
#else
#define CheckIndexSparseMat(v,i,j) {}                          
#define CheckSparseMatMatch(lhs, rhs) {}
#define CheckSparseMatVectIsCompatible(mat, vect){}
#endif /* PNL_RANGE_CHECK_OFF */

/**
 * \defgroup PnlMorseMat Morse Matrix 
 *
 * Morse matrices are only used for for construction 
 * 
 * m,n :  size of matrix
 *
 * D[i]-> points to the ith column
 *
 * D[i]-> size : number of non zero elements in column i
 *
 * D[i]->Index : array of the index  of the elements of line i
 *
 * D[i]->Value : array of the real part of the elements of line i
 *
 * D[i]->C : array of the imaginary part of the elements of line i
 *
 */
/*@{*/


/* used to store a matlab compatible representation */
typedef struct _sprow  SpRow ;

struct _sprow {
  int size; /*!< size of a row */
  int Max_size; /*!< max size allocation of a row */
  int    *Index; /*!< pointer to an int array giving the columns or row i */
  double *Value; /*!< Pointer on values */
};



typedef struct PnlMorseMat{
  int m; /*!< nb rows */ 
  int n; /*!< nb columns */ 
  SpRow * array; 
  /*!< pointer in each row or col to store no nul coefficients */
  int RC; /*!> 0 if store row, 1 if store col*/ 
} PnlMorseMat;


extern PnlMorseMat * pnl_morse_mat_create(int m, int n,int Max_size_row,int RC);
extern PnlMorseMat * pnl_morse_mat_create_fromfull(PnlMat * FM,int RC);
extern void pnl_morse_mat_free(PnlMorseMat ** M);

extern double pnl_morse_mat_get(PnlMorseMat* M, int i, int j);
extern int pnl_morse_mat_set(PnlMorseMat* M, int i, int j,double Val);
extern double* pnl_morse_mat_lget(PnlMorseMat* M, int i, int j);
extern int pnl_morse_mat_freeze(PnlMorseMat* M);
extern void pnl_morse_mat_mult_vect_inplace(PnlVect *lhs, const PnlMorseMat *M, const PnlVect *rhs);
extern PnlVect* pnl_morse_mat_mult_vect(const PnlMorseMat *M, const PnlVect *vec);
extern void pnl_morse_mat_print (const PnlMorseMat *M);
extern PnlMat * pnl_morse_mat_full(PnlMorseMat * M);

typedef struct _spwaverow  SpWaveRow ;

struct _spwaverow {
  int level_size;
  SpRow * L_Value;
};


typedef struct PnlMorseWaveMat{
  int m; /*!< nb rows */ 
  int n; /*!< nb columns */ 
  SpWaveRow * array; 
  /*!< pointer in each row and each column level to store no nul coefficients */
} PnlMorseWaveMat;



typedef struct PnlMorseLUFact
/* LU, factorization */
{
    PnlMorseMat *L ; /* L for LU and Cholesky, V for QR */
    PnlMorseMat *U ; /* U for LU, R for QR, not used for Cholesky */
    PnlPermutation *pinv; /* partial pivoting for LU */
}PnlMorseLUFact;

extern PnlMorseLUFact *pnl_morsefact_empty_create(int m,int n,int max_size);
extern PnlMorseLUFact *pnl_morse_lu_fact_create(const PnlMorseMat *A,double tol);
extern void pnl_morse_lu_fact_solve_inplace(const PnlMorseLUFact *F, PnlVect *x);
extern void pnl_morsefact_free(PnlMorseLUFact ** Fatc);

/*@}*/

/**
 * \defgroup PnlSparseMat Sparse Matrix 
 *
 * This is the cs struct of Csparse library written by Timothy A.Davis
 *
 * For convenience, we have renamed somes functions and structures.
 * We have also reduced then numbers of function parameters for non
 * expert in sparse matrix.
 *
 * In the following, we only use the LU factorization for sparse system.  The
 * factorization can be stored for PDE-parabolic problem with the same
 * operator at each time discretisation step.
 *
 */
/*@{*/

typedef cs PnlSparseMat;

extern PnlSparseMat *pnl_sparse_mat_create_fromfull(PnlMat * M);
extern PnlSparseMat *pnl_sparse_mat_create_frommorse(PnlMorseMat * M);
extern void pnl_sparse_mat_free(PnlSparseMat **M);
extern int pnl_sparse_mat_gaxpby(PnlVect *lhs, const PnlSparseMat *M, const PnlVect *rhs);
extern int pnl_sparse_mat_mult_vect_inplace(PnlVect *lhs, const PnlSparseMat *M, const PnlVect *rhs);
extern void pnl_sparse_mat_print(PnlSparseMat *A);
extern void pnl_sparse_mat_map_inplace(PnlSparseMat *lhs,double(*f)(double ));
extern void pnl_sparse_mat_plus_double(PnlSparseMat *lhs , double x);
extern void pnl_sparse_mat_minus_double(PnlSparseMat *lhs , double x);
extern void pnl_sparse_mat_mult_double(PnlSparseMat *lhs , double x);
extern void pnl_sparse_mat_div_double(PnlSparseMat *lhs , double x);
extern void pnl_sparse_mat_plus_mat(PnlSparseMat *lhs, const PnlSparseMat *rhs);
extern void pnl_sparse_mat_minus_mat(PnlSparseMat *lhs, const PnlSparseMat *rhs);
extern void pnl_sparse_mat_inv_term(PnlSparseMat *lhs);
extern void pnl_sparse_mat_div_mat_term(PnlSparseMat *lhs, const PnlSparseMat *rhs);
extern void pnl_sparse_mat_mult_mat_term(PnlSparseMat *lhs, const PnlSparseMat *rhs);

typedef struct PnlSparseFactorization
{
  css *S ;
  csn *N ;
}PnlSparseFactorization;

extern PnlSparseFactorization * pnl_sparse_factorization_lu_create (const PnlSparseMat*A, double tol);
extern void pnl_sparse_factorization_free(PnlSparseFactorization ** F);
extern void pnl_sparse_factorization_lu_syslin(const PnlSparseFactorization * N, PnlVect *x);

/*@}*/


#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* SPECIAL_MATRIX_H */
