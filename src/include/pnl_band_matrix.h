#ifndef BAND_MATRIX_H
#define BAND_MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl_matrix.h"
#include "pnl_vector.h"

#ifndef PNL_RANGE_CHECK_OFF
#define CheckIndexBandMat(v,i,j){             \
    if (i>=v->size || abs(j)>1  || i<0 )       \
      {perror("index out of range"); abort();\
      }}
#define CheckBandMatMatch(lhs, rhs){          \
    if ((lhs)->size != (rhs)->size )             \
      {perror("non compatible dimensions"); abort();\
}}
#define CheckBandMatVectIsCompatible(mat, vect){            \
    if((mat)->size != (vect)->size)                        \
      {perror("non compatible dimensions"); abort();}}
#else
#define CheckIndexBandMat(v,i,j) {}                          
#define CheckBandMatMatch(lhs, rhs) {}
#define CheckBandMatVectIsCompatible(mat, vect){}
#endif /* PNL_RANGE_CHECK_OFF */

/** \defgroup PnlBandMat Band Matrix
 *
 * A standard Lapack band storage is used.
 *
 * An m x n band matrix with nl lowerdiagonals and nu upperdiagonals is stored
 * in a two-dimensional array with kl+ku+1  rows and n columns. Columns of the
 * matrix are stored in corresponding columns of the array, and diagonals of the
 * matrix are stored in rows of the array. 
 *
 * To be precise, M(i,j) is stored in B(nu+i-j,j) for max(0,j-nu) <= i <=
 * min(m,j+nl)
 */
/*@{*/
typedef struct
{
  int m; /*!< nb rows */ 
  int n; /*!< nb columns */ 
  int nu; /*!< nb of upperdiagonals */
  int nl; /*!< nb of lowerdiagonals */
  int m_band; /*!< nb rows of the band storage */
  int n_band; /*!< nb columns of the band storage */
  double *array;  /*!< a block to store the bands */  
} PnlBandMat;

extern PnlBandMat* pnl_bandmat_create (int m, int n, int nl, int nu);
extern PnlBandMat* pnl_bandmat_create_from_mat (const PnlMat *M, int nl, int nu);
extern void pnl_bandmat_free (PnlBandMat **BM);
extern int pnl_bandmat_resize (PnlBandMat *BM, int m, int n, int nl, int nu);
extern void pnl_bandmat_clone (PnlBandMat * Bclone, const PnlBandMat * BM);
extern PnlBandMat* pnl_bandmat_copy (const PnlBandMat * BM);
extern PnlMat* pnl_bandmat_to_mat (const PnlBandMat *BM);
extern void pnl_bandmat_print_as_full (const PnlBandMat *BM);
extern void pnl_bandmat_map (PnlBandMat *lhs, const PnlBandMat *rhs, double(*f)(double));
extern void pnl_bandmat_map_inplace (PnlBandMat *BM, double(*f)(double));
extern void pnl_bandmat_map_bandmat (PnlBandMat *BA, const PnlBandMat *BB, double(*f)(double, double));

extern void pnl_bandmat_plus_double (PnlBandMat *BM, double x);
extern void pnl_bandmat_minus_double (PnlBandMat *BM, double x);
extern void pnl_bandmat_mult_double (PnlBandMat *BM, double x);
extern void pnl_bandmat_div_double (PnlBandMat *BM, double x);

extern void pnl_bandmat_plus_bandmat (PnlBandMat *lhs, const PnlBandMat *rhs);
extern void pnl_bandmat_minus_bandmat (PnlBandMat *lhs, const PnlBandMat *rhs);
extern void pnl_bandmat_div_bandmat_term (PnlBandMat *lhs, const PnlBandMat *rhs);
extern void pnl_bandmat_mult_bandmat_term (PnlBandMat *lhs, const PnlBandMat *rhs);

extern double pnl_bandmat_get (PnlBandMat * M, int i, int j);
extern double* pnl_bandmat_lget (PnlBandMat * M, int i, int j);
extern void pnl_bandmat_set (PnlBandMat * M, int i, int j, double x);
extern void pnl_bandmat_set_double (PnlBandMat* BM, double x);

extern void pnl_bandmat_mult_vect_inplace (PnlVect *lhs, const PnlBandMat *mat, const 
PnlVect *rhs);
extern void pnl_bandmat_lAxpby (double l, const PnlBandMat *A, const PnlVect *x, double b, PnlVect * y);

extern void pnl_bandmat_lu (PnlBandMat *BM, PnlVectInt *p);
extern void pnl_bandmat_syslin_inplace (PnlBandMat *BM, PnlVect *b);
extern void pnl_bandmat_syslin (PnlVect *x, PnlBandMat *BM, const PnlVect *b);
extern void pnl_bandmat_lu_syslin_inplace (const PnlBandMat *LU, const PnlVectInt *p, PnlVect *b);
extern void pnl_bandmat_lu_syslin (PnlVect *x, const PnlBandMat *LU, const PnlVectInt *p, const PnlVect *b);

/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* BAND_MATRIX_H */
