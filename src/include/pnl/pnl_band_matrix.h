#ifndef _PNL_BAND_MATRIX_H
#define _PNL_BAND_MATRIX_H

#include "pnl/pnl_matrix.h"
#include "pnl/pnl_vector.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


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

typedef struct _PnlBandMatObject PnlBandMatObject;
typedef struct _PnlBandMat PnlBandMat;

struct _PnlBandMat
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows any PnlVectXXX pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int m; /*!< nb rows */ 
  int n; /*!< nb columns */ 
  int nu; /*!< nb of upperdiagonals */
  int nl; /*!< nb of lowerdiagonals */
  int m_band; /*!< nb rows of the band storage */
  int n_band; /*!< nb columns of the band storage */
  double *array;  /*!< a block to store the bands */  
} ;

struct _PnlBandMatObject
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows any PnlVectXXX pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int m; /*!< nb rows */ 
  int n; /*!< nb columns */ 
  int nu; /*!< nb of upperdiagonals */
  int nl; /*!< nb of lowerdiagonals */
  int m_band; /*!< nb rows of the band storage */
  int n_band; /*!< nb columns of the band storage */
  void *array;  /*!< a block to store the bands */  
} ;

extern PnlBandMatObject* pnl_band_mat_object_new ();
extern void pnl_band_mat_object_free (PnlBandMatObject **BM);
extern PnlBandMat* pnl_band_mat_new ();
extern PnlBandMat* pnl_band_mat_create (int m, int n, int nl, int nu);
extern PnlBandMat* pnl_band_mat_create_from_mat (const PnlMat *M, int nl, int nu);
extern void pnl_band_mat_free (PnlBandMat **BM);
extern int pnl_band_mat_resize (PnlBandMat *BM, int m, int n, int nl, int nu);
extern void pnl_band_mat_clone (PnlBandMat * Bclone, const PnlBandMat * BM);
extern PnlBandMat* pnl_band_mat_copy (const PnlBandMat * BM);
extern PnlMat* pnl_band_mat_to_mat (const PnlBandMat *BM);
extern void pnl_band_mat_print_as_full (const PnlBandMat *BM);
extern void pnl_band_mat_map (PnlBandMat *lhs, const PnlBandMat *rhs, double(*f)(double));
extern void pnl_band_mat_map_inplace (PnlBandMat *BM, double(*f)(double));
extern void pnl_band_mat_map_band_mat_inplace (PnlBandMat *BA, const PnlBandMat *BB, double(*f)(double, double));

extern void pnl_band_mat_plus_scalar (PnlBandMat *BM, double x);
extern void pnl_band_mat_minus_scalar (PnlBandMat *BM, double x);
extern void pnl_band_mat_mult_scalar (PnlBandMat *BM, double x);
extern void pnl_band_mat_div_scalar (PnlBandMat *BM, double x);

extern void pnl_band_mat_plus_band_mat (PnlBandMat *lhs, const PnlBandMat *rhs);
extern void pnl_band_mat_minus_band_mat (PnlBandMat *lhs, const PnlBandMat *rhs);
extern void pnl_band_mat_div_band_mat_term (PnlBandMat *lhs, const PnlBandMat *rhs);
extern void pnl_band_mat_mult_band_mat_term (PnlBandMat *lhs, const PnlBandMat *rhs);

extern double pnl_band_mat_get (PnlBandMat * M, int i, int j);
extern double* pnl_band_mat_lget (PnlBandMat * M, int i, int j);
extern void pnl_band_mat_set (PnlBandMat * M, int i, int j, double x);
extern void pnl_band_mat_set_all (PnlBandMat* BM, double x);

extern void pnl_band_mat_mult_vect_inplace (PnlVect *lhs, const PnlBandMat *mat, const 
PnlVect *rhs);
extern void pnl_band_mat_lAxpby (double l, const PnlBandMat *A, const PnlVect *x, double b, PnlVect * y);

extern void pnl_band_mat_lu (PnlBandMat *BM, PnlVectInt *p);
extern void pnl_band_mat_syslin_inplace (PnlBandMat *BM, PnlVect *b);
extern void pnl_band_mat_syslin (PnlVect *x, PnlBandMat *BM, const PnlVect *b);
extern void pnl_band_mat_lu_syslin_inplace (const PnlBandMat *LU, const PnlVectInt *p, PnlVect *b);
extern void pnl_band_mat_lu_syslin (PnlVect *x, const PnlBandMat *LU, const PnlVectInt *p, const PnlVect *b);

/*@}*/

/*
 * Somes deprecated names
 */

#include "pnl/pnl_deprecated.h"


#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_BAND_MATRIX_H */
