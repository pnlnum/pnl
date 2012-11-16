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
 * \defgroup PnlMatrices a Matrix object
 *
 * Matrix are stored row-wise in a one dimensional array.
 * The element (i,j) of a matrix is stored in array[i*n+j]
 */
/*@{*/
#define PNL_MGET(v,i,j) (v)->array[(i)*(v)->n+(j)]
#define PNL_MSET(v,i,j, x) (v)->array[(i)*(v)->n+(j)] = (x)
#define PNL_MLET(v,i,j) (v)->array[(i)*(v)->n+(j)]

#include "pnl/pnl_matvect.h"
#include "pnl/pnl_matrix_double.h"
#include "pnl/pnl_matrix_complex.h"
#include "pnl/pnl_matrix_int.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

extern PnlMatObject* pnl_mat_object_new ();
extern void pnl_mat_object_free (PnlMatObject **);
extern int pnl_mat_object_resize(PnlMatObject *M, int m, int n);
/*@}*/

/**
 * \defgroup PnlHmatrices Hyper Matrix object
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


#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_MATRIX_H */
