#ifndef _PNL_MATRIX_H
#define _PNL_MATRIX_H


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


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
/*@}*/

#include "pnl_object.h"
/**
 * \defgroup PnlHmatrices Hyper Matrix object
 */

typedef struct _PnlMatObject PnlMatObject;
typedef struct _PnlHmatObject PnlHmatObject;
typedef struct _PnlMat PnlMat;
typedef struct _PnlMatComplex PnlMatComplex;
typedef struct _PnlMatInt PnlMatInt;

/**
 * \private 
 * This structure is only used internally and should never be accessed directly.
 * It is only useful for handling the different types of vectors together
 */
struct _PnlMatObject
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows any PnlMatXXX pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int m; /*!< nb rows */ 
  int n; /*!< nb columns */ 
  int mn; /*!< product m*n */
  int mem_size; /*!< size of the memory block allocated for array */
  void *array; /*!< pointer to store the data row-wise */
  int owner; /*!< 1 if the object owns its array member, 0 otherwise */
}; 

struct _PnlHmatObject
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows any PnlHmatXXX pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int ndim; /*!< nb dimensions */ 
  int *dims; /*!< pointer to store the value of the ndim dimensions */ 
  int mn; /*!< product dim_1 *...*dim_ndim */
  void *array; /*!< pointer to store */
} ;

extern PnlMatObject* pnl_mat_object_new ();
extern int pnl_mat_object_resize(PnlMatObject *M, int m, int n);
extern PnlHmatObject* pnl_hmat_object_new ();
extern int pnl_hmat_object_resize(PnlHmatObject *H, int ndim, const int *dims);

#include "pnl_matrix_double.h"
#include "pnl_matrix_complex.h"
#include "pnl_matrix_int.h"

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_MATRIX_H */
