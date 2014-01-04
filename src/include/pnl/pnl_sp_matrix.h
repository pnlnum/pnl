#ifndef _PNL_SP_MATRIX_H
#define _PNL_SP_MATRIX_H 

#include "pnl/pnl_object.h"
#include "pnl/pnl_complex.h"
#include "pnl/pnl_matrix.h"

#ifndef PNL_RANGE_CHECK_OFF
#define CheckIndexSpMat(M,i,j) \
  if ( (i)>((M)->m-1) || (j)>((M)->n-1) || (i)<0 || (j)<0 ) { perror("index of range"); abort(); }

#endif /* PNL_RANGE_CHECK_OFF */

#include "pnl/pnl_config.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/**
 * \defgroup PnlSpMatrices a sparse matrix object
 *
 * Sparse matrices are stored using compressed storage row format.
 * Full[i,j] = Sparse[I[i] + j] 
 */
/*@{*/


typedef struct _PnlSpMatObject PnlSpMatObject;
typedef struct _PnlSpMat PnlSpMat;
typedef struct _PnlSpMatInt PnlSpMatInt;
typedef struct _PnlSpMatComplex PnlSpMatComplex;

/**
 * \private 
 * This structure is only used internally and should never be accessed directly.
 * It is only useful for handling the different types of sparse matrices together.
 */
struct _PnlSpMatObject
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows any PnlSpMatXXX pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int m; /*!< number of rows */
  int n; /*!< number of columns */
  int nz; /*!< number of non-zero elements */
  int *J; /*!< column indices, vector of size nzmax */
  int *I; /*!< row offset integer vector, 
            array[I[i]] is the first element of row i.
            Vector of size m */ 
  void *array; /*!< pointer to store the data of size nzmax */
  int nzmax; /*!< size of the memory block allocated for array */
};

extern PnlSpMatObject* pnl_sp_mat_object_new ();
extern void pnl_sp_mat_object_free (PnlSpMatObject **);
extern int pnl_sp_mat_object_resize(PnlSpMatObject *M, int m, int n, int nzmax);


/**
 * \ingroup PnlSpMat
 */
/*@{*/
struct _PnlSpMat
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows a PnlSpMat pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int m; /*!< number of rows */
  int n; /*!< number of columns */
  int nz; /*!< number of non-zero elements */
  int *J; /*!< column indices, vector of size nzmax */
  int *I; /*!< row offset integer vector, 
            array[I[i]] is the first element of row i.
            Vector of size (m+1) */ 
  double *array; /*!< pointer to store the data of size nzmax */
  int nzmax; /*!< size of the memory block allocated for array */
};

extern PnlSpMat* pnl_sp_mat_new ();
extern void pnl_sp_mat_free (PnlSpMat **);
extern PnlSpMat* pnl_sp_mat_copy (const PnlSpMat *M);
extern void pnl_sp_mat_clone (PnlSpMat *dest, const PnlSpMat *src);
extern PnlSpMat* pnl_sp_mat_create(int m, int n, int nzmax);
extern int pnl_sp_mat_resize(PnlSpMat *M, int m, int n, int nzmax);
extern void pnl_sp_mat_set(PnlSpMat *M, int i, int j, double x);
extern double pnl_sp_mat_get(const PnlSpMat *M, int i, int j);
extern PnlMat* pnl_mat_create_from_sp_mat(const PnlSpMat *M);
/*@}*/


/**
 * \ingroup PnlSpMatInt
 */
/*@{*/
struct _PnlSpMatInt
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows a PnlSpMatInt pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int m; /*!< number of rows */
  int n; /*!< number of columns */
  int nz; /*!< number of non-zero elements */
  int *J; /*!< column indices, vector of size nzmax */
  int *I; /*!< row offset integer vector, 
            array[I[i]] is the first element of row i.
            Vector of size (m+1) */ 
  int *array; /*!< pointer to store the data of size nzmax */
  int nzmax; /*!< size of the memory block allocated for array */
}; 

extern PnlSpMatInt* pnl_sp_mat_int_new ();
extern void pnl_sp_mat_int_free (PnlSpMatInt **);
extern PnlSpMatInt* pnl_sp_mat_int_copy (const PnlSpMatInt *M);
extern void pnl_sp_mat_int_clone (PnlSpMatInt *dest, const PnlSpMatInt *src);
extern PnlSpMatInt* pnl_sp_mat_int_create(int m, int n, int nzmax);
extern int pnl_sp_mat_int_resize(PnlSpMatInt *M, int m, int n, int nzmax);
extern void pnl_sp_mat_int_set(PnlSpMatInt *M, int i, int j, int x);
extern int pnl_sp_mat_int_get(const PnlSpMatInt *M, int i, int j);
extern PnlMatInt* pnl_mat_int_create_from_sp_mat (const PnlSpMatInt *M);
/*@}*/


/**
 * \ingroup PnlSpMatComplex
 */
/*@{*/
struct _PnlSpMatComplex
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows a PnlSpMatComplex pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int m; /*!< number of rows */
  int n; /*!< number of columns */
  int nz; /*!< number of non-zero elements */
  int *J; /*!< column indices, vector of size nzmax */
  int *I; /*!< row offset integer vector, 
            array[I[i]] is the first element of row i.
            Vector of size (m+1) */ 
  dcomplex *array; /*!< pointer to store the data of size nzmax */
  int nzmax; /*!< size of the memory block allocated for array */
};

extern PnlSpMatComplex* pnl_sp_mat_complex_new ();
extern void pnl_sp_mat_complex_free (PnlSpMatComplex **);
extern PnlSpMatComplex* pnl_sp_mat_complex_copy (const PnlSpMatComplex *M);
extern void pnl_sp_mat_complex_clone (PnlSpMatComplex *dest, const PnlSpMatComplex *src);
extern PnlSpMatComplex* pnl_sp_mat_complex_create(int m, int n, int nzmax);
extern int pnl_sp_mat_complex_resize(PnlSpMatComplex *M, int m, int n, int nzmax);
extern void pnl_sp_mat_complex_set(PnlSpMatComplex *M, int i, int j, dcomplex x);
extern dcomplex pnl_sp_mat_complex_get(const PnlSpMatComplex *M, int i, int j);
extern PnlMatComplex* pnl_mat_complex_create_from_sp_mat (const PnlSpMatComplex *M);
/*@}*/
/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _PNL_SP_MATRIX_H */