#ifndef _PNL_MATVECT_H
#define _PNL_MATVECT_H

#include "pnl/pnl_object.h"
#include "pnl/pnl_complex.h"

/**
 * \ingroup PnlVectors
 */
/*@{*/


typedef struct _PnlVectObject PnlVectObject;
typedef struct _PnlVect PnlVect;
typedef struct _PnlVectInt PnlVectInt;
typedef struct _PnlVectComplex PnlVectComplex;

/**
 * \private 
 * This structure is only used internally and should never be accessed directly.
 * It is only useful for handling the different types of vectors together
 */
struct _PnlVectObject
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows any PnlVectXXX pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int size;/*!< size of the vector */ 
  void *array;/*!< pointer to store the data */
  int mem_size; /*!< size of the memory block allocated for array */
  int owner; /*!< 1 if the object owns its array member, 0 otherwise */
};

/**
 * \ingroup PnlVect
 */
struct _PnlVect
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows a PnlVect pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int size;/*!< size of the vector */ 
  double *array;/*!< pointer to store the data */
  int mem_size; /*!< size of the memory block allocated for array */
  int owner; /*!< 1 if the object owns its array member, 0 otherwise */
};
/*@}*/


/**
 * \ingroup PnlVectInt
 */
/*@{*/
struct _PnlVectInt
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows a PnlVectInt pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int size;/*!< size of the vector */ 
  int *array;/*!< pointer to store the data */
  int mem_size; /*!< size of the memory block allocated for array */
  int owner; /*!< 1 if the object owns its array member, 0 otherwise */
}; 
/*@}*/


/**
 * \ingroup PnlVectComplex
 */
/*@{*/
struct _PnlVectComplex
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows a PnlVectComplex pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int size;/*!< size of the vector */ 
  dcomplex *array;/*!< pointer to store the data */
  int mem_size; /*!< size of the memory block allocated for array */
  int owner; /*!< 1 if the object owns its array member, 0 otherwise */
};
/*@}*/
/*@}*/

/**
 * \ingroup PnlMatrices
 */
/*@{*/

typedef struct _PnlMatObject PnlMatObject;
typedef struct _PnlMat PnlMat;
typedef struct _PnlMatInt PnlMatInt;
typedef struct _PnlMatComplex PnlMatComplex;

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


/**
 * \ingroup PnlMat
 */
/*@{*/

struct _PnlMat
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
  double *array; /*!< pointer to store the data row-wise */
  int owner; /*!< 1 if the object owns its array member, 0 otherwise */
};
/*@}*/

/**
 * \ingroup PnlMatInt
 */
/*@{*/
struct _PnlMatInt
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
  int *array; /*!< pointer to store the data row-wise */
  int owner; /*!< 1 if the object owns its array member, 0 otherwise */
};
/*@}*/

/**
 * \ingroup PnlMatComplex
 */
/*@{*/
struct _PnlMatComplex
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
  dcomplex *array; /*!< pointer to store the data row-wise */
  int owner; /*!< 1 if the object owns its array member, 0 otherwise */
};

/*@}*/
/*@}*/

/**
 * \ingroup PnlHmatrices
 */
/*@{*/

typedef struct _PnlHmatObject PnlHmatObject;
typedef struct _PnlHmat PnlHmat;
typedef struct _PnlHmatInt PnlHmatInt;
typedef struct _PnlHmatComplex PnlHmatComplex;


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


/**
 * \ingroup PnlHMat
 */
/*@{*/

struct _PnlHmat
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows any PnlMatXXX pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int ndim; /*!< nb dimensions */ 
  int *dims; /*!< pointer to store the value of the ndim dimensions */ 
  int mn; /*!< product dim_1 *...*dim_ndim */
  double *array; /*!< pointer to store */
};

/**
 * \ingroup PnlHMatInt
 */
/*@{*/
struct _PnlHmatInt
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows any PnlMatXXX pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int ndim; /*!< nb dimensions */ 
  int *dims; /*!< pointer to store the value of the ndim dimensions */ 
  int mn; /*!< product dim_1 *...*dim_ndim */
  int *array; /*!< pointer to store */
};


/**
 * \ingroup PnlHMatComplex
 */
/*@{*/
struct _PnlHmatComplex
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows any PnlMatXXX pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int ndim; /*!< nb dimensions */ 
  int *dims; /*!< pointer to store the value of the ndim dimensions */ 
  int mn; /*!< product dim_1 *...*dim_ndim */
  dcomplex *array; /*!< pointer to store */
};
/*@}*/
/*@}*/



#endif /* _PNL_MATVECT_H */