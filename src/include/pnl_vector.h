#ifndef _PNL_VECTOR_H
#define _PNL_VECTOR_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>
#include <string.h>

#ifndef PNL_RANGE_CHECK_OFF
#define CheckIndexVect(v,i) {                                           \
    if (i>=v->size || i<0) {perror("index out of range"); abort();}}
#define CheckVectMatch(lhs, rhs) { if ((lhs)->size != (rhs)->size)  \
      {perror("non compatible dimensions"); abort();}}
#else
#define CheckIndexVect(v,i) {}
#define CheckVectMatch(lhs, rhs){}
#endif /* PNL_RANGE_CHECK_OFF */

#include "pnl_object.h"
/**
 * \defgroup PnlVectors  a Vector object
 */

#define PNL_GET(v,i) (v)->array[i]
#define PNL_LET(v,i) (v)->array[i]
#define PNL_SET(v,i,x) (v)->array[i]=(x)

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
  
extern PnlVectObject* pnl_vect_object_new ();

#include "pnl_vector_double.h"
#include "pnl_vector_complex.h"
#include "pnl_vector_int.h"

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_VECTOR_H */
