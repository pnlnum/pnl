#ifndef __VECTOR_H__
#define __VECTOR_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>

#ifndef PNL_RANGE_CHECK_OFF
#define CheckIndexVect(v,i) {                                           \
    if (i>=v->size || i<0) {perror("index out of range"); abort();}}
#define CheckVectMatch(lhs, rhs) { if ((lhs)->size != (rhs)->size)  \
      {perror("non compatible dimensions"); abort();}}
#else
#define CheckIndexVect(v,i) {}
#define CheckVectMatch(lhs, rhs){}
#endif /* PNL_RANGE_CHECK_OFF */

#define PNL_GET(v,i) (v)->array[i]
#define PNL_LET(v,i) (v)->array[i]
#define PNL_SET(v,i,x) (v)->array[i]=(x)


#include "pnl_vector_double.h"
#include "pnl_vector_complex.h"
#include "pnl_vector_int.h"
#include "pnl_vector_uint.h"

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* __VECTOR_H__ */
