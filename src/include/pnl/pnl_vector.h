#ifndef _PNL_VECTOR_H
#define _PNL_VECTOR_H

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

#include "pnl/pnl_config.h"

/**
 * \defgroup PnlVectors  a Vector object
 */

#define PNL_GET(v,i) (v)->array[i]
#define PNL_LET(v,i) (v)->array[i]
#define PNL_SET(v,i,x) (v)->array[i]=(x)

#include "pnl/pnl_matvect.h"
#include "pnl/pnl_vector_double.h"
#include "pnl/pnl_vector_complex.h"
#include "pnl/pnl_vector_int.h"
  
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

extern PnlVectObject* pnl_vect_object_new ();
extern void pnl_vect_object_free (PnlVectObject **);
extern int pnl_vect_object_resize(PnlVectObject * v, int size);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_VECTOR_H */
