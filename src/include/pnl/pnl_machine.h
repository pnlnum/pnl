#ifndef _PNL_MACHINE_H 
#define _PNL_MACHINE_H 


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl/pnl_config.h"

#ifdef PNL_HAVE_FORTRAN_COMPILER
#  include "pnl/FC.h"
#  ifdef FC_GLOBAL
#    define C2F(name) (FC_GLOBAL(name,name))
#  else
#    define C2F(name) name##_
#  endif
#else
#  define C2F(name) name##_
#endif

#ifdef USE_INTERNAL_BLAS
#  define C2F(name) name##_
#endif


#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_MACHINE_H  */
