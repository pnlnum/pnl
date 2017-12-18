#ifndef _PNL_EXTERN_H
#define _PNL_EXTERN_H

#include "pnl/pnl_config.h"

#ifdef _WIN32
  #if (!defined(_COMPILING_PNL)) && defined(PNL_DLL)
    #define DECLEXTERN extern __declspec(dllimport)
  #else
    #define DECLEXTERN extern
  #endif
#else
  #define DECLEXTERN extern
#endif

#endif /* _PNL_EXTERN_H */
