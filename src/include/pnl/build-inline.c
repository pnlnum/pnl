
/* To prevent from loading pnl/pnl_config.h */
#define _PNL_CONFIG_H 

#ifdef PNL_HAVE_INLINE
#undef PNL_HAVE_INLINE
#endif
#define PNL_HAVE_INLINE

#ifdef PNL_INLINE_DECL
#undef PNL_INLINE_DECL
#endif
#define PNL_INLINE_DECL

#ifdef PNL_INLINE_FUNC
#undef PNL_INLINE_FUNC
#endif
#define PNL_INLINE_FUNC

#include "pnl_matrix.h"
#include "pnl_vector.h"
#include "pnl_random.h"
