#include <stdlib.h>

#include "config.h"
#include "pnl_vector.h"
#include "pnl_matrix.h"

typedef int(*cmp_func)(const void *, const void *);

#define BASE_UINT
#include "pnl_templates_on.h"
#include "pnl_TMatrix_source.c"
#include "dgemm_source.c"
#include "pnl_templates_off.h"
#undef  BASE_UINT


#define BASE_INT
#include "pnl_templates_on.h"
#include "pnl_TMatrix_source.c"
#include "dgemm_source.c"
#include "pnl_templates_off.h"
#undef  BASE_INT

#define BASE_DOUBLE
#include "pnl_templates_on.h"
#include "pnl_TMatrix_source.c"
#include "dgemm_source.c"
#include "pnl_templates_off.h"
#undef  BASE_DOUBLE


#define BASE_PNL_COMPLEX
#include "pnl_templates_on.h"
#include "pnl_TMatrix_source.c"
#include "dgemm_source.c"
#include "pnl_templates_off.h"
#undef  BASE_PNL_COMPLEX

