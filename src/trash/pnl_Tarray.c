#include <stdlib.h>

#include "pnl/pnl_config.h"
#include "pnl/pnl_mathtools.h"


#define BASE_DOUBLE
#include "pnl/pnl_templates_on.h"
#include "pnl/pnl_Tarray_source.c"
#include "pnl/pnl_templates_off.h"
#undef BASE_DOUBLE


#define BASE_PNL_COMPLEX
#include "pnl/pnl_templates_on.h"
#include "pnl/pnl_Tarray_source.c"
#include "pnl/pnl_templates_off.h"
#undef BASE_PNL_COMPLEX


#define BASE_INT
#include "pnl/pnl_templates_on.h"
#include "pnl/pnl_Tarray_source.c"
#include "pnl/pnl_templates_off.h"
#undef  BASE_INT

#define BASE_UINT
#include "pnl/pnl_templates_on.h"
#include "pnl/pnl_Tarray_source.c"
#include "pnl/pnl_templates_off.h"
#undef  BASE_UINT


