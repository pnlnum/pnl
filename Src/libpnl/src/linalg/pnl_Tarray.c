#include <stdlib.h>

#include "config.h"
#include "pnl_mathtools.h"


#define BASE_DOUBLE
#include "pnl_templates_on.h"
#include "pnl_Tarray_source.c"
#include "pnl_templates_off.h"
#undef BASE_DOUBLE


#define BASE_PNL_COMPLEX
#include "pnl_templates_on.h"
#include "pnl_Tarray_source.c"
#include "pnl_templates_off.h"
#undef BASE_PNL_COMPLEX


#define BASE_INT
#include "pnl_templates_on.h"
#include "pnl_Tarray_source.c"
#include "pnl_templates_off.h"
#undef  BASE_INT

#define BASE_UINT
#include "pnl_templates_on.h"
#include "pnl_Tarray_source.c"
#include "pnl_templates_off.h"
#undef  BASE_UINT


