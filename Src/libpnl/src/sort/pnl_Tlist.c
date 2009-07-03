#include <stdlib.h>

#include "config.h"
#include "pnl_list.h"

#define BASE_CONTAIN
#include "pnl_templates_on.h"
#include "pnl_Tlist_source.c"
#include "pnl_templates_off.h"
#undef BASE_CONTAIN

#define BASE_SPARSE_POINT
#include "pnl_templates_on.h"
#include "pnl_Tlist_source.c"
#include "pnl_templates_off.h"
#undef BASE_SPARSE_POINT


