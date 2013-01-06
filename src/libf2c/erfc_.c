#include "pnl/pnl_f2c.h"
#include "pnl/pnl_specfun.h"

double erfc_(real *x)
{
return( pnl_sf_erfc(*x) );
}
