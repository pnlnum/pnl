#include "pnl/pnl_f2c.h"
#include "pnl/pnl_specfun.h"

double derfc_(doublereal *x)
{
return( pnl_sf_erfc(*x) );
}
