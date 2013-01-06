#include "pnl/pnl_f2c.h"
#include "pnl/pnl_specfun.h"

double derf_(doublereal *x)
{
return( pnl_sf_erf(*x) );
}
