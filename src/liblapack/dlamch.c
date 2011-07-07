#include "pnl/pnl_f2c.h"

/*
 * Turn off compiler optimization for Win32.
 */
#if defined _WIN32
#pragma optimize ( "", off )
#endif
#include "pnl/pnl_internals.h"

/*
 * The code of Lapack's dlamch is already in libamos, so we it if needed 
 */
double dlamch_(char *cmach)
{
  return pnl_dlamch(cmach);
}


