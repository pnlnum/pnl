#ifndef _PNL_ROOT_H
#define _PNL_ROOT_H

#include "pnl_mathtools.h"

/**
 * \defgroup Root Searching for zeros
 */

/*@{*/ 

extern double pnl_root_brent (PnlFunc * F, double  x1, double  x2, double *tol);
extern int pnl_find_root (PnlFuncDFunc *Func, double x_min, double x_max,
                         double Tol,int N_Max,double * res);
extern int pnl_root_newton (PnlFuncDFunc *Func, double x0, double epsrel,
                            double epsabs, int N_max, double * res);
extern int pnl_root_bisection (PnlFunc *Func, double xmin, double xmax,
                               double epsrel, double espabs, int N_max, double *res);

/*@}*/

#endif /* _PNL_ROOT_H */
