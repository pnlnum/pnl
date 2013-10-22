#ifndef _PNL_ROOT_H
#define _PNL_ROOT_H

#include "pnl/pnl_mathtools.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/**
 * \defgroup Root Searching for zeros
 */

/*@{*/ 

extern double pnl_root_brent (PnlFunc * F, double  x1, double  x2, double *tol);
extern int pnl_root_newton_bisection (PnlFuncDFunc *Func, double x_min, double x_max,
                          double Tol,int N_Max,double * res);
extern int pnl_root_newton (PnlFuncDFunc *Func, double x0, double epsrel,
                            double epsabs, int N_max, double * res);
extern int pnl_root_bisection (PnlFunc *Func, double xmin, double xmax,
                               double epsrel, double espabs, int N_max, double *res);
extern int pnl_multiroot_newton (PnlRnFuncRnDFunc *func, const PnlVect *x0, double x_eps, 
                                 double fx_eps, int max_iter, int verbose, PnlVect *res);
extern int pnl_root_fsolve (PnlRnFuncRnDFunc *f, PnlVect *x, PnlVect *fx,
                            double xtol, int maxfev, int *nfev, PnlVect
                            *scale, int error_msg);
extern int pnl_root_fsolve_lsq (PnlRnFuncRmDFunc *f, PnlVect *x, int m,
                                PnlVect *fx,  double xtol, double ftol,
                                double gtol, int maxfev, int *nfev, PnlVect
                                *scale, int error_msg);

/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_ROOT_H */
