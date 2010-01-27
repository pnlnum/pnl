#ifndef  _SPECIAL_FUNC_H
#define _SPECIAL_FUNS_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

extern double pnl_sf_gamma_inc(double a,double x);
extern double pnl_sf_gamma_inc_P(double a,double x);
extern double pnl_sf_gamma_inc_Q(double a,double x);
extern double pnl_sf_expint_En(int n,double x);

extern double pnl_sf_gamma (double);
extern double pnl_sf_log_gamma (double);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif
