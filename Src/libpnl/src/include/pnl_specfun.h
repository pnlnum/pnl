#ifndef  _SPECIAL_FUNC_H
#define _SPECIAL_FUNS_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "pnl_bessel.h"

extern double pnl_sf_gamma_inc(double a,double x);
extern double pnl_sf_gamma_inc_P(double a,double x);
extern double pnl_sf_gamma_inc_Q(double a,double x);
extern double pnl_sf_expint_En(int n,double x);
/* extern double pnl_sf_expint_Ei(double x); */

extern double pnl_sf_gamma (double);
extern double pnl_sf_log_gamma (double);
extern double pnl_sf_log_erf (double x);
extern double pnl_sf_log_erfc (double x);
extern double pnl_sf_erf (double x);
extern double pnl_sf_erfc (double x);

extern double pnl_sf_hyperg_2F1 (double a, double b, double c, double x);
extern double pnl_sf_hyperg_1F1 (double a, double b, double x);
extern double pnl_sf_hyperg_2F0 (double a, double b, double x);
extern double pnl_sf_hyperg_0F1 (double c, double x);


#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif
