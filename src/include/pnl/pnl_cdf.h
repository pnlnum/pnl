#ifndef _PNL_CDF_H
#define _PNL_CDF_H

#include "pnl/pnl_mathtools.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**
 * \defgroup PnlCDF Cumulative Distribution Functions
 */
/*@{*/


void pnl_cdf_bet(int *which,double *p,double *q,double *x,double *y,
                 double *a,double *b,int *status,double *bound);
void pnl_cdf_bin(int *which,double *p,double *q,double *s,double *xn,
                 double *pr,double *ompr,int *status,double *bound);
void pnl_cdf_chi(int *which,double *p,double *q,double *x,double *df,
                 int *status,double *bound);
void pnl_cdf_chn(int *which,double *p,double *q,double *x,double *df,
                 double *pnonc,int *status,double *bound);
void pnl_cdf_f(int *which,double *p,double *q,double *f,double *dfn,
               double *dfd,int *status,double *bound);
void pnl_cdf_fnc(int *which,double *p,double *q,double *f,double *dfn,
                 double *dfd,double *phonc,int *status,double *bound);
void pnl_cdf_gam(int *which,double *p,double *q,double *x,double *shape,
                 double *scale,int *status,double *bound);
void pnl_cdf_nbn(int *which,double *p,double *q,double *s,double *xn,
                 double *pr,double *ompr,int *status,double *bound);
void pnl_cdf_nor(int *which,double *p,double *q,double *x,double *mean,
                 double *sd,int *status,double *bound);
void pnl_cdf_poi(int *which,double *p,double *q,double *s,double *xlam,
                 int *status,double *bound);
void pnl_cdf_t(int *which,double *p,double *q,double *t,double *df,
               int *status,double *bound);

/*
 * Simple useable cumulative functions
 */ 

double pnl_cdfchi2n(double x, double df, double ncparam);
void pnl_cdfbchi2n(double x, double nu, double lambda, double beta, double *P);
double cdf_nor(double);
double pnl_cdfnor(double);
double pnl_cdf2nor( double a, double b, double r);
double pnl_normal_density(double x);
double pnl_inv_cdfnor (double u);
  
/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_CDF_H */
