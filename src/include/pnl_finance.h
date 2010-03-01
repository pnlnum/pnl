#ifndef _FINANCE_FUNCTION_H
#define _FINANCE_FUNCTION_H


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl_matrix_uint.h"
#include "pnl_cdf.h"

/**
 * \defgroup Finance A few financial functions
 */
/*@{*/

extern int pnl_cf_call_bs(double s,double k,double t,double r,double divid,double sigma,
                          double *ptprice,double *ptdelta);
extern int pnl_cf_put_bs(double s,double k,double t,double r,double divid,double sigma,
                         double *ptprice,double *ptdelta);
extern double pnl_bs_call(double s,double k,double t,double r,double divid,double sigma);
extern double pnl_bs_put(double s,double k,double t,double r,double divid,double sigma);
extern double pnl_bs_call_put(int iscall, double s,double k,double t,double r,double divid,double sigma);
extern double pnl_bs_gamma(double s,double k,double t,double r,double divid,double sigma);
extern double pnl_bs_vega (double s, double k, double t, double r, double divid, double sigma);
double pnl_bs_implicit_vol (int is_call, double Price, double spot, double Strike,
                            double T, double r, double divid, int *error);
extern int pnl_bs_matrix_implicit_vol(const PnlMatInt * Is_Call, const PnlMat * Price,
                                      double spot,double rate, double divid,
                                      const PnlVect * Strike,const PnlVect * Maturity,PnlMat * Vol);
/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif
