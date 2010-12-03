#ifndef _PNL_FINANCE_H
#define _PNL_FINANCE_H


#include "pnl/pnl_matrix.h"
#include "pnl/pnl_cdf.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

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


#endif /* _PNL_FINANCE_H */
