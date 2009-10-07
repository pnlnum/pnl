#ifndef  _FINANCE_FUNCTION_H
#define _FINANCE_FUNCTION_H


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl_matrix_uint.h"
#include "pnl_cdf.h"

/*  Compute delta forward because this quantity  is sold/bought to hedge option with forward*/
extern double pnl_forward_price(double Spot,double r, double divid, double Maturity);
extern double pnl_bs_call (double Vol,double Bond, double Forward, double Strike, double Maturity);
extern double pnl_bs_put  (double Vol,double Bond, double Forward, double Strike, double Maturity);
extern double pnl_bs_call_delta_forward (double Vol,double Bond, double Forward,
                                         double Strike, double Maturity); 
extern double pnl_bs_put_delta_forward  (double Vol,double Bond, double Forward,
                                         double Strike, double Maturity); 
extern double pnl_bs_call_put (int Is_Call, double Vol,double Bond, double Forward,
                               double Strike, double Maturity);
extern double pnl_bs_call_put_delta_forward (int Is_Call, double Vol,double Bond,
                                             double Forward, double Strike, double Maturity);
extern double pnl_bs_vega(double Vol,double Bond, double Forward, double Strike, double Maturity);
extern double pnl_bs_gamma(double Vol,double Bond, double Forward, double Strike, double Maturity);
extern double pnl_bs_s_square_gamma (double Vol,double Bond, double Forward, double Strike,double Maturity);
extern double pnl_bs_implicit_vol(int Is_Call, double Price,double Bond,
                                  double Forward, double Strike, double Maturity);  
extern int pnl_bs_matrix_implicit_vol(const PnlMatUint * Is_Call, const PnlMat * Price,
                                      double spot,double rate, double divid,
                                      const PnlVect * Strike,const PnlVect * Maturity,PnlMat * Vol);

extern int pnl_cf_call_bs(double s,double k,double t,double r,double divid,double sigma,
                          double *ptprice,double *ptdelta);
extern int pnl_cf_put_bs(double s,double k,double t,double r,double divid,double sigma,
                         double *ptprice,double *ptdelta);
extern double pnl_cf_bs_gamma(double s,double k,double t,double r,double divid,double sigma);


#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif
