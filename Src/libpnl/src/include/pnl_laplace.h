#ifndef PNL_LAPLACE
#define PNL_LAPLACE


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl_mathtools.h"
#include "pnl_complex.h"
#include "pnl_vector.h"

extern double pnl_ilap_cdf_euler(PnlCmplxFunc *f, double t, double h, int N, int M);
extern double pnl_ilap_euler(PnlCmplxFunc *f, double t, int N, int M);
extern void pnl_ilap_fft(PnlVect *res, PnlCmplxFunc *f, double T, double eps);
extern double pnl_ilap_gs_basic (PnlFunc *fhat, double t, int n);
extern double pnl_ilap_gs (PnlFunc *fhat, double t, int n);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* PNL_LAPLACE */
