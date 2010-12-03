#ifndef _PNL_LAPLACE_H
#define _PNL_LAPLACE_H


#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_complex.h"
#include "pnl/pnl_vector.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**
 * \defgroup Laplace Laplace transforms
 *
 * Numercial inversion of Laplace transforms either on the real axis or
 * in the complex plane
 */
/*@{*/
extern double pnl_ilap_cdf_euler(PnlCmplxFunc *f, double t, double h, int N, int M);
extern double pnl_ilap_euler(PnlCmplxFunc *f, double t, int N, int M);
extern void pnl_ilap_fft(PnlVect *res, PnlCmplxFunc *f, double T, double eps);
extern double pnl_ilap_gs_basic (PnlFunc *fhat, double t, int n);
extern double pnl_ilap_gs (PnlFunc *fhat, double t, int n);
/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_LAPLACE_H */
