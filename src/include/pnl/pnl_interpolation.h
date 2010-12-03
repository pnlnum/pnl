#ifndef _PNL_INTERPOLATION_H
#define _PNL_INTERPOLATION_H

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

enum 
  {
    NOT_A_KNOT, 
    NATURAL, 
    CLAMPED, 
    PERIODIC, 
    FAST, 
    FAST_PERIODIC,
    MONOTONE, 
    BY_ZERO, 
    C0, 
    LINEAR, 
    BY_NAN, 
    UNDEFINED
  };

int pnl_bicubic_spline(PnlVect *x, PnlVect *y, PnlMat *u, double *C, int type);

void pnl_eval_bicubic(PnlVect *x, PnlVect *y, double *C, PnlVect *x_eval,
                      PnlVect *y_eval, PnlMat *z_eval, PnlMat *dzdx_eval,
                      PnlMat *dzdy_eval, int outmode);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _PNL_INTERPOLATION_H */
