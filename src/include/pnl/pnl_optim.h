#ifndef PNL_INTPOINTS_BFGS_H
#define PNL_INTPOINTS_BFGS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


extern int pnl_optim_intpoints_bfgs_solve(PnlRnFuncR * func, PnlRnFuncRm *grad_func, PnlRnFuncRm * nl_constraints, PnlVect* lower_bounds, PnlVect* upper_bounds, PnlVect * x_input , double tolerance, int iter_max, int print_inner_steps, PnlVect *output);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* PNL_INTPOINTS_BFGS_H */
