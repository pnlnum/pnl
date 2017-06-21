#ifndef PNL_INTPOINTS_BFGS_H
#define PNL_INTPOINTS_BFGS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_sp_matrix.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


extern int pnl_optim_intpoints_bfgs_solve(PnlRnFuncR * func, PnlRnFuncRm *grad_func, PnlRnFuncRm * nl_constraints, PnlVect* lower_bounds, PnlVect* upper_bounds, PnlVect * x_input , double tolerance, int iter_max, int print_inner_steps, PnlVect *output);
extern int pnl_optim_linprog(const PnlVect *C, const PnlMat *A_ineq, const PnlVect *B_ineq, const PnlMat *A_eq, const PnlVect *B_eq, const PnlVect *x_min, const PnlVect *x_max, int debug, PnlVect *xopt, double *fobj_opt);
extern int pnl_optim_linprog_sp(const PnlSpMat *C, const PnlSpMat *A_ineq, const PnlVect *B_ineq, const PnlSpMat *A_eq, const PnlVect *B_eq, const PnlVectInt *index_min, const PnlVect *x_min, const PnlVectInt *index_max, const PnlVect *x_max, int debug, PnlVect *xopt, double *fobj_opt);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* PNL_INTPOINTS_BFGS_H */
