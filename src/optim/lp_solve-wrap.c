#include <stdio.h>
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_sp_matrix.h"
#include "pnl/pnl_mathtools.h"
#include "lp_solve/lp_lib.h"

/**
 * Solve a linear programming optimization problem calling LPSolve.
 *
 *      min C^T x
 *      s.t.
 *          A_ineq x <= B_ineq
 *          A_eq x = B_eq
 *          x_min <= x <= x_max
 *
 * @param C defines the linear objective fucntion
 * @param A_ineq The matrix of inequality constraints. Can be NULL, if there is no such constraint.
 * @param B_ineq Right hand side member of the inequality constraints. The length of B_ineq must
 * equal the number of rows of A_ineq.
 * @param A_eq The matrix of equality constraints. Can be NULL, if there is no such constraint.
 * @param B_eq Right hand side member of the equality constraints. The length of B_eq must
 * equal the number of rows of A_eq.
 * @param x_min The lower bounds on the solution. If NULL, all the variales are supposed to be
 * non negative. Use -PNL_INF for every coordinate with no lower bound.
 * @param x_max The upper bounds on the solution. If NULL, it means +Infinity.
 * Use PNL_INF for every coordinate with no upper bound.
 * @param debug TRUE or FALSE. If TRUE, it prints all intermediate results and the branch-and-bound decisions
 * @param[out] xopt contains the optimal x on output.
 * @param[out] fobj_opt contains the value of the objective function at xopt
 *
 * @return FAIL or OK
 */
int pnl_optim_linprog(const PnlVect *C, const PnlMat *A_ineq, const PnlVect *B_ineq, const PnlMat *A_eq, const PnlVect *B_eq,
                      const PnlVect *x_min, const PnlVect *x_max, int debug, PnlVect *xopt, double *fobj_opt)
{
  lprec *lp;
  int i, ret;
  int nVar;
  int nIneqConstraint, nEqConstraint, nConstraint;
#ifndef PNL_RANGE_CHECK_OFF
  if (C == NULL)
    {
      perror("Objective function not defined.\n");
      abort();
    }
  if (((A_ineq != NULL) && (B_ineq != NULL) && ((A_ineq->n != C->size) || (A_ineq->m != B_ineq->size )))
      || ((A_ineq == NULL) && (B_ineq != NULL)) || ((A_ineq == NULL) && (B_ineq != NULL)))
    {
      perror("Non-compatible inequality constraints.\n");
      abort();
    }
  if (((A_eq != NULL) && (B_eq != NULL) && ((A_eq->n != C->size) || (A_eq->m != B_eq->size )))
      || ((A_eq == NULL) && (B_eq != NULL)) || ((A_eq == NULL) && (B_eq != NULL)))
    {
      perror("Non-compatible equality constraints.\n");
      abort();
    }
  if (((x_min != NULL) && (x_min->size != C->size)) || ((x_max != NULL) && (x_max->size != C->size)))
    {
      perror("Non-compatible bounds.\n");
      abort();
    }
#endif

  nVar = C->size;
  nIneqConstraint = nEqConstraint = 0;
  if (A_ineq != NULL) nIneqConstraint = A_ineq->m;
  if (A_eq != NULL) nEqConstraint = A_eq->m;
  nConstraint = nIneqConstraint + nEqConstraint;
  if ((lp = make_lp(nConstraint, nVar)) == NULL)
    {
      printf("Cannot allocate a linear programming problem.\n");
      return FAIL;
    }
  set_verbose(lp, IMPORTANT);
  if (debug == TRUE) set_debug(lp, TRUE);
  /* Set the objective function */
  set_obj_fnex(lp, C->size, C->array - 1, NULL);
  /* Set the inequality constraints */
  for (i = 0; i < nIneqConstraint; i++)
    {
      set_rowex(lp, i + 1, A_ineq->n, A_ineq->array + i * A_ineq->n - 1, NULL);
      set_constr_type(lp, i + 1, LE);
      set_rh(lp, i + 1, GET(B_ineq, i));
    }
  /* Set the equality constraints */
  for (i = 0; i < nEqConstraint; i++)
    {
      set_rowex(lp, nIneqConstraint + i + 1, A_eq->n,  A_eq->array + i * A_eq->n - 1, NULL);
      set_constr_type(lp, nIneqConstraint + i + 1, EQ);
      set_rh(lp, nIneqConstraint + i + 1, GET(B_eq, i));
    }
  /* Set the lower bounds. The default lower bound is 0 */
  if (x_min != NULL)
    {
      for (i = 0; i < x_min->size; i++)
        {
          double bound_i = GET(x_min, i);
          if (bound_i == 0) continue;
          if (pnl_isinf(bound_i)) bound_i = -DEF_INFINITE;
          set_lowbo(lp, i + 1, bound_i);
        }
    }
  /* Set the upper bounds. The Default upper bound is +Inf */
  if (x_max != NULL)
    {
      for (i = 0; i < x_max->size; i++)
        {
          double bound_i = GET(x_max, i);
          if (pnl_isinf(bound_i)) continue;
          set_upbo(lp, i + 1, bound_i);
        }
    }

  if (debug) print_lp(lp);
  ret = solve(lp);
  switch(ret)
    {
    case NOMEMORY:
      perror("LPSolve ran out of memory.\n"); break;
    case INFEASIBLE:
      perror("The model is infeasible.\n"); break;
    case UNBOUNDED:
      perror("The model is unbounded.\n"); break;
    case DEGENERATE:
      perror("The model is degenerative.\n"); break;
    case NUMFAILURE:
      perror("A numerical failure was encoutered.\n"); break;
    case OPTIMAL:
      pnl_vect_resize(xopt, nVar);
      get_variables(lp, xopt->array);
      *fobj_opt = get_objective(lp);
      break;
    }

  delete_lp(lp);
  return (ret == OPTIMAL) ? OK : FAIL;
}

int * incr(int n, int *T)
{
  int i;
  int *ret = malloc(n * sizeof(int));
  for (i = 0; i < n; i++) { ret[i] = T[i] + 1; }
  return ret;
}


/**
 * Solve a linear programming optimization problem calling LPSolve.
 *
 *      min C^T x
 *      s.t.
 *          A_ineq x <= B_ineq
 *          A_eq x = B_eq
 *          x_min <= x <= x_max
 *
 * @param C defines the linear objective fucntion
 * @param A_ineq The matrix of inequality constraints. Can be NULL, if there is no such constraint.
 * @param B_ineq Right hand side member of the inequality constraints. The length of B_ineq must
 * equal the number of rows of A_ineq.
 * @param A_eq The matrix of equality constraints. Can be NULL, if there is no such constraint.
 * @param B_eq Right hand side member of the equality constraints. The length of B_eq must
 * equal the number of rows of A_eq.
 * @param x_min The lower bounds on the solution. If NULL, all the variales are supposed to be
 * non negative. Use -PNL_INF for every coordinate with no lower bound.
 * @param x_max The upper bounds on the solution. If NULL, it means +Infinity.
 * Use PNL_INF for every coordinate with no upper bound.
 * @param index_min the indices of the variables with lower bounds specified by x_min. Must be NULL
 * when x_min is NULL.
 * @param index_max the indices of the variables with lower bounds specified by x_min. Must be NULL
 * when x_max is NULL.
 * @param debug TRUE or FALSE. If TRUE, it prints all intermediate results and the branch-and-bound decisions
 * @param[out] xopt contains the optimal x on output.
 * @param[out] fobj_opt contains the value of the objective function at xopt
 *
 * @return FAIL or OK
 */
int pnl_optim_linprog_sp(const PnlSpMat *C, const PnlSpMat *A_ineq, const PnlVect *B_ineq, const PnlSpMat *A_eq, const PnlVect *B_eq,
                         const PnlVectInt *index_min, const PnlVect *x_min, const PnlVectInt *index_max, const PnlVect *x_max,
                         int debug, PnlVect *xopt, double *fobj_opt)
{
  lprec *lp;
  int i, ret;
  int nVar;
  int *colC, *colA_ineq, *colA_eq;
  int nIneqConstraint, nEqConstraint, nConstraint;
#ifndef PNL_RANGE_CHECK_OFF
  if (C == NULL)
    {
      perror("Objective function not defined.\n");
      abort();
    }
  if (((A_ineq != NULL) && (B_ineq != NULL) && ((A_ineq->n != C->m) || (A_ineq->m != B_ineq->size )))
      || ((A_ineq == NULL) && (B_ineq != NULL)) || ((A_ineq == NULL) && (B_ineq != NULL)))
    {
      perror("Non-compatible inequality constraints.\n");
      abort();
    }
  if (((A_eq != NULL) && (B_eq != NULL) && ((A_eq->n != C->m) || (A_eq->m != B_eq->size )))
      || ((A_eq == NULL) && (B_eq != NULL)) || ((A_eq == NULL) && (B_eq != NULL)))
    {
      perror("Non-compatible equality constraints.\n");
      abort();
    }
  if (((x_min != NULL) && ((x_min->size > C->m) || (x_min->size != index_min->size))) || ((x_max != NULL) && ((x_max->size > C->m) || (x_max->size != index_max->size))))
    {
      perror("Non-compatible bounds.\n");
      abort();
    }
#endif

  nVar = C->m;
  nIneqConstraint = nEqConstraint = 0;
  if (A_ineq != NULL) nIneqConstraint = A_ineq->m;
  if (A_eq != NULL) nEqConstraint = A_eq->m;
  nConstraint = nIneqConstraint + nEqConstraint;
  colC = incr(C->nz, C->I);
  colA_ineq = incr(A_ineq->nz, A_ineq->J);
  colA_eq = incr(A_eq->nz, A_eq->J);
  if ((lp = make_lp(nConstraint, nVar)) == NULL)
    {
      printf("Cannot allocate a linear programming problem.\n");
      return FAIL;
    }
  set_verbose(lp, IMPORTANT);
  if (debug == TRUE) set_debug(lp, TRUE);
  /* Set the objective function */
  set_obj_fnex(lp, C->nz, C->array, colC);
  /* Set the inequality constraints */
  for (i = 0; i < nIneqConstraint; i++)
    {
      int offset = A_ineq->I[i];
      int nElements = A_ineq->I[i + 1] - offset;
      set_rowex(lp, i + 1, nElements, A_ineq->array + offset, colA_ineq + offset);
      set_constr_type(lp, i + 1, LE);
      set_rh(lp, i + 1, GET(B_ineq, i));
    }
  /* Set the equality constraints */
  for (i = 0; i < nEqConstraint; i++)
    {
      int offset = A_eq->I[i];
      int nElements = A_eq->I[i + 1] - offset;
      set_rowex(lp, nIneqConstraint + i + 1, nElements, A_eq->array + offset, colA_eq + offset);
      set_constr_type(lp, nIneqConstraint + i + 1, EQ);
      set_rh(lp, nIneqConstraint + i + 1, GET(B_eq, i));
    }
  /* Set the lower bounds. The default lower bound is 0 */
  if (x_min != NULL)
    {
      for (i = 0; i < x_min->size; i++)
        {
          double bound_i = GET(x_min, i);
          if (bound_i == 0) continue;
          int index_i = GET_INT(index_min, i);
          if (pnl_isinf(bound_i)) bound_i = -DEF_INFINITE;
          set_lowbo(lp, index_i + 1, bound_i);
        }
    }
  /* Set the upper bounds. The Default upper bound is +Inf */
  if (x_max != NULL)
    {
      for (i = 0; i < x_max->size; i++)
        {
          double bound_i = GET(x_max, i);
          if (bound_i == 0) continue;
          int index_i = GET_INT(index_max, i);
          set_upbo(lp, index_i + 1, bound_i);
        }
    }

  if (debug) print_lp(lp);
  ret = solve(lp);
  switch(ret)
    {
    case NOMEMORY:
      perror("LPSolve ran out of memory.\n"); break;
    case INFEASIBLE:
      perror("The model is infeasible.\n"); break;
    case UNBOUNDED:
      perror("The model is unbounded.\n"); break;
    case DEGENERATE:
      perror("The model is degenerative.\n"); break;
    case NUMFAILURE:
      perror("A numerical failure was encoutered.\n"); break;
    case OPTIMAL:
      pnl_vect_resize(xopt, nVar);
      get_variables(lp, xopt->array);
      *fobj_opt = get_objective(lp);
      break;
    }

  delete_lp(lp);
  free(colC);
  free(colA_eq);
  free(colA_ineq);
  return (ret == OPTIMAL) ? OK : FAIL;
}
