/************************************************************************/
/* Copyright Ismail Laachir <ismail.laachir@gmail.fr> 2009.		          */
/*                                                                      */
/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as       */
/* published by the Free Software Foundation, either version 3 of the   */
/* License, or (at your option) any later version.                      */
/*                                                                      */
/* This program is distributed in the hope that it will be useful, but  */
/* WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    */
/* Lesser General Public License for more details.                      */
/*                                                                      */
/* You should have received a copy of the GNU Lesser General Public     */
/* License along with this program.  If not, see                        */
/* <http://www.gnu.org/licenses/>.                                      */
/************************************************************************/

/************************************************************************/
/* This program was originally written in C#, it was translated into C  */
/* by David Pommier <pommier.david@gmail.com> 2009.			                */
/************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

#include "pnl/pnl_optim.h"


/**
 * Structure containing the inequality constraints.
 */
typedef struct AllConstraints
{
  PnlRnFuncRm* NL_Constraints; // non linear constraints
  PnlVect* LowerBounds; // lower bound constraints
  PnlVect* UpperBounds; // upper bound constraints

  PnlVectInt* LowerBoundsIndex; // index of variable where a bound constraints is applied
  PnlVectInt* UpperBoundsIndex; // index of variable where a bound constraints is applied

  int nbr_nl_constraints;
  int nbr_lower_bounds;
  int nbr_upper_bounds;

} AllConstraints;

/**
 * Frees the structure AllConstraints
 *
 * @param all_constraints
 */
static void AllConstraints_Free(AllConstraints *all_constraints)
{
  pnl_vect_free(&(all_constraints->LowerBounds));
  pnl_vect_free(&(all_constraints->UpperBounds));

  pnl_vect_int_free(&(all_constraints->LowerBoundsIndex));
  pnl_vect_int_free(&(all_constraints->UpperBoundsIndex));
}

/**
 * Initialize all_constraints with:
 *
 * @param all_constraints
 * @param NL_Constraints: the PnlRnFuncRm that computes a vector of nonlinear constraints
 * @param LowerBounds: the lower bound constraints. LowerBounds[i] may be
 * -inf if no lower bound constraint is applied on x[i]. If LowerBound ==
 *  NULL, it means there is no lower bound constraint
 * @param UpperBounds: the upper bound constraints. UpperBounds[i] may be
 * +inf if no upper bound constraint is applied on x[i].  If UpperBound ==
 *  NULL, it means there is no upper bound constraint.
 * @param xinit: the initial value of x. Is should be strictly feasible. ie
 * constraint(xinit)>0.
 * @return xinit_feasible: integer such that (xinit_feasible=0) if
 * (c(xinit)>0), 1 otherwise.
 */
static int AllConstraints_Init(AllConstraints *all_constraints,
                               PnlRnFuncRm* NL_Constraints, PnlVect*
                               LowerBounds, PnlVect* UpperBounds, PnlVect*
                               xinit)
{
  int i, j, k, nbr_var, xinit_feasible;
  double lb, ub;

  PnlVect* dummy;

  dummy = pnl_vect_create(0);
  all_constraints->NL_Constraints = NL_Constraints;

  nbr_var = xinit->size;
  all_constraints->NL_Constraints = NL_Constraints;
  all_constraints->LowerBoundsIndex = pnl_vect_int_create(nbr_var);
  all_constraints->UpperBoundsIndex = pnl_vect_int_create(nbr_var);
  all_constraints->LowerBounds = pnl_vect_create(nbr_var);
  all_constraints->UpperBounds = pnl_vect_create(nbr_var);

  xinit_feasible=0;
  if (all_constraints->NL_Constraints!=NULL)
    {
      PNL_EVAL_RNFUNCRM(NL_Constraints, xinit, dummy);
      all_constraints->nbr_nl_constraints = dummy->size;
      for (i=0; i<all_constraints->nbr_nl_constraints; i++)
        {
          if (GET(dummy, i)<=0)
            {
              xinit_feasible = 1;
              break;
            }
        }
    }
  else all_constraints->nbr_nl_constraints = 0;

  k=0;
  j=0;
  if ( LowerBounds != NULL )
    {
      for (i=0; i<nbr_var; i++)
        {
          lb = GET(LowerBounds, i);
          if (pnl_isinf(lb)==0)
            {
              pnl_vect_int_set(all_constraints->LowerBoundsIndex, j, i);
              LET(all_constraints->LowerBounds, j) = lb;
              j++;

              if (GET(xinit, i)<=lb) xinit_feasible = 1;
            }
        } 
    }
  if ( UpperBounds != NULL )
    {
      for (i=0; i<nbr_var; i++)
        {
          ub = GET(UpperBounds, i);
          if (pnl_isinf(ub)==0)
            {
              pnl_vect_int_set(all_constraints->UpperBoundsIndex, k, i);
              LET(all_constraints->UpperBounds, k) = ub;
              k++;

              if (GET(xinit, i)>=ub) xinit_feasible = 1;
            }
        }
    }
  all_constraints->nbr_lower_bounds = j;
  all_constraints->nbr_upper_bounds = k;

  pnl_vect_int_resize(all_constraints->LowerBoundsIndex, j);
  pnl_vect_int_resize(all_constraints->UpperBoundsIndex, k);
  pnl_vect_resize(all_constraints->LowerBounds, j);
  pnl_vect_resize(all_constraints->UpperBounds, k);

  pnl_vect_free(&dummy);
  return xinit_feasible;
}

/**
 * Computes the value of constraints at vector x
 *
 * @param constraints_x: vector containting value of constraints valued at x
 * @param x: vector where constraints are valued
 * @param all_constraints: AllConstraints structure
 */
static void AllConstraints_Value(PnlVect *constraints_x, const PnlVect *x,
                                 AllConstraints *all_constraints)
{
  int i, index;

  if (all_constraints->NL_Constraints!=NULL) 
    PNL_EVAL_RNFUNCRM(all_constraints->NL_Constraints, x, constraints_x);

  pnl_vect_resize(constraints_x, all_constraints->nbr_nl_constraints 
                                  + all_constraints->nbr_lower_bounds 
                                  + all_constraints->nbr_upper_bounds);

  for (i=0; i<(all_constraints->nbr_lower_bounds); i++)
    {
      index = PNL_GET(all_constraints->LowerBoundsIndex, i);
      LET(constraints_x, all_constraints->nbr_nl_constraints + i) = GET(x, index)-GET(all_constraints->LowerBounds, i);
    }

  for (i=0; i<(all_constraints->nbr_upper_bounds); i++)
    {
      index = PNL_GET(all_constraints->UpperBoundsIndex, i);
      LET(constraints_x, all_constraints->nbr_nl_constraints 
                        + all_constraints->nbr_lower_bounds + i) 
           = GET(all_constraints->UpperBounds, i) - GET(x, index);
    }
}


/**
 * Computes the value of gradient of the function f.
 *
 * @param grad_func: function computing the gradient. It may be NULL, in
 * this case the gradient of objective function is estimated with finites
 * differences.
 * @param func: function to be differentiated. It's callled in finites
 * differences method when grad_func is not given.
 * @param x: vector where gradient is computed
 * @param res: vector containting value of gradient valued at x
 */
static void gradient_f(PnlRnFuncRm *grad_func, PnlRnFuncR * func, PnlVect * x, PnlVect *res)
{
  int i=0,nbr_var = x->size;
  double h = 0.00001;
  double dx=0;

  /* If the gradient of objective function is not provided, we estimate it
   * using finites differences. */
  if (grad_func==NULL)
    {

      PnlVect * x1=pnl_vect_copy(x);
      PnlVect * x2=pnl_vect_copy(x);

      pnl_vect_resize(res, nbr_var);

      for (i = 0; i < nbr_var; i++)
        {
          if (fabs(GET(x,i))<h) dx=h;
          else dx = h*GET(x,i);

          /* Shift x[i] */
          LET(x1,i)=GET(x,i) + dx;
          LET(x2,i)=GET(x,i) - dx;
          LET(res,i) = (PNL_EVAL_RNFUNCR(func, x1) - PNL_EVAL_RNFUNCR(func, x2)) / (2 * dx);
          /* Restore x[i] */
          LET(x1,i)=GET(x,i);
          LET(x2,i)=GET(x,i);

        }
      pnl_vect_free(&x1);
      pnl_vect_free(&x2);
    }
  else
    {
      PNL_EVAL_RNFUNCRM(grad_func, x, res);
    }
}


/**
 * Computes the value of gradient of constraints.
 *
 * @param all_constraints: AllConstraints structure
 * @param x: vector where gradient is computed
 * @param grad_c: vector containting value of gradient valued at x
 */
static void gradient_c(AllConstraints* all_constraints, PnlVect* x, PnlMat* grad_c)
{
  int i, j, nbr_nl_constraints, nbr_lower_bounds, nbr_upper_bounds, nbr_var = x->size;
  double h = 0.0001;
  double dx;
  PnlVect * x1 = pnl_vect_copy(x);
  PnlVect * x2 = pnl_vect_copy(x);
  PnlVect *constraintsmd = pnl_vect_create(0);
  PnlVect *constraintspd = pnl_vect_create(0);

  nbr_nl_constraints = all_constraints->nbr_nl_constraints;
  nbr_lower_bounds = all_constraints->nbr_lower_bounds;
  nbr_upper_bounds = all_constraints->nbr_upper_bounds;

  pnl_mat_resize(grad_c, nbr_var, nbr_nl_constraints+nbr_lower_bounds+nbr_upper_bounds);
  pnl_mat_set_all(grad_c, 0.);

  if (all_constraints->NL_Constraints!=NULL)
    {
      for (i=0; i<nbr_var; i++)
        {
          if (fabs(GET(x,i))<h) dx=h;
          else dx = h*GET(x,i);

          LET(x1,i)=GET(x,i) + dx;
          LET(x2,i)=GET(x,i) - dx;

          PNL_EVAL_RNFUNCRM(all_constraints->NL_Constraints, x1, constraintspd);
          PNL_EVAL_RNFUNCRM(all_constraints->NL_Constraints, x2, constraintsmd);
          pnl_vect_axpby(-1/(2*dx), constraintsmd, 1/(2*dx), constraintspd);

          for (j=0; j<nbr_nl_constraints; j++)
            {
              MLET(grad_c, i, j) = GET(constraintspd, j);
            }

          LET(x1,i)=GET(x,i);
          LET(x2,i)=GET(x,i);
        }
    }

  for (i=0; i<nbr_lower_bounds; i++)
    {
      j = pnl_vect_int_get(all_constraints->LowerBoundsIndex, i);
      MLET(grad_c, j, nbr_nl_constraints+i) = 1.;
    }

  for (i=0; i<nbr_upper_bounds; i++)
    {
      j = pnl_vect_int_get(all_constraints->UpperBoundsIndex, i);
      MLET(grad_c, j, nbr_nl_constraints+nbr_lower_bounds+i) = -1.;
    }

  pnl_vect_free(&constraintsmd);
  pnl_vect_free(&constraintspd);
  pnl_vect_free(&x1);
  pnl_vect_free(&x2);
}


/**********************************************************************************
  Solving inequality constrained optimization problems

 ***********************
 * min {f(x) ; c(x)>=0}*    (P)
 ***********************

 This Algorithm was proposed by P.ARMAND, J-C.GILBERT, and S.JAN-JEGOU in their article
 "A feasible BFGS interior point algorithm for solving strongly convex minimization problems."
 SIAM J. Optim. (2000) 11, 199 222"
 For sake of simplicity, we used the same names for the variables in this implentation.

 Input parameters:
 func              ->  Function to minimize
 grad_func         ->  Gradient of func
 nl_constraints    ->  Nonlinear inequality constraints
 lower_bounds      ->  Lower bound constraints
 upper_bounds      ->  Upper bound constraints
 x_input           ->  Initial point wher the algorithm starts
 tolerance         ->  Precision required in solving (P)
 iter_max          ->  Maximum number of iterations in the algorithm.
 print_algo_steps  ->  Flag to decide to print information.
 =0 : we print no information
 =1 : we print final results
 >1 : we print information at each iteration

 Output parameters:
 x_output         ->   Output point
 CONVERGENCE ->   Convergence state

 !!!!!!!!!!!! IMPORTANT : Initial point must be strictly feasible. ie constraint(x) > 0.  !!!!!!!!!!!!
 !!!!! Do not choose x as constraint(x)=0 !!!!!

 CONVERGENCE=0 : Output Status: Failure: Initial point is not strictly feasible
 CONVERGENCE=1 : Output Status: Step too small, we stop the algorithm.
 CONVERGENCE=2 : Output Status: Maximum iteration reached.
 CONVERGENCE=3 : Output Status: A solution has been found up to the required accuracy.
 ****************************************************************************************/
/** 
 * 
 * min {f(x) ; c(x)>=0} 
 * 
 * @param func the function to minimize
 * @param grad_func the gradient of func. It may be NULL, in this case a
 * finite difference approach is used
 * @param nl_constraints the constraint function with the convention c(x) * >= 0
 * @param lower_bounds Lower bound constraints. If NULL, no lower bound
 * constraint. Some components can be -Inf
 * @param upper_bounds Upper bound constraints. If NULL, no upper ound
 * constraint. Some components can be +Inf
 * @param x_input Initial starting point of the algorithm
 * @param tolerance Precision required in solving the minimization problem
 * @param iter_max Maximum number of iterations 
 * @param print_algo_steps Flag to print extra information
 *               =0 : we print no information
 *               =1 : we print final results
 *               >1 : we print information at each iteration
 * @param x_output the solution
 * 
 * @return an integer error code
 *    =0 : Output Status: Failure: Initial point is not strictly feasible
 *    =1 : Output Status: Step too small, we stop the algorithm.
 *    =2 : Output Status: Maximum iteration reached.
 *    =3 : Output Status: A solution has been found up to the required accuracy.
 */
int pnl_optim_intpoints_bfgs_solve(PnlRnFuncR * func, PnlRnFuncRm
                                   *grad_func, PnlRnFuncRm *
                                   nl_constraints, PnlVect* lower_bounds,
                                   PnlVect* upper_bounds, PnlVect *
                                   x_input, double tolerance, int iter_max,
                                   int print_algo_steps, PnlVect *
                                   x_output)
{
  int i, j, k, nbr_func_eval=0;
  int CONVERGENCE; //Output Status
  double norm_grad_lagrangien=0.; // Norm of the gradient of langrange function.
  double norm_comp_cond=0.; // Norm of complmentarity condition.
  double norm_delta_x_min = 1e-15; // Minimal step that the algorithm can make.
  double norm_delta_x = 0.0; // Step between successive iterates
  double norm_x, norm_d_x; // Norm of current iterate. Norm of the descent direction.

  int nbr_var = x_input->size; // Number of variables
  int nbr_nl_constraints, nbr_lower_bounds, nbr_upper_bounds, nbr_constraints; // nbr of constraints.
  int nbr_iterations=0; // Current number of iterations
  int TestInfeasibility=0; // Test of feasibility of (x, lambda)
  int OptimalityCriterion=0., inner_iter=0, inner_iter_max=5;

  double f;     // Current value of objective function
  double f_min=0.; // Minimum value reached during the algorithm
  double mu = 1., mu_min=1e-30; // Perturbation parameter. Minimum value for the perturbation parameter.

  //****** Parameters of Line Search ******//
  double alpha = 1.0;  // Descent step size. x=x+alpha*d_x and lambda=lambda+alpha*d_lambda
  double decrease_alpha = 0.5; // Decreasing rate of alpha in Line Search step
  double omega = 1e-4; // omega in [0,1/2], used in Line Search
  double grad_psi_d; //Scalar product between gradient of merit function psi and primal descent direction d_x: <grad_psi,d_x>
  double psi;          // psi : value of merit function used in Line Search
  double psi_alpha;    // psi_alpha : value of merit function used in Line Search

  //****** Parameters of BFGS formula ******//
  double kappa = 0.2;     // Parameter used in Powell's correction
  double eta = 1e-50;     // Parameter used in Powell's correction
  double theta;           // Parameter used in Powell's correction
  double delta_M_delta;   // <delta,M*delta>
  double delta_gamma;     // <delta,gamma>
  double norm_gamma;      // <gamma,gamma>
  double SelfScalingBFGS = 1;

  double sum = 0.0, temps_cpu;

  PnlVect *comp_cond_vect, *grad_lagrangien, *vect_constraints, *lambda,*lambda_old, *x, *x_old,* d_x,* d_lambda, * delta,*M_delta,* delta_M,* gamma,*grad_f,* grad_f_old;
  PnlMat * M,* M_old,* A,* grad_c,* grad_c_old;
  PnlVect *x_min, *lambda_min;

  AllConstraints all_constraints;

  clock_t temps_initial, temps_final;
  temps_initial = clock();

  TestInfeasibility = AllConstraints_Init(&all_constraints, nl_constraints, lower_bounds, upper_bounds, x_input);
  nbr_nl_constraints = (all_constraints.nbr_nl_constraints);
  nbr_lower_bounds = (all_constraints.nbr_lower_bounds);
  nbr_upper_bounds = (all_constraints.nbr_upper_bounds);
  nbr_constraints = nbr_nl_constraints + nbr_lower_bounds + nbr_upper_bounds; // Total number of constraints

  //Create vectors and matrices used in the algorithm
  x = pnl_vect_copy(x_input); // current primal variable
  x_old = pnl_vect_create(nbr_var);// previous primal variable
  x_min = pnl_vect_create(nbr_var); // primal variable corresponding to minimal value reached during the algorithm
  grad_lagrangien = pnl_vect_create(nbr_var); // gradient of the lagrangian
  vect_constraints = pnl_vect_create(nbr_constraints); // current value of the constraints
  comp_cond_vect = pnl_vect_create(nbr_constraints); // term by term multiplication of constraints(x) and lambda
  lambda = pnl_vect_create(nbr_constraints); // current dual variable
  lambda_old = pnl_vect_create(nbr_constraints); // previous dual variable
  lambda_min = pnl_vect_create(nbr_constraints); // dual variable corresponding to minimal value reached during the algorithm
  d_x = pnl_vect_create(nbr_var); // Primal descent direction
  d_lambda = pnl_vect_create(nbr_constraints); // Dual descent direction
  delta = pnl_vect_create(nbr_var); // x(k+1) - x(k) = alpha*d_x
  M_delta = pnl_vect_create(nbr_var); // matrix*vector
  delta_M = pnl_vect_create(nbr_var); // vector*matrix
  gamma = pnl_vect_create(nbr_var); // grad_lagrangian(x(k+1), lambda)-grad_lagrangian(x(k), lambda)
  grad_f = pnl_vect_create(nbr_var); // grandient of objective function
  grad_f_old = pnl_vect_create(nbr_var); // previous value of grandient of objective function
  M = pnl_mat_create(nbr_var, nbr_var); // BFGS matrix
  M_old = pnl_mat_create(nbr_var, nbr_var); // previous BFGS matrix
  A = pnl_mat_create(nbr_var,nbr_var); // Matrix used to compute the descent direction d_x, by solving a linear system
  grad_c = pnl_mat_create(nbr_var,nbr_constraints); // gradient of constraints
  grad_c_old = pnl_mat_create(nbr_var,nbr_constraints); //previous value of gradient of constraints

  // Initial point must be strictly admissible. ie constraint(x) > 0, otherwise we stop the algorithm.
  if (TestInfeasibility==1)
    {
      printf("Initial point must be strictly admissible. ie constraint(x) > 0.\n");
      CONVERGENCE=0;
      goto Sortie;
    }

  AllConstraints_Value(vect_constraints, x, &all_constraints);
  f = PNL_EVAL_RNFUNCR(func,x);
  nbr_func_eval++;

  // We choose the initial value of mu in accordance with the objective function and the constraints
  // with some max and min values
  mu=0.;
  for (j=0; j<nbr_constraints; j++) mu += (log(GET(vect_constraints,j)));
  mu = MAX(fabs(f/mu), 1e-1);
  mu = MIN(mu, 1e1);

  // Initialization of lagrange variable lambda
  for (j=0; j<nbr_constraints; j++) LET(lambda,j) = mu / GET(vect_constraints,j);

  // Initialization of BFGS matrix using identity.
  // M should be positive definite.
  pnl_mat_set_id(M);

  // We keep in memory the best iteration
  // Since all teh iteration of the algorithm are feasible (even strictly), x_min is feasible
  f_min = f;
  pnl_vect_clone(x_min, x);
  pnl_vect_clone(lambda_min, lambda);

  //****************** Start of Global Algorithm (A) *******************//
  do
    {
      inner_iter=0;
      // We start the inner algorithm at the best iteration.
      pnl_vect_clone(x, x_min);
      pnl_vect_clone(lambda, lambda_min);
      AllConstraints_Value(vect_constraints, x, &all_constraints);
      f = f_min;

      //****************** Start of Inner Algorithm (A_mu) *******************//
      do
        {
          nbr_iterations++;
          inner_iter++;

          // Objective function gradient, Constraint gradient
          gradient_f(grad_func, func, x, grad_f);
          if (grad_func==NULL) nbr_func_eval+=2;
          gradient_c(&all_constraints, x, grad_c);

          for (i = 0; i < nbr_var; i++)
            {
              for (j = 0; j < nbr_var; j++)
                {
                  sum = 0;
                  for (k = 0; k < nbr_constraints; k++)
                    {
                      sum += MGET(grad_c,i,k) * GET(lambda, k) * MGET(grad_c,j,k) / GET(vect_constraints,k);
                    }
                  MLET(A,i,j) = MGET(M,i,j) + sum;
                  MLET(M_old,i,j) = MGET(M,i,j); // Save old BFGS matrix M
                }
            }

          for (i = 0; i < nbr_var; i++)
            {
              sum = 0;
              for (k = 0; k < nbr_constraints; k++)
                {
                  sum += MGET(grad_c,i,k)/GET(vect_constraints,k);
                }
              LET(d_x,i) = -GET(grad_f,i) + mu*sum;

            }

          //Compute primal descent direction d_x by solving linear system
          pnl_mat_ls (A, d_x);

          //Compute dual descent direction d_lambda
          for (i = 0; i < nbr_constraints; i++)
            {
              sum = 0;
              for (k = 0; k < nbr_var; k++)
                {
                  sum += MLET(grad_c,k,i) * GET(d_x,k);
                }
              LET(d_lambda,i) = (mu - GET(lambda,i)*sum)/GET(vect_constraints,i) - GET(lambda,i);
            }

          // Compute the value of merit function psi, used in the Line Search :
          psi = f;
          for (i = 0; i < nbr_constraints; i++)
            {
              psi += GET(lambda,i)*GET(vect_constraints,i) - mu*log(GET(lambda,i)*SQR(GET(vect_constraints,i)));
            }

          // Scalar product between gradient of merit function psi and primal descent direction d_x : <grad_psi, d_x>
          grad_psi_d = 0;
          for (i = 0; i < nbr_var; i++)
            {
              sum = 0;
              for (k = 0; k < nbr_constraints; k++)
                {
                  sum += MGET(grad_c,i,k)*(GET(lambda,k)-2*mu/GET(vect_constraints,k));
                }
              grad_psi_d += (GET(grad_f,i)+ sum) * GET(d_x,i);
            }

          for (i = 0; i < nbr_constraints; i++)
            {
              grad_psi_d += (GET(vect_constraints,i)-mu/GET(lambda,i))*GET(d_lambda,i);
            }

          // Save old values of x, lambda and gradients
          pnl_vect_clone(x_old, x);
          pnl_vect_clone(lambda_old, lambda);
          pnl_vect_clone(grad_f_old, grad_f);
          pnl_mat_clone(grad_c_old, grad_c);
          norm_x = pnl_vect_norm_one(x);
          norm_d_x = pnl_vect_norm_one(d_x);

          // Line Search : find a good step size alpha, verifying Armijo rule.
          alpha = 1. / (decrease_alpha); // Initial value for coefficient alpha (= 1 ?)
          do
            {
              alpha *= decrease_alpha; // decrease

              // To prevent d_x = 0 or alpha = 0, in this case we end the inner loop.
              norm_delta_x = alpha*norm_d_x;
              if (norm_delta_x/norm_x < norm_delta_x_min)
                {
                  CONVERGENCE=1;
                  pnl_vect_clone(x, x_min);
                  pnl_vect_clone(lambda, lambda_min);
                  goto End_Inner_Algo;
                }

              for (i = 0; i < nbr_var; i++) LET(x,i) = GET(x_old,i) + alpha*GET(d_x,i);
              for (i = 0; i < nbr_constraints; i++) LET(lambda,i) = GET(lambda_old,i) + alpha*GET(d_lambda,i);


              // Compute Psi(x+alpha*d_x, lambda+alpha*d_lambda)
              AllConstraints_Value(vect_constraints, x, &all_constraints);
              psi_alpha = 0;
              TestInfeasibility = 0; // If (x, lambda) is not feasible we go back to the decrease of alpha.
              for (i = 0; i < nbr_constraints; i++)
                {
                  if ((GET(lambda,i) < 0)|| (GET(vect_constraints,i) < 0))
                    {
                      TestInfeasibility = 1;
                      break;
                    }
                  psi_alpha += GET(lambda,i)*GET(vect_constraints,i) - mu * log(GET(lambda,i)*SQR(GET(vect_constraints,i)));
                }

              if (TestInfeasibility==0)
                {
                  // We compute f(x) only if it's strictly feasible.
                  f = PNL_EVAL_RNFUNCR(func,x);
                  nbr_func_eval++;
                  psi_alpha += f;
                  if (f<f_min)// We keep in memory the best iteration.
                    {
                      f_min = f;
                      pnl_vect_clone(x_min, x);
                      pnl_vect_clone(lambda_min, lambda);
                    }
                }
            }
          while ((TestInfeasibility==1)||((psi_alpha >= psi + omega * alpha * grad_psi_d)));

          /***************************************************************************
           * Update of matrix M by the BFGS formula with the use of Powell correction *
           ****************************************************************************/
          gradient_f(grad_func, func, x, grad_f);
          if (grad_func==NULL) nbr_func_eval+=2;
          gradient_c(&all_constraints, x, grad_c);

          norm_gamma = 0;
          delta_gamma = 0;
          // Compute delta_k and gamma_k (notation of the article):
          // delta_k = x_{k+1}-x_k
          // gamma_k = grad_lagrangian_x(x_{k+1}, lambda_k) - grad_lagrangian_x(l)(x_k, lambda_k)
          for (i = 0; i < nbr_var; i++)
            {
              LET(delta, i) = GET(x, i) - GET(x_old, i);
              LET(gamma, i) = GET(grad_f, i) - GET(grad_f_old, i);

              sum = 0; // sum = grad_c * lambda
              for (j = 0; j < nbr_constraints; j++)
                {
                  sum += (MGET(grad_c,i,j) - MGET(grad_c_old,i,j)) * GET(lambda, j);
                }

              LET(gamma, i) = GET(gamma, i) - sum;

              delta_gamma += GET(delta, i) * GET(gamma, i); // delta_gamma=<gamma,delta>
              norm_gamma += SQR(GET(gamma, i)); // norm_gamma=norm(gamma)^2
            }

          delta_M_delta=pnl_mat_scalar_prod(M_old,delta,delta); //transpose(delta)*M*delta
          // Up date of M only if BFGS formula is well defined.
          if ((fabs(delta_gamma) > eta) && (fabs(delta_M_delta) > eta))
            {
              // Powell's Correction, to ensure that the BFGS matrix is def.positive
              // We replace gamma by a vector to have <delta,gamma> > 0
              if (delta_gamma < kappa * delta_M_delta)
                {
                  theta = (1 - kappa) * delta_M_delta / (delta_M_delta - delta_gamma);
                  pnl_vect_axpby((1 - theta), M_delta, theta, gamma);
                  delta_gamma = kappa * delta_M_delta;
                }

              // Compute M*delta and transpose(delta)*M
              pnl_mat_mult_vect_inplace(M_delta,M_old,delta);
              pnl_mat_mult_vect_transpose_inplace(delta_M,M_old,delta);
              // Self-Scaling coefficient as proposed by Al-Baali
              SelfScalingBFGS = 1.;
              if (delta_gamma < delta_M_delta)
                {
                  SelfScalingBFGS = delta_gamma / delta_M_delta;
                }

              // M_k+1 = BFGS(SelfScalingBFGS * M_k, gamma_k, delta_k)
              for (i = 0; i < nbr_var; i++)
                {
                  for (j = 0; j < nbr_var; j++)
                    {
                      MLET(M,i,j) = SelfScalingBFGS * (MGET(M_old,i,j) - GET(M_delta,i) * GET(delta_M,j) / delta_M_delta) + GET(gamma,i) * GET(gamma,j)/delta_gamma;
                    }
                }

            }
          ///*** END UP-DATE OF BFGS MATRIX M

          // Errors used in stopping test of global algorithm
          pnl_mat_mult_vect_inplace(grad_lagrangien, grad_c, lambda);
          pnl_vect_axpby(1., grad_f, -1., grad_lagrangien);
          norm_grad_lagrangien = pnl_vect_norm_infty(grad_lagrangien);

          pnl_vect_clone(comp_cond_vect, lambda);
          pnl_vect_mult_vect_term(comp_cond_vect, vect_constraints);
          norm_comp_cond = pnl_vect_norm_infty(comp_cond_vect);
          OptimalityCriterion = ((norm_grad_lagrangien < tolerance) && (norm_comp_cond < tolerance));

          if (print_algo_steps>1)
            {
              printf("iter=%3i  f_eval=%3i   mu=%.3e   f=%.10e   ", nbr_iterations, nbr_func_eval, mu, f);
              printf("|grad_lagr|=%e   |comp_cond|=%e |Step|=%e\n", norm_grad_lagrangien, norm_comp_cond, norm_delta_x);
            }
        }
      while ((!OptimalityCriterion) && (inner_iter < inner_iter_max));

End_Inner_Algo:

      // Decrease of perturbation parameter
      if (inner_iter<MAX(3,inner_iter_max/2)) mu *= 0.1;
      else mu *= 0.5;
    }
  while ((!OptimalityCriterion) && (nbr_iterations < iter_max) && mu>mu_min);

Sortie:

  if (nbr_iterations >= iter_max) CONVERGENCE=2;
  if (OptimalityCriterion) CONVERGENCE=3;

  pnl_vect_clone(x_output, x_min);
  pnl_vect_clone(lambda, lambda_min);
  temps_final = clock();
  temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;

  if (print_algo_steps>0 && CONVERGENCE!=0)
    {
      printf("\n**********************************************************************************\n");
      if (CONVERGENCE==0) printf("************ Output Status: Failure: Initial point is not strictly feasible. \n");
      if (CONVERGENCE==1) printf("************ Output Status: Step too small, we stop the algorithm. \n");
      if (CONVERGENCE==2) printf("************ Output Status: Maximum iteration reached. \n");
      if (CONVERGENCE==3) printf("************ Output Status: A solution has been found up to the required accuracy. \n");
      printf("\n****** Initial point= ");
      pnl_vect_print(x_input);
      printf("****** Objective_function(Initial point)  = %e \n\n", PNL_EVAL_RNFUNCR(func, x_input));

      printf("******** Output primal-point= ");
      pnl_vect_print(x_output);
      printf("******** Objective_function(Output point)  = %.10e \n", f_min);
      printf("******** Output dual-point  = ");
      pnl_vect_print(lambda_min);
      printf("******** Norm_Grad_Lagrangien(Output point)= %.10e \n", norm_grad_lagrangien);
      printf("******** Norm_Comp_Condition(Output point) = %.10e \n", norm_comp_cond);
      printf("******** Nbr of iterations  = %i \n", nbr_iterations);
      printf("******** Nbr of evaluations = %i \n", nbr_func_eval);
      printf("******** Time elapsed       = %f seconds", temps_cpu);
      printf("\n**********************************************************************************\n");
    }

  AllConstraints_Free(&all_constraints);
  pnl_vect_free(&comp_cond_vect);
  pnl_vect_free(&grad_lagrangien);
  pnl_vect_free(&vect_constraints);
  pnl_vect_free(&lambda);
  pnl_vect_free(&lambda_old);
  pnl_vect_free(&lambda_min);
  pnl_vect_free(&x);
  pnl_vect_free(&x_old);
  pnl_vect_free(&x_min);
  pnl_vect_free(&d_x);
  pnl_vect_free(&d_lambda );
  pnl_vect_free(&delta);
  pnl_vect_free(&M_delta);
  pnl_vect_free(&delta_M);
  pnl_vect_free(&gamma);
  pnl_vect_free(&grad_f);
  pnl_vect_free(&grad_f_old);
  pnl_mat_free(&M);
  pnl_mat_free(&M_old);
  pnl_mat_free(&A);
  pnl_mat_free(&grad_c);
  pnl_mat_free(&grad_c_old);

  return CONVERGENCE;
}
