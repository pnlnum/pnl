/************************************************************************/
/* Copyright Ismail Laachir <ismail.laachir@gmail.fr> 2009.             */
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_random.h"
#include "pnl/pnl_optim.h"
#include "tests_utils.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 *  This file contains routines to implement test problem HS15 from
 *  the Hock & Schittkowski collection.
 *
 *  min   100 (x[1] - x[0]^2)^2 + (1 - x[0])^2
 *  s.t.  x[0] x[1] >= 1
 *        x[0] + x[1]^2 >= 0
 *        x[0] <= 0.5
 *
 *  The standard minimum is (0.5, 2.0), with final objective = 306.5.
 *  Sometimes the solver converges to another local minimum
 *  at (-0.79212, -1.26243), with final objective = 360.4.
 *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// Constraints on variables
// Here we implement the two constraints "x1*x2>=1" and "x1+x2^2>=0"
static void NL_Constraints_HS15(const PnlVect *x, PnlVect *res, void *UpperLower)
{
  pnl_vect_resize(res, 2);

  LET(res, 0) = GET(x,0) * GET(x,1)-1.;
  LET(res, 1) = GET(x,0) + SQR(GET(x,1));
}

// Objective function. f(x1,x2)=100*(x2-x1^2)^2 + (1-x1)^2
static double func_HS15(const PnlVect *x, void *pt)
{
  return 100.0 * SQR(GET(x,1) - SQR(GET(x,0))) + SQR(1.0 - GET(x,0));
}

// Gradient of objective function
// grad[i] = derivative of f with respect to x[i]
static void grad_func_HS15(const PnlVect *x, PnlVect *grad, void *pt)
{
  double dTmp = GET(x,1) - SQR(GET(x,0));
  pnl_vect_resize(grad, 2);

  LET(grad, 0) = (-400.0 * dTmp * GET(x,0)) - (2.0 * (1.0 - GET(x,0)));
  LET(grad, 1) = 200.0 * dTmp;
}

// In order to test the algorithm, the starting point is chosen randomly.
// Of course, in practical cases, one chooses a starting point in accordance with the problem.
static void choose_random_input_point(int type_generator, PnlVect* x_input)
{
  double theta, x_min, x_max;
  theta = pnl_rand_uni(type_generator);
  x_min = 0;
  x_max = 0.5;
  LET(x_input, 0) = theta*x_min + (1-theta)*x_max;

  theta = pnl_rand_uni(type_generator);
  x_min = 1./GET(x_input, 0);
  x_max = x_min+10.;
  LET(x_input, 1) = theta*x_min + (1-theta)*x_max;
}

static void minimize_func_HS15()
{
  double tolerance;
  int CONVERGENCE, type_generator, print_algo_steps, iter_max;

  PnlVect *x_input, *x_output;
  PnlVect *lower_bounds, *upper_bounds;

  PnlRnFuncR FuncToMinimize; // Structure that contains the objective function
  PnlRnFuncRm GradFuncToMinimize; // Structure that contains the gradient objective function
  PnlRnFuncRm NL_Constraints; // Structure that contains the constraints

  type_generator = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init(type_generator, 2, 1);

  tolerance = 1e-9; // Optimality termination tolerance
  iter_max = 500; // Maximum number of iteration

  if ( pnl_test_is_verbose () )
    print_algo_steps = 2; 
  else
    print_algo_steps = 0; 

  /*print_algo_steps=0 : we print no information
    print_algo_steps=1 : we print final results
    print_algo_steps>1 : we print information at each iteration*/

  x_input = pnl_vect_create(2); //starting point
  x_output = pnl_vect_create(2); //Output point

  // Lower and upper bound of the variable x
  lower_bounds = pnl_vect_create_from_list(2, PNL_NEGINF, PNL_NEGINF);
  upper_bounds = pnl_vect_create_from_list(2, 0.5, PNL_POSINF);

  // Objective function
  FuncToMinimize.F = &func_HS15;
  FuncToMinimize.params = NULL;

  // Gradient of objective function
  GradFuncToMinimize.F = &grad_func_HS15;
  GradFuncToMinimize.params = NULL;

  // Constraints
  NL_Constraints.F = &NL_Constraints_HS15;
  NL_Constraints.params = NULL;

  // We use a random starting point
  // Of course, in practical cases, one chooses a starting point in accordance with the problem.
  // !! Initial point must be strictly admissible. ie constraint(x) > 0. !!
  choose_random_input_point(type_generator, x_input);

  // Call of the optimization routine.
  CONVERGENCE = pnl_optim_intpoints_bfgs_solve(&FuncToMinimize, &GradFuncToMinimize, &NL_Constraints, lower_bounds, upper_bounds, x_input, tolerance, iter_max, print_algo_steps, x_output);


  if (CONVERGENCE!=0)
    {
    /* "\n\nThe standard minimum is (%f, %f), Sometimes the solver
     * converges to another local minimum at (-0.79212, -1.26243)\n", 0.5,
     * 2.0 */
      PnlVect *res = pnl_vect_create_from_list (2, 0.5, 2.0);
      pnl_test_vect_eq_abs (x_output, res, 1E-8, "intpoints_bfgs", "An other solution may be (-0.79212, -1.26243)");
      pnl_vect_free (&res);
    }

  pnl_vect_free(&x_input);
  pnl_vect_free(&x_output);
  pnl_vect_free(&lower_bounds);
  pnl_vect_free(&upper_bounds);
}


int main (int argc, char **argv)
{
  pnl_test_init (argc, argv);
  minimize_func_HS15();
  exit (pnl_test_finalize ("Optim"));
}
