
/************************************************************************/
/* Copyright Jérôme Lelong <jerome.lelong@gmail.com>                    */
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
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_basis.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_random.h"
#include "tests_utils.h"

#define function exp

static int PRINT_COEFF = 0;

/*
 * example of how to use  pnl_basis_fit_ls to regress on a basis.
 * regression of the exponential function on the grid [0:0.05:5]
 */
static void exp_regression2()
{
  int i, n, basis_name, basis_dim, space_dim;
  double a, b, h, err;
  PnlMat *t;
  PnlVect *y, *alpha;
  PnlBasis *basis;

  /* creating the grid */
  a = 0.0;
  b = 5.0;
  n = 100;
  h = (b - a) / n;
  t = pnl_mat_create_from_scalar(n + 1, 1, h);
  pnl_mat_set(t, 0, 0, 0.0);
  pnl_mat_cumsum(t, 'r');

  /* creating the values of exp on the grid */
  y = pnl_vect_create(n + 1);
  for (i = 0 ; i < n + 1 ; i++)
    {
      pnl_vect_set(y, i, function(pnl_mat_get(t, i, 0)));
    }


  basis_name = PNL_BASIS_HERMITIAN; /* PNL_BASIS_TCHEBYCHEV; */
  space_dim = 1; /* real valued basis */
  basis_dim = 5; /* number of elements in the basis */
  basis = pnl_basis_create(basis_name, basis_dim, space_dim);
  alpha = pnl_vect_new();
  pnl_basis_fit_ls(basis, alpha, t, y);
  if (PRINT_COEFF)
    {
      printf("\tcoefficients of the decomposition : ");
      pnl_vect_print_asrow(alpha);
    }
  /* computing the infinity norm of the error */
  err = 0.;
  for (i = 0; i < t->m; i++)
    {
      double tmp = function(pnl_mat_get(t, i, 0)) -
                   pnl_basis_eval(basis, alpha, pnl_mat_lget(t, i, 0));
      if (fabs(tmp) > err) err = fabs(tmp);
    }

  pnl_test_eq_abs(err, 1.812972, 1E-5, "pnl_basis_eval", "exponential function on [0:0.05:5]");

  pnl_basis_free(&basis);
  pnl_mat_free(&t);
  pnl_vect_free(&y);
  pnl_vect_free(&alpha);
}

static double function2d_1(double *x)
{
  return (x[0]*x[0] + x[1]*x[1] + x[0] * x[1]) * (x[0] - x[1]);
}

static double function2d_2(double *x)
{
  return log(1 + 2 * (x[0] * x[0] + x[1] * x[1])) + (x[0] * x[0] + x[1] * x[1]);
}

static double myfunction(const PnlVect *x, void *p)
{
  double a = ((double *) p)[0];
  return log(1 + a * SQR(pnl_vect_norm_two(x)));
}

static double regression_multid_aux(double(*f)(double *), PnlBasis *basis)
{
  int i, j, n;
  double a, b, h, err;
  PnlMat *t;
  PnlVect *y, *alpha;

  /* creating the grid */
  a = -2.0;
  b = 2.0;
  n = 100;
  h = (b - a) / n;
  t = pnl_mat_create((n + 1) * (n + 1), 2);

  /* creating the values of exp on the grid */
  y = pnl_vect_create((n + 1) * (n + 1));
  for (i = 0; i < n + 1; i++)
    {
      for (j = 0; j < n + 1; j++)
        {
          MLET(t, i * (n + 1) + j, 0) =  a + i * h;
          MLET(t, i * (n + 1) + j, 1) =  a + j * h;
          LET(y, i * (n + 1) + j) = (*f)(pnl_mat_lget(t, i * (n + 1) + j, 0));
        }
    }

  alpha = pnl_vect_new();
  pnl_basis_fit_ls(basis, alpha, t, y);
  if (PRINT_COEFF)
    {
      printf("\tcoefficients of the decomposition : ");
      pnl_vect_print_asrow(alpha);
    }

  /* Compute the infinity norm of the error */
  err = 0.;
  for (i = 0; i < t->m; i++)
    {
      double tmp = (*f)(pnl_mat_lget(t, i, 0)) -
                   pnl_basis_eval(basis, alpha, pnl_mat_lget(t, i, 0));
      if (fabs(tmp) > err) err = fabs(tmp);
    }
  return err;
}

static void regression_multid()
{
    {
      double err;
      int basis_name = PNL_BASIS_HERMITE;
      int nb_variates = 2; /* functions with values in R^2 */
      int degree = 3; /* total sum degree */
      PnlBasis *basis = pnl_basis_create_from_degree(basis_name, degree, nb_variates);
      err = regression_multid_aux(function2d_1, basis);
      pnl_test_eq_abs(err, 0., 1E-5, "pnl_basis_eval (sum degree)", "polynomial with deg <= 3");
      pnl_basis_free(&basis);
    }

    {
      double err;
      int basis_name = PNL_BASIS_HERMITE;
      int nb_variates = 2; /* functions with values in R^2 */
      int degree = 3; /* total product degree */
      PnlBasis *basis = pnl_basis_create_from_prod_degree(basis_name, degree, nb_variates);
      err = regression_multid_aux(function2d_1, basis);
      pnl_test_eq_abs(err, 0., 1E-5, "pnl_basis_eval (prod degree)", "polynomial with deg <= 3");
      pnl_basis_free(&basis);
    }

    {
      double err;
      PnlRnFuncR F;
      double a = 2;
      int basis_name = PNL_BASIS_HERMITE;
      int nb_variates = 2; /* functions with values in R^2 */
      int degree = 3; /* total product degree */
      PnlBasis *basis = pnl_basis_create_from_prod_degree(basis_name, degree, nb_variates);
      F.F = myfunction;
      F.params = (void *) &a;
      pnl_basis_add_function(basis, &F);
      err = regression_multid_aux(function2d_2, basis);
      pnl_test_eq_abs(err, 0., 1E-5, "pnl_basis_eval (extra functions)", "log (1+x[0]*x[0] + x[1]*x[1]) + (x[0]*x[0] + x[1]*x[1]) on [-2,2]^2");
      pnl_basis_free(&basis);
    }
}

static double fonction_a_retrouver(double t, double x)
{
  return (1 - t) * x * x / 5;
}

static double derive_x_fonction_a_retrouver(double t, double x)
{
  return (1 - t) * 2 * x / 5;
}

static double derive_xt_fonction_a_retrouver(double t, double x)
{
  return -2 * x / 5;
}

static double derive_xx_fonction_a_retrouver(double t, double x)
{
  return (1 - t) * 2 / 5;
}

static double derive_t_fonction_a_retrouver(double t, double x)
{
  return -x * x / 5;
}

static double derive_x_approx_fonction(PnlBasis *B, PnlVect *alpha, double t, double x)
{
  double arg[2];
  arg[0] = t;
  arg[1] = x;
  return pnl_basis_eval_D(B, alpha, arg, 1);
}

static void derive_approx_fonction2(PnlBasis *B, PnlVect *D, PnlVect *alpha, double t, double x)
{
  double arg[2];
  arg[0] = t;
  arg[1] = x;

  LET(D, 3) = pnl_basis_eval_D(B, alpha, arg, 0); // GET(grad,0);
  LET(D, 1) = pnl_basis_eval_D(B, alpha, arg, 1); //GET(grad,1);
  LET(D, 0) = pnl_basis_eval(B, alpha, arg);
  LET(D, 2) = pnl_basis_eval_D2(B, alpha, arg, 1, 1); // PNL_MGET(Hes, 1, 1);
  LET(D, 4) = pnl_basis_eval_D2(B, alpha, arg, 0, 1); // PNL_MGET(Hes, 0, 1);
}

static void derive_approx_fonction(PnlBasis *B, PnlVect *D, PnlVect *alpha, double t, double x)
{
  double sum0;
  double arg[2];
  PnlMat *Hes = pnl_mat_new();
  PnlVect *grad = pnl_vect_new();
  arg[0] = t;
  arg[1] = x;
  sum0 = 0.0; //calcule la valeur de la fonction

  pnl_basis_eval_derivs(B, alpha, arg, &sum0, grad, Hes);
  LET(D, 3) = GET(grad, 0);
  LET(D, 1) = GET(grad, 1);
  LET(D, 0) = sum0;
  LET(D, 2) = PNL_MGET(Hes, 1, 1);
  LET(D, 4) = PNL_MGET(Hes, 0, 1);

  pnl_mat_free(&Hes);
  pnl_vect_free(&grad);

}

static void pnl_basis_eval_test()
{
  PnlMat    *X;
  PnlVect   *V, *x, *t, *D, *alpha, *lower, *upper;
  PnlRng    *rng;
  PnlBasis  *basis;
  int        j, deg, n;
  double     t0, x0, tol;

  tol = 1E-5;
  deg = 5; //total degree
  n = 50;
  D = pnl_vect_create(5);
  x = pnl_vect_create(n);
  t = pnl_vect_create(n);
  t0 = 0.5;
  x0 = 2.5;
  rng = pnl_rng_create(PNL_RNG_MERSENNE);
  pnl_rng_sseed(rng, 0);

  /*
   * Random points where the function will be evaluted
   */
  pnl_vect_rng_uni(x, n, -5, 4, rng);
  pnl_vect_rng_uni(t, n, 0, 1, rng);
  basis = pnl_basis_create_from_degree(PNL_BASIS_HERMITIAN, deg, 2);
  alpha = pnl_vect_create(basis->nb_func);
  X = pnl_mat_create(n, 2);
  for (j = 0; j < n; j++)
    {
      MLET(X, j, 0) = GET(t, j);
      MLET(X, j, 1) = GET(x, j);
    }
  V = pnl_vect_create(n);
  /*
   * Vector of values for the function to recover
   */
  for (j = 0; j < n; j++)
    {
      LET(V, j) = fonction_a_retrouver(GET(t, j), GET(x, j));
    }
  pnl_basis_fit_ls(basis, alpha, X, V);


  /*
   * Test pnl_basis_eval_derivs
   */
  derive_approx_fonction(basis, D, alpha, t0, x0);
  pnl_test_eq_abs(pnl_vect_get(D, 0), fonction_a_retrouver(t0, x0), tol,
                  "pnl_basis_eval_derivs", "derivative 0");
  pnl_test_eq_abs(derive_x_approx_fonction(basis, alpha, t0, x0),
                  derive_x_fonction_a_retrouver(t0, x0), tol,
                  "pnl_basis_eval_derivs", "derivative %% x");
  pnl_test_eq_abs(pnl_vect_get(D, 2), derive_xx_fonction_a_retrouver(t0, x0),
                  tol, "pnl_basis_eval_derivs", "derivative %% xx");
  pnl_test_eq_abs(pnl_vect_get(D, 3), derive_t_fonction_a_retrouver(t0, x0),
                  tol, "pnl_basis_eval_derivs", "derivative %% t");
  pnl_test_eq_abs(pnl_vect_get(D, 4), derive_xt_fonction_a_retrouver(t0, x0),
                  tol, "pnl_basis_eval_derivs", "derivative %% tx");

  pnl_basis_free(&basis);

  /* reduced basis */
  basis = pnl_basis_create_from_degree(PNL_BASIS_HERMITIAN, deg, 2);
  lower = pnl_vect_create_from_list(2, 0., -5.);
  upper = pnl_vect_create_from_list(2, 1., 4.);
  pnl_basis_set_domain(basis, lower, upper);
  pnl_basis_fit_ls(basis, alpha, X, V);

  derive_approx_fonction(basis, D, alpha, t0, x0);
  pnl_test_eq_abs(pnl_vect_get(D, 0), fonction_a_retrouver(t0, x0), tol,
                  "pnl_basis_eval_derivs (reduced)", "derivative 0");
  pnl_test_eq_abs(derive_x_approx_fonction(basis, alpha, t0, x0),
                  derive_x_fonction_a_retrouver(t0, x0), tol,
                  "pnl_basis_eval_derivs (reduced)", "derivative %% x");
  pnl_test_eq_abs(pnl_vect_get(D, 2), derive_xx_fonction_a_retrouver(t0, x0),
                  tol, "pnl_basis_eval_derivs (reduced)",  "derivative %% xx");
  pnl_test_eq_abs(pnl_vect_get(D, 3), derive_t_fonction_a_retrouver(t0, x0),
                  tol, "pnl_basis_eval_derivs (reduced)", "derivative %% t");
  pnl_test_eq_abs(pnl_vect_get(D, 4), derive_xt_fonction_a_retrouver(t0, x0),
                  tol, "pnl_basis_eval_derivs (reduced)",  "derivative %% tx");

  /*
   * Test pnl_basis_eval_D and pnl_basis_eval_D2
   */
  derive_approx_fonction2(basis, D, alpha, t0, x0);
  pnl_test_eq_abs(pnl_vect_get(D, 0), fonction_a_retrouver(t0, x0), tol,
                  "pnl_basis_eval", "derivative 0");
  pnl_test_eq_abs(derive_x_approx_fonction(basis, alpha, t0, x0),
                  derive_x_fonction_a_retrouver(t0, x0), tol,
                  "pnl_basis_eval_D", "derivative %% x");
  pnl_test_eq_abs(pnl_vect_get(D, 2), derive_xx_fonction_a_retrouver(t0, x0),
                  tol, "pnl_basis_eval_D2", "derivative %% xx");
  pnl_test_eq_abs(pnl_vect_get(D, 3), derive_t_fonction_a_retrouver(t0, x0),
                  tol, "pnl_basis_eval_D", "derivative %% t");
  pnl_test_eq_abs(pnl_vect_get(D, 4), derive_xt_fonction_a_retrouver(t0, x0),
                  tol, "pnl_basis_eval_D2", "derivative %% tx");

  pnl_vect_free(&lower);
  pnl_vect_free(&upper);
  pnl_basis_free(&basis);

  /* reduced basis */
  basis = pnl_basis_create_from_degree(PNL_BASIS_HERMITIAN, deg, 2);
  lower = pnl_vect_create_from_list(2, 0., -5.);
  upper = pnl_vect_create_from_list(2, 1., 4.);
  pnl_basis_set_domain(basis, lower, upper);
  pnl_basis_fit_ls(basis, alpha, X, V);

  derive_approx_fonction2(basis, D, alpha, t0, x0);
  pnl_test_eq_abs(pnl_vect_get(D, 0), fonction_a_retrouver(t0, x0), tol,
                  "pnl_basis_eval (reduced)", "derivative 0");
  pnl_test_eq_abs(derive_x_approx_fonction(basis, alpha, t0, x0),
                  derive_x_fonction_a_retrouver(t0, x0), tol,
                  "pnl_basis_eval_D (reduced)", "derivative %% x");
  pnl_test_eq_abs(pnl_vect_get(D, 2), derive_xx_fonction_a_retrouver(t0, x0),
                  tol, "pnl_basis_eval_D2 (reduced)",  "derivative %% xx");
  pnl_test_eq_abs(pnl_vect_get(D, 3), derive_t_fonction_a_retrouver(t0, x0),
                  tol, "pnl_basis_eval_D (reduced)", "derivative %% t");
  pnl_test_eq_abs(pnl_vect_get(D, 4), derive_xt_fonction_a_retrouver(t0, x0),
                  tol, "pnl_basis_eval_D2 (reduced)",  "derivative %% tx");


  pnl_basis_free(&basis);
  pnl_rng_free(&rng);
  pnl_vect_free(&alpha);
  pnl_vect_free(&x);
  pnl_vect_free(&t);
  pnl_vect_free(&V);
  pnl_vect_free(&D);
  pnl_vect_free(&lower);
  pnl_vect_free(&upper);
  pnl_mat_free(&X);
}

int main(int argc, char **argv)
{
  pnl_test_init(argc, argv);
  if (pnl_test_is_verbose()) PRINT_COEFF = 1;
  exp_regression2();
  regression_multid();
  pnl_basis_eval_test();
  pnl_test_finalize("Basis functions");
  return OK;
}
