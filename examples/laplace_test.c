
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

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_laplace.h"
#include "tests_utils.h"

const double tol = 1E-4;

static double exp_law_density (double x, void *p)
{
  double mu = *((double *)p);
  if (x < 0) return 0.;
  return mu * exp ( -mu * x);
}

static double exp_law_cdf (double x, void *p)
{
  double mu = *((double *)p);
  if (x < 0) return 0.;
  return 1 - exp ( -mu * x);
}

static dcomplex exp_law_laplace (dcomplex l, void *p)
{
  double mu = *((double *)p);
  return RCdiv (mu, RCadd (mu, l));
}

static double exp_law_real_laplace (double l, void *p)
{
  double mu = *((double *)p);
  return mu / (mu + l);
}



static void euler_test()
{
  PnlCmplxFunc lap;
  PnlFunc density, cdf;
  PnlVect *res;
  double t;
  double mu;
  int M, N;

  lap.F = exp_law_laplace;
  lap.params = &mu;
  density.F = exp_law_density;
  density.params = &mu;
  cdf.F = exp_law_cdf;
  cdf.params = &mu;
  res = pnl_vect_create (0);

  M = N = 15;
  mu = 2;
  t = .5;
  pnl_test_eq(pnl_ilap_euler (&lap, t, N, M), PNL_EVAL_FUNC (&density, t),
              tol, "pnl_ilap_euler", "");
  pnl_ilap_fft (res, &lap, t, 1E-6);
  pnl_test_eq (pnl_vect_get (res, res->size - 1), PNL_EVAL_FUNC (&density, t),
               tol, "pnl_ilap_ftt", "");
  pnl_test_eq (pnl_ilap_cdf_euler (&lap, t, 0.05, 10000, 1), PNL_EVAL_FUNC (&cdf, t),
               tol, "pnl_ilap_cdf_euler", "");
  pnl_vect_free (&res);
}


static void gs_test ()
{
  PnlFunc lap;
  PnlFunc density;
  double t;
  double mu;
  int n;

  lap.F = exp_law_real_laplace;
  lap.params = &mu;
  density.F = exp_law_density;
  density.params = &mu;

  n = 10;
  mu = 2;
  t = .5;
  pnl_test_eq (pnl_ilap_gs_basic (&lap, t, n), PNL_EVAL_FUNC (&density, t),
               1E-1, "pnl_ilap_gs_basic", "");
  pnl_test_eq (pnl_ilap_gs (&lap, t, n), PNL_EVAL_FUNC (&density, t),
               tol, "pnl_ilap_gs", "");

}


int main (int argc, char **argv)
{
  pnl_test_init (argc, argv);
  euler_test ();
  gs_test ();
  exit(pnl_test_finalize ("Laplace transform"));
}
