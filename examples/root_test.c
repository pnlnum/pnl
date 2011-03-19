
/*
 * Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as 
 * published by  the Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License  along with this program.  If not, see
 * <http://www.gnu.org/licenses/>. 
 */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pnl/pnl_root.h"
#include "pnl/pnl_matrix.h"

#include "tests_utils.h"

static double epsabs = 0.0001;
static double epsrel = 0.00001;
static int N_max = 1000;

#define FUNC  f_cos
#define FDF_FUNC fdf_cos 

static double x_square(double x, void *p) { return x*x -2.; }

static void x_fdf_square(double x, double *f, double *df, void *p)
{
  *f = x*x -2.;
  *df = 2*x;
}

static double f_cos (double x, void *p) {return cos (x);}
static void fdf_cos (double x, double *f, double *df, void *p) {*f = cos (x); *df = -sin(x);}

static void bisection_test ()
{
  double x1, x2, r;
  PnlFunc func;
  int status;
  
  x1 = 0;
  x2 = 3;
  func.function = FUNC;
  func.params = NULL;

  status = pnl_root_bisection(&func, x1, x2, epsrel, epsabs, N_max, &r);
  printf("bisection : root is %f, status==OK %d\n", r, status==OK);

}

static void newton_test ()
{
  double x0, r;
  PnlFuncDFunc func;
  int status;
  
  x0 = 0.5;
  func.function = FDF_FUNC;
  func.params = NULL;

  status = pnl_root_newton(&func, x0, epsrel, epsabs, N_max, &r);
  printf("newton : root is %f, status==OK %d\n", r, status==OK);
}

static void brent_test()
{
  double x1, x2, tol, r;
  PnlFunc func;

  x1 = -1;
  x2 = 3;
  tol = 0.001;
  func.function = FUNC;
  func.params = NULL;

  r = pnl_root_brent(&func, x1, x2, &tol);
  printf("brent : root is %f (tol = %f) \n", r, tol);
}

static void find_root_test ()
{
  double x1, x2, tol, r;
  PnlFuncDFunc func;
  int status;
  
  x1 = -1;
  x2 = 3;
  tol = 0.001;
  func.function = FDF_FUNC;
  func.params = NULL;

  status = pnl_find_root(&func, x1, x2, tol, N_max, &r);
  printf("find_root : root is %f, status==OK %d\n", r, status==OK);
}

/*
 * Test of minpack routines
 *
 * These examples were include in the C/C++ Minpack package and
 * have translated into the Pnl format
 */

static void Dfcn_fsolve(const PnlVect *x, PnlMat *fjac, void *p)
{
  int k, j, n;
  double one=1, four=4, three=3, two=2, zero=0;

  n = x->size;
  pnl_mat_resize (fjac, n, n);
  for (k=0; k<n; k++)
    {
      for (j=0; j < n; j++)
        {
          MLET(fjac, k, j) = zero;
        }
      MLET(fjac, k, k) = three - four * GET(x,k);
      if (k != 0) MLET (fjac, k-1, k) = -one;
      if (k != n-1) MLET(fjac, k+1, k) = -two;
    }      
}

static void fcn_fsolve(const PnlVect *x, PnlVect *fvec, void *p)
{
  int k, n;
  double one=1, temp, temp1, temp2, three=3, two=2, zero=0;

  n = x->size;
  pnl_vect_resize (fvec, n);
  for (k=0; k<n; k++)
    {
      
      temp = (three - two*GET(x,k))*GET(x,k);
      temp1 = zero;
      if (k != 0) temp1 = GET(x,k-1);
      temp2 = zero;
      if (k != n-1) temp2 = GET(x,k+1);
      LET(fvec,k) = temp - temp1 - two*temp2 + one;
    }
}

static void fcn_lsq(const PnlVect *x, PnlVect *fvec, void *p)
{
  int i;
  double tmp1, tmp2, tmp3;
  double y[15] = {1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1, 3.5e-1,
		  3.9e-1, 3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1, 1.34, 2.1, 4.39};

  for (i = 0; i < 15; i++)
    {
      tmp1 = i+1;
      tmp2 = 15 - i;
      tmp3 = tmp1;
      if (i > 7) tmp3 = tmp2;
      LET(fvec,i) = y[i] - (GET(x,0) + tmp1/(GET(x,1)*tmp2 + GET(x,2)*tmp3));
    }
}

static void Dfcn_lsq(const PnlVect *x, PnlMat *fjac, void *p)
{
  int i;
  double tmp1, tmp2, tmp3, tmp4;
  for (i = 0; i < 15; i++)
    {
      tmp1 = i+1;
      tmp2 = 15 - i;
      tmp3 = tmp1;
      if (i > 7) tmp3 = tmp2;
      tmp4 = (GET(x,1)*tmp2 + GET(x,2)*tmp3); tmp4 = tmp4*tmp4;
      MLET(fjac, i, 0) = -1.;
      MLET(fjac, i, 1) = tmp1*tmp2/tmp4;
      MLET(fjac, i, 2) = tmp1*tmp3/tmp4;
    }
}

static void test_hybrX ()
{
  int j, n, maxfev, info, nfev;
  double xtol, fnorm;
  PnlVect *x, *fvec, *diag;
  PnlRnFuncRnDFunc f;

  n = 9;
  x = pnl_vect_create (n);
  fvec = pnl_vect_create (n);
  diag = pnl_vect_create (n);

  /* the following starting values provide a rough solution. */
  pnl_vect_set_double (x, -1);

  /* default value for xtol */
  xtol = 0;

  maxfev = 2000;
  pnl_vect_set_double (diag, 1);


  /*
   * Test without Jacobian
   */
  printf ("Test of pnl_root_fsolve without user supplied Jacobian.\n\n");
  f.function = fcn_fsolve;
  f.Dfunction = NULL;
  f.params = NULL;
  info = pnl_root_fsolve (&f, x, fvec, xtol, maxfev, &nfev, diag, FALSE);
  fnorm = pnl_vect_norm_two(fvec);
  printf("     final l2 norm of the residuals %15.7g\n\n", fnorm);
  printf("     number of function evaluations  %10i\n\n", nfev);
  printf("     exit parameter                  %10i\n\n", info);
  printf("     final approximate solution\n");
  for (j=1; j<=n; j++) printf("%s%15.7g", j%3==1?"\n     ":"", GET(x,j-1));
  printf("\n\n");
  
  /*
   * Test with Jacobian
   */
  printf ("Test of pnl_root_fsolve without user supplied Jacobian.\n\n");
  f.function = fcn_fsolve;
  f.Dfunction = Dfcn_fsolve;
  f.params = NULL;
  info = pnl_root_fsolve (&f, x, fvec, xtol, maxfev, &nfev, diag, FALSE);
  fnorm = pnl_vect_norm_two(fvec);
  printf("     final l2 norm of the residuals %15.7g\n\n", fnorm);
  printf("     number of function evaluations  %10i\n\n", nfev);
  printf("     exit parameter                  %10i\n\n", info);
  printf("     final approximate solution\n");
  for (j=1; j<=n; j++) printf("%s%15.7g", j%3==1?"\n     ":"", GET(x,j-1));
  printf("\n\n");

  pnl_vect_free (&x);
  pnl_vect_free (&fvec);
  pnl_vect_free (&diag);
}


static void test_lmdif ()
{
  int m, n, info, nfev, maxfev;
  double tol, fnorm;
  PnlVect *x, *fvec;
  PnlRnFuncRmDFunc f;

  m = 15;
  n = 3;

  x = pnl_vect_create (n);
  fvec = pnl_vect_create (m);

  /* the following starting values provide a rough fit. */

  pnl_vect_set_double (x, 1.);
  /* default vlaues */
  tol = 0;
  maxfev = 0;

  /*
   * Test without user supplied Jacobian
   */
  printf ("Test of pnl_root_fsolve_lsq without user supplied Jacobian.\n\n");
  f.function = fcn_lsq;
  f.Dfunction = NULL;
  f.params = NULL;
  info = pnl_root_fsolve_lsq(&f, x, m, fvec, tol, tol, 0., maxfev, &nfev, NULL, TRUE);

  fnorm = pnl_vect_norm_two(fvec);

  printf("      final l2 norm of the residuals%15.7f\n\n",fnorm);
  printf("      exit parameter                %10i\n\n", info);
  printf("      final approximate solution\n\n %15.7f%15.7f%15.7f\n\n",
	 GET(x,0), GET(x,1), GET(x,2));

  /*
   * Test with user supplied Jacobian
   */
  printf ("Test of pnl_root_fsolve_lsq with user supplied Jacobian.\n\n");
  f.function = fcn_lsq;
  f.Dfunction = Dfcn_lsq;
  f.params = NULL;
  info = pnl_root_fsolve_lsq(&f, x, m, fvec, tol, tol, 0., maxfev, &nfev, NULL, TRUE);

  fnorm = pnl_vect_norm_two(fvec);

  printf("      final l2 norm of the residuals%15.7f\n\n",fnorm);
  printf("      exit parameter                %10i\n\n", info);
  printf("      final approximate solution\n\n %15.7f%15.7f%15.7f\n\n",
	 GET(x,0), GET(x,1), GET(x,2));


  pnl_vect_free (&x);
  pnl_vect_free (&fvec);
}


int main ()
{
  brent_test () ;
  bisection_test ();
  newton_test ();
  find_root_test ();
  test_hybrX ();
  test_lmdif ();
  return OK;
}
