#include "pnl/pnl_mathtools.h"


extern int pnl_root_fsolve (PnlRnFuncRnDFunc *f, PnlVect *x, PnlVect *fx,
                       double xtol, int maxfev, int *nfev, PnlVect *scale,
                       int error_msg);
extern int pnl_root_fsolve_lsq (PnlRnFuncRmDFunc *f, PnlVect *x, int m, PnlVect
                           *fx,  double xtol, double ftol, double gtol, int
                           maxfev, int *nfev, PnlVect *scale, int
                           error_msg);

static void Dfcn_fsolve(const PnlVect *x, PnlMat *fjac, void *p)
{
  /*      subroutine fcn for hybrd example. */

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
  /*      subroutine fcn for hybrd example. */

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



int main()
{
  test_hybrX ();
  test_lmdif ();
  return 0;
}


