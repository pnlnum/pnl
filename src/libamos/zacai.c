#include "amos.h"


/* Copyright, Donald E. Amos: sandia national laboratories
 *            from slatec library or amos library.
 *
 * issued by sandia laboratories, a prime contractor to the 
 * united states department of energy 
 * notice:
 * this report was prepared as an account of work sponsored by the 
 * united states government.  neither the united states nor the 
 * united states department of energy, nor any of their 
 * employees, nor any of their contractors, subcontractors, or their 
 * employees, makes any warranty, express or implied, or assumes any 
 * legal liability or responsibility for the accuracy, completeness 
 * or usefulness of any information, apparatus, product or process 
 * disclosed, or represents that its use would not infringe 
 * privately owned rights. 
 *
 *
 * this code has been approved for unlimited release. 
 */


/*
 *    zacai applies the analytic continuation formula 
 * 
 *        k(fnu,zn*exp(mp))=k(fnu,zn)*exp(-mp*fnu) - mp*i(fnu,zn) 
 *                mp=pi*mr*cmplx(0.0,1.0) 
 * 
 *    to continue the k function from the right half to the left 
 *    half z plane for use with zairy where fnu=1/3 or 2/3 and n=1. 
 *    zacai is the same as zacon with the parts for larger orders and 
 *    recurrence removed. a recursive call to zacon can result if zacon 
 *    is called from zairy. 
 * 
 */

static const int c__1 = 1;

int amos_zacai (double *zr, double *zi, double *fnu,const int *kode, int *mr,const int *n,
		double *yr, double *yi, int *nz, double *rl, double *tol,
		double *elim, double *alim)
{
  static const double pi = 3.14159265358979324;
  /* Local variables */
  double dfnu;
  double ascle;
  double csgni, csgnr, cspni, cspnr;
  double az;
  int nn, nw;
  double yy, c1i, c2i, c1r, c2r, arg;
  int iuf;
  double cyi[2], fmr, sgn;
  int inu;
  double cyr[2], zni, znr;

  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  *nz = 0;
  znr = -(*zr);
  zni = -(*zi);
  az = amos_azabs (zr, zi);
  nn = *n;
  dfnu = *fnu + (double) (*n - 1);
  if (az <= 2.)
    {
      goto L10;
    }
  if (az * az * .25 > dfnu + 1.)
    {
      goto L20;
    }
 L10:
  /*----------------------------------------------------------------------- 
   *    POWER SERIES FOR THE I FUNCTION 
   *----------------------------------------------------------------------- 
   */
  amos_zseri (&znr, &zni, fnu, kode, &nn, &yr[1], &yi[1], &nw, tol, elim,
	      alim);
  goto L40;
 L20:
  if (az < *rl)
    {
      goto L30;
    }
  /*----------------------------------------------------------------------- 
   *    ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION 
   *----------------------------------------------------------------------- 
   */
  amos_zasyi (&znr, &zni, fnu, kode, &nn, &yr[1], &yi[1], &nw, rl, tol, elim,
	      alim);
  if (nw < 0)
    {
      goto L80;
    }
  goto L40;
 L30:
  /*----------------------------------------------------------------------- 
   *    MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION 
   *----------------------------------------------------------------------- 
   */
  amos_zmlri (&znr, &zni, fnu, kode, &nn, &yr[1], &yi[1], &nw, tol);
  if (nw < 0)
    {
      goto L80;
    }
 L40:
  /*----------------------------------------------------------------------- 
   *    ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION 
   *----------------------------------------------------------------------- 
   */
  amos_zbknu (&znr, &zni, fnu, kode, &c__1, cyr, cyi, &nw, tol, elim, alim);
  if (nw != 0)
    {
      goto L80;
    }
  fmr = (double) (*mr);
  sgn = - D_SIGN (pi, fmr);
  csgnr = 0.;
  csgni = sgn;
  if (*kode == 1)
    {
      goto L50;
    }
  yy = -zni;
  csgnr = -csgni * sin (yy);
  csgni *= cos (yy);
 L50:
  /*----------------------------------------------------------------------- 
   *    CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE 
   *    WHEN FNU IS LARGE 
   *----------------------------------------------------------------------- 
   */
  inu = (int) (*fnu);
  arg = (*fnu - (double) inu) * sgn;
  cspnr = cos (arg);
  cspni = sin (arg);
  if (inu % 2 == 0)
    {
      goto L60;
    }
  cspnr = -cspnr;
  cspni = -cspni;
 L60:
  c1r = cyr[0];
  c1i = cyi[0];
  c2r = yr[1];
  c2i = yi[1];
  if (*kode == 1)
    {
      goto L70;
    }
  iuf = 0;
  ascle = pnl_d1mach (1) * 1e3 / *tol;
  amos_zs1s2 (&znr, &zni, &c1r, &c1i, &c2r, &c2i, &nw, &ascle, alim, &iuf);
  *nz += nw;
 L70:
  yr[1] = cspnr * c1r - cspni * c1i + csgnr * c2r - csgni * c2i;
  yi[1] = cspnr * c1i + cspni * c1r + csgnr * c2i + csgni * c2r;
  return 0;
 L80:
  *nz = -1;
  if (nw == -2)
    {
      *nz = -2;
    }
  return 0;
}

