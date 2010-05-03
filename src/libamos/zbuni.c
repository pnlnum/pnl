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



static const int c__2 = 2;

/*
 *    zbuni computes the i bessel function for large cabs(z).gt. 
 *    fnul and fnu+n-1.lt.fnul. the order is increased from 
 *    fnu+n-1 greater than fnul by adding nui and computing 
 *    according to the uniform asymptotic expansion for i(fnu,z) 
 *    on iform=1 and the expansion for j(fnu,z) on iform=2 
 * 
 */

int amos_zbuni (double *zr, double *zi, double *fnu, int *kode,const int *n,
		double *yr, double *yi, int *nz, int *nui, int *nlast,
		double *fnul, double *tol, double *elim, double *alim)
{
  /* System generated locals */
  int i__1;

  /* Local variables */
  double dfnu, fnui;
  int i__, k, iflag;
  double ascle;
  double csclr, cscrr;
  int iform;
  double ax, ay;
  int nl, nw;
  double c1i, c1m, c1r, s1i, s2i, s1r, s2r, cyi[2], gnu, raz, cyr[2], sti,
    bry[3], rzi, str, rzr;

  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  *nz = 0;
  ax = fabs (*zr) * 1.7321;
  ay = fabs (*zi);
  iform = 1;
  if (ay > ax)
    {
      iform = 2;
    }
  if (*nui == 0)
    {
      goto L60;
    }
  fnui = (double) (*nui);
  dfnu = *fnu + (double) (*n - 1);
  gnu = dfnu + fnui;
  if (iform == 2)
    {
      goto L10;
    }
  /*----------------------------------------------------------------------- 
   *    ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN 
   *    -PI/3.LE.ARG(Z).LE.PI/3 
   *----------------------------------------------------------------------- 
   */
  amos_zuni1 (zr, zi, &gnu, kode, &c__2, cyr, cyi, &nw, nlast, fnul, tol,
	      elim, alim);
  goto L20;
 L10:
  /*----------------------------------------------------------------------- 
   *    ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU 
   *    APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I 
   *    AND HPI=PI/2 
   *----------------------------------------------------------------------- 
   */
  pnl_zuni2 (zr, zi, &gnu, kode, &c__2, cyr, cyi, &nw, nlast, fnul, tol,
	      elim, alim);
 L20:
  if (nw < 0)
    {
      goto L50;
    }
  if (nw != 0)
    {
      goto L90;
    }
  str = amos_azabs (cyr, cyi);
  /*---------------------------------------------------------------------- 
   *    SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED 
   *---------------------------------------------------------------------- 
   */
  bry[0] = pnl_d1mach (1) * 1e3 / *tol;
  bry[1] = 1. / bry[0];
  bry[2] = bry[1];
  iflag = 2;
  ascle = bry[1];
  csclr = 1.;
  if (str > bry[0])
    {
      goto L21;
    }
  iflag = 1;
  ascle = bry[0];
  csclr = 1. / *tol;
  goto L25;
 L21:
  if (str < bry[1])
    {
      goto L25;
    }
  iflag = 3;
  ascle = bry[2];
  csclr = *tol;
 L25:
  cscrr = 1. / csclr;
  s1r = cyr[1] * csclr;
  s1i = cyi[1] * csclr;
  s2r = cyr[0] * csclr;
  s2i = cyi[0] * csclr;
  raz = 1. / amos_azabs (zr, zi);
  str = *zr * raz;
  sti = -(*zi) * raz;
  rzr = (str + str) * raz;
  rzi = (sti + sti) * raz;
  i__1 = *nui;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      str = s2r;
      sti = s2i;
      s2r = (dfnu + fnui) * (rzr * str - rzi * sti) + s1r;
      s2i = (dfnu + fnui) * (rzr * sti + rzi * str) + s1i;
      s1r = str;
      s1i = sti;
      fnui += -1.;
      if (iflag >= 3)
	{
	  goto L30;
	}
      str = s2r * cscrr;
      sti = s2i * cscrr;
      c1r = fabs (str);
      c1i = fabs (sti);
      c1m = MAX (c1r, c1i);
      if (c1m <= ascle)
	{
	  goto L30;
	}
      ++iflag;
      ascle = bry[iflag - 1];
      s1r *= cscrr;
      s1i *= cscrr;
      s2r = str;
      s2i = sti;
      csclr *= *tol;
      cscrr = 1. / csclr;
      s1r *= csclr;
      s1i *= csclr;
      s2r *= csclr;
      s2i *= csclr;
    L30:
      ;
    }
  yr[*n] = s2r * cscrr;
  yi[*n] = s2i * cscrr;
  if (*n == 1)
    {
      return 0;
    }
  nl = *n - 1;
  fnui = (double) nl;
  k = nl;
  i__1 = nl;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      str = s2r;
      sti = s2i;
      s2r = (*fnu + fnui) * (rzr * str - rzi * sti) + s1r;
      s2i = (*fnu + fnui) * (rzr * sti + rzi * str) + s1i;
      s1r = str;
      s1i = sti;
      str = s2r * cscrr;
      sti = s2i * cscrr;
      yr[k] = str;
      yi[k] = sti;
      fnui += -1.;
      --k;
      if (iflag >= 3)
	{
	  goto L40;
	}
      c1r = fabs (str);
      c1i = fabs (sti);
      c1m = MAX (c1r, c1i);
      if (c1m <= ascle)
	{
	  goto L40;
	}
      ++iflag;
      ascle = bry[iflag - 1];
      s1r *= cscrr;
      s1i *= cscrr;
      s2r = str;
      s2i = sti;
      csclr *= *tol;
      cscrr = 1. / csclr;
      s1r *= csclr;
      s1i *= csclr;
      s2r *= csclr;
      s2i *= csclr;
    L40:
      ;
    }
  return 0;
 L50:
  *nz = -1;
  if (nw == -2)
    {
      *nz = -2;
    }
  return 0;
 L60:
  if (iform == 2)
    {
      goto L70;
    }
  /*----------------------------------------------------------------------- 
   *    ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN 
   *    -PI/3.LE.ARG(Z).LE.PI/3 
   *----------------------------------------------------------------------- 
   */
  amos_zuni1 (zr, zi, fnu, kode, n, &yr[1], &yi[1], &nw, nlast, fnul, tol,
	      elim, alim);
  goto L80;
 L70:
  /*----------------------------------------------------------------------- 
   *    ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU 
   *    APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I 
   *    AND HPI=PI/2 
   *----------------------------------------------------------------------- 
   */
  pnl_zuni2 (zr, zi, fnu, kode, n, &yr[1], &yi[1], &nw, nlast, fnul, tol,
	      elim, alim);
 L80:
  if (nw < 0)
    {
      goto L50;
    }
  *nz = nw;
  return 0;
 L90:
  *nlast = *n;
  return 0;
}	
