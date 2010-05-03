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
 *    zseri computes the i bessel function for real(z).ge.0.0 by 
 *    means of the power series for large cabs(z) in the 
 *    region cabs(z).le.2*sqrt(fnu+1). nz=0 is a normal return. 
 *    nz.gt.0 means that the last nz components were set to zero 
 *    due to underflow. nz.lt.0 means underflow occurred, but the 
 *    condition cabs(z).le.2*sqrt(fnu+1) was violated and the 
 *    computation must be completed in another routine with n=n-abs(nz). 
 *
 */

int
amos_zseri (double *zr, double *zi, double *fnu,const int *kode,const int *n,
	    double *yr, double *yi, int *nz, double *tol, double *elim,
	    double *alim)
{
  static const double zeror = 0.;
  static const double zeroi = 0.;
  static const double coner = 1.;
  static const double conei = 0.;

  int i__1;
  double dfnu;
  int idum;
  double atol, fnup;
  int i__, k, l, m, iflag;
  double s, coefi, ascle, coefr;
  double crscr;
  double aa;
  int ib;
  double ak;
  int il;
  double az;
  int nn;
  double wi[2];
  double rs, ss=0;
  int nw;
  double wr[2], s1i, s2i, s1r, s2r, cki, acz, arm, ckr, czi, hzi, raz, czr,
    sti, hzr, rzi, str, rzr, ak1i, ak1r, rtr1;

  --yi;
  --yr;

  *nz = 0;
  az = amos_azabs (zr, zi);
  if (az == 0.)
    {
      goto L160;
    }
  arm = pnl_d1mach (1) * 1e3;
  rtr1 = sqrt (arm);
  crscr = 1.;
  iflag = 0;
  if (az < arm)
    {
      goto L150;
    }
  hzr = *zr * .5;
  hzi = *zi * .5;
  czr = zeror;
  czi = zeroi;
  if (az <= rtr1)
    {
      goto L10;
    }
  amos_zmlt (&hzr, &hzi, &hzr, &hzi, &czr, &czi);
 L10:
  acz = amos_azabs (&czr, &czi);
  nn = *n;
  amos_azlog (&hzr, &hzi, &ckr, &cki, &idum);
 L20:
  dfnu = *fnu + (double) (nn - 1);
  fnup = dfnu + 1.;
  /*----------------------------------------------------------------------- 
   *    UNDERFLOW TEST 
   *----------------------------------------------------------------------- 
   */
  ak1r = ckr * dfnu;
  ak1i = cki * dfnu;
  ak = amos_dgamln (&fnup, &idum);
  ak1r -= ak;
  if (*kode == 2)
    {
      ak1r -= *zr;
    }
  if (ak1r > -(*elim))
    {
      goto L40;
    }
 L30:
  ++(*nz);
  yr[nn] = zeror;
  yi[nn] = zeroi;
  if (acz > dfnu)
    {
      goto L190;
    }
  --nn;
  if (nn == 0)
    {
      return 0;
    }
  goto L20;
 L40:
  if (ak1r > -(*alim))
    {
      goto L50;
    }
  iflag = 1;
  ss = 1. / *tol;
  crscr = *tol;
  ascle = arm * ss;
 L50:
  aa = exp (ak1r);
  if (iflag == 1)
    {
      aa *= ss;
    }
  coefr = aa * cos (ak1i);
  coefi = aa * sin (ak1i);
  atol = *tol * acz / fnup;
  il = MIN (2, nn);
  i__1 = il;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      dfnu = *fnu + (double) (nn - i__);
      fnup = dfnu + 1.;
      s1r = coner;
      s1i = conei;
      if (acz < *tol * fnup)
	{
	  goto L70;
	}
      ak1r = coner;
      ak1i = conei;
      ak = fnup + 2.;
      s = fnup;
      aa = 2.;
    L60:
      rs = 1. / s;
      str = ak1r * czr - ak1i * czi;
      sti = ak1r * czi + ak1i * czr;
      ak1r = str * rs;
      ak1i = sti * rs;
      s1r += ak1r;
      s1i += ak1i;
      s += ak;
      ak += 2.;
      aa = aa * acz * rs;
      if (aa > atol)
	{
	  goto L60;
	}
    L70:
      s2r = s1r * coefr - s1i * coefi;
      s2i = s1r * coefi + s1i * coefr;
      wr[i__ - 1] = s2r;
      wi[i__ - 1] = s2i;
      if (iflag == 0)
	{
	  goto L80;
	}
      amos_zuchk (&s2r, &s2i, &nw, &ascle, tol);
      if (nw != 0)
	{
	  goto L30;
	}
    L80:
      m = nn - i__ + 1;
      yr[m] = s2r * crscr;
      yi[m] = s2i * crscr;
      if (i__ == il)
	{
	  goto L90;
	}
      amos_zdiv (&coefr, &coefi, &hzr, &hzi, &str, &sti);
      coefr = str * dfnu;
      coefi = sti * dfnu;
    L90:
      ;
    }
  if (nn <= 2)
    {
      return 0;
    }
  k = nn - 2;
  ak = (double) k;
  raz = 1. / az;
  str = *zr * raz;
  sti = -(*zi) * raz;
  rzr = (str + str) * raz;
  rzi = (sti + sti) * raz;
  if (iflag == 1)
    {
      goto L120;
    }
  ib = 3;
 L100:
  i__1 = nn;
  for (i__ = ib; i__ <= i__1; ++i__)
    {
      yr[k] = (ak + *fnu) * (rzr * yr[k + 1] - rzi * yi[k + 1]) + yr[k + 2];
      yi[k] = (ak + *fnu) * (rzr * yi[k + 1] + rzi * yr[k + 1]) + yi[k + 2];
      ak += -1.;
      --k;
      /* L110: */
    }
  return 0;
  /*----------------------------------------------------------------------- 
   *    RECUR BACKWARD WITH SCALED VALUES 
   *----------------------------------------------------------------------- 
   */
 L120:
  /*----------------------------------------------------------------------- 
   *    EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE 
   *    UNDERFLOW LIMIT = ASCLE = D1MACH(1)*SS*1.0D+3 
   *----------------------------------------------------------------------- 
   */
  s1r = wr[0];
  s1i = wi[0];
  s2r = wr[1];
  s2i = wi[1];
  i__1 = nn;
  for (l = 3; l <= i__1; ++l)
    {
      ckr = s2r;
      cki = s2i;
      s2r = s1r + (ak + *fnu) * (rzr * ckr - rzi * cki);
      s2i = s1i + (ak + *fnu) * (rzr * cki + rzi * ckr);
      s1r = ckr;
      s1i = cki;
      ckr = s2r * crscr;
      cki = s2i * crscr;
      yr[k] = ckr;
      yi[k] = cki;
      ak += -1.;
      --k;
      if (amos_azabs (&ckr, &cki) > ascle)
	{
	  goto L140;
	}
      /* L130: */
    }
  return 0;
 L140:
  ib = l + 1;
  if (ib > nn)
    {
      return 0;
    }
  goto L100;
 L150:
  *nz = *n;
  if (*fnu == 0.)
    {
      --(*nz);
    }
 L160:
  yr[1] = zeror;
  yi[1] = zeroi;
  if (*fnu != 0.)
    {
      goto L170;
    }
  yr[1] = coner;
  yi[1] = conei;
 L170:
  if (*n == 1)
    {
      return 0;
    }
  i__1 = *n;
  for (i__ = 2; i__ <= i__1; ++i__)
    {
      yr[i__] = zeror;
      yi[i__] = zeroi;
    }
  return 0;
  /*----------------------------------------------------------------------- 
   *    RETURN WITH NZ.LT.0 IF CABS(Z*Z/4).GT.FNU+N-NZ-1 COMPLETE 
   *    THE CALCULATION IN CBINU WITH N=N-IABS(NZ) 
   *----------------------------------------------------------------------- 
   */
 L190:
  *nz = -(*nz);
  return 0;
}	
