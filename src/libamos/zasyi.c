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
 *    zasyi computes the i bessel function for real(z).ge.0.0 by 
 *    means of the asymptotic expansion for large cabs(z) in the 
 *    region cabs(z).gt.max(rl,fnu*fnu/2). nz=0 is a normal return. 
 *    nz.lt.0 indicates an overflow on kode=1. 
 *
 */

int amos_zasyi (double *zr, double *zi, double *fnu,const int *kode,const  int *n,
		double *yr, double *yi, int *nz, double *rl, double *tol,
		double *elim, double *alim)
{
  static const double pi = 3.14159265358979324;
  static const double rtpi = .159154943091895336;
  static const double zeror = 0.;
  static const double zeroi = 0.;
  static const double coner = 1.;
  static const double conei = 0.;

  /* System generated locals */
  int i__1, i__2;
  double d__1, d__2;

  /* Local variables */
  double dfnu, atol;
  int i__, j, k, m;
  double s;
  int koded;
  double aa, bb;
  int ib;
  double ak, bk;
  int il, jl;
  double az;
  int nn;
  double p1i, s2i, p1r;
  double s2r, cki, dki, fdn, arg, aez, arm, ckr, dkr, czi, ezi, sgn;
  int inu;
  double raz, czr, ezr, sqk, sti, rzi, tzi, str, rzr, tzr, ak1i, ak1r, cs1i,
    cs2i, cs1r, cs2r, dnu2, rtr1;

  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  /* 
   */
  *nz = 0;
  az = amos_azabs (zr, zi);
  arm = pnl_d1mach (1) * 1e3;
  rtr1 = sqrt (arm);
  il = MIN (2, *n);
  dfnu = *fnu + (double) (*n - il);
  /*----------------------------------------------------------------------- 
   *    OVERFLOW TEST 
   *----------------------------------------------------------------------- 
   */
  raz = 1. / az;
  str = *zr * raz;
  sti = -(*zi) * raz;
  ak1r = rtpi * str * raz;
  ak1i = rtpi * sti * raz;
  amos_azsqrt (&ak1r, &ak1i, &ak1r, &ak1i);
  czr = *zr;
  czi = *zi;
  if (*kode != 2)
    {
      goto L10;
    }
  czr = zeror;
  czi = *zi;
 L10:
  if (fabs (czr) > *elim)
    {
      goto L100;
    }
  dnu2 = dfnu + dfnu;
  koded = 1;
  if (fabs (czr) > *alim && *n > 2)
    {
      goto L20;
    }
  koded = 0;
  amos_azexp (&czr, &czi, &str, &sti);
  amos_zmlt (&ak1r, &ak1i, &str, &sti, &ak1r, &ak1i);
 L20:
  fdn = 0.;
  if (dnu2 > rtr1)
    {
      fdn = dnu2 * dnu2;
    }
  ezr = *zr * 8.;
  ezi = *zi * 8.;
  /*----------------------------------------------------------------------- 
   *    WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE 
   *    FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE 
   *    EXPANSION FOR THE IMAGINARY PART. 
   *----------------------------------------------------------------------- 
   */
  aez = az * 8.;
  s = *tol / aez;
  jl = (int) (*rl + *rl) + 2;
  p1r = zeror;
  p1i = zeroi;
  if (*zi == 0.)
    {
      goto L30;
    }
  /*----------------------------------------------------------------------- 
   *    CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF 
   *    SIGNIFICANCE WHEN FNU OR N IS LARGE 
   *----------------------------------------------------------------------- 
   */
  inu = (int) (*fnu);
  arg = (*fnu - (double) inu) * pi;
  inu = inu + *n - il;
  ak = -sin (arg);
  bk = cos (arg);
  if (*zi < 0.)
    {
      bk = -bk;
    }
  p1r = ak;
  p1i = bk;
  if (inu % 2 == 0)
    {
      goto L30;
    }
  p1r = -p1r;
  p1i = -p1i;
 L30:
  i__1 = il;
  for (k = 1; k <= i__1; ++k)
    {
      sqk = fdn - 1.;
      atol = s * fabs (sqk);
      sgn = 1.;
      cs1r = coner;
      cs1i = conei;
      cs2r = coner;
      cs2i = conei;
      ckr = coner;
      cki = conei;
      ak = 0.;
      aa = 1.;
      bb = aez;
      dkr = ezr;
      dki = ezi;
      i__2 = jl;
      for (j = 1; j <= i__2; ++j)
	{
	  amos_zdiv (&ckr, &cki, &dkr, &dki, &str, &sti);
	  ckr = str * sqk;
	  cki = sti * sqk;
	  cs2r += ckr;
	  cs2i += cki;
	  sgn = -sgn;
	  cs1r += ckr * sgn;
	  cs1i += cki * sgn;
	  dkr += ezr;
	  dki += ezi;
	  aa = aa * fabs (sqk) / bb;
	  bb += aez;
	  ak += 8.;
	  sqk -= ak;
	  if (aa <= atol)
	    {
	      goto L50;
	    }
	  /* L40: */
	}
      goto L110;
    L50:
      s2r = cs1r;
      s2i = cs1i;
      if (*zr + *zr >= *elim)
	{
	  goto L60;
	}
      tzr = *zr + *zr;
      tzi = *zi + *zi;
      d__1 = -tzr;
      d__2 = -tzi;
      amos_azexp (&d__1, &d__2, &str, &sti);
      amos_zmlt (&str, &sti, &p1r, &p1i, &str, &sti);
      amos_zmlt (&str, &sti, &cs2r, &cs2i, &str, &sti);
      s2r += str;
      s2i += sti;
    L60:
      fdn = fdn + dfnu * 8. + 4.;
      p1r = -p1r;
      p1i = -p1i;
      m = *n - il + k;
      yr[m] = s2r * ak1r - s2i * ak1i;
      yi[m] = s2r * ak1i + s2i * ak1r;
      /* L70: */
    }
  if (*n <= 2)
    {
      return 0;
    }
  nn = *n;
  k = nn - 2;
  ak = (double) k;
  str = *zr * raz;
  sti = -(*zi) * raz;
  rzr = (str + str) * raz;
  rzi = (sti + sti) * raz;
  ib = 3;
  i__1 = nn;
  for (i__ = ib; i__ <= i__1; ++i__)
    {
      yr[k] = (ak + *fnu) * (rzr * yr[k + 1] - rzi * yi[k + 1]) + yr[k + 2];
      yi[k] = (ak + *fnu) * (rzr * yi[k + 1] + rzi * yr[k + 1]) + yi[k + 2];
      ak += -1.;
      --k;
      /* L80: */
    }
  if (koded == 0)
    {
      return 0;
    }
  amos_azexp (&czr, &czi, &ckr, &cki);
  i__1 = nn;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      str = yr[i__] * ckr - yi[i__] * cki;
      yi[i__] = yr[i__] * cki + yi[i__] * ckr;
      yr[i__] = str;
      /* L90: */
    }
  return 0;
 L100:
  *nz = -1;
  return 0;
 L110:
  *nz = -2;
  return 0;
}
