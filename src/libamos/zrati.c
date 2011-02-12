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
 *    zrati computes ratios of i bessel functions by backward 
 *    recurrence.  the starting index is determined by forward 
 *    recurrence as described in j. res. of nat. bur. of standards-b, 
 *    mathematical sciences, vol 77b, p111-114, september, 1973, 
 *    bessel functions i and j of complex argument and int order, 
 *    by d. j. sookne. 
 *
 *        input      zr,zi,fnu are double precision 
 *          zr,zi  - z=cmplx(zr,zi),  -pi.lt.arg(z).le.pi 
 *          fnu    - order of initial j function, fnu.ge.0.0d0 
 *          n      - number of members of the sequence, n.ge.1 
 * 
 *        output     cyr,cyi are double precision 
 *          cyr,cyi- double precision vectors whose first n components 
 *                   contain real and imaginary parts for the sequence 
 *                   cy(i)=R(fnu+i-1,z)  i=1,...,n 
 */
int pnl_zrati (double *zr, double *zi, double *fnu,const int *n, double *cyr,
	    double *cyi, double *tol)
{
  /* Initialized data */

  static const double czeror = 0.;
  static const double czeroi = 0.;
  static const double coner = 1.;
  static const double conei = 0.;
  static const double rt2 = 1.41421356237309505;

  /* System generated locals */
  int i__1;
  double d__1;

  /* Local variables */
  double flam, dfnu, fdnu;
  int magz, idnu;
  double fnup;
  double test, test1;
  int i__, k;
  double amagz;
  int itime;
  double ak;
  int id, kk;
  double az, cdfnui, cdfnur, ap1, ap2, p1i, p2i, t1i, p1r, p2r, t1r, arg, rak,
    rho;
  int inu;
  double pti, tti, rzi, ptr, ttr, rzr, rap1;

  /* Parameter adjustments */
  --cyi;
  --cyr;

  /* Function Body */
  az = amos_azabs (zr, zi);
  inu = (int) (*fnu);
  idnu = inu + *n - 1;
  magz = (int) az;
  amagz = (double) (magz + 1);
  fdnu = (double) idnu;
  fnup = MAX (amagz, fdnu);
  id = idnu - magz - 1;
  itime = 1;
  k = 1;
  ptr = 1. / az;
  rzr = ptr * (*zr + *zr) * ptr;
  rzi = -ptr * (*zi + *zi) * ptr;
  t1r = rzr * fnup;
  t1i = rzi * fnup;
  p2r = -t1r;
  p2i = -t1i;
  p1r = coner;
  p1i = conei;
  t1r += rzr;
  t1i += rzi;
  if (id > 0)
    {
      id = 0;
    }
  ap2 = amos_azabs (&p2r, &p2i);
  ap1 = amos_azabs (&p1r, &p1i);
  /*----------------------------------------------------------------------- 
   *    THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNU 
   *    GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT 
   *    P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR 
   *    PREMATURELY. 
   *----------------------------------------------------------------------- 
   */
  arg = (ap2 + ap2) / (ap1 * *tol);
  test1 = sqrt (arg);
  test = test1;
  rap1 = 1. / ap1;
  p1r *= rap1;
  p1i *= rap1;
  p2r *= rap1;
  p2i *= rap1;
  ap2 *= rap1;
 L10:
  ++k;
  ap1 = ap2;
  ptr = p2r;
  pti = p2i;
  p2r = p1r - (t1r * ptr - t1i * pti);
  p2i = p1i - (t1r * pti + t1i * ptr);
  p1r = ptr;
  p1i = pti;
  t1r += rzr;
  t1i += rzi;
  ap2 = amos_azabs (&p2r, &p2i);
  if (ap1 <= test)
    {
      goto L10;
    }
  if (itime == 2)
    {
      goto L20;
    }
  ak = amos_azabs (&t1r, &t1i) * .5;
  flam = ak + sqrt (ak * ak - 1.);
  /*Computing MIN 
   */
  d__1 = ap2 / ap1;
  rho = MIN (d__1, flam);
  test = test1 * sqrt (rho / (rho * rho - 1.));
  itime = 2;
  goto L10;
 L20:
  kk = k + 1 - id;
  ak = (double) kk;
  t1r = ak;
  t1i = czeroi;
  dfnu = *fnu + (double) (*n - 1);
  p1r = 1. / ap2;
  p1i = czeroi;
  p2r = czeror;
  p2i = czeroi;
  i__1 = kk;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      ptr = p1r;
      pti = p1i;
      rap1 = dfnu + t1r;
      ttr = rzr * rap1;
      tti = rzi * rap1;
      p1r = ptr * ttr - pti * tti + p2r;
      p1i = ptr * tti + pti * ttr + p2i;
      p2r = ptr;
      p2i = pti;
      t1r -= coner;
      /* L30: */
    }
  if (p1r != czeror || p1i != czeroi)
    {
      goto L40;
    }
  p1r = *tol;
  p1i = *tol;
 L40:
  amos_zdiv (&p2r, &p2i, &p1r, &p1i, &cyr[*n], &cyi[*n]);
  if (*n == 1)
    {
      return 0;
    }
  k = *n - 1;
  ak = (double) k;
  t1r = ak;
  t1i = czeroi;
  cdfnur = *fnu * rzr;
  cdfnui = *fnu * rzi;
  i__1 = *n;
  for (i__ = 2; i__ <= i__1; ++i__)
    {
      ptr = cdfnur + (t1r * rzr - t1i * rzi) + cyr[k + 1];
      pti = cdfnui + (t1r * rzi + t1i * rzr) + cyi[k + 1];
      ak = amos_azabs (&ptr, &pti);
      if (ak != czeror)
	{
	  goto L50;
	}
      ptr = *tol;
      pti = *tol;
      ak = *tol * rt2;
    L50:
      rak = coner / ak;
      cyr[k] = rak * ptr * rak;
      cyi[k] = -rak * pti * rak;
      t1r -= coner;
      --k;
      /* L60: */
    }
  return 0;
}
