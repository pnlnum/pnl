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
 * 
 *    zbknu computes the k bessel function in the right half z plane. 
 * 
 */

/* Table of constant values */

static const double c_b10 = .5;
static const double c_b11 = 0.;

/* Subroutine */ int
amos_zbknu (double *zr, double *zi, double *fnu,const int *kode,const int *n,
	    double *yr, double *yi, int *nz, double *tol, double *elim,
	    double *alim)
{
  /* Initialized data */

  static const int kmax = 30;
  static const double spi = 1.90985931710274403;
  static const double hpi = 1.57079632679489662;
  static const double fpi = 1.89769999331517738;
  static const double tth = .666666666666666666;
  static const double cc[8] =
    { .577215664901532861, -.0420026350340952355, -.0421977345555443367,
      .00721894324666309954, -2.15241674114950973e-4, -2.01348547807882387e-5,
      1.13302723198169588e-6, 6.11609510448141582e-9 };
  static const double czeror = 0.;
  static const double czeroi = 0.;
  static const double coner = 1.;
  static const double conei = 0.;
  static const double ctwor = 2.;
  static const double r1 = 2.;
  static const double dpi = 3.14159265358979324;
  static const double rthpi = 1.25331413731550025;

  /* System generated locals */
  int i__1;
  double d__1;

  /* Local variables */
  double cchi, cchr, alas, cshi;
  int inub, idum;
  double cshr, fmui, rcaz, csrr[3], cssr[3], fmur;
  double smui, smur;
  int i__, j, k, iflag;
  double s;
  int kflag;
  double coefi;
  int koded;
  double ascle, coefr, helim;
  double celmr, csclr, crscr;
  double a1, a2, etest;
  double g1, g2;
  double t1, t2, aa, bb, fc, ak, bk;
  int ic;
  double fi, fk, as;
  int kk;
  double fr, pi, qi, tm, pr, qr;
  int nw;
  double p1i, p2i, s1i, s2i, p2m, p1r, p2r;
  double s1r, s2r, cbi, cbr, cki=0, caz, csi, ckr=0, fhs, fks, rak, czi, dnu, csr,
    elm, zdi, bry[3], pti, czr, sti, zdr, cyr[2], rzi, ptr, cyi[2];
  int inu;
  double str, rzr, dnu2=0.0;

  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  /* 
   */
  caz = amos_azabs (zr, zi);
  csclr = 1. / *tol;
  crscr = *tol;
  cssr[0] = csclr;
  cssr[1] = 1.;
  cssr[2] = crscr;
  csrr[0] = crscr;
  csrr[1] = 1.;
  csrr[2] = csclr;
  bry[0] = pnl_d1mach (1) * 1e3 / *tol;
  bry[1] = 1. / bry[0];
  bry[2] = pnl_d1mach (2);
  *nz = 0;
  iflag = 0;
  koded = *kode;
  rcaz = 1. / caz;
  str = *zr * rcaz;
  sti = -(*zi) * rcaz;
  rzr = (str + str) * rcaz;
  rzi = (sti + sti) * rcaz;
  inu = (int) (*fnu + .5);
  dnu = *fnu - (double) inu;
  if (fabs (dnu) == .5)
    {
      goto L110;
    }
  dnu2 = 0.;
  if (fabs (dnu) > *tol)
    {
      dnu2 = dnu * dnu;
    }
  if (caz > r1)
    {
      goto L110;
    }
  /*----------------------------------------------------------------------- 
   *    SERIES FOR CABS(Z).LE.R1 
   *----------------------------------------------------------------------- 
   */
  fc = 1.;
  amos_azlog (&rzr, &rzi, &smur, &smui, &idum);
  fmur = smur * dnu;
  fmui = smui * dnu;
  amos_zshch (&fmur, &fmui, &cshr, &cshi, &cchr, &cchi);
  if (dnu == 0.)
    {
      goto L10;
    }
  fc = dnu * dpi;
  fc /= sin (fc);
  smur = cshr / dnu;
  smui = cshi / dnu;
 L10:
  a2 = dnu + 1.;
  /*----------------------------------------------------------------------- 
   *    GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU) 
   *----------------------------------------------------------------------- 
   */
  t2 = exp (-amos_dgamln (&a2, &idum));
  t1 = 1. / (t2 * fc);
  if (fabs (dnu) > .1)
    {
      goto L40;
    }
  /*----------------------------------------------------------------------- 
   *    SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU) 
   *----------------------------------------------------------------------- 
   */
  ak = 1.;
  s = cc[0];
  for (k = 2; k <= 8; ++k)
    {
      ak *= dnu2;
      tm = cc[k - 1] * ak;
      s += tm;
      if (fabs (tm) < *tol)
	{
	  goto L30;
	}
      /* L20: */
    }
 L30:
  g1 = -s;
  goto L50;
 L40:
  g1 = (t1 - t2) / (dnu + dnu);
 L50:
  g2 = (t1 + t2) * .5;
  fr = fc * (cchr * g1 + smur * g2);
  fi = fc * (cchi * g1 + smui * g2);
  amos_azexp (&fmur, &fmui, &str, &sti);
  pr = str * .5 / t2;
  pi = sti * .5 / t2;
  amos_zdiv (&c_b10, &c_b11, &str, &sti, &ptr, &pti);
  qr = ptr / t1;
  qi = pti / t1;
  s1r = fr;
  s1i = fi;
  s2r = pr;
  s2i = pi;
  ak = 1.;
  a1 = 1.;
  ckr = coner;
  cki = conei;
  bk = 1. - dnu2;
  if (inu > 0 || *n > 1)
    {
      goto L80;
    }
  /*----------------------------------------------------------------------- 
   *    GENERATE K(FNU,Z), 0.0D0 .LE. FNU .LT. 0.5D0 AND N=1 
   *----------------------------------------------------------------------- 
   */
  if (caz < *tol)
    {
      goto L70;
    }
  amos_zmlt (zr, zi, zr, zi, &czr, &czi);
  czr *= .25;
  czi *= .25;
  t1 = caz * .25 * caz;
 L60:
  fr = (fr * ak + pr + qr) / bk;
  fi = (fi * ak + pi + qi) / bk;
  str = 1. / (ak - dnu);
  pr *= str;
  pi *= str;
  str = 1. / (ak + dnu);
  qr *= str;
  qi *= str;
  str = ckr * czr - cki * czi;
  rak = 1. / ak;
  cki = (ckr * czi + cki * czr) * rak;
  ckr = str * rak;
  s1r = ckr * fr - cki * fi + s1r;
  s1i = ckr * fi + cki * fr + s1i;
  a1 = a1 * t1 * rak;
  bk = bk + ak + ak + 1.;
  ak += 1.;
  if (a1 > *tol)
    {
      goto L60;
    }
 L70:
  yr[1] = s1r;
  yi[1] = s1i;
  if (koded == 1)
    {
      return 0;
    }
  amos_azexp (zr, zi, &str, &sti);
  amos_zmlt (&s1r, &s1i, &str, &sti, &yr[1], &yi[1]);
  return 0;
  /*----------------------------------------------------------------------- 
   *    GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE 
   *----------------------------------------------------------------------- 
   */
 L80:
  if (caz < *tol)
    {
      goto L100;
    }
  amos_zmlt (zr, zi, zr, zi, &czr, &czi);
  czr *= .25;
  czi *= .25;
  t1 = caz * .25 * caz;
 L90:
  fr = (fr * ak + pr + qr) / bk;
  fi = (fi * ak + pi + qi) / bk;
  str = 1. / (ak - dnu);
  pr *= str;
  pi *= str;
  str = 1. / (ak + dnu);
  qr *= str;
  qi *= str;
  str = ckr * czr - cki * czi;
  rak = 1. / ak;
  cki = (ckr * czi + cki * czr) * rak;
  ckr = str * rak;
  s1r = ckr * fr - cki * fi + s1r;
  s1i = ckr * fi + cki * fr + s1i;
  str = pr - fr * ak;
  sti = pi - fi * ak;
  s2r = ckr * str - cki * sti + s2r;
  s2i = ckr * sti + cki * str + s2i;
  a1 = a1 * t1 * rak;
  bk = bk + ak + ak + 1.;
  ak += 1.;
  if (a1 > *tol)
    {
      goto L90;
    }
 L100:
  kflag = 2;
  a1 = *fnu + 1.;
  ak = a1 * fabs (smur);
  if (ak > *alim)
    {
      kflag = 3;
    }
  str = cssr[kflag - 1];
  p2r = s2r * str;
  p2i = s2i * str;
  amos_zmlt (&p2r, &p2i, &rzr, &rzi, &s2r, &s2i);
  s1r *= str;
  s1i *= str;
  if (koded == 1)
    {
      goto L210;
    }
  amos_azexp (zr, zi, &fr, &fi);
  amos_zmlt (&s1r, &s1i, &fr, &fi, &s1r, &s1i);
  amos_zmlt (&s2r, &s2i, &fr, &fi, &s2r, &s2i);
  goto L210;
  /*----------------------------------------------------------------------- 
   *    IFLAG=0 MEANS NO UNDERFLOW OCCURRED 
   *    IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH 
   *    KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD 
   *    RECURSION 
   *----------------------------------------------------------------------- 
   */
 L110:
  amos_azsqrt (zr, zi, &str, &sti);
  amos_zdiv (&rthpi, &czeroi, &str, &sti, &coefr, &coefi);
  kflag = 2;
  if (koded == 2)
    {
      goto L120;
    }
  if (*zr > *alim)
    {
      goto L290;
    }
  /*    BLANK LINE 
   */
  str = exp (-(*zr)) * cssr[kflag - 1];
  sti = -str * sin (*zi);
  str *= cos (*zi);
  amos_zmlt (&coefr, &coefi, &str, &sti, &coefr, &coefi);
 L120:
  if (fabs (dnu) == .5)
    {
      goto L300;
    }
  /*----------------------------------------------------------------------- 
   *    MILLER ALGORITHM FOR CABS(Z).GT.R1 
   *----------------------------------------------------------------------- 
   */
  ak = cos (dpi * dnu);
  ak = fabs (ak);
  if (ak == czeror)
    {
      goto L300;
    }
  fhs = (d__1 = .25 - dnu2, fabs (d__1));
  if (fhs == czeror)
    {
      goto L300;
    }
  /*----------------------------------------------------------------------- 
   *    COMPUTE R2=F(E). IF CABS(Z).GE.R2, USE FORWARD RECURRENCE TO 
   *    DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON 
   *    12.LE.E.LE.60. E IS COMPUTED FROM 2**(-E)=B**(1-I1MACH(14))= 
   *    TOL WHERE B IS THE BASE OF THE ARITHMETIC. 
   *----------------------------------------------------------------------- 
   */
  t1 = (double) (amos_i1mach (14) - 1);
  t1 = t1 * pnl_d1mach (5) * 3.321928094;
  t1 = MAX (t1, 12.);
  t1 = MIN (t1, 60.);
  t2 = tth * t1 - 6.;
  if (*zr != 0.)
    {
      goto L130;
    }
  t1 = hpi;
  goto L140;
 L130:
  t1 = atan (*zi / *zr);
  t1 = fabs (t1);
 L140:
  if (t2 > caz)
    {
      goto L170;
    }
  /*----------------------------------------------------------------------- 
   *    FORWARD RECURRENCE LOOP WHEN CABS(Z).GE.R2 
   *----------------------------------------------------------------------- 
   */
  etest = ak / (dpi * caz * *tol);
  fk = coner;
  if (etest < coner)
    {
      goto L180;
    }
  fks = ctwor;
  ckr = caz + caz + ctwor;
  p1r = czeror;
  p2r = coner;
  i__1 = kmax;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      ak = fhs / fks;
      cbr = ckr / (fk + coner);
      ptr = p2r;
      p2r = cbr * p2r - p1r * ak;
      p1r = ptr;
      ckr += ctwor;
      fks = fks + fk + fk + ctwor;
      fhs = fhs + fk + fk;
      fk += coner;
      str = fabs (p2r) * fk;
      if (etest < str)
	{
	  goto L160;
	}
      /* L150: */
    }
  goto L310;
 L160:
  fk += spi * t1 * sqrt (t2 / caz);
  fhs = (d__1 = .25 - dnu2, fabs (d__1));
  goto L180;
 L170:
  /*----------------------------------------------------------------------- 
   *    COMPUTE BACKWARD INDEX K FOR CABS(Z).LT.R2 
   *----------------------------------------------------------------------- 
   */
  a2 = sqrt (caz);
  ak = fpi * ak / (*tol * sqrt (a2));
  aa = t1 * 3. / (caz + 1.);
  bb = t1 * 14.7 / (caz + 28.);
  ak = (log (ak) + caz * cos (aa) / (caz * .008 + 1.)) / cos (bb);
  fk = ak * .12125 * ak / caz + 1.5;
 L180:
  /*----------------------------------------------------------------------- 
   *    BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM 
   *----------------------------------------------------------------------- 
   */
  k = (int) fk;
  fk = (double) k;
  fks = fk * fk;
  p1r = czeror;
  p1i = czeroi;
  p2r = *tol;
  p2i = czeroi;
  csr = p2r;
  csi = p2i;
  i__1 = k;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      a1 = fks - fk;
      ak = (fks + fk) / (a1 + fhs);
      rak = 2. / (fk + coner);
      cbr = (fk + *zr) * rak;
      cbi = *zi * rak;
      ptr = p2r;
      pti = p2i;
      p2r = (ptr * cbr - pti * cbi - p1r) * ak;
      p2i = (pti * cbr + ptr * cbi - p1i) * ak;
      p1r = ptr;
      p1i = pti;
      csr += p2r;
      csi += p2i;
      fks = a1 - fk + coner;
      fk -= coner;
      /* L190: */
    }
  /*----------------------------------------------------------------------- 
   *    COMPUTE (P2/CS)=(P2/CABS(CS))*(CONJG(CS)/CABS(CS)) FOR BETTER 
   *    SCALING 
   *----------------------------------------------------------------------- 
   */
  tm = amos_azabs (&csr, &csi);
  ptr = 1. / tm;
  s1r = p2r * ptr;
  s1i = p2i * ptr;
  csr *= ptr;
  csi = -csi * ptr;
  amos_zmlt (&coefr, &coefi, &s1r, &s1i, &str, &sti);
  amos_zmlt (&str, &sti, &csr, &csi, &s1r, &s1i);
  if (inu > 0 || *n > 1)
    {
      goto L200;
    }
  zdr = *zr;
  zdi = *zi;
  if (iflag == 1)
    {
      goto L270;
    }
  goto L240;
 L200:
  /*----------------------------------------------------------------------- 
   *    COMPUTE P1/P2=(P1/CABS(P2)*CONJG(P2)/CABS(P2) FOR SCALING 
   *----------------------------------------------------------------------- 
   */
  tm = amos_azabs (&p2r, &p2i);
  ptr = 1. / tm;
  p1r *= ptr;
  p1i *= ptr;
  p2r *= ptr;
  p2i = -p2i * ptr;
  amos_zmlt (&p1r, &p1i, &p2r, &p2i, &ptr, &pti);
  str = dnu + .5 - ptr;
  sti = -pti;
  amos_zdiv (&str, &sti, zr, zi, &str, &sti);
  str += 1.;
  amos_zmlt (&str, &sti, &s1r, &s1i, &s2r, &s2i);
  /*----------------------------------------------------------------------- 
   *    FORWARD RECURSION ON THE THREE TERM RECURSION WITH RELATION WITH 
   *    SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3 
   *----------------------------------------------------------------------- 
   */
 L210:
  str = dnu + 1.;
  ckr = str * rzr;
  cki = str * rzi;
  if (*n == 1)
    {
      --inu;
    }
  if (inu > 0)
    {
      goto L220;
    }
  if (*n > 1)
    {
      goto L215;
    }
  s1r = s2r;
  s1i = s2i;
 L215:
  zdr = *zr;
  zdi = *zi;
  if (iflag == 1)
    {
      goto L270;
    }
  goto L240;
 L220:
  inub = 1;
  if (iflag == 1)
    {
      goto L261;
    }
 L225:
  p1r = csrr[kflag - 1];
  ascle = bry[kflag - 1];
  i__1 = inu;
  for (i__ = inub; i__ <= i__1; ++i__)
    {
      str = s2r;
      sti = s2i;
      s2r = ckr * str - cki * sti + s1r;
      s2i = ckr * sti + cki * str + s1i;
      s1r = str;
      s1i = sti;
      ckr += rzr;
      cki += rzi;
      if (kflag >= 3)
	{
	  goto L230;
	}
      p2r = s2r * p1r;
      p2i = s2i * p1r;
      str = fabs (p2r);
      sti = fabs (p2i);
      p2m = MAX (str, sti);
      if (p2m <= ascle)
	{
	  goto L230;
	}
      ++kflag;
      ascle = bry[kflag - 1];
      s1r *= p1r;
      s1i *= p1r;
      s2r = p2r;
      s2i = p2i;
      str = cssr[kflag - 1];
      s1r *= str;
      s1i *= str;
      s2r *= str;
      s2i *= str;
      p1r = csrr[kflag - 1];
    L230:
      ;
    }
  if (*n != 1)
    {
      goto L240;
    }
  s1r = s2r;
  s1i = s2i;
 L240:
  str = csrr[kflag - 1];
  yr[1] = s1r * str;
  yi[1] = s1i * str;
  if (*n == 1)
    {
      return 0;
    }
  yr[2] = s2r * str;
  yi[2] = s2i * str;
  if (*n == 2)
    {
      return 0;
    }
  kk = 2;
 L250:
  ++kk;
  if (kk > *n)
    {
      return 0;
    }
  p1r = csrr[kflag - 1];
  ascle = bry[kflag - 1];
  i__1 = *n;
  for (i__ = kk; i__ <= i__1; ++i__)
    {
      p2r = s2r;
      p2i = s2i;
      s2r = ckr * p2r - cki * p2i + s1r;
      s2i = cki * p2r + ckr * p2i + s1i;
      s1r = p2r;
      s1i = p2i;
      ckr += rzr;
      cki += rzi;
      p2r = s2r * p1r;
      p2i = s2i * p1r;
      yr[i__] = p2r;
      yi[i__] = p2i;
      if (kflag >= 3)
	{
	  goto L260;
	}
      str = fabs (p2r);
      sti = fabs (p2i);
      p2m = MAX (str, sti);
      if (p2m <= ascle)
	{
	  goto L260;
	}
      ++kflag;
      ascle = bry[kflag - 1];
      s1r *= p1r;
      s1i *= p1r;
      s2r = p2r;
      s2i = p2i;
      str = cssr[kflag - 1];
      s1r *= str;
      s1i *= str;
      s2r *= str;
      s2i *= str;
      p1r = csrr[kflag - 1];
    L260:
      ;
    }
  return 0;
  /*----------------------------------------------------------------------- 
   *    IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW 
   *----------------------------------------------------------------------- 
   */
 L261:
  helim = *elim * .5;
  elm = exp (-(*elim));
  celmr = elm;
  ascle = bry[0];
  zdr = *zr;
  zdi = *zi;
  ic = -1;
  j = 2;
  i__1 = inu;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      str = s2r;
      sti = s2i;
      s2r = str * ckr - sti * cki + s1r;
      s2i = sti * ckr + str * cki + s1i;
      s1r = str;
      s1i = sti;
      ckr += rzr;
      cki += rzi;
      as = amos_azabs (&s2r, &s2i);
      alas = log (as);
      p2r = -zdr + alas;
      if (p2r < -(*elim))
	{
	  goto L263;
	}
      amos_azlog (&s2r, &s2i, &str, &sti, &idum);
      p2r = -zdr + str;
      p2i = -zdi + sti;
      p2m = exp (p2r) / *tol;
      p1r = p2m * cos (p2i);
      p1i = p2m * sin (p2i);
      amos_zuchk (&p1r, &p1i, &nw, &ascle, tol);
      if (nw != 0)
	{
	  goto L263;
	}
      j = 3 - j;
      cyr[j - 1] = p1r;
      cyi[j - 1] = p1i;
      if (ic == i__ - 1)
	{
	  goto L264;
	}
      ic = i__;
      goto L262;
    L263:
      if (alas < helim)
	{
	  goto L262;
	}
      zdr -= *elim;
      s1r *= celmr;
      s1i *= celmr;
      s2r *= celmr;
      s2i *= celmr;
    L262:
      ;
    }
  if (*n != 1)
    {
      goto L270;
    }
  s1r = s2r;
  s1i = s2i;
  goto L270;
 L264:
  kflag = 1;
  inub = i__ + 1;
  s2r = cyr[j - 1];
  s2i = cyi[j - 1];
  j = 3 - j;
  s1r = cyr[j - 1];
  s1i = cyi[j - 1];
  if (inub <= inu)
    {
      goto L225;
    }
  if (*n != 1)
    {
      goto L240;
    }
  s1r = s2r;
  s1i = s2i;
  goto L240;
 L270:
  yr[1] = s1r;
  yi[1] = s1i;
  if (*n == 1)
    {
      goto L280;
    }
  yr[2] = s2r;
  yi[2] = s2i;
 L280:
  ascle = bry[0];
  amos_zkscl (&zdr, &zdi, fnu, n, &yr[1], &yi[1], nz, &rzr, &rzi, &ascle, tol,  elim);
  inu = *n - *nz;
  if (inu <= 0)
    {
      return 0;
    }
  kk = *nz + 1;
  s1r = yr[kk];
  s1i = yi[kk];
  yr[kk] = s1r * csrr[0];
  yi[kk] = s1i * csrr[0];
  if (inu == 1)
    {
      return 0;
    }
  kk = *nz + 2;
  s2r = yr[kk];
  s2i = yi[kk];
  yr[kk] = s2r * csrr[0];
  yi[kk] = s2i * csrr[0];
  if (inu == 2)
    {
      return 0;
    }
  t2 = *fnu + (double) (kk - 1);
  ckr = t2 * rzr;
  cki = t2 * rzi;
  kflag = 1;
  goto L250;
 L290:
  /*----------------------------------------------------------------------- 
   *    SCALE BY DEXP(Z), IFLAG = 1 CASES 
   *----------------------------------------------------------------------- 
   */
  koded = 2;
  iflag = 1;
  kflag = 2;
  goto L120;
  /*----------------------------------------------------------------------- 
   *    FNU=HALF ODD INT CASE, DNU=-0.5 
   *----------------------------------------------------------------------- 
   */
 L300:
  s1r = coefr;
  s1i = coefi;
  s2r = coefr;
  s2i = coefi;
  goto L210;
  /* 
   * 
   */
 L310:
  *nz = -2;
  return 0;
}				/* zbknu_ */
