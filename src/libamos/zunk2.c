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



static const int c__1 = 1;
static const int c__2 = 2;
static const int c__0 = 0;

/*
 *    zunk2 computes k(fnu,z) and its analytic continuation from the 
 *    right half plane to the left half plane by means of the 
 *    uniform asymptotic expansions for h(kind,fnu,zn) and j(fnu,zn) 
 *    where zn is in the right half plane, kind=(3-mr)/2, mr=+1 or 
 *    -1. here zn=zr*i or -zr*i where zr=z if z is in the right 
 *    half plane or zr=-z if z is in the left half plane. mr indic- 
 *    ates the direction of rotation for analytic continuation. 
 *    nz=-1 means an overflow will occur 
 *
 */

int amos_zunk2 (double *zr, double *zi, double *fnu, int *kode, int *mr,const int *n,
		double *yr, double *yi, int *nz, double *tol, double *elim,
		double *alim)
{

  static const double zeror = 0.;
  static const double aic = 1.26551212348464539;
  static const double cipr[4] = { 1., 0., -1., 0. };
  static const double cipi[4] = { 0., -1., 0., 1. };
  static const double zeroi = 0.;
  static const double coner = 1.;
  static const double cr1r = 1.;
  static const double cr1i = 1.73205080756887729;
  static const double cr2r = -.5;
  static const double cr2i = -.866025403784438647;
  static const double hpi = 1.57079632679489662;
  static const double pi = 3.14159265358979324;

  /* System generated locals */
  int i__1;

  /* Local variables */
  double daii, aarg;
  int ndai;
  double dair, aphi, argi[2], cscl, phii[2], crsc, argr[2];
  int idum;
  double phir[2], csrr[3], cssr[3], rast, razr;
  int i__, k, j, iflag=1, kflag=1;
  double argdi, ascle;
  int kdflg;
  double phidi, argdr;
  int ipard;
  double csgni, phidr, cspni, asumi[2], bsumi[2];
  double cspnr, asumr[2], bsumr[2];
  double zeta1i[2], zeta2i[2], zet1di, zet2di, zeta1r[2], zeta2r[2], zet1dr,
    zet2dr;
  int ib, ic;
  double fn;
  int il, kk, in, nw;
  double asumdi, bsumdi, yy, asumdr, bsumdr, c1i, c2i, c2m, c1r, c2r, s1i,
    s2i, rs1, s1r, s2r, aii, ang, asc, car, cki, fnf;
  int nai;
  double air;
  int ifn;
  double csi, ckr;
  int iuf;
  double cyi[2], fmr, sar, csr, sgn, zbi;
  int inu;
  double bry[3], cyr[2], pti, sti, zbr, zni, rzi, ptr, zri, str, znr, rzr,
    zrr;

  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  /* 
   */
  kdflg = 1;
  *nz = 0;
  /*----------------------------------------------------------------------- 
   *    EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN 
   *    THE UNDERFLOW LIMIT 
   *----------------------------------------------------------------------- 
   */
  cscl = 1. / *tol;
  crsc = *tol;
  cssr[0] = cscl;
  cssr[1] = coner;
  cssr[2] = crsc;
  csrr[0] = crsc;
  csrr[1] = coner;
  csrr[2] = cscl;
  bry[0] = pnl_d1mach (1) * 1e3 / *tol;
  bry[1] = 1. / bry[0];
  bry[2] = pnl_d1mach (2);
  zrr = *zr;
  zri = *zi;
  if (*zr >= 0.)
    {
      goto L10;
    }
  zrr = -(*zr);
  zri = -(*zi);
 L10:
  yy = zri;
  znr = zri;
  zni = -zrr;
  zbr = zrr;
  zbi = zri;
  inu = (int) (*fnu);
  fnf = *fnu - (double) inu;
  ang = -hpi * fnf;
  car = cos (ang);
  sar = sin (ang);
  c2r = hpi * sar;
  c2i = -hpi * car;
  kk = inu % 4 + 1;
  str = c2r * cipr[kk - 1] - c2i * cipi[kk - 1];
  sti = c2r * cipi[kk - 1] + c2i * cipr[kk - 1];
  csr = cr1r * str - cr1i * sti;
  csi = cr1r * sti + cr1i * str;
  if (yy > 0.)
    {
      goto L20;
    }
  znr = -znr;
  zbi = -zbi;
 L20:
  /*----------------------------------------------------------------------- 
   *    K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST 
   *    QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY 
   *    CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS 
   *----------------------------------------------------------------------- 
   */
  j = 2;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /*----------------------------------------------------------------------- 
       *    J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J 
       *----------------------------------------------------------------------- 
       */
      j = 3 - j;
      fn = *fnu + (double) (i__ - 1);
      amos_zunhj (&znr, &zni, &fn, &c__0, tol, &phir[j - 1], &phii[j - 1],
		  &argr[j - 1], &argi[j - 1], &zeta1r[j - 1], &zeta1i[j - 1],
		  &zeta2r[j - 1], &zeta2i[j - 1], &asumr[j - 1],
		  &asumi[j - 1], &bsumr[j - 1], &bsumi[j - 1]);
      if (*kode == 1)
	{
	  goto L30;
	}
      str = zbr + zeta2r[j - 1];
      sti = zbi + zeta2i[j - 1];
      rast = fn / amos_azabs (&str, &sti);
      str = str * rast * rast;
      sti = -sti * rast * rast;
      s1r = zeta1r[j - 1] - str;
      s1i = zeta1i[j - 1] - sti;
      goto L40;
    L30:
      s1r = zeta1r[j - 1] - zeta2r[j - 1];
      s1i = zeta1i[j - 1] - zeta2i[j - 1];
    L40:
      /*----------------------------------------------------------------------- 
       *    TEST FOR UNDERFLOW AND OVERFLOW 
       *----------------------------------------------------------------------- 
       */
      rs1 = s1r;
      if (fabs (rs1) > *elim)
	{
	  goto L70;
	}
      if (kdflg == 1)
	{
	  kflag = 2;
	}
      if (fabs (rs1) < *alim)
	{
	  goto L50;
	}
      /*----------------------------------------------------------------------- 
       *    REFINE  TEST AND SCALE 
       *----------------------------------------------------------------------- 
       */
      aphi = amos_azabs (&phir[j - 1], &phii[j - 1]);
      aarg = amos_azabs (&argr[j - 1], &argi[j - 1]);
      rs1 = rs1 + log (aphi) - log (aarg) * .25 - aic;
      if (fabs (rs1) > *elim)
	{
	  goto L70;
	}
      if (kdflg == 1)
	{
	  kflag = 1;
	}
      if (rs1 < 0.)
	{
	  goto L50;
	}
      if (kdflg == 1)
	{
	  kflag = 3;
	}
    L50:
      /*----------------------------------------------------------------------- 
       *    SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR 
       *    EXPONENT EXTREMES 
       *----------------------------------------------------------------------- 
       */
      c2r = argr[j - 1] * cr2r - argi[j - 1] * cr2i;
      c2i = argr[j - 1] * cr2i + argi[j - 1] * cr2r;
      pnl_zairy (&c2r, &c2i, &c__0, &c__2, &air, &aii, &nai, &idum);
      pnl_zairy (&c2r, &c2i, &c__1, &c__2, &dair, &daii, &ndai, &idum);
      str = dair * bsumr[j - 1] - daii * bsumi[j - 1];
      sti = dair * bsumi[j - 1] + daii * bsumr[j - 1];
      ptr = str * cr2r - sti * cr2i;
      pti = str * cr2i + sti * cr2r;
      str = ptr + (air * asumr[j - 1] - aii * asumi[j - 1]);
      sti = pti + (air * asumi[j - 1] + aii * asumr[j - 1]);
      ptr = str * phir[j - 1] - sti * phii[j - 1];
      pti = str * phii[j - 1] + sti * phir[j - 1];
      s2r = ptr * csr - pti * csi;
      s2i = ptr * csi + pti * csr;
      str = exp (s1r) * cssr[kflag - 1];
      s1r = str * cos (s1i);
      s1i = str * sin (s1i);
      str = s2r * s1r - s2i * s1i;
      s2i = s1r * s2i + s2r * s1i;
      s2r = str;
      if (kflag != 1)
	{
	  goto L60;
	}
      amos_zuchk (&s2r, &s2i, &nw, bry, tol);
      if (nw != 0)
	{
	  goto L70;
	}
    L60:
      if (yy <= 0.)
	{
	  s2i = -s2i;
	}
      cyr[kdflg - 1] = s2r;
      cyi[kdflg - 1] = s2i;
      yr[i__] = s2r * csrr[kflag - 1];
      yi[i__] = s2i * csrr[kflag - 1];
      str = csi;
      csi = -csr;
      csr = str;
      if (kdflg == 2)
	{
	  goto L85;
	}
      kdflg = 2;
      goto L80;
    L70:
      if (rs1 > 0.)
	{
	  goto L320;
	}
      /*----------------------------------------------------------------------- 
       *    FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW 
       *----------------------------------------------------------------------- 
       */
      if (*zr < 0.)
	{
	  goto L320;
	}
      kdflg = 1;
      yr[i__] = zeror;
      yi[i__] = zeroi;
      ++(*nz);
      str = csi;
      csi = -csr;
      csr = str;
      if (i__ == 1)
	{
	  goto L80;
	}
      if (yr[i__ - 1] == zeror && yi[i__ - 1] == zeroi)
	{
	  goto L80;
	}
      yr[i__ - 1] = zeror;
      yi[i__ - 1] = zeroi;
      ++(*nz);
    L80:
      ;
    }
  i__ = *n;
 L85:
  razr = 1. / amos_azabs (&zrr, &zri);
  str = zrr * razr;
  sti = -zri * razr;
  rzr = (str + str) * razr;
  rzi = (sti + sti) * razr;
  ckr = fn * rzr;
  cki = fn * rzi;
  ib = i__ + 1;
  if (*n < ib)
    {
      goto L180;
    }
  /*----------------------------------------------------------------------- 
   *    TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO 
   *    ON UNDERFLOW. 
   *----------------------------------------------------------------------- 
   */
  fn = *fnu + (double) (*n - 1);
  ipard = 1;
  if (*mr != 0)
    {
      ipard = 0;
    }
  amos_zunhj (&znr, &zni, &fn, &ipard, tol, &phidr, &phidi, &argdr, &argdi,
	      &zet1dr, &zet1di, &zet2dr, &zet2di, &asumdr, &asumdi, &bsumdr,
	      &bsumdi);
  if (*kode == 1)
    {
      goto L90;
    }
  str = zbr + zet2dr;
  sti = zbi + zet2di;
  rast = fn / amos_azabs (&str, &sti);
  str = str * rast * rast;
  sti = -sti * rast * rast;
  s1r = zet1dr - str;
  s1i = zet1di - sti;
  goto L100;
 L90:
  s1r = zet1dr - zet2dr;
  s1i = zet1di - zet2di;
 L100:
  rs1 = s1r;
  if (fabs (rs1) > *elim)
    {
      goto L105;
    }
  if (fabs (rs1) < *alim)
    {
      goto L120;
    }
  /*---------------------------------------------------------------------------- 
   *    REFINE ESTIMATE AND TEST 
   *------------------------------------------------------------------------- 
   */
  aphi = amos_azabs (&phidr, &phidi);
  rs1 += log (aphi);
  if (fabs (rs1) < *elim)
    {
      goto L120;
    }
 L105:
  if (rs1 > 0.)
    {
      goto L320;
    }
  /*----------------------------------------------------------------------- 
   *    FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW 
   *----------------------------------------------------------------------- 
   */
  if (*zr < 0.)
    {
      goto L320;
    }
  *nz = *n;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      yr[i__] = zeror;
      yi[i__] = zeroi;
      /* L106: */
    }
  return 0;
 L120:
  s1r = cyr[0];
  s1i = cyi[0];
  s2r = cyr[1];
  s2i = cyi[1];
  c1r = csrr[kflag - 1];
  ascle = bry[kflag - 1];
  i__1 = *n;
  for (i__ = ib; i__ <= i__1; ++i__)
    {
      c2r = s2r;
      c2i = s2i;
      s2r = ckr * c2r - cki * c2i + s1r;
      s2i = ckr * c2i + cki * c2r + s1i;
      s1r = c2r;
      s1i = c2i;
      ckr += rzr;
      cki += rzi;
      c2r = s2r * c1r;
      c2i = s2i * c1r;
      yr[i__] = c2r;
      yi[i__] = c2i;
      if (kflag >= 3)
	{
	  goto L130;
	}
      str = fabs (c2r);
      sti = fabs (c2i);
      c2m = MAX (str, sti);
      if (c2m <= ascle)
	{
	  goto L130;
	}
      ++kflag;
      ascle = bry[kflag - 1];
      s1r *= c1r;
      s1i *= c1r;
      s2r = c2r;
      s2i = c2i;
      s1r *= cssr[kflag - 1];
      s1i *= cssr[kflag - 1];
      s2r *= cssr[kflag - 1];
      s2i *= cssr[kflag - 1];
      c1r = csrr[kflag - 1];
    L130:
      ;
    }
 L180:
  if (*mr == 0)
    {
      return 0;
    }
  /*----------------------------------------------------------------------- 
   *    ANALYTIC CONTINUATION FOR RE(Z).LT.0.0D0 
   *----------------------------------------------------------------------- 
   */
  *nz = 0;
  fmr = (double) (*mr);
  sgn = -D_SIGN  (pi, fmr);
  /*----------------------------------------------------------------------- 
   *    CSPN AND CSGN ARE COEFF OF K AND I FUNCIONS RESP. 
   *----------------------------------------------------------------------- 
   */
  csgni = sgn;
  if (yy <= 0.)
    {
      csgni = -csgni;
    }
  ifn = inu + *n - 1;
  ang = fnf * sgn;
  cspnr = cos (ang);
  cspni = sin (ang);
  if (ifn % 2 == 0)
    {
      goto L190;
    }
  cspnr = -cspnr;
  cspni = -cspni;
 L190:
  /*----------------------------------------------------------------------- 
   *    CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION. I(FNU,Z) IS 
   *    COMPUTED FROM EXP(I*FNU*HPI)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST 
   *    QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY 
   *    CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS 
   *----------------------------------------------------------------------- 
   */
  csr = sar * csgni;
  csi = car * csgni;
  in = ifn % 4 + 1;
  c2r = cipr[in - 1];
  c2i = cipi[in - 1];
  str = csr * c2r + csi * c2i;
  csi = -csr * c2i + csi * c2r;
  csr = str;
  asc = bry[0];
  iuf = 0;
  kk = *n;
  kdflg = 1;
  --ib;
  ic = ib - 1;
  i__1 = *n;
  for (k = 1; k <= i__1; ++k)
    {
      fn = *fnu + (double) (kk - 1);
      /*----------------------------------------------------------------------- 
       *    LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K 
       *    FUNCTION ABOVE 
       *----------------------------------------------------------------------- 
       */
      if (*n > 2)
	{
	  goto L175;
	}
    L172:
      phidr = phir[j - 1];
      phidi = phii[j - 1];
      argdr = argr[j - 1];
      argdi = argi[j - 1];
      zet1dr = zeta1r[j - 1];
      zet1di = zeta1i[j - 1];
      zet2dr = zeta2r[j - 1];
      zet2di = zeta2i[j - 1];
      asumdr = asumr[j - 1];
      asumdi = asumi[j - 1];
      bsumdr = bsumr[j - 1];
      bsumdi = bsumi[j - 1];
      j = 3 - j;
      goto L210;
    L175:
      if (kk == *n && ib < *n)
	{
	  goto L210;
	}
      if (kk == ib || kk == ic)
	{
	  goto L172;
	}
      amos_zunhj (&znr, &zni, &fn, &c__0, tol, &phidr, &phidi, &argdr, &argdi,
		  &zet1dr, &zet1di, &zet2dr, &zet2di, &asumdr, &asumdi,
		  &bsumdr, &bsumdi);
    L210:
      if (*kode == 1)
	{
	  goto L220;
	}
      str = zbr + zet2dr;
      sti = zbi + zet2di;
      rast = fn / amos_azabs (&str, &sti);
      str = str * rast * rast;
      sti = -sti * rast * rast;
      s1r = -zet1dr + str;
      s1i = -zet1di + sti;
      goto L230;
    L220:
      s1r = -zet1dr + zet2dr;
      s1i = -zet1di + zet2di;
    L230:
      /*----------------------------------------------------------------------- 
       *    TEST FOR UNDERFLOW AND OVERFLOW 
       *----------------------------------------------------------------------- 
       */
      rs1 = s1r;
      if (fabs (rs1) > *elim)
	{
	  goto L280;
	}
      if (kdflg == 1)
	{
	  iflag = 2;
	}
      if (fabs (rs1) < *alim)
	{
	  goto L240;
	}
      /*----------------------------------------------------------------------- 
       *    REFINE  TEST AND SCALE 
       *----------------------------------------------------------------------- 
       */
      aphi = amos_azabs (&phidr, &phidi);
      aarg = amos_azabs (&argdr, &argdi);
      rs1 = rs1 + log (aphi) - log (aarg) * .25 - aic;
      if (fabs (rs1) > *elim)
	{
	  goto L280;
	}
      if (kdflg == 1)
	{
	  iflag = 1;
	}
      if (rs1 < 0.)
	{
	  goto L240;
	}
      if (kdflg == 1)
	{
	  iflag = 3;
	}
    L240:
      pnl_zairy (&argdr, &argdi, &c__0, &c__2, &air, &aii, &nai, &idum);
      pnl_zairy (&argdr, &argdi, &c__1, &c__2, &dair, &daii, &ndai, &idum);
      str = dair * bsumdr - daii * bsumdi;
      sti = dair * bsumdi + daii * bsumdr;
      str += air * asumdr - aii * asumdi;
      sti += air * asumdi + aii * asumdr;
      ptr = str * phidr - sti * phidi;
      pti = str * phidi + sti * phidr;
      s2r = ptr * csr - pti * csi;
      s2i = ptr * csi + pti * csr;
      str = exp (s1r) * cssr[iflag - 1];
      s1r = str * cos (s1i);
      s1i = str * sin (s1i);
      str = s2r * s1r - s2i * s1i;
      s2i = s2r * s1i + s2i * s1r;
      s2r = str;
      if (iflag != 1)
	{
	  goto L250;
	}
      amos_zuchk (&s2r, &s2i, &nw, bry, tol);
      if (nw == 0)
	{
	  goto L250;
	}
      s2r = zeror;
      s2i = zeroi;
    L250:
      if (yy <= 0.)
	{
	  s2i = -s2i;
	}
      cyr[kdflg - 1] = s2r;
      cyi[kdflg - 1] = s2i;
      c2r = s2r;
      c2i = s2i;
      s2r *= csrr[iflag - 1];
      s2i *= csrr[iflag - 1];
      /*----------------------------------------------------------------------- 
       *    ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N 
       *----------------------------------------------------------------------- 
       */
      s1r = yr[kk];
      s1i = yi[kk];
      if (*kode == 1)
	{
	  goto L270;
	}
      amos_zs1s2 (&zrr, &zri, &s1r, &s1i, &s2r, &s2i, &nw, &asc, alim, &iuf);
      *nz += nw;
    L270:
      yr[kk] = s1r * cspnr - s1i * cspni + s2r;
      yi[kk] = s1r * cspni + s1i * cspnr + s2i;
      --kk;
      cspnr = -cspnr;
      cspni = -cspni;
      str = csi;
      csi = -csr;
      csr = str;
      if (c2r != 0. || c2i != 0.)
	{
	  goto L255;
	}
      kdflg = 1;
      goto L290;
    L255:
      if (kdflg == 2)
	{
	  goto L295;
	}
      kdflg = 2;
      goto L290;
    L280:
      if (rs1 > 0.)
	{
	  goto L320;
	}
      s2r = zeror;
      s2i = zeroi;
      goto L250;
    L290:
      ;
    }
  k = *n;
 L295:
  il = *n - k;
  if (il == 0)
    {
      return 0;
    }
  /*----------------------------------------------------------------------- 
   *    RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE 
   *    K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP 
   *    INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES. 
   *----------------------------------------------------------------------- 
   */
  s1r = cyr[0];
  s1i = cyi[0];
  s2r = cyr[1];
  s2i = cyi[1];
  csr = csrr[iflag - 1];
  ascle = bry[iflag - 1];
  fn = (double) (inu + il);
  i__1 = il;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      c2r = s2r;
      c2i = s2i;
      s2r = s1r + (fn + fnf) * (rzr * c2r - rzi * c2i);
      s2i = s1i + (fn + fnf) * (rzr * c2i + rzi * c2r);
      s1r = c2r;
      s1i = c2i;
      fn += -1.;
      c2r = s2r * csr;
      c2i = s2i * csr;
      ckr = c2r;
      cki = c2i;
      c1r = yr[kk];
      c1i = yi[kk];
      if (*kode == 1)
	{
	  goto L300;
	}
      amos_zs1s2 (&zrr, &zri, &c1r, &c1i, &c2r, &c2i, &nw, &asc, alim, &iuf);
      *nz += nw;
    L300:
      yr[kk] = c1r * cspnr - c1i * cspni + c2r;
      yi[kk] = c1r * cspni + c1i * cspnr + c2i;
      --kk;
      cspnr = -cspnr;
      cspni = -cspni;
      if (iflag >= 3)
	{
	  goto L310;
	}
      c2r = fabs (ckr);
      c2i = fabs (cki);
      c2m = MAX (c2r, c2i);
      if (c2m <= ascle)
	{
	  goto L310;
	}
      ++iflag;
      ascle = bry[iflag - 1];
      s1r *= csr;
      s1i *= csr;
      s2r = ckr;
      s2i = cki;
      s1r *= cssr[iflag - 1];
      s1i *= cssr[iflag - 1];
      s2r *= cssr[iflag - 1];
      s2i *= cssr[iflag - 1];
      csr = csrr[iflag - 1];
    L310:
      ;
    }
  return 0;
 L320:
  *nz = -1;
  return 0;
}

