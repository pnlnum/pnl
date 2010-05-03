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



/* Table of constant values */

static const int c__1 = 1;
static const int c__0 = 0;

/*
 *    zuni1 computes i(fnu,z)  by means of the uniform asymptotic 
 *    expansion for i(fnu,z) in -pi/3.le.arg z.le.pi/3. 
 * 
 *    fnul is the smallest order permitted for the asymptotic 
 *    expansion. nlast=0 means all of the y values were set. 
 *    nlast.ne.0 is the number left to be computed by another 
 *    formula for orders fnu to fnu+nlast-1 because fnu+nlast-1.lt.fnul. 
 *    y(i)=czero for i=nlast+1,n 
 *
 */

int amos_zuni1 (double *zr, double *zi, double *fnu, int *kode,const int *n,
		double *yr, double *yi, int *nz, int *nlast, double *fnul,
		double *tol, double *elim, double *alim)
{
  /* Initialized data */

  static const double zeror = 0.;
  static const double zeroi = 0.;
  static const double coner = 1.;

  /* System generated locals */
  int i__1;

  /* Local variables */
  double aphi, cscl, phii, crsc, phir;
  int init;
  double csrr[3], cssr[3], rast, sumi, sumr;
  int i__, k, m, iflag=1;
  double ascle;
  double cwrki[16];
  double cwrkr[16];
  double zeta1i, zeta2i, zeta1r, zeta2r;
  int nd;
  double fn;
  int nn, nw;
  double c2i, c2m, c1r, c2r, s1i, s2i, rs1, s1r, s2r, cyi[2];
  int nuf;
  double bry[3], cyr[2], sti, rzi, str, rzr;

  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  /* 
   */
  *nz = 0;
  nd = *n;
  *nlast = 0;
  /*----------------------------------------------------------------------- 
   *    COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG- 
   *    NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE, 
   *    EXP(ALIM)=EXP(ELIM)*TOL 
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
  /*----------------------------------------------------------------------- 
   *    CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER 
   *----------------------------------------------------------------------- 
   */
  fn = MAX (*fnu, 1.);
  init = 0;
  amos_zunik (zr, zi, &fn, &c__1, &c__1, tol, &init, &phir, &phii, &zeta1r,
	      &zeta1i, &zeta2r, &zeta2i, &sumr, &sumi, cwrkr, cwrki);
  if (*kode == 1)
    {
      goto L10;
    }
  str = *zr + zeta2r;
  sti = *zi + zeta2i;
  rast = fn / amos_azabs (&str, &sti);
  str = str * rast * rast;
  sti = -sti * rast * rast;
  s1r = -zeta1r + str;
  s1i = -zeta1i + sti;
  goto L20;
 L10:
  s1r = -zeta1r + zeta2r;
  s1i = -zeta1i + zeta2i;
 L20:
  rs1 = s1r;
  if (fabs (rs1) > *elim)
    {
      goto L130;
    }
 L30:
  nn = MIN (2, nd);
  i__1 = nn;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      fn = *fnu + (double) (nd - i__);
      init = 0;
      amos_zunik (zr, zi, &fn, &c__1, &c__0, tol, &init, &phir, &phii,
		  &zeta1r, &zeta1i, &zeta2r, &zeta2i, &sumr, &sumi, cwrkr,
		  cwrki);
      if (*kode == 1)
	{
	  goto L40;
	}
      str = *zr + zeta2r;
      sti = *zi + zeta2i;
      rast = fn / amos_azabs (&str, &sti);
      str = str * rast * rast;
      sti = -sti * rast * rast;
      s1r = -zeta1r + str;
      s1i = -zeta1i + sti + *zi;
      goto L50;
    L40:
      s1r = -zeta1r + zeta2r;
      s1i = -zeta1i + zeta2i;
    L50:
      /*----------------------------------------------------------------------- 
       *    TEST FOR UNDERFLOW AND OVERFLOW 
       *----------------------------------------------------------------------- 
       */
      rs1 = s1r;
      if (fabs (rs1) > *elim)
	{
	  goto L110;
	}
      if (i__ == 1)
	{
	  iflag = 2;
	}
      if (fabs (rs1) < *alim)
	{
	  goto L60;
	}
      /*----------------------------------------------------------------------- 
       *    REFINE  TEST AND SCALE 
       *----------------------------------------------------------------------- 
       */
      aphi = amos_azabs (&phir, &phii);
      rs1 += log (aphi);
      if (fabs (rs1) > *elim)
	{
	  goto L110;
	}
      if (i__ == 1)
	{
	  iflag = 1;
	}
      if (rs1 < 0.)
	{
	  goto L60;
	}
      if (i__ == 1)
	{
	  iflag = 3;
	}
    L60:
      /*----------------------------------------------------------------------- 
       *    SCALE S1 IF CABS(S1).LT.ASCLE 
       *----------------------------------------------------------------------- 
       */
      s2r = phir * sumr - phii * sumi;
      s2i = phir * sumi + phii * sumr;
      str = exp (s1r) * cssr[iflag - 1];
      s1r = str * cos (s1i);
      s1i = str * sin (s1i);
      str = s2r * s1r - s2i * s1i;
      s2i = s2r * s1i + s2i * s1r;
      s2r = str;
      if (iflag != 1)
	{
	  goto L70;
	}
      amos_zuchk (&s2r, &s2i, &nw, bry, tol);
      if (nw != 0)
	{
	  goto L110;
	}
    L70:
      cyr[i__ - 1] = s2r;
      cyi[i__ - 1] = s2i;
      m = nd - i__ + 1;
      yr[m] = s2r * csrr[iflag - 1];
      yi[m] = s2i * csrr[iflag - 1];
      /* L80: */
    }
  if (nd <= 2)
    {
      goto L100;
    }
  rast = 1. / amos_azabs (zr, zi);
  str = *zr * rast;
  sti = -(*zi) * rast;
  rzr = (str + str) * rast;
  rzi = (sti + sti) * rast;
  bry[1] = 1. / bry[0];
  bry[2] = pnl_d1mach (2);
  s1r = cyr[0];
  s1i = cyi[0];
  s2r = cyr[1];
  s2i = cyi[1];
  c1r = csrr[iflag - 1];
  ascle = bry[iflag - 1];
  k = nd - 2;
  fn = (double) k;
  i__1 = nd;
  for (i__ = 3; i__ <= i__1; ++i__)
    {
      c2r = s2r;
      c2i = s2i;
      s2r = s1r + (*fnu + fn) * (rzr * c2r - rzi * c2i);
      s2i = s1i + (*fnu + fn) * (rzr * c2i + rzi * c2r);
      s1r = c2r;
      s1i = c2i;
      c2r = s2r * c1r;
      c2i = s2i * c1r;
      yr[k] = c2r;
      yi[k] = c2i;
      --k;
      fn += -1.;
      if (iflag >= 3)
	{
	  goto L90;
	}
      str = fabs (c2r);
      sti = fabs (c2i);
      c2m = MAX (str, sti);
      if (c2m <= ascle)
	{
	  goto L90;
	}
      ++iflag;
      ascle = bry[iflag - 1];
      s1r *= c1r;
      s1i *= c1r;
      s2r = c2r;
      s2i = c2i;
      s1r *= cssr[iflag - 1];
      s1i *= cssr[iflag - 1];
      s2r *= cssr[iflag - 1];
      s2i *= cssr[iflag - 1];
      c1r = csrr[iflag - 1];
    L90:
      ;
    }
 L100:
  return 0;
  /*----------------------------------------------------------------------- 
   *    SET UNDERFLOW AND UPDATE PARAMETERS 
   *----------------------------------------------------------------------- 
   */
 L110:
  if (rs1 > 0.)
    {
      goto L120;
    }
  yr[nd] = zeror;
  yi[nd] = zeroi;
  ++(*nz);
  --nd;
  if (nd == 0)
    {
      goto L100;
    }
  amos_zuoik (zr, zi, fnu, kode, &c__1, &nd, &yr[1], &yi[1], &nuf, tol, elim, alim);
  if (nuf < 0)
    {
      goto L120;
    }
  nd -= nuf;
  *nz += nuf;
  if (nd == 0)
    {
      goto L100;
    }
  fn = *fnu + (double) (nd - 1);
  if (fn >= *fnul)
    {
      goto L30;
    }
  *nlast = nd;
  return 0;
 L120:
  *nz = -1;
  return 0;
 L130:
  if (rs1 > 0.)
    {
      goto L120;
    }
  *nz = *n;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      yr[i__] = zeror;
      yi[i__] = zeroi;
      /* L140: */
    }
  return 0;
}				/* zuni1_ */
