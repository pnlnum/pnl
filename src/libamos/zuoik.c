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

/* 
 *    zuoik computes the leading terms of the uniform asymptotic 
 *    expansions for the i and k functions and compares them 
 *    (in logarithmic form) to alim and elim for over and underflow 
 *    where alim.lt.elim. if the magnitude, based on the leading 
 *    exponential, is less than alim or greater than -alim, then 
 *    the result is on scale. if not, then a refined test using other 
 *    multipliers (in logarithmic form) is made based on elim. here 
 *    exp(-elim)=smallest machine number*1.0e+3 and exp(-alim)= 
 *    exp(-elim)/tol 
 * 
 *    ikflg=1 means the i sequence is tested 
 *         =2 means the k sequence is tested 
 *    nuf = 0 means the last member of the sequence is on scale 
 *        =-1 means an overflow would occur 
 *    ikflg=1 and nuf.gt.0 means the last nuf y values were set to zero 
 *            the first n-nuf values must be set by another routine 
 *    ikflg=2 and nuf.eq.n means all y values were set to zero 
 *    ikflg=2 and 0.lt.nuf.lt.n not considered. y must be set by 
 *            another routine 
 */

int amos_zuoik (double *zr, double *zi, double *fnu,const int *kode,const int *ikflg,
		const int *n, double *yr, double *yi, int *nuf, double *tol,
		double *elim, double *alim)
{

  static const double zeror = 0.;
  static const double zeroi = 0.;
  static const double aic = 1.265512123484645396;

  /* System generated locals */
  int i__1;

  /* Local variables */
  double aarg=0.0, aphi, argi, phii, argr;
  int idum;
  double phir;
  int init;
  double sumi, sumr;
  int i__;
  double ascle;
  int iform;
  double asumi, bsumi, cwrki[16];
  double asumr, bsumr, cwrkr[16];
  double zeta1i, zeta2i, zeta1r, zeta2r, ax, ay;
  int nn, nw;
  double fnn, gnn, zbi, czi, gnu, zbr, czr, rcz, sti, zni, zri, str, znr, zrr;

  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  *nuf = 0;
  nn = *n;
  zrr = *zr;
  zri = *zi;
  if (*zr >= 0.)
    {
      goto L10;
    }
  zrr = -(*zr);
  zri = -(*zi);
 L10:
  zbr = zrr;
  zbi = zri;
  ax = fabs (*zr) * 1.7321;
  ay = fabs (*zi);
  iform = 1;
  if (ay > ax)
    {
      iform = 2;
    }
  gnu = MAX (*fnu, 1.);
  if (*ikflg == 1)
    {
      goto L20;
    }
  fnn = (double) nn;
  gnn = *fnu + fnn - 1.;
  gnu = MAX (gnn, fnn);
 L20:
  /*----------------------------------------------------------------------- 
   *    ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE 
   *    REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET 
   *    THE SIGN OF THE IMAGINARY PART CORRECT. 
   *----------------------------------------------------------------------- 
   */
  if (iform == 2)
    {
      goto L30;
    }
  init = 0;
  amos_zunik (&zrr, &zri, &gnu, ikflg, &c__1, tol, &init, &phir, &phii,
	      &zeta1r, &zeta1i, &zeta2r, &zeta2i, &sumr, &sumi, cwrkr, cwrki);
  czr = -zeta1r + zeta2r;
  czi = -zeta1i + zeta2i;
  goto L50;
 L30:
  znr = zri;
  zni = -zrr;
  if (*zi > 0.)
    {
      goto L40;
    }
  znr = -znr;
 L40:
  amos_zunhj (&znr, &zni, &gnu, &c__1, tol, &phir, &phii, &argr, &argi,
	      &zeta1r, &zeta1i, &zeta2r, &zeta2i, &asumr, &asumi, &bsumr,
	      &bsumi);
  czr = -zeta1r + zeta2r;
  czi = -zeta1i + zeta2i;
  aarg = amos_azabs (&argr, &argi);
 L50:
  if (*kode == 1)
    {
      goto L60;
    }
  czr -= zbr;
  czi -= zbi;
 L60:
  if (*ikflg == 1)
    {
      goto L70;
    }
  czr = -czr;
  czi = -czi;
 L70:
  aphi = amos_azabs (&phir, &phii);
  rcz = czr;
  /*----------------------------------------------------------------------- 
   *    OVERFLOW TEST 
   *----------------------------------------------------------------------- 
   */
  if (rcz > *elim)
    {
      goto L210;
    }
  if (rcz < *alim)
    {
      goto L80;
    }
  rcz += log (aphi);
  if (iform == 2)
    {
      rcz = rcz - log (aarg) * .25 - aic;
    }
  if (rcz > *elim)
    {
      goto L210;
    }
  goto L130;
 L80:
  /*----------------------------------------------------------------------- 
   *    UNDERFLOW TEST 
   *----------------------------------------------------------------------- 
   */
  if (rcz < -(*elim))
    {
      goto L90;
    }
  if (rcz > -(*alim))
    {
      goto L130;
    }
  rcz += log (aphi);
  if (iform == 2)
    {
      rcz = rcz - log (aarg) * .25 - aic;
    }
  if (rcz > -(*elim))
    {
      goto L110;
    }
 L90:
  i__1 = nn;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      yr[i__] = zeror;
      yi[i__] = zeroi;
      /* L100: */
    }
  *nuf = nn;
  return 0;
 L110:
  ascle = pnl_d1mach (1) * 1e3 / *tol;
  amos_azlog (&phir, &phii, &str, &sti, &idum);
  czr += str;
  czi += sti;
  if (iform == 1)
    {
      goto L120;
    }
  amos_azlog (&argr, &argi, &str, &sti, &idum);
  czr = czr - str * .25 - aic;
  czi -= sti * .25;
 L120:
  ax = exp (rcz) / *tol;
  ay = czi;
  czr = ax * cos (ay);
  czi = ax * sin (ay);
  amos_zuchk (&czr, &czi, &nw, &ascle, tol);
  if (nw != 0)
    {
      goto L90;
    }
 L130:
  if (*ikflg == 2)
    {
      return 0;
    }
  if (*n == 1)
    {
      return 0;
    }
  /*----------------------------------------------------------------------- 
   *    SET UNDERFLOWS ON I SEQUENCE 
   *----------------------------------------------------------------------- 
   */
 L140:
  gnu = *fnu + (double) (nn - 1);
  if (iform == 2)
    {
      goto L150;
    }
  init = 0;
  amos_zunik (&zrr, &zri, &gnu, ikflg, &c__1, tol, &init, &phir, &phii,
	      &zeta1r, &zeta1i, &zeta2r, &zeta2i, &sumr, &sumi, cwrkr, cwrki);
  czr = -zeta1r + zeta2r;
  czi = -zeta1i + zeta2i;
  goto L160;
 L150:
  amos_zunhj (&znr, &zni, &gnu, &c__1, tol, &phir, &phii, &argr, &argi,
	      &zeta1r, &zeta1i, &zeta2r, &zeta2i, &asumr, &asumi, &bsumr,
	      &bsumi);
  czr = -zeta1r + zeta2r;
  czi = -zeta1i + zeta2i;
  aarg = amos_azabs (&argr, &argi);
 L160:
  if (*kode == 1)
    {
      goto L170;
    }
  czr -= zbr;
  czi -= zbi;
 L170:
  aphi = amos_azabs (&phir, &phii);
  rcz = czr;
  if (rcz < -(*elim))
    {
      goto L180;
    }
  if (rcz > -(*alim))
    {
      return 0;
    }
  rcz += log (aphi);
  if (iform == 2)
    {
      rcz = rcz - log (aarg) * .25 - aic;
    }
  if (rcz > -(*elim))
    {
      goto L190;
    }
 L180:
  yr[nn] = zeror;
  yi[nn] = zeroi;
  --nn;
  ++(*nuf);
  if (nn == 0)
    {
      return 0;
    }
  goto L140;
 L190:
  ascle = pnl_d1mach (1) * 1e3 / *tol;
  amos_azlog (&phir, &phii, &str, &sti, &idum);
  czr += str;
  czi += sti;
  if (iform == 1)
    {
      goto L200;
    }
  amos_azlog (&argr, &argi, &str, &sti, &idum);
  czr = czr - str * .25 - aic;
  czi -= sti * .25;
 L200:
  ax = exp (rcz) / *tol;
  ay = czi;
  czr = ax * cos (ay);
  czi = ax * sin (ay);
  amos_zuchk (&czr, &czi, &nw, &ascle, tol);
  if (nw != 0)
    {
      goto L180;
    }
  return 0;
 L210:
  *nuf = -1;
  return 0;
}				/* zuoik_ */
