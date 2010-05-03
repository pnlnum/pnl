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
 *    zmlri computes the i bessel function for re(z).ge.0.0 by the 
 *    miller algorithm normalized by a neumann series. 
 *
 */


/* Subroutine */ int
amos_zmlri (double *zr, double *zi, double *fnu,const int *kode,const int *n,
	    double *yr, double *yi, int *nz, double *tol)
{
  /* Initialized data */

  static const double zeror = 0.;
  static const double zeroi = 0.;
  static const double coner = 1.;
  static const double conei = 0.;

  /* System generated locals */
  int i__1, i__2;
  double d__1, d__2, d__3;

  /* Local variables */
  double flam, fkap, scle, tfnf;
  int idum, ifnu;
  double sumi, sumr;
  int i__, k, m;
  int itime;
  double ak, bk, ap, at;
  int kk, km;
  double az;
  double cnormi, cnormr, p1i, p2i, p1r, p2r, ack, cki, fnf, fkk, ckr;
  int iaz;
  double rho;
  int inu;
  double pti, raz, sti, rzi, ptr, str, tst, rzr, rho2;

  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  scle = pnl_d1mach (1) / *tol;
  *nz = 0;
  az = amos_azabs (zr, zi);
  iaz = (int) az;
  ifnu = (int) (*fnu);
  inu = ifnu + *n - 1;
  at = (double) iaz + 1.;
  raz = 1. / az;
  str = *zr * raz;
  sti = -(*zi) * raz;
  ckr = str * at * raz;
  cki = sti * at * raz;
  rzr = (str + str) * raz;
  rzi = (sti + sti) * raz;
  p1r = zeror;
  p1i = zeroi;
  p2r = coner;
  p2i = conei;
  ack = (at + 1.) * raz;
  rho = ack + sqrt (ack * ack - 1.);
  rho2 = rho * rho;
  tst = (rho2 + rho2) / ((rho2 - 1.) * (rho - 1.));
  tst /= *tol;
  /*----------------------------------------------------------------------- 
   *    COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES 
   *----------------------------------------------------------------------- 
   */
  ak = at;
  for (i__ = 1; i__ <= 80; ++i__)
    {
      ptr = p2r;
      pti = p2i;
      p2r = p1r - (ckr * ptr - cki * pti);
      p2i = p1i - (cki * ptr + ckr * pti);
      p1r = ptr;
      p1i = pti;
      ckr += rzr;
      cki += rzi;
      ap = amos_azabs (&p2r, &p2i);
      if (ap > tst * ak * ak)
	{
	  goto L20;
	}
      ak += 1.;
      /* L10: */
    }
  goto L110;
 L20:
  ++i__;
  k = 0;
  if (inu < iaz)
    {
      goto L40;
    }
  /*----------------------------------------------------------------------- 
   *    COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS 
   *----------------------------------------------------------------------- 
   */
  p1r = zeror;
  p1i = zeroi;
  p2r = coner;
  p2i = conei;
  at = (double) inu + 1.;
  str = *zr * raz;
  sti = -(*zi) * raz;
  ckr = str * at * raz;
  cki = sti * at * raz;
  ack = at * raz;
  tst = sqrt (ack / *tol);
  itime = 1;
  for (k = 1; k <= 80; ++k)
    {
      ptr = p2r;
      pti = p2i;
      p2r = p1r - (ckr * ptr - cki * pti);
      p2i = p1i - (ckr * pti + cki * ptr);
      p1r = ptr;
      p1i = pti;
      ckr += rzr;
      cki += rzi;
      ap = amos_azabs (&p2r, &p2i);
      if (ap < tst)
	{
	  goto L30;
	}
      if (itime == 2)
	{
	  goto L40;
	}
      ack = amos_azabs (&ckr, &cki);
      flam = ack + sqrt (ack * ack - 1.);
      fkap = ap / amos_azabs (&p1r, &p1i);
      rho = MIN (flam, fkap);
      tst *= sqrt (rho / (rho * rho - 1.));
      itime = 2;
    L30:
      ;
    }
  goto L110;
 L40:
  /*----------------------------------------------------------------------- 
   *    BACKWARD RECURRENCE AND SUM NORMALIZING RELATION 
   *----------------------------------------------------------------------- 
   */
  ++k;
  /*Computing MAX 
   */
  i__1 = i__ + iaz, i__2 = k + inu;
  kk = MAX (i__1, i__2);
  fkk = (double) kk;
  p1r = zeror;
  p1i = zeroi;
  /*----------------------------------------------------------------------- 
   *    SCALE P2 AND SUM BY SCLE 
   *----------------------------------------------------------------------- 
   */
  p2r = scle;
  p2i = zeroi;
  fnf = *fnu - (double) ifnu;
  tfnf = fnf + fnf;
  d__1 = fkk + tfnf + 1.;
  d__2 = fkk + 1.;
  d__3 = tfnf + 1.;
  bk =
    amos_dgamln (&d__1, &idum) - amos_dgamln (&d__2,
					      &idum) - amos_dgamln (&d__3,
								    &idum);
  bk = exp (bk);
  sumr = zeror;
  sumi = zeroi;
  km = kk - inu;
  i__1 = km;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      ptr = p2r;
      pti = p2i;
      p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
      p2i = p1i + (fkk + fnf) * (rzi * ptr + rzr * pti);
      p1r = ptr;
      p1i = pti;
      ak = 1. - tfnf / (fkk + tfnf);
      ack = bk * ak;
      sumr += (ack + bk) * p1r;
      sumi += (ack + bk) * p1i;
      bk = ack;
      fkk += -1.;
      /* L50: */
    }
  yr[*n] = p2r;
  yi[*n] = p2i;
  if (*n == 1)
    {
      goto L70;
    }
  i__1 = *n;
  for (i__ = 2; i__ <= i__1; ++i__)
    {
      ptr = p2r;
      pti = p2i;
      p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
      p2i = p1i + (fkk + fnf) * (rzi * ptr + rzr * pti);
      p1r = ptr;
      p1i = pti;
      ak = 1. - tfnf / (fkk + tfnf);
      ack = bk * ak;
      sumr += (ack + bk) * p1r;
      sumi += (ack + bk) * p1i;
      bk = ack;
      fkk += -1.;
      m = *n - i__ + 1;
      yr[m] = p2r;
      yi[m] = p2i;
      /* L60: */
    }
 L70:
  if (ifnu <= 0)
    {
      goto L90;
    }
  i__1 = ifnu;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      ptr = p2r;
      pti = p2i;
      p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
      p2i = p1i + (fkk + fnf) * (rzr * pti + rzi * ptr);
      p1r = ptr;
      p1i = pti;
      ak = 1. - tfnf / (fkk + tfnf);
      ack = bk * ak;
      sumr += (ack + bk) * p1r;
      sumi += (ack + bk) * p1i;
      bk = ack;
      fkk += -1.;
      /* L80: */
    }
 L90:
  ptr = *zr;
  pti = *zi;
  if (*kode == 2)
    {
      ptr = zeror;
    }
  amos_azlog (&rzr, &rzi, &str, &sti, &idum);
  p1r = -fnf * str + ptr;
  p1i = -fnf * sti + pti;
  d__1 = fnf + 1.;
  ap = amos_dgamln (&d__1, &idum);
  ptr = p1r - ap;
  pti = p1i;
  /*----------------------------------------------------------------------- 
   *    THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW 
   *    IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES 
   *----------------------------------------------------------------------- 
   */
  p2r += sumr;
  p2i += sumi;
  ap = amos_azabs (&p2r, &p2i);
  p1r = 1. / ap;
  amos_azexp (&ptr, &pti, &str, &sti);
  ckr = str * p1r;
  cki = sti * p1r;
  ptr = p2r * p1r;
  pti = -p2i * p1r;
  amos_zmlt (&ckr, &cki, &ptr, &pti, &cnormr, &cnormi);
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      str = yr[i__] * cnormr - yi[i__] * cnormi;
      yi[i__] = yr[i__] * cnormi + yi[i__] * cnormr;
      yr[i__] = str;
      /* L100: */
    }
  return 0;
 L110:
  *nz = -2;
  return 0;
}				/* zmlri_ */
