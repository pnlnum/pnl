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
 *    set k functions to zero on underflow, continue recurrence 
 *    on scaled functions until two members come on scale, then 
 *    return with min(nz+2,n) values scaled by 1/tol. 
 *
 */

/* Subroutine */ int
amos_zkscl (double *zrr, double *zri, double *fnu,const int *n, double *yr,
	    double *yi, int *nz, double *rzr, double *rzi, double *ascle,
	    double *tol, double *elim)
{
  /* Initialized data */

  static const double zeror = 0.;
  static const double zeroi = 0.;

  /* System generated locals */
  int i__1;

  /* Local variables */
  double alas;
  int idum, i__;
  double helim;
  double celmr;
  int ic;
  double as, fn;
  int kk, nn, nw;
  double s1i, s2i, s1r, s2r, acs, cki, elm, csi, ckr, cyi[2], zdi, csr,
    cyr[2], zdr, str;

  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  /* 
   */
  *nz = 0;
  ic = 0;
  nn = MIN (2, *n);
  i__1 = nn;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      s1r = yr[i__];
      s1i = yi[i__];
      cyr[i__ - 1] = s1r;
      cyi[i__ - 1] = s1i;
      as = amos_azabs (&s1r, &s1i);
      acs = -(*zrr) + log (as);
      ++(*nz);
      yr[i__] = zeror;
      yi[i__] = zeroi;
      if (acs < -(*elim))
	{
	  goto L10;
	}
      amos_azlog (&s1r, &s1i, &csr, &csi, &idum);
      csr -= *zrr;
      csi -= *zri;
      str = exp (csr) / *tol;
      csr = str * cos (csi);
      csi = str * sin (csi);
      amos_zuchk (&csr, &csi, &nw, ascle, tol);
      if (nw != 0)
	{
	  goto L10;
	}
      yr[i__] = csr;
      yi[i__] = csi;
      ic = i__;
      --(*nz);
    L10:
      ;
    }
  if (*n == 1)
    {
      return 0;
    }
  if (ic > 1)
    {
      goto L20;
    }
  yr[1] = zeror;
  yi[1] = zeroi;
  *nz = 2;
 L20:
  if (*n == 2)
    {
      return 0;
    }
  if (*nz == 0)
    {
      return 0;
    }
  fn = *fnu + 1.;
  ckr = fn * *rzr;
  cki = fn * *rzi;
  s1r = cyr[0];
  s1i = cyi[0];
  s2r = cyr[1];
  s2i = cyi[1];
  helim = *elim * .5;
  elm = exp (-(*elim));
  celmr = elm;
  zdr = *zrr;
  zdi = *zri;
  /* 
   *    FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF 
   *    S2 GETS LARGER THAN EXP(ELIM/2) 
   * 
   */
  i__1 = *n;
  for (i__ = 3; i__ <= i__1; ++i__)
    {
      kk = i__;
      csr = s2r;
      csi = s2i;
      s2r = ckr * csr - cki * csi + s1r;
      s2i = cki * csr + ckr * csi + s1i;
      s1r = csr;
      s1i = csi;
      ckr += *rzr;
      cki += *rzi;
      as = amos_azabs (&s2r, &s2i);
      alas = log (as);
      acs = -zdr + alas;
      ++(*nz);
      yr[i__] = zeror;
      yi[i__] = zeroi;
      if (acs < -(*elim))
	{
	  goto L25;
	}
      amos_azlog (&s2r, &s2i, &csr, &csi, &idum);
      csr -= zdr;
      csi -= zdi;
      str = exp (csr) / *tol;
      csr = str * cos (csi);
      csi = str * sin (csi);
      amos_zuchk (&csr, &csi, &nw, ascle, tol);
      if (nw != 0)
	{
	  goto L25;
	}
      yr[i__] = csr;
      yi[i__] = csi;
      --(*nz);
      if (ic == kk - 1)
	{
	  goto L40;
	}
      ic = kk;
      goto L30;
    L25:
      if (alas < helim)
	{
	  goto L30;
	}
      zdr -= *elim;
      s1r *= celmr;
      s1i *= celmr;
      s2r *= celmr;
      s2i *= celmr;
    L30:
      ;
    }
  *nz = *n;
  if (ic == *n)
    {
      *nz = *n - 1;
    }
  goto L45;
 L40:
  *nz = kk - 2;
 L45:
  i__1 = *nz;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      yr[i__] = zeror;
      yi[i__] = zeroi;
      /* L50: */
    }
  return 0;
}				/* zkscl_ */
