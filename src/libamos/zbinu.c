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
static const int c__2 = 2;

/*
 * 
 *    zbinu computes the i function in the right half z plane 
 *
 */

int amos_zbinu (double *zr, double *zi, double *fnu, int *kode,const int *n,
		double *cyr, double *cyi, int *nz, double *rl, double *fnul,
		double *tol, double *elim, double *alim)
{
  /* Initialized data */

  static const double zeror = 0.;
  static const double zeroi = 0.;

  /* System generated locals */
  int i__1;

  /* Local variables */
  double dfnu;
  int i__;
  int nlast;
  double az;
  int nn, nw;
  double cwi[2], cwr[2];
  int nui, inw;

  /* Parameter adjustments */
  --cyi;
  --cyr;

  /* Function Body */
  /* 
   */
  *nz = 0;
  az = amos_azabs (zr, zi);
  nn = *n;
  dfnu = *fnu + (double) (*n - 1);
  if (az <= 2.)
    {
      goto L10;
    }
  if (az * az * .25 > dfnu + 1.)
    {
      goto L20;
    }
 L10:
  /*----------------------------------------------------------------------- 
   *    POWER SERIES 
   *----------------------------------------------------------------------- 
   */
  amos_zseri (zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, tol, elim, alim);
  inw = abs (nw);
  *nz += inw;
  nn -= inw;
  if (nn == 0)
    {
      return 0;
    }
  if (nw >= 0)
    {
      goto L120;
    }
  dfnu = *fnu + (double) (nn - 1);
 L20:
  if (az < *rl)
    {
      goto L40;
    }
  if (dfnu <= 1.)
    {
      goto L30;
    }
  if (az + az < dfnu * dfnu)
    {
      goto L50;
    }
  /*----------------------------------------------------------------------- 
   *    ASYMPTOTIC EXPANSION FOR LARGE Z 
   *----------------------------------------------------------------------- 
   */
 L30:
  amos_zasyi (zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, rl, tol, elim,
	      alim);
  if (nw < 0)
    {
      goto L130;
    }
  goto L120;
 L40:
  if (dfnu <= 1.)
    {
      goto L70;
    }
 L50:
  /*----------------------------------------------------------------------- 
   *    OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM 
   *----------------------------------------------------------------------- 
   */
  amos_zuoik (zr, zi, fnu, kode, &c__1, &nn, &cyr[1], &cyi[1], &nw, tol, elim,
	      alim);
  if (nw < 0)
    {
      goto L130;
    }
  *nz += nw;
  nn -= nw;
  if (nn == 0)
    {
      return 0;
    }
  dfnu = *fnu + (double) (nn - 1);
  if (dfnu > *fnul)
    {
      goto L110;
    }
  if (az > *fnul)
    {
      goto L110;
    }
 L60:
  if (az > *rl)
    {
      goto L80;
    }
 L70:
  /*----------------------------------------------------------------------- 
   *    MILLER ALGORITHM NORMALIZED BY THE SERIES 
   *----------------------------------------------------------------------- 
   */
  amos_zmlri (zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, tol);
  if (nw < 0)
    {
      goto L130;
    }
  goto L120;
 L80:
  /*----------------------------------------------------------------------- 
   *    MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN 
   *----------------------------------------------------------------------- 
   *----------------------------------------------------------------------- 
   *    OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN 
   *----------------------------------------------------------------------- 
   */
  amos_zuoik (zr, zi, fnu, kode, &c__2, &c__2, cwr, cwi, &nw, tol, elim,
	      alim);
  if (nw >= 0)
    {
      goto L100;
    }
  *nz = nn;
  i__1 = nn;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      cyr[i__] = zeror;
      cyi[i__] = zeroi;
      /* L90: */
    }
  return 0;
 L100:
  if (nw > 0)
    {
      goto L130;
    }
  amos_zwrsk (zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, cwr, cwi, tol,
	      elim, alim);
  if (nw < 0)
    {
      goto L130;
    }
  goto L120;
 L110:
  /*----------------------------------------------------------------------- 
   *    INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD 
   *----------------------------------------------------------------------- 
   */
  nui = (int) (*fnul - dfnu) + 1;
  nui = MAX (nui, 0);
  amos_zbuni (zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, &nui, &nlast,
	      fnul, tol, elim, alim);
  if (nw < 0)
    {
      goto L130;
    }
  *nz += nw;
  if (nlast == 0)
    {
      goto L120;
    }
  nn = nlast;
  goto L60;
 L120:
  return 0;
 L130:
  *nz = -1;
  if (nw == -2)
    {
      *nz = -2;
    }
  return 0;
}				/* zbinu_ */
