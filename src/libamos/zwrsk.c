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



static const int c__2 = 2;

/*
 *    zwrsk computes the i bessel function for re(z).ge.0.0 by 
 *    normalizing the i function ratios from zrati by the wronskian 
 * 
 *    i(fnu+i-1,z) by backward recurrence for ratios 
 *    y(i)=i(fnu+i,z)/i(fnu+i-1,z) from crati normalized by the 
 *    wronskian with k(fnu,z) and k(fnu+1,z) from cbknu. 
 */

int amos_zwrsk (double *zrr, double *zri, double *fnu, int *kode,const int *n,
		double *yr, double *yi, int *nz, double *cwr, double *cwi,
		double *tol, double *elim, double *alim)
{
  int i__1;
  double ract;
  int i__;
  double ascle;
  double csclr, cinui, cinur;
  int nw;
  double c1i, c2i, c1r, c2r, act, acw, cti, ctr, pti, sti, ptr, str;

  /* Parameter adjustments */
  --yi;
  --yr;
  --cwr;
  --cwi;

  /* Function Body */
  *nz = 0;
  amos_zbknu (zrr, zri, fnu, kode, &c__2, &cwr[1], &cwi[1], &nw, tol, elim, alim);
  if (nw != 0)
    {
      goto L50;
    }
  pnl_zrati (zrr, zri, fnu, n, &yr[1], &yi[1], tol);
  /*----------------------------------------------------------------------- 
   *    RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z), 
   *    R(FNU+J-1,Z)=Y(J),  J=1,...,N 
   *----------------------------------------------------------------------- 
   */
  cinur = 1.;
  cinui = 0.;
  if (*kode == 1)
    {
      goto L10;
    }
  cinur = cos (*zri);
  cinui = sin (*zri);
 L10:
  /*----------------------------------------------------------------------- 
   *    ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH 
   *    THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE 
   *    SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT 
   *    THE RESULT IS ON SCALE. 
   *----------------------------------------------------------------------- 
   */
  acw = amos_azabs (&cwr[2], &cwi[2]);
  ascle = pnl_d1mach (1) * 1e3 / *tol;
  csclr = 1.;
  if (acw > ascle)
    {
      goto L20;
    }
  csclr = 1. / *tol;
  goto L30;
 L20:
  ascle = 1. / ascle;
  if (acw < ascle)
    {
      goto L30;
    }
  csclr = *tol;
 L30:
  c1r = cwr[1] * csclr;
  c1i = cwi[1] * csclr;
  c2r = cwr[2] * csclr;
  c2i = cwi[2] * csclr;
  str = yr[1];
  sti = yi[1];
  /*----------------------------------------------------------------------- 
   *    CINU=CINU*(CONJG(CT)/CABS(CT))*(1.0D0/CABS(CT) PREVENTS 
   *    UNDER- OR OVERFLOW PREMATURELY BY SQUARING CABS(CT) 
   *----------------------------------------------------------------------- 
   */
  ptr = str * c1r - sti * c1i;
  pti = str * c1i + sti * c1r;
  ptr += c2r;
  pti += c2i;
  ctr = *zrr * ptr - *zri * pti;
  cti = *zrr * pti + *zri * ptr;
  act = amos_azabs (&ctr, &cti);
  ract = 1. / act;
  ctr *= ract;
  cti = -cti * ract;
  ptr = cinur * ract;
  pti = cinui * ract;
  cinur = ptr * ctr - pti * cti;
  cinui = ptr * cti + pti * ctr;
  yr[1] = cinur * csclr;
  yi[1] = cinui * csclr;
  if (*n == 1)
    {
      return 0;
    }
  i__1 = *n;
  for (i__ = 2; i__ <= i__1; ++i__)
    {
      ptr = str * cinur - sti * cinui;
      cinui = str * cinui + sti * cinur;
      cinur = ptr;
      str = yr[i__];
      sti = yi[i__];
      yr[i__] = cinur * csclr;
      yi[i__] = cinui * csclr;
      /* L40: */
    }
  return 0;
 L50:
  *nz = -1;
  if (nw == -2)
    {
      *nz = -2;
    }
  return 0;
}

