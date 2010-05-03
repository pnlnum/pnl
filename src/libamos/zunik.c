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
 *       zunik computes parameters for the uniform asymptotic 
 *       expansions of the i and k functions on ikflg= 1 or 2 
 *       respectively by 
 * 
 *       w(fnu,zr) = phi*exp(zeta)*sum 
 * 
 *       where       zeta=-zeta1 + zeta2       or 
 *                         zeta1 - zeta2 
 * 
 *       the first call must have init=0. subsequent calls with the 
 *       same zr and fnu will return the i or k function on ikflg= 
 *       1 or 2 with no change in init. cwrk is a complex work 
 *       array. ipmtr=0 computes all parameters. ipmtr=1 computes phi, 
 *       zeta1,zeta2. 
 *
 */


/* Subroutine */ int
amos_zunik (double *zrr, double *zri, double *fnu,const int *ikflg,const int *ipmtr,
	    double *tol, int *init, double *phir, double *phii,
	    double *zeta1r, double *zeta1i, double *zeta2r, double *zeta2i,
	    double *sumr, double *sumi, double *cwrkr, double *cwrki)
{
  /* Initialized data */

  static const double zeror = 0.;
  static const double zeroi = 0.;
  static const double coner = 1.;
  static const double conei = 0.;
  static const double con[2] = { .398942280401432678, 1.25331413731550025 };
  static const double c__[120] =
    { 1., -.208333333333333333, .125, .334201388888888889,
      -.401041666666666667, .0703125, -1.02581259645061728, 1.84646267361111111, -.8912109375,
      .0732421875, 4.66958442342624743, -11.2070026162229938, 8.78912353515625, -2.3640869140625,
      .112152099609375, -28.2120725582002449, 84.6362176746007346, -91.8182415432400174,
      42.5349987453884549, -7.3687943594796317, .227108001708984375, 212.570130039217123,
      -765.252468141181642, 1059.99045252799988, -699.579627376132541, 218.19051174421159,
      -26.4914304869515555, .572501420974731445, -1919.457662318407, 8061.72218173730938,
      -13586.5500064341374, 11655.3933368645332, -5305.64697861340311, 1200.90291321635246,
      -108.090919788394656, 1.7277275025844574, 20204.2913309661486, -96980.5983886375135,
      192547.001232531532, -203400.177280415534, 122200.46498301746, -41192.6549688975513,
      7109.51430248936372, -493.915304773088012, 6.07404200127348304, -242919.187900551333,
      1311763.6146629772, -2998015.91853810675, 3763271.297656404, -2813563.22658653411,
      1268365.27332162478, -331645.172484563578, 45218.7689813627263, -2499.83048181120962,
      24.3805296995560639, 3284469.85307203782, -19706819.1184322269, 50952602.4926646422,
      -74105148.2115326577, 66344512.2747290267, -37567176.6607633513, 13288767.1664218183,
      -2785618.12808645469, 308186.404612662398, -13886.0897537170405, 110.017140269246738,
      -49329253.664509962, 325573074.185765749, -939462359.681578403, 1553596899.57058006,
      -1621080552.10833708, 1106842816.82301447, -495889784.275030309, 142062907.797533095,
      -24474062.7257387285, 2243768.17792244943, -84005.4336030240853, 551.335896122020586,
      814789096.118312115, -5866481492.05184723, 18688207509.2958249, -34632043388.1587779,
      41280185579.753974, -33026599749.8007231, 17954213731.1556001, -6563293792.61928433,
      1559279864.87925751, -225105661.889415278, 17395107.5539781645, -549842.327572288687,
      3038.09051092238427, -14679261247.6956167, 114498237732.02581, -399096175224.466498,
      819218669548.577329, -1098375156081.22331, 1008158106865.38209, -645364869245.376503,
      287900649906.150589, -87867072178.0232657, 17634730606.8349694, -2167164983.22379509,
      143157876.718888981, -3871833.44257261262, 18257.7554742931747, 286464035717.679043,
      -2406297900028.50396, 9109341185239.89896, -20516899410934.4374, 30565125519935.3206,
      -31667088584785.1584, 23348364044581.8409, -12320491305598.2872, 4612725780849.13197,
      -1196552880196.1816, 205914503232.410016, -21822927757.5292237, 1247009293.51271032,
      -29188388.1222208134, 118838.426256783253 };

  /* System generated locals */
  int i__1;
  double d__1, d__2;

  /* Local variables */
  int idum;
  double test;
  int i__, j, k, l;
  double crfni, crfnr;
  double ac, si, ti, sr, tr, t2i;
  double t2r, rfn, sri, sti, zni, srr, str, znr;

  /* Parameter adjustments */
  --cwrki;
  --cwrkr;

  /* Function Body */
  /* 
   */
  if (*init != 0)
    {
      goto L40;
    }
  /*----------------------------------------------------------------------- 
   *    INITIALIZE ALL VARIABLES 
   *----------------------------------------------------------------------- 
   */
  rfn = 1. / *fnu;
  /*----------------------------------------------------------------------- 
   *    OVERFLOW TEST (ZR/FNU TOO SMALL) 
   *----------------------------------------------------------------------- 
   */
  test = pnl_d1mach (1) * 1e3;
  ac = *fnu * test;
  if (fabs (*zrr) > ac || fabs (*zri) > ac)
    {
      goto L15;
    }
  *zeta1r = (d__1 = log (test), fabs (d__1)) * 2. + *fnu;
  *zeta1i = 0.;
  *zeta2r = *fnu;
  *zeta2i = 0.;
  *phir = 1.;
  *phii = 0.;
  return 0;
 L15:
  tr = *zrr * rfn;
  ti = *zri * rfn;
  sr = coner + (tr * tr - ti * ti);
  si = conei + (tr * ti + ti * tr);
  amos_azsqrt (&sr, &si, &srr, &sri);
  str = coner + srr;
  sti = conei + sri;
  amos_zdiv (&str, &sti, &tr, &ti, &znr, &zni);
  amos_azlog (&znr, &zni, &str, &sti, &idum);
  *zeta1r = *fnu * str;
  *zeta1i = *fnu * sti;
  *zeta2r = *fnu * srr;
  *zeta2i = *fnu * sri;
  amos_zdiv (&coner, &conei, &srr, &sri, &tr, &ti);
  srr = tr * rfn;
  sri = ti * rfn;
  amos_azsqrt (&srr, &sri, &cwrkr[16], &cwrki[16]);
  *phir = cwrkr[16] * con[*ikflg - 1];
  *phii = cwrki[16] * con[*ikflg - 1];
  if (*ipmtr != 0)
    {
      return 0;
    }
  amos_zdiv (&coner, &conei, &sr, &si, &t2r, &t2i);
  cwrkr[1] = coner;
  cwrki[1] = conei;
  crfnr = coner;
  crfni = conei;
  ac = 1.;
  l = 1;
  for (k = 2; k <= 15; ++k)
    {
      sr = zeror;
      si = zeroi;
      i__1 = k;
      for (j = 1; j <= i__1; ++j)
	{
	  ++l;
	  str = sr * t2r - si * t2i + c__[l - 1];
	  si = sr * t2i + si * t2r;
	  sr = str;
	  /* L10: */
	}
      str = crfnr * srr - crfni * sri;
      crfni = crfnr * sri + crfni * srr;
      crfnr = str;
      cwrkr[k] = crfnr * sr - crfni * si;
      cwrki[k] = crfnr * si + crfni * sr;
      ac *= rfn;
      test = (d__1 = cwrkr[k], fabs (d__1)) + (d__2 = cwrki[k], fabs (d__2));
      if (ac < *tol && test < *tol)
	{
	  goto L30;
	}
      /* L20: */
    }
  k = 15;
 L30:
  *init = k;
 L40:
  if (*ikflg == 2)
    {
      goto L60;
    }
  /*----------------------------------------------------------------------- 
   *    COMPUTE SUM FOR THE I FUNCTION 
   *----------------------------------------------------------------------- 
   */
  sr = zeror;
  si = zeroi;
  i__1 = *init;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      sr += cwrkr[i__];
      si += cwrki[i__];
      /* L50: */
    }
  *sumr = sr;
  *sumi = si;
  *phir = cwrkr[16] * con[0];
  *phii = cwrki[16] * con[0];
  return 0;
 L60:
  /*----------------------------------------------------------------------- 
   *    COMPUTE SUM FOR THE K FUNCTION 
   *----------------------------------------------------------------------- 
   */
  sr = zeror;
  si = zeroi;
  tr = coner;
  i__1 = *init;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      sr += tr * cwrkr[i__];
      si += tr * cwrki[i__];
      tr = -tr;
      /* L70: */
    }
  *sumr = sr;
  *sumi = si;
  *phir = cwrkr[16] * con[1];
  *phii = cwrki[16] * con[1];
  return 0;
}				/* zunik_ */
