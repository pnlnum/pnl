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
 * gamma function,logarithm of gamma function 
 ****author  amos, donald e., sandia national laboratories 
 ****purpose  to compute the logarithm of the gamma function 
 ****description 
 * 
 *              **** a double precision routine **** 
 *        dgamln computes the natural log of the gamma function for 
 *        z.gt.0.  the asymptotic expansion is used to generate values 
 *        greater than zmin which are adjusted by the recursion 
 *        g(z+1)=z*g(z) for z.le.zmin.  the function was made as 
 *        portable as possible by computimg zmin from the number of base 
 *        10 digits in a word, rln=amax1(-alog10(r1mach(4)),0.5e-18) 
 *        limited to 18 digits of (relative) accuracy. 
 * 
 *        since int arguments are common, a table look up on 100 
 *        values is used for speed of execution. 
 * 
 *    description of arguments 
 * 
 *        input      z is d0uble precision 
 *          z      - argument, z.gt.0.0d0 
 * 
 *        output      dgamln is double precision 
 *          dgamln  - natural log of the gamma function at z.ne.0.0d0 
 *          ierr    - error flag 
 *                    ierr=0, normal return, computation completed 
 *                    ierr=1, z.le.0.0d0,    no computation 
 * 
 * 
 ****references  computation of bessel functions of complex argument 
 *                by d. e. amos, sand83-0083, may, 1983. 
 ****routines called  i1mach,d1mach 
 ****end prologue  dgamln 
 *          lngamma(n), n=1,100 
 */

double amos_dgamln (double *z__, int *ierr)
{
  static const double gln[100] =
    { 0., 0., .693147180559945309, 1.791759469228055, 3.17805383034794562,
      4.78749174278204599, 6.579251212010101, 8.5251613610654143, 10.6046029027452502,
      12.8018274800814696, 15.1044125730755153, 17.5023078458738858, 19.9872144956618861,
      22.5521638531234229, 25.1912211827386815, 27.8992713838408916, 30.6718601060806728,
      33.5050734501368889, 36.3954452080330536, 39.339884187199494, 42.335616460753485,
      45.380138898476908, 48.4711813518352239, 51.6066755677643736, 54.7847293981123192,
      58.0036052229805199, 61.261701761002002, 64.5575386270063311, 67.889743137181535,
      71.257038967168009, 74.6582363488301644, 78.0922235533153106, 81.5579594561150372,
      85.0544670175815174, 88.5808275421976788, 92.1361756036870925, 95.7196945421432025,
      99.3306124547874269, 102.968198614513813, 106.631760260643459, 110.320639714757395,
      114.034211781461703, 117.771881399745072, 121.533081515438634, 125.317271149356895,
      129.123933639127215, 132.95257503561631, 136.802722637326368, 140.673923648234259,
      144.565743946344886, 148.477766951773032, 152.409592584497358, 156.360836303078785,
      160.331128216630907, 164.320112263195181, 168.327445448427652, 172.352797139162802,
      176.395848406997352, 180.456291417543771, 184.533828861449491, 188.628173423671591,
      192.739047287844902, 196.866181672889994, 201.009316399281527, 205.168199482641199,
      209.342586752536836, 213.532241494563261, 217.736934113954227, 221.956441819130334,
      226.190548323727593, 230.439043565776952, 234.701723442818268, 238.978389561834323,
      243.268849002982714, 247.572914096186884, 251.890402209723194, 256.221135550009525,
      260.564940971863209, 264.921649798552801, 269.291097651019823, 273.673124285693704,
      278.067573440366143, 282.474292687630396, 286.893133295426994, 291.323950094270308,
      295.766601350760624, 300.220948647014132, 304.686856765668715, 309.164193580146922,
      313.652829949879062, 318.152639620209327, 322.663499126726177, 327.185287703775217,
      331.717887196928473, 336.261181979198477, 340.815058870799018, 345.379407062266854,
      349.954118040770237, 354.539085519440809, 359.134205369575399 };
  static const double cf[22] =
    { .0833333333333333333, -.00277777777777777778, 7.93650793650793651e-4,
      -5.95238095238095238e-4, 8.41750841750841751e-4, -.00191752691752691753,
      .00641025641025641026, -.0295506535947712418, .179644372368830573, -1.39243221690590112,
      13.402864044168392, -156.848284626002017, 2193.10333333333333, -36108.7712537249894,
      691472.268851313067, -15238221.5394074162, 382900751.391414141, -10882266035.7843911,
      347320283765.002252, -12369602142269.2745, 488788064793079.335, -21320333960919373.9 };
  static const double con = 1.83787706640934548;

  int i__1;
  double zinc, zmin, zdmy;
  int i__, k;
  double s, wdtol;
  double t1, fz, zm;
  int mz, nz=0;
  double zp;
  int i1m;
  double fln, tlg, rln, trm, tst, zsq;

  *ierr = 0;
  if (*z__ <= 0.)
    {
      goto L70;
    }
  if (*z__ > 101.)
    {
      goto L10;
    }
  nz = (int) (*z__);
  fz = *z__ - (double) nz;
  if (fz > 0.)
    {
      goto L10;
    }
  if (nz > 100)
    {
      goto L10;
    }
  return  gln[nz - 1];
 L10:
  wdtol = pnl_d1mach (4);
  wdtol = MAX (wdtol, 5e-19);
  i1m = amos_i1mach (14);
  rln = pnl_d1mach (5) * (double) i1m;
  fln = MIN (rln, 20.);
  fln = MAX (fln, 3.);
  fln += -3.;
  zm = fln * .3875 + 1.8;
  mz = (int) zm + 1;
  zmin = (double) mz;
  zdmy = *z__;
  zinc = 0.;
  if (*z__ >= zmin)
    {
      goto L20;
    }
  zinc = zmin - (double) nz;
  zdmy = *z__ + zinc;
 L20:
  zp = 1. / zdmy;
  t1 = cf[0] * zp;
  s = t1;
  if (zp < wdtol)
    {
      goto L40;
    }
  zsq = zp * zp;
  tst = t1 * wdtol;
  for (k = 2; k <= 22; ++k)
    {
      zp *= zsq;
      trm = cf[k - 1] * zp;
      if (fabs (trm) < tst)
	{
	  goto L40;
	}
      s += trm;
      /* L30: */
    }
 L40:
  if (zinc != 0.)
    {
      goto L50;
    }
  tlg = log (*z__);
  return  *z__ * (tlg - 1.) + (con - tlg) * .5 + s;
 L50:
  zp = 1.;
  nz = (int) zinc;
  i__1 = nz;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      zp *= *z__ + (double) (i__ - 1);
      /* L60: */
    }
  tlg = log (zdmy);
  return zdmy * (tlg - 1.) - log (zp) + (con - tlg) * .5 + s;
  /* 
   * 
   */
 L70:
  *ierr = 1;
  return 0.0;
}

