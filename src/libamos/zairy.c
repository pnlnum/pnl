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

/*
 *  airy function,bessel functions of order one third 
 *  amos, donald e., sandia national laboratories 
 * purpose  to compute airy functions ai(z) and dai(z) for complex z 
 * 
 *                     ***a double precision routine*** 
 *        on kode=1, zairy computes the complex airy function ai(z) or 
 *        its derivative dai(z)/dz on id=0 or id=1 respectively. on 
 *        kode=2, a scaling option cexp(zta)*ai(z) or cexp(zta)* 
 *        dai(z)/dz is provided to remove the exponential decay in 
 *        -pi/3.lt.arg(z).lt.pi/3 and the exponential growth in 
 *        pi/3.lt.abs(arg(z)).lt.pi where zta=(2/3)*z*csqrt(z). 
 * 
 *        while the airy functions ai(z) and dai(z)/dz are analytic in 
 *        the whole z plane, the corresponding scaled functions defined 
 *        for kode=2 have a cut along the negative real axis. 
 *        defintions and notation are found in the nbs handbook of 
 *        mathematical functions (ref. 1). 
 * 
 *        input      zr,zi are double precision 
 *          zr,zi  - z=cmplx(zr,zi) 
 *          id     - order of derivative, id=0 or id=1 
 *          kode   - a parameter to indicate the scaling option 
 *                   kode= 1  returns 
 *                            ai=ai(z)                on id=0 or 
 *                            ai=dai(z)/dz            on id=1 
 *                       = 2  returns 
 *                            ai=cexp(zta)*ai(z)       on id=0 or 
 *                            ai=cexp(zta)*dai(z)/dz   on id=1 where 
 *                            zta=(2/3)*z*csqrt(z) 
 * 
 *        output     air,aii are double precision 
 *          air,aii- complex answer depending on the choices for id and 
 *                   kode 
 *          nz     - underflow indicator 
 *                   nz= 0   , normal return 
 *                   nz= 1   , ai=cmplx(0.0d0,0.0d0) due to underflow in 
 *                             -pi/3.lt.arg(z).lt.pi/3 on kode=1 
 *          ierr   - error flag 
 *                   ierr=0, normal return - computation completed 
 *                   ierr=1, input error   - no computation 
 *                   ierr=2, overflow      - no computation, real(zta) 
 *                           too large on kode=1 
 *                   ierr=3, cabs(z) large      - computation completed 
 *                           losses of signifcance by argument reduction 
 *                           produce less than half of machine accuracy 
 *                   ierr=4, cabs(z) too large  - no computation 
 *                           complete loss of accuracy by argument 
 *                           reduction 
 *                   ierr=5, error              - no computation, 
 *                           algorithm termination condition not met 
 * 
 ****long description 
 * 
 *        ai and dai are computed for cabs(z).gt.1.0 from the k bessel 
 *        functions by 
 * 
 *           ai(z)=c*sqrt(z)*k(1/3,zta) , dai(z)=-c*z*k(2/3,zta) 
 *                          c=1.0/(pi*sqrt(3.0)) 
 *                           zta=(2/3)*z**(3/2) 
 * 
 *        with the power series for cabs(z).le.1.0. 
 * 
 *        in most complex variable computation, one must evaluate ele- 
 *        mentary functions. when the magnitude of z is large, losses 
 *        of significance by argument reduction occur. consequently, if 
 *        the magnitude of zeta=(2/3)*z**1.5 exceeds u1=sqrt(0.5/ur), 
 *        then losses exceeding half precision are likely and an error 
 *        flag ierr=3 is triggered where ur=dmax1(d1mach(4),1.0d-18) is 
 *        double precision unit roundoff limited to 18 digits precision. 
 *        also, if the magnitude of zeta is larger than u2=0.5/ur, then 
 *        all significance is lost and ierr=4. in order to use the int 
 *        function, zeta must be further restricted not to exceed the 
 *        largest int, u3=i1mach(9). thus, the magnitude of zeta 
 *        must be restricted by min(u2,u3). on 32 bit machines, u1,u2, 
 *        and u3 are approximately 2.0e+3, 4.2e+6, 2.1e+9 in single 
 *        precision arithmetic and 1.3e+8, 1.8e+16, 2.1e+9 in double 
 *        precision arithmetic respectively. this makes u2 and u3 limit- 
 *        ing in their respective arithmetics. this means that the mag- 
 *        nitude of z cannot exceed 3.1e+4 in single and 2.1e+6 in 
 *        double precision arithmetic. this also means that one can 
 *        expect to retain, in the worst cases on 32 bit machines, 
 *        no digits in single precision and only 7 digits in double 
 *        precision arithmetic. similar considerations hold for other 
 *        machines. 
 * 
 *        the approximate relative error in the magnitude of a complex 
 *        bessel function can be expressed by p*10**s where p=max(unit 
 *        roundoff,1.0e-18) is the nominal precision and 10**s repre- 
 *        sents the increase in error due to argument reduction in the 
 *        elementary functions. here, s=max(1,abs(log10(cabs(z))), 
 *        abs(log10(fnu))) approximately (i.e. s=max(1,abs(exponent of 
 *        cabs(z),abs(exponent of fnu)) ). however, the phase angle may 
 *        have only absolute accuracy. this is most likely to occur when 
 *        one component (in absolute value) is larger than the other by 
 *        several orders of magnitude. if one component is 10**k larger 
 *        than the other, then one can expect only max(abs(log10(p))-k, 
 *        0) significant digits; or, stated another way, when k exceeds 
 *        the exponent of p, no significant digits remain in the smaller 
 *        component. however, the phase angle retains absolute accuracy 
 *        because, in complex arithmetic with precision p, the smaller 
 *        component will not (as a rule) decrease below p times the 
 *        magnitude of the larger component. in these extreme cases, 
 *        the principal phase angle is on the order of +p, -p, pi/2-p, 
 *        or -pi/2+p. 
 * 
 * references:  handbook of mathematical functions by m. abramowitz 
 *                and i. a. stegun, nbs ams series 55, u.s. dept. of 
 *                commerce, 1955. 
 * 
 *              computation of bessel functions of complex argument 
 *                and large order by d. e. amos, sand83-0643, may, 1983 
 * 
 *              a subroutine package for bessel functions of a complex 
 *                argument and nonnegative order by d. e. amos, sand85- 
 *                1018, may, 1985 
 * 
 *              a portable package for bessel functions of a complex 
 *                argument and nonnegative order by d. e. amos, trans. 
 *                math. software, 1986 
 * 
 */

int pnl_zairy (double *zr, double *zi,const int *id,const int *kode, double *air,
		double *aii, int *nz, int *ierr)
{

  static const double tth = .666666666666666667;
  static const double c1 = .35502805388781724;
  static const double c2 = .258819403792806799;
  static const double coef = .183776298473930683;
  static const double zeror = 0.;
  static const double zeroi = 0.;
  static const double coner = 1.;
  static const double conei = 0.;

  /* System generated locals */
  int i__1, i__2;
  double d__1;

  /* Local variables */
  double sfac, alim, elim, alaz, csqi, atrm, ztai, csqr, ztar, trm1i, trm2i,
    trm1r, trm2r;
  int k, iflag;
  double d1, d2;
  int k1;
  int k2;
  double aa, bb, ad, cc, ak, bk, ck, dk, az;
  int nn;
  double rl;
  int mr;
  double s1i, az3, s2i;
  double s1r, s2r, z3i, z3r, dig, fid, cyi[1], r1m5, fnu, cyr[1], tol, sti,
    ptr, str;

  *ierr = 0;
  *nz = 0;
  if (*id < 0 || *id > 1)
    {
      *ierr = 1;
    }
  if (*kode < 1 || *kode > 2)
    {
      *ierr = 1;
    }
  if (*ierr != 0)
    {
      return 0;
    }
  az = amos_azabs (zr, zi);
  /*Computing MAX 
   */
  d__1 = pnl_d1mach (4);	/*  */
  tol = MAX (d__1, 1e-18);
  fid = (double) (*id);
  if (az > 1.)
    {
      goto L70;
    }
  /*----------------------------------------------------------------------- 
   *    POWER SERIES FOR CABS(Z).LE.1. 
   *----------------------------------------------------------------------- 
   */
  s1r = coner;
  s1i = conei;
  s2r = coner;
  s2i = conei;
  if (az < tol)
    {
      goto L170;
    }
  aa = az * az;
  if (aa < tol / az)
    {
      goto L40;
    }
  trm1r = coner;
  trm1i = conei;
  trm2r = coner;
  trm2i = conei;
  atrm = 1.;
  str = *zr * *zr - *zi * *zi;
  sti = *zr * *zi + *zi * *zr;
  z3r = str * *zr - sti * *zi;
  z3i = str * *zi + sti * *zr;
  az3 = az * aa;
  ak = fid + 2.;
  bk = 3. - fid - fid;
  ck = 4. - fid;
  dk = fid + 3. + fid;
  d1 = ak * dk;
  d2 = bk * ck;
  ad = MIN (d1, d2);
  ak = fid * 9. + 24.;
  bk = 30. - fid * 9.;
  for (k = 1; k <= 25; ++k)
    {
      str = (trm1r * z3r - trm1i * z3i) / d1;
      trm1i = (trm1r * z3i + trm1i * z3r) / d1;
      trm1r = str;
      s1r += trm1r;
      s1i += trm1i;
      str = (trm2r * z3r - trm2i * z3i) / d2;
      trm2i = (trm2r * z3i + trm2i * z3r) / d2;
      trm2r = str;
      s2r += trm2r;
      s2i += trm2i;
      atrm = atrm * az3 / ad;
      d1 += ak;
      d2 += bk;
      ad = MIN (d1, d2);
      if (atrm < tol * ad)
	{
	  goto L40;
	}
      ak += 18.;
      bk += 18.;
      /* L30: */
    }
 L40:
  if (*id == 1)
    {
      goto L50;
    }
  *air = s1r * c1 - c2 * (*zr * s2r - *zi * s2i);
  *aii = s1i * c1 - c2 * (*zr * s2i + *zi * s2r);
  if (*kode == 1)
    {
      return 0;
    }
  amos_azsqrt (zr, zi, &str, &sti);
  ztar = tth * (*zr * str - *zi * sti);
  ztai = tth * (*zr * sti + *zi * str);
  amos_azexp (&ztar, &ztai, &str, &sti);
  ptr = *air * str - *aii * sti;
  *aii = *air * sti + *aii * str;
  *air = ptr;
  return 0;
 L50:
  *air = -s2r * c2;
  *aii = -s2i * c2;
  if (az <= tol)
    {
      goto L60;
    }
  str = *zr * s1r - *zi * s1i;
  sti = *zr * s1i + *zi * s1r;
  cc = c1 / (fid + 1.);
  *air += cc * (str * *zr - sti * *zi);
  *aii += cc * (str * *zi + sti * *zr);
 L60:
  if (*kode == 1)
    {
      return 0;
    }
  amos_azsqrt (zr, zi, &str, &sti);
  ztar = tth * (*zr * str - *zi * sti);
  ztai = tth * (*zr * sti + *zi * str);
  amos_azexp (&ztar, &ztai, &str, &sti);
  ptr = str * *air - sti * *aii;
  *aii = str * *aii + sti * *air;
  *air = ptr;
  return 0;
  /*----------------------------------------------------------------------- 
   *    CASE FOR CABS(Z).GT.1.0 
   *----------------------------------------------------------------------- 
   */
 L70:
  fnu = (fid + 1.) / 3.;
  /*----------------------------------------------------------------------- 
   *    SET PARAMETERS RELATED TO MACHINE CONSTANTS. 
   *    TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0D-18. 
   *    ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT. 
   *    EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND 
   *    EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR 
   *    UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE. 
   *    RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z. 
   *    DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG). 
   *----------------------------------------------------------------------- 
   */
  k1 = amos_i1mach (15);
  k2 = amos_i1mach (16);
  r1m5 = pnl_d1mach (5);
  /*Computing MIN 
   */
  i__1 = abs (k1), i__2 = abs (k2);
  k = MIN (i__1, i__2);
  elim = ((double) k * r1m5 - 3.) * 2.303;
  k1 = amos_i1mach (14) - 1;
  aa = r1m5 * (double) k1;
  dig = MIN (aa, 18.);
  aa *= 2.303;
  /*Computing MAX 
   */
  d__1 = -aa;
  alim = elim + MAX (d__1, -41.45);
  rl = dig * 1.2 + 3.;
  alaz = log (az);
  /*-------------------------------------------------------------------------- 
   *    TEST FOR PROPER RANGE 
   *----------------------------------------------------------------------- 
   */
  aa = .5 / tol;
  bb = (double) amos_i1mach (9) * .5;
  aa = MIN (aa, bb);
  aa = pow (aa, tth);
  if (az > aa)
    {
      goto L260;
    }
  aa = sqrt (aa);
  if (az > aa)
    {
      *ierr = 3;
    }
  amos_azsqrt (zr, zi, &csqr, &csqi);
  ztar = tth * (*zr * csqr - *zi * csqi);
  ztai = tth * (*zr * csqi + *zi * csqr);
  /*----------------------------------------------------------------------- 
   *    RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL 
   *----------------------------------------------------------------------- 
   */
  iflag = 0;
  sfac = 1.;
  ak = ztai;
  if (*zr >= 0.)
    {
      goto L80;
    }
  bk = ztar;
  ck = - fabs(bk);
  ztar = ck;
  ztai = ak;
 L80:
  if (*zi != 0.)
    {
      goto L90;
    }
  if (*zr > 0.)
    {
      goto L90;
    }
  ztar = 0.;
  ztai = ak;
 L90:
  aa = ztar;
  if (aa >= 0. && *zr > 0.)
    {
      goto L110;
    }
  if (*kode == 2)
    {
      goto L100;
    }
  /*----------------------------------------------------------------------- 
   *    OVERFLOW TEST 
   *----------------------------------------------------------------------- 
   */
  if (aa > -alim)
    {
      goto L100;
    }
  aa = -aa + alaz * .25;
  iflag = 1;
  sfac = tol;
  if (aa > elim)
    {
      goto L270;
    }
 L100:
  /*----------------------------------------------------------------------- 
   *    CBKNU AND CACON RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2 
   *----------------------------------------------------------------------- 
   */
  mr = 1;
  if (*zi < 0.)
    {
      mr = -1;
    }
  amos_zacai (&ztar, &ztai, &fnu, kode, &mr, &c__1, cyr, cyi, &nn, &rl, &tol,  &elim, &alim);
  if (nn < 0)
    {
      goto L280;
    }
  *nz += nn;
  goto L130;
 L110:
  if (*kode == 2)
    {
      goto L120;
    }
  /*----------------------------------------------------------------------- 
   *    UNDERFLOW TEST 
   *----------------------------------------------------------------------- 
   */
  if (aa < alim)
    {
      goto L120;
    }
  aa = -aa - alaz * .25;
  iflag = 2;
  sfac = 1. / tol;
  if (aa < -elim)
    {
      goto L210;
    }
 L120:
  amos_zbknu (&ztar, &ztai, &fnu, kode, &c__1, cyr, cyi, nz, &tol, &elim,
	      &alim);
 L130:
  s1r = cyr[0] * coef;
  s1i = cyi[0] * coef;
  if (iflag != 0)
    {
      goto L150;
    }
  if (*id == 1)
    {
      goto L140;
    }
  *air = csqr * s1r - csqi * s1i;
  *aii = csqr * s1i + csqi * s1r;
  return 0;
 L140:
  *air = -(*zr * s1r - *zi * s1i);
  *aii = -(*zr * s1i + *zi * s1r);
  return 0;
 L150:
  s1r *= sfac;
  s1i *= sfac;
  if (*id == 1)
    {
      goto L160;
    }
  str = s1r * csqr - s1i * csqi;
  s1i = s1r * csqi + s1i * csqr;
  s1r = str;
  *air = s1r / sfac;
  *aii = s1i / sfac;
  return 0;
 L160:
  str = -(s1r * *zr - s1i * *zi);
  s1i = -(s1r * *zi + s1i * *zr);
  s1r = str;
  *air = s1r / sfac;
  *aii = s1i / sfac;
  return 0;
 L170:
  aa = pnl_d1mach (1) * 1e3;
  s1r = zeror;
  s1i = zeroi;
  if (*id == 1)
    {
      goto L190;
    }
  if (az <= aa)
    {
      goto L180;
    }
  s1r = c2 * *zr;
  s1i = c2 * *zi;
 L180:
  *air = c1 - s1r;
  *aii = -s1i;
  return 0;
 L190:
  *air = -c2;
  *aii = 0.;
  aa = sqrt (aa);
  if (az <= aa)
    {
      goto L200;
    }
  s1r = (*zr * *zr - *zi * *zi) * .5;
  s1i = *zr * *zi;
 L200:
  *air += c1 * s1r;
  *aii += c1 * s1i;
  return 0;
 L210:
  *nz = 1;
  *air = zeror;
  *aii = zeroi;
  return 0;
 L270:
  *nz = 0;
  *ierr = 2;
  return 0;
 L280:
  if (nn == -1)
    {
      goto L270;
    }
  *nz = 0;
  *ierr = 5;
  return 0;
 L260:
  *ierr = 4;
  *nz = 0;
  return 0;
}				/* zairy_ */
