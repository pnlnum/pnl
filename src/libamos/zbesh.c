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
 * h-bessel functions,bessel functions of complex argument, 
 *          bessel functions of third kind,hankel functions 
 * Amos, Donald E., sandia national laboratories 
 *
 * description 
 * 
 *                     ***a double precision routine*** 
 *        on kode=1, zbesh computes an n member sequence of complex 
 *        hankel (bessel) functions cy(j)=h(m,fnu+j-1,z) for kinds m=1 
 *        or 2, real, nonnegative orders fnu+j-1, j=1,...,n, and complex 
 *        z.ne.cmplx(0.0,0.0) in the cut plane -pi.lt.arg(z).le.pi. 
 *        on kode=2, zbesh returns the scaled hankel functions 
 * 
 *        cy(i)=exp(-mm*z*i)*h(m,fnu+j-1,z)       mm=3-2*m,   i**2=-1. 
 * 
 *        which removes the exponential behavior in both the upper and 
 *        lower half planes. definitions and notation are found in the 
 *        nbs handbook of mathematical functions (ref. 1). 
 * 
 *        input      zr,zi,fnu are double precision 
 *          zr,zi  - z=cmplx(zr,zi), z.ne.cmplx(0.0d0,0.0d0), 
 *                   -pt.lt.arg(z).le.pi 
 *          fnu    - order of initial h function, fnu.ge.0.0d0 
 *          kode   - a parameter to indicate the scaling option 
 *                   kode= 1  returns 
 *                            cy(j)=h(m,fnu+j-1,z),   j=1,...,n 
 *                       = 2  returns 
 *                            cy(j)=h(m,fnu+j-1,z)*exp(-i*z*(3-2m)) 
 *                                 j=1,...,n  ,  i**2=-1 
 *          m      - kind of hankel function, m=1 or 2 
 *          n      - number of members in the sequence, n.ge.1 
 * 
 *        output     cyr,cyi are double precision 
 *          cyr,cyi- double precision vectors whose first n components 
 *                   contain real and imaginary parts for the sequence 
 *                   cy(j)=h(m,fnu+j-1,z)  or 
 *                   cy(j)=h(m,fnu+j-1,z)*exp(-i*z*(3-2m))  j=1,...,n 
 *                   depending on kode, i**2=-1. 
 *          nz     - number of components set to zero due to underflow, 
 *                   nz= 0   , normal return 
 *                   nz.gt.0 , first nz components of cy set to zero due 
 *                             to underflow, cy(j)=cmplx(0.0d0,0.0d0) 
 *                             j=1,...,nz when y.gt.0.0 and m=1 or 
 *                             y.lt.0.0 and m=2. for the complmentary 
 *                             half planes, nz states only the number 
 *                             of underflows. 
 *          ierr   - error flag 
 *                   ierr=0, normal return - computation completed 
 *                   ierr=1, input error   - no computation 
 *                   ierr=2, overflow      - no computation, fnu too 
 *                           large or cabs(z) too small or both 
 *                   ierr=3, cabs(z) or fnu+n-1 large - computation done 
 *                           but losses of signifcance by argument 
 *                           reduction produce less than half of machine 
 *                           accuracy 
 *                   ierr=4, cabs(z) or fnu+n-1 too large - no computa- 
 *                           tion because of complete losses of signifi- 
 *                           cance by argument reduction 
 *                   ierr=5, error              - no computation, 
 *                           algorithm termination condition not met 
 * 
 ****long description 
 * 
 *        the computation is carried out by the relation 
 * 
 *        h(m,fnu,z)=(1/mp)*exp(-mp*fnu)*k(fnu,z*exp(-mp)) 
 *            mp=mm*hpi*i,  mm=3-2*m,  hpi=pi/2,  i**2=-1 
 * 
 *        for m=1 or 2 where the k bessel function is computed for the 
 *        right half plane re(z).ge.0.0. the k function is continued 
 *        to the left half plane by the relation 
 * 
 *        k(fnu,z*exp(mp)) = exp(-mp*fnu)*k(fnu,z)-mp*i(fnu,z) 
 *        mp=mr*pi*i, mr=+1 or -1, re(z).gt.0, i**2=-1 
 * 
 *        where i(fnu,z) is the i bessel function. 
 * 
 *        exponential decay of h(m,fnu,z) occurs in the upper half z 
 *        plane for m=1 and the lower half z plane for m=2.  exponential 
 *        growth occurs in the complementary half planes.  scaling 
 *        by exp(-mm*z*i) removes the exponential behavior in the 
 *        whole z plane for z to infinity. 
 * 
 *        for negative orders,the formulae 
 * 
 *              h(1,-fnu,z) = h(1,fnu,z)*cexp( pi*fnu*i) 
 *              h(2,-fnu,z) = h(2,fnu,z)*cexp(-pi*fnu*i) 
 *                        i**2=-1 
 * 
 *        can be used. 
 * 
 *        in most complex variable computation, one must evaluate ele- 
 *        mentary functions. when the magnitude of z or fnu+n-1 is 
 *        large, losses of significance by argument reduction occur. 
 *        consequently, if either one exceeds u1=sqrt(0.5/ur), then 
 *        losses exceeding half precision are likely and an error flag 
 *        ierr=3 is triggered where ur=dmax1(d1mach(4),1.0d-18) is 
 *        double precision unit roundoff limited to 18 digits precision. 
 *        if either is larger than u2=0.5/ur, then all significance is 
 *        lost and ierr=4. in order to use the int function, arguments 
 *        must be further restricted not to exceed the largest machine 
 *        int, u3=i1mach(9). thus, the magnitude of z and fnu+n-1 is 
 *        restricted by min(u2,u3). on 32 bit machines, u1,u2, and u3 
 *        are approximately 2.0e+3, 4.2e+6, 2.1e+9 in single precision 
 *        arithmetic and 1.3e+8, 1.8e+16, 2.1e+9 in double precision 
 *        arithmetic respectively. this makes u2 and u3 limiting in 
 *        their respective arithmetics. this means that one can expect 
 *        to retain, in the worst cases on 32 bit machines, no digits 
 *        in single and only 7 digits in double precision arithmetic. 
 *        similar considerations hold for other machines. 
 * 
 *        the approximate relative error in the magnitude of a complex 
 *        bessel function can be expressed by p*10**s where p=max(unit 
 *        roundoff,1.0d-18) is the nominal precision and 10**s repre- 
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
 ****references  handbook of mathematical functions by m. abramowitz 
 *                and i. a. stegun, nbs ams series 55, u.s. dept. of 
 *                commerce, 1955. 
 * 
 *              computation of bessel functions of complex argument 
 *                by d. e. amos, sand83-0083, may, 1983. 
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

int pnl_zbesh (double *zr, double *zi, double *fnu, int *kode,const int *m,const int *n,
		double *cyr, double *cyi, int *nz, int *ierr)
{
  static const double hpi = 1.57079632679489662;

  int i__1, i__2;
  double d__1, d__2;

  /* Local variables */
  double alim, elim, atol, rhpi;
  int inuh;
  double fnul, rtol;
  int i__, k;
  double ascle;
  double csgni;
  double csgnr;
  int k1;
  int k2;
  double aa, bb, fn;
  int mm;
  double az;
  int ir, nn;
  double rl;
  int mr, nw;
  double dig, arg, aln, fmm, r1m5, ufl, sgn;
  int nuf, inu;
  double tol, sti, zni, zti, str, znr;

  /* Parameter adjustments */
  --cyi;
  --cyr;

  /* Function Body */
  /* 
****FIRST EXECUTABLE STATEMENT  ZBESH 
*/
  *ierr = 0;
  *nz = 0;
  if (*zr == 0. && *zi == 0.)
    {
      *ierr = 1;
    }
  if (*fnu < 0.)
    {
      *ierr = 1;
    }
  if (*m < 1 || *m > 2)
    {
      *ierr = 1;
    }
  if (*kode < 1 || *kode > 2)
    {
      *ierr = 1;
    }
  if (*n < 1)
    {
      *ierr = 1;
    }
  if (*ierr != 0)
    {
      return 0;
    }
  nn = *n;
  /*----------------------------------------------------------------------- 
   *    SET PARAMETERS RELATED TO MACHINE CONSTANTS. 
   *    TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18. 
   *    ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT. 
   *    EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND 
   *    EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR 
   *    UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE. 
   *    RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z. 
   *    DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG). 
   *    FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU 
   *----------------------------------------------------------------------- 
   *Computing MAX 
   */
  d__1 = pnl_d1mach (4);
  tol = MAX (d__1, 1e-18);
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
  fnul = (dig - 3.) * 6. + 10.;
  rl = dig * 1.2 + 3.;
  fn = *fnu + (double) (nn - 1);
  mm = 3 - *m - *m;
  fmm = (double) mm;
  znr = fmm * *zi;
  zni = -fmm * *zr;
  /*----------------------------------------------------------------------- 
   *    TEST FOR PROPER RANGE 
   *----------------------------------------------------------------------- 
   */
  az = amos_azabs (zr, zi);
  aa = .5 / tol;
  bb = (double) amos_i1mach (9) * .5;
  aa = MIN (aa, bb);
  if (az > aa)
    {
      goto L260;
    }
  if (fn > aa)
    {
      goto L260;
    }
  aa = sqrt (aa);
  if (az > aa)
    {
      *ierr = 3;
    }
  if (fn > aa)
    {
      *ierr = 3;
    }
  /*----------------------------------------------------------------------- 
   *    OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE 
   *----------------------------------------------------------------------- 
   */
  ufl = pnl_d1mach (1) * 1e3;
  if (az < ufl)
    {
      goto L230;
    }
  if (*fnu > fnul)
    {
      goto L90;
    }
  if (fn <= 1.)
    {
      goto L70;
    }
  if (fn > 2.)
    {
      goto L60;
    }
  if (az > tol)
    {
      goto L70;
    }
  arg = az * .5;
  aln = -fn * log (arg);
  if (aln > elim)
    {
      goto L230;
    }
  goto L70;
 L60:
  amos_zuoik (&znr, &zni, fnu, kode, &c__2, &nn, &cyr[1], &cyi[1], &nuf, &tol,
	      &elim, &alim);
  if (nuf < 0)
    {
      goto L230;
    }
  *nz += nuf;
  nn -= nuf;
  /*----------------------------------------------------------------------- 
   *    HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK 
   *    IF NUF=NN, THEN CY(I)=CZERO FOR ALL I 
   *----------------------------------------------------------------------- 
   */
  if (nn == 0)
    {
      goto L140;
    }
 L70:
  if (znr < 0. || (znr == 0. && zni < 0. && *m == 2) )
    {
      goto L80;
    }
  /*----------------------------------------------------------------------- 
   *    RIGHT HALF PLANE COMPUTATION, XN.GE.0. .AND. (XN.NE.0. .OR. 
   *    YN.GE.0. .OR. M=1) 
   *----------------------------------------------------------------------- 
   */
  amos_zbknu (&znr, &zni, fnu, kode, &nn, &cyr[1], &cyi[1], nz, &tol, &elim,
	      &alim);
  goto L110;
  /*----------------------------------------------------------------------- 
   *    LEFT HALF PLANE COMPUTATION 
   *----------------------------------------------------------------------- 
   */
 L80:
  mr = -mm;
  amos_zacon (&znr, &zni, fnu, kode, &mr, &nn, &cyr[1], &cyi[1], &nw, &rl,
	      &fnul, &tol, &elim, &alim);
  if (nw < 0)
    {
      goto L240;
    }
  *nz = nw;
  goto L110;
 L90:
  /*----------------------------------------------------------------------- 
   *    UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL 
   *----------------------------------------------------------------------- 
   */
  mr = 0;
  if (znr >= 0. && (znr != 0. || zni >= 0. || *m != 2))
    {
      goto L100;
    }
  mr = -mm;
  if (znr != 0. || zni >= 0.)
    {
      goto L100;
    }
  znr = -znr;
  zni = -zni;
 L100:
  amos_zbunk (&znr, &zni, fnu, kode, &mr, &nn, &cyr[1], &cyi[1], &nw, &tol,
	      &elim, &alim);
  if (nw < 0)
    {
      goto L240;
    }
  *nz += nw;
 L110:
  /*----------------------------------------------------------------------- 
   *    H(M,FNU,Z) = -FMM*(I/HPI)*(ZT**FNU)*K(FNU,-Z*ZT) 
   * 
   *    ZT=EXP(-FMM*HPI*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2 
   *----------------------------------------------------------------------- 
   */
  d__1 = -fmm;
  sgn = D_SIGN (hpi, d__1);	/*  */
  /*----------------------------------------------------------------------- 
   *    CALCULATE EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE 
   *    WHEN FNU IS LARGE 
   *----------------------------------------------------------------------- 
   */
  inu = (int) (*fnu);
  inuh = inu / 2;
  ir = inu - (inuh << 1);
  arg = (*fnu - (double) (inu - ir)) * sgn;
  rhpi = 1. / sgn;
  /*    ZNI = RHPI*DCOS(ARG) 
   *    ZNR = -RHPI*DSIN(ARG) 
   */
  csgni = rhpi * cos (arg);
  csgnr = -rhpi * sin (arg);
  if (inuh % 2 == 0)
    {
      goto L120;
    }
  /*    ZNR = -ZNR 
   *    ZNI = -ZNI 
   */
  csgnr = -csgnr;
  csgni = -csgni;
 L120:
  zti = -fmm;
  rtol = 1. / tol;
  ascle = ufl * rtol;
  i__1 = nn;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /*      STR = CYR(I)*ZNR - CYI(I)*ZNI 
       *      CYI(I) = CYR(I)*ZNI + CYI(I)*ZNR 
       *      CYR(I) = STR 
       *      STR = -ZNI*ZTI 
       *      ZNI = ZNR*ZTI 
       *      ZNR = STR 
       */
      aa = cyr[i__];
      bb = cyi[i__];
      atol = 1.;
      /*Computing MAX 
       */
      d__1 = fabs (aa), d__2 = fabs (bb);
      if (MAX (d__1, d__2) > ascle)
	{
	  goto L135;
	}
      aa *= rtol;
      bb *= rtol;
      atol = tol;
    L135:
      str = aa * csgnr - bb * csgni;
      sti = aa * csgni + bb * csgnr;
      cyr[i__] = str * atol;
      cyi[i__] = sti * atol;
      str = -csgni * zti;
      csgni = csgnr * zti;
      csgnr = str;
      /* L130: */
    }
  return 0;
 L140:
  if (znr < 0.)
    {
      goto L230;
    }
  return 0;
 L230:
  *nz = 0;
  *ierr = 2;
  return 0;
 L240:
  if (nw == -1)
    {
      goto L230;
    }
  *nz = 0;
  *ierr = 5;
  return 0;
 L260:
  *nz = 0;
  *ierr = 4;
  return 0;
}				/* zbesh_ */
