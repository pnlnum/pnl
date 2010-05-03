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
 * i-bessel function,complex bessel function, 
 *            modified bessel function of the first kind 
 ****author  amos, donald e., sandia national laboratories 
 ****purpose  to compute i-bessel functions of complex argument 
 ****description 
 * 
 *                   ***a double precision routine*** 
 *        on kode=1, zbesi computes an n member sequence of complex 
 *        bessel functions cy(j)=i(fnu+j-1,z) for real, nonnegative 
 *        orders fnu+j-1, j=1,...,n and complex z in the cut plane 
 *        -pi.lt.arg(z).le.pi. on kode=2, zbesi returns the scaled 
 *        functions 
 * 
 *        cy(j)=exp(-abs(x))*i(fnu+j-1,z)   j = 1,...,n , x=real(z) 
 * 
 *        with the exponential growth removed in both the left and 
 *        right half planes for z to infinity. definitions and notation 
 *        are found in the nbs handbook of mathematical functions 
 *        (ref. 1). 
 * 
 *        input      zr,zi,fnu are double precision 
 *          zr,zi  - z=cmplx(zr,zi),  -pi.lt.arg(z).le.pi 
 *          fnu    - order of initial i function, fnu.ge.0.0d0 
 *          kode   - a parameter to indicate the scaling option 
 *                   kode= 1  returns 
 *                            cy(j)=i(fnu+j-1,z), j=1,...,n 
 *                       = 2  returns 
 *                            cy(j)=i(fnu+j-1,z)*exp(-abs(x)), j=1,...,n 
 *          n      - number of members of the sequence, n.ge.1 
 * 
 *        output     cyr,cyi are double precision 
 *          cyr,cyi- double precision vectors whose first n components 
 *                   contain real and imaginary parts for the sequence 
 *                   cy(j)=i(fnu+j-1,z)  or 
 *                   cy(j)=i(fnu+j-1,z)*exp(-abs(x))  j=1,...,n 
 *                   depending on kode, x=real(z) 
 *          nz     - number of components set to zero due to underflow, 
 *                   nz= 0   , normal return 
 *                   nz.gt.0 , last nz components of cy set to zero 
 *                             to underflow, cy(j)=cmplx(0.0d0,0.0d0) 
 *                             j = n-nz+1,...,n 
 *          ierr   - error flag 
 *                   ierr=0, normal return - computation completed 
 *                   ierr=1, input error   - no computation 
 *                   ierr=2, overflow      - no computation, real(z) too 
 *                           large on kode=1 
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
 *        the computation is carried out by the power series for 
 *        small cabs(z), the asymptotic expansion for large cabs(z), 
 *        the miller algorithm normalized by the wronskian and a 
 *        neumann series for imtermediate magnitudes, and the 
 *        uniform asymptotic expansions for i(fnu,z) and j(fnu,z) 
 *        for large orders. backward recurrence is used to generate 
 *        sequences or reduce orders when necessary. 
 * 
 *        the calculations above are done in the right half plane and 
 *        continued into the left half plane by the formula 
 * 
 *        i(fnu,z*exp(m*pi)) = exp(m*pi*fnu)*i(fnu,z)  real(z).gt.0.0 
 *                      m = +i or -i,  i**2=-1 
 * 
 *        for negative orders,the formula 
 * 
 *             i(-fnu,z) = i(fnu,z) + (2/pi)*sin(pi*fnu)*k(fnu,z) 
 * 
 *        can be used. however,for large orders close to ints, the 
 *        the function changes radically. when fnu is a large positive 
 *        int,the magnitude of i(-fnu,z)=i(fnu,z) is a large 
 *        negative power of ten. but when fnu is not an int, 
 *        k(fnu,z) dominates in magnitude with a large positive power of 
 *        ten and the most that the second term can be reduced is by 
 *        unit roundoff from the coefficient. thus, wide changes can 
 *        occur within unit roundoff of a large int for fnu. here, 
 *        large means fnu.gt.cabs(z). 
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

int pnl_zbesi (double *zr, double *zi, double *fnu, int *kode,const int *n,
		double *cyr, double *cyi, int *nz, int *ierr)
{
  static const double pi = 3.14159265358979324;
  static const double coner = 1.;
  static const double conei = 0.;

  int i__1, i__2;
  double d__1, d__2;

  /* Local variables */
  double alim, elim, atol, fnul, rtol;
  int i__, k;
  double ascle;
  double csgni, csgnr;
  int k1;
  int k2;
  double aa, bb, fn, az;
  int nn;
  double rl, dig, arg, r1m5;
  int inu;
  double tol, sti, zni, str, znr;

  /* Parameter adjustments */
  --cyi;
  --cyr;

  /* Function Body */
  *ierr = 0;
  *nz = 0;
  if (*fnu < 0.)
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
  /*----------------------------------------------------------------------- 
   *    SET PARAMETERS RELATED TO MACHINE CONSTANTS. 
   *    TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18. 
   *    ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT. 
   *    EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND 
   *    EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR 
   *    UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE. 
   *    RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z. 
   *    DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG). 
   *    FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU. 
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
  rl = dig * 1.2 + 3.;
  fnul = (dig - 3.) * 6. + 10.;
  /*----------------------------------------------------------------------------- 
   *    TEST FOR PROPER RANGE 
   *----------------------------------------------------------------------- 
   */
  az = amos_azabs (zr, zi);
  fn = *fnu + (double) (*n - 1);
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
  znr = *zr;
  zni = *zi;
  csgnr = coner;
  csgni = conei;
  if (*zr >= 0.)
    {
      goto L40;
    }
  znr = -(*zr);
  zni = -(*zi);
  /*----------------------------------------------------------------------- 
   *    CALCULATE CSGN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE 
   *    WHEN FNU IS LARGE 
   *----------------------------------------------------------------------- 
   */
  inu = (int) (*fnu);
  arg = (*fnu - (double) inu) * pi;
  if (*zi < 0.)
    {
      arg = -arg;
    }
  csgnr = cos (arg);
  csgni = sin (arg);
  if (inu % 2 == 0)
    {
      goto L40;
    }
  csgnr = -csgnr;
  csgni = -csgni;
 L40:
  amos_zbinu (&znr, &zni, fnu, kode, n, &cyr[1], &cyi[1], nz, &rl, &fnul,
	      &tol, &elim, &alim);
  if (*nz < 0)
    {
      goto L120;
    }
  if (*zr >= 0.)
    {
      return 0;
    }
  /*----------------------------------------------------------------------- 
   *    ANALYTIC CONTINUATION TO THE LEFT HALF PLANE 
   *----------------------------------------------------------------------- 
   */
  nn = *n - *nz;
  if (nn == 0)
    {
      return 0;
    }
  rtol = 1. / tol;
  ascle = pnl_d1mach (1) * rtol * 1e3;
  i__1 = nn;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /*      STR = CYR(I)*CSGNR - CYI(I)*CSGNI 
       *      CYI(I) = CYR(I)*CSGNI + CYI(I)*CSGNR 
       *      CYR(I) = STR 
       */
      aa = cyr[i__];
      bb = cyi[i__];
      atol = 1.;
      /*Computing MAX 
       */
      d__1 = fabs (aa), d__2 = fabs (bb);
      if (MAX (d__1, d__2) > ascle)
	{
	  goto L55;
	}
      aa *= rtol;
      bb *= rtol;
      atol = tol;
    L55:
      str = aa * csgnr - bb * csgni;
      sti = aa * csgni + bb * csgnr;
      cyr[i__] = str * atol;
      cyi[i__] = sti * atol;
      csgnr = -csgnr;
      csgni = -csgni;
      /* L50: */
    }
  return 0;
 L120:
  if (*nz == -2)
    {
      goto L130;
    }
  *nz = 0;
  *ierr = 2;
  return 0;
 L130:
  *nz = 0;
  *ierr = 5;
  return 0;
 L260:
  *nz = 0;
  *ierr = 4;
  return 0;
}				/* zbesi_ */
