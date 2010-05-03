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
static const int c__2 = 2;

/*
 * y-bessel function,bessel function of complex argument, 
 *            bessel function of second kind 
 ****author  amos, donald e., sandia national laboratories 
 ****purpose  to compute the y-bessel function of a complex argument 
 ****description 
 * 
 *                     ***a double precision routine*** 
 * 
 *        on kode=1, cbesy computes an n member sequence of complex 
 *        bessel functions cy(i)=y(fnu+i-1,z) for real, nonnegative 
 *        orders fnu+i-1, i=1,...,n and complex z in the cut plane 
 *        -pi.lt.arg(z).le.pi. on kode=2, cbesy returns the scaled 
 *        functions 
 * 
 *        cy(i)=exp(-abs(y))*y(fnu+i-1,z)   i = 1,...,n , y=aimag(z) 
 * 
 *        which remove the exponential growth in both the upper and 
 *        lower half planes for z to infinity. definitions and notation 
 *        are found in the nbs handbook of mathematical functions 
 *        (ref. 1). 
 * 
 *        input      zr,zi,fnu are double precision 
 *          zr,zi  - z=cmplx(zr,zi), z.ne.cmplx(0.0d0,0.0d0), 
 *                   -pi.lt.arg(z).le.pi 
 *          fnu    - order of initial y function, fnu.ge.0.0d0 
 *          kode   - a parameter to indicate the scaling option 
 *                   kode= 1  returns 
 *                            cy(i)=y(fnu+i-1,z), i=1,...,n 
 *                       = 2  returns 
 *                            cy(i)=y(fnu+i-1,z)*exp(-abs(y)), i=1,...,n 
 *                            where y=aimag(z) 
 *          n      - number of members of the sequence, n.ge.1 
 *          cwrkr, - double precision work vectors of dimension at 
 *          cwrki    at least n 
 * 
 *        output     cyr,cyi are double precision 
 *          cyr,cyi- double precision vectors whose first n components 
 *                   contain real and imaginary parts for the sequence 
 *                   cy(i)=y(fnu+i-1,z)  or 
 *                   cy(i)=y(fnu+i-1,z)*exp(-abs(y))  i=1,...,n 
 *                   depending on kode. 
 *          nz     - nz=0 , a normal return 
 *                   nz.gt.0 , nz components of cy set to zero due to 
 *                   underflow (generally on kode=2) 
 *          ierr   - error flag 
 *                   ierr=0, normal return - computation completed 
 *                   ierr=1, input error   - no computation 
 *                   ierr=2, overflow      - no computation, fnu is 
 *                           too large or cabs(z) is too small or both 
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
 *        the computation is carried out by the formula 
 * 
 *        y(fnu,z)=0.5*(h(1,fnu,z)-h(2,fnu,z))/i 
 * 
 *        where i**2 = -1 and the hankel bessel functions h(1,fnu,z) 
 *        and h(2,fnu,z) are calculated in cbesh. 
 * 
 *        for negative orders,the formula 
 * 
 *             y(-fnu,z) = y(fnu,z)*cos(pi*fnu) + j(fnu,z)*sin(pi*fnu) 
 * 
 *        can be used. however,for large orders close to half odd 
 *        ints the function changes radically. when fnu is a large 
 *        positive half odd int,the magnitude of y(-fnu,z)=j(fnu,z)* 
 *        sin(pi*fnu) is a large negative power of ten. but when fnu is 
 *        not a half odd int, y(fnu,z) dominates in magnitude with a 
 *        large positive power of ten and the most that the second term 
 *        can be reduced is by unit roundoff from the coefficient. thus, 
 *        wide changes can occur within unit roundoff of a large half 
 *        odd int. here, large means fnu.gt.cabs(z). 
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


int pnl_zbesy (double *zr, double *zi, double *fnu, int *kode,const int *n,
		double *cyr, double *cyi, int *nz, double *cwrkr, double *cwrki,
		int *ierr)
{
  int i__1, i__2;
  double d__1, d__2;

  double hcii, elim, atol, rtol;
  int i__, k;
  double ascle;
  int k1, k2;
  double aa, bb, ey, c1i, c2i, c1r, c2r;
  int nz1, nz2;
  double exi, r1m5, exr, sti, tay, tol, str;

  /* Parameter adjustments */
  --cwrki;
  --cwrkr;
  --cyi;
  --cyr;

  /* Function Body */
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
  hcii = .5;
  pnl_zbesh (zr, zi, fnu, kode, &c__1, n, &cyr[1], &cyi[1], &nz1, ierr);
  if (*ierr != 0 && *ierr != 3)
    {
      goto L170;
    }
  pnl_zbesh (zr, zi, fnu, kode, &c__2, n, &cwrkr[1], &cwrki[1], &nz2, ierr);
  if (*ierr != 0 && *ierr != 3)
    {
      goto L170;
    }
  *nz = MIN (nz1, nz2);
  if (*kode == 2)
    {
      goto L60;
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      str = cwrkr[i__] - cyr[i__];
      sti = cwrki[i__] - cyi[i__];
      cyr[i__] = -sti * hcii;
      cyi[i__] = str * hcii;
      /* L50: */
    }
  return 0;
 L60:
  /*Computing MAX 
   */
  d__1 = pnl_d1mach (4);
  tol = MAX (d__1, 1e-18);
  k1 = amos_i1mach (15);
  k2 = amos_i1mach (16);
  /*Computing MIN 
   */
  i__1 = abs (k1), i__2 = abs (k2);
  k = MIN (i__1, i__2);
  r1m5 = pnl_d1mach (5);
  /*----------------------------------------------------------------------- 
   *    ELIM IS THE APPROXIMATE EXPONENTIAL UNDER- AND OVERFLOW LIMIT 
   *----------------------------------------------------------------------- 
   */
  elim = ((double) k * r1m5 - 3.) * 2.303;
  exr = cos (*zr);
  exi = sin (*zr);
  ey = 0.;
  tay = (d__1 = *zi + *zi, fabs (d__1));
  if (tay < elim)
    {
      ey = exp (-tay);
    }
  if (*zi < 0.)
    {
      goto L90;
    }
  c1r = exr * ey;
  c1i = exi * ey;
  c2r = exr;
  c2i = -exi;
 L70:
  *nz = 0;
  rtol = 1. / tol;
  ascle = pnl_d1mach (1) * rtol * 1e3;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /*      STR = C1R*CYR(I) - C1I*CYI(I) 
       *      STI = C1R*CYI(I) + C1I*CYR(I) 
       *      STR = -STR + C2R*CWRKR(I) - C2I*CWRKI(I) 
       *      STI = -STI + C2R*CWRKI(I) + C2I*CWRKR(I) 
       *      CYR(I) = -STI*HCII 
       *      CYI(I) = STR*HCII 
       */
      aa = cwrkr[i__];
      bb = cwrki[i__];
      atol = 1.;
      /*Computing MAX 
       */
      d__1 = fabs (aa), d__2 = fabs (bb);
      if (MAX (d__1, d__2) > ascle)
	{
	  goto L75;
	}
      aa *= rtol;
      bb *= rtol;
      atol = tol;
    L75:
      str = (aa * c2r - bb * c2i) * atol;
      sti = (aa * c2i + bb * c2r) * atol;
      aa = cyr[i__];
      bb = cyi[i__];
      atol = 1.;
      /*Computing MAX 
       */
      d__1 = fabs (aa), d__2 = fabs (bb);
      if (MAX (d__1, d__2) > ascle)
	{
	  goto L85;
	}
      aa *= rtol;
      bb *= rtol;
      atol = tol;
    L85:
      str -= (aa * c1r - bb * c1i) * atol;
      sti -= (aa * c1i + bb * c1r) * atol;
      cyr[i__] = -sti * hcii;
      cyi[i__] = str * hcii;
      if (str == 0. && sti == 0. && ey == 0.)
	{
	  ++(*nz);
	}
      /* L80: */
    }
  return 0;
 L90:
  c1r = exr;
  c1i = exi;
  c2r = exr * ey;
  c2i = -exi * ey;
  goto L70;
 L170:
  *nz = 0;
  return 0;
}				/* zbesy_ */
