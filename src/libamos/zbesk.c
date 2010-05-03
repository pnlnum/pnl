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
 * k-bessel function,complex bessel function, 
 *            modified bessel function of the second kind, 
 *            bessel function of the third kind 
 ****author  amos, donald e., sandia national laboratories 
 ****purpose  to compute k-bessel functions of complex argument 
 ****description 
 * 
 *                     ***a double precision routine*** 
 * 
 *        on kode=1, cbesk computes an n member sequence of complex 
 *        bessel functions cy(j)=k(fnu+j-1,z) for real, nonnegative 
 *        orders fnu+j-1, j=1,...,n and complex z.ne.cmplx(0.0,0.0) 
 *        in the cut plane -pi.lt.arg(z).le.pi. on kode=2, cbesk 
 *        returns the scaled k functions, 
 * 
 *        cy(j)=exp(z)*k(fnu+j-1,z) , j=1,...,n, 
 * 
 *        which remove the exponential behavior in both the left and 
 *        right half planes for z to infinity. definitions and 
 *        notation are found in the nbs handbook of mathematical 
 *        functions (ref. 1). 
 * 
 *        input      zr,zi,fnu are double precision 
 *          zr,zi  - z=cmplx(zr,zi), z.ne.cmplx(0.0d0,0.0d0), 
 *                   -pi.lt.arg(z).le.pi 
 *          fnu    - order of initial k function, fnu.ge.0.0d0 
 *          n      - number of members of the sequence, n.ge.1 
 *          kode   - a parameter to indicate the scaling option 
 *                   kode= 1  returns 
 *                            cy(i)=k(fnu+i-1,z), i=1,...,n 
 *                       = 2  returns 
 *                            cy(i)=k(fnu+i-1,z)*exp(z), i=1,...,n 
 * 
 *        output     cyr,cyi are double precision 
 *          cyr,cyi- double precision vectors whose first n components 
 *                   contain real and imaginary parts for the sequence 
 *                   cy(i)=k(fnu+i-1,z), i=1,...,n or 
 *                   cy(i)=k(fnu+i-1,z)*exp(z), i=1,...,n 
 *                   depending on kode 
 *          nz     - number of components set to zero due to underflow. 
 *                   nz= 0   , normal return 
 *                   nz.gt.0 , first nz components of cy set to zero due 
 *                             to underflow, cy(i)=cmplx(0.0d0,0.0d0), 
 *                             i=1,...,n when x.ge.0.0. when x.lt.0.0 
 *                             nz states only the number of underflows 
 *                             in the sequence. 
 * 
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
 *        equations of the reference are implemented for small orders 
 *        dnu and dnu+1.0 in the right half plane x.ge.0.0. forward 
 *        recurrence generates higher orders. k is continued to the left 
 *        half plane by the relation 
 * 
 *        k(fnu,z*exp(mp)) = exp(-mp*fnu)*k(fnu,z)-mp*i(fnu,z) 
 *        mp=mr*pi*i, mr=+1 or -1, re(z).gt.0, i**2=-1 
 * 
 *        where i(fnu,z) is the i bessel function. 
 * 
 *        for large orders, fnu.gt.fnul, the k function is computed 
 *        by means of its uniform asymptotic expansions. 
 * 
 *        for negative orders, the formula 
 * 
 *                      k(-fnu,z) = k(fnu,z) 
 * 
 *        can be used. 
 * 
 *        cbesk assumes that a significant digit sinh(x) function is 
 *        available. 
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
 *                and large order by d. e. amos, sand83-0643, may, 1983. 
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

int pnl_zbesk (double *zr, double *zi, double *fnu, int *kode,const int *n,
		double *cyr, double *cyi, int *nz, int *ierr)
{
  int i__1, i__2;
  double d__1;

  /* Local variables */
  double alim, elim, fnul;
  int k;
  int k1;
  int k2;
  double aa, bb, fn, az;
  int nn;
  double rl;
  int mr, nw;
  double dig, arg, aln, r1m5, ufl;
  int nuf;
  double tol;

  /* Parameter adjustments */
  --cyi;
  --cyr;

  /* Function Body */
  *ierr = 0;
  *nz = 0;
  if (*zi == 0. && *zr == 0.)
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
  /*----------------------------------------------------------------------------- 
   *    TEST FOR PROPER RANGE 
   *----------------------------------------------------------------------- 
   */
  az = amos_azabs (zr, zi);
  fn = *fnu + (double) (nn - 1);
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
   *    UFL = DEXP(-ELIM) 
   */
  ufl = pnl_d1mach (1) * 1e3;
  if (az < ufl)
    {
      goto L180;
    }
  if (*fnu > fnul)
    {
      goto L80;
    }
  if (fn <= 1.)
    {
      goto L60;
    }
  if (fn > 2.)
    {
      goto L50;
    }
  if (az > tol)
    {
      goto L60;
    }
  arg = az * .5;
  aln = -fn * log (arg);
  if (aln > elim)
    {
      goto L180;
    }
  goto L60;
 L50:
  amos_zuoik (zr, zi, fnu, kode, &c__2, &nn, &cyr[1], &cyi[1], &nuf, &tol,
	      &elim, &alim);
  if (nuf < 0)
    {
      goto L180;
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
      goto L100;
    }
 L60:
  if (*zr < 0.)
    {
      goto L70;
    }
  /*----------------------------------------------------------------------- 
   *    RIGHT HALF PLANE COMPUTATION, REAL(Z).GE.0. 
   *----------------------------------------------------------------------- 
   */
  amos_zbknu (zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, &tol, &elim,
	      &alim);
  if (nw < 0)
    {
      goto L200;
    }
  *nz = nw;
  return 0;
  /*----------------------------------------------------------------------- 
   *    LEFT HALF PLANE COMPUTATION 
   *    PI/2.LT.ARG(Z).LE.PI AND -PI.LT.ARG(Z).LT.-PI/2. 
   *----------------------------------------------------------------------- 
   */
 L70:
  if (*nz != 0)
    {
      goto L180;
    }
  mr = 1;
  if (*zi < 0.)
    {
      mr = -1;
    }
  amos_zacon (zr, zi, fnu, kode, &mr, &nn, &cyr[1], &cyi[1], &nw, &rl, &fnul,
	      &tol, &elim, &alim);
  if (nw < 0)
    {
      goto L200;
    }
  *nz = nw;
  return 0;
  /*----------------------------------------------------------------------- 
   *    UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL 
   *----------------------------------------------------------------------- 
   */
 L80:
  mr = 0;
  if (*zr >= 0.)
    {
      goto L90;
    }
  mr = 1;
  if (*zi < 0.)
    {
      mr = -1;
    }
 L90:
  amos_zbunk (zr, zi, fnu, kode, &mr, &nn, &cyr[1], &cyi[1], &nw, &tol, &elim,
	      &alim);
  if (nw < 0)
    {
      goto L200;
    }
  *nz += nw;
  return 0;
 L100:
  if (*zr < 0.)
    {
      goto L180;
    }
  return 0;
 L180:
  *nz = 0;
  *ierr = 2;
  return 0;
 L200:
  if (nw == -1)
    {
      goto L180;
    }
  *nz = 0;
  *ierr = 5;
  return 0;
 L260:
  *nz = 0;
  *ierr = 4;
  return 0;
}				/* zbesk_ */
