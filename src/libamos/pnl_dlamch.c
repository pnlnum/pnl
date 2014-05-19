#include "stdio.h"

#ifndef _MSC_VER
#include <strings.h>            /* usefull for the function strncasecmp */
#else
/* _strnicmp is used instead of strncasecmp
   because the later doesn't exist in Windows' standard C libraries */
#define strncasecmp _strnicmp   
#endif
#include "pnl/pnl_mathtools.h"


static int pnl_lsame (const char *ca,const  char *cb);
static int pnl_dlamc1 (int *beta, int *t, int *rnd, int *ieee1);
static int pnl_dlamc2 (int *beta, int *t, int *rnd, double *eps,int *emin, double *rmin,int *emax, double *rmax);
double pnl_dlamc3 (double *a, double *b);
static  int pnl_dlamc4 (int *emin, double *start, int *base);
static  int pnl_dlamc5 (int *beta, int *p, int *emin, int *ieee, int *emax, double *rmax);


/*--------------------------------------------------------------------
 *     LAPACK auxiliary routine (version 2.0) -- 
 *     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., 
 *     Courant Institute, Argonne National Lab, and Rice University 
 *     October 31, 1992 
 *  DLAMCH determines double precision machine parameters. 
 *  CMACH   (input) CHARACTER*1 
 *          Specifies the value to be returned by DLAMCH: 
 *          = 'E' or 'e',   DLAMCH := eps 
 *          = 'S' or 's ,   DLAMCH := sfmin 
 *          = 'B' or 'b',   DLAMCH := base 
 *          = 'P' or 'p',   DLAMCH := eps*base 
 *          = 'N' or 'n',   DLAMCH := t 
 *          = 'R' or 'r',   DLAMCH := rnd 
 *          = 'M' or 'm',   DLAMCH := emin 
 *          = 'U' or 'u',   DLAMCH := rmin 
 *          = 'L' or 'l',   DLAMCH := emax 
 *          = 'O' or 'o',   DLAMCH := rmax 
 *          where 
 *          eps   = relative machine precision 
 *          sfmin = safe minimum, such that 1/sfmin does not overflow 
 *          base  = base of the machine 
 *          prec  = eps*base 
 *          t     = number of (base) digits in the mantissa 
 *          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise 
 *          emin  = minimum exponent before (gradual) underflow 
 *          rmin  = underflow threshold - base**(emin-1) 
 *          emax  = largest exponent before overflow 
 *          rmax  = overflow threshold  - (base**emax)*(1-eps) 
 *------------------------------------------------------------------*/

double pnl_dlamch (char *cmach)
{
  static int first = TRUE;
  static double base, emin, prec, emax, rmin, rmax,t, sfmin, rnd, eps,  rmach=0.0, small;
  const double zero= 0.0,one=1.0;
  int i_1, beta, imin, imax, lrnd,  it;

  if (first)
    {
      first = FALSE;
      pnl_dlamc2 (&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
      base = (double) beta;
      t = (double) it;
      if (lrnd)
	{
	  rnd = one;
	  i_1 = 1 - it;
	  eps = pnl_pow_i(base, i_1) / 2; /* pow_di */
	}
      else
	{
	  rnd = zero;
	  i_1 = 1 - it;
	  eps = pnl_pow_i(base, i_1);/* pow_di */
	}
      prec = eps * base;
      emin = (double) imin;
      emax = (double) imax;
      sfmin = rmin;
      small = one / rmax;
      if (small >= sfmin)
	{
	  /*           Use SMALL plus a bit, to avoid the possibility of rounding 
	   *           causing overflow when computing  1/sfmin. 
	   */
	  sfmin = small * (one + eps);
	}
    }

  if (pnl_lsame (cmach, "E"))
    {
      rmach = eps;
    }
  else if (pnl_lsame (cmach, "S"))
    {
      rmach = sfmin;
    }
  else if (pnl_lsame (cmach, "B"))
    {
      rmach = base;
    }
  else if (pnl_lsame (cmach, "P"))
    {
      rmach = prec;
    }
  else if (pnl_lsame (cmach, "N"))
    {
      rmach = t;
    }
  else if (pnl_lsame (cmach, "R"))
    {
      rmach = rnd;
    }
  else if (pnl_lsame (cmach, "M"))
    {
      rmach = emin;
    }
  else if (pnl_lsame (cmach, "U"))
    {
      rmach = rmin;
    }
  else if (pnl_lsame (cmach, "L"))
    {
      rmach = emax;
    }
  else if (pnl_lsame (cmach, "O"))
    {
      rmach = rmax;
    }
  return rmach;
}


/*  -- LAPACK auxiliary routine (version 2.0) -- 
 *     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., 
 *     Courant Institute, Argonne National Lab, and Rice University 
 *     October 31, 1992 
 *  DLAMC1 determines the machine parameters given by BETA, T, RND, and  IEEE1. 
 *  Arguments 
 *  ========= 
 *  BETA    (output) INT 
 *          The base of the machine. 
 *  T       (output) INT 
 *          The number of ( BETA ) digits in the mantissa. 
 *  RND     (output) INT 
 *          Specifies whether proper rounding  ( RND = .TRUE. )  or 
 *          chopping  ( RND = .FALSE. )  occurs in addition. This may not 
 *          be a reliable guide to the way in which the machine performs 
 *          its arithmetic. 
 *  IEEE1   (output) INT 
 *          Specifies whether rounding appears to be done in the IEEE 
 *          'round to nearest' style. 
 *  Further Details 
 *  =============== 
 *  The routine is based on the routine  ENVRON  by Malcolm and 
 *  incorporates suggestions by Gentleman and Marovich. See 
 *     Malcolm M. A. (1972) Algorithms to reveal properties of 
 *        floating-point arithmetic. Comms. of the ACM, 15, 949-951. 
 *     Gentleman W. M. and Marovich S. B. (1974) More on algorithms 
 *        that reveal properties of floating point arithmetic units. 
 *        Comms. of the ACM, 17, 276-277. 
 */

static int pnl_dlamc1 (int *beta, int *t, int *rnd, int *ieee1)
{
  static int first = TRUE;
  double d_1, d_2;
  /* Local variables */
  static int lrnd;
  double a, b, c_, f;
  static int lbeta;
  double savec;
  static int lieee1;
  double t1, t2;
  static int lt;
  double one, qtr;
  if (first)
    {
      first = FALSE;
      one = 1.;
      /*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA, 
       *        IEEE1, T and RND. 
       *        Throughout this routine  we use the function  DLAMC3  to ensure 
       *        that relevant values are  stored and not held in registers,  or 
       *        are not affected by optimizers. 
       *        Compute  a = 2.0**m  with the  smallest positive int m such 
       *        that 
       *           fl( a + 1.0 ) = a. 
       */
      a = 1.;
      c_ = 1.;
      /* +       WHILE( C.EQ.ONE )LOOP */
    L10:
      if (c_ == one)
	{
	  a *= 2;
	  c_ = pnl_dlamc3 (&a, &one);
	  d_1 = -a;
	  c_ = pnl_dlamc3 (&c_, &d_1);
	  goto L10;
	}
      /* +       END WHILE 
       *        Now compute  b = 2.0**m  with the smallest positive int m 
       *        such that 
       *           fl( a + b ) .gt. a. 
       */
      b = 1.;
      c_ = pnl_dlamc3 (&a, &b);

      /* +       WHILE( C.EQ.A )LOOP */
    L20:
      if (c_ == a)
	{
	  b *= 2;
	  c_ = pnl_dlamc3 (&a, &b);
	  goto L20;
	}
      /* +       END WHILE 
       *        Now compute the base.  a and c  are neighbouring floating point 
       *        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so 
       *        their difference is beta. Adding 0.25 to c is to ensure that it 
       *        is truncated to beta and not ( beta - 1 ). 
       */
      qtr = one / 4;
      savec = c_;
      d_1 = -a;
      c_ = pnl_dlamc3 (&c_, &d_1);
      lbeta = (int) (c_ + qtr);
      /*        Now determine whether rounding or chopping occurs,  by adding a 
       *        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a. 
       */
      b = (double) lbeta;
      d_1 = b / 2;
      d_2 = -b / 100;
      f = pnl_dlamc3 (&d_1, &d_2);
      c_ = pnl_dlamc3 (&f, &a);
      if (c_ == a)
	{
	  lrnd = TRUE;
	}
      else
	{
	  lrnd = FALSE;
	}
      d_1 = b / 2;
      d_2 = b / 100;
      f = pnl_dlamc3 (&d_1, &d_2);
      c_ = pnl_dlamc3 (&f, &a);
      if (lrnd && c_ == a)
	{
	  lrnd = FALSE;
	}
      /*        Try and decide whether rounding is done in the  IEEE  'round to 
       *        nearest' style. B/2 is half a unit in the last place of the two 
       *        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit 
       *        zero, and SAVEC is odd. Thus adding B/2 to A should not  change 
       *        A, but adding B/2 to SAVEC should change SAVEC. 
       */
      d_1 = b / 2;
      t1 = pnl_dlamc3 (&d_1, &a);
      d_1 = b / 2;
      t2 = pnl_dlamc3 (&d_1, &savec);
      lieee1 = t1 == a && t2 > savec && lrnd;
      /*        Now find  the  mantissa, t.  It should  be the  int part of 
       *        log to the base beta of a,  however it is safer to determine  t 
       *        by powering.  So we find t as the smallest positive int for 
       *        which 
       *           fl( beta**t + 1.0 ) = 1.0. 
       */
      lt = 0;
      a = 1.;
      c_ = 1.;
      /* +       WHILE( C.EQ.ONE )LOOP */
    L30:
      if (c_ == one)
	{
	  ++lt;
	  a *= lbeta;
	  c_ = pnl_dlamc3 (&a, &one);
	  d_1 = -a;
	  c_ = pnl_dlamc3 (&c_, &d_1);
	  goto L30;
	}
      /* +       END WHILE 
       */
    }

  *beta = lbeta;
  *t = lt;
  *rnd = lrnd;
  *ieee1 = lieee1;
  return 0;
}


/*  -- LAPACK auxiliary routine (version 2.0) -- 
 *     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., 
 *     Courant Institute, Argonne National Lab, and Rice University 
 *     October 31, 1992 
 *  Purpose 
 *  ======= 
 *  DLAMC2 determines the machine parameters specified in its argument 
 *  list. 
 *  Arguments 
 *  ========= 
 *  BETA    (output) INT 
 *          The base of the machine. 
 *  T       (output) INT 
 *          The number of ( BETA ) digits in the mantissa. 
 *  RND     (output) INT 
 *          Specifies whether proper rounding  ( RND = .TRUE. )  or 
 *          chopping  ( RND = .FALSE. )  occurs in addition. This may not 
 *          be a reliable guide to the way in which the machine performs 
 *          its arithmetic. 
 *  EPS     (output) DOUBLE PRECISION 
 *          The smallest positive number such that 
 *             fl( 1.0 - EPS ) .LT. 1.0, 
 *          where fl denotes the computed value. 
 *  EMIN    (output) INT 
 *          The minimum exponent before (gradual) underflow occurs. 
 *  RMIN    (output) DOUBLE PRECISION 
 *          The smallest normalized number for the machine, given by 
 *          BASE**( EMIN - 1 ), where  BASE  is the floating point value 
 *          of BETA. 
 *  EMAX    (output) INT 
 *          The maximum exponent before overflow occurs. 
 *  RMAX    (output) DOUBLE PRECISION 
 *          The largest positive number for the machine, given by 
 *          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point 
 *          value of BETA. 
 *  Further Details 
 *  =============== 
 *  The computation of  EPS  is based on a routine PARANOIA by 
 *  W. Kahan of the University of California at Berkeley.
 */

static int pnl_dlamc2 (int *beta, int *t, int *rnd, double *eps, 
		       int *emin, double *rmin,
		       int *emax, double *rmax)
{
  static int first = TRUE;
  static int iwarn = FALSE;
  /* System generated locals */
  int i_1;
  double d_1, d_2, d_3, d_4, d_5;

  int ieee;
  double half;
  static int lrnd;
  static double leps;
  double zero, a, b, c_;
  int i_;
  static int lbeta;
  double rbase;
  static int lemin, lemax;
  int gnmin;
  double small;
  int gpmin;
  double third;
  static double lrmin, lrmax;
  double sixth;
  int lieee1;
  static int lt;
  int ngnmin, ngpmin;
  double one, two;


  if (first)
    {
      first = FALSE;
      zero = 0.;
      one = 1.;
      two = 2.;
      /*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of 
       *        BETA, T, RND, EPS, EMIN and RMIN. 
       *        Throughout this routine  we use the function  DLAMC3  to ensure 
       *        that relevant values are stored  and not held in registers,  or 
       *        are not affected by optimizers. 
       *        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1. 
       */

      pnl_dlamc1 (&lbeta, &lt, &lrnd, &lieee1);
      
      /*        Start to find EPS. */

      b = (double) lbeta;
      i_1 = -lt;
      a = pnl_pow_i(b, i_1); /* pow_di */
      leps = a;
      /*        Try some tricks to see whether or not this is the correct  EPS. */
      b = two / 3;
      half = one / 2;
      d_1 = -half;
      sixth = pnl_dlamc3 (&b, &d_1);
      third = pnl_dlamc3 (&sixth, &sixth);
      d_1 = -half;
      b = pnl_dlamc3 (&third, &d_1);
      b = pnl_dlamc3 (&b, &sixth);
      b = fabs (b);
      if (b < leps)
	{
	  b = leps;
	}

      leps = 1.;
      /* +       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP */
    L10:
      if (leps > b && b > zero)
	{
	  leps = b;
	  d_1 = half * leps;
	  /* Computing 5th power */
	  d_3 = two, d_4 = d_3, d_3 *= d_3;
	  /* Computing 2nd power */
	  d_5 = leps;
	  d_2 = d_4 * (d_3 * d_3) * (d_5 * d_5);
	  c_ = pnl_dlamc3 (&d_1, &d_2);
	  d_1 = -c_;
	  c_ = pnl_dlamc3 (&half, &d_1);
	  b = pnl_dlamc3 (&half, &c_);
	  d_1 = -b;
	  c_ = pnl_dlamc3 (&half, &d_1);
	  b = pnl_dlamc3 (&half, &c_);
	  goto L10;
	}
      /* +       END WHILE */
      if (a < leps)
	{
	  leps = a;
	}
      /*
       *        Computation of EPS complete. 
       *        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)). 
       *        Keep dividing  A by BETA until (gradual) underflow occurs. This 
       *        is detected when we cannot recover the previous A. 
       */

      rbase = one / lbeta;
      small = one;
      for (i_ = 1; i_ <= 3; ++i_)
	{
	  d_1 = small * rbase;
	  small = pnl_dlamc3 (&d_1, &zero);
	  /* L20: */
	}
      a = pnl_dlamc3 (&one, &small);
      pnl_dlamc4 (&ngpmin, &one, &lbeta);
      d_1 = -one;
      pnl_dlamc4 (&ngnmin, &d_1, &lbeta);
      pnl_dlamc4 (&gpmin, &a, &lbeta);
      d_1 = -a;
      pnl_dlamc4 (&gnmin, &d_1, &lbeta);
      ieee = FALSE;

      if (ngpmin == ngnmin && gpmin == gnmin)
	{
	  if (ngpmin == gpmin)
	    {
	      lemin = ngpmin;
	      /*            ( Non twos-complement machines, no gradual underflow; 
	       *              e.g.,  VAX ) 
	       */
	    }
	  else if (gpmin - ngpmin == 3)
	    {
	      lemin = ngpmin - 1 + lt;
	      ieee = TRUE;
	      /*            ( Non twos-complement machines, with gradual underflow; 
	       *              e.g., IEEE standard followers ) 
	       */
	    }
	  else
	    {
	      lemin = MIN (ngpmin, gpmin);
	      /*            ( A guess; no known machine ) */
	      iwarn = TRUE;
	    }

	}
      else if (ngpmin == gpmin && ngnmin == gnmin)
	{
	  if ((i_1 = ngpmin - ngnmin, fabs (i_1)) == 1)
	    {
	      lemin = MAX (ngpmin, ngnmin);
	      /*            ( Twos-complement machines, no gradual underflow; 
	       *              e.g., CYBER 205 ) 
	       */
	    }
	  else
	    {
	      lemin = MIN (ngpmin, ngnmin);
	      /*            ( A guess; no known machine ) */
	      iwarn = TRUE;
	    }

	}
      else if ((i_1 = ngpmin - ngnmin, fabs (i_1)) == 1 && gpmin == gnmin)
	{
	  if (gpmin - MIN (ngpmin, ngnmin) == 3)
	    {
	      lemin = MAX (ngpmin, ngnmin) - 1 + lt;
	      /*            ( Twos-complement machines with gradual underflow; 
	       *              no known machine ) 
	       */
	    }
	  else
	    {
	      lemin = MIN (ngpmin, ngnmin);
	      /*            ( A guess; no known machine ) */
	      iwarn = TRUE;
	    }

	}
      else
	{
	  /* Computing MIN */
	  i_1 = MIN (ngpmin, ngnmin), i_1 = MIN (i_1, gpmin);
	  lemin = MIN (i_1, gnmin);
	  /*         ( A guess; no known machine ) */
	  iwarn = TRUE;
	}
      /* ** */
      /* Comment out this if block if EMIN is ok */
      if (iwarn)
	{
	  printf("WARNING. The value EMIN may be incorrect: Emin=%d\n",lemin);
	  printf(" If, after ipnlection, the value EMIN looks\n");
	  printf(" acceptable please comment out\n");
	  printf(" the IF block as marked within the code of routine\n");
	  printf(" dlamc2");
	}
      /* ** 
       *        Assume IEEE arithmetic if we found denormalised  numbers above, 
       *        or if arithmetic seems to round in the  IEEE style,  determined 
       *        in routine DLAMC1. A true IEEE machine should have both  things 
       *        true; however, faulty machines may have one or the other. 
       */
      ieee = ieee || lieee1;
      /*        Compute  RMIN by successive division by  BETA. We could compute 
       *        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during 
       *        this computation. 
       */
      lrmin = 1.;
      i_1 = 1 - lemin;
      for (i_ = 1; i_ <= i_1; ++i_)
	{
	  d_1 = lrmin * rbase;
	  lrmin = pnl_dlamc3 (&d_1, &zero);
	  /* L30: */
	}
      /*        Finally, call DLAMC5 to compute EMAX and RMAX. */
      pnl_dlamc5 (&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
    }

  *beta = lbeta;
  *t = lt;
  *rnd = lrnd;
  *eps = leps;
  *emin = lemin;
  *rmin = lrmin;
  *emax = lemax;
  *rmax = lrmax;

  return 0;
}


/*  -- LAPACK auxiliary routine (version 2.0) -- 
 *     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., 
 *     Courant Institute, Argonne National Lab, and Rice University 
 *     October 31, 1992 
 *  Purpose 
 *  ======= 
 *  DLAMC3  is intended to force  A  and  B  to be stored prior to doing 
 *  the addition of  A  and  B ,  for use in situations where optimizers 
 *  might hold one of these in a register. 
 *  Arguments 
 *  ========= 
 *  A, B    (input) DOUBLE PRECISION 
 *          The values A and B. 
 */


double pnl_dlamc3 (double *a, double *b)
{
  double ret_val;
  ret_val = *a + *b;
  return ret_val;
}

/*  -- LAPACK auxiliary routine (version 2.0) -- 
 *     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., 
 *     Courant Institute, Argonne National Lab, and Rice University 
 *     October 31, 1992 
 *     .. Scalar Arguments .. 
 *     .. 
 *  Purpose 
 *  ======= 
 *  DLAMC4 is a service routine for DLAMC2. 
 *  Arguments 
 *  ========= 
 *  EMIN    (output) EMIN 
 *          The minimum exponent before (gradual) underflow, computed by 
 *          setting A = START and dividing by BASE until the previous A 
 *          can not be recovered. 
 *  START   (input) DOUBLE PRECISION 
 *          The starting point for determining EMIN. 
 *  BASE    (input) INT 
 *          The base of the machine. 
 */

static  int pnl_dlamc4 (int *emin, double *start, int *base)
{
  int i_1;
  double d_1;

  /* Local variables */
  double zero, a;
  int i_;
  double rbase, b1, b2, c1, c2, d1, d2;
  double one;


  a = *start;
  one = 1.;
  rbase = one / *base;
  zero = 0.;
  *emin = 1;
  d_1 = a * rbase;
  b1 = pnl_dlamc3 (&d_1, &zero);
  c1 = a;
  c2 = a;
  d1 = a;
  d2 = a;
  /* +    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND. */
  /*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP */
 L10:
  if (c1 == a && c2 == a && d1 == a && d2 == a)
    {
      --(*emin);
      a = b1;
      d_1 = a / *base;
      b1 = pnl_dlamc3 (&d_1, &zero);
      d_1 = b1 * *base;
      c1 = pnl_dlamc3 (&d_1, &zero);
      d1 = zero;
      i_1 = *base;
      for (i_ = 1; i_ <= i_1; ++i_)
	{
	  d1 += b1;
	  /* L20: */
	}
      d_1 = a * rbase;
      b2 = pnl_dlamc3 (&d_1, &zero);
      d_1 = b2 / rbase;
      c2 = pnl_dlamc3 (&d_1, &zero);
      d2 = zero;
      i_1 = *base;
      for (i_ = 1; i_ <= i_1; ++i_)
	{
	  d2 += b2;
	  /* L30: */
	}
      goto L10;
    }
  /* +    END WHILE */

  return 0;
}


/*  -- LAPACK auxiliary routine (version 2.0) -- 
 *     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., 
 *     Courant Institute, Argonne National Lab, and Rice University 
 *     October 31, 1992 
 *  Purpose 
 *  ======= 
 *  DLAMC5 attempts to compute RMAX, the largest machine floating-point 
 *  number, without overflow.  It assumes that EMAX + fabs(EMIN) sum 
 *  approximately to a power of 2.  It will fail on machines where this 
 *  assumption does not hold, for example, the Cyber 205 (EMIN = -28625, 
 *  EMAX = 28718).  It will also fail if the value supplied for EMIN is 
 *  too large (i.e. too close to zero), probably with overflow. 
 *  Arguments 
 *  ========= 
 *  BETA    (input) INT 
 *          The base of floating-point arithmetic. 
 *  P       (input) INT 
 *          The number of base BETA digits in the mantissa of a 
 *          floating-point value. 
 *  EMIN    (input) INT 
 *          The minimum exponent before (gradual) underflow. 
 *  IEEE    (input) INT 
 *          A int flag specifying whether or not the arithmetic 
 *          system is thought to comply with the IEEE standard. 
 *  EMAX    (output) INT 
 *          The largest exponent before overflow 
 *  RMAX    (output) DOUBLE PRECISION 
 *          The largest machine floating-point number. 
 *     First compute LEXP and UEXP, two powers of 2 that bound 
 *     fabs(EMIN). We then assume that EMAX + fabs(EMIN) will sum 
 *     approximately to the bound that is closest to fabs(EMIN). 
 *     (EMAX is the exponent of the required number RMAX). 
 */

static  int pnl_dlamc5 (int *beta, int *p, int *emin, int *ieee, int *emax, double *rmax)
{
  /*  ZERO = 0.0D0, ONE = 1.0D0 */
  int i_1;
  double d_1;
  int lexp;
  double oldy=0.0, zero=0.0;
  int uexp, i_;
  double y, z_;
  int nbits;
  double recbas;
  int exbits, expsum;
  double one=1.0;
  int try_;


  lexp = 1;
  exbits = 1;
 L10:
  try_ = lexp << 1;
  if (try_ <= -(*emin))
    {
      lexp = try_;
      ++exbits;
      goto L10;
    }
  if (lexp == -(*emin))
    {
      uexp = lexp;
    }
  else
    {
      uexp = try_;
      ++exbits;
    }

  /*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
   *     than or equal to EMIN. EXBITS is the number of bits needed to
   *     store the exponent. 
   */

  if (uexp + *emin > -lexp - *emin)
    {
      expsum = lexp << 1;
    }
  else
    {
      expsum = uexp << 1;
    }
  /*     EXPSUM is the exponent range, approximately equal to 
   *     EMAX - EMIN + 1 . 
   */
  *emax = expsum + *emin - 1;
  nbits = exbits + 1 + *p;
  /*     NBITS is the total number of bits needed to store a 
   *     floating-point number. 
   */

  if (nbits % 2 == 1 && *beta == 2)
    {
      /*        Either there are an odd number of bits used to store a 
       *        floating-point number, which is unlikely, or some bits are 
       *        not used in the representation of numbers, which is possible, 
       *        (e.g. Cray machines) or the mantissa has an implicit bit, 
       *        (e.g. IEEE machines, Dec Vax machines), which is perhaps the 
       *        most likely. We have to assume the last alternative. 
       *        If this is true, then we need to reduce EMAX by one because 
       *        there must be some way of representing zero in an implicit-bit 
       *        system. On machines like Cray, we are reducing EMAX by one 
       *        unnecessarily. 
       */

      --(*emax);
    }

  if (*ieee)
    {
      /*        Assume we are on an IEEE machine which reserves one exponent 
       *        for infinity and NaN. 
       */
      --(*emax);
    }
  /*     Now create RMAX, the largest machine number, which should 
   *     be equal to (1.0 - BETA**(-P)) * BETA**EMAX . 
   *     First compute 1.0 - BETA**(-P), being careful that the 
   *     result is less than 1.0 . 
   */
  recbas = one / *beta;
  z_ = *beta - one;
  y = zero;
  i_1 = *p;
  for (i_ = 1; i_ <= i_1; ++i_)
    {
      z_ *= recbas;
      if (y < one)
	{
	  oldy = y;
	}
      y = pnl_dlamc3 (&y, &z_);
      /* L20: */
    }
  if (y >= one)
    {
      y = oldy;
    }
  /*     Now multiply by BETA**EMAX to get RMAX. */
  i_1 = *emax;
  for (i_ = 1; i_ <= i_1; ++i_)
    {
      d_1 = y * *beta;
      y = pnl_dlamc3 (&d_1, &zero);
      /* L30: */
    }
  *rmax = y;
  return 0;
}


/*  LSAME returns .TRUE. if CA is the same letter as CB regardless of case. */

static int pnl_lsame (const char *ca,const  char *cb)
{
  return strncasecmp(ca,cb,1)==0;
}
