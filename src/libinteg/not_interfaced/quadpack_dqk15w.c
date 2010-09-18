/*
 * This file contains the routine dqk15w from the quadpack distribution
 * available on www.netlib.org
 *
 * This routine has been translated from Fortran using f2c and then slightly
 * modified by Jérôme Lelong for integration into Pnl. All Quadpack names have
 * been prefixed by pnl to avoid possible clash in names. Calls to d1mach have
 * been replaced by calls to dlamch coming from Lapack
 */

#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"
typedef double(*D_weight)();

int pnl_dqk15w(PnlFunc * f, D_weight w, double *p1, double *p2, 
	double *p3, double *p4, int *kp, double *a, 
	double *b, double *result, double *abserr, double *
	resabs, double *resasc)
{
    /* Initialized data */

    double xgk[8] = { .9914553711208126,.9491079123427585,
	    .8648644233597691,.7415311855993944,.5860872354676911,
	    .4058451513773972,.2077849550078985,0. };
    double wgk[8] = { .02293532201052922,.06309209262997855,
	    .1047900103222502,.1406532597155259,.1690047266392679,
	    .1903505780647854,.2044329400752989,.2094821410847278 };
    double wg[4] = { .1294849661688697,.2797053914892767,
	    .3818300505051889,.4179591836734694 };

    double c_b7 = 1.5;

    /* System generated locals */
    double d__1, d__2, d__3;

    /* Local variables */
    int j;
    double fc, fv1[7], fv2[7];
    int jtw;
    double absc, resg, resk, fsum, absc1, absc2, fval1, fval2;
    int jtwm1;
    double hlgth, centr, reskh, uflow;
    
    double epmach, dhlgth;

/* ***begin prologue  dqk15w */
/* ***date written   810101   (yymmdd) */
/* ***revision date  830518   (mmddyy) */
/* ***category no.  h2a2a2 */
/* ***keywords  15-point gauss-kronrod rules */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  to compute i = integral of f*w over (a,b), with error */
/*                           estimate */
/*                       j = integral of ABS(f*w) over (a,b) */
/* ***description */

/*           integration rules */
/*           standard fortran subroutine */
/*           double precision version */

/*           parameters */
/*             on entry */
/*              f      - double precision */
/*                       function subprogram defining the integrand */
/*                       function f(x). the actual name for f needs to be */
/*                       declared e x t e r n a l in the driver program. */

/*              w      - double precision */
/*                       function subprogram defining the integrand */
/*                       weight function w(x). the actual name for w */
/*                       needs to be declared e x t e r n a l in the */
/*                       calling program. */

/*              p1, p2, p3, p4 - double precision */
/*                       parameters in the weight function */

/*              kp     - int */
/*                       key for indicating the type of weight function */

/*              a      - double precision */
/*                       lower limit of integration */

/*              b      - double precision */
/*                       upper limit of integration */

/*            on return */
/*              result - double precision */
/*                       approximation to the integral i */
/*                       result is computed by applying the 15-point */
/*                       kronrod rule (resk) obtained by optimal addition */
/*                       of abscissae to the 7-point gauss rule (resg). */

/*              abserr - double precision */
/*                       estimate of the modulus of the absolute error, */
/*                       which should equal or exceed ABS(i-result) */

/*              resabs - double precision */
/*                       approximation to the integral of ABS(f) */

/*              resasc - double precision */
/*                       approximation to the integral of ABS(f-i/(b-a)) */


/* ***references  (none) */
/* ***end prologue  dqk15w */



/*           the abscissae and weights are given for the interval (-1,1). */
/*           because of symmetry only the positive abscissae and their */
/*           corresponding weights are given. */

/*           xgk    - abscissae of the 15-point gauss-kronrod rule */
/*                    xgk(2), xgk(4), ... abscissae of the 7-point */
/*                    gauss rule */
/*                    xgk(1), xgk(3), ... abscissae which are optimally */
/*                    added to the 7-point gauss rule */

/*           wgk    - weights of the 15-point gauss-kronrod rule */

/*           wg     - weights of the 7-point gauss rule */





/*           list of major variables */
/*           ----------------------- */

/*           centr  - mid point of the interval */
/*           hlgth  - half-length of the interval */
/*           absc*  - abscissa */
/*           fval*  - function value */
/*           resg   - result of the 7-point gauss formula */
/*           resk   - result of the 15-point kronrod formula */
/*           reskh  - approximation to the mean value of f*w over (a,b), */
/*                    i.e. to i/(b-a) */

/*           machine dependent constants */
/*           --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */

/* ***first executable statement  dqk15w */
    epmach = pnl_d1mach(4);
    uflow = pnl_d1mach(1);

    centr = (*a + *b) * .5;
    hlgth = (*b - *a) * .5;
    dhlgth = ABS(hlgth);

/*           compute the 15-point kronrod approximation to the */
/*           integral, and estimate the error. */

    fc = PNL_EVAL_FUNC(f,centr) * (*w)(&centr, p1, p2, p3, p4, kp);
    resg = wg[3] * fc;
    resk = wgk[7] * fc;
    *resabs = ABS(resk);
    for (j = 1; j <= 3; ++j) {
	jtw = j << 1;
	absc = hlgth * xgk[jtw - 1];
	absc1 = centr - absc;
	absc2 = centr + absc;
	fval1 = PNL_EVAL_FUNC(f,absc1) * (*w)(&absc1, p1, p2, p3, p4, kp);
	fval2 = PNL_EVAL_FUNC(f,absc2) * (*w)(&absc2, p1, p2, p3, p4, kp);
	fv1[jtw - 1] = fval1;
	fv2[jtw - 1] = fval2;
	fsum = fval1 + fval2;
	resg += wg[j - 1] * fsum;
	resk += wgk[jtw - 1] * fsum;
	*resabs += wgk[jtw - 1] * (ABS(fval1) + ABS(fval2));
/* L10: */
    }
    for (j = 1; j <= 4; ++j) {
	jtwm1 = (j << 1) - 1;
	absc = hlgth * xgk[jtwm1 - 1];
	absc1 = centr - absc;
	absc2 = centr + absc;
	fval1 = PNL_EVAL_FUNC(f,absc1) * (*w)(&absc1, p1, p2, p3, p4, kp);
	fval2 = PNL_EVAL_FUNC(f,absc2) * (*w)(&absc2, p1, p2, p3, p4, kp);
	fv1[jtwm1 - 1] = fval1;
	fv2[jtwm1 - 1] = fval2;
	fsum = fval1 + fval2;
	resk += wgk[jtwm1 - 1] * fsum;
	*resabs += wgk[jtwm1 - 1] * (ABS(fval1) + ABS(fval2));
/* L15: */
    }
    reskh = resk * .5;
    *resasc = wgk[7] * (d__1 = fc - reskh, ABS(d__1));
    for (j = 1; j <= 7; ++j) {
	*resasc += wgk[j - 1] * ((d__1 = fv1[j - 1] - reskh, ABS(d__1)) + (
		d__2 = fv2[j - 1] - reskh, ABS(d__2)));
/* L20: */
    }
    *result = resk * hlgth;
    *resabs *= dhlgth;
    *resasc *= dhlgth;
    *abserr = (d__1 = (resk - resg) * hlgth, ABS(d__1));
    if (*resasc != 0. && *abserr != 0.) {
/* Computing MIN */
	d__3 = *abserr * 200. / *resasc;
	d__1 = 1., d__2 = pow(d__3, c_b7);
	*abserr = *resasc * MIN(d__1,d__2);
    }
    if (*resabs > uflow / (epmach * 50.)) {
/* Computing MAX */
	d__1 = epmach * 50. * *resabs;
	*abserr = MAX(d__1,*abserr);
    }
    return 0;
} /* pnl_dqk15w */

