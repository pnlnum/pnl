/*
 * This file contains the routine dqk21 from the quadpack distribution
 * available on www.netlib.org
 *
 * This routine has been translated from Fortran using f2c and then slightly
 * modified by Jérôme Lelong for integration into Pnl. All Quadpack names have
 * been prefixed by pnl to avoid possible clash in names. Calls to d1mach have
 * been replaced by calls to dlamch coming from Lapack
 */

#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"


int pnl_dqk21(PnlFunc * f, double *a, double *b, double *
	result, double *abserr, double *resabs, double *resasc)
{
    /* Initialized data */

    double wg[5] = { .066671344308688137593568809893332,
	    .149451349150580593145776339657697,
	    .219086362515982043995534934228163,
	    .269266719309996355091226921569469,
	    .295524224714752870173892994651338 };
    double xgk[11] = { .995657163025808080735527280689003,
	    .973906528517171720077964012084452,
	    .930157491355708226001207180059508,
	    .865063366688984510732096688423493,
	    .780817726586416897063717578345042,
	    .679409568299024406234327365114874,
	    .562757134668604683339000099272694,
	    .433395394129247190799265943165784,
	    .294392862701460198131126603103866,
	    .14887433898163121088482600112972,0. };
    double wgk[11] = { .011694638867371874278064396062192,
	    .03255816230796472747881897245939,
	    .05475589657435199603138130024458,
	    .07503967481091995276704314091619,
	    .093125454583697605535065465083366,
	    .109387158802297641899210590325805,
	    .123491976262065851077958109831074,
	    .134709217311473325928054001771707,
	    .142775938577060080797094273138717,
	    .147739104901338491374841515972068,
	    .149445554002916905664936468389821 };

    double c_b7 = 1.5;

    /* System generated locals */
    double d__1, d__2, d__3;

    /* Local variables */
    int j;
    double fc, fv1[10], fv2[10];
    int jtw;
    double absc, resg, resk, fsum, fval1, fval2;
    int jtwm1;
    double hlgth, centr, reskh, uflow;
    
    double epmach, dhlgth;

/* ***begin prologue  dqk21 */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a1a2 */
/* ***keywords  21-point gauss-kronrod rules */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  to compute i = integral of f over (a,b), with error */
/*                           estimate */
/*                       j = integral of ABS(f) over (a,b) */
/* ***description */

/*           integration rules */
/*           standard fortran subroutine */
/*           double precision version */

/*           parameters */
/*            on entry */
/*              f      - double precision */
/*                       function subprogram defining the integrand */
/*                       function f(x). the actual name for f needs to be */
/*                       declared e x t e r n a l in the driver program. */

/*              a      - double precision */
/*                       lower limit of integration */

/*              b      - double precision */
/*                       upper limit of integration */

/*            on return */
/*              result - double precision */
/*                       approximation to the integral i */
/*                       result is computed by applying the 21-point */
/*                       kronrod rule (resk) obtained by optimal addition */
/*                       of abscissae to the 10-point gauss rule (resg). */

/*              abserr - double precision */
/*                       estimate of the modulus of the absolute error, */
/*                       which should not exceed ABS(i-result) */

/*              resabs - double precision */
/*                       approximation to the integral j */

/*              resasc - double precision */
/*                       approximation to the integral of ABS(f-i/(b-a)) */
/*                       over (a,b) */

/* ***references  (none) */
/* ***end prologue  dqk21 */



/*           the abscissae and weights are given for the interval (-1,1). */
/*           because of symmetry only the positive abscissae and their */
/*           corresponding weights are given. */

/*           xgk    - abscissae of the 21-point kronrod rule */
/*                    xgk(2), xgk(4), ...  abscissae of the 10-point */
/*                    gauss rule */
/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
/*                    added to the 10-point gauss rule */

/*           wgk    - weights of the 21-point kronrod rule */

/*           wg     - weights of the 10-point gauss rule */


/* gauss quadrature weights and kronron quadrature abscissae and weights */
/* as evaluated with 80 decimal digit arithmetic by l. w. fullerton, */
/* bell labs, nov. 1981. */





/*           list of major variables */
/*           ----------------------- */

/*           centr  - mid point of the interval */
/*           hlgth  - half-length of the interval */
/*           absc   - abscissa */
/*           fval*  - function value */
/*           resg   - result of the 10-point gauss formula */
/*           resk   - result of the 21-point kronrod formula */
/*           reskh  - approximation to the mean value of f over (a,b), */
/*                    i.e. to i/(b-a) */


/*           machine dependent constants */
/*           --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */

/* ***first executable statement  dqk21 */
    epmach = pnl_d1mach(4);
    uflow = pnl_d1mach(1);

    centr = (*a + *b) * .5;
    hlgth = (*b - *a) * .5;
    dhlgth = ABS(hlgth);

/*           compute the 21-point kronrod approximation to */
/*           the integral, and estimate the absolute error. */

    resg = 0.;
    fc = PNL_EVAL_FUNC(f,centr);
    resk = wgk[10] * fc;
    *resabs = ABS(resk);
    for (j = 1; j <= 5; ++j) {
	jtw = j << 1;
	absc = hlgth * xgk[jtw - 1];
	d__1 = centr - absc;
	fval1 = PNL_EVAL_FUNC(f,d__1);
	d__1 = centr + absc;
	fval2 = PNL_EVAL_FUNC(f,d__1);
	fv1[jtw - 1] = fval1;
	fv2[jtw - 1] = fval2;
	fsum = fval1 + fval2;
	resg += wg[j - 1] * fsum;
	resk += wgk[jtw - 1] * fsum;
	*resabs += wgk[jtw - 1] * (ABS(fval1) + ABS(fval2));
/* L10: */
    }
    for (j = 1; j <= 5; ++j) {
	jtwm1 = (j << 1) - 1;
	absc = hlgth * xgk[jtwm1 - 1];
	d__1 = centr - absc;
	fval1 = PNL_EVAL_FUNC(f,d__1);
	d__1 = centr + absc;
	fval2 = PNL_EVAL_FUNC(f,d__1);
	fv1[jtwm1 - 1] = fval1;
	fv2[jtwm1 - 1] = fval2;
	fsum = fval1 + fval2;
	resk += wgk[jtwm1 - 1] * fsum;
	*resabs += wgk[jtwm1 - 1] * (ABS(fval1) + ABS(fval2));
/* L15: */
    }
    reskh = resk * .5;
    *resasc = wgk[10] * (d__1 = fc - reskh, ABS(d__1));
    for (j = 1; j <= 10; ++j) {
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
} /* pnl_dqk21 */

