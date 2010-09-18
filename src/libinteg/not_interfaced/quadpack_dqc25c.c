/*
 * This file contains the routine dqc25c from the quadpack distribution
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

extern int pnl_dqk15w(PnlFunc *, D_weight, double *, double *
                   , double *, double *, int *, double *, double 
                   *, double *, double *, double *, double *), 
       pnl_dqcheb(double *, double *, double *, double *);
extern double pnl_dqwgtc();

int pnl_dqc25c(PnlFunc * f, double *a, double *b, double 
	*c__, double *result, double *abserr, int *krul, int *
	neval)
{
    /* Initialized data */

    double x[11] = { .991444861373810411144557526928563,
	    .965925826289068286749743199728897,
	    .923879532511286756128183189396788,
	    .866025403784438646763723170752936,
	    .793353340291235164579776961501299,
	    .707106781186547524400844362104849,
	    .608761429008720639416097542898164,.5,
	    .382683432365089771728459984030399,
	    .258819045102520762348898837624048,
	    .130526192220051591548406227895489 };

    /* System generated locals */
    double d__1;

    /* Local variables */
    int i__, k;
    double u, p2, p3, p4, cc;
    int kp;
    double ak22, fval[25], res12, res24;
    int isym;
    double amom0, amom1, amom2, cheb12[13], cheb24[25], hlgth, centr;
    double resabs, resasc;

/* ***begin prologue  dqc25c */
/* ***date written   810101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a2a2,j4 */
/* ***keywords  25-point clenshaw-curtis integration */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  to compute i = integral of f*w over (a,b) with */
/*            error estimate, where w(x) = 1/(x-c) */
/* ***description */

/*        integration rules for the computation of cauchy */
/*        principal value integrals */
/*        standard fortran subroutine */
/*        double precision version */

/*        parameters */
/*           f      - double precision */
/*                    function subprogram defining the integrand function */
/*                    f(x). the actual name for f needs to be declared */
/*                    e x t e r n a l  in the driver program. */

/*           a      - double precision */
/*                    left end point of the integration interval */

/*           b      - double precision */
/*                    right end point of the integration interval, b.gt.a */

/*           c      - double precision */
/*                    parameter in the weight function */

/*           result - double precision */
/*                    approximation to the integral */
/*                    result is computed by using a generalized */
/*                    clenshaw-curtis method if c lies within ten percent */
/*                    of the integration interval. in the other case the */
/*                    15-point kronrod rule obtained by optimal addition */
/*                    of abscissae to the 7-point gauss rule, is applied. */

/*           abserr - double precision */
/*                    estimate of the modulus of the absolute error, */
/*                    which should equal or exceed ABS(i-result) */

/*           krul   - int */
/*                    key which is decreased by 1 if the 15-point */
/*                    gauss-kronrod scheme has been used */

/*           neval  - int */
/*                    number of integrand evaluations */

/* ....................................................................... */
/* ***references  (none) */
/* ***routines called  dqcheb,dqk15w,dqwgtc */
/* ***end prologue  dqc25c */




/*           the vector x contains the values cos(k*pi/24), */
/*           k = 1, ..., 11, to be used for the chebyshev series */
/*           expansion of f */


/*           list of major variables */
/*           ---------------------- */
/*           fval   - value of the function f at the points */
/*                    cos(k*pi/24),  k = 0, ..., 24 */
/*           cheb12 - chebyshev series expansion coefficients, */
/*                    for the function f, of degree 12 */
/*           cheb24 - chebyshev series expansion coefficients, */
/*                    for the function f, of degree 24 */
/*           res12  - approximation to the integral corresponding */
/*                    to the use of cheb12 */
/*           res24  - approximation to the integral corresponding */
/*                    to the use of cheb24 */
/*           dqwgtc - external function subprogram defining */
/*                    the weight function */
/*           hlgth  - half-length of the interval */
/*           centr  - mid point of the interval */


/*           check the position of c. */

/* ***first executable statement  dqc25c */
    cc = (2. * *c__ - *b - *a) / (*b - *a);
    if (ABS(cc) < 1.1) {
	goto L10;
    }

/*           apply the 15-point gauss-kronrod scheme. */

    --(*krul);
    pnl_dqk15w((PnlFunc *)f, (D_weight)pnl_dqwgtc, c__, &p2, &p3, &p4, &kp, a, b, result, 
	    abserr, &resabs, &resasc);
    *neval = 15;
    if (resasc == *abserr) {
	++(*krul);
    }
    goto L50;

/*           use the generalized clenshaw-curtis method. */

L10:
    hlgth = (*b - *a) * .5;
    centr = (*b + *a) * .5;
    *neval = 25;
    d__1 = hlgth + centr;
    fval[0] = PNL_EVAL_FUNC(f,d__1) * .5;
    fval[12] = PNL_EVAL_FUNC(f,centr);
    d__1 = centr - hlgth;
    fval[24] = PNL_EVAL_FUNC(f,d__1) * .5;
    for (i__ = 2; i__ <= 12; ++i__) {
	u = hlgth * x[i__ - 2];
	isym = 26 - i__;
	d__1 = u + centr;
	fval[i__ - 1] = PNL_EVAL_FUNC(f,d__1);
	d__1 = centr - u;
	fval[isym - 1] = PNL_EVAL_FUNC(f,d__1);
/* L20: */
    }

/*           compute the chebyshev series expansion. */

    pnl_dqcheb(x, fval, cheb12, cheb24);

/*           the modified chebyshev moments are computed by forward */
/*           recursion, using amom0 and amom1 as starting values. */

    amom0 = log((d__1 = (1. - cc) / (cc + 1.), ABS(d__1)));
    amom1 = cc * amom0 + 2.;
    res12 = cheb12[0] * amom0 + cheb12[1] * amom1;
    res24 = cheb24[0] * amom0 + cheb24[1] * amom1;
    for (k = 3; k <= 13; ++k) {
	amom2 = cc * 2. * amom1 - amom0;
	ak22 = (double) ((k - 2) * (k - 2));
	if (k / 2 << 1 == k) {
	    amom2 -= 4. / (ak22 - 1.);
	}
	res12 += cheb12[k - 1] * amom2;
	res24 += cheb24[k - 1] * amom2;
	amom0 = amom1;
	amom1 = amom2;
/* L30: */
    }
    for (k = 14; k <= 25; ++k) {
	amom2 = cc * 2. * amom1 - amom0;
	ak22 = (double) ((k - 2) * (k - 2));
	if (k / 2 << 1 == k) {
	    amom2 -= 4. / (ak22 - 1.);
	}
	res24 += cheb24[k - 1] * amom2;
	amom0 = amom1;
	amom1 = amom2;
/* L40: */
    }
    *result = res24;
    *abserr = (d__1 = res24 - res12, ABS(d__1));
L50:
    return 0;
} /* pnl_dqc25c */

