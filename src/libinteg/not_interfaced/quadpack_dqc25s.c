/*
 * This file contains the routine dqc25s from the quadpack distribution
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
                   *, double *, double *, double *, double *);
extern int pnl_dqcheb(double *, double *, double *, double *);
extern double pnl_dqwgts();

int pnl_dqc25s(PnlFunc * f, double *a, double *b, double 
	*bl, double *br, double *alfa, double *beta, double *
	ri, double *rj, double *rg, double *rh, double *
	result, double *abserr, double *resasc, int *integr, 
	int *nev)
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
    double d__1, d__2;

    /* Local variables */
    int i__;
    double u, dc, fix, fval[25], res12, res24;
    int isym;
    double cheb12[13], cheb24[25], hlgth, centr;
    double factor, resabs;

/* ***begin prologue  dqc25s */
/* ***date written   810101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a2a2 */
/* ***keywords  25-point clenshaw-curtis integration */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  to compute i = integral of f*w over (bl,br), with error */
/*            estimate, where the weight function w has a singular */
/*            behaviour of algebraico-logarithmic type at the points */
/*            a and/or b. (bl,br) is a part of (a,b). */
/* ***description */

/*        integration rules for integrands having algebraico-logarithmic */
/*        end point singularities */
/*        standard fortran subroutine */
/*        double precision version */

/*        parameters */
/*           f      - double precision */
/*                    function subprogram defining the integrand */
/*                    f(x). the actual name for f needs to be declared */
/*                    e x t e r n a l  in the driver program. */

/*           a      - double precision */
/*                    left end point of the original interval */

/*           b      - double precision */
/*                    right end point of the original interval, b.gt.a */

/*           bl     - double precision */
/*                    lower limit of integration, bl.ge.a */

/*           br     - double precision */
/*                    upper limit of integration, br.le.b */

/*           alfa   - double precision */
/*                    parameter in the weight function */

/*           beta   - double precision */
/*                    parameter in the weight function */

/*           ri,rj,rg,rh - double precision */
/*                    modified chebyshev moments for the application */
/*                    of the generalized clenshaw-curtis */
/*                    method (computed in subroutine dqmomo) */

/*           result - double precision */
/*                    approximation to the integral */
/*                    result is computed by using a generalized */
/*                    clenshaw-curtis method if b1 = a or br = b. */
/*                    in all other cases the 15-point kronrod */
/*                    rule is applied, obtained by optimal addition of */
/*                    abscissae to the 7-point gauss rule. */

/*           abserr - double precision */
/*                    estimate of the modulus of the absolute error, */
/*                    which should equal or exceed ABS(i-result) */

/*           resasc - double precision */
/*                    approximation to the integral of ABS(f*w-i/(b-a)) */

/*           integr - int */
/*                    which deterMINes the weight function */
/*                    = 1   w(x) = (x-a)**alfa*(b-x)**beta */
/*                    = 2   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a) */
/*                    = 3   w(x) = (x-a)**alfa*(b-x)**beta*log(b-x) */
/*                    = 4   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)* */
/*                                 log(b-x) */

/*           nev    - int */
/*                    number of integrand evaluations */
/* ***references  (none) */
/* ***routines called  dqcheb,dqk15w */
/* ***end prologue  dqc25s */




/*           the vector x contains the values cos(k*pi/24) */
/*           k = 1, ..., 11, to be used for the computation of the */
/*           chebyshev series expansion of f. */

    /* Parameter adjustments */
    --rh;
    --rg;
    --rj;
    --ri;

    /* Function Body */

/*           list of major variables */
/*           ----------------------- */

/*           fval   - value of the function f at the points */
/*                    (br-bl)*0.5*cos(k*pi/24)+(br+bl)*0.5 */
/*                    k = 0, ..., 24 */
/*           cheb12 - coefficients of the chebyshev series expansion */
/*                    of degree 12, for the function f, in the */
/*                    interval (bl,br) */
/*           cheb24 - coefficients of the chebyshev series expansion */
/*                    of degree 24, for the function f, in the */
/*                    interval (bl,br) */
/*           res12  - approximation to the integral obtained from cheb12 */
/*           res24  - approximation to the integral obtained from cheb24 */
/*           dqwgts - external function subprogram defining */
/*                    the four possible weight functions */
/*           hlgth  - half-length of the interval (bl,br) */
/*           centr  - mid point of the interval (bl,br) */

/* ***first executable statement  dqc25s */
    *nev = 25;
    if (*bl == *a && (*alfa != 0. || *integr == 2 || *integr == 4)) {
	goto L10;
    }
    if (*br == *b && (*beta != 0. || *integr == 3 || *integr == 4)) {
	goto L140;
    }

/*           if a.gt.bl and b.lt.br, apply the 15-point gauss-kronrod */
/*           scheme. */


    pnl_dqk15w((PnlFunc *)f, (D_weight)pnl_dqwgts, a, b, alfa, beta, integr, bl, br, result, 
	    abserr, &resabs, resasc);
    *nev = 15;
    goto L270;

/*           this part of the program is executed only if a = bl. */
/*           ---------------------------------------------------- */

/*           compute the chebyshev series expansion of the */
/*           following function */
/*           f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta */
/*                  *f(0.5*(br-a)*x+0.5*(br+a)) */

L10:
    hlgth = (*br - *bl) * .5;
    centr = (*br + *bl) * .5;
    fix = *b - centr;
    d__1 = hlgth + centr;
    d__2 = fix - hlgth;
    fval[0] = PNL_EVAL_FUNC(f,d__1) * .5 * pow(d__2, *beta);
    fval[12] = PNL_EVAL_FUNC(f,centr) * pow(fix, *beta);
    d__1 = centr - hlgth;
    d__2 = fix + hlgth;
    fval[24] = PNL_EVAL_FUNC(f,d__1) * .5 * pow(d__2, *beta);
    for (i__ = 2; i__ <= 12; ++i__) {
	u = hlgth * x[i__ - 2];
	isym = 26 - i__;
	d__1 = u + centr;
	d__2 = fix - u;
	fval[i__ - 1] = PNL_EVAL_FUNC(f,d__1) * pow(d__2, *beta);
	d__1 = centr - u;
	d__2 = fix + u;
	fval[isym - 1] = PNL_EVAL_FUNC(f,d__1) * pow(d__2, *beta);
/* L20: */
    }
    d__1 = *alfa + 1.;
    factor = pow(hlgth, d__1);
    *result = 0.;
    *abserr = 0.;
    res12 = 0.;
    res24 = 0.;
    if (*integr > 2) {
	goto L70;
    }
    pnl_dqcheb(x, fval, cheb12, cheb24);

/*           integr = 1  (or 2) */

    for (i__ = 1; i__ <= 13; ++i__) {
	res12 += cheb12[i__ - 1] * ri[i__];
	res24 += cheb24[i__ - 1] * ri[i__];
/* L30: */
    }
    for (i__ = 14; i__ <= 25; ++i__) {
	res24 += cheb24[i__ - 1] * ri[i__];
/* L40: */
    }
    if (*integr == 1) {
	goto L130;
    }

/*           integr = 2 */

    dc = log(*br - *bl);
    *result = res24 * dc;
    *abserr = (d__1 = (res24 - res12) * dc, ABS(d__1));
    res12 = 0.;
    res24 = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
	res12 += cheb12[i__ - 1] * rg[i__];
	res24 = res12 + cheb24[i__ - 1] * rg[i__];
/* L50: */
    }
    for (i__ = 14; i__ <= 25; ++i__) {
	res24 += cheb24[i__ - 1] * rg[i__];
/* L60: */
    }
    goto L130;

/*           compute the chebyshev series expansion of the */
/*           following function */
/*           f4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x) */

L70:
    fval[0] *= log(fix - hlgth);
    fval[12] *= log(fix);
    fval[24] *= log(fix + hlgth);
    for (i__ = 2; i__ <= 12; ++i__) {
	u = hlgth * x[i__ - 2];
	isym = 26 - i__;
	fval[i__ - 1] *= log(fix - u);
	fval[isym - 1] *= log(fix + u);
/* L80: */
    }
    pnl_dqcheb(x, fval, cheb12, cheb24);

/*           integr = 3  (or 4) */

    for (i__ = 1; i__ <= 13; ++i__) {
	res12 += cheb12[i__ - 1] * ri[i__];
	res24 += cheb24[i__ - 1] * ri[i__];
/* L90: */
    }
    for (i__ = 14; i__ <= 25; ++i__) {
	res24 += cheb24[i__ - 1] * ri[i__];
/* L100: */
    }
    if (*integr == 3) {
	goto L130;
    }

/*           integr = 4 */

    dc = log(*br - *bl);
    *result = res24 * dc;
    *abserr = (d__1 = (res24 - res12) * dc, ABS(d__1));
    res12 = 0.;
    res24 = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
	res12 += cheb12[i__ - 1] * rg[i__];
	res24 += cheb24[i__ - 1] * rg[i__];
/* L110: */
    }
    for (i__ = 14; i__ <= 25; ++i__) {
	res24 += cheb24[i__ - 1] * rg[i__];
/* L120: */
    }
L130:
    *result = (*result + res24) * factor;
    *abserr = (*abserr + (d__1 = res24 - res12, ABS(d__1))) * factor;
    goto L270;

/*           this part of the program is executed only if b = br. */
/*           ---------------------------------------------------- */

/*           compute the chebyshev series expansion of the */
/*           following function */
/*           f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa */
/*                *f(0.5*(b-bl)*x+0.5*(b+bl)) */

L140:
    hlgth = (*br - *bl) * .5;
    centr = (*br + *bl) * .5;
    fix = centr - *a;
    d__1 = hlgth + centr;
    d__2 = fix + hlgth;
    fval[0] = PNL_EVAL_FUNC(f,d__1) * .5 * pow(d__2, *alfa);
    fval[12] = PNL_EVAL_FUNC(f,centr) * pow(fix, *alfa);
    d__1 = centr - hlgth;
    d__2 = fix - hlgth;
    fval[24] = PNL_EVAL_FUNC(f,d__1) * .5 * pow(d__2, *alfa);
    for (i__ = 2; i__ <= 12; ++i__) {
	u = hlgth * x[i__ - 2];
	isym = 26 - i__;
	d__1 = u + centr;
	d__2 = fix + u;
	fval[i__ - 1] = PNL_EVAL_FUNC(f,d__1) * pow(d__2, *alfa);
	d__1 = centr - u;
	d__2 = fix - u;
	fval[isym - 1] = PNL_EVAL_FUNC(f,d__1) * pow(d__2, *alfa);
/* L150: */
    }
    d__1 = *beta + 1.;
    factor = pow(hlgth, d__1);
    *result = 0.;
    *abserr = 0.;
    res12 = 0.;
    res24 = 0.;
    if (*integr == 2 || *integr == 4) {
	goto L200;
    }

/*           integr = 1  (or 3) */

    pnl_dqcheb(x, fval, cheb12, cheb24);
    for (i__ = 1; i__ <= 13; ++i__) {
	res12 += cheb12[i__ - 1] * rj[i__];
	res24 += cheb24[i__ - 1] * rj[i__];
/* L160: */
    }
    for (i__ = 14; i__ <= 25; ++i__) {
	res24 += cheb24[i__ - 1] * rj[i__];
/* L170: */
    }
    if (*integr == 1) {
	goto L260;
    }

/*           integr = 3 */

    dc = log(*br - *bl);
    *result = res24 * dc;
    *abserr = (d__1 = (res24 - res12) * dc, ABS(d__1));
    res12 = 0.;
    res24 = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
	res12 += cheb12[i__ - 1] * rh[i__];
	res24 += cheb24[i__ - 1] * rh[i__];
/* L180: */
    }
    for (i__ = 14; i__ <= 25; ++i__) {
	res24 += cheb24[i__ - 1] * rh[i__];
/* L190: */
    }
    goto L260;

/*           compute the chebyshev series expansion of the */
/*           following function */
/*           f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a)) */

L200:
    fval[0] *= log(hlgth + fix);
    fval[12] *= log(fix);
    fval[24] *= log(fix - hlgth);
    for (i__ = 2; i__ <= 12; ++i__) {
	u = hlgth * x[i__ - 2];
	isym = 26 - i__;
	fval[i__ - 1] *= log(u + fix);
	fval[isym - 1] *= log(fix - u);
/* L210: */
    }
    pnl_dqcheb(x, fval, cheb12, cheb24);

/*           integr = 2  (or 4) */

    for (i__ = 1; i__ <= 13; ++i__) {
	res12 += cheb12[i__ - 1] * rj[i__];
	res24 += cheb24[i__ - 1] * rj[i__];
/* L220: */
    }
    for (i__ = 14; i__ <= 25; ++i__) {
	res24 += cheb24[i__ - 1] * rj[i__];
/* L230: */
    }
    if (*integr == 2) {
	goto L260;
    }
    dc = log(*br - *bl);
    *result = res24 * dc;
    *abserr = (d__1 = (res24 - res12) * dc, ABS(d__1));
    res12 = 0.;
    res24 = 0.;

/*           integr = 4 */

    for (i__ = 1; i__ <= 13; ++i__) {
	res12 += cheb12[i__ - 1] * rh[i__];
	res24 += cheb24[i__ - 1] * rh[i__];
/* L240: */
    }
    for (i__ = 14; i__ <= 25; ++i__) {
	res24 += cheb24[i__ - 1] * rh[i__];
/* L250: */
    }
L260:
    *result = (*result + res24) * factor;
    *abserr = (*abserr + (d__1 = res24 - res12, ABS(d__1))) * factor;
L270:
    return 0;
} /* pnl_dqc25s */

