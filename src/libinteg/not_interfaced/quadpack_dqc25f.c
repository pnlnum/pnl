/*
 * This file contains the routine dqc25f from the quadpack distribution
 * available on www.netlib.org
 *
 * This routine has been translated from Fortran using f2c and then slightly
 * modified by Jérôme Lelong for integration into Pnl. All Quadpack names have
 * been prefixed by pnl to avoid possible clash in names. Calls to d1mach have
 * been replaced by calls to dlamch coming from Lapack
 */

#include "pnl/pnl_machine.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/clapack.h"

typedef double(*D_weight)();

extern int pnl_dqk15w(PnlFunc *, D_weight, double *, double *
                   , double *, double *, int *, double *, double 
                   *, double *, double *, double *, double *);
extern int pnl_dqcheb(double *, double *, 
                   double *, double *);
extern double pnl_dqwgtf();

int pnl_dqc25f(PnlFunc * f, double *a, double *b, double 
	*omega, int *integr, int *nrmom, int *MAXp1, int *
	ksave, double *result, double *abserr, int *neval, 
	double *resabs, double *resasc, int *momcom, double *
	chebmo)
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
    int chebmo_dim1, chebmo_offset, i__1, nrhs;
    double d__1, d__2;

    /* Local variables */
    double d__[25];
    int i__, j, k, m;
    double v[28], d1[25], d2[25], p2, p3, p4, ac, an, as, an2, ass, par2, 
	    conc, asap, par22, fval[25], estc, cons;
    int iers;
    double ests;
    int isym, noeq1;
    double cheb12[13], cheb24[25], resc12, resc24, hlgth, centr;
    double ress12, ress24, oflow;
    int noequ;
    
    double cospar;
    double parint, sinpar;

/* ***begin prologue  dqc25f */
/* ***date written   810101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a2a2 */
/* ***keywords  integration rules for functions with cos or sin */
/*             factor, clenshaw-curtis, gauss-kronrod */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  to compute the integral i=integral of f(x) over (a,b) */
/*            where w(x) = cos(omega*x) or w(x)=sin(omega*x) and to */
/*            compute j = integral of ABS(f) over (a,b). for small value */
/*            of omega or small intervals (a,b) the 15-point gauss-kronro */
/*            rule is used. otherwise a generalized clenshaw-curtis */
/*            method is used. */
/* ***description */

/*        integration rules for functions with cos or sin factor */
/*        standard fortran subroutine */
/*        double precision version */

/*        parameters */
/*         on entry */
/*           f      - double precision */
/*                    function subprogram defining the integrand */
/*                    function f(x). the actual name for f needs to */
/*                    be declared e x t e r n a l in the calling program. */

/*           a      - double precision */
/*                    lower limit of integration */

/*           b      - double precision */
/*                    upper limit of integration */

/*           omega  - double precision */
/*                    parameter in the weight function */

/*           integr - int */
/*                    indicates which weight function is to be used */
/*                       integr = 1   w(x) = cos(omega*x) */
/*                       integr = 2   w(x) = sin(omega*x) */

/*           nrmom  - int */
/*                    the length of interval (a,b) is equal to the length */
/*                    of the original integration interval divided by */
/*                    2**nrmom (we suppose that the routine is used in an */
/*                    adaptive integration process, otherwise set */
/*                    nrmom = 0). nrmom must be zero at the first call. */

/*           MAXp1  - int */
/*                    gives an upper bound on the number of chebyshev */
/*                    moments which can be stored, i.e. for the */
/*                    intervals of lengths ABS(bb-aa)*2**(-l), */
/*                    l = 0,1,2, ..., MAXp1-2. */

/*           ksave  - int */
/*                    key which is one when the moments for the */
/*                    current interval have been computed */

/*         on return */
/*           result - double precision */
/*                    approximation to the integral i */

/*           abserr - double precision */
/*                    estimate of the modulus of the absolute */
/*                    error, which should equal or exceed ABS(i-result) */

/*           neval  - int */
/*                    number of integrand evaluations */

/*           resabs - double precision */
/*                    approximation to the integral j */

/*           resasc - double precision */
/*                    approximation to the integral of ABS(f-i/(b-a)) */

/*         on entry and return */
/*           momcom - int */
/*                    for each interval length we need to compute the */
/*                    chebyshev moments. momcom counts the number of */
/*                    intervals for which these moments have already been */
/*                    computed. if nrmom.lt.momcom or ksave = 1, the */
/*                    chebyshev moments for the interval (a,b) have */
/*                    already been computed and stored, otherwise we */
/*                    compute them and we increase momcom. */

/*           chebmo - double precision */
/*                    array of dimension at least (MAXp1,25) containing */
/*                    the modified chebyshev moments for the first momcom */
/*                    momcom interval lengths */

/* ...................................................................... */
/* ***references  (none) */
/* ***routines called  dqcheb,dqk15w,dqwgtf */
/* ***end prologue  dqc25f */




/*           the vector x contains the values cos(k*pi/24) */
/*           k = 1, ...,11, to be used for the chebyshev expansion of f */

    /* Parameter adjustments */
    chebmo_dim1 = *MAXp1;
    chebmo_offset = 1 + chebmo_dim1;
    chebmo -= chebmo_offset;

    /* Function Body */

/*           list of major variables */
/*           ----------------------- */

/*           centr  - mid point of the integration interval */
/*           hlgth  - half-length of the integration interval */
/*           fval   - value of the function f at the points */
/*                    (b-a)*0.5*cos(k*pi/12) + (b+a)*0.5, k = 0, ..., 24 */
/*           cheb12 - coefficients of the chebyshev series expansion */
/*                    of degree 12, for the function f, in the */
/*                    interval (a,b) */
/*           cheb24 - coefficients of the chebyshev series expansion */
/*                    of degree 24, for the function f, in the */
/*                    interval (a,b) */
/*           resc12 - approximation to the integral of */
/*                    cos(0.5*(b-a)*omega*x)*f(0.5*(b-a)*x+0.5*(b+a)) */
/*                    over (-1,+1), using the chebyshev series */
/*                    expansion of degree 12 */
/*           resc24 - approximation to the same integral, using the */
/*                    chebyshev series expansion of degree 24 */
/*           ress12 - the analogue of resc12 for the sine */
/*           ress24 - the analogue of resc24 for the sine */


/*           machine dependent constant */
/*           -------------------------- */

/*           oflow is the largest positive magnitude. */

/* ***first executable statement  dqc25f */
    oflow = pnl_d1mach(2);

    centr = (*b + *a) * .5;
    hlgth = (*b - *a) * .5;
    parint = *omega * hlgth;

/*           compute the integral using the 15-point gauss-kronrod */
/*           formula if the value of the parameter in the integrand */
/*           is small. */

    if (ABS(parint) > 2.) {
	goto L10;
    }
    pnl_dqk15w((PnlFunc *)f, (D_weight)pnl_dqwgtf, omega, &p2, &p3, &p4, integr, a, b, 
	    result, abserr, resabs, resasc);
    *neval = 15;
    goto L170;

/*           compute the integral using the generalized clenshaw- */
/*           curtis method. */

L10:
    conc = hlgth * cos(centr * *omega);
    cons = hlgth * sin(centr * *omega);
    *resasc = oflow;
    *neval = 25;

/*           check whether the chebyshev moments for this interval */
/*           have already been computed. */

    if (*nrmom < *momcom || *ksave == 1) {
	goto L120;
    }

/*           compute a new set of chebyshev moments. */

    m = *momcom + 1;
    par2 = parint * parint;
    par22 = par2 + 2.;
    sinpar = sin(parint);
    cospar = cos(parint);

/*           compute the chebyshev moments with respect to cosine. */

    v[0] = sinpar * 2. / parint;
    v[1] = (cospar * 8. + (par2 + par2 - 8.) * sinpar / parint) / par2;
    v[2] = ((par2 - 12.) * 32. * cospar + ((par2 - 80.) * par2 + 192.) * 2. * 
	    sinpar / parint) / (par2 * par2);
    ac = cospar * 8.;
    as = parint * 24. * sinpar;
    if (ABS(parint) > 24.) {
	goto L30;
    }

/*           compute the chebyshev moments as the solutions of a */
/*           boundary value problem with 1 initial value (v(3)) and 1 */
/*           end value (computed using an asymptotic formula). */

    noequ = 25;
    noeq1 = noequ - 1;
    an = 6.;
    i__1 = noeq1;
    for (k = 1; k <= i__1; ++k) {
	an2 = an * an;
	d__[k - 1] = (an2 - 4.) * -2. * (par22 - an2 - an2);
	d2[k - 1] = (an - 1.) * (an - 2.) * par2;
	d1[k] = (an + 3.) * (an + 4.) * par2;
	v[k + 2] = as - (an2 - 4.) * ac;
	an += 2.;
/* L20: */
    }
    an2 = an * an;
    d__[noequ - 1] = (an2 - 4.) * -2. * (par22 - an2 - an2);
    v[noequ + 2] = as - (an2 - 4.) * ac;
    v[3] -= par2 * 56. * v[2];
    ass = parint * sinpar;
    asap = (((((par2 * 210. - 1.) * cospar - (par2 * 105. - 63.) * ass) / an2 
	    - (1. - par2 * 15.) * cospar + ass * 15.) / an2 - cospar + ass * 
	    3.) / an2 - cospar) / an2;
    v[noequ + 2] -= asap * 2. * par2 * (an - 1.) * (an - 2.);

/*           solve the tridiagonal system by means of gaussian */
/*           eliMINation with partial pivoting. */

    /*
     * Call to dgtsl from linpack replaced by call to dgtsv.
     * The storage convention of the lower diagonal is different in linpack and
     * in lapack, this is why we use d1 - 1 instead of d1
     */
    nrhs = 1;
    C2F(dgtsv)(&noequ, &nrhs, d1 - 1, d__, d2, &v[3], &noequ, &iers);
    goto L50;

/*           compute the chebyshev moments by means of forward */
/*           recursion. */

L30:
    an = 4.;
    for (i__ = 4; i__ <= 13; ++i__) {
	an2 = an * an;
	v[i__ - 1] = ((an2 - 4.) * ((par22 - an2 - an2) * 2. * v[i__ - 2] - 
		ac) + as - par2 * (an + 1.) * (an + 2.) * v[i__ - 3]) / (par2 
		* (an - 1.) * (an - 2.));
	an += 2.;
/* L40: */
    }
L50:
    for (j = 1; j <= 13; ++j) {
	chebmo[m + ((j << 1) - 1) * chebmo_dim1] = v[j - 1];
/* L60: */
    }

/*           compute the chebyshev moments with respect to sine. */

    v[0] = (sinpar - parint * cospar) * 2. / par2;
    v[1] = (18. - 48. / par2) * sinpar / par2 + (48. / par2 - 2.) * cospar / 
	    parint;
    ac = parint * -24. * cospar;
    as = sinpar * -8.;
    if (ABS(parint) > 24.) {
	goto L80;
    }

/*           compute the chebyshev moments as the solutions of a boundary */
/*           value problem with 1 initial value (v(2)) and 1 end value */
/*           (computed using an asymptotic formula). */

    an = 5.;
    i__1 = noeq1;
    for (k = 1; k <= i__1; ++k) {
	an2 = an * an;
	d__[k - 1] = (an2 - 4.) * -2. * (par22 - an2 - an2);
	d2[k - 1] = (an - 1.) * (an - 2.) * par2;
	d1[k] = (an + 3.) * (an + 4.) * par2;
	v[k + 1] = ac + (an2 - 4.) * as;
	an += 2.;
/* L70: */
    }
    an2 = an * an;
    d__[noequ - 1] = (an2 - 4.) * -2. * (par22 - an2 - an2);
    v[noequ + 1] = ac + (an2 - 4.) * as;
    v[2] -= par2 * 42. * v[1];
    ass = parint * cospar;
    asap = (((((par2 * 105. - 63.) * ass + (par2 * 210. - 1.) * sinpar) / an2 
	    + (par2 * 15. - 1.) * sinpar - ass * 15.) / an2 - ass * 3. - 
	    sinpar) / an2 - sinpar) / an2;
    v[noequ + 1] -= asap * 2. * par2 * (an - 1.) * (an - 2.);

/*           solve the tridiagonal system by means of gaussian */
/*           eliMINation with partial pivoting. */

    /*
     * Call to dgtsl from linpack replaced by call to dgtsv.
     * The storage convention of the lower diagonal is different in linpack and
     * in lapack, this is why we use d1 - 1 instead of d1
     */
    nrhs = 1;
    C2F(dgtsv)(&noequ, &nrhs, d1 - 1, d__, d2, &v[2], &noequ, &iers);
    goto L100;

/*           compute the chebyshev moments by means of forward recursion. */

L80:
    an = 3.;
    for (i__ = 3; i__ <= 12; ++i__) {
	an2 = an * an;
	v[i__ - 1] = ((an2 - 4.) * ((par22 - an2 - an2) * 2. * v[i__ - 2] + 
		as) + ac - par2 * (an + 1.) * (an + 2.) * v[i__ - 3]) / (par2 
		* (an - 1.) * (an - 2.));
	an += 2.;
/* L90: */
    }
L100:
    for (j = 1; j <= 12; ++j) {
	chebmo[m + (j << 1) * chebmo_dim1] = v[j - 1];
/* L110: */
    }
L120:
    if (*nrmom < *momcom) {
	m = *nrmom + 1;
    }
    if (*momcom < *MAXp1 - 1 && *nrmom >= *momcom) {
	++(*momcom);
    }

/*           compute the coefficients of the chebyshev expansions */
/*           of degrees 12 and 24 of the function f. */

    d__1 = centr + hlgth;
    fval[0] = PNL_EVAL_FUNC(f,d__1) * .5;
    fval[12] = PNL_EVAL_FUNC(f,centr);
    d__1 = centr - hlgth;
    fval[24] = PNL_EVAL_FUNC(f,d__1) * .5;
    for (i__ = 2; i__ <= 12; ++i__) {
	isym = 26 - i__;
	d__1 = hlgth * x[i__ - 2] + centr;
	fval[i__ - 1] = PNL_EVAL_FUNC(f,d__1);
	d__1 = centr - hlgth * x[i__ - 2];
	fval[isym - 1] = PNL_EVAL_FUNC(f,d__1);
/* L130: */
    }
    pnl_dqcheb(x, fval, cheb12, cheb24);

/*           compute the integral and error estimates. */

    resc12 = cheb12[12] * chebmo[m + chebmo_dim1 * 13];
    ress12 = 0.;
    k = 11;
    for (j = 1; j <= 6; ++j) {
	resc12 += cheb12[k - 1] * chebmo[m + k * chebmo_dim1];
	ress12 += cheb12[k] * chebmo[m + (k + 1) * chebmo_dim1];
	k += -2;
/* L140: */
    }
    resc24 = cheb24[24] * chebmo[m + chebmo_dim1 * 25];
    ress24 = 0.;
    *resabs = ABS(cheb24[24]);
    k = 23;
    for (j = 1; j <= 12; ++j) {
	resc24 += cheb24[k - 1] * chebmo[m + k * chebmo_dim1];
	ress24 += cheb24[k] * chebmo[m + (k + 1) * chebmo_dim1];
	*resabs = (d__1 = cheb24[k - 1], ABS(d__1)) + (d__2 = cheb24[k], ABS(
		d__2));
	k += -2;
/* L150: */
    }
    estc = (d__1 = resc24 - resc12, ABS(d__1));
    ests = (d__1 = ress24 - ress12, ABS(d__1));
    *resabs *= ABS(hlgth);
    if (*integr == 2) {
	goto L160;
    }
    *result = conc * resc24 - cons * ress24;
    *abserr = (d__1 = conc * estc, ABS(d__1)) + (d__2 = cons * ests, ABS(d__2)
	    );
    goto L170;
L160:
    *result = conc * ress24 + cons * resc24;
    *abserr = (d__1 = conc * ests, ABS(d__1)) + (d__2 = cons * estc, ABS(d__2)
	    );
L170:
    return 0;
} /* pnl_dqc25f */

