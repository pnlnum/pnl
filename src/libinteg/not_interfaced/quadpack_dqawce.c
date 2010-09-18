/*
 * This file contains the routine dqawce from the quadpack distribution
 * available on www.netlib.org
 *
 * This routine has been translated from Fortran using f2c and then slightly
 * modified by Jérôme Lelong for integration into Pnl. All Quadpack names have
 * been prefixed by pnl to avoid possible clash in names. Calls to d1mach have
 * been replaced by calls to dlamch coming from Lapack
 */

#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"

extern int pnl_dqc25c(PnlFunc *, double *, double *, double *, double *, double *, int *, int *);
extern int pnl_dqpsrt(int *, int *, int *, double *, double *, int *, int *);

int pnl_dqawce(PnlFunc * f, double *a, double *b, double 
	*c__, double *epsabs, double *epsrel, int *limit, 
	double *result, double *abserr, int *neval, int *ier, 
	double *alist__, double *blist, double *rlist, double 
	*elist, int *iord, int *last)
{
    /* System generated locals */
    int i__1;
    double d__1, d__2;

    /* Local variables */
    int k;
    double a1, a2, b1, b2, aa, bb;
    int nev;
    double area, area1, area2, area12;
    double erro12;
    int krule, nrMAX;
    double uflow;
    
    int iroff1, iroff2;
    double error1, error2, epmach, errbnd, errMAX;
    int MAXerr;
    double errsum;

/* ***begin prologue  dqawce */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a2a1,j4 */
/* ***keywords  automatic integrator, special-purpose, */
/*             cauchy principal value, clenshaw-curtis method */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***  purpose  the routine calculates an approximation result to a */
/*              cauchy principal value i = integral of f*w over (a,b) */
/*              (w(x) = 1/(x-c), (c.ne.a, c.ne.b), hopefully satisfying */
/*              following claim for accuracy */
/*              ABS(i-result).le.MAX(epsabs,epsrel*ABS(i)) */
/* ***description */

/*        computation of a cauchy principal value */
/*        standard fortran subroutine */
/*        double precision version */

/*        parameters */
/*         on entry */
/*            f      - double precision */
/*                     function subprogram defining the integrand */
/*                     function f(x). the actual name for f needs to be */
/*                     declared e x t e r n a l in the driver program. */

/*            a      - double precision */
/*                     lower limit of integration */

/*            b      - double precision */
/*                     upper limit of integration */

/*            c      - double precision */
/*                     parameter in the weight function, c.ne.a, c.ne.b */
/*                     if c = a or c = b, the routine will end with */
/*                     ier = 6. */

/*            epsabs - double precision */
/*                     absolute accuracy requested */
/*            epsrel - double precision */
/*                     relative accuracy requested */
/*                     if  epsabs.le.0 */
/*                     and epsrel.lt.MAX(50*rel.mach.acc.,0.5d-28), */
/*                     the routine will end with ier = 6. */

/*            limit  - int */
/*                     gives an upper bound on the number of subintervals */
/*                     in the partition of (a,b), limit.ge.1 */

/*         on return */
/*            result - double precision */
/*                     approximation to the integral */

/*            abserr - double precision */
/*                     estimate of the modulus of the absolute error, */
/*                     which should equal or exceed ABS(i-result) */

/*            neval  - int */
/*                     number of integrand evaluations */

/*            ier    - int */
/*                     ier = 0 normal and reliable terMINation of the */
/*                             routine. it is assumed that the requested */
/*                             accuracy has been achieved. */
/*                     ier.gt.0 abnormal terMINation of the routine */
/*                             the estimates for integral and error are */
/*                             less reliable. it is assumed that the */
/*                             requested accuracy has not been achieved. */
/*            error messages */
/*                     ier = 1 MAXimum number of subdivisions allowed */
/*                             has been achieved. one can allow more sub- */
/*                             divisions by increasing the value of */
/*                             limit. however, if this yields no */
/*                             improvement it is advised to analyze the */
/*                             the integrand, in order to deterMINe the */
/*                             the integration difficulties. if the */
/*                             position of a local difficulty can be */
/*                             deterMINed (e.g. singularity, */
/*                             discontinuity within the interval) one */
/*                             will probably gain from splitting up the */
/*                             interval at this point and calling */
/*                             appropriate integrators on the subranges. */
/*                         = 2 the occurrence of roundoff error is detec- */
/*                             ted, which prevents the requested */
/*                             tolerance from being achieved. */
/*                         = 3 extremely bad integrand behaviour */
/*                             occurs at some interior points of */
/*                             the integration interval. */
/*                         = 6 the input is invalid, because */
/*                             c = a or c = b or */
/*                             (epsabs.le.0 and */
/*                              epsrel.lt.MAX(50*rel.mach.acc.,0.5d-28)) */
/*                             or limit.lt.1. */
/*                             result, abserr, neval, rlist(1), elist(1), */
/*                             iord(1) and last are set to zero. alist(1) */
/*                             and blist(1) are set to a and b */
/*                             respectively. */

/*            alist   - double precision */
/*                      vector of dimension at least limit, the first */
/*                       last  elements of which are the left */
/*                      end points of the subintervals in the partition */
/*                      of the given integration range (a,b) */

/*            blist   - double precision */
/*                      vector of dimension at least limit, the first */
/*                       last  elements of which are the right */
/*                      end points of the subintervals in the partition */
/*                      of the given integration range (a,b) */

/*            rlist   - double precision */
/*                      vector of dimension at least limit, the first */
/*                       last  elements of which are the integral */
/*                      approximations on the subintervals */

/*            elist   - double precision */
/*                      vector of dimension limit, the first  last */
/*                      elements of which are the moduli of the absolute */
/*                      error estimates on the subintervals */

/*            iord    - int */
/*                      vector of dimension at least limit, the first k */
/*                      elements of which are pointers to the error */
/*                      estimates over the subintervals, so that */
/*                      elist(iord(1)), ..., elist(iord(k)) with k = last */
/*                      if last.le.(limit/2+2), and k = limit+1-last */
/*                      otherwise, form a decreasing sequence */

/*            last    - int */
/*                      number of subintervals actually produced in */
/*                      the subdivision process */

/* ***references  (none) */
/* ***routines called  dqc25c,dqpsrt */
/* ***end prologue  dqawce */




/*            list of major variables */
/*            ----------------------- */

/*           alist     - list of left end points of all subintervals */
/*                       considered up to now */
/*           blist     - list of right end points of all subintervals */
/*                       considered up to now */
/*           rlist(i)  - approximation to the integral over */
/*                       (alist(i),blist(i)) */
/*           elist(i)  - error estimate applying to rlist(i) */
/*           MAXerr    - pointer to the interval with largest */
/*                       error estimate */
/*           errMAX    - elist(MAXerr) */
/*           area      - sum of the integrals over the subintervals */
/*           errsum    - sum of the errors over the subintervals */
/*           errbnd    - requested accuracy MAX(epsabs,epsrel* */
/*                       ABS(result)) */
/*           *****1    - variable for the left subinterval */
/*           *****2    - variable for the right subinterval */
/*           last      - index for subdivision */


/*            machine dependent constants */
/*            --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */

/* ***first executable statement  dqawce */
    /* Parameter adjustments */
    --iord;
    --elist;
    --rlist;
    --blist;
    --alist__;

    /* Function Body */
    epmach = pnl_d1mach(4);
    uflow = pnl_d1mach(1);


/*           test on validity of parameters */
/*           ------------------------------ */

    *ier = 6;
    *neval = 0;
    *last = 0;
    alist__[1] = *a;
    blist[1] = *b;
    rlist[1] = 0.;
    elist[1] = 0.;
    iord[1] = 0;
    *result = 0.;
    *abserr = 0.;
/* Computing MAX */
    d__1 = epmach * 50.;
    if (*c__ == *a || *c__ == *b || (*epsabs <= 0. && *epsrel < MAX(d__1,5e-29))) {
	goto L999;
    }

/*           first approximation to the integral */
/*           ----------------------------------- */

    aa = *a;
    bb = *b;
    if (*a <= *b) {
	goto L10;
    }
    aa = *b;
    bb = *a;
L10:
    *ier = 0;
    krule = 1;
    pnl_dqc25c((PnlFunc *)f, &aa, &bb, c__, result, abserr, &krule, neval);
    *last = 1;
    rlist[1] = *result;
    elist[1] = *abserr;
    iord[1] = 1;
    alist__[1] = *a;
    blist[1] = *b;

/*           test on accuracy */

/* Computing MAX */
    d__1 = *epsabs, d__2 = *epsrel * ABS(*result);
    errbnd = MAX(d__1,d__2);
    if (*limit == 1) {
	*ier = 1;
    }
/* Computing MIN */
    d__1 = ABS(*result) * .01;
    if (*abserr < MIN(d__1,errbnd) || *ier == 1) {
	goto L70;
    }

/*           initialization */
/*           -------------- */

    alist__[1] = aa;
    blist[1] = bb;
    rlist[1] = *result;
    errMAX = *abserr;
    MAXerr = 1;
    area = *result;
    errsum = *abserr;
    nrMAX = 1;
    iroff1 = 0;
    iroff2 = 0;

/*           main do-loop */
/*           ------------ */

    i__1 = *limit;
    for (*last = 2; *last <= i__1; ++(*last)) {

/*           bisect the subinterval with nrMAX-th largest */
/*           error estimate. */

	a1 = alist__[MAXerr];
	b1 = (alist__[MAXerr] + blist[MAXerr]) * .5;
	b2 = blist[MAXerr];
	if (*c__ <= b1 && *c__ > a1) {
	    b1 = (*c__ + b2) * .5;
	}
	if (*c__ > b1 && *c__ < b2) {
	    b1 = (a1 + *c__) * .5;
	}
	a2 = b1;
	krule = 2;
	pnl_dqc25c((PnlFunc *)f, &a1, &b1, c__, &area1, &error1, &krule, &nev);
	*neval += nev;
	pnl_dqc25c((PnlFunc *)f, &a2, &b2, c__, &area2, &error2, &krule, &nev);
	*neval += nev;

/*           improve previous approximations to integral */
/*           and error and test for accuracy. */

	area12 = area1 + area2;
	erro12 = error1 + error2;
	errsum = errsum + erro12 - errMAX;
	area = area + area12 - rlist[MAXerr];
	if ((d__1 = rlist[MAXerr] - area12, ABS(d__1)) < ABS(area12) * 1e-5 &&
		 erro12 >= errMAX * .99 && krule == 0) {
	    ++iroff1;
	}
	if (*last > 10 && erro12 > errMAX && krule == 0) {
	    ++iroff2;
	}
	rlist[MAXerr] = area1;
	rlist[*last] = area2;
/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * ABS(area);
	errbnd = MAX(d__1,d__2);
	if (errsum <= errbnd) {
	    goto L15;
	}

/*           test for roundoff error and eventually set error flag. */

	if (iroff1 >= 6 && iroff2 > 20) {
	    *ier = 2;
	}

/*           set error flag in the case that number of interval */
/*           bisections exceeds limit. */

	if (*last == *limit) {
	    *ier = 1;
	}

/*           set error flag in the case of bad integrand behaviour */
/*           at a point of the integration range. */

/* Computing MAX */
	d__1 = ABS(a1), d__2 = ABS(b2);
	if (MAX(d__1,d__2) <= (epmach * 100. + 1.) * (ABS(a2) + uflow * 1e3)) 
		{
	    *ier = 3;
	}

/*           append the newly-created intervals to the list. */

L15:
	if (error2 > error1) {
	    goto L20;
	}
	alist__[*last] = a2;
	blist[MAXerr] = b1;
	blist[*last] = b2;
	elist[MAXerr] = error1;
	elist[*last] = error2;
	goto L30;
L20:
	alist__[MAXerr] = a2;
	alist__[*last] = a1;
	blist[*last] = b1;
	rlist[MAXerr] = area2;
	rlist[*last] = area1;
	elist[MAXerr] = error2;
	elist[*last] = error1;

/*           call subroutine dqpsrt to maintain the descending ordering */
/*           in the list of error estimates and select the subinterval */
/*           with nrMAX-th largest error estimate (to be bisected next). */

L30:
	pnl_dqpsrt(limit, last, &MAXerr, &errMAX, &elist[1], &iord[1], &nrMAX);
/* ***jump out of do-loop */
	if (*ier != 0 || errsum <= errbnd) {
	    goto L50;
	}
/* L40: */
    }

/*           compute final result. */
/*           --------------------- */

L50:
    *result = 0.;
    i__1 = *last;
    for (k = 1; k <= i__1; ++k) {
	*result += rlist[k];
/* L60: */
    }
    *abserr = errsum;
L70:
    if (aa == *b) {
	*result = -(*result);
    }
L999:
    return 0;
} /* pnl_dqawce */

