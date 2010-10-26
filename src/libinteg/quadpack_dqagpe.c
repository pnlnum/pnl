/*
 * This file contains the routine dqagpe from the quadpack distribution
 * available on www.netlib.org
 *
 * This routine has been translated from Fortran using f2c and then slightly
 * modified by Jérôme Lelong for integration into Pnl. All Quadpack names have
 * been prefixed by pnl to avoid possible clash in names. Calls to d1mach have
 * been replaced by calls to dlamch coming from Lapack
 */

#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"

extern int pnl_dqk21(PnlFunc *, double *, double *, double *, double *, double *, double *);
extern int pnl_dqelg(int *, double *, double *, double *, double *, int *);
extern int pnl_dqpsrt(int *, int *, int *, double *, double *, int *, int *);

int pnl_dqagpe(PnlFunc * f, double *a, double *b, int *
	npts2, double *points, double *epsabs, double *epsrel,
	int *limit, double *result, double *abserr, int *
	neval, int *ier, double *alist__, double *blist,
	double *rlist, double *elist, double *pts, int *iord,
	int *level, int *ndin, int *last)
{
    /* System generated locals */
    int i__1, i__2;
    double d__1, d__2;

    /* Local variables */
    int i__, j, k;
    double a1, a2, b1, b2;
    int id, ip1, ind1, ind2;
    double area;
    double resa, dres, sign;
    int ksgn;
    double temp;
    int nres, nint, jlow, npts;
    double area1, area2, area12;
    double erro12;
    int ierro;
    double defab1, defab2;
    int ktMIN, nrMAX;
    double oflow, uflow;

    int noext;
    int iroff1, iroff2, iroff3;
    double res3la[3];
    int nintp1;
    double error1, error2, rlist2[52];
    int numrl2;
    double defabs, epmach, erlarg, abseps, correc=0, errbnd, resabs;
    int jupbnd;
    double erlast;
    int levMAX;
    double errMAX;
    int MAXerr, levcur;
    double reseps;
    int extrap;
    double ertest, errsum;

/* ***begin prologue  dqagpe */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a2a1 */
/* ***keywords  automatic integrator, general-purpose, */
/*             singularities at user specified points, */
/*             extrapolation, globally adaptive. */
/* ***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  the routine calculates an approximation result to a given */
/*            definite integral i = integral of f over (a,b), hopefully */
/*            satisfying following claim for accuracy ABS(i-result).le. */
/*            MAX(epsabs,epsrel*ABS(i)). break points of the integration */
/*            interval, where local difficulties of the integrand may */
/*            occur(e.g. singularities,discontinuities),provided by user. */
/* ***description */

/*        computation of a definite integral */
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

/*            npts2  - int */
/*                     number equal to two more than the number of */
/*                     user-supplied break points within the integration */
/*                     range, npts2.ge.2. */
/*                     if npts2.lt.2, the routine will end with ier = 6. */

/*            points - double precision */
/*                     vector of dimension npts2, the first (npts2-2) */
/*                     elements of which are the user provided break */
/*                     points. if these points do not constitute an */
/*                     ascending sequence there will be an automatic */
/*                     sorting. */

/*            epsabs - double precision */
/*                     absolute accuracy requested */
/*            epsrel - double precision */
/*                     relative accuracy requested */
/*                     if  epsabs.le.0 */
/*                     and epsrel.lt.MAX(50*rel.mach.acc.,0.5d-28), */
/*                     the routine will end with ier = 6. */

/*            limit  - int */
/*                     gives an upper bound on the number of subintervals */
/*                     in the partition of (a,b), limit.ge.npts2 */
/*                     if limit.lt.npts2, the routine will end with */
/*                     ier = 6. */

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
/*                     ier.gt.0 abnormal terMINation of the routine. */
/*                             the estimates for integral and error are */
/*                             less reliable. it is assumed that the */
/*                             requested accuracy has not been achieved. */
/*            error messages */
/*                     ier = 1 MAXimum number of subdivisions allowed */
/*                             has been achieved. one can allow more */
/*                             subdivisions by increasing the value of */
/*                             limit (and taking the according dimension */
/*                             adjustments into account). however, if */
/*                             this yields no improvement it is advised */
/*                             to analyze the integrand in order to */
/*                             deterMINe the integration difficulties. if */
/*                             the position of a local difficulty can be */
/*                             deterMINed (i.e. singularity, */
/*                             discontinuity within the interval), it */
/*                             should be supplied to the routine as an */
/*                             element of the vector points. if necessary */
/*                             an appropriate special-purpose integrator */
/*                             must be used, which is designed for */
/*                             handling the type of difficulty involved. */
/*                         = 2 the occurrence of roundoff error is */
/*                             detected, which prevents the requested */
/*                             tolerance from being achieved. */
/*                             the error may be under-estimated. */
/*                         = 3 extremely bad integrand behaviour occurs */
/*                             at some points of the integration */
/*                             interval. */
/*                         = 4 the algorithm does not converge. */
/*                             roundoff error is detected in the */
/*                             extrapolation table. it is presumed that */
/*                             the requested tolerance cannot be */
/*                             achieved, and that the returned result is */
/*                             the best which can be obtained. */
/*                         = 5 the integral is probably divergent, or */
/*                             slowly convergent. it must be noted that */
/*                             divergence can occur with any other value */
/*                             of ier.gt.0. */
/*                         = 6 the input is invalid because */
/*                             npts2.lt.2 or */
/*                             break points are specified outside */
/*                             the integration range or */
/*                             (epsabs.le.0 and */
/*                              epsrel.lt.MAX(50*rel.mach.acc.,0.5d-28)) */
/*                             or limit.lt.npts2. */
/*                             result, abserr, neval, last, rlist(1), */
/*                             and elist(1) are set to zero. alist(1) and */
/*                             blist(1) are set to a and b respectively. */

/*            alist  - double precision */
/*                     vector of dimension at least limit, the first */
/*                      last  elements of which are the left end points */
/*                     of the subintervals in the partition of the given */
/*                     integration range (a,b) */

/*            blist  - double precision */
/*                     vector of dimension at least limit, the first */
/*                      last  elements of which are the right end points */
/*                     of the subintervals in the partition of the given */
/*                     integration range (a,b) */

/*            rlist  - double precision */
/*                     vector of dimension at least limit, the first */
/*                      last  elements of which are the integral */
/*                     approximations on the subintervals */

/*            elist  - double precision */
/*                     vector of dimension at least limit, the first */
/*                      last  elements of which are the moduli of the */
/*                     absolute error estimates on the subintervals */

/*            pts    - double precision */
/*                     vector of dimension at least npts2, containing the */
/*                     integration limits and the break points of the */
/*                     interval in ascending sequence. */

/*            level  - int */
/*                     vector of dimension at least limit, containing the */
/*                     subdivision levels of the subinterval, i.e. if */
/*                     (aa,bb) is a subinterval of (p1,p2) where p1 as */
/*                     well as p2 is a user-provided break point or */
/*                     integration limit, then (aa,bb) has level l if */
/*                     ABS(bb-aa) = ABS(p2-p1)*2**(-l). */

/*            ndin   - int */
/*                     vector of dimension at least npts2, after first */
/*                     integration over the intervals (pts(i)),pts(i+1), */
/*                     i = 0,1, ..., npts2-2, the error estimates over */
/*                     some of the intervals may have been increased */
/*                     artificially, in order to put their subdivision */
/*                     forward. if this happens for the subinterval */
/*                     numbered k, ndin(k) is put to 1, otherwise */
/*                     ndin(k) = 0. */

/*            iord   - int */
/*                     vector of dimension at least limit, the first k */
/*                     elements of which are pointers to the */
/*                     error estimates over the subintervals, */
/*                     such that elist(iord(1)), ..., elist(iord(k)) */
/*                     form a decreasing sequence, with k = last */
/*                     if last.le.(limit/2+2), and k = limit+1-last */
/*                     otherwise */

/*            last   - int */
/*                     number of subintervals actually produced in the */
/*                     subdivisions process */

/* ***references  (none) */
/* ***routines called  dqelg,dqk21,dqpsrt */
/* ***end prologue  dqagpe */




/*            the dimension of rlist2 is deterMINed by the value of */
/*            limexp in subroutine epsalg (rlist2 should be of dimension */
/*            (limexp+2) at least). */


/*            list of major variables */
/*            ----------------------- */

/*           alist     - list of left end points of all subintervals */
/*                       considered up to now */
/*           blist     - list of right end points of all subintervals */
/*                       considered up to now */
/*           rlist(i)  - approximation to the integral over */
/*                       (alist(i),blist(i)) */
/*           rlist2    - array of dimension at least limexp+2 */
/*                       containing the part of the epsilon table which */
/*                       is still needed for further computations */
/*           elist(i)  - error estimate applying to rlist(i) */
/*           MAXerr    - pointer to the interval with largest error */
/*                       estimate */
/*           errMAX    - elist(MAXerr) */
/*           erlast    - error on the interval currently subdivided */
/*                       (before that subdivision has taken place) */
/*           area      - sum of the integrals over the subintervals */
/*           errsum    - sum of the errors over the subintervals */
/*           errbnd    - requested accuracy MAX(epsabs,epsrel* */
/*                       ABS(result)) */
/*           *****1    - variable for the left subinterval */
/*           *****2    - variable for the right subinterval */
/*           last      - index for subdivision */
/*           nres      - number of calls to the extrapolation routine */
/*           numrl2    - number of elements in rlist2. if an appropriate */
/*                       approximation to the compounded integral has */
/*                       been obtained, it is put in rlist2(numrl2) after */
/*                       numrl2 has been increased by one. */
/*           erlarg    - sum of the errors over the intervals larger */
/*                       than the smallest interval considered up to now */
/*           extrap    - int variable denoting that the routine */
/*                       is attempting to perform extrapolation. i.e. */
/*                       before subdividing the smallest interval we */
/*                       try to decrease the value of erlarg. */
/*           noext     - int variable denoting that extrapolation is */
/*                       no longer allowed (true-value) */

/*            machine dependent constants */
/*            --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */
/*           oflow is the largest positive magnitude. */

/* ***first executable statement  dqagpe */
    /* Parameter adjustments */
    --ndin;
    --pts;
    --points;
    --level;
    --iord;
    --elist;
    --rlist;
    --blist;
    --alist__;
    k = 0 ; /* avoid warning */

    /* Function Body */
    epmach = pnl_d1mach(4);

/*            test on validity of parameters */
/*            ----------------------------- */

    *ier = 0;
    *neval = 0;
    *last = 0;
    *result = 0.;
    *abserr = 0.;
    alist__[1] = *a;
    blist[1] = *b;
    rlist[1] = 0.;
    elist[1] = 0.;
    iord[1] = 0;
    level[1] = 0;
    npts = *npts2 - 2;
/* Computing MAX */
    d__1 = epmach * 50.;
    if (*npts2 < 2 || *limit <= npts || (*epsabs <= 0. && *epsrel < MAX(d__1, 5e-29))) {
	*ier = 6;
    }
    if (*ier == 6) {
	goto L999;
    }

/*            if any break points are provided, sort them into an */
/*            ascending sequence. */

    sign = 1.;
    if (*a > *b) {
	sign = -1.;
    }
    pts[1] = MIN(*a,*b);
    if (npts == 0) {
	goto L15;
    }
    i__1 = npts;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pts[i__ + 1] = points[i__];
/* L10: */
    }
L15:
    pts[npts + 2] = MAX(*a,*b);
    nint = npts + 1;
    a1 = pts[1];
    if (npts == 0) {
	goto L40;
    }
    nintp1 = nint + 1;
    i__1 = nint;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ip1 = i__ + 1;
	i__2 = nintp1;
	for (j = ip1; j <= i__2; ++j) {
	    if (pts[i__] <= pts[j]) {
		goto L20;
	    }
	    temp = pts[i__];
	    pts[i__] = pts[j];
	    pts[j] = temp;
L20:
	    ;
	}
    }
    if (pts[1] != MIN(*a,*b) || pts[nintp1] != MAX(*a,*b)) {
	*ier = 6;
    }
    if (*ier == 6) {
	goto L999;
    }

/*            compute first integral and error approximations. */
/*            ------------------------------------------------ */

L40:
    resabs = 0.;
    i__2 = nint;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b1 = pts[i__ + 1];
	pnl_dqk21((PnlFunc *)f, &a1, &b1, &area1, &error1, &defabs, &resa);
	*abserr += error1;
	*result += area1;
	ndin[i__] = 0;
	if (error1 == resa && error1 != 0.) {
	    ndin[i__] = 1;
	}
	resabs += defabs;
	level[i__] = 0;
	elist[i__] = error1;
	alist__[i__] = a1;
	blist[i__] = b1;
	rlist[i__] = area1;
	iord[i__] = i__;
	a1 = b1;
/* L50: */
    }
    errsum = 0.;
    i__2 = nint;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (ndin[i__] == 1) {
	    elist[i__] = *abserr;
	}
	errsum += elist[i__];
/* L55: */
    }

/*           test on accuracy. */

    *last = nint;
    *neval = nint * 21;
    dres = ABS(*result);
/* Computing MAX */
    d__1 = *epsabs, d__2 = *epsrel * dres;
    errbnd = MAX(d__1,d__2);
    if (*abserr <= epmach * 100. * resabs && *abserr > errbnd) {
	*ier = 2;
    }
    if (nint == 1) {
	goto L80;
    }
    i__2 = npts;
    for (i__ = 1; i__ <= i__2; ++i__) {
	jlow = i__ + 1;
	ind1 = iord[i__];
	i__1 = nint;
	for (j = jlow; j <= i__1; ++j) {
	    ind2 = iord[j];
	    if (elist[ind1] > elist[ind2]) {
		goto L60;
	    }
	    ind1 = ind2;
	    k = j;
L60:
	    ;
	}
	if (ind1 == iord[i__]) {
	    goto L70;
	}
	iord[k] = iord[i__];
	iord[i__] = ind1;
L70:
	;
    }
    if (*limit < *npts2) {
	*ier = 1;
    }
L80:
    if (*ier != 0 || *abserr <= errbnd) {
	goto L210;
    }

/*           initialization */
/*           -------------- */

    rlist2[0] = *result;
    MAXerr = iord[1];
    errMAX = elist[MAXerr];
    area = *result;
    nrMAX = 1;
    nres = 0;
    numrl2 = 1;
    ktMIN = 0;
    extrap = FALSE;
    noext = FALSE;
    erlarg = errsum;
    ertest = errbnd;
    levMAX = 1;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ierro = 0;
    uflow = pnl_d1mach(1);
    oflow = pnl_d1mach(2);
    *abserr = oflow;
    ksgn = -1;
    if (dres >= (1. - epmach * 50.) * resabs) {
	ksgn = 1;
    }

/*           main do-loop */
/*           ------------ */

    i__2 = *limit;
    for (*last = *npts2; *last <= i__2; ++(*last)) {

/*           bisect the subinterval with the nrMAX-th largest error */
/*           estimate. */

	levcur = level[MAXerr] + 1;
	a1 = alist__[MAXerr];
	b1 = (alist__[MAXerr] + blist[MAXerr]) * .5;
	a2 = b1;
	b2 = blist[MAXerr];
	erlast = errMAX;
	pnl_dqk21((PnlFunc *)f, &a1, &b1, &area1, &error1, &resa, &defab1);
	pnl_dqk21((PnlFunc *)f, &a2, &b2, &area2, &error2, &resa, &defab2);

/*           improve previous approximations to integral */
/*           and error and test for accuracy. */

	*neval += 42;
	area12 = area1 + area2;
	erro12 = error1 + error2;
	errsum = errsum + erro12 - errMAX;
	area = area + area12 - rlist[MAXerr];
	if (defab1 == error1 || defab2 == error2) {
	    goto L95;
	}
	if ((d__1 = rlist[MAXerr] - area12, ABS(d__1)) > ABS(area12) * 1e-5 ||
		 erro12 < errMAX * .99) {
	    goto L90;
	}
	if (extrap) {
	    ++iroff2;
	}
	if (! extrap) {
	    ++iroff1;
	}
L90:
	if (*last > 10 && erro12 > errMAX) {
	    ++iroff3;
	}
L95:
	level[MAXerr] = levcur;
	level[*last] = levcur;
	rlist[MAXerr] = area1;
	rlist[*last] = area2;
/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * ABS(area);
	errbnd = MAX(d__1,d__2);

/*           test for roundoff error and eventually set error flag. */

	if (iroff1 + iroff2 >= 10 || iroff3 >= 20) {
	    *ier = 2;
	}
	if (iroff2 >= 5) {
	    ierro = 3;
	}

/*           set error flag in the case that the number of */
/*           subintervals equals limit. */

	if (*last == *limit) {
	    *ier = 1;
	}

/*           set error flag in the case of bad integrand behaviour */
/*           at a point of the integration range */

/* Computing MAX */
	d__1 = ABS(a1), d__2 = ABS(b2);
	if (MAX(d__1,d__2) <= (epmach * 100. + 1.) * (ABS(a2) + uflow * 1e3))
		{
	    *ier = 4;
	}

/*           append the newly-created intervals to the list. */

	if (error2 > error1) {
	    goto L100;
	}
	alist__[*last] = a2;
	blist[MAXerr] = b1;
	blist[*last] = b2;
	elist[MAXerr] = error1;
	elist[*last] = error2;
	goto L110;
L100:
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

L110:
	pnl_dqpsrt(limit, last, &MAXerr, &errMAX, &elist[1], &iord[1], &nrMAX);
/* ***jump out of do-loop */
	if (errsum <= errbnd) {
	    goto L190;
	}
/* ***jump out of do-loop */
	if (*ier != 0) {
	    goto L170;
	}
	if (noext) {
	    goto L160;
	}
	erlarg -= erlast;
	if (levcur + 1 <= levMAX) {
	    erlarg += erro12;
	}
	if (extrap) {
	    goto L120;
	}

/*           test whether the interval to be bisected next is the */
/*           smallest interval. */

	if (level[MAXerr] + 1 <= levMAX) {
	    goto L160;
	}
	extrap = TRUE;
	nrMAX = 2;
L120:
	if (ierro == 3 || erlarg <= ertest) {
	    goto L140;
	}

/*           the smallest interval has the largest error. */
/*           before bisecting decrease the sum of the errors over */
/*           the larger intervals (erlarg) and perform extrapolation. */

	id = nrMAX;
	jupbnd = *last;
	if (*last > *limit / 2 + 2) {
	    jupbnd = *limit + 3 - *last;
	}
	i__1 = jupbnd;
	for (k = id; k <= i__1; ++k) {
	    MAXerr = iord[nrMAX];
	    errMAX = elist[MAXerr];
/* ***jump out of do-loop */
	    if (level[MAXerr] + 1 <= levMAX) {
		goto L160;
	    }
	    ++nrMAX;
/* L130: */
	}

/*           perform extrapolation. */

L140:
	++numrl2;
	rlist2[numrl2 - 1] = area;
	if (numrl2 <= 2) {
	    goto L155;
	}
	pnl_dqelg(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
	++ktMIN;
	if (ktMIN > 5 && *abserr < errsum * .001) {
	    *ier = 5;
	}
	if (abseps >= *abserr) {
	    goto L150;
	}
	ktMIN = 0;
	*abserr = abseps;
	*result = reseps;
	correc = erlarg;
/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * ABS(reseps);
	ertest = MAX(d__1,d__2);
/* ***jump out of do-loop */
	if (*abserr < ertest) {
	    goto L170;
	}

/*           prepare bisection of the smallest interval. */

L150:
	if (numrl2 == 1) {
	    noext = TRUE;
	}
	if (*ier >= 5) {
	    goto L170;
	}
L155:
	MAXerr = iord[1];
	errMAX = elist[MAXerr];
	nrMAX = 1;
	extrap = FALSE;
	++levMAX;
	erlarg = errsum;
L160:
	;
    }

/*           set the final result. */
/*           --------------------- */


L170:
    if (*abserr == oflow) {
	goto L190;
    }
    if (*ier + ierro == 0) {
	goto L180;
    }
    if (ierro == 3) {
	*abserr += correc;
    }
    if (*ier == 0) {
	*ier = 3;
    }
    if (*result != 0. && area != 0.) {
	goto L175;
    }
    if (*abserr > errsum) {
	goto L190;
    }
    if (area == 0.) {
	goto L210;
    }
    goto L180;
L175:
    if (*abserr / ABS(*result) > errsum / ABS(area)) {
	goto L190;
    }

/*           test on divergence. */

L180:
/* Computing MAX */
    d__1 = ABS(*result), d__2 = ABS(area);
    if (ksgn == -1 && MAX(d__1,d__2) <= resabs * .01) {
	goto L210;
    }
    if (.01 > *result / area || *result / area > 100. || errsum > ABS(area)) {
	*ier = 6;
    }
    goto L210;

/*           compute global integral sum. */

L190:
    *result = 0.;
    i__2 = *last;
    for (k = 1; k <= i__2; ++k) {
	*result += rlist[k];
/* L200: */
    }
    *abserr = errsum;
L210:
    if (*ier > 2) {
	--(*ier);
    }
    *result *= sign;
L999:
    return 0;
} /* pnl_dqagpe */

