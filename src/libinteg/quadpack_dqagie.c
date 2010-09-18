
/*
 * This file contains the routine dqagie from the quadpack distribution
 * available on www.netlib.org
 *
 * This routine has been translated from Fortran using f2c and then slightly
 * modified by Jérôme Lelong for integration into Pnl. All Quadpack names have
 * been prefixed by pnl to avoid possible clash in names. Calls to d1mach have
 * been replaced by calls to dlamch coming from Lapack
 */

#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"

extern int pnl_dqelg(int *, double *, double *, double *, double *, int *);
extern int pnl_dqk15i(PnlFunc *, double * , int *, double *, double *, double *, double *, double *, double *);
extern int pnl_dqpsrt(int *, int *, int *, double *, double *, int *, int *);

int pnl_dqagie(PnlFunc *f, double *bound, int *inf, 
	double *epsabs, double *epsrel, int *limit, double *
	result, double *abserr, int *neval, int *ier, double *
	alist__, double *blist, double *rlist, double *elist, 
	int *iord, int *last)
{
    /* System generated locals */
    int i__1, i__2;
    double d__1, d__2, c_b4, c_b5;

    /* Local variables */
    int k;
    double a1, a2, b1, b2;
    int id;
    double area, dres;
    int ksgn;
    double boun;
    int nres;
    double area1, area2, area12;
    double small=0, erro12;
    int ierro;
    double defab1, defab2;
    int ktMIN, nrMAX;
    double oflow, uflow;
    
    int noext;
    int iroff1, iroff2, iroff3;
    double res3la[3], error1, error2, rlist2[52];
    int numrl2;
    double defabs, epmach, erlarg=0, abseps, correc=0, errbnd, resabs;
    int jupbnd;
    double erlast, errMAX;
    int MAXerr;
    double reseps;
    int extrap;
    double ertest=0, errsum;

/* ***begin prologue  dqagie */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a3a1,h2a4a1 */
/* ***keywords  automatic integrator, infinite intervals, */
/*             general-purpose, transformation, extrapolation, */
/*             globally adaptive */
/* ***author  piessens,robert,appl. math & progr. div - k.u.leuven */
/*           de doncker,elise,appl. math & progr. div - k.u.leuven */
/* ***purpose  the routine calculates an approximation result to a given */
/*            integral   i = integral of f over (bound,+infinity) */
/*            or i = integral of f over (-infinity,bound) */
/*            or i = integral of f over (-infinity,+infinity), */
/*            hopefully satisfying following claim for accuracy */
/*            ABS(i-result).le.MAX(epsabs,epsrel*ABS(i)) */
/* ***description */

/* integration over infinite intervals */
/* standard fortran subroutine */

/*            f      - double precision */
/*                     function subprogram defining the integrand */
/*                     function f(x). the actual name for f needs to be */
/*                     declared e x t e r n a l in the driver program. */

/*            bound  - double precision */
/*                     finite bound of integration range */
/*                     (has no meaning if interval is doubly-infinite) */

/*            inf    - double precision */
/*                     indicating the kind of integration range involved */
/*                     inf = 1 corresponds to  (bound,+infinity), */
/*                     inf = -1            to  (-infinity,bound), */
/*                     inf = 2             to (-infinity,+infinity). */

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
/*                   - ier.gt.0 abnormal terMINation of the routine. the */
/*                             estimates for result and error are less */
/*                             reliable. it is assumed that the requested */
/*                             accuracy has not been achieved. */
/*            error messages */
/*                     ier = 1 MAXimum number of subdivisions allowed */
/*                             has been achieved. one can allow more */
/*                             subdivisions by increasing the value of */
/*                             limit (and taking the according dimension */
/*                             adjustments into account). however,if */
/*                             this yields no improvement it is advised */
/*                             to analyze the integrand in order to */
/*                             deterMINe the integration difficulties. */
/*                             if the position of a local difficulty can */
/*                             be deterMINed (e.g. singularity, */
/*                             discontinuity within the interval) one */
/*                             will probably gain from splitting up the */
/*                             interval at this point and calling the */
/*                             integrator on the subranges. if possible, */
/*                             an appropriate special-purpose integrator */
/*                             should be used, which is designed for */
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
/*                             extrapolation table. */
/*                             it is assumed that the requested tolerance */
/*                             cannot be achieved, and that the returned */
/*                             result is the best which can be obtained. */
/*                         = 5 the integral is probably divergent, or */
/*                             slowly convergent. it must be noted that */
/*                             divergence can occur with any other value */
/*                             of ier. */
/*                         = 6 the input is invalid, because */
/*                             (epsabs.le.0 and */
/*                              epsrel.lt.MAX(50*rel.mach.acc.,0.5d-28), */
/*                             result, abserr, neval, last, rlist(1), */
/*                             elist(1) and iord(1) are set to zero. */
/*                             alist(1) and blist(1) are set to 0 */
/*                             and 1 respectively. */

/*            alist  - double precision */
/*                     vector of dimension at least limit, the first */
/*                      last  elements of which are the left */
/*                     end points of the subintervals in the partition */
/*                     of the transformed integration range (0,1). */

/*            blist  - double precision */
/*                     vector of dimension at least limit, the first */
/*                      last  elements of which are the right */
/*                     end points of the subintervals in the partition */
/*                     of the transformed integration range (0,1). */

/*            rlist  - double precision */
/*                     vector of dimension at least limit, the first */
/*                      last  elements of which are the integral */
/*                     approximations on the subintervals */

/*            elist  - double precision */
/*                     vector of dimension at least limit,  the first */
/*                     last elements of which are the moduli of the */
/*                     absolute error estimates on the subintervals */

/*            iord   - int */
/*                     vector of dimension limit, the first k */
/*                     elements of which are pointers to the */
/*                     error estimates over the subintervals, */
/*                     such that elist(iord(1)), ..., elist(iord(k)) */
/*                     form a decreasing sequence, with k = last */
/*                     if last.le.(limit/2+2), and k = limit+1-last */
/*                     otherwise */

/*            last   - int */
/*                     number of subintervals actually produced */
/*                     in the subdivision process */

/* ***references  (none) */
/* ***routines called  dqelg,dqk15i,dqpsrt */
/* ***end prologue  dqagie */



/*            the dimension of rlist2 is deterMINed by the value of */
/*            limexp in subroutine dqelg. */


/*            list of major variables */
/*            ----------------------- */

/*           alist     - list of left end points of all subintervals */
/*                       considered up to now */
/*           blist     - list of right end points of all subintervals */
/*                       considered up to now */
/*           rlist(i)  - approximation to the integral over */
/*                       (alist(i),blist(i)) */
/*           rlist2    - array of dimension at least (limexp+2), */
/*                       containing the part of the epsilon table */
/*                       wich is still needed for further computations */
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
/*           numrl2    - number of elements currently in rlist2. if an */
/*                       appropriate approximation to the compounded */
/*                       integral has been obtained, it is put in */
/*                       rlist2(numrl2) after numrl2 has been increased */
/*                       by one. */
/*           small     - length of the smallest interval considered up */
/*                       to now, multiplied by 1.5 */
/*           erlarg    - sum of the errors over the intervals larger */
/*                       than the smallest interval considered up to now */
/*           extrap    - int variable denoting that the routine */
/*                       is attempting to perform extrapolation. i.e. */
/*                       before subdividing the smallest interval we */
/*                       try to decrease the value of erlarg. */
/*           noext     - int variable denoting that extrapolation */
/*                       is no longer allowed (true-value) */

/*            machine dependent constants */
/*            --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */
/*           oflow is the largest positive magnitude. */

/* ***first executable statement  dqagie */
    /* Parameter adjustments */
    --iord;
    --elist;
    --rlist;
    --blist;
    --alist__;
    c_b4 = 0;
    c_b5 = 1;

    /* Function Body */
    epmach = pnl_d1mach(4);

/*           test on validity of parameters */
/*           ----------------------------- */

    *ier = 0;
    *neval = 0;
    *last = 0;
    *result = 0.;
    *abserr = 0.;
    alist__[1] = 0.;
    blist[1] = 1.;
    rlist[1] = 0.;
    elist[1] = 0.;
    iord[1] = 0;
/* Computing MAX */
    d__1 = epmach * 50.;
    if (*epsabs <= 0. && *epsrel < MAX(d__1,5e-29)) {
	*ier = 6;
    }
    if (*ier == 6) {
	goto L999;
    }


/*           first approximation to the integral */
/*           ----------------------------------- */

/*           deterMINe the interval to be mapped onto (0,1). */
/*           if inf = 2 the integral is computed as i = i1+i2, where */
/*           i1 = integral of f over (-infinity,0), */
/*           i2 = integral of f over (0,+infinity). */

    boun = *bound;
    if (*inf == 2) {
	boun = 0.;
    }
    pnl_dqk15i((PnlFunc *)f, &boun, inf, &c_b4, &c_b5, result, abserr, &defabs, &
	    resabs);

/*           test on accuracy */

    *last = 1;
    rlist[1] = *result;
    elist[1] = *abserr;
    iord[1] = 1;
    dres = ABS(*result);
/* Computing MAX */
    d__1 = *epsabs, d__2 = *epsrel * dres;
    errbnd = MAX(d__1,d__2);
    if (*abserr <= epmach * 100. * defabs && *abserr > errbnd) {
	*ier = 2;
    }
    if (*limit == 1) {
	*ier = 1;
    }
    if (*ier != 0 || (*abserr <= errbnd && *abserr != resabs) || *abserr == 0.) 
	    {
	goto L130;
    }

/*           initialization */
/*           -------------- */

    uflow = pnl_d1mach(1);
    oflow = pnl_d1mach(2);
    rlist2[0] = *result;
    errMAX = *abserr;
    MAXerr = 1;
    area = *result;
    errsum = *abserr;
    *abserr = oflow;
    nrMAX = 1;
    nres = 0;
    ktMIN = 0;
    numrl2 = 2;
    extrap = FALSE;
    noext = FALSE;
    ierro = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ksgn = -1;
    if (dres >= (1. - epmach * 50.) * defabs) {
	ksgn = 1;
    }

/*           main do-loop */
/*           ------------ */

    i__1 = *limit;
    for (*last = 2; *last <= i__1; ++(*last)) {

/*           bisect the subinterval with nrMAX-th largest error estimate. */

	a1 = alist__[MAXerr];
	b1 = (alist__[MAXerr] + blist[MAXerr]) * .5;
	a2 = b1;
	b2 = blist[MAXerr];
	erlast = errMAX;
	pnl_dqk15i((PnlFunc *)f, &boun, inf, &a1, &b1, &area1, &error1, &resabs, &
		defab1);
	pnl_dqk15i((PnlFunc *)f, &boun, inf, &a2, &b2, &area2, &error2, &resabs, &
		defab2);

/*           improve previous approximations to integral */
/*           and error and test for accuracy. */

	area12 = area1 + area2;
	erro12 = error1 + error2;
	errsum = errsum + erro12 - errMAX;
	area = area + area12 - rlist[MAXerr];
	if (defab1 == error1 || defab2 == error2) {
	    goto L15;
	}
	if ((d__1 = rlist[MAXerr] - area12, ABS(d__1)) > ABS(area12) * 1e-5 ||
		 erro12 < errMAX * .99) {
	    goto L10;
	}
	if (extrap) {
	    ++iroff2;
	}
	if (! extrap) {
	    ++iroff1;
	}
L10:
	if (*last > 10 && erro12 > errMAX) {
	    ++iroff3;
	}
L15:
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
/*           at some points of the integration range. */

/* Computing MAX */
	d__1 = ABS(a1), d__2 = ABS(b2);
	if (MAX(d__1,d__2) <= (epmach * 100. + 1.) * (ABS(a2) + uflow * 1e3)) 
		{
	    *ier = 4;
	}

/*           append the newly-created intervals to the list. */

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
	if (errsum <= errbnd) {
	    goto L115;
	}
	if (*ier != 0) {
	    goto L100;
	}
	if (*last == 2) {
	    goto L80;
	}
	if (noext) {
	    goto L90;
	}
	erlarg -= erlast;
	if ((d__1 = b1 - a1, ABS(d__1)) > small) {
	    erlarg += erro12;
	}
	if (extrap) {
	    goto L40;
	}

/*           test whether the interval to be bisected next is the */
/*           smallest interval. */

	if ((d__1 = blist[MAXerr] - alist__[MAXerr], ABS(d__1)) > small) {
	    goto L90;
	}
	extrap = TRUE;
	nrMAX = 2;
L40:
	if (ierro == 3 || erlarg <= ertest) {
	    goto L60;
	}

/*           the smallest interval has the largest error. */
/*           before bisecting decrease the sum of the errors over the */
/*           larger intervals (erlarg) and perform extrapolation. */

	id = nrMAX;
	jupbnd = *last;
	if (*last > *limit / 2 + 2) {
	    jupbnd = *limit + 3 - *last;
	}
	i__2 = jupbnd;
	for (k = id; k <= i__2; ++k) {
	    MAXerr = iord[nrMAX];
	    errMAX = elist[MAXerr];
	    if ((d__1 = blist[MAXerr] - alist__[MAXerr], ABS(d__1)) > small) {
		goto L90;
	    }
	    ++nrMAX;
/* L50: */
	}

/*           perform extrapolation. */

L60:
	++numrl2;
	rlist2[numrl2 - 1] = area;
	pnl_dqelg(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
	++ktMIN;
	if (ktMIN > 5 && *abserr < errsum * .001) {
	    *ier = 5;
	}
	if (abseps >= *abserr) {
	    goto L70;
	}
	ktMIN = 0;
	*abserr = abseps;
	*result = reseps;
	correc = erlarg;
/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * ABS(reseps);
	ertest = MAX(d__1,d__2);
	if (*abserr <= ertest) {
	    goto L100;
	}

/*            prepare bisection of the smallest interval. */

L70:
	if (numrl2 == 1) {
	    noext = TRUE;
	}
	if (*ier == 5) {
	    goto L100;
	}
	MAXerr = iord[1];
	errMAX = elist[MAXerr];
	nrMAX = 1;
	extrap = FALSE;
	small *= .5;
	erlarg = errsum;
	goto L90;
L80:
	small = .375;
	erlarg = errsum;
	ertest = errbnd;
	rlist2[1] = area;
L90:
	;
    }

/*           set final result and error estimate. */
/*           ------------------------------------ */

L100:
    if (*abserr == oflow) {
	goto L115;
    }
    if (*ier + ierro == 0) {
	goto L110;
    }
    if (ierro == 3) {
	*abserr += correc;
    }
    if (*ier == 0) {
	*ier = 3;
    }
    if (*result != 0. && area != 0.) {
	goto L105;
    }
    if (*abserr > errsum) {
	goto L115;
    }
    if (area == 0.) {
	goto L130;
    }
    goto L110;
L105:
    if (*abserr / ABS(*result) > errsum / ABS(area)) {
	goto L115;
    }

/*           test on divergence */

L110:
/* Computing MAX */
    d__1 = ABS(*result), d__2 = ABS(area);
    if (ksgn == -1 && MAX(d__1,d__2) <= defabs * .01) {
	goto L130;
    }
    if (.01 > *result / area || *result / area > 100. || errsum > ABS(area)) {
	*ier = 6;
    }
    goto L130;

/*           compute global integral sum. */

L115:
    *result = 0.;
    i__1 = *last;
    for (k = 1; k <= i__1; ++k) {
	*result += rlist[k];
/* L120: */
    }
    *abserr = errsum;
L130:
    *neval = *last * 30 - 15;
    if (*inf == 2) {
	*neval <<= 1;
    }
    if (*ier > 2) {
	--(*ier);
    }
L999:
    return 0;
} /* dqagie */

