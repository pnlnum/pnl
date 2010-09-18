/*
 * This file contains the routine dqawoe from the quadpack distribution
 * available on www.netlib.org
 *
 * This routine has been translated from Fortran using f2c and then slightly
 * modified by Jérôme Lelong for integration into Pnl. All Quadpack names have
 * been prefixed by pnl to avoid possible clash in names. Calls to d1mach have
 * been replaced by calls to dlamch coming from Lapack
 */

#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"

extern  int pnl_dqc25f(PnlFunc *, double *, double *, 
                    double *, int *, int *, int *, int *, 
                    double *, double *, int *, double *, double *,
                    int *, double *), pnl_dqelg(int *, double *, 
                                             double *, double *, double *, int *);
extern  int pnl_dqpsrt(int *, int *, int *, 
                    double *, double *, int *, int *);

int pnl_dqawoe(PnlFunc * f, double *a, double *b, double 
	*omega, int *integr, double *epsabs, double *epsrel, 
	int *limit, int *icall, int *MAXp1, double *result, 
	double *abserr, int *neval, int *ier, int *last, 
	double *alist__, double *blist, double *rlist, double 
	*elist, int *iord, int *nnlog, int *momcom, double *
	chebmo)
{
    /* System generated locals */
    int chebmo_dim1, chebmo_offset, i__1, i__2;
    double d__1, d__2;
    int c__1, c__0;

    /* Local variables */
    int k;
    double a1, a2, b1, b2;
    int id, nev;
    double area, dres;
    int ksgn, nres;
    double area1, area2, area12;
    double small, erro12, width, defab1, defab2;
    int ierro, ktMIN;
    double oflow;
    int nrMAX, nrmom;
    double uflow;
    
    int noext;
    int iroff1, iroff2, iroff3;
    double res3la[3], error1, error2, rlist2[52];
    int numrl2;
    double defabs, domega, epmach, erlarg, abseps, correc, errbnd, resabs;
    int jupbnd;
    int extall;
    double erlast, errMAX;
    int MAXerr;
    double reseps;
    int extrap;
    double ertest, errsum;

/* ***begin prologue  dqawoe */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a2a1 */
/* ***keywords  automatic integrator, special-purpose, */
/*             integrand with oscillatory cos or sin factor, */
/*             clenshaw-curtis method, (end point) singularities, */
/*             extrapolation, globally adaptive */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  the routine calculates an approximation result to a given */
/*            definite integral */
/*            i = integral of f(x)*w(x) over (a,b) */
/*            where w(x) = cos(omega*x) or w(x)=sin(omega*x), */
/*            hopefully satisfying following claim for accuracy */
/*            ABS(i-result).le.MAX(epsabs,epsrel*ABS(i)). */
/* ***description */

/*        computation of oscillatory integrals */
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

/*            omega  - double precision */
/*                     parameter in the integrand weight function */

/*            integr - int */
/*                     indicates which of the weight functions is to be */
/*                     used */
/*                     integr = 1      w(x) = cos(omega*x) */
/*                     integr = 2      w(x) = sin(omega*x) */
/*                     if integr.ne.1 and integr.ne.2, the routine */
/*                     will end with ier = 6. */

/*            epsabs - double precision */
/*                     absolute accuracy requested */
/*            epsrel - double precision */
/*                     relative accuracy requested */
/*                     if  epsabs.le.0 */
/*                     and epsrel.lt.MAX(50*rel.mach.acc.,0.5d-28), */
/*                     the routine will end with ier = 6. */

/*            limit  - int */
/*                     gives an upper bound on the number of subdivisions */
/*                     in the partition of (a,b), limit.ge.1. */

/*            icall  - int */
/*                     if dqawoe is to be used only once, icall must */
/*                     be set to 1.  assume that during this call, the */
/*                     chebyshev moments (for clenshaw-curtis integration */
/*                     of degree 24) have been computed for intervals of */
/*                     lenghts (ABS(b-a))*2**(-l), l=0,1,2,...momcom-1. */
/*                     if icall.gt.1 this means that dqawoe has been */
/*                     called twice or more on intervals of the same */
/*                     length ABS(b-a). the chebyshev moments already */
/*                     computed are then re-used in subsequent calls. */
/*                     if icall.lt.1, the routine will end with ier = 6. */

/*            MAXp1  - int */
/*                     gives an upper bound on the number of chebyshev */
/*                     moments which can be stored, i.e. for the */
/*                     intervals of lenghts ABS(b-a)*2**(-l), */
/*                     l=0,1, ..., MAXp1-2, MAXp1.ge.1. */
/*                     if MAXp1.lt.1, the routine will end with ier = 6. */

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
/*                             routine. it is assumed that the */
/*                             requested accuracy has been achieved. */
/*                   - ier.gt.0 abnormal terMINation of the routine. */
/*                             the estimates for integral and error are */
/*                             less reliable. it is assumed that the */
/*                             requested accuracy has not been achieved. */
/*            error messages */
/*                     ier = 1 MAXimum number of subdivisions allowed */
/*                             has been achieved. one can allow more */
/*                             subdivisions by increasing the value of */
/*                             limit (and taking according dimension */
/*                             adjustments into account). however, if */
/*                             this yields no improvement it is advised */
/*                             to analyze the integrand, in order to */
/*                             deterMINe the integration difficulties. */
/*                             if the position of a local difficulty can */
/*                             be deterMINed (e.g. singularity, */
/*                             discontinuity within the interval) one */
/*                             will probably gain from splitting up the */
/*                             interval at this point and calling the */
/*                             integrator on the subranges. if possible, */
/*                             an appropriate special-purpose integrator */
/*                             should be used which is designed for */
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
/*                             it is presumed that the requested */
/*                             tolerance cannot be achieved due to */
/*                             roundoff in the extrapolation table, */
/*                             and that the returned result is the */
/*                             best which can be obtained. */
/*                         = 5 the integral is probably divergent, or */
/*                             slowly convergent. it must be noted that */
/*                             divergence can occur with any other value */
/*                             of ier.gt.0. */
/*                         = 6 the input is invalid, because */
/*                             (epsabs.le.0 and */
/*                              epsrel.lt.MAX(50*rel.mach.acc.,0.5d-28)) */
/*                             or (integr.ne.1 and integr.ne.2) or */
/*                             icall.lt.1 or MAXp1.lt.1. */
/*                             result, abserr, neval, last, rlist(1), */
/*                             elist(1), iord(1) and nnlog(1) are set */
/*                             to zero. alist(1) and blist(1) are set */
/*                             to a and b respectively. */

/*            last  -  int */
/*                     on return, last equals the number of */
/*                     subintervals produces in the subdivision */
/*                     process, which deterMINes the number of */
/*                     significant elements actually in the */
/*                     work arrays. */
/*            alist  - double precision */
/*                     vector of dimension at least limit, the first */
/*                      last  elements of which are the left */
/*                     end points of the subintervals in the partition */
/*                     of the given integration range (a,b) */

/*            blist  - double precision */
/*                     vector of dimension at least limit, the first */
/*                      last  elements of which are the right */
/*                     end points of the subintervals in the partition */
/*                     of the given integration range (a,b) */

/*            rlist  - double precision */
/*                     vector of dimension at least limit, the first */
/*                      last  elements of which are the integral */
/*                     approximations on the subintervals */

/*            elist  - double precision */
/*                     vector of dimension at least limit, the first */
/*                      last  elements of which are the moduli of the */
/*                     absolute error estimates on the subintervals */

/*            iord   - int */
/*                     vector of dimension at least limit, the first k */
/*                     elements of which are pointers to the error */
/*                     estimates over the subintervals, */
/*                     such that elist(iord(1)), ..., */
/*                     elist(iord(k)) form a decreasing sequence, with */
/*                     k = last if last.le.(limit/2+2), and */
/*                     k = limit+1-last otherwise. */

/*            nnlog  - int */
/*                     vector of dimension at least limit, containing the */
/*                     subdivision levels of the subintervals, i.e. */
/*                     iwork(i) = l means that the subinterval */
/*                     numbered i is of length ABS(b-a)*2**(1-l) */

/*         on entry and return */
/*            momcom - int */
/*                     indicating that the chebyshev moments */
/*                     have been computed for intervals of lengths */
/*                     (ABS(b-a))*2**(-l), l=0,1,2, ..., momcom-1, */
/*                     momcom.lt.MAXp1 */

/*            chebmo - double precision */
/*                     array of dimension (MAXp1,25) containing the */
/*                     chebyshev moments */

/* ***references  (none) */
/* ***routines called  dqc25f,dqelg,dqpsrt */
/* ***end prologue  dqawoe */




/*            the dimension of rlist2 is deterMINed by  the value of */
/*            limexp in subroutine dqelg (rlist2 should be of */
/*            dimension (limexp+2) at least). */

/*            list of major variables */
/*            ----------------------- */

/*           alist     - list of left end points of all subintervals */
/*                       considered up to now */
/*           blist     - list of right end points of all subintervals */
/*                       considered up to now */
/*           rlist(i)  - approximation to the integral over */
/*                       (alist(i),blist(i)) */
/*           rlist2    - array of dimension at least limexp+2 */
/*                       containing the part of the epsilon table */
/*                       which is still needed for further computations */
/*           elist(i)  - error estimate applying to rlist(i) */
/*           MAXerr    - pointer to the interval with largest */
/*                       error estimate */
/*           errMAX    - elist(MAXerr) */
/*           erlast    - error on the interval currently subdivided */
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
/*                       been obtained it is put in rlist2(numrl2) after */
/*                       numrl2 has been increased by one */
/*           small     - length of the smallest interval considered */
/*                       up to now, multiplied by 1.5 */
/*           erlarg    - sum of the errors over the intervals larger */
/*                       than the smallest interval considered up to now */
/*           extrap    - int variable denoting that the routine is */
/*                       attempting to perform extrapolation, i.e. before */
/*                       subdividing the smallest interval we try to */
/*                       decrease the value of erlarg */
/*           noext     - int variable denoting that extrapolation */
/*                       is no longer allowed (true  value) */

/*            machine dependent constants */
/*            --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */
/*           oflow is the largest positive magnitude. */

/* ***first executable statement  dqawoe */
    /* Parameter adjustments */
    --nnlog;
    --iord;
    --elist;
    --rlist;
    --blist;
    --alist__;
    chebmo_dim1 = *MAXp1;
    chebmo_offset = 1 + chebmo_dim1;
    chebmo -= chebmo_offset;
    c__1 = 1;
    c__0 = 0;

    /* Function Body */
    epmach = pnl_d1mach(4);

/*         test on validity of parameters */
/*         ------------------------------ */

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
    nnlog[1] = 0;
/* Computing MAX */
    d__1 = epmach * 50.;
    if ((*integr != 1 && *integr != 2) || (*epsabs <= 0. && *epsrel < MAX(d__1, 5e-29))
        || *icall < 1 || *MAXp1 < 1) {
	*ier = 6;
    }
    if (*ier == 6) {
	goto L999;
    }

/*           first approximation to the integral */
/*           ----------------------------------- */

    domega = ABS(*omega);
    nrmom = 0;
    if (*icall > 1) {
	goto L5;
    }
    *momcom = 0;
L5:
    pnl_dqc25f((PnlFunc *)f, a, b, &domega, integr, &nrmom, MAXp1, &c__0, result, 
	    abserr, neval, &defabs, &resabs, momcom, &chebmo[chebmo_offset]);

/*           test on accuracy. */

    dres = ABS(*result);
/* Computing MAX */
    d__1 = *epsabs, d__2 = *epsrel * dres;
    errbnd = MAX(d__1,d__2);
    rlist[1] = *result;
    elist[1] = *abserr;
    iord[1] = 1;
    if (*abserr <= epmach * 100. * defabs && *abserr > errbnd) {
	*ier = 2;
    }
    if (*limit == 1) {
	*ier = 1;
    }
    if (*ier != 0 || *abserr <= errbnd) {
	goto L200;
    }

/*           initializations */
/*           --------------- */

    uflow = pnl_d1mach(1);
    oflow = pnl_d1mach(2);
    errMAX = *abserr;
    MAXerr = 1;
    area = *result;
    errsum = *abserr;
    *abserr = oflow;
    nrMAX = 1;
    extrap = FALSE;
    noext = FALSE;
    ierro = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ktMIN = 0;
    small = (d__1 = *b - *a, ABS(d__1)) * .75;
    nres = 0;
    numrl2 = 0;
    extall = FALSE;
    if ((d__1 = *b - *a, ABS(d__1)) * .5 * domega > 2.) {
	goto L10;
    }
    numrl2 = 1;
    extall = TRUE;
    rlist2[0] = *result;
L10:
    if ((d__1 = *b - *a, ABS(d__1)) * .25 * domega <= 2.) {
	extall = TRUE;
    }
    ksgn = -1;
    if (dres >= (1. - epmach * 50.) * defabs) {
	ksgn = 1;
    }

/*           main do-loop */
/*           ------------ */

    i__1 = *limit;
    for (*last = 2; *last <= i__1; ++(*last)) {

/*           bisect the subinterval with the nrMAX-th largest */
/*           error estimate. */

	nrmom = nnlog[MAXerr] + 1;
	a1 = alist__[MAXerr];
	b1 = (alist__[MAXerr] + blist[MAXerr]) * .5;
	a2 = b1;
	b2 = blist[MAXerr];
	erlast = errMAX;
	pnl_dqc25f((PnlFunc *)f, &a1, &b1, &domega, integr, &nrmom, MAXp1, &c__0, &
		area1, &error1, &nev, &resabs, &defab1, momcom, &chebmo[
		chebmo_offset]);
	*neval += nev;
	pnl_dqc25f((PnlFunc *)f, &a2, &b2, &domega, integr, &nrmom, MAXp1, &c__1, &
		area2, &error2, &nev, &resabs, &defab2, momcom, &chebmo[
		chebmo_offset]);
	*neval += nev;

/*           improve previous approximations to integral */
/*           and error and test for accuracy. */

	area12 = area1 + area2;
	erro12 = error1 + error2;
	errsum = errsum + erro12 - errMAX;
	area = area + area12 - rlist[MAXerr];
	if (defab1 == error1 || defab2 == error2) {
	    goto L25;
	}
	if ((d__1 = rlist[MAXerr] - area12, ABS(d__1)) > ABS(area12) * 1e-5 ||
		 erro12 < errMAX * .99) {
	    goto L20;
	}
	if (extrap) {
	    ++iroff2;
	}
	if (! extrap) {
	    ++iroff1;
	}
L20:
	if (*last > 10 && erro12 > errMAX) {
	    ++iroff3;
	}
L25:
	rlist[MAXerr] = area1;
	rlist[*last] = area2;
	nnlog[MAXerr] = nrmom;
	nnlog[*last] = nrmom;
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
/*           at a point of the integration range. */

/* Computing MAX */
	d__1 = ABS(a1), d__2 = ABS(b2);
	if (MAX(d__1,d__2) <= (epmach * 100. + 1.) * (ABS(a2) + uflow * 1e3)) 
		{
	    *ier = 4;
	}

/*           append the newly-created intervals to the list. */

	if (error2 > error1) {
	    goto L30;
	}
	alist__[*last] = a2;
	blist[MAXerr] = b1;
	blist[*last] = b2;
	elist[MAXerr] = error1;
	elist[*last] = error2;
	goto L40;
L30:
	alist__[MAXerr] = a2;
	alist__[*last] = a1;
	blist[*last] = b1;
	rlist[MAXerr] = area2;
	rlist[*last] = area1;
	elist[MAXerr] = error2;
	elist[*last] = error1;

/*           call subroutine dqpsrt to maintain the descending ordering */
/*           in the list of error estimates and select the subinterval */
/*           with nrMAX-th largest error estimate (to bisected next). */

L40:
	pnl_dqpsrt(limit, last, &MAXerr, &errMAX, &elist[1], &iord[1], &nrMAX);
/* ***jump out of do-loop */
	if (errsum <= errbnd) {
	    goto L170;
	}
	if (*ier != 0) {
	    goto L150;
	}
	if (*last == 2 && extall) {
	    goto L120;
	}
	if (noext) {
	    goto L140;
	}
	if (! extall) {
	    goto L50;
	}
	erlarg -= erlast;
	if ((d__1 = b1 - a1, ABS(d__1)) > small) {
	    erlarg += erro12;
	}
	if (extrap) {
	    goto L70;
	}

/*           test whether the interval to be bisected next is the */
/*           smallest interval. */

L50:
	width = (d__1 = blist[MAXerr] - alist__[MAXerr], ABS(d__1));
	if (width > small) {
	    goto L140;
	}
	if (extall) {
	    goto L60;
	}

/*           test whether we can start with the extrapolation procedure */
/*           (we do this if we integrate over the next interval with */
/*           use of a gauss-kronrod rule - see subroutine dqc25f). */

	small *= .5;
	if (width * .25 * domega > 2.) {
	    goto L140;
	}
	extall = TRUE;
	goto L130;
L60:
	extrap = TRUE;
	nrMAX = 2;
L70:
	if (ierro == 3 || erlarg <= ertest) {
	    goto L90;
	}

/*           the smallest interval has the largest error. */
/*           before bisecting decrease the sum of the errors over */
/*           the larger intervals (erlarg) and perform extrapolation. */

	jupbnd = *last;
	if (*last > *limit / 2 + 2) {
	    jupbnd = *limit + 3 - *last;
	}
	id = nrMAX;
	i__2 = jupbnd;
	for (k = id; k <= i__2; ++k) {
	    MAXerr = iord[nrMAX];
	    errMAX = elist[MAXerr];
	    if ((d__1 = blist[MAXerr] - alist__[MAXerr], ABS(d__1)) > small) {
		goto L140;
	    }
	    ++nrMAX;
/* L80: */
	}

/*           perform extrapolation. */

L90:
	++numrl2;
	rlist2[numrl2 - 1] = area;
	if (numrl2 < 3) {
	    goto L110;
	}
	pnl_dqelg(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
	++ktMIN;
	if (ktMIN > 5 && *abserr < errsum * .001) {
	    *ier = 5;
	}
	if (abseps >= *abserr) {
	    goto L100;
	}
	ktMIN = 0;
	*abserr = abseps;
	*result = reseps;
	correc = erlarg;
/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * ABS(reseps);
	ertest = MAX(d__1,d__2);
/* ***jump out of do-loop */
	if (*abserr <= ertest) {
	    goto L150;
	}

/*           prepare bisection of the smallest interval. */

L100:
	if (numrl2 == 1) {
	    noext = TRUE;
	}
	if (*ier == 5) {
	    goto L150;
	}
L110:
	MAXerr = iord[1];
	errMAX = elist[MAXerr];
	nrMAX = 1;
	extrap = FALSE;
	small *= .5;
	erlarg = errsum;
	goto L140;
L120:
	small *= .5;
	++numrl2;
	rlist2[numrl2 - 1] = area;
L130:
	ertest = errbnd;
	erlarg = errsum;
L140:
	;
    }

/*           set the final result. */
/*           --------------------- */

L150:
    if (*abserr == oflow || nres == 0) {
	goto L170;
    }
    if (*ier + ierro == 0) {
	goto L165;
    }
    if (ierro == 3) {
	*abserr += correc;
    }
    if (*ier == 0) {
	*ier = 3;
    }
    if (*result != 0. && area != 0.) {
	goto L160;
    }
    if (*abserr > errsum) {
	goto L170;
    }
    if (area == 0.) {
	goto L190;
    }
    goto L165;
L160:
    if (*abserr / ABS(*result) > errsum / ABS(area)) {
	goto L170;
    }

/*           test on divergence. */

L165:
/* Computing MAX */
    d__1 = ABS(*result), d__2 = ABS(area);
    if (ksgn == -1 && MAX(d__1,d__2) <= defabs * .01) {
	goto L190;
    }
    if (.01 > *result / area || *result / area > 100. || errsum >= ABS(area)) 
	    {
	*ier = 6;
    }
    goto L190;

/*           compute global integral sum. */

L170:
    *result = 0.;
    i__1 = *last;
    for (k = 1; k <= i__1; ++k) {
	*result += rlist[k];
/* L180: */
    }
    *abserr = errsum;
L190:
    if (*ier > 2) {
	--(*ier);
    }
L200:
    if (*integr == 2 && *omega < 0.) {
	*result = -(*result);
    }
L999:
    return 0;
} /* pnl_dqawoe */

