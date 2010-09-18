/*
 * This file contains the routine dqawfe from the quadpack distribution
 * available on www.netlib.org
 *
 * This routine has been translated from Fortran using f2c and then slightly
 * modified by Jérôme Lelong for integration into Pnl. All Quadpack names have
 * been prefixed by pnl to avoid possible clash in names. Calls to d1mach have
 * been replaced by calls to dlamch coming from Lapack
 */

#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"

extern int pnl_dqelg(int *, double *, double *, 
                  double *, double *, int *);
extern int pnl_dqagie(PnlFunc *, double *, int *, 
                   double *, double *, int *, double *, double *,
                   int *, int *, double *, double *, double *, 
                   double *, int *, int *);
extern int pnl_dqawoe(PnlFunc *, double *, double *, 
                   double *, int *, double *, double *, int *, 
                   int *, int *, double *, double *, int *, 
                   int *, int *, double *, double *, double *, 
                   double *, int *, int *, int *, double *);


int pnl_dqawfe(PnlFunc * f, double *a, double *omega, 
	int *integr, double *epsabs, int *limlst, int *limit, 
	int *MAXp1, double *result, double *abserr, int *
	neval, int *ier, double *rslst, double *erlst, int *
	ierlst, int *lst, double *alist__, double *blist, 
	double *rlist, double *elist, int *iord, int *nnlog, 
	double *chebmo)
{
    /* Initialized data */

    double p = .9;

    /* System generated locals */
    int chebmo_dim1, chebmo_offset, i__1;
    double d__1, d__2, c_b4;
    int c__1;

    /* Local variables */
    int l;
    double c1, c2, p1, dl, ep;
    int ll;
    double dla, drl, eps;
    int nev;
    double fact, epsa;
    int last, nres;
    double psum[52];
    double cycle;
    int ktMIN;
    double uflow;
    
    double res3la[3];
    int numrl2;
    double abseps, correc;
    int momcom;
    double reseps, errsum;

/* ***begin prologue  dqawfe */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a3a1 */
/* ***keywords  automatic integrator, special-purpose, */
/*             fourier integrals, */
/*             integration between zeros with dqawoe, */
/*             convergence acceleration with dqelg */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           dedoncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  the routine calculates an approximation result to a */
/*            given fourier integal */
/*            i = integral of f(x)*w(x) over (a,infinity) */
/*            where w(x)=cos(omega*x) or w(x)=sin(omega*x), */
/*            hopefully satisfying following claim for accuracy */
/*            ABS(i-result).le.epsabs. */
/* ***description */

/*        computation of fourier integrals */
/*        standard fortran subroutine */
/*        double precision version */

/*        parameters */
/*         on entry */
/*            f      - double precision */
/*                     function subprogram defining the integrand */
/*                     function f(x). the actual name for f needs to */
/*                     be declared e x t e r n a l in the driver program. */

/*            a      - double precision */
/*                     lower limit of integration */

/*            omega  - double precision */
/*                     parameter in the weight function */

/*            integr - int */
/*                     indicates which weight function is used */
/*                     integr = 1      w(x) = cos(omega*x) */
/*                     integr = 2      w(x) = sin(omega*x) */
/*                     if integr.ne.1.and.integr.ne.2, the routine will */
/*                     end with ier = 6. */

/*            epsabs - double precision */
/*                     absolute accuracy requested, epsabs.gt.0 */
/*                     if epsabs.le.0, the routine will end with ier = 6. */

/*            limlst - int */
/*                     limlst gives an upper bound on the number of */
/*                     cycles, limlst.ge.1. */
/*                     if limlst.lt.3, the routine will end with ier = 6. */

/*            limit  - int */
/*                     gives an upper bound on the number of subintervals */
/*                     allowed in the partition of each cycle, limit.ge.1 */
/*                     each cycle, limit.ge.1. */

/*            MAXp1  - int */
/*                     gives an upper bound on the number of */
/*                     chebyshev moments which can be stored, i.e. */
/*                     for the intervals of lengths ABS(b-a)*2**(-l), */
/*                     l=0,1, ..., MAXp1-2, MAXp1.ge.1 */

/*         on return */
/*            result - double precision */
/*                     approximation to the integral x */

/*            abserr - double precision */
/*                     estimate of the modulus of the absolute error, */
/*                     which should equal or exceed ABS(i-result) */

/*            neval  - int */
/*                     number of integrand evaluations */

/*            ier    - ier = 0 normal and reliable terMINation of */
/*                             the routine. it is assumed that the */
/*                             requested accuracy has been achieved. */
/*                     ier.gt.0 abnormal terMINation of the routine. the */
/*                             estimates for integral and error are less */
/*                             reliable. it is assumed that the requested */
/*                             accuracy has not been achieved. */
/*            error messages */
/*                    if omega.ne.0 */
/*                     ier = 1 MAXimum number of  cycles  allowed */
/*                             has been achieved., i.e. of subintervals */
/*                             (a+(k-1)c,a+kc) where */
/*                             c = (2*int(ABS(omega))+1)*pi/ABS(omega), */
/*                             for k = 1, 2, ..., lst. */
/*                             one can allow more cycles by increasing */
/*                             the value of limlst (and taking the */
/*                             according dimension adjustments into */
/*                             account). */
/*                             exaMINe the array iwork which contains */
/*                             the error flags on the cycles, in order to */
/*                             look for eventual local integration */
/*                             difficulties. if the position of a local */
/*                             difficulty can be deterMINed (e.g. */
/*                             singularity, discontinuity within the */
/*                             interval) one will probably gain from */
/*                             splitting up the interval at this point */
/*                             and calling appropriate integrators on */
/*                             the subranges. */
/*                         = 4 the extrapolation table constructed for */
/*                             convergence acceleration of the series */
/*                             formed by the integral contributions over */
/*                             the cycles, does not converge to within */
/*                             the requested accuracy. as in the case of */
/*                             ier = 1, it is advised to exaMINe the */
/*                             array iwork which contains the error */
/*                             flags on the cycles. */
/*                         = 6 the input is invalid because */
/*                             (integr.ne.1 and integr.ne.2) or */
/*                              epsabs.le.0 or limlst.lt.3. */
/*                              result, abserr, neval, lst are set */
/*                              to zero. */
/*                         = 7 bad integrand behaviour occurs within one */
/*                             or more of the cycles. location and type */
/*                             of the difficulty involved can be */
/*                             deterMINed from the vector ierlst. here */
/*                             lst is the number of cycles actually */
/*                             needed (see below). */
/*                             ierlst(k) = 1 the MAXimum number of */
/*                                           subdivisions (= limit) has */
/*                                           been achieved on the k th */
/*                                           cycle. */
/*                                       = 2 occurrence of roundoff error */
/*                                           is detected and prevents the */
/*                                           tolerance imposed on the */
/*                                           k th cycle, from being */
/*                                           achieved. */
/*                                       = 3 extremely bad integrand */
/*                                           behaviour occurs at some */
/*                                           points of the k th cycle. */
/*                                       = 4 the integration procedure */
/*                                           over the k th cycle does */
/*                                           not converge (to within the */
/*                                           required accuracy) due to */
/*                                           roundoff in the */
/*                                           extrapolation procedure */
/*                                           invoked on this cycle. it */
/*                                           is assumed that the result */
/*                                           on this interval is the */
/*                                           best which can be obtained. */
/*                                       = 5 the integral over the k th */
/*                                           cycle is probably divergent */
/*                                           or slowly convergent. it */
/*                                           must be noted that */
/*                                           divergence can occur with */
/*                                           any other value of */
/*                                           ierlst(k). */
/*                    if omega = 0 and integr = 1, */
/*                    the integral is calculated by means of dqagie */
/*                    and ier = ierlst(1) (with meaning as described */
/*                    for ierlst(k), k = 1). */

/*            rslst  - double precision */
/*                     vector of dimension at least limlst */
/*                     rslst(k) contains the integral contribution */
/*                     over the interval (a+(k-1)c,a+kc) where */
/*                     c = (2*int(ABS(omega))+1)*pi/ABS(omega), */
/*                     k = 1, 2, ..., lst. */
/*                     note that, if omega = 0, rslst(1) contains */
/*                     the value of the integral over (a,infinity). */

/*            erlst  - double precision */
/*                     vector of dimension at least limlst */
/*                     erlst(k) contains the error estimate corresponding */
/*                     with rslst(k). */

/*            ierlst - int */
/*                     vector of dimension at least limlst */
/*                     ierlst(k) contains the error flag corresponding */
/*                     with rslst(k). for the meaning of the local error */
/*                     flags see description of output parameter ier. */

/*            lst    - int */
/*                     number of subintervals needed for the integration */
/*                     if omega = 0 then lst is set to 1. */

/*            alist, blist, rlist, elist - double precision */
/*                     vector of dimension at least limit, */

/*            iord, nnlog - int */
/*                     vector of dimension at least limit, providing */
/*                     space for the quantities needed in the subdivision */
/*                     process of each cycle */

/*            chebmo - double precision */
/*                     array of dimension at least (MAXp1,25), providing */
/*                     space for the chebyshev moments needed within the */
/*                     cycles */

/* ***references  (none) */
/* ***routines called  dqagie,dqawoe,dqelg */
/* ***end prologue  dqawfe */





/*            the dimension of  psum  is deterMINed by the value of */
/*            limexp in subroutine dqelg (psum must be of dimension */
/*            (limexp+2) at least). */

/*           list of major variables */
/*           ----------------------- */

/*           c1, c2    - end points of subinterval (of length cycle) */
/*           cycle     - (2*int(ABS(omega))+1)*pi/ABS(omega) */
/*           psum      - vector of dimension at least (limexp+2) */
/*                       (see routine dqelg) */
/*                       psum contains the part of the epsilon table */
/*                       which is still needed for further computations. */
/*                       each element of psum is a partial sum of the */
/*                       series which should sum to the value of the */
/*                       integral. */
/*           errsum    - sum of error estimates over the subintervals, */
/*                       calculated cumulatively */
/*           epsa      - absolute tolerance requested over current */
/*                       subinterval */
/*           chebmo    - array containing the modified chebyshev */
/*                       moments (see also routine dqc25f) */

    /* Parameter adjustments */
    --ierlst;
    --erlst;
    --rslst;
    --nnlog;
    --iord;
    --elist;
    --rlist;
    --blist;
    --alist__;
    chebmo_dim1 = *MAXp1;
    chebmo_offset = 1 + chebmo_dim1;
    chebmo -= chebmo_offset;

    /* Function Body */

/*           test on validity of parameters */
/*           ------------------------------ */

/* ***first executable statement  dqawfe */
    *result = 0.;
    *abserr = 0.;
    *neval = 0;
    *lst = 0;
    *ier = 0;
    if ((*integr != 1 && *integr != 2) || *epsabs <= 0. || *limlst < 3) {
	*ier = 6;
    }
    if (*ier == 6) {
	goto L999;
    }
    if (*omega != 0.) {
	goto L10;
    }

/*           integration by dqagie if omega is zero */
/*           -------------------------------------- */

    if (*integr == 1) {
	pnl_dqagie((PnlFunc *)f, &c_b4, &c__1, epsabs, &c_b4, limit, result, abserr, 
		neval, ier, &alist__[1], &blist[1], &rlist[1], &elist[1], &
		iord[1], &last);
    }
    rslst[1] = *result;
    erlst[1] = *abserr;
    ierlst[1] = *ier;
    *lst = 1;
    goto L999;

/*           initializations */
/*           --------------- */

L10:
    l = (int) ABS(*omega);
    dl = (double) ((l << 1) + 1);
    cycle = dl * M_PI / ABS(*omega);
    *ier = 0;
    ktMIN = 0;
    *neval = 0;
    numrl2 = 0;
    nres = 0;
    c1 = *a;
    c2 = cycle + *a;
    p1 = 1. - p;
    uflow = pnl_d1mach(1);
    eps = *epsabs;
    if (*epsabs > uflow / p1) {
	eps = *epsabs * p1;
    }
    ep = eps;
    fact = 1.;
    correc = 0.;
    *abserr = 0.;
    errsum = 0.;
    c__1 = 1;
    c_b4 = 0;

/*           main do-loop */
/*           ------------ */

    i__1 = *limlst;
    for (*lst = 1; *lst <= i__1; ++(*lst)) {

/*           integrate over current subinterval. */

	dla = (double) (*lst);
	epsa = eps * fact;
	pnl_dqawoe((PnlFunc *)f, &c1, &c2, omega, integr, &epsa, &c_b4, limit, lst, 
		MAXp1, &rslst[*lst], &erlst[*lst], &nev, &ierlst[*lst], &last,
		 &alist__[1], &blist[1], &rlist[1], &elist[1], &iord[1], &
		nnlog[1], &momcom, &chebmo[chebmo_offset]);
	*neval += nev;
	fact *= p;
	errsum += erlst[*lst];
	drl = (d__1 = rslst[*lst], ABS(d__1)) * 50.;

/*           test on accuracy with partial sum */

	if (errsum + drl <= *epsabs && *lst >= 6) {
	    goto L80;
	}
/* Computing MAX */
	d__1 = correc, d__2 = erlst[*lst];
	correc = MAX(d__1,d__2);
	if (ierlst[*lst] != 0) {
/* Computing MAX */
	    d__1 = ep, d__2 = correc * p1;
	    eps = MAX(d__1,d__2);
	}
	if (ierlst[*lst] != 0) {
	    *ier = 7;
	}
	if (*ier == 7 && errsum + drl <= correc * 10. && *lst > 5) {
	    goto L80;
	}
	++numrl2;
	if (*lst > 1) {
	    goto L20;
	}
	psum[0] = rslst[1];
	goto L40;
L20:
	psum[numrl2 - 1] = psum[ll - 1] + rslst[*lst];
	if (*lst == 2) {
	    goto L40;
	}

/*           test on MAXimum number of subintervals */

	if (*lst == *limlst) {
	    *ier = 1;
	}

/*           perform new extrapolation */

	pnl_dqelg(&numrl2, psum, &reseps, &abseps, res3la, &nres);

/*           test whether extrapolated result is influenced by roundoff */

	++ktMIN;
	if (ktMIN >= 15 && *abserr <= (errsum + drl) * .001) {
	    *ier = 4;
	}
	if (abseps > *abserr && *lst != 3) {
	    goto L30;
	}
	*abserr = abseps;
	*result = reseps;
	ktMIN = 0;

/*           if ier is not 0, check whether direct result (partial sum) */
/*           or extrapolated result yields the best integral */
/*           approximation */

	if (*abserr + correc * 10. <= *epsabs || (*abserr <= *epsabs && correc * 10. >= *epsabs)) {
	    goto L60;
	}
L30:
	if (*ier != 0 && *ier != 7) {
	    goto L60;
	}
L40:
	ll = numrl2;
	c1 = c2;
	c2 += cycle;
/* L50: */
    }

/*         set final result and error estimate */
/*         ----------------------------------- */

L60:
    *abserr += correc * 10.;
    if (*ier == 0) {
	goto L999;
    }
    if (*result != 0. && psum[numrl2 - 1] != 0.) {
	goto L70;
    }
    if (*abserr > errsum) {
	goto L80;
    }
    if (psum[numrl2 - 1] == 0.) {
	goto L999;
    }
L70:
    if (*abserr / ABS(*result) > (errsum + drl) / (d__1 = psum[numrl2 - 1], 
	    ABS(d__1))) {
	goto L80;
    }
    if (*ier >= 1 && *ier != 7) {
	*abserr += drl;
    }
    goto L999;
L80:
    *result = psum[numrl2 - 1];
    *abserr = errsum + drl;
L999:
    return 0;
} /* pnl_dqawfe */

