/*
 * This file contains the routine dqpsrt from the quadpack distribution
 * available on www.netlib.org
 *
 * This routine has been translated from Fortran using f2c and then slightly
 * modified by Jérôme Lelong for integration into Pnl. All Quadpack names have
 * been prefixed by pnl to avoid possible clash in names. Calls to d1mach have
 * been replaced by calls to dlamch coming from Lapack
 */

#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"

int pnl_dqpsrt(int *limit, int *last, int *MAXerr, 
	double *erMAX, double *elist, int *iord, int *nrMAX)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    int i__, j, k, ido, ibeg, jbnd, isucc, jupbn;
    double errMIN, errMAX;

/* ***begin prologue  dqpsrt */
/* ***refer to  dqage,dqagie,dqagpe,dqawse */
/* ***routines called  (none) */
/* ***revision date  810101   (yymmdd) */
/* ***keywords  sequential sorting */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  this routine maintains the descending ordering in the */
/*            list of the local error estimated resulting from the */
/*            interval subdivision process. at each call two error */
/*            estimates are inserted using the sequential search */
/*            method, top-down for the largest error estimate and */
/*            bottom-up for the smallest error estimate. */
/* ***description */

/*           ordering routine */
/*           standard fortran subroutine */
/*           double precision version */

/*           parameters (meaning at output) */
/*              limit  - int */
/*                       MAXimum number of error estimates the list */
/*                       can contain */

/*              last   - int */
/*                       number of error estimates currently in the list */

/*              MAXerr - int */
/*                       MAXerr points to the nrMAX-th largest error */
/*                       estimate currently in the list */

/*              erMAX  - double precision */
/*                       nrMAX-th largest error estimate */
/*                       erMAX = elist(MAXerr) */

/*              elist  - double precision */
/*                       vector of dimension last containing */
/*                       the error estimates */

/*              iord   - int */
/*                       vector of dimension last, the first k elements */
/*                       of which contain pointers to the error */
/*                       estimates, such that */
/*                       elist(iord(1)),...,  elist(iord(k)) */
/*                       form a decreasing sequence, with */
/*                       k = last if last.le.(limit/2+2), and */
/*                       k = limit+1-last otherwise */

/*              nrMAX  - int */
/*                       MAXerr = iord(nrMAX) */

/* ***end prologue  dqpsrt */


/*           check whether the list contains more than */
/*           two error estimates. */

/* ***first executable statement  dqpsrt */
    /* Parameter adjustments */
    --iord;
    --elist;

    /* Function Body */
    if (*last > 2) {
	goto L10;
    }
    iord[1] = 1;
    iord[2] = 2;
    goto L90;

/*           this part of the routine is only executed if, due to a */
/*           difficult integrand, subdivision increased the error */
/*           estimate. in the normal case the insert procedure should */
/*           start after the nrMAX-th largest error estimate. */

L10:
    errMAX = elist[*MAXerr];
    if (*nrMAX == 1) {
	goto L30;
    }
    ido = *nrMAX - 1;
    i__1 = ido;
    for (i__ = 1; i__ <= i__1; ++i__) {
	isucc = iord[*nrMAX - 1];
/* ***jump out of do-loop */
	if (errMAX <= elist[isucc]) {
	    goto L30;
	}
	iord[*nrMAX] = isucc;
	--(*nrMAX);
/* L20: */
    }

/*           compute the number of elements in the list to be maintained */
/*           in descending order. this number depends on the number of */
/*           subdivisions still allowed. */

L30:
    jupbn = *last;
    if (*last > *limit / 2 + 2) {
	jupbn = *limit + 3 - *last;
    }
    errMIN = elist[*last];

/*           insert errMAX by traversing the list top-down, */
/*           starting comparison from the element elist(iord(nrMAX+1)). */

    jbnd = jupbn - 1;
    ibeg = *nrMAX + 1;
    if (ibeg > jbnd) {
	goto L50;
    }
    i__1 = jbnd;
    for (i__ = ibeg; i__ <= i__1; ++i__) {
	isucc = iord[i__];
/* ***jump out of do-loop */
	if (errMAX >= elist[isucc]) {
	    goto L60;
	}
	iord[i__ - 1] = isucc;
/* L40: */
    }
L50:
    iord[jbnd] = *MAXerr;
    iord[jupbn] = *last;
    goto L90;

/*           insert errMIN by traversing the list bottom-up. */

L60:
    iord[i__ - 1] = *MAXerr;
    k = jbnd;
    i__1 = jbnd;
    for (j = i__; j <= i__1; ++j) {
	isucc = iord[k];
/* ***jump out of do-loop */
	if (errMIN < elist[isucc]) {
	    goto L80;
	}
	iord[k + 1] = isucc;
	--k;
/* L70: */
    }
    iord[i__] = *last;
    goto L90;
L80:
    iord[k + 1] = *last;

/*           set MAXerr and erMAX. */

L90:
    *MAXerr = iord[*nrMAX];
    *erMAX = elist[*MAXerr];
    return 0;
} /* pnl_dqpsrt */

