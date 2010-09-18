/*
 * This file contains the routine dqelg from the quadpack distribution
 * available on www.netlib.org
 *
 * This routine has been translated from Fortran using f2c and then slightly
 * modified by Jérôme Lelong for integration into Pnl. All Quadpack names have
 * been prefixed by pnl to avoid possible clash in names. Calls to d1mach have
 * been replaced by calls to dlamch coming from Lapack
 */

#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"

int pnl_dqelg(int *n, double *epstab, double *
	result, double *abserr, double *res3la, int *nres)
{
    /* System generated locals */
    int i__1;
    double d__1, d__2, d__3;

    /* Local variables */
    int i__;
    double e0, e1, e2, e3;
    int k1, k2, k3, ib, ie;
    double ss;
    int ib2;
    double res;
    int num;
    double err1, err2, err3, tol1, tol2, tol3;
    int indx;
    double e1abs, oflow, error;
    
    double delta1, delta2, delta3, epmach, epsinf;
    int newelm, limexp;

/* ***begin prologue  dqelg */
/* ***refer to  dqagie,dqagoe,dqagpe,dqagse */
/* ***revision date  830518   (yymmdd) */
/* ***keywords  epsilon algorithm, convergence acceleration, */
/*             extrapolation */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math & progr. div. - k.u.leuven */
/* ***purpose  the routine deterMINes the limit of a given sequence of */
/*            approximations, by means of the epsilon algorithm of */
/*            p.wynn. an estimate of the absolute error is also given. */
/*            the condensed epsilon table is computed. only those */
/*            elements needed for the computation of the next diagonal */
/*            are preserved. */
/* ***description */

/*           epsilon algorithm */
/*           standard fortran subroutine */
/*           double precision version */

/*           parameters */
/*              n      - int */
/*                       epstab(n) contains the new element in the */
/*                       first column of the epsilon table. */

/*              epstab - double precision */
/*                       vector of dimension 52 containing the elements */
/*                       of the two lower diagonals of the triangular */
/*                       epsilon table. the elements are numbered */
/*                       starting at the right-hand corner of the */
/*                       triangle. */

/*              result - double precision */
/*                       resulting approximation to the integral */

/*              abserr - double precision */
/*                       estimate of the absolute error computed from */
/*                       result and the 3 previous results */

/*              res3la - double precision */
/*                       vector of dimension 3 containing the last 3 */
/*                       results */

/*              nres   - int */
/*                       number of calls to the routine */
/*                       (should be zero at first call) */

/* ***end prologue  dqelg */


/*           list of major variables */
/*           ----------------------- */

/*           e0     - the 4 elements on which the computation of a new */
/*           e1       element in the epsilon table is based */
/*           e2 */
/*           e3                 e0 */
/*                        e3    e1    new */
/*                              e2 */
/*           newelm - number of elements to be computed in the new */
/*                    diagonal */
/*           error  - error = ABS(e1-e0)+ABS(e2-e1)+ABS(new-e2) */
/*           result - the element in the new diagonal with least value */
/*                    of error */

/*           machine dependent constants */
/*           --------------------------- */

/*           epmach is the largest relative spacing. */
/*           oflow is the largest positive magnitude. */
/*           limexp is the MAXimum number of elements the epsilon */
/*           table can contain. if this number is reached, the upper */
/*           diagonal of the epsilon table is deleted. */

/* ***first executable statement  dqelg */
    /* Parameter adjustments */
    --res3la;
    --epstab;

    /* Function Body */
    epmach = pnl_d1mach(4);
    oflow = pnl_d1mach(2);
    ++(*nres);
    *abserr = oflow;
    *result = epstab[*n];
    if (*n < 3) {
	goto L100;
    }
    limexp = 50;
    epstab[*n + 2] = epstab[*n];
    newelm = (*n - 1) / 2;
    epstab[*n] = oflow;
    num = *n;
    k1 = *n;
    i__1 = newelm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k2 = k1 - 1;
	k3 = k1 - 2;
	res = epstab[k1 + 2];
	e0 = epstab[k3];
	e1 = epstab[k2];
	e2 = res;
	e1abs = ABS(e1);
	delta2 = e2 - e1;
	err2 = ABS(delta2);
/* Computing MAX */
	d__1 = ABS(e2);
	tol2 = MAX(d__1,e1abs) * epmach;
	delta3 = e1 - e0;
	err3 = ABS(delta3);
/* Computing MAX */
	d__1 = e1abs, d__2 = ABS(e0);
	tol3 = MAX(d__1,d__2) * epmach;
	if (err2 > tol2 || err3 > tol3) {
	    goto L10;
	}

/*           if e0, e1 and e2 are equal to within machine */
/*           accuracy, convergence is assumed. */
/*           result = e2 */
/*           abserr = ABS(e1-e0)+ABS(e2-e1) */

	*result = res;
	*abserr = err2 + err3;
/* ***jump out of do-loop */
	goto L100;
L10:
	e3 = epstab[k1];
	epstab[k1] = e1;
	delta1 = e1 - e3;
	err1 = ABS(delta1);
/* Computing MAX */
	d__1 = e1abs, d__2 = ABS(e3);
	tol1 = MAX(d__1,d__2) * epmach;

/*           if two elements are very close to each other, omit */
/*           a part of the table by adjusting the value of n */

	if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3) {
	    goto L20;
	}
	ss = 1. / delta1 + 1. / delta2 - 1. / delta3;
	epsinf = (d__1 = ss * e1, ABS(d__1));

/*           test to detect irregular behaviour in the table, and */
/*           eventually omit a part of the table adjusting the value */
/*           of n. */

	if (epsinf > 1e-4) {
	    goto L30;
	}
L20:
	*n = i__ + i__ - 1;
/* ***jump out of do-loop */
	goto L50;

/*           compute a new element and eventually adjust */
/*           the value of result. */

L30:
	res = e1 + 1. / ss;
	epstab[k1] = res;
	k1 += -2;
	error = err2 + (d__1 = res - e2, ABS(d__1)) + err3;
	if (error > *abserr) {
	    goto L40;
	}
	*abserr = error;
	*result = res;
L40:
	;
    }

/*           shift the table. */

L50:
    if (*n == limexp) {
	*n = (limexp / 2 << 1) - 1;
    }
    ib = 1;
    if (num / 2 << 1 == num) {
	ib = 2;
    }
    ie = newelm + 1;
    i__1 = ie;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ib2 = ib + 2;
	epstab[ib] = epstab[ib2];
	ib = ib2;
/* L60: */
    }
    if (num == *n) {
	goto L80;
    }
    indx = num - *n + 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	epstab[i__] = epstab[indx];
	++indx;
/* L70: */
    }
L80:
    if (*nres >= 4) {
	goto L90;
    }
    res3la[*nres] = *result;
    *abserr = oflow;
    goto L100;

/*           compute error estimate */

L90:
    *abserr = (d__1 = *result - res3la[3], ABS(d__1)) + (d__2 = *result - 
	    res3la[2], ABS(d__2)) + (d__3 = *result - res3la[1], ABS(d__3));
    res3la[1] = res3la[2];
    res3la[2] = res3la[3];
    res3la[3] = *result;
L100:
/* Computing MAX */
    d__1 = *abserr, d__2 = epmach * 5. * ABS(*result);
    *abserr = MAX(d__1,d__2);
    return 0;
} /* pnl_dqelg */

