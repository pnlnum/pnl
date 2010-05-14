/* Minpack Copyright Notice (1999) University of Chicago.  All rights reserved
 * 
 * Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the
 * following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above
 * copyright notice, this list of conditions and the following
 * disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above
 * copyright notice, this list of conditions and the following
 * disclaimer in the documentation and/or other materials
 * provided with the distribution.
 * 
 * 3. The end-user documentation included with the
 * redistribution, if any, must include the following
 * acknowledgment:
 * 
 *    "This product includes software developed by the
 *    University of Chicago, as Operator of Argonne National
 *    Laboratory.
 * 
 * Alternately, this acknowledgment may appear in the software
 * itself, if and wherever such third-party acknowledgments
 * normally appear.
 * 
 * 4. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
 * WITHOUT WARRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
 * UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
 * THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
 * OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
 * OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
 * USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
 * THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
 * DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
 * UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL
 * BE CORRECTED.
 * 
 * 5. LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT
 * HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
 * ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
 * INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
 * ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
 * PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
 * SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
 * (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE,
 * EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
 * POSSIBILITY OF SUCH LOSS OR DAMAGES.
 */

#include <math.h>
#include "cminpack.h"



/* Table of constant values */

 void pnl_minpack_dogleg(int n, const double *r__, int lr, 
	const double *diag, const double *qtb, double delta, double *x, 
	double *wa1, double *wa2)
{
    /* System generated locals */
    int i__1, i__2;
    double d__1, d__2, d__3, d__4;

    /* Local variables */
    int i__, j, k, l, jj, jp1;
    double sum, temp, alpha, bnorm;
    double gnorm, qnorm, epsmch;
    double sgnorm;

/*     ********** */

/*     subroutine dogleg */

/*     given an m by n matrix a, an n by n nonsingular diagonal */
/*     matrix d, an m-vector b, and a positive number delta, the */
/*     problem is to determine the convex combination x of the */
/*     gauss-newton and scaled gradient directions that minimizes */
/*     (a*x - b) in the least squares sense, subject to the */
/*     restriction that the euclidean norm of d*x be at most delta. */

/*     this subroutine completes the solution of the problem */
/*     if it is provided with the necessary information from the */
/*     qr factorization of a. that is, if a = q*r, where q has */
/*     orthogonal columns and r is an upper triangular matrix, */
/*     then dogleg expects the full upper triangle of r and */
/*     the first n components of (q transpose)*b. */

/*     the subroutine statement is */

/*       subroutine pnl_minpack_dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2) */

/*     where */

/*       n is a positive integer input variable set to the order of r. */

/*       r is an input array of length lr which must contain the upper */
/*         triangular matrix r stored by rows. */

/*       lr is a positive integer input variable not less than */
/*         (n*(n+1))/2. */

/*       diag is an input array of length n which must contain the */
/*         diagonal elements of the matrix d. */

/*       qtb is an input array of length n which must contain the first */
/*         n elements of the vector (q transpose)*b. */

/*       delta is a positive input variable which specifies an upper */
/*         bound on the euclidean norm of d*x. */

/*       x is an output array of length n which contains the desired */
/*         convex combination of the gauss-newton direction and the */
/*         scaled gradient direction. */

/*       wa1 and wa2 are work arrays of length n. */

/*     subprograms called */

/*       minpack-supplied ... dpmpar,enorm */

/*       fortran-supplied ... dabs,dmax1,dmin1,dsqrt */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

/*     ********** */
    /* Parameter adjustments */
    --wa2;
    --wa1;
    --x;
    --qtb;
    --diag;
    --r__;

    /* Function Body */

/*     epsmch is the machine precision. */

    epsmch = pnl_minpack_dpmpar(1);

/*     first, calculate the gauss-newton direction. */

    jj = n * (n + 1) / 2 + 1;
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	j = n - k + 1;
	jp1 = j + 1;
	jj -= k;
	l = jj + 1;
	sum = 0.;
	if (n < jp1) {
	    goto L20;
	}
	i__2 = n;
	for (i__ = jp1; i__ <= i__2; ++i__) {
	    sum += r__[l] * x[i__];
	    ++l;
/* L10: */
	}
L20:
	temp = r__[jj];
	if (temp != 0.) {
	    goto L40;
	}
	l = j;
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = temp, d__3 = fabs(r__[l]);
	    temp = MAX(d__2,d__3);
	    l = l + n - i__;
/* L30: */
	}
	temp = epsmch * temp;
	if (temp == 0.) {
	    temp = epsmch;
	}
L40:
	x[j] = (qtb[j] - sum) / temp;
/* L50: */
    }

/*     test whether the gauss-newton direction is acceptable. */

    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	wa1[j] = 0.;
	wa2[j] = diag[j] * x[j];
/* L60: */
    }
    qnorm = pnl_minpack_enorm(n, &wa2[1]);
    if (qnorm <= delta) {
	/* goto L140; */
        return;
    }

/*     the gauss-newton direction is not acceptable. */
/*     next, calculate the scaled gradient direction. */

    l = 1;
    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	temp = qtb[j];
	i__2 = n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    wa1[i__] += r__[l] * temp;
	    ++l;
/* L70: */
	}
	wa1[j] /= diag[j];
/* L80: */
    }

/*     calculate the norm of the scaled gradient and test for */
/*     the special case in which the scaled gradient is zero. */

    gnorm = pnl_minpack_enorm(n, &wa1[1]);
    sgnorm = 0.;
    alpha = delta / qnorm;
    if (gnorm == 0.) {
	goto L120;
    }

/*     calculate the point along the scaled gradient */
/*     at which the quadratic is minimized. */

    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	wa1[j] = wa1[j] / gnorm / diag[j];
/* L90: */
    }
    l = 1;
    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	sum = 0.;
	i__2 = n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    sum += r__[l] * wa1[i__];
	    ++l;
/* L100: */
	}
	wa2[j] = sum;
/* L110: */
    }
    temp = pnl_minpack_enorm(n, &wa2[1]);
    sgnorm = gnorm / temp / temp;

/*     test whether the scaled gradient direction is acceptable. */

    alpha = 0.;
    if (sgnorm >= delta) {
	goto L120;
    }

/*     the scaled gradient direction is not acceptable. */
/*     finally, calculate the point along the dogleg */
/*     at which the quadratic is minimized. */

    bnorm = pnl_minpack_enorm(n, &qtb[1]);
    temp = bnorm / gnorm * (bnorm / qnorm) * (sgnorm / delta);
/* Computing 2nd power */
    d__1 = sgnorm / delta;
/* Computing 2nd power */
    d__2 = temp - delta / qnorm;
/* Computing 2nd power */
    d__3 = delta / qnorm;
/* Computing 2nd power */
    d__4 = sgnorm / delta;
    temp = temp - delta / qnorm * (d__1 * d__1) + sqrt(d__2 * d__2 + (1. - 
	    d__3 * d__3) * (1. - d__4 * d__4));
/* Computing 2nd power */
    d__1 = sgnorm / delta;
    alpha = delta / qnorm * (1. - d__1 * d__1) / temp;
L120:

/*     form appropriate convex combination of the gauss-newton */
/*     direction and the scaled gradient direction. */

    temp = (1. - alpha) * MIN(sgnorm,delta);
    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	x[j] = temp * wa1[j] + alpha * x[j];
/* L130: */
    }
/* L140: */
    return;

/*     last card of subroutine dogleg. */

} /* dogleg_ */

