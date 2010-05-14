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

/* subroutine covar
 * 
 * given an m by n matrix a, the problem is to determine
 * the covariance matrix corresponding to a, defined as
 * 
 *                t
 *       inverse(a *a) .
 * 
 * this subroutine completes the solution of the problem
 * if it is provided with the necessary information from the
 * qr factorization, with column pivoting, of a. that is, if
 * a*p = q*r, where p is a permutation matrix, q has orthogonal
 * columns, and r is an upper triangular matrix with diagonal
 * elements of nonincreasing magnitude, then covar expects
 * the full upper triangle of r and the permutation matrix p.
 * the covariance matrix is then computed as
 * 
 *                  t     t
 *       p*inverse(r *r)*p  .
 * 
 * if a is nearly rank deficient, it may be desirable to compute
 * the covariance matrix corresponding to the linearly independent
 * columns of a. to define the numerical rank of a, covar uses
 * the tolerance tol. if l is the largest integer such that
 * 
 *       ABS(r(l,l)) .gt. tol*ABS(r(1,1)) ,
 * 
 * then covar computes the covariance matrix corresponding to
 * the first l columns of r. for k greater than l, column
 * and row ipvt(k) of the covariance matrix are set to zero.
 * 
 * the subroutine statement is
 * 
 *   pnl_minpack_covar(n,r,ldr,ipvt,tol,wa)
 * 
 * where
 * 
 *   n is a positive integer input variable set to the order of r.
 * 
 *   r is an n by n array. on input the full upper triangle must
 *     contain the full upper triangle of the matrix r. on output
 *     r contains the square symmetric covariance matrix.
 * 
 *   ldr is a positive integer input variable not less than n
 *     which specifies the leading dimension of the array r.
 * 
 *   ipvt is an integer input array of length n which defines the
 *     permutation matrix p such that a*p = q*r. column j of p
 *     is column ipvt(j) of the identity matrix.
 * 
 *   tol is a nonnegative input variable used to define the
 *     numerical rank of a in the manner described above.
 * 
 *   wa is a work array of length n.
 * 
 * subprograms called
 * 
 *   fortran-supplied ... dabs
 * 
 * argonne national laboratory. minpack project. august 1980.
 * burton s. garbow, kenneth e. hillstrom, jorge j. more
 */

void pnl_minpack_covar(int n, double *r__, int ldr, 
	const int *ipvt, double tol, double *wa)
{
    /* System generated locals */
    int r_dim1, r_offset, i__1, i__2, i__3;

    /* Local variables */
    int i__, j, k, l, ii, jj, km1;
    int sing;
    double temp, tolr;

    /* Parameter adjustments */
    --wa;
    --ipvt;
    tolr = tol * fabs(r__[0]);
    r_dim1 = ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;

    /* Function Body */

/*     form the inverse of r in the full upper triangle of r. */

    l = 0;
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	if (fabs(r__[k + k * r_dim1]) <= tolr) {
	    goto L50;
	}
	r__[k + k * r_dim1] = 1. / r__[k + k * r_dim1];
	km1 = k - 1;
	if (km1 < 1) {
	    goto L30;
	}
	i__2 = km1;
	for (j = 1; j <= i__2; ++j) {
	    temp = r__[k + k * r_dim1] * r__[j + k * r_dim1];
	    r__[j + k * r_dim1] = 0.;
	    i__3 = j;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		r__[i__ + k * r_dim1] -= temp * r__[i__ + j * r_dim1];
/* L10: */
	    }
/* L20: */
	}
L30:
	l = k;
/* L40: */
    }
L50:

/*     form the full upper triangle of the inverse of (r transpose)*r */
/*     in the full upper triangle of r. */

    if (l < 1) {
	goto L110;
    }
    i__1 = l;
    for (k = 1; k <= i__1; ++k) {
	km1 = k - 1;
	if (km1 < 1) {
	    goto L80;
	}
	i__2 = km1;
	for (j = 1; j <= i__2; ++j) {
	    temp = r__[j + k * r_dim1];
	    i__3 = j;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		r__[i__ + j * r_dim1] += temp * r__[i__ + k * r_dim1];
/* L60: */
	    }
/* L70: */
	}
L80:
	temp = r__[k + k * r_dim1];
	i__2 = k;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    r__[i__ + k * r_dim1] = temp * r__[i__ + k * r_dim1];
/* L90: */
	}
/* L100: */
    }
L110:

/*     form the full lower triangle of the covariance matrix */
/*     in the strict lower triangle of r and in wa. */

    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	jj = ipvt[j];
	sing = j > l;
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (sing) {
		r__[i__ + j * r_dim1] = 0.;
	    }
	    ii = ipvt[i__];
	    if (ii > jj) {
		r__[ii + jj * r_dim1] = r__[i__ + j * r_dim1];
	    }
	    if (ii < jj) {
		r__[jj + ii * r_dim1] = r__[i__ + j * r_dim1];
	    }
/* L120: */
	}
	wa[jj] = r__[j + j * r_dim1];
/* L130: */
    }

/*     symmetrize the covariance matrix in r. */

    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    r__[i__ + j * r_dim1] = r__[j + i__ * r_dim1];
/* L140: */
	}
	r__[j + j * r_dim1] = wa[j];
/* L150: */
    }
    /*return 0;*/

/*     last card of subroutine covar. */

} /* covar_ */

