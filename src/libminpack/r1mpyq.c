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

/* subroutine r1mpyq
 * 
 * given an m by n matrix a, this subroutine computes a*q where
 * q is the product of 2*(n - 1) transformations
 * 
 *       gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
 * 
 * and gv(i), gw(i) are givens rotations in the (i,n) plane which
 * eliminate elements in the i-th and n-th planes, respectively.
 * q itself is not given, rather the information to recover the
 * gv, gw rotations is supplied.
 * 
 * the subroutine statement is
 * 
 *   subroutine pnl_minpack_r1mpyq(m,n,a,lda,v,w)
 * 
 * where
 * 
 *   m is a positive integer input variable set to the number
 *     of rows of a.
 * 
 *   n is a positive integer input variable set to the number
 *     of columns of a.
 * 
 *   a is an m by n array. on input a must contain the matrix
 *     to be postmultiplied by the orthogonal matrix q
 *     described above. on output a*q has replaced a.
 * 
 *   lda is a positive integer input variable not less than m
 *     which specifies the leading dimension of the array a.
 * 
 *   v is an input array of length n. v(i) must contain the
 *     information necessary to recover the givens rotation gv(i)
 *     described above.
 * 
 *   w is an input array of length n. w(i) must contain the
 *     information necessary to recover the givens rotation gw(i)
 *     described above.
 * 
 * subroutines called
 * 
 *   fortran-supplied ... dabs,dsqrt
 * 
 * argonne national laboratory. minpack project. march 1980.
 * burton s. garbow, kenneth e. hillstrom, jorge j. more
 * 
 */

void pnl_minpack_r1mpyq(int m, int n, double *a, int
	lda, const double *v, const double *w)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;
    double d__1, d__2;

    /* Local variables */
    int i__, j, nm1, nmj;
    double cos__=0.0, sin__=0.0, temp;


    /* Parameter adjustments */
    --w;
    --v;
    a_dim1 = lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */

/*     apply the first set of givens rotations to a. */

    nm1 = n - 1;
    if (nm1 < 1) {
	/* goto L50; */
        return;
    }
    i__1 = nm1;
    for (nmj = 1; nmj <= i__1; ++nmj) {
	j = n - nmj;
	if ((d__1 = v[j], ABS(d__1)) > 1.) {
	    cos__ = 1. / v[j];
	}
	if ((d__1 = v[j], ABS(d__1)) > 1.) {
/* Computing 2nd power */
	    d__2 = cos__;
	    sin__ = sqrt(1. - d__2 * d__2);
	}
	if ((d__1 = v[j], ABS(d__1)) <= 1.) {
	    sin__ = v[j];
	}
	if ((d__1 = v[j], ABS(d__1)) <= 1.) {
/* Computing 2nd power */
	    d__2 = sin__;
	    cos__ = sqrt(1. - d__2 * d__2);
	}
	i__2 = m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    temp = cos__ * a[i__ + j * a_dim1] - sin__ * a[i__ + n * a_dim1];
	    a[i__ + n * a_dim1] = sin__ * a[i__ + j * a_dim1] + cos__ * a[
		    i__ + n * a_dim1];
	    a[i__ + j * a_dim1] = temp;
/* L10: */
	}
/* L20: */
    }

/*     apply the second set of givens rotations to a. */

    i__1 = nm1;
    for (j = 1; j <= i__1; ++j) {
	if ((d__1 = w[j], ABS(d__1)) > 1.) {
	    cos__ = 1. / w[j];
	}
	if ((d__1 = w[j], ABS(d__1)) > 1.) {
/* Computing 2nd power */
	    d__2 = cos__;
	    sin__ = sqrt(1. - d__2 * d__2);
	}
	if ((d__1 = w[j], ABS(d__1)) <= 1.) {
	    sin__ = w[j];
	}
	if ((d__1 = w[j], ABS(d__1)) <= 1.) {
/* Computing 2nd power */
	    d__2 = sin__;
	    cos__ = sqrt(1. - d__2 * d__2);
	}
	i__2 = m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    temp = cos__ * a[i__ + j * a_dim1] + sin__ * a[i__ + n * a_dim1];
	    a[i__ + n * a_dim1] = -sin__ * a[i__ + j * a_dim1] + cos__ * a[
		    i__ + n * a_dim1];
	    a[i__ + j * a_dim1] = temp;
/* L30: */
	}
/* L40: */
    }
/* L50: */
    return;

/*     last card of subroutine r1mpyq. */

} /* r1mpyq_ */

