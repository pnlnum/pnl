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

/* subroutine rwupdt
 * 
 * given an n by n upper triangular matrix r, this subroutine
 * computes the qr decomposition of the matrix formed when a row
 * is added to r. if the row is specified by the vector w, then
 * rwupdt determines an orthogonal matrix q such that when the
 * n+1 by n matrix composed of r augmented by w is premultiplied
 * by (q transpose), the resulting matrix is upper trapezoidal.
 * the matrix (q transpose) is the product of n transformations
 * 
 *       g(n)*g(n-1)* ... *g(1)
 * 
 * where g(i) is a givens rotation in the (i,n+1) plane which
 * eliminates elements in the (n+1)-st plane. rwupdt also
 * computes the product (q transpose)*c where c is the
 * (n+1)-vector (b,alpha). q itself is not accumulated, rather
 * the information to recover the g rotations is supplied.
 * 
 * the subroutine statement is
 * 
 *   subroutine pnl_minpack_rwupdt(n,r,ldr,w,b,alpha,cos,sin)
 * 
 * where
 * 
 *   n is a positive integer input variable set to the order of r.
 * 
 *   r is an n by n array. on input the upper triangular part of
 *     r must contain the matrix to be updated. on output r
 *     contains the updated triangular matrix.
 * 
 *   ldr is a positive integer input variable not less than n
 *     which specifies the leading dimension of the array r.
 * 
 *   w is an input array of length n which must contain the row
 *     vector to be added to r.
 * 
 *   b is an array of length n. on input b must contain the
 *     first n elements of the vector c. on output b contains
 *     the first n elements of the vector (q transpose)*c.
 * 
 *   alpha is a variable. on input alpha must contain the
 *     (n+1)-st element of the vector c. on output alpha contains
 *     the (n+1)-st element of the vector (q transpose)*c.
 * 
 *   cos is an output array of length n which contains the
 *     cosines of the transforming givens rotations.
 * 
 *   sin is an output array of length n which contains the
 *     sines of the transforming givens rotations.
 * 
 * subprograms called
 * 
 *   fortran-supplied ... dabs,dsqrt
 * 
 * argonne national laboratory. minpack project. march 1980.
 * burton s. garbow, dudley v. goetschel, kenneth e. hillstrom,
 * jorge j. more
 */

void pnl_minpack_rwupdt(int n, double *r__, int ldr, 
	const double *w, double *b, double *alpha, double *cos__, 
	double *sin__)
{

double p5 = .5;
double p25 = .25;

    /* System generated locals */
    int r_dim1, r_offset, i__1, i__2;
    double d__1;

    /* Local variables */
    int i__, j, jm1;
    double tan__, temp, rowj, cotan;

    /* Parameter adjustments */
    --sin__;
    --cos__;
    --b;
    --w;
    r_dim1 = ldr;
    r_offset = 1 + r_dim1 * 1;
    r__ -= r_offset;

    /* Function Body */

    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	rowj = w[j];
	jm1 = j - 1;

/*        apply the previous transformations to */
/*        r(i,j), i=1,2,...,j-1, and to w(j). */

	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    temp = cos__[i__] * r__[i__ + j * r_dim1] + sin__[i__] * rowj;
	    rowj = -sin__[i__] * r__[i__ + j * r_dim1] + cos__[i__] * rowj;
	    r__[i__ + j * r_dim1] = temp;
/* L10: */
	}
L20:

/*        determine a givens rotation which eliminates w(j). */

	cos__[j] = 1.;
	sin__[j] = 0.;
	if (rowj == 0.) {
	    goto L50;
	}
	if ((d__1 = r__[j + j * r_dim1], ABS(d__1)) >= ABS(rowj)) {
	    goto L30;
	}
	cotan = r__[j + j * r_dim1] / rowj;
/* Computing 2nd power */
	d__1 = cotan;
	sin__[j] = p5 / sqrt(p25 + p25 * (d__1 * d__1));
	cos__[j] = sin__[j] * cotan;
	goto L40;
L30:
	tan__ = rowj / r__[j + j * r_dim1];
/* Computing 2nd power */
	d__1 = tan__;
	cos__[j] = p5 / sqrt(p25 + p25 * (d__1 * d__1));
	sin__[j] = cos__[j] * tan__;
L40:

/*        apply the current transformation to r(j,j), b(j), and alpha. */

	r__[j + j * r_dim1] = cos__[j] * r__[j + j * r_dim1] + sin__[j] * 
		rowj;
	temp = cos__[j] * b[j] + sin__[j] * *alpha;
	*alpha = -sin__[j] * b[j] + cos__[j] * *alpha;
	b[j] = temp;
L50:
/* L60: */
	;
    }
    return;

/*     last card of subroutine rwupdt. */

} /* rwupdt_ */

