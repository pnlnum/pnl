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

/* subroutine r1updt
 * 
 * given an m by n lower trapezoidal matrix s, an m-vector u,
 * and an n-vector v, the problem is to determine an
 * orthogonal matrix q such that
 * 
 *               t
 *       (s + u*v )*q
 * 
 * is again lower trapezoidal.
 * 
 * this subroutine determines q as the product of 2*(n - 1)
 * transformations
 * 
 *       gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
 * 
 * where gv(i), gw(i) are givens rotations in the (i,n) plane
 * which eliminate elements in the i-th and n-th planes,
 * respectively. q itself is not accumulated, rather the
 * information to recover the gv, gw rotations is returned.
 * 
 * the subroutine statement is
 * 
 *   subroutine pnl_minpack_r1updt(m,n,s,ls,u,v,w,sing)
 * 
 * where
 * 
 *   m is a positive integer input variable set to the number
 *     of rows of s.
 * 
 *   n is a positive integer input variable set to the number
 *     of columns of s. n must not exceed m.
 * 
 *   s is an array of length ls. on input s must contain the lower
 *     trapezoidal matrix s stored by columns. on output s contains
 *     the lower trapezoidal matrix produced as described above.
 * 
 *   ls is a positive integer input variable not less than
 *     (n*(2*m-n+1))/2.
 * 
 *   u is an input array of length m which must contain the
 *     vector u.
 * 
 *   v is an array of length n. on input v must contain the vector
 *     v. on output v(i) contains the information necessary to
 *     recover the givens rotation gv(i) described above.
 * 
 *   w is an output array of length m. w(i) contains information
 *     necessary to recover the givens rotation gw(i) described
 *     above.
 * 
 *   sing is a logical output variable. sing is set true if any
 *     of the diagonal elements of the output s are zero. otherwise
 *     sing is set false.
 * 
 * subprograms called
 * 
 *   minpack-supplied ... dpmpar
 * 
 *   fortran-supplied ... dabs,dsqrt
 * 
 * argonne national laboratory. minpack project. march 1980.
 * burton s. garbow, kenneth e. hillstrom, jorge j. more,
 * john l. nazareth
 */


void pnl_minpack_r1updt(int m, int n, double *s, int
	ls, const double *u, double *v, double *w, int *sing)
{

    double p5 = .5;
    double p25 = .25;

    /* System generated locals */
    int i__1, i__2;
    double d__1, d__2;

    /* Local variables */
    int i__, j, l, jj, nm1;
    double tan__;
    int nmj;
    double cos__, sin__, tau, temp, giant, cotan;


/*     ********** */
    /* Parameter adjustments */
    --w;
    --u;
    --v;
    --s;

    /* Function Body */

/*     giant is the largest magnitude. */

    giant = pnl_minpack_dpmpar(3);

/*     initialize the diagonal element pointer. */

    jj = n * ((m << 1) - n + 1) / 2 - (m - n);

/*     move the nontrivial part of the last column of s into w. */

    l = jj;
    i__1 = m;
    for (i__ = n; i__ <= i__1; ++i__) {
	w[i__] = s[l];
	++l;
/* L10: */
    }

/*     rotate the vector v into a multiple of the n-th unit vector */
/*     in such a way that a spike is introduced into w. */

    nm1 = n - 1;
    if (nm1 < 1) {
	goto L70;
    }
    i__1 = nm1;
    for (nmj = 1; nmj <= i__1; ++nmj) {
	j = n - nmj;
	jj -= m - j + 1;
	w[j] = 0.;
	if (v[j] == 0.) {
	    goto L50;
	}

/*        determine a givens rotation which eliminates the */
/*        j-th element of v. */

	if ((d__1 = v[n], ABS(d__1)) >= (d__2 = v[j], ABS(d__2))) {
	    goto L20;
	}
	cotan = v[n] / v[j];
/* Computing 2nd power */
	d__1 = cotan;
	sin__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
	cos__ = sin__ * cotan;
	tau = 1.;
	if (ABS(cos__) * giant > 1.) {
	    tau = 1. / cos__;
	}
	goto L30;
L20:
	tan__ = v[j] / v[n];
/* Computing 2nd power */
	d__1 = tan__;
	cos__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
	sin__ = cos__ * tan__;
	tau = sin__;
L30:

/*        apply the transformation to v and store the information */
/*        necessary to recover the givens rotation. */

	v[n] = sin__ * v[j] + cos__ * v[n];
	v[j] = tau;

/*        apply the transformation to s and extend the spike in w. */

	l = jj;
	i__2 = m;
	for (i__ = j; i__ <= i__2; ++i__) {
	    temp = cos__ * s[l] - sin__ * w[i__];
	    w[i__] = sin__ * s[l] + cos__ * w[i__];
	    s[l] = temp;
	    ++l;
/* L40: */
	}
L50:
/* L60: */
	;
    }
L70:

/*     add the spike from the rank 1 update to w. */

    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__] += v[n] * u[i__];
/* L80: */
    }

/*     eliminate the spike. */

    *sing = FALSE;
    if (nm1 < 1) {
	goto L140;
    }
    i__1 = nm1;
    for (j = 1; j <= i__1; ++j) {
	if (w[j] == 0.) {
	    goto L120;
	}

/*        determine a givens rotation which eliminates the */
/*        j-th element of the spike. */

	if ((d__1 = s[jj], ABS(d__1)) >= (d__2 = w[j], ABS(d__2))) {
	    goto L90;
	}
	cotan = s[jj] / w[j];
/* Computing 2nd power */
	d__1 = cotan;
	sin__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
	cos__ = sin__ * cotan;
	tau = 1.;
	if (ABS(cos__) * giant > 1.) {
	    tau = 1. / cos__;
	}
	goto L100;
L90:
	tan__ = w[j] / s[jj];
/* Computing 2nd power */
	d__1 = tan__;
	cos__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
	sin__ = cos__ * tan__;
	tau = sin__;
L100:

/*        apply the transformation to s and reduce the spike in w. */

	l = jj;
	i__2 = m;
	for (i__ = j; i__ <= i__2; ++i__) {
	    temp = cos__ * s[l] + sin__ * w[i__];
	    w[i__] = -sin__ * s[l] + cos__ * w[i__];
	    s[l] = temp;
	    ++l;
/* L110: */
	}

/*        store the information necessary to recover the */
/*        givens rotation. */

	w[j] = tau;
L120:

/*        test for zero diagonal elements in the output s. */

	if (s[jj] == 0.) {
	    *sing = TRUE;
	}
	jj += m - j + 1;
/* L130: */
    }
L140:

/*     move w back into the last column of the output s. */

    l = jj;
    i__1 = m;
    for (i__ = n; i__ <= i__1; ++i__) {
	s[l] = w[i__];
	++l;
/* L150: */
    }
    if (s[jj] == 0.) {
	*sing = TRUE;
    }
    return;

/*     last card of subroutine r1updt. */

} /* r1updt_ */

