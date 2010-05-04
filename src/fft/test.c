/* 
 * This file is based largely on the following software distribution:
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 *                              FFTPACK
 * 
 * Reference                                                                                                                        
 *    P.N. Swarztrauber, Vectorizing the FFTs, in Parallel Computations
 *    (G. Rodrigue, ed.), Academic Press, 1982, pp. 51--83.                                                                                                                   
 * 
 *     http://www.netlib.org/fftpack/
 * 
 * Updated to single, double, and extended precision,
 * and translated to ISO-Standard C/C++ (without aliasing)
 * on 10 October 2005 by Andrew Fernandes <andrew_AT_fernandes.org>
 * 
 *                   Version 4  April 1985
 * 
 *      A Package of Fortran Subprograms for the Fast Fourier
 *       Transform of Periodic and other Symmetric Sequences
 * 
 *                          by
 * 
 *                   Paul N Swarztrauber
 * 
 *   National Center for Atmospheric Research, Boulder, Colorado 80307,
 * 
 *    which is sponsored by the National Science Foundation
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * There appears to be no explicit license for FFTPACK. However, the
 * package has been incorporated verbatim into a large number of software
 * systems over the years with numerous types of license without complaint
 * from the original author; therefore it would appear
 * that the code is effectively public domain. If you are in doubt,
 * however, you will need to contact the author or the  National Center
 * for Atmospheric Research to be sure.
 * 
 * All the changes from the original FFTPACK to the current file
 * fall under the following BSD-style open-source license:
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * Copyright (c) 2005, Andrew Fernandes (andrew@fernandes.org);
 * All rights reserved.
 *  
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  
 * - Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 
 * - Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 
 * - Neither the name of the North Carolina State University nor the
 * names of its contributors may be used to endorse or promote products
 * derived from this software without specific prior written permission.
 *  
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */
#include "f2c.h"

/* test.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* Table of constant values */

static integer c_n1 = -1;
static integer c__1 = 1;

/* Main program */ int MAIN__(void)
{
    /* Initialized data */

    static integer nd[10] = { 120,54,49,32,4,3,2 };

    /* Format strings */
    static char fmt_1001[] = "(\0020N\002,i5,\002 RFFTF  \002,e10.3,\002 RFF"
	    "TB  \002,e10.3,\002 RFFTFB \002,e10.3,\002 SINT   \002,e10.3,"
	    "\002 SINTFB \002,e10.3,\002 COST   \002,e10.3/7x,\002 COSTFB "
	    "\002,e10.3,\002 SINQF  \002,e10.3,\002 SINQB  \002,e10.3,\002 SI"
	    "NQFB \002,e10.3,\002 COSQF  \002,e10.3,\002 COSQB  \002,e10.3/7x,"
	    "\002 COSQFB \002,e10.3,\002 DEZF   \002,e10.3,\002 DEZB   \002,e"
	    "10.3,\002 DEZFB  \002,e10.3,\002 CFFTF  \002,e10.3,\002 CFFTB "
	    " \002,e10.3/7x,\002 CFFTFB \002,e10.3)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double sqrt(doublereal), sin(doublereal), cos(doublereal);
    integer pow_ii(integer *, integer *);
    double atan(doublereal), z_abs(doublecomplex *);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    doublereal a[100], b[100];
    integer i__, j, k, n;
    doublereal w[2000], x[200], y[200], ah[100], bh[100], cf, fn, dt, pi;
    doublecomplex cx[200], cy[200];
    doublereal xh[200];
    integer nz, nm1, np1, ns2;
    doublereal arg, tfn, tpi;
    integer nns;
    doublereal sum, arg1, arg2;
    integer ns2m;
    doublereal sum1, sum2, dcfb;
    integer ifac[64], modn;
    doublereal rftb, rftf;
    extern /* Subroutine */ void cost(integer *, doublereal *, doublereal *, 
	    integer *), sint(integer *, doublereal *, doublereal *, integer *
	    );
    doublereal dezb1, dezf1, sqrt2;
    extern /* Subroutine */ void cfftb(integer *, doublecomplex *, doublereal 
	    *, integer *), cfftf(integer *, doublecomplex *, doublereal *, 
	    integer *);
    doublereal dezfb;
    extern /* Subroutine */ void cffti(integer *, doublereal *, integer *), 
	    rfftb(integer *, doublereal *, doublereal *, integer *);
    doublereal rftfb;
    extern /* Subroutine */ void rfftf(integer *, doublereal *, doublereal *, 
	    integer *), cosqb(integer *, doublereal *, doublereal *, integer 
	    *), rffti(integer *, doublereal *, integer *), cosqf(integer *, 
	    doublereal *, doublereal *, integer *), sinqb(integer *, 
	    doublereal *, doublereal *, integer *), cosqi(integer *, 
	    doublereal *, integer *), sinqf(integer *, doublereal *, 
	    doublereal *, integer *), costi(integer *, doublereal *, integer 
	    *);
    doublereal azero;
    extern /* Subroutine */ void sinqi(integer *, doublereal *, integer *), 
	    sinti(integer *, doublereal *, integer *);
    doublereal costt, sintt, dcfftb, dcfftf, cosqfb, costfb;
    extern /* Subroutine */ void ezfftb(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *);
    doublereal sinqfb;
    extern /* Subroutine */ void ezfftf(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *);
    doublereal sintfb;
    extern /* Subroutine */ void ezffti(integer *, doublereal *, integer *);
    doublereal azeroh, cosqbt, cosqft, sinqbt, sinqft;

    /* Fortran I/O blocks */
    static cilist io___58 = { 0, 6, 0, fmt_1001, 0 };



/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*                       VERSION 4  APRIL 1985 */

/*                         A TEST DRIVER FOR */
/*          A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE FAST FOURIER */
/*           TRANSFORM OF PERIODIC AND OTHER SYMMETRIC SEQUENCES */

/*                              BY */

/*                       PAUL N SWARZTRAUBER */

/*       NATIONAL CENTER FOR ATMOSPHERIC RESEARCH  BOULDER,COLORADO 80307 */

/*        WHICH IS SPONSORED BY THE NATIONAL SCIENCE FOUNDATION */

/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/*             THIS PROGRAM TESTS THE PACKAGE OF FAST FOURIER */
/*     TRANSFORMS FOR BOTH COMPLEX AND REAL PERIODIC SEQUENCES AND */
/*     CERTIAN OTHER SYMMETRIC SEQUENCES THAT ARE LISTED BELOW. */

/*     1.   RFFTI     INITIALIZE  RFFTF AND RFFTB */
/*     2.   RFFTF     FORWARD TRANSFORM OF A REAL PERIODIC SEQUENCE */
/*     3.   RFFTB     BACKWARD TRANSFORM OF A REAL COEFFICIENT ARRAY */

/*     4.   EZFFTI    INITIALIZE EZFFTF AND EZFFTB */
/*     5.   EZFFTF    A SIMPLIFIED REAL PERIODIC FORWARD TRANSFORM */
/*     6.   EZFFTB    A SIMPLIFIED REAL PERIODIC BACKWARD TRANSFORM */

/*     7.   SINTI     INITIALIZE SINT */
/*     8.   SINT      SINE TRANSFORM OF A REAL ODD SEQUENCE */

/*     9.   COSTI     INITIALIZE COST */
/*     10.  COST      COSINE TRANSFORM OF A REAL EVEN SEQUENCE */

/*     11.  SINQI     INITIALIZE SINQF AND SINQB */
/*     12.  SINQF     FORWARD SINE TRANSFORM WITH ODD WAVE NUMBERS */
/*     13.  SINQB     UNNORMALIZED INVERSE OF SINQF */

/*     14.  COSQI     INITIALIZE COSQF AND COSQB */
/*     15.  COSQF     FORWARD COSINE TRANSFORM WITH ODD WAVE NUMBERS */
/*     16.  COSQB     UNNORMALIZED INVERSE OF COSQF */

/*     17.  CFFTI     INITIALIZE CFFTF AND CFFTB */
/*     18.  CFFTF     FORWARD TRANSFORM OF A COMPLEX PERIODIC SEQUENCE */
/*     19.  CFFTB     UNNORMALIZED INVERSE OF CFFTF */


    sqrt2 = sqrt(2.0);
    nns = 7;
    i__1 = nns;
    for (nz = 1; nz <= i__1; ++nz) {
	n = nd[nz - 1];
	modn = n % 2;
	fn = (real) n;
	tfn = fn + fn;
	np1 = n + 1;
	nm1 = n - 1;
	i__2 = np1;
	for (j = 1; j <= i__2; ++j) {
	    x[j - 1] = sin((real) j * sqrt2);
	    y[j - 1] = x[j - 1];
	    xh[j - 1] = x[j - 1];
/* L101: */
	}

/*     TEST SUBROUTINES RFFTI,RFFTF AND RFFTB */

	rffti(&n, w, ifac);
	pi = 3.141592653589793238462643383279502884197169399375108209749445923;
	dt = (pi + pi) / fn;
	ns2 = (n + 1) / 2;
	if (ns2 < 2) {
	    goto L104;
	}
	i__2 = ns2;
	for (k = 2; k <= i__2; ++k) {
	    sum1 = 0.0;
	    sum2 = 0.0;
	    arg = (real) (k - 1) * dt;
	    i__3 = n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		arg1 = (real) (i__ - 1) * arg;
		sum1 += x[i__ - 1] * cos(arg1);
		sum2 += x[i__ - 1] * sin(arg1);
/* L102: */
	    }
	    y[(k << 1) - 3] = sum1;
	    y[(k << 1) - 2] = -sum2;
/* L103: */
	}
L104:
	sum1 = 0.0;
	sum2 = 0.0;
	i__2 = nm1;
	for (i__ = 1; i__ <= i__2; i__ += 2) {
	    sum1 += x[i__ - 1];
	    sum2 += x[i__];
/* L105: */
	}
	if (modn == 1) {
	    sum1 += x[n - 1];
	}
	y[0] = sum1 + sum2;
	if (modn == 0) {
	    y[n - 1] = sum1 - sum2;
	}
	rfftf(&n, x, w, ifac);
	rftf = 0.0;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = rftf, d__3 = (d__1 = x[i__ - 1] - y[i__ - 1], abs(d__1));
	    rftf = max(d__2,d__3);
	    x[i__ - 1] = xh[i__ - 1];
/* L106: */
	}
	rftf /= fn;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sum = x[0] * 0.5;
	    arg = (real) (i__ - 1) * dt;
	    if (ns2 < 2) {
		goto L108;
	    }
	    i__3 = ns2;
	    for (k = 2; k <= i__3; ++k) {
		arg1 = (real) (k - 1) * arg;
		sum = sum + x[(k << 1) - 3] * cos(arg1) - x[(k << 1) - 2] * 
			sin(arg1);
/* L107: */
	    }
L108:
	    if (modn == 0) {
		i__3 = i__ - 1;
		sum += (real) pow_ii(&c_n1, &i__3) * 0.5 * x[n - 1];
	    }
	    y[i__ - 1] = sum + sum;
/* L109: */
	}
	rfftb(&n, x, w, ifac);
	rftb = 0.0;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = rftb, d__3 = (d__1 = x[i__ - 1] - y[i__ - 1], abs(d__1));
	    rftb = max(d__2,d__3);
	    x[i__ - 1] = xh[i__ - 1];
	    y[i__ - 1] = xh[i__ - 1];
/* L110: */
	}
	rfftb(&n, y, w, ifac);
	rfftf(&n, y, w, ifac);
	cf = 1.0 / fn;
	rftfb = 0.0;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = rftfb, d__3 = (d__1 = cf * y[i__ - 1] - x[i__ - 1], abs(
		    d__1));
	    rftfb = max(d__2,d__3);
/* L111: */
	}

/*     TEST SUBROUTINES SINTI AND SINT */

	dt = pi / fn;
	i__2 = nm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x[i__ - 1] = xh[i__ - 1];
/* L112: */
	}
	i__2 = nm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__ - 1] = 0.0;
	    arg1 = (real) i__ * dt;
	    i__3 = nm1;
	    for (k = 1; k <= i__3; ++k) {
		y[i__ - 1] += x[k - 1] * sin((real) k * arg1);
/* L113: */
	    }
	    y[i__ - 1] += y[i__ - 1];
/* L114: */
	}
	sinti(&nm1, w, ifac);
	sint(&nm1, x, w, ifac);
	cf = 0.5 / fn;
	sintt = 0.0;
	i__2 = nm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = sintt, d__3 = (d__1 = x[i__ - 1] - y[i__ - 1], abs(d__1));
	    sintt = max(d__2,d__3);
	    x[i__ - 1] = xh[i__ - 1];
	    y[i__ - 1] = x[i__ - 1];
/* L115: */
	}
	sintt = cf * sintt;
	sint(&nm1, x, w, ifac);
	sint(&nm1, x, w, ifac);
	sintfb = 0.0;
	i__2 = nm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = sintfb, d__3 = (d__1 = cf * x[i__ - 1] - y[i__ - 1], abs(
		    d__1));
	    sintfb = max(d__2,d__3);
/* L116: */
	}

/*     TEST SUBROUTINES COSTI AND COST */

	i__2 = np1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x[i__ - 1] = xh[i__ - 1];
/* L117: */
	}
	i__2 = np1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + 1;
	    y[i__ - 1] = (x[0] + (real) pow_ii(&c_n1, &i__3) * x[n]) * 0.5;
	    arg = (real) (i__ - 1) * dt;
	    i__3 = n;
	    for (k = 2; k <= i__3; ++k) {
		y[i__ - 1] += x[k - 1] * cos((real) (k - 1) * arg);
/* L118: */
	    }
	    y[i__ - 1] += y[i__ - 1];
/* L119: */
	}
	costi(&np1, w, ifac);
	cost(&np1, x, w, ifac);
	costt = 0.0;
	i__2 = np1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = costt, d__3 = (d__1 = x[i__ - 1] - y[i__ - 1], abs(d__1));
	    costt = max(d__2,d__3);
	    x[i__ - 1] = xh[i__ - 1];
	    y[i__ - 1] = xh[i__ - 1];
/* L120: */
	}
	costt = cf * costt;
	cost(&np1, x, w, ifac);
	cost(&np1, x, w, ifac);
	costfb = 0.0;
	i__2 = np1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = costfb, d__3 = (d__1 = cf * x[i__ - 1] - y[i__ - 1], abs(
		    d__1));
	    costfb = max(d__2,d__3);
/* L121: */
	}

/*     TEST SUBROUTINES SINQI,SINQF AND SINQB */

	cf = 0.25 / fn;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__ - 1] = xh[i__ - 1];
/* L122: */
	}
	dt = pi / (fn + fn);
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x[i__ - 1] = 0.0;
	    arg = dt * (real) i__;
	    i__3 = n;
	    for (k = 1; k <= i__3; ++k) {
		x[i__ - 1] += y[k - 1] * sin((real) (k + k - 1) * arg);
/* L123: */
	    }
	    x[i__ - 1] *= 4.0;
/* L124: */
	}
	sinqi(&n, w, ifac);
	sinqb(&n, y, w, ifac);
	sinqbt = 0.0;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = sinqbt, d__3 = (d__1 = y[i__ - 1] - x[i__ - 1], abs(d__1));
	    sinqbt = max(d__2,d__3);
	    x[i__ - 1] = xh[i__ - 1];
/* L125: */
	}
	sinqbt = cf * sinqbt;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    arg = (real) (i__ + i__ - 1) * dt;
	    i__3 = i__ + 1;
	    y[i__ - 1] = (real) pow_ii(&c_n1, &i__3) * 0.5 * x[n - 1];
	    i__3 = nm1;
	    for (k = 1; k <= i__3; ++k) {
		y[i__ - 1] += x[k - 1] * sin((real) k * arg);
/* L126: */
	    }
	    y[i__ - 1] += y[i__ - 1];
/* L127: */
	}
	sinqf(&n, x, w, ifac);
	sinqft = 0.0;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = sinqft, d__3 = (d__1 = x[i__ - 1] - y[i__ - 1], abs(d__1));
	    sinqft = max(d__2,d__3);
	    y[i__ - 1] = xh[i__ - 1];
	    x[i__ - 1] = xh[i__ - 1];
/* L128: */
	}
	sinqf(&n, y, w, ifac);
	sinqb(&n, y, w, ifac);
	sinqfb = 0.0;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = sinqfb, d__3 = (d__1 = cf * y[i__ - 1] - x[i__ - 1], abs(
		    d__1));
	    sinqfb = max(d__2,d__3);
/* L129: */
	}

/*     TEST SUBROUTINES COSQI,COSQF AND COSQB */

	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__ - 1] = xh[i__ - 1];
/* L130: */
	}
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x[i__ - 1] = 0.0;
	    arg = (real) (i__ - 1) * dt;
	    i__3 = n;
	    for (k = 1; k <= i__3; ++k) {
		x[i__ - 1] += y[k - 1] * cos((real) (k + k - 1) * arg);
/* L131: */
	    }
	    x[i__ - 1] *= 4.0;
/* L132: */
	}
	cosqi(&n, w, ifac);
	cosqb(&n, y, w, ifac);
	cosqbt = 0.0;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = cosqbt, d__3 = (d__1 = x[i__ - 1] - y[i__ - 1], abs(d__1));
	    cosqbt = max(d__2,d__3);
	    x[i__ - 1] = xh[i__ - 1];
/* L133: */
	}
	cosqbt = cf * cosqbt;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__ - 1] = x[0] * 0.5;
	    arg = (real) (i__ + i__ - 1) * dt;
	    i__3 = n;
	    for (k = 2; k <= i__3; ++k) {
		y[i__ - 1] += x[k - 1] * cos((real) (k - 1) * arg);
/* L134: */
	    }
	    y[i__ - 1] += y[i__ - 1];
/* L135: */
	}
	cosqf(&n, x, w, ifac);
	cosqft = 0.0;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = cosqft, d__3 = (d__1 = y[i__ - 1] - x[i__ - 1], abs(d__1));
	    cosqft = max(d__2,d__3);
	    x[i__ - 1] = xh[i__ - 1];
	    y[i__ - 1] = xh[i__ - 1];
/* L136: */
	}
	cosqft = cf * cosqft;
	cosqb(&n, x, w, ifac);
	cosqf(&n, x, w, ifac);
	cosqfb = 0.0;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = cosqfb, d__3 = (d__1 = cf * x[i__ - 1] - y[i__ - 1], abs(
		    d__1));
	    cosqfb = max(d__2,d__3);
/* L137: */
	}

/*     TEST PROGRAMS EZFFTI,EZFFTF,EZFFTB */

	ezffti(&n, w, ifac);
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x[i__ - 1] = xh[i__ - 1];
/* L138: */
	}
	tpi = atan(1.0) * 8.0;
	dt = tpi / (real) n;
	ns2 = (n + 1) / 2;
	cf = 2.0 / (real) n;
	ns2m = ns2 - 1;
	if (ns2m <= 0) {
	    goto L141;
	}
	i__2 = ns2m;
	for (k = 1; k <= i__2; ++k) {
	    sum1 = 0.0;
	    sum2 = 0.0;
	    arg = (real) k * dt;
	    i__3 = n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		arg1 = (real) (i__ - 1) * arg;
		sum1 += x[i__ - 1] * cos(arg1);
		sum2 += x[i__ - 1] * sin(arg1);
/* L139: */
	    }
	    a[k - 1] = cf * sum1;
	    b[k - 1] = cf * sum2;
/* L140: */
	}
L141:
	nm1 = n - 1;
	sum1 = 0.0;
	sum2 = 0.0;
	i__2 = nm1;
	for (i__ = 1; i__ <= i__2; i__ += 2) {
	    sum1 += x[i__ - 1];
	    sum2 += x[i__];
/* L142: */
	}
	if (modn == 1) {
	    sum1 += x[n - 1];
	}
	azero = cf * 0.5 * (sum1 + sum2);
	if (modn == 0) {
	    a[ns2 - 1] = cf * 0.5 * (sum1 - sum2);
	}
	ezfftf(&n, x, &azeroh, ah, bh, w, ifac);
	dezf1 = (d__1 = azeroh - azero, abs(d__1));
	if (modn == 0) {
/* Computing MAX */
	    d__2 = dezf1, d__3 = (d__1 = a[ns2 - 1] - ah[ns2 - 1], abs(d__1));
	    dezf1 = max(d__2,d__3);
	}
	if (ns2m <= 0) {
	    goto L144;
	}
	i__2 = ns2m;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__3 = dezf1, d__4 = (d__1 = ah[i__ - 1] - a[i__ - 1], abs(d__1)),
		     d__3 = max(d__3,d__4), d__4 = (d__2 = bh[i__ - 1] - b[
		    i__ - 1], abs(d__2));
	    dezf1 = max(d__3,d__4);
/* L143: */
	}
L144:
	ns2 = n / 2;
	if (modn == 0) {
	    b[ns2 - 1] = 0.0;
	}
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sum = azero;
	    arg1 = (real) (i__ - 1) * dt;
	    i__3 = ns2;
	    for (k = 1; k <= i__3; ++k) {
		arg2 = (real) k * arg1;
		sum = sum + a[k - 1] * cos(arg2) + b[k - 1] * sin(arg2);
/* L145: */
	    }
	    x[i__ - 1] = sum;
/* L146: */
	}
	ezfftb(&n, y, &azero, a, b, w, ifac);
	dezb1 = 0.0;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = dezb1, d__3 = (d__1 = x[i__ - 1] - y[i__ - 1], abs(d__1));
	    dezb1 = max(d__2,d__3);
	    x[i__ - 1] = xh[i__ - 1];
/* L147: */
	}
	ezfftf(&n, x, &azero, a, b, w, ifac);
	ezfftb(&n, y, &azero, a, b, w, ifac);
	dezfb = 0.0;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = dezfb, d__3 = (d__1 = x[i__ - 1] - y[i__ - 1], abs(d__1));
	    dezfb = max(d__2,d__3);
/* L148: */
	}

/*     TEST  CFFTI,CFFTF,CFFTB */

	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ - 1;
	    d__1 = cos(sqrt2 * (real) i__);
	    d__2 = sin(sqrt2 * (real) (i__ * i__));
	    z__1.r = d__1, z__1.i = d__2;
	    cx[i__3].r = z__1.r, cx[i__3].i = z__1.i;
/* L149: */
	}
	dt = (pi + pi) / fn;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    arg1 = -((real) (i__ - 1)) * dt;
	    i__3 = i__ - 1;
	    cy[i__3].r = 0.0, cy[i__3].i = 0.0;
	    i__3 = n;
	    for (k = 1; k <= i__3; ++k) {
		arg2 = (real) (k - 1) * arg1;
		i__4 = i__ - 1;
		i__5 = i__ - 1;
		d__1 = cos(arg2);
		d__2 = sin(arg2);
		z__3.r = d__1, z__3.i = d__2;
		i__6 = k - 1;
		z__2.r = z__3.r * cx[i__6].r - z__3.i * cx[i__6].i, z__2.i = 
			z__3.r * cx[i__6].i + z__3.i * cx[i__6].r;
		z__1.r = cy[i__5].r + z__2.r, z__1.i = cy[i__5].i + z__2.i;
		cy[i__4].r = z__1.r, cy[i__4].i = z__1.i;
/* L150: */
	    }
/* L151: */
	}
	cffti(&n, w, ifac);
	cfftf(&n, cx, w, ifac);
	dcfftf = 0.0;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    i__3 = i__ - 1;
	    i__4 = i__ - 1;
	    z__1.r = cx[i__3].r - cy[i__4].r, z__1.i = cx[i__3].i - cy[i__4]
		    .i;
	    d__1 = dcfftf, d__2 = z_abs(&z__1);
	    dcfftf = max(d__1,d__2);
	    i__3 = i__ - 1;
	    i__4 = i__ - 1;
	    z__1.r = cx[i__4].r / fn, z__1.i = cx[i__4].i / fn;
	    cx[i__3].r = z__1.r, cx[i__3].i = z__1.i;
/* L152: */
	}
	dcfftf /= fn;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    arg1 = (real) (i__ - 1) * dt;
	    i__3 = i__ - 1;
	    cy[i__3].r = 0.0, cy[i__3].i = 0.0;
	    i__3 = n;
	    for (k = 1; k <= i__3; ++k) {
		arg2 = (real) (k - 1) * arg1;
		i__4 = i__ - 1;
		i__5 = i__ - 1;
		d__1 = cos(arg2);
		d__2 = sin(arg2);
		z__3.r = d__1, z__3.i = d__2;
		i__6 = k - 1;
		z__2.r = z__3.r * cx[i__6].r - z__3.i * cx[i__6].i, z__2.i = 
			z__3.r * cx[i__6].i + z__3.i * cx[i__6].r;
		z__1.r = cy[i__5].r + z__2.r, z__1.i = cy[i__5].i + z__2.i;
		cy[i__4].r = z__1.r, cy[i__4].i = z__1.i;
/* L153: */
	    }
/* L154: */
	}
	cfftb(&n, cx, w, ifac);
	dcfftb = 0.0;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    i__3 = i__ - 1;
	    i__4 = i__ - 1;
	    z__1.r = cx[i__3].r - cy[i__4].r, z__1.i = cx[i__3].i - cy[i__4]
		    .i;
	    d__1 = dcfftb, d__2 = z_abs(&z__1);
	    dcfftb = max(d__1,d__2);
	    i__3 = i__ - 1;
	    i__4 = i__ - 1;
	    cx[i__3].r = cy[i__4].r, cx[i__3].i = cy[i__4].i;
/* L155: */
	}
	cf = 1.0 / fn;
	cfftf(&n, cx, w, ifac);
	cfftb(&n, cx, w, ifac);
	dcfb = 0.0;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    i__3 = i__ - 1;
	    z__2.r = cf * cx[i__3].r, z__2.i = cf * cx[i__3].i;
	    i__4 = i__ - 1;
	    z__1.r = z__2.r - cy[i__4].r, z__1.i = z__2.i - cy[i__4].i;
	    d__1 = dcfb, d__2 = z_abs(&z__1);
	    dcfb = max(d__1,d__2);
/* L156: */
	}
	s_wsfe(&io___58);
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&rftf, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&rftb, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&rftfb, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&sintt, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&sintfb, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&costt, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&costfb, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&sinqft, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&sinqbt, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&sinqfb, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&cosqft, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&cosqbt, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&cosqfb, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&dezf1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&dezb1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&dezfb, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&dcfftf, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&dcfftb, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&dcfb, (ftnlen)sizeof(doublereal));
	e_wsfe();
/* L157: */
    }




    return 0;
} /* MAIN__ */

/* Main program alias */ int tstfft_ () { MAIN__ (); return 0; }
