/* ssyequb.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "pnl/pnl_f2c.h"

/* Table of constant values */

static int c__1 = 1;

 int ssyequb_(char *uplo, int *n, float *a, int *lda, 
	float *s, float *scond, float *amax, float *work, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;
    float r__1, r__2, r__3;

    /* Builtin functions */
    double sqrt(double), log(double), pow_ri(float *, int *);

    /* Local variables */
    float d__;
    int i__, j;
    float t, u, c0, c1, c2, si;
    int up;
    float avg, std, tol, base;
    int iter;
    float smin, smax, scale;
    extern int lsame_(char *, char *);
    float sumsq;
    extern double slamch_(char *);
    extern  int xerbla_(char *, int *);
    float bignum;
    extern  int slassq_(int *, float *, int *, float *, 
	    float *);
    float smlnum;


/*     -- LAPACK routine (version 3.2)                                 -- */
/*     -- Contributed by James Demmel, Deaglan Halligan, Yozo Hida and -- */
/*     -- Jason Riedy of Univ. of California Berkeley.                 -- */
/*     -- November 2008                                                -- */

/*     -- LAPACK is a software package provided by Univ. of Tennessee, -- */
/*     -- Univ. of California Berkeley and NAG Ltd.                    -- */

/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SSYEQUB computes row and column scalings intended to equilibrate a */
/*  symmetric matrix A and reduce its condition number */
/*  (with respect to the two-norm).  S contains the scale factors, */
/*  S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with */
/*  elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This */
/*  choice of S puts the condition number of B within a factor N of the */
/*  smallest possible condition number over all possible diagonal */
/*  scalings. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0. */

/*  A       (input) REAL array, dimension (LDA,N) */
/*          The N-by-N symmetric matrix whose scaling */
/*          factors are to be computed.  Only the diagonal elements of A */
/*          are referenced. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= MAX(1,N). */

/*  S       (output) REAL array, dimension (N) */
/*          If INFO = 0, S contains the scale factors for A. */

/*  SCOND   (output) REAL */
/*          If INFO = 0, S contains the ratio of the smallest S(i) to */
/*          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too */
/*          large nor too small, it is not worth scaling by S. */

/*  AMAX    (output) REAL */
/*          Absolute value of largest matrix element.  If AMAX is very */
/*          close to overflow or very close to underflow, the matrix */
/*          should be scaled. */
/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
/*          > 0:  if INFO = i, the i-th diagonal element is nonpositive. */

/*  Further Details */
/*  ======= ======= */

/*  Reference: Livne, O.E. and Golub, G.H., "Scaling by Binormalization", */
/*  Numerical Algorithms, vol. 35, no. 1, pp. 97-120, January 2004. */
/*  DOI 10.1023/B:NUMA.0000016606.32820.69 */
/*  Tech report version: http://ruready.utah.edu/archive/papers/bin.pdf */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --s;
    --work;

    /* Function Body */
    *info = 0;
    if (! (lsame_(uplo, "U") || lsame_(uplo, "L"))) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < MAX(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SSYEQUB", &i__1);
	return 0;
    }
    up = lsame_(uplo, "U");
    *amax = 0.f;

/*     Quick return if possible. */

    if (*n == 0) {
	*scond = 1.f;
	return 0;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s[i__] = 0.f;
    }
    *amax = 0.f;
    if (up) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = j - 1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
		r__2 = s[i__], r__3 = (r__1 = a[i__ + j * a_dim1], ABS(r__1))
			;
		s[i__] = MAX(r__2,r__3);
/* Computing MAX */
		r__2 = s[j], r__3 = (r__1 = a[i__ + j * a_dim1], ABS(r__1));
		s[j] = MAX(r__2,r__3);
/* Computing MAX */
		r__2 = *amax, r__3 = (r__1 = a[i__ + j * a_dim1], ABS(r__1));
		*amax = MAX(r__2,r__3);
	    }
/* Computing MAX */
	    r__2 = s[j], r__3 = (r__1 = a[j + j * a_dim1], ABS(r__1));
	    s[j] = MAX(r__2,r__3);
/* Computing MAX */
	    r__2 = *amax, r__3 = (r__1 = a[j + j * a_dim1], ABS(r__1));
	    *amax = MAX(r__2,r__3);
	}
    } else {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	    r__2 = s[j], r__3 = (r__1 = a[j + j * a_dim1], ABS(r__1));
	    s[j] = MAX(r__2,r__3);
/* Computing MAX */
	    r__2 = *amax, r__3 = (r__1 = a[j + j * a_dim1], ABS(r__1));
	    *amax = MAX(r__2,r__3);
	    i__2 = *n;
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
/* Computing MAX */
		r__2 = s[i__], r__3 = (r__1 = a[i__ + j * a_dim1], ABS(r__1))
			;
		s[i__] = MAX(r__2,r__3);
/* Computing MAX */
		r__2 = s[j], r__3 = (r__1 = a[i__ + j * a_dim1], ABS(r__1));
		s[j] = MAX(r__2,r__3);
/* Computing MAX */
		r__2 = *amax, r__3 = (r__1 = a[i__ + j * a_dim1], ABS(r__1));
		*amax = MAX(r__2,r__3);
	    }
	}
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	s[j] = 1.f / s[j];
    }
    tol = 1.f / sqrt(*n * 2.f);
    for (iter = 1; iter <= 100; ++iter) {
	scale = 0.f;
	sumsq = 0.f;
/*       BETA = |A|S */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    work[i__] = 0.f;
	}
	if (up) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    t = (r__1 = a[i__ + j * a_dim1], ABS(r__1));
		    work[i__] += (r__1 = a[i__ + j * a_dim1], ABS(r__1)) * s[
			    j];
		    work[j] += (r__1 = a[i__ + j * a_dim1], ABS(r__1)) * s[
			    i__];
		}
		work[j] += (r__1 = a[j + j * a_dim1], ABS(r__1)) * s[j];
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		work[j] += (r__1 = a[j + j * a_dim1], ABS(r__1)) * s[j];
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    t = (r__1 = a[i__ + j * a_dim1], ABS(r__1));
		    work[i__] += (r__1 = a[i__ + j * a_dim1], ABS(r__1)) * s[
			    j];
		    work[j] += (r__1 = a[i__ + j * a_dim1], ABS(r__1)) * s[
			    i__];
		}
	    }
	}
/*       avg = s^T beta / n */
	avg = 0.f;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    avg += s[i__] * work[i__];
	}
	avg /= *n;
	std = 0.f;
	i__1 = *n * 3;
	for (i__ = (*n << 1) + 1; i__ <= i__1; ++i__) {
	    work[i__] = s[i__ - (*n << 1)] * work[i__ - (*n << 1)] - avg;
	}
	slassq_(n, &work[(*n << 1) + 1], &c__1, &scale, &sumsq);
	std = scale * sqrt(sumsq / *n);
	if (std < tol * avg) {
	    goto L999;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    t = (r__1 = a[i__ + i__ * a_dim1], ABS(r__1));
	    si = s[i__];
	    c2 = (*n - 1) * t;
	    c1 = (*n - 2) * (work[i__] - t * si);
	    c0 = -(t * si) * si + work[i__] * 2 * si - *n * avg;
	    d__ = c1 * c1 - c0 * 4 * c2;
	    if (d__ <= 0.f) {
		*info = -1;
		return 0;
	    }
	    si = c0 * -2 / (c1 + sqrt(d__));
	    d__ = si - s[i__];
	    u = 0.f;
	    if (up) {
		i__2 = i__;
		for (j = 1; j <= i__2; ++j) {
		    t = (r__1 = a[j + i__ * a_dim1], ABS(r__1));
		    u += s[j] * t;
		    work[j] += d__ * t;
		}
		i__2 = *n;
		for (j = i__ + 1; j <= i__2; ++j) {
		    t = (r__1 = a[i__ + j * a_dim1], ABS(r__1));
		    u += s[j] * t;
		    work[j] += d__ * t;
		}
	    } else {
		i__2 = i__;
		for (j = 1; j <= i__2; ++j) {
		    t = (r__1 = a[i__ + j * a_dim1], ABS(r__1));
		    u += s[j] * t;
		    work[j] += d__ * t;
		}
		i__2 = *n;
		for (j = i__ + 1; j <= i__2; ++j) {
		    t = (r__1 = a[j + i__ * a_dim1], ABS(r__1));
		    u += s[j] * t;
		    work[j] += d__ * t;
		}
	    }
	    avg += (u + work[i__]) * d__ / *n;
	    s[i__] = si;
	}
    }
L999:
    smlnum = slamch_("SAFEMIN");
    bignum = 1.f / smlnum;
    smin = bignum;
    smax = 0.f;
    t = 1.f / sqrt(avg);
    base = slamch_("B");
    u = 1.f / log(base);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = (int) (u * log(s[i__] * t));
	s[i__] = pow_ri(&base, &i__2);
/* Computing MIN */
	r__1 = smin, r__2 = s[i__];
	smin = MIN(r__1,r__2);
/* Computing MAX */
	r__1 = smax, r__2 = s[i__];
	smax = MAX(r__1,r__2);
    }
    *scond = MAX(smin,smlnum) / MIN(smax,bignum);

    return 0;
} /* ssyequb_ */
