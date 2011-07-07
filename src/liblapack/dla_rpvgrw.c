/* dla_rpvgrw.f -- translated by f2c (version 20061008).
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

double dla_rpvgrw__(int *n, int *ncols, double *a, int *
	lda, double *af, int *ldaf)
{
    /* System generated locals */
    int a_dim1, a_offset, af_dim1, af_offset, i__1, i__2;
    double ret_val, d__1, d__2;

    /* Local variables */
    int i__, j;
    double amax, umax, rpvgrw;


/*     -- LAPACK routine (version 3.2.1)                                 -- */
/*     -- Contributed by James Demmel, Deaglan Halligan, Yozo Hida and -- */
/*     -- Jason Riedy of Univ. of California Berkeley.                 -- */
/*     -- April 2009                                                   -- */

/*     -- LAPACK is a software package provided by Univ. of Tennessee, -- */
/*     -- Univ. of California Berkeley and NAG Ltd.                    -- */

/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLA_RPVGRW computes the reciprocal pivot growth factor */
/*  norm(A)/norm(U). The "max absolute element" norm is used. If this is */
/*  much less than 1, the stability of the LU factorization of the */
/*  (equilibrated) matrix A could be poor. This also means that the */
/*  solution X, estimated condition numbers, and error bounds could be */
/*  unreliable. */

/*  Arguments */
/*  ========= */

/*     N       (input) INTEGER */
/*     The number of linear equations, i.e., the order of the */
/*     matrix A.  N >= 0. */

/*     NCOLS   (input) INTEGER */
/*     The number of columns of the matrix A. NCOLS >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*     On entry, the N-by-N matrix A. */

/*     LDA     (input) INTEGER */
/*     The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     AF      (input) DOUBLE PRECISION array, dimension (LDAF,N) */
/*     The factors L and U from the factorization */
/*     A = P*L*U as computed by DGETRF. */

/*     LDAF    (input) INTEGER */
/*     The leading dimension of the array AF.  LDAF >= MAX(1,N). */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    af_dim1 = *ldaf;
    af_offset = 1 + af_dim1;
    af -= af_offset;

    /* Function Body */
    rpvgrw = 1.;
    i__1 = *ncols;
    for (j = 1; j <= i__1; ++j) {
	amax = 0.;
	umax = 0.;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = (d__1 = a[i__ + j * a_dim1], ABS(d__1));
	    amax = MAX(d__2,amax);
	}
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = (d__1 = af[i__ + j * af_dim1], ABS(d__1));
	    umax = MAX(d__2,umax);
	}
	if (umax != 0.) {
/* Computing MIN */
	    d__1 = amax / umax;
	    rpvgrw = MIN(d__1,rpvgrw);
	}
    }
    ret_val = rpvgrw;
    return ret_val;
} /* dla_rpvgrw__ */
