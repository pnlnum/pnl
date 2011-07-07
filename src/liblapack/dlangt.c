/* dlangt.f -- translated by f2c (version 20061008).
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

double dlangt_(char *norm, int *n, double *dl, double *d__, 
	double *du)
{
    /* System generated locals */
    int i__1;
    double ret_val, d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    int i__;
    double sum, scale;
    extern int lsame_(char *, char *);
    double anorm;
    extern  int dlassq_(int *, double *, int *, 
	    double *, double *);


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLANGT  returns the value of the one norm,  or the Frobenius norm, or */
/*  the  infinity norm,  or the  element of  largest absolute value  of a */
/*  float tridiagonal matrix A. */

/*  Description */
/*  =========== */

/*  DLANGT returns the value */

/*     DLANGT = ( MAX(ABS(A(i,j))), NORM = 'M' or 'm' */
/*              ( */
/*              ( norm1(A),         NORM = '1', 'O' or 'o' */
/*              ( */
/*              ( normI(A),         NORM = 'I' or 'i' */
/*              ( */
/*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e' */

/*  where  norm1  denotes the  one norm of a matrix (maximum column sum), */
/*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and */
/*  normF  denotes the  Frobenius norm of a matrix (square root of sum of */
/*  squares).  Note that  MAX(ABS(A(i,j)))  is not a consistent matrix norm. */

/*  Arguments */
/*  ========= */

/*  NORM    (input) CHARACTER*1 */
/*          Specifies the value to be returned in DLANGT as described */
/*          above. */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0.  When N = 0, DLANGT is */
/*          set to zero. */

/*  DL      (input) DOUBLE PRECISION array, dimension (N-1) */
/*          The (n-1) sub-diagonal elements of A. */

/*  D       (input) DOUBLE PRECISION array, dimension (N) */
/*          The diagonal elements of A. */

/*  DU      (input) DOUBLE PRECISION array, dimension (N-1) */
/*          The (n-1) super-diagonal elements of A. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --du;
    --d__;
    --dl;

    /* Function Body */
    if (*n <= 0) {
	anorm = 0.;
    } else if (lsame_(norm, "M")) {

/*        Find MAX(ABS(A(i,j))). */

	anorm = (d__1 = d__[*n], ABS(d__1));
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__2 = anorm, d__3 = (d__1 = dl[i__], ABS(d__1));
	    anorm = MAX(d__2,d__3);
/* Computing MAX */
	    d__2 = anorm, d__3 = (d__1 = d__[i__], ABS(d__1));
	    anorm = MAX(d__2,d__3);
/* Computing MAX */
	    d__2 = anorm, d__3 = (d__1 = du[i__], ABS(d__1));
	    anorm = MAX(d__2,d__3);
/* L10: */
	}
    } else if (lsame_(norm, "O") || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

	if (*n == 1) {
	    anorm = ABS(d__[1]);
	} else {
/* Computing MAX */
	    d__3 = ABS(d__[1]) + ABS(dl[1]), d__4 = (d__1 = d__[*n], ABS(d__1)
		    ) + (d__2 = du[*n - 1], ABS(d__2));
	    anorm = MAX(d__3,d__4);
	    i__1 = *n - 1;
	    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing MAX */
		d__4 = anorm, d__5 = (d__1 = d__[i__], ABS(d__1)) + (d__2 = 
			dl[i__], ABS(d__2)) + (d__3 = du[i__ - 1], ABS(d__3));
		anorm = MAX(d__4,d__5);
/* L20: */
	    }
	}
    } else if (lsame_(norm, "I")) {

/*        Find normI(A). */

	if (*n == 1) {
	    anorm = ABS(d__[1]);
	} else {
/* Computing MAX */
	    d__3 = ABS(d__[1]) + ABS(du[1]), d__4 = (d__1 = d__[*n], ABS(d__1)
		    ) + (d__2 = dl[*n - 1], ABS(d__2));
	    anorm = MAX(d__3,d__4);
	    i__1 = *n - 1;
	    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing MAX */
		d__4 = anorm, d__5 = (d__1 = d__[i__], ABS(d__1)) + (d__2 = 
			du[i__], ABS(d__2)) + (d__3 = dl[i__ - 1], ABS(d__3));
		anorm = MAX(d__4,d__5);
/* L30: */
	    }
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.;
	sum = 1.;
	dlassq_(n, &d__[1], &c__1, &scale, &sum);
	if (*n > 1) {
	    i__1 = *n - 1;
	    dlassq_(&i__1, &dl[1], &c__1, &scale, &sum);
	    i__1 = *n - 1;
	    dlassq_(&i__1, &du[1], &c__1, &scale, &sum);
	}
	anorm = scale * sqrt(sum);
    }

    ret_val = anorm;
    return ret_val;

/*     End of DLANGT */

} /* dlangt_ */
