/* sptrfs.f -- translated by f2c (version 20061008).
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
static float c_b11 = 1.f;

 int sptrfs_(int *n, int *nrhs, float *d__, float *e, 
	float *df, float *ef, float *b, int *ldb, float *x, int *ldx, 
	float *ferr, float *berr, float *work, int *info)
{
    /* System generated locals */
    int b_dim1, b_offset, x_dim1, x_offset, i__1, i__2;
    float r__1, r__2, r__3;

    /* Local variables */
    int i__, j;
    float s, bi, cx, dx, ex;
    int ix, nz;
    float eps, safe1, safe2;
    int count;
    extern  int saxpy_(int *, float *, float *, int *, 
	    float *, int *);
    extern double slamch_(char *);
    float safmin;
    extern  int xerbla_(char *, int *);
    extern int isamax_(int *, float *, int *);
    float lstres;
    extern  int spttrs_(int *, int *, float *, float *, 
	    float *, int *, int *);


/*  -- LAPACK routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SPTRFS improves the computed solution to a system of linear */
/*  equations when the coefficient matrix is symmetric positive definite */
/*  and tridiagonal, and provides error bounds and backward error */
/*  estimates for the solution. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0. */

/*  NRHS    (input) INTEGER */
/*          The number of right hand sides, i.e., the number of columns */
/*          of the matrix B.  NRHS >= 0. */

/*  D       (input) REAL array, dimension (N) */
/*          The n diagonal elements of the tridiagonal matrix A. */

/*  E       (input) REAL array, dimension (N-1) */
/*          The (n-1) subdiagonal elements of the tridiagonal matrix A. */

/*  DF      (input) REAL array, dimension (N) */
/*          The n diagonal elements of the diagonal matrix D from the */
/*          factorization computed by SPTTRF. */

/*  EF      (input) REAL array, dimension (N-1) */
/*          The (n-1) subdiagonal elements of the unit bidiagonal factor */
/*          L from the factorization computed by SPTTRF. */

/*  B       (input) REAL array, dimension (LDB,NRHS) */
/*          The right hand side matrix B. */

/*  LDB     (input) INTEGER */
/*          The leading dimension of the array B.  LDB >= MAX(1,N). */

/*  X       (input/output) REAL array, dimension (LDX,NRHS) */
/*          On entry, the solution matrix X, as computed by SPTTRS. */
/*          On exit, the improved solution matrix X. */

/*  LDX     (input) INTEGER */
/*          The leading dimension of the array X.  LDX >= MAX(1,N). */

/*  FERR    (output) REAL array, dimension (NRHS) */
/*          The forward error bound for each solution vector */
/*          X(j) (the j-th column of the solution matrix X). */
/*          If XTRUE is the true solution corresponding to X(j), FERR(j) */
/*          is an estimated upper bound for the magnitude of the largest */
/*          element in (X(j) - XTRUE) divided by the magnitude of the */
/*          largest element in X(j). */

/*  BERR    (output) REAL array, dimension (NRHS) */
/*          The componentwise relative backward error of each solution */
/*          vector X(j) (i.e., the smallest relative change in */
/*          any element of A or B that makes X(j) an exact solution). */

/*  WORK    (workspace) REAL array, dimension (2*N) */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */

/*  Internal Parameters */
/*  =================== */

/*  ITMAX is the maximum number of steps of iterative refinement. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --d__;
    --e;
    --df;
    --ef;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --ferr;
    --berr;
    --work;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*nrhs < 0) {
	*info = -2;
    } else if (*ldb < MAX(1,*n)) {
	*info = -8;
    } else if (*ldx < MAX(1,*n)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SPTRFS", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	i__1 = *nrhs;
	for (j = 1; j <= i__1; ++j) {
	    ferr[j] = 0.f;
	    berr[j] = 0.f;
/* L10: */
	}
	return 0;
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

    nz = 4;
    eps = slamch_("Epsilon");
    safmin = slamch_("Safe minimum");
    safe1 = nz * safmin;
    safe2 = safe1 / eps;

/*     Do for each right hand side */

    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {

	count = 1;
	lstres = 3.f;
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - A * X.  Also compute */
/*        ABS(A)*ABS(x) + ABS(b) for use in the backward error bound. */

	if (*n == 1) {
	    bi = b[j * b_dim1 + 1];
	    dx = d__[1] * x[j * x_dim1 + 1];
	    work[*n + 1] = bi - dx;
	    work[1] = ABS(bi) + ABS(dx);
	} else {
	    bi = b[j * b_dim1 + 1];
	    dx = d__[1] * x[j * x_dim1 + 1];
	    ex = e[1] * x[j * x_dim1 + 2];
	    work[*n + 1] = bi - dx - ex;
	    work[1] = ABS(bi) + ABS(dx) + ABS(ex);
	    i__2 = *n - 1;
	    for (i__ = 2; i__ <= i__2; ++i__) {
		bi = b[i__ + j * b_dim1];
		cx = e[i__ - 1] * x[i__ - 1 + j * x_dim1];
		dx = d__[i__] * x[i__ + j * x_dim1];
		ex = e[i__] * x[i__ + 1 + j * x_dim1];
		work[*n + i__] = bi - cx - dx - ex;
		work[i__] = ABS(bi) + ABS(cx) + ABS(dx) + ABS(ex);
/* L30: */
	    }
	    bi = b[*n + j * b_dim1];
	    cx = e[*n - 1] * x[*n - 1 + j * x_dim1];
	    dx = d__[*n] * x[*n + j * x_dim1];
	    work[*n + *n] = bi - cx - dx;
	    work[*n] = ABS(bi) + ABS(cx) + ABS(dx);
	}

/*        Compute componentwise relative backward error from formula */

/*        MAX(i) ( ABS(R(i)) / ( ABS(A)*ABS(X) + ABS(B) )(i) ) */

/*        where ABS(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

	s = 0.f;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (work[i__] > safe2) {
/* Computing MAX */
		r__2 = s, r__3 = (r__1 = work[*n + i__], ABS(r__1)) / work[
			i__];
		s = MAX(r__2,r__3);
	    } else {
/* Computing MAX */
		r__2 = s, r__3 = ((r__1 = work[*n + i__], ABS(r__1)) + safe1)
			 / (work[i__] + safe1);
		s = MAX(r__2,r__3);
	    }
/* L40: */
	}
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

	if (berr[j] > eps && berr[j] * 2.f <= lstres && count <= 5) {

/*           Update solution and try again. */

	    spttrs_(n, &c__1, &df[1], &ef[1], &work[*n + 1], n, info);
	    saxpy_(n, &c_b11, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1)
		    ;
	    lstres = berr[j];
	    ++count;
	    goto L20;
	}

/*        Bound error from formula */

/*        norm(X - XTRUE) / norm(X) .le. FERR = */
/*        norm( ABS(inv(A))* */
/*           ( ABS(R) + NZ*EPS*( ABS(A)*ABS(X)+ABS(B) ))) / norm(X) */

/*        where */
/*          norm(Z) is the magnitude of the largest component of Z */
/*          inv(A) is the inverse of A */
/*          ABS(Z) is the componentwise absolute value of the matrix or */
/*             vector Z */
/*          NZ is the maximum number of nonzeros in any row of A, plus 1 */
/*          EPS is machine epsilon */

/*        The i-th component of ABS(R)+NZ*EPS*(ABS(A)*ABS(X)+ABS(B)) */
/*        is incremented by SAFE1 if the i-th component of */
/*        ABS(A)*ABS(X) + ABS(B) is less than SAFE2. */

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (work[i__] > safe2) {
		work[i__] = (r__1 = work[*n + i__], ABS(r__1)) + nz * eps * 
			work[i__];
	    } else {
		work[i__] = (r__1 = work[*n + i__], ABS(r__1)) + nz * eps * 
			work[i__] + safe1;
	    }
/* L50: */
	}
	ix = isamax_(n, &work[1], &c__1);
	ferr[j] = work[ix];

/*        Estimate the norm of inv(A). */

/*        Solve M(A) * x = e, where M(A) = (m(i,j)) is given by */

/*           m(i,j) =  ABS(A(i,j)), i = j, */
/*           m(i,j) = -ABS(A(i,j)), i .ne. j, */

/*        and e = [ 1, 1, ..., 1 ]'.  Note M(A) = M(L)*D*M(L)'. */

/*        Solve M(L) * x = e. */

	work[1] = 1.f;
	i__2 = *n;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    work[i__] = work[i__ - 1] * (r__1 = ef[i__ - 1], ABS(r__1)) + 
		    1.f;
/* L60: */
	}

/*        Solve D * M(L)' * x = b. */

	work[*n] /= df[*n];
	for (i__ = *n - 1; i__ >= 1; --i__) {
	    work[i__] = work[i__] / df[i__] + work[i__ + 1] * (r__1 = ef[i__],
		     ABS(r__1));
/* L70: */
	}

/*        Compute norm(inv(A)) = MAX(x(i)), 1<=i<=n. */

	ix = isamax_(n, &work[1], &c__1);
	ferr[j] *= (r__1 = work[ix], ABS(r__1));

/*        Normalize error. */

	lstres = 0.f;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__2 = lstres, r__3 = (r__1 = x[i__ + j * x_dim1], ABS(r__1));
	    lstres = MAX(r__2,r__3);
/* L80: */
	}
	if (lstres != 0.f) {
	    ferr[j] /= lstres;
	}

/* L90: */
    }

    return 0;

/*     End of SPTRFS */

} /* sptrfs_ */
