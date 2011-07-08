/* cla_gbrcond_x.f -- translated by f2c (version 20061008).
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

double cla_gbrcond_x__(char *trans, int *n, int *kl, int *ku, 
	complex *ab, int *ldab, complex *afb, int *ldafb, int *
	ipiv, complex *x, int *info, complex *work, float *rwork, unsigned long 
	trans_len)
{
    /* System generated locals */
    int ab_dim1, ab_offset, afb_dim1, afb_offset, i__1, i__2, i__3, i__4;
    float ret_val, r__1, r__2;
    complex q__1, q__2;

    /* Builtin functions */
    double r_imag(complex *);
    void c_div(complex *, complex *, complex *);

    /* Local variables */
    int i__, j, kd, ke;
    float tmp;
    int kase;
    extern int lsame_(char *, char *);
    int isave[3];
    float anorm;
    extern  int clacn2_(int *, complex *, complex *, float 
	    *, int *, int *), xerbla_(char *, int *), 
	    cgbtrs_(char *, int *, int *, int *, int *, 
	    complex *, int *, int *, complex *, int *, int *);
    float ainvnm;
    int notrans;


/*     -- LAPACK routine (version 3.2.1)                               -- */
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

/*     CLA_GBRCOND_X Computes the infinity norm condition number of */
/*     op(A) * diag(X) where X is a COMPLEX vector. */

/*  Arguments */
/*  ========= */

/*     TRANS   (input) CHARACTER*1 */
/*     Specifies the form of the system of equations: */
/*       = 'N':  A * X = B     (No transpose) */
/*       = 'T':  A**T * X = B  (Transpose) */
/*       = 'C':  A**H * X = B  (Conjugate Transpose = Transpose) */

/*     N       (input) INTEGER */
/*     The number of linear equations, i.e., the order of the */
/*     matrix A.  N >= 0. */

/*     KL      (input) INTEGER */
/*     The number of subdiagonals within the band of A.  KL >= 0. */

/*     KU      (input) INTEGER */
/*     The number of superdiagonals within the band of A.  KU >= 0. */

/*     AB      (input) COMPLEX array, dimension (LDAB,N) */
/*     On entry, the matrix A in band storage, in rows 1 to KL+KU+1. */
/*     The j-th column of A is stored in the j-th column of the */
/*     array AB as follows: */
/*     AB(KU+1+i-j,j) = A(i,j) for MAX(1,j-KU)<=i<=MIN(N,j+kl) */

/*     LDAB    (input) INTEGER */
/*     The leading dimension of the array AB.  LDAB >= KL+KU+1. */

/*     AFB     (input) COMPLEX array, dimension (LDAFB,N) */
/*     Details of the LU factorization of the band matrix A, as */
/*     computed by CGBTRF.  U is stored as an upper triangular */
/*     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, */
/*     and the multipliers used during the factorization are stored */
/*     in rows KL+KU+2 to 2*KL+KU+1. */

/*     LDAFB   (input) INTEGER */
/*     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1. */

/*     IPIV    (input) INTEGER array, dimension (N) */
/*     The pivot indices from the factorization A = P*L*U */
/*     as computed by CGBTRF; row i of the matrix was interchanged */
/*     with row IPIV(i). */

/*     X       (input) COMPLEX array, dimension (N) */
/*     The vector X in the formula op(A) * diag(X). */

/*     INFO    (output) INTEGER */
/*       = 0:  Successful exit. */
/*     i > 0:  The ith argument is invalid. */

/*     WORK    (input) COMPLEX array, dimension (2*N). */
/*     Workspace. */

/*     RWORK   (input) REAL array, dimension (N). */
/*     Workspace. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function Definitions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    afb_dim1 = *ldafb;
    afb_offset = 1 + afb_dim1;
    afb -= afb_offset;
    --ipiv;
    --x;
    --work;
    --rwork;

    /* Function Body */
    ret_val = 0.f;

    *info = 0;
    notrans = lsame_(trans, "N");
    if (! notrans && ! lsame_(trans, "T") && ! lsame_(
	    trans, "C")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*kl < 0 || *kl > *n - 1) {
	*info = -3;
    } else if (*ku < 0 || *ku > *n - 1) {
	*info = -4;
    } else if (*ldab < *kl + *ku + 1) {
	*info = -6;
    } else if (*ldafb < (*kl << 1) + *ku + 1) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CLA_GBRCOND_X", &i__1);
	return ret_val;
    }

/*     Compute norm of op(A)*op2(C). */

    kd = *ku + 1;
    ke = *kl + 1;
    anorm = 0.f;
    if (notrans) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tmp = 0.f;
/* Computing MAX */
	    i__2 = i__ - *kl;
/* Computing MIN */
	    i__4 = i__ + *ku;
	    i__3 = MIN(i__4,*n);
	    for (j = MAX(i__2,1); j <= i__3; ++j) {
		i__2 = kd + i__ - j + j * ab_dim1;
		i__4 = j;
		q__2.r = ab[i__2].r * x[i__4].r - ab[i__2].i * x[i__4].i, 
			q__2.i = ab[i__2].r * x[i__4].i + ab[i__2].i * x[i__4]
			.r;
		q__1.r = q__2.r, q__1.i = q__2.i;
		tmp += (r__1 = q__1.r, ABS(r__1)) + (r__2 = r_imag(&q__1), 
			ABS(r__2));
	    }
	    rwork[i__] = tmp;
	    anorm = MAX(anorm,tmp);
	}
    } else {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tmp = 0.f;
/* Computing MAX */
	    i__3 = i__ - *kl;
/* Computing MIN */
	    i__4 = i__ + *ku;
	    i__2 = MIN(i__4,*n);
	    for (j = MAX(i__3,1); j <= i__2; ++j) {
		i__3 = ke - i__ + j + i__ * ab_dim1;
		i__4 = j;
		q__2.r = ab[i__3].r * x[i__4].r - ab[i__3].i * x[i__4].i, 
			q__2.i = ab[i__3].r * x[i__4].i + ab[i__3].i * x[i__4]
			.r;
		q__1.r = q__2.r, q__1.i = q__2.i;
		tmp += (r__1 = q__1.r, ABS(r__1)) + (r__2 = r_imag(&q__1), 
			ABS(r__2));
	    }
	    rwork[i__] = tmp;
	    anorm = MAX(anorm,tmp);
	}
    }

/*     Quick return if possible. */

    if (*n == 0) {
	ret_val = 1.f;
	return ret_val;
    } else if (anorm == 0.f) {
	return ret_val;
    }

/*     Estimate the norm of inv(op(A)). */

    ainvnm = 0.f;

    kase = 0;
L10:
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
    if (kase != 0) {
	if (kase == 2) {

/*           Multiply by R. */

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = i__;
		i__3 = i__;
		i__4 = i__;
		q__1.r = rwork[i__4] * work[i__3].r, q__1.i = rwork[i__4] * 
			work[i__3].i;
		work[i__2].r = q__1.r, work[i__2].i = q__1.i;
	    }

	    if (notrans) {
		cgbtrs_("No transpose", n, kl, ku, &c__1, &afb[afb_offset], 
			ldafb, &ipiv[1], &work[1], n, info);
	    } else {
		cgbtrs_("Conjugate transpose", n, kl, ku, &c__1, &afb[
			afb_offset], ldafb, &ipiv[1], &work[1], n, info);
	    }

/*           Multiply by inv(X). */

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = i__;
		c_div(&q__1, &work[i__], &x[i__]);
		work[i__2].r = q__1.r, work[i__2].i = q__1.i;
	    }
	} else {

/*           Multiply by inv(X'). */

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = i__;
		c_div(&q__1, &work[i__], &x[i__]);
		work[i__2].r = q__1.r, work[i__2].i = q__1.i;
	    }

	    if (notrans) {
		cgbtrs_("Conjugate transpose", n, kl, ku, &c__1, &afb[
			afb_offset], ldafb, &ipiv[1], &work[1], n, info);
	    } else {
		cgbtrs_("No transpose", n, kl, ku, &c__1, &afb[afb_offset], 
			ldafb, &ipiv[1], &work[1], n, info);
	    }

/*           Multiply by R. */

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = i__;
		i__3 = i__;
		i__4 = i__;
		q__1.r = rwork[i__4] * work[i__3].r, q__1.i = rwork[i__4] * 
			work[i__3].i;
		work[i__2].r = q__1.r, work[i__2].i = q__1.i;
	    }
	}
	goto L10;
    }

/*     Compute the estimate of the reciprocal condition number. */

    if (ainvnm != 0.f) {
	ret_val = 1.f / ainvnm;
    }

    return ret_val;

} /* cla_gbrcond_x__ */