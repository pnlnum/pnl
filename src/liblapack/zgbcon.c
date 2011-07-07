/* zgbcon.f -- translated by f2c (version 20061008).
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

 int zgbcon_(char *norm, int *n, int *kl, int *ku, 
	 doublecomplex *ab, int *ldab, int *ipiv, double *anorm, 
	double *rcond, doublecomplex *work, double *rwork, int *
	info)
{
    /* System generated locals */
    int ab_dim1, ab_offset, i__1, i__2, i__3;
    double d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    int j;
    doublecomplex t;
    int kd, lm, jp, ix, kase, kase1;
    double scale;
    extern int lsame_(char *, char *);
    int isave[3];
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, int *, 
	    doublecomplex *, int *, doublecomplex *, int *);
    int lnoti;
    extern  int zaxpy_(int *, doublecomplex *, 
	    doublecomplex *, int *, doublecomplex *, int *), zlacn2_(
	    int *, doublecomplex *, doublecomplex *, double *, 
	    int *, int *);
    extern double dlamch_(char *);
    extern  int xerbla_(char *, int *);
    double ainvnm;
    extern int izamax_(int *, doublecomplex *, int *);
    int onenrm;
    extern  int zlatbs_(char *, char *, char *, char *, 
	    int *, int *, doublecomplex *, int *, doublecomplex *, 
	     double *, double *, int *), zdrscl_(int *, double *, doublecomplex *, 
	    int *);
    char normin[1];
    double smlnum;


/*  -- LAPACK routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     Modified to call ZLACN2 in place of ZLACON, 10 Feb 03, SJH. */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  ZGBCON estimates the reciprocal of the condition number of a complex */
/*  general band matrix A, in either the 1-norm or the infinity-norm, */
/*  using the LU factorization computed by ZGBTRF. */

/*  An estimate is obtained for norm(inv(A)), and the reciprocal of the */
/*  condition number is computed as */
/*     RCOND = 1 / ( norm(A) * norm(inv(A)) ). */

/*  Arguments */
/*  ========= */

/*  NORM    (input) CHARACTER*1 */
/*          Specifies whether the 1-norm condition number or the */
/*          infinity-norm condition number is required: */
/*          = '1' or 'O':  1-norm; */
/*          = 'I':         Infinity-norm. */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0. */

/*  KL      (input) INTEGER */
/*          The number of subdiagonals within the band of A.  KL >= 0. */

/*  KU      (input) INTEGER */
/*          The number of superdiagonals within the band of A.  KU >= 0. */

/*  AB      (input) COMPLEX*16 array, dimension (LDAB,N) */
/*          Details of the LU factorization of the band matrix A, as */
/*          computed by ZGBTRF.  U is stored as an upper triangular band */
/*          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and */
/*          the multipliers used during the factorization are stored in */
/*          rows KL+KU+2 to 2*KL+KU+1. */

/*  LDAB    (input) INTEGER */
/*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1. */

/*  IPIV    (input) INTEGER array, dimension (N) */
/*          The pivot indices; for 1 <= i <= N, row i of the matrix was */
/*          interchanged with row IPIV(i). */

/*  ANORM   (input) DOUBLE PRECISION */
/*          If NORM = '1' or 'O', the 1-norm of the original matrix A. */
/*          If NORM = 'I', the infinity-norm of the original matrix A. */

/*  RCOND   (output) DOUBLE PRECISION */
/*          The reciprocal of the condition number of the matrix A, */
/*          computed as RCOND = 1/(norm(A) * norm(inv(A))). */

/*  WORK    (workspace) COMPLEX*16 array, dimension (2*N) */

/*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N) */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0: if INFO = -i, the i-th argument had an illegal value */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
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
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --ipiv;
    --work;
    --rwork;

    /* Function Body */
    *info = 0;
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O");
    if (! onenrm && ! lsame_(norm, "I")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*kl < 0) {
	*info = -3;
    } else if (*ku < 0) {
	*info = -4;
    } else if (*ldab < (*kl << 1) + *ku + 1) {
	*info = -6;
    } else if (*anorm < 0.) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZGBCON", &i__1);
	return 0;
    }

/*     Quick return if possible */

    *rcond = 0.;
    if (*n == 0) {
	*rcond = 1.;
	return 0;
    } else if (*anorm == 0.) {
	return 0;
    }

    smlnum = dlamch_("Safe minimum");

/*     Estimate the norm of inv(A). */

    ainvnm = 0.;
    *(unsigned char *)normin = 'N';
    if (onenrm) {
	kase1 = 1;
    } else {
	kase1 = 2;
    }
    kd = *kl + *ku + 1;
    lnoti = *kl > 0;
    kase = 0;
L10:
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
    if (kase != 0) {
	if (kase == kase1) {

/*           Multiply by inv(L). */

	    if (lnoti) {
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
		    i__2 = *kl, i__3 = *n - j;
		    lm = MIN(i__2,i__3);
		    jp = ipiv[j];
		    i__2 = jp;
		    t.r = work[i__2].r, t.i = work[i__2].i;
		    if (jp != j) {
			i__2 = jp;
			i__3 = j;
			work[i__2].r = work[i__3].r, work[i__2].i = work[i__3]
				.i;
			i__2 = j;
			work[i__2].r = t.r, work[i__2].i = t.i;
		    }
		    z__1.r = -t.r, z__1.i = -t.i;
		    zaxpy_(&lm, &z__1, &ab[kd + 1 + j * ab_dim1], &c__1, &
			    work[j + 1], &c__1);
/* L20: */
		}
	    }

/*           Multiply by inv(U). */

	    i__1 = *kl + *ku;
	    zlatbs_("Upper", "No transpose", "Non-unit", normin, n, &i__1, &
		    ab[ab_offset], ldab, &work[1], &scale, &rwork[1], info);
	} else {

/*           Multiply by inv(U'). */

	    i__1 = *kl + *ku;
	    zlatbs_("Upper", "Conjugate transpose", "Non-unit", normin, n, &
		    i__1, &ab[ab_offset], ldab, &work[1], &scale, &rwork[1], 
		    info);

/*           Multiply by inv(L'). */

	    if (lnoti) {
		for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
		    i__1 = *kl, i__2 = *n - j;
		    lm = MIN(i__1,i__2);
		    i__1 = j;
		    i__2 = j;
		    zdotc_(&z__2, &lm, &ab[kd + 1 + j * ab_dim1], &c__1, &
			    work[j + 1], &c__1);
		    z__1.r = work[i__2].r - z__2.r, z__1.i = work[i__2].i - 
			    z__2.i;
		    work[i__1].r = z__1.r, work[i__1].i = z__1.i;
		    jp = ipiv[j];
		    if (jp != j) {
			i__1 = jp;
			t.r = work[i__1].r, t.i = work[i__1].i;
			i__1 = jp;
			i__2 = j;
			work[i__1].r = work[i__2].r, work[i__1].i = work[i__2]
				.i;
			i__1 = j;
			work[i__1].r = t.r, work[i__1].i = t.i;
		    }
/* L30: */
		}
	    }
	}

/*        Divide X by 1/SCALE if doing so will not cause overflow. */

	*(unsigned char *)normin = 'Y';
	if (scale != 1.) {
	    ix = izamax_(n, &work[1], &c__1);
	    i__1 = ix;
	    if (scale < ((d__1 = work[i__1].r, ABS(d__1)) + (d__2 = d_imag(&
		    work[ix]), ABS(d__2))) * smlnum || scale == 0.) {
		goto L40;
	    }
	    zdrscl_(n, &scale, &work[1], &c__1);
	}
	goto L10;
    }

/*     Compute the estimate of the reciprocal condition number. */

    if (ainvnm != 0.) {
	*rcond = 1. / ainvnm / *anorm;
    }

L40:
    return 0;

/*     End of ZGBCON */

} /* zgbcon_ */
