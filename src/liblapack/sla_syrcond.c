/* sla_syrcond.f -- translated by f2c (version 20061008).
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

double sla_syrcond__(char *uplo, int *n, float *a, int *lda, float *
	af, int *ldaf, int *ipiv, int *cmode, float *c__, int *
	info, float *work, int *iwork, unsigned long uplo_len)
{
    /* System generated locals */
    int a_dim1, a_offset, af_dim1, af_offset, i__1, i__2;
    float ret_val, r__1;

    /* Local variables */
    int i__, j;
    int up;
    float tmp;
    int kase;
    extern int lsame_(char *, char *);
    int isave[3];
    extern  int slacn2_(int *, float *, float *, int *, 
	    float *, int *, int *);
    extern double slamch_(char *);
    extern  int xerbla_(char *, int *);
    float ainvnm;
    char normin[1];
    float smlnum;
    extern  int ssytrs_(char *, int *, int *, float *, 
	    int *, int *, float *, int *, int *);


/*     -- LAPACK routine (version 3.2.1)                                 -- */
/*     -- Contributed by James Demmel, Deaglan Halligan, Yozo Hida and -- */
/*     -- Jason Riedy of Univ. of California Berkeley.                 -- */
/*     -- April 2009                                                   -- */

/*     -- LAPACK is a software package provided by Univ. of Tennessee, -- */
/*     -- Univ. of California Berkeley and NAG Ltd.                    -- */

/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments */
/*     .. */

/*  Purpose */
/*  ======= */

/*     SLA_SYRCOND estimates the Skeel condition number of  op(A) * op2(C) */
/*     where op2 is determined by CMODE as follows */
/*     CMODE =  1    op2(C) = C */
/*     CMODE =  0    op2(C) = I */
/*     CMODE = -1    op2(C) = inv(C) */
/*     The Skeel condition number cond(A) = norminf( |inv(A)||A| ) */
/*     is computed by computing scaling factors R such that */
/*     diag(R)*A*op2(C) is row equilibrated and computing the standard */
/*     infinity-norm condition number. */

/*  Arguments */
/*  ========== */

/*     UPLO    (input) CHARACTER*1 */
/*       = 'U':  Upper triangle of A is stored; */
/*       = 'L':  Lower triangle of A is stored. */

/*     N       (input) INTEGER */
/*     The number of linear equations, i.e., the order of the */
/*     matrix A.  N >= 0. */

/*     A       (input) REAL array, dimension (LDA,N) */
/*     On entry, the N-by-N matrix A. */

/*     LDA     (input) INTEGER */
/*     The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     AF      (input) REAL array, dimension (LDAF,N) */
/*     The block diagonal matrix D and the multipliers used to */
/*     obtain the factor U or L as computed by SSYTRF. */

/*     LDAF    (input) INTEGER */
/*     The leading dimension of the array AF.  LDAF >= MAX(1,N). */

/*     IPIV    (input) INTEGER array, dimension (N) */
/*     Details of the interchanges and the block structure of D */
/*     as determined by SSYTRF. */

/*     CMODE   (input) INTEGER */
/*     Determines op2(C) in the formula op(A) * op2(C) as follows: */
/*     CMODE =  1    op2(C) = C */
/*     CMODE =  0    op2(C) = I */
/*     CMODE = -1    op2(C) = inv(C) */

/*     C       (input) REAL array, dimension (N) */
/*     The vector C in the formula op(A) * op2(C). */

/*     INFO    (output) INTEGER */
/*       = 0:  Successful exit. */
/*     i > 0:  The ith argument is invalid. */

/*     WORK    (input) REAL array, dimension (3*N). */
/*     Workspace. */

/*     IWORK   (input) INTEGER array, dimension (N). */
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
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    af_dim1 = *ldaf;
    af_offset = 1 + af_dim1;
    af -= af_offset;
    --ipiv;
    --c__;
    --work;
    --iwork;

    /* Function Body */
    ret_val = 0.f;

    *info = 0;
    if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SLA_SYRCOND", &i__1);
	return ret_val;
    }
    if (*n == 0) {
	ret_val = 1.f;
	return ret_val;
    }
    up = FALSE;
    if (lsame_(uplo, "U")) {
	up = TRUE;
    }

/*     Compute the equilibration matrix R such that */
/*     inv(R)*A*C has unit 1-norm. */

    if (up) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tmp = 0.f;
	    if (*cmode == 1) {
		i__2 = i__;
		for (j = 1; j <= i__2; ++j) {
		    tmp += (r__1 = a[j + i__ * a_dim1] * c__[j], ABS(r__1));
		}
		i__2 = *n;
		for (j = i__ + 1; j <= i__2; ++j) {
		    tmp += (r__1 = a[i__ + j * a_dim1] * c__[j], ABS(r__1));
		}
	    } else if (*cmode == 0) {
		i__2 = i__;
		for (j = 1; j <= i__2; ++j) {
		    tmp += (r__1 = a[j + i__ * a_dim1], ABS(r__1));
		}
		i__2 = *n;
		for (j = i__ + 1; j <= i__2; ++j) {
		    tmp += (r__1 = a[i__ + j * a_dim1], ABS(r__1));
		}
	    } else {
		i__2 = i__;
		for (j = 1; j <= i__2; ++j) {
		    tmp += (r__1 = a[j + i__ * a_dim1] / c__[j], ABS(r__1));
		}
		i__2 = *n;
		for (j = i__ + 1; j <= i__2; ++j) {
		    tmp += (r__1 = a[i__ + j * a_dim1] / c__[j], ABS(r__1));
		}
	    }
	    work[(*n << 1) + i__] = tmp;
	}
    } else {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tmp = 0.f;
	    if (*cmode == 1) {
		i__2 = i__;
		for (j = 1; j <= i__2; ++j) {
		    tmp += (r__1 = a[i__ + j * a_dim1] * c__[j], ABS(r__1));
		}
		i__2 = *n;
		for (j = i__ + 1; j <= i__2; ++j) {
		    tmp += (r__1 = a[j + i__ * a_dim1] * c__[j], ABS(r__1));
		}
	    } else if (*cmode == 0) {
		i__2 = i__;
		for (j = 1; j <= i__2; ++j) {
		    tmp += (r__1 = a[i__ + j * a_dim1], ABS(r__1));
		}
		i__2 = *n;
		for (j = i__ + 1; j <= i__2; ++j) {
		    tmp += (r__1 = a[j + i__ * a_dim1], ABS(r__1));
		}
	    } else {
		i__2 = i__;
		for (j = 1; j <= i__2; ++j) {
		    tmp += (r__1 = a[i__ + j * a_dim1] / c__[j], ABS(r__1));
		}
		i__2 = *n;
		for (j = i__ + 1; j <= i__2; ++j) {
		    tmp += (r__1 = a[j + i__ * a_dim1] / c__[j], ABS(r__1));
		}
	    }
	    work[(*n << 1) + i__] = tmp;
	}
    }

/*     Estimate the norm of inv(op(A)). */

    smlnum = slamch_("Safe minimum");
    ainvnm = 0.f;
    *(unsigned char *)normin = 'N';
    kase = 0;
L10:
    slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
    if (kase != 0) {
	if (kase == 2) {

/*           Multiply by R. */

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		work[i__] *= work[(*n << 1) + i__];
	    }
	    if (up) {
		ssytrs_("U", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info);
	    } else {
		ssytrs_("L", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info);
	    }

/*           Multiply by inv(C). */

	    if (*cmode == 1) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    work[i__] /= c__[i__];
		}
	    } else if (*cmode == -1) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    work[i__] *= c__[i__];
		}
	    }
	} else {

/*           Multiply by inv(C'). */

	    if (*cmode == 1) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    work[i__] /= c__[i__];
		}
	    } else if (*cmode == -1) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    work[i__] *= c__[i__];
		}
	    }
	    if (up) {
		ssytrs_("U", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info);
	    } else {
		ssytrs_("L", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info);
	    }

/*           Multiply by R. */

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		work[i__] *= work[(*n << 1) + i__];
	    }
	}

	goto L10;
    }

/*     Compute the estimate of the reciprocal condition number. */

    if (ainvnm != 0.f) {
	ret_val = 1.f / ainvnm;
    }

    return ret_val;

} /* sla_syrcond__ */
