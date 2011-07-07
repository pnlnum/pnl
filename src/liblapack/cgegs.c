/* cgegs.f -- translated by f2c (version 20061008).
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

static complex c_b1 = {0.f,0.f};
static complex c_b2 = {1.f,0.f};
static int c__1 = 1;
static int c_n1 = -1;

 int cgegs_(char *jobvsl, char *jobvsr, int *n, complex *
	a, int *lda, complex *b, int *ldb, complex *alpha, complex *
	beta, complex *vsl, int *ldvsl, complex *vsr, int *ldvsr, 
	complex *work, int *lwork, float *rwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, vsl_dim1, vsl_offset, 
	    vsr_dim1, vsr_offset, i__1, i__2, i__3;

    /* Local variables */
    int nb, nb1, nb2, nb3, ihi, ilo;
    float eps, anrm, bnrm;
    int itau, lopt;
    extern int lsame_(char *, char *);
    int ileft, iinfo, icols;
    int ilvsl;
    int iwork;
    int ilvsr;
    int irows;
    extern  int cggbak_(char *, char *, int *, int *, 
	    int *, float *, float *, int *, complex *, int *, 
	    int *), cggbal_(char *, int *, complex *, 
	    int *, complex *, int *, int *, int *, float *, 
	    float *, float *, int *);
    extern double clange_(char *, int *, int *, complex *, 
	    int *, float *);
    extern  int cgghrd_(char *, char *, int *, int *, 
	    int *, complex *, int *, complex *, int *, complex *, 
	    int *, complex *, int *, int *), 
	    clascl_(char *, int *, int *, float *, float *, int *, 
	    int *, complex *, int *, int *);
    int ilascl, ilbscl;
    extern  int cgeqrf_(int *, int *, complex *, 
	    int *, complex *, complex *, int *, int *);
    extern double slamch_(char *);
    extern  int clacpy_(char *, int *, int *, complex 
	    *, int *, complex *, int *), claset_(char *, 
	    int *, int *, complex *, complex *, complex *, int *);
    float safmin;
    extern  int xerbla_(char *, int *);
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *);
    float bignum;
    extern  int chgeqz_(char *, char *, char *, int *, 
	    int *, int *, complex *, int *, complex *, int *, 
	    complex *, complex *, complex *, int *, complex *, int *, 
	    complex *, int *, float *, int *);
    int ijobvl, iright, ijobvr;
    float anrmto;
    int lwkmin;
    float bnrmto;
    extern  int cungqr_(int *, int *, int *, 
	    complex *, int *, complex *, complex *, int *, int *),
	     cunmqr_(char *, char *, int *, int *, int *, complex 
	    *, int *, complex *, complex *, int *, complex *, int 
	    *, int *);
    float smlnum;
    int irwork, lwkopt;
    int lquery;


/*  -- LAPACK driver routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  This routine is deprecated and has been replaced by routine CGGES. */

/*  CGEGS computes the eigenvalues, Schur form, and, optionally, the */
/*  left and or/right Schur vectors of a complex matrix pair (A,B). */
/*  Given two square matrices A and B, the generalized Schur */
/*  factorization has the form */

/*     A = Q*S*Z**H,  B = Q*T*Z**H */

/*  where Q and Z are unitary matrices and S and T are upper triangular. */
/*  The columns of Q are the left Schur vectors */
/*  and the columns of Z are the right Schur vectors. */

/*  If only the eigenvalues of (A,B) are needed, the driver routine */
/*  CGEGV should be used instead.  See CGEGV for a description of the */
/*  eigenvalues of the generalized nonsymmetric eigenvalue problem */
/*  (GNEP). */

/*  Arguments */
/*  ========= */

/*  JOBVSL   (input) CHARACTER*1 */
/*          = 'N':  do not compute the left Schur vectors; */
/*          = 'V':  compute the left Schur vectors (returned in VSL). */

/*  JOBVSR   (input) CHARACTER*1 */
/*          = 'N':  do not compute the right Schur vectors; */
/*          = 'V':  compute the right Schur vectors (returned in VSR). */

/*  N       (input) INTEGER */
/*          The order of the matrices A, B, VSL, and VSR.  N >= 0. */

/*  A       (input/output) COMPLEX array, dimension (LDA, N) */
/*          On entry, the matrix A. */
/*          On exit, the upper triangular matrix S from the generalized */
/*          Schur factorization. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of A.  LDA >= MAX(1,N). */

/*  B       (input/output) COMPLEX array, dimension (LDB, N) */
/*          On entry, the matrix B. */
/*          On exit, the upper triangular matrix T from the generalized */
/*          Schur factorization. */

/*  LDB     (input) INTEGER */
/*          The leading dimension of B.  LDB >= MAX(1,N). */

/*  ALPHA   (output) COMPLEX array, dimension (N) */
/*          The complex scalars alpha that define the eigenvalues of */
/*          GNEP.  ALPHA(j) = S(j,j), the diagonal element of the Schur */
/*          form of A. */

/*  BETA    (output) COMPLEX array, dimension (N) */
/*          The non-negative float scalars beta that define the */
/*          eigenvalues of GNEP.  BETA(j) = T(j,j), the diagonal element */
/*          of the triangular factor T. */

/*          Together, the quantities alpha = ALPHA(j) and beta = BETA(j) */
/*          represent the j-th eigenvalue of the matrix pair (A,B), in */
/*          one of the forms lambda = alpha/beta or mu = beta/alpha. */
/*          Since either lambda or mu may overflow, they should not, */
/*          in general, be computed. */

/*  VSL     (output) COMPLEX array, dimension (LDVSL,N) */
/*          If JOBVSL = 'V', the matrix of left Schur vectors Q. */
/*          Not referenced if JOBVSL = 'N'. */

/*  LDVSL   (input) INTEGER */
/*          The leading dimension of the matrix VSL. LDVSL >= 1, and */
/*          if JOBVSL = 'V', LDVSL >= N. */

/*  VSR     (output) COMPLEX array, dimension (LDVSR,N) */
/*          If JOBVSR = 'V', the matrix of right Schur vectors Z. */
/*          Not referenced if JOBVSR = 'N'. */

/*  LDVSR   (input) INTEGER */
/*          The leading dimension of the matrix VSR. LDVSR >= 1, and */
/*          if JOBVSR = 'V', LDVSR >= N. */

/*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK)) */
/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

/*  LWORK   (input) INTEGER */
/*          The dimension of the array WORK.  LWORK >= MAX(1,2*N). */
/*          For good performance, LWORK must generally be larger. */
/*          To compute the optimal value of LWORK, call ILAENV to get */
/*          blocksizes (for CGEQRF, CUNMQR, and CUNGQR.)  Then compute: */
/*          NB  -- MAX of the blocksizes for CGEQRF, CUNMQR, and CUNGQR; */
/*          the optimal LWORK is N*(NB+1). */

/*          If LWORK = -1, then a workspace query is assumed; the routine */
/*          only calculates the optimal size of the WORK array, returns */
/*          this value as the first entry of the WORK array, and no error */
/*          message related to LWORK is issued by XERBLA. */

/*  RWORK   (workspace) REAL array, dimension (3*N) */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/*          =1,...,N: */
/*                The QZ iteration failed.  (A,B) are not in Schur */
/*                form, but ALPHA(j) and BETA(j) should be correct for */
/*                j=INFO+1,...,N. */
/*          > N:  errors that usually indicate LAPACK problems: */
/*                =N+1: error return from CGGBAL */
/*                =N+2: error return from CGEQRF */
/*                =N+3: error return from CUNMQR */
/*                =N+4: error return from CUNGQR */
/*                =N+5: error return from CGGHRD */
/*                =N+6: error return from CHGEQZ (other than failed */
/*                                               iteration) */
/*                =N+7: error return from CGGBAK (computing VSL) */
/*                =N+8: error return from CGGBAK (computing VSR) */
/*                =N+9: error return from CLASCL (various places) */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --alpha;
    --beta;
    vsl_dim1 = *ldvsl;
    vsl_offset = 1 + vsl_dim1;
    vsl -= vsl_offset;
    vsr_dim1 = *ldvsr;
    vsr_offset = 1 + vsr_dim1;
    vsr -= vsr_offset;
    --work;
    --rwork;

    /* Function Body */
    if (lsame_(jobvsl, "N")) {
	ijobvl = 1;
	ilvsl = FALSE;
    } else if (lsame_(jobvsl, "V")) {
	ijobvl = 2;
	ilvsl = TRUE;
    } else {
	ijobvl = -1;
	ilvsl = FALSE;
    }

    if (lsame_(jobvsr, "N")) {
	ijobvr = 1;
	ilvsr = FALSE;
    } else if (lsame_(jobvsr, "V")) {
	ijobvr = 2;
	ilvsr = TRUE;
    } else {
	ijobvr = -1;
	ilvsr = FALSE;
    }

/*     Test the input arguments */

/* Computing MAX */
    i__1 = *n << 1;
    lwkmin = MAX(i__1,1);
    lwkopt = lwkmin;
    work[1].r = (float) lwkopt, work[1].i = 0.f;
    lquery = *lwork == -1;
    *info = 0;
    if (ijobvl <= 0) {
	*info = -1;
    } else if (ijobvr <= 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < MAX(1,*n)) {
	*info = -5;
    } else if (*ldb < MAX(1,*n)) {
	*info = -7;
    } else if (*ldvsl < 1 || ilvsl && *ldvsl < *n) {
	*info = -11;
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
	*info = -13;
    } else if (*lwork < lwkmin && ! lquery) {
	*info = -15;
    }

    if (*info == 0) {
	nb1 = ilaenv_(&c__1, "CGEQRF", " ", n, n, &c_n1, &c_n1);
	nb2 = ilaenv_(&c__1, "CUNMQR", " ", n, n, n, &c_n1);
	nb3 = ilaenv_(&c__1, "CUNGQR", " ", n, n, n, &c_n1);
/* Computing MAX */
	i__1 = MAX(nb1,nb2);
	nb = MAX(i__1,nb3);
	lopt = *n * (nb + 1);
	work[1].r = (float) lopt, work[1].i = 0.f;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CGEGS ", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Get machine constants */

    eps = slamch_("E") * slamch_("B");
    safmin = slamch_("S");
    smlnum = *n * safmin / eps;
    bignum = 1.f / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

    anrm = clange_("M", n, n, &a[a_offset], lda, &rwork[1]);
    ilascl = FALSE;
    if (anrm > 0.f && anrm < smlnum) {
	anrmto = smlnum;
	ilascl = TRUE;
    } else if (anrm > bignum) {
	anrmto = bignum;
	ilascl = TRUE;
    }

    if (ilascl) {
	clascl_("G", &c_n1, &c_n1, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

    bnrm = clange_("M", n, n, &b[b_offset], ldb, &rwork[1]);
    ilbscl = FALSE;
    if (bnrm > 0.f && bnrm < smlnum) {
	bnrmto = smlnum;
	ilbscl = TRUE;
    } else if (bnrm > bignum) {
	bnrmto = bignum;
	ilbscl = TRUE;
    }

    if (ilbscl) {
	clascl_("G", &c_n1, &c_n1, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
    }

/*     Permute the matrix to make it more nearly triangular */

    ileft = 1;
    iright = *n + 1;
    irwork = iright + *n;
    iwork = 1;
    cggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[
	    ileft], &rwork[iright], &rwork[irwork], &iinfo);
    if (iinfo != 0) {
	*info = *n + 1;
	goto L10;
    }

/*     Reduce B to triangular form, and initialize VSL and/or VSR */

    irows = ihi + 1 - ilo;
    icols = *n + 1 - ilo;
    itau = iwork;
    iwork = itau + irows;
    i__1 = *lwork + 1 - iwork;
    cgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwork], &i__1, &iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__3 = iwork;
	i__1 = lwkopt, i__2 = (int) work[i__3].r + iwork - 1;
	lwkopt = MAX(i__1,i__2);
    }
    if (iinfo != 0) {
	*info = *n + 2;
	goto L10;
    }

    i__1 = *lwork + 1 - iwork;
    cunmqr_("L", "C", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwork], &i__1, &
	    iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__3 = iwork;
	i__1 = lwkopt, i__2 = (int) work[i__3].r + iwork - 1;
	lwkopt = MAX(i__1,i__2);
    }
    if (iinfo != 0) {
	*info = *n + 3;
	goto L10;
    }

    if (ilvsl) {
	claset_("Full", n, n, &c_b1, &c_b2, &vsl[vsl_offset], ldvsl);
	i__1 = irows - 1;
	i__2 = irows - 1;
	clacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[ilo 
		+ 1 + ilo * vsl_dim1], ldvsl);
	i__1 = *lwork + 1 - iwork;
	cungqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &
		work[itau], &work[iwork], &i__1, &iinfo);
	if (iinfo >= 0) {
/* Computing MAX */
	    i__3 = iwork;
	    i__1 = lwkopt, i__2 = (int) work[i__3].r + iwork - 1;
	    lwkopt = MAX(i__1,i__2);
	}
	if (iinfo != 0) {
	    *info = *n + 4;
	    goto L10;
	}
    }

    if (ilvsr) {
	claset_("Full", n, n, &c_b1, &c_b2, &vsr[vsr_offset], ldvsr);
    }

/*     Reduce to generalized Hessenberg form */

    cgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
	    ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &iinfo);
    if (iinfo != 0) {
	*info = *n + 5;
	goto L10;
    }

/*     Perform QZ algorithm, computing Schur vectors if desired */

    iwork = itau;
    i__1 = *lwork + 1 - iwork;
    chgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, &
	    vsr[vsr_offset], ldvsr, &work[iwork], &i__1, &rwork[irwork], &
	    iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__3 = iwork;
	i__1 = lwkopt, i__2 = (int) work[i__3].r + iwork - 1;
	lwkopt = MAX(i__1,i__2);
    }
    if (iinfo != 0) {
	if (iinfo > 0 && iinfo <= *n) {
	    *info = iinfo;
	} else if (iinfo > *n && iinfo <= *n << 1) {
	    *info = iinfo - *n;
	} else {
	    *info = *n + 6;
	}
	goto L10;
    }

/*     Apply permutation to VSL and VSR */

    if (ilvsl) {
	cggbak_("P", "L", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsl[vsl_offset], ldvsl, &iinfo);
	if (iinfo != 0) {
	    *info = *n + 7;
	    goto L10;
	}
    }
    if (ilvsr) {
	cggbak_("P", "R", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsr[vsr_offset], ldvsr, &iinfo);
	if (iinfo != 0) {
	    *info = *n + 8;
	    goto L10;
	}
    }

/*     Undo scaling */

    if (ilascl) {
	clascl_("U", &c_n1, &c_n1, &anrmto, &anrm, n, n, &a[a_offset], lda, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
	clascl_("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &alpha[1], n, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
    }

    if (ilbscl) {
	clascl_("U", &c_n1, &c_n1, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
	clascl_("G", &c_n1, &c_n1, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
    }

L10:
    work[1].r = (float) lwkopt, work[1].i = 0.f;

    return 0;

/*     End of CGEGS */

} /* cgegs_ */
