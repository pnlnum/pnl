/* dgegs.f -- translated by f2c (version 20061008).
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
static int c_n1 = -1;
static double c_b36 = 0.;
static double c_b37 = 1.;

 int dgegs_(char *jobvsl, char *jobvsr, int *n, 
	double *a, int *lda, double *b, int *ldb, double *
	alphar, double *alphai, double *beta, double *vsl, 
	int *ldvsl, double *vsr, int *ldvsr, double *work, 
	int *lwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, vsl_dim1, vsl_offset, 
	    vsr_dim1, vsr_offset, i__1, i__2;

    /* Local variables */
    int nb, nb1, nb2, nb3, ihi, ilo;
    double eps, anrm, bnrm;
    int itau, lopt;
    extern int lsame_(char *, char *);
    int ileft, iinfo, icols;
    int ilvsl;
    int iwork;
    int ilvsr;
    int irows;
    extern  int dggbak_(char *, char *, int *, int *, 
	    int *, double *, double *, int *, double *, 
	    int *, int *), dggbal_(char *, int *, 
	    double *, int *, double *, int *, int *, 
	    int *, double *, double *, double *, int *);
    extern double dlamch_(char *), dlange_(char *, int *, 
	    int *, double *, int *, double *);
    extern  int dgghrd_(char *, char *, int *, int *, 
	    int *, double *, int *, double *, int *, 
	    double *, int *, double *, int *, int *), dlascl_(char *, int *, int *, double 
	    *, double *, int *, int *, double *, int *, 
	    int *);
    int ilascl, ilbscl;
    extern  int dgeqrf_(int *, int *, double *, 
	    int *, double *, double *, int *, int *), 
	    dlacpy_(char *, int *, int *, double *, int *, 
	    double *, int *);
    double safmin;
    extern  int dlaset_(char *, int *, int *, 
	    double *, double *, double *, int *), 
	    xerbla_(char *, int *);
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *);
    double bignum;
    extern  int dhgeqz_(char *, char *, char *, int *, 
	    int *, int *, double *, int *, double *, 
	    int *, double *, double *, double *, double *, 
	     int *, double *, int *, double *, int *, 
	    int *);
    int ijobvl, iright, ijobvr;
    extern  int dorgqr_(int *, int *, int *, 
	    double *, int *, double *, double *, int *, 
	    int *);
    double anrmto;
    int lwkmin;
    double bnrmto;
    extern  int dormqr_(char *, char *, int *, int *, 
	    int *, double *, int *, double *, double *, 
	    int *, double *, int *, int *);
    double smlnum;
    int lwkopt;
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

/*  This routine is deprecated and has been replaced by routine DGGES. */

/*  DGEGS computes the eigenvalues, float Schur form, and, optionally, */
/*  left and or/right Schur vectors of a float matrix pair (A,B). */
/*  Given two square matrices A and B, the generalized float Schur */
/*  factorization has the form */

/*    A = Q*S*Z**T,  B = Q*T*Z**T */

/*  where Q and Z are orthogonal matrices, T is upper triangular, and S */
/*  is an upper quasi-triangular matrix with 1-by-1 and 2-by-2 diagonal */
/*  blocks, the 2-by-2 blocks corresponding to complex conjugate pairs */
/*  of eigenvalues of (A,B).  The columns of Q are the left Schur vectors */
/*  and the columns of Z are the right Schur vectors. */

/*  If only the eigenvalues of (A,B) are needed, the driver routine */
/*  DGEGV should be used instead.  See DGEGV for a description of the */
/*  eigenvalues of the generalized nonsymmetric eigenvalue problem */
/*  (GNEP). */

/*  Arguments */
/*  ========= */

/*  JOBVSL  (input) CHARACTER*1 */
/*          = 'N':  do not compute the left Schur vectors; */
/*          = 'V':  compute the left Schur vectors (returned in VSL). */

/*  JOBVSR  (input) CHARACTER*1 */
/*          = 'N':  do not compute the right Schur vectors; */
/*          = 'V':  compute the right Schur vectors (returned in VSR). */

/*  N       (input) INTEGER */
/*          The order of the matrices A, B, VSL, and VSR.  N >= 0. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N) */
/*          On entry, the matrix A. */
/*          On exit, the upper quasi-triangular matrix S from the */
/*          generalized float Schur factorization. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of A.  LDA >= MAX(1,N). */

/*  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N) */
/*          On entry, the matrix B. */
/*          On exit, the upper triangular matrix T from the generalized */
/*          float Schur factorization. */

/*  LDB     (input) INTEGER */
/*          The leading dimension of B.  LDB >= MAX(1,N). */

/*  ALPHAR  (output) DOUBLE PRECISION array, dimension (N) */
/*          The float parts of each scalar alpha defining an eigenvalue */
/*          of GNEP. */

/*  ALPHAI  (output) DOUBLE PRECISION array, dimension (N) */
/*          The imaginary parts of each scalar alpha defining an */
/*          eigenvalue of GNEP.  If ALPHAI(j) is zero, then the j-th */
/*          eigenvalue is float; if positive, then the j-th and (j+1)-st */
/*          eigenvalues are a complex conjugate pair, with */
/*          ALPHAI(j+1) = -ALPHAI(j). */

/*  BETA    (output) DOUBLE PRECISION array, dimension (N) */
/*          The scalars beta that define the eigenvalues of GNEP. */
/*          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and */
/*          beta = BETA(j) represent the j-th eigenvalue of the matrix */
/*          pair (A,B), in one of the forms lambda = alpha/beta or */
/*          mu = beta/alpha.  Since either lambda or mu may overflow, */
/*          they should not, in general, be computed. */

/*  VSL     (output) DOUBLE PRECISION array, dimension (LDVSL,N) */
/*          If JOBVSL = 'V', the matrix of left Schur vectors Q. */
/*          Not referenced if JOBVSL = 'N'. */

/*  LDVSL   (input) INTEGER */
/*          The leading dimension of the matrix VSL. LDVSL >=1, and */
/*          if JOBVSL = 'V', LDVSL >= N. */

/*  VSR     (output) DOUBLE PRECISION array, dimension (LDVSR,N) */
/*          If JOBVSR = 'V', the matrix of right Schur vectors Z. */
/*          Not referenced if JOBVSR = 'N'. */

/*  LDVSR   (input) INTEGER */
/*          The leading dimension of the matrix VSR. LDVSR >= 1, and */
/*          if JOBVSR = 'V', LDVSR >= N. */

/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

/*  LWORK   (input) INTEGER */
/*          The dimension of the array WORK.  LWORK >= MAX(1,4*N). */
/*          For good performance, LWORK must generally be larger. */
/*          To compute the optimal value of LWORK, call ILAENV to get */
/*          blocksizes (for DGEQRF, DORMQR, and DORGQR.)  Then compute: */
/*          NB  -- MAX of the blocksizes for DGEQRF, DORMQR, and DORGQR */
/*          The optimal LWORK is  2*N + N*(NB+1). */

/*          If LWORK = -1, then a workspace query is assumed; the routine */
/*          only calculates the optimal size of the WORK array, returns */
/*          this value as the first entry of the WORK array, and no error */
/*          message related to LWORK is issued by XERBLA. */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/*          = 1,...,N: */
/*                The QZ iteration failed.  (A,B) are not in Schur */
/*                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should */
/*                be correct for j=INFO+1,...,N. */
/*          > N:  errors that usually indicate LAPACK problems: */
/*                =N+1: error return from DGGBAL */
/*                =N+2: error return from DGEQRF */
/*                =N+3: error return from DORMQR */
/*                =N+4: error return from DORGQR */
/*                =N+5: error return from DGGHRD */
/*                =N+6: error return from DHGEQZ (other than failed */
/*                                                iteration) */
/*                =N+7: error return from DGGBAK (computing VSL) */
/*                =N+8: error return from DGGBAK (computing VSR) */
/*                =N+9: error return from DLASCL (various places) */

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
    --alphar;
    --alphai;
    --beta;
    vsl_dim1 = *ldvsl;
    vsl_offset = 1 + vsl_dim1;
    vsl -= vsl_offset;
    vsr_dim1 = *ldvsr;
    vsr_offset = 1 + vsr_dim1;
    vsr -= vsr_offset;
    --work;

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
    i__1 = *n << 2;
    lwkmin = MAX(i__1,1);
    lwkopt = lwkmin;
    work[1] = (double) lwkopt;
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
	*info = -12;
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
	*info = -14;
    } else if (*lwork < lwkmin && ! lquery) {
	*info = -16;
    }

    if (*info == 0) {
	nb1 = ilaenv_(&c__1, "DGEQRF", " ", n, n, &c_n1, &c_n1);
	nb2 = ilaenv_(&c__1, "DORMQR", " ", n, n, n, &c_n1);
	nb3 = ilaenv_(&c__1, "DORGQR", " ", n, n, n, &c_n1);
/* Computing MAX */
	i__1 = MAX(nb1,nb2);
	nb = MAX(i__1,nb3);
	lopt = (*n << 1) + *n * (nb + 1);
	work[1] = (double) lopt;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEGS ", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Get machine constants */

    eps = dlamch_("E") * dlamch_("B");
    safmin = dlamch_("S");
    smlnum = *n * safmin / eps;
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

    anrm = dlange_("M", n, n, &a[a_offset], lda, &work[1]);
    ilascl = FALSE;
    if (anrm > 0. && anrm < smlnum) {
	anrmto = smlnum;
	ilascl = TRUE;
    } else if (anrm > bignum) {
	anrmto = bignum;
	ilascl = TRUE;
    }

    if (ilascl) {
	dlascl_("G", &c_n1, &c_n1, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

    bnrm = dlange_("M", n, n, &b[b_offset], ldb, &work[1]);
    ilbscl = FALSE;
    if (bnrm > 0. && bnrm < smlnum) {
	bnrmto = smlnum;
	ilbscl = TRUE;
    } else if (bnrm > bignum) {
	bnrmto = bignum;
	ilbscl = TRUE;
    }

    if (ilbscl) {
	dlascl_("G", &c_n1, &c_n1, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
    }

/*     Permute the matrix to make it more nearly triangular */
/*     Workspace layout:  (2*N words -- "work..." not actually used) */
/*        left_permutation, right_permutation, work... */

    ileft = 1;
    iright = *n + 1;
    iwork = iright + *n;
    dggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[
	    ileft], &work[iright], &work[iwork], &iinfo);
    if (iinfo != 0) {
	*info = *n + 1;
	goto L10;
    }

/*     Reduce B to triangular form, and initialize VSL and/or VSR */
/*     Workspace layout:  ("work..." must have at least N words) */
/*        left_permutation, right_permutation, tau, work... */

    irows = ihi + 1 - ilo;
    icols = *n + 1 - ilo;
    itau = iwork;
    iwork = itau + irows;
    i__1 = *lwork + 1 - iwork;
    dgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwork], &i__1, &iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__1 = lwkopt, i__2 = (int) work[iwork] + iwork - 1;
	lwkopt = MAX(i__1,i__2);
    }
    if (iinfo != 0) {
	*info = *n + 2;
	goto L10;
    }

    i__1 = *lwork + 1 - iwork;
    dormqr_("L", "T", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwork], &i__1, &
	    iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__1 = lwkopt, i__2 = (int) work[iwork] + iwork - 1;
	lwkopt = MAX(i__1,i__2);
    }
    if (iinfo != 0) {
	*info = *n + 3;
	goto L10;
    }

    if (ilvsl) {
	dlaset_("Full", n, n, &c_b36, &c_b37, &vsl[vsl_offset], ldvsl);
	i__1 = irows - 1;
	i__2 = irows - 1;
	dlacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[ilo 
		+ 1 + ilo * vsl_dim1], ldvsl);
	i__1 = *lwork + 1 - iwork;
	dorgqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &
		work[itau], &work[iwork], &i__1, &iinfo);
	if (iinfo >= 0) {
/* Computing MAX */
	    i__1 = lwkopt, i__2 = (int) work[iwork] + iwork - 1;
	    lwkopt = MAX(i__1,i__2);
	}
	if (iinfo != 0) {
	    *info = *n + 4;
	    goto L10;
	}
    }

    if (ilvsr) {
	dlaset_("Full", n, n, &c_b36, &c_b37, &vsr[vsr_offset], ldvsr);
    }

/*     Reduce to generalized Hessenberg form */

    dgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
	    ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &iinfo);
    if (iinfo != 0) {
	*info = *n + 5;
	goto L10;
    }

/*     Perform QZ algorithm, computing Schur vectors if desired */
/*     Workspace layout:  ("work..." must have at least 1 word) */
/*        left_permutation, right_permutation, work... */

    iwork = itau;
    i__1 = *lwork + 1 - iwork;
    dhgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[vsl_offset]
, ldvsl, &vsr[vsr_offset], ldvsr, &work[iwork], &i__1, &iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__1 = lwkopt, i__2 = (int) work[iwork] + iwork - 1;
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
	dggbak_("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vsl[
		vsl_offset], ldvsl, &iinfo);
	if (iinfo != 0) {
	    *info = *n + 7;
	    goto L10;
	}
    }
    if (ilvsr) {
	dggbak_("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vsr[
		vsr_offset], ldvsr, &iinfo);
	if (iinfo != 0) {
	    *info = *n + 8;
	    goto L10;
	}
    }

/*     Undo scaling */

    if (ilascl) {
	dlascl_("H", &c_n1, &c_n1, &anrmto, &anrm, n, n, &a[a_offset], lda, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
	dlascl_("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &alphar[1], n, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
	dlascl_("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &alphai[1], n, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
    }

    if (ilbscl) {
	dlascl_("U", &c_n1, &c_n1, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
	dlascl_("G", &c_n1, &c_n1, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return 0;
	}
    }

L10:
    work[1] = (double) lwkopt;

    return 0;

/*     End of DGEGS */

} /* dgegs_ */
