/* cgges.f -- translated by f2c (version 20061008).
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
static int c__0 = 0;
static int c_n1 = -1;

 int cgges_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, int *n, complex *a, int *lda, complex *b, int *
	ldb, int *sdim, complex *alpha, complex *beta, complex *vsl, 
	int *ldvsl, complex *vsr, int *ldvsr, complex *work, int *
	lwork, float *rwork, int *bwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, vsl_dim1, vsl_offset, 
	    vsr_dim1, vsr_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    int i__;
    float dif[2];
    int ihi, ilo;
    float eps, anrm, bnrm;
    int idum[1], ierr, itau, iwrk;
    float pvsl, pvsr;
    extern int lsame_(char *, char *);
    int ileft, icols;
    int cursl, ilvsl, ilvsr;
    int irwrk, irows;
    extern  int cggbak_(char *, char *, int *, int *, 
	    int *, float *, float *, int *, complex *, int *, 
	    int *), cggbal_(char *, int *, complex *, 
	    int *, complex *, int *, int *, int *, float *, 
	    float *, float *, int *), slabad_(float *, float *);
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
	    int *, int *, complex *, complex *, complex *, int *), xerbla_(char *, int *);
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *);
    float bignum;
    extern  int chgeqz_(char *, char *, char *, int *, 
	    int *, int *, complex *, int *, complex *, int *, 
	    complex *, complex *, complex *, int *, complex *, int *, 
	    complex *, int *, float *, int *), 
	    ctgsen_(int *, int *, int *, int *, int *, 
	    complex *, int *, complex *, int *, complex *, complex *, 
	    complex *, int *, complex *, int *, int *, float *, 
	    float *, float *, complex *, int *, int *, int *, 
	    int *);
    int ijobvl, iright, ijobvr;
    float anrmto;
    int lwkmin;
    int lastsl;
    float bnrmto;
    extern  int cungqr_(int *, int *, int *, 
	    complex *, int *, complex *, complex *, int *, int *),
	     cunmqr_(char *, char *, int *, int *, int *, complex 
	    *, int *, complex *, complex *, int *, complex *, int 
	    *, int *);
    float smlnum;
    int wantst, lquery;
    int lwkopt;


/*  -- LAPACK driver routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Function Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  CGGES computes for a pair of N-by-N complex nonsymmetric matrices */
/*  (A,B), the generalized eigenvalues, the generalized complex Schur */
/*  form (S, T), and optionally left and/or right Schur vectors (VSL */
/*  and VSR). This gives the generalized Schur factorization */

/*          (A,B) = ( (VSL)*S*(VSR)**H, (VSL)*T*(VSR)**H ) */

/*  where (VSR)**H is the conjugate-transpose of VSR. */

/*  Optionally, it also orders the eigenvalues so that a selected cluster */
/*  of eigenvalues appears in the leading diagonal blocks of the upper */
/*  triangular matrix S and the upper triangular matrix T. The leading */
/*  columns of VSL and VSR then form an unitary basis for the */
/*  corresponding left and right eigenspaces (deflating subspaces). */

/*  (If only the generalized eigenvalues are needed, use the driver */
/*  CGGEV instead, which is faster.) */

/*  A generalized eigenvalue for a pair of matrices (A,B) is a scalar w */
/*  or a ratio alpha/beta = w, such that  A - w*B is singular.  It is */
/*  usually represented as the pair (alpha,beta), as there is a */
/*  reasonable interpretation for beta=0, and even for both being zero. */

/*  A pair of matrices (S,T) is in generalized complex Schur form if S */
/*  and T are upper triangular and, in addition, the diagonal elements */
/*  of T are non-negative float numbers. */

/*  Arguments */
/*  ========= */

/*  JOBVSL  (input) CHARACTER*1 */
/*          = 'N':  do not compute the left Schur vectors; */
/*          = 'V':  compute the left Schur vectors. */

/*  JOBVSR  (input) CHARACTER*1 */
/*          = 'N':  do not compute the right Schur vectors; */
/*          = 'V':  compute the right Schur vectors. */

/*  SORT    (input) CHARACTER*1 */
/*          Specifies whether or not to order the eigenvalues on the */
/*          diagonal of the generalized Schur form. */
/*          = 'N':  Eigenvalues are not ordered; */
/*          = 'S':  Eigenvalues are ordered (see SELCTG). */

/*  SELCTG  (external procedure) LOGICAL FUNCTION of two COMPLEX arguments */
/*          SELCTG must be declared EXTERNAL in the calling subroutine. */
/*          If SORT = 'N', SELCTG is not referenced. */
/*          If SORT = 'S', SELCTG is used to select eigenvalues to sort */
/*          to the top left of the Schur form. */
/*          An eigenvalue ALPHA(j)/BETA(j) is selected if */
/*          SELCTG(ALPHA(j),BETA(j)) is true. */

/*          Note that a selected complex eigenvalue may no longer satisfy */
/*          SELCTG(ALPHA(j),BETA(j)) = .TRUE. after ordering, since */
/*          ordering may change the value of complex eigenvalues */
/*          (especially if the eigenvalue is ill-conditioned), in this */
/*          case INFO is set to N+2 (See INFO below). */

/*  N       (input) INTEGER */
/*          The order of the matrices A, B, VSL, and VSR.  N >= 0. */

/*  A       (input/output) COMPLEX array, dimension (LDA, N) */
/*          On entry, the first of the pair of matrices. */
/*          On exit, A has been overwritten by its generalized Schur */
/*          form S. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of A.  LDA >= MAX(1,N). */

/*  B       (input/output) COMPLEX array, dimension (LDB, N) */
/*          On entry, the second of the pair of matrices. */
/*          On exit, B has been overwritten by its generalized Schur */
/*          form T. */

/*  LDB     (input) INTEGER */
/*          The leading dimension of B.  LDB >= MAX(1,N). */

/*  SDIM    (output) INTEGER */
/*          If SORT = 'N', SDIM = 0. */
/*          If SORT = 'S', SDIM = number of eigenvalues (after sorting) */
/*          for which SELCTG is true. */

/*  ALPHA   (output) COMPLEX array, dimension (N) */
/*  BETA    (output) COMPLEX array, dimension (N) */
/*          On exit,  ALPHA(j)/BETA(j), j=1,...,N, will be the */
/*          generalized eigenvalues.  ALPHA(j), j=1,...,N  and  BETA(j), */
/*          j=1,...,N  are the diagonals of the complex Schur form (A,B) */
/*          output by CGGES. The  BETA(j) will be non-negative float. */

/*          Note: the quotients ALPHA(j)/BETA(j) may easily over- or */
/*          underflow, and BETA(j) may even be zero.  Thus, the user */
/*          should avoid naively computing the ratio alpha/beta. */
/*          However, ALPHA will be always less than and usually */
/*          comparable with norm(A) in magnitude, and BETA always less */
/*          than and usually comparable with norm(B). */

/*  VSL     (output) COMPLEX array, dimension (LDVSL,N) */
/*          If JOBVSL = 'V', VSL will contain the left Schur vectors. */
/*          Not referenced if JOBVSL = 'N'. */

/*  LDVSL   (input) INTEGER */
/*          The leading dimension of the matrix VSL. LDVSL >= 1, and */
/*          if JOBVSL = 'V', LDVSL >= N. */

/*  VSR     (output) COMPLEX array, dimension (LDVSR,N) */
/*          If JOBVSR = 'V', VSR will contain the right Schur vectors. */
/*          Not referenced if JOBVSR = 'N'. */

/*  LDVSR   (input) INTEGER */
/*          The leading dimension of the matrix VSR. LDVSR >= 1, and */
/*          if JOBVSR = 'V', LDVSR >= N. */

/*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK)) */
/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

/*  LWORK   (input) INTEGER */
/*          The dimension of the array WORK.  LWORK >= MAX(1,2*N). */
/*          For good performance, LWORK must generally be larger. */

/*          If LWORK = -1, then a workspace query is assumed; the routine */
/*          only calculates the optimal size of the WORK array, returns */
/*          this value as the first entry of the WORK array, and no error */
/*          message related to LWORK is issued by XERBLA. */

/*  RWORK   (workspace) REAL array, dimension (8*N) */

/*  BWORK   (workspace) LOGICAL array, dimension (N) */
/*          Not referenced if SORT = 'N'. */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/*          =1,...,N: */
/*                The QZ iteration failed.  (A,B) are not in Schur */
/*                form, but ALPHA(j) and BETA(j) should be correct for */
/*                j=INFO+1,...,N. */
/*          > N:  =N+1: other than QZ iteration failed in CHGEQZ */
/*                =N+2: after reordering, roundoff changed values of */
/*                      some complex eigenvalues so that leading */
/*                      eigenvalues in the Generalized Schur form no */
/*                      longer satisfy SELCTG=.TRUE.  This could also */
/*                      be caused due to scaling. */
/*                =N+3: reordering falied in CTGSEN. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
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
    --bwork;

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

    wantst = lsame_(sort, "S");

/*     Test the input arguments */

    *info = 0;
    lquery = *lwork == -1;
    if (ijobvl <= 0) {
	*info = -1;
    } else if (ijobvr <= 0) {
	*info = -2;
    } else if (! wantst && ! lsame_(sort, "N")) {
	*info = -3;
    } else if (*n < 0) {
	*info = -5;
    } else if (*lda < MAX(1,*n)) {
	*info = -7;
    } else if (*ldb < MAX(1,*n)) {
	*info = -9;
    } else if (*ldvsl < 1 || ilvsl && *ldvsl < *n) {
	*info = -14;
    } else if (*ldvsr < 1 || ilvsr && *ldvsr < *n) {
	*info = -16;
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

    if (*info == 0) {
/* Computing MAX */
	i__1 = 1, i__2 = *n << 1;
	lwkmin = MAX(i__1,i__2);
/* Computing MAX */
	i__1 = 1, i__2 = *n + *n * ilaenv_(&c__1, "CGEQRF", " ", n, &c__1, n, 
		&c__0);
	lwkopt = MAX(i__1,i__2);
/* Computing MAX */
	i__1 = lwkopt, i__2 = *n + *n * ilaenv_(&c__1, "CUNMQR", " ", n, &
		c__1, n, &c_n1);
	lwkopt = MAX(i__1,i__2);
	if (ilvsl) {
/* Computing MAX */
	    i__1 = lwkopt, i__2 = *n + *n * ilaenv_(&c__1, "CUNGQR", " ", n, &
		    c__1, n, &c_n1);
	    lwkopt = MAX(i__1,i__2);
	}
	work[1].r = (float) lwkopt, work[1].i = 0.f;

	if (*lwork < lwkmin && ! lquery) {
	    *info = -18;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CGGES ", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	*sdim = 0;
	return 0;
    }

/*     Get machine constants */

    eps = slamch_("P");
    smlnum = slamch_("S");
    bignum = 1.f / smlnum;
    slabad_(&smlnum, &bignum);
    smlnum = sqrt(smlnum) / eps;
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
	clascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr);
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
	clascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr);
    }

/*     Permute the matrix to make it more nearly triangular */
/*     (Real Workspace: need 6*N) */

    ileft = 1;
    iright = *n + 1;
    irwrk = iright + *n;
    cggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[
	    ileft], &rwork[iright], &rwork[irwrk], &ierr);

/*     Reduce B to triangular form (QR decomposition of B) */
/*     (Complex Workspace: need N, prefer N*NB) */

    irows = ihi + 1 - ilo;
    icols = *n + 1 - ilo;
    itau = 1;
    iwrk = itau + irows;
    i__1 = *lwork + 1 - iwrk;
    cgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwrk], &i__1, &ierr);

/*     Apply the orthogonal transformation to matrix A */
/*     (Complex Workspace: need N, prefer N*NB) */

    i__1 = *lwork + 1 - iwrk;
    cunmqr_("L", "C", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &
	    ierr);

/*     Initialize VSL */
/*     (Complex Workspace: need N, prefer N*NB) */

    if (ilvsl) {
	claset_("Full", n, n, &c_b1, &c_b2, &vsl[vsl_offset], ldvsl);
	if (irows > 1) {
	    i__1 = irows - 1;
	    i__2 = irows - 1;
	    clacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[
		    ilo + 1 + ilo * vsl_dim1], ldvsl);
	}
	i__1 = *lwork + 1 - iwrk;
	cungqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &
		work[itau], &work[iwrk], &i__1, &ierr);
    }

/*     Initialize VSR */

    if (ilvsr) {
	claset_("Full", n, n, &c_b1, &c_b2, &vsr[vsr_offset], ldvsr);
    }

/*     Reduce to generalized Hessenberg form */
/*     (Workspace: none needed) */

    cgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
	    ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &ierr);

    *sdim = 0;

/*     Perform QZ algorithm, computing Schur vectors if desired */
/*     (Complex Workspace: need N) */
/*     (Real Workspace: need N) */

    iwrk = itau;
    i__1 = *lwork + 1 - iwrk;
    chgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, &
	    vsr[vsr_offset], ldvsr, &work[iwrk], &i__1, &rwork[irwrk], &ierr);
    if (ierr != 0) {
	if (ierr > 0 && ierr <= *n) {
	    *info = ierr;
	} else if (ierr > *n && ierr <= *n << 1) {
	    *info = ierr - *n;
	} else {
	    *info = *n + 1;
	}
	goto L30;
    }

/*     Sort eigenvalues ALPHA/BETA if desired */
/*     (Workspace: none needed) */

    if (wantst) {

/*        Undo scaling on eigenvalues before selecting */

	if (ilascl) {
	    clascl_("G", &c__0, &c__0, &anrm, &anrmto, n, &c__1, &alpha[1], n, 
		     &ierr);
	}
	if (ilbscl) {
	    clascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, &c__1, &beta[1], n, 
		    &ierr);
	}

/*        Select eigenvalues */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    bwork[i__] = (*selctg)(&alpha[i__], &beta[i__]);
/* L10: */
	}

	i__1 = *lwork - iwrk + 1;
	ctgsen_(&c__0, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &b[
		b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, 
		&vsr[vsr_offset], ldvsr, sdim, &pvsl, &pvsr, dif, &work[iwrk], 
		 &i__1, idum, &c__1, &ierr);
	if (ierr == 1) {
	    *info = *n + 3;
	}

    }

/*     Apply back-permutation to VSL and VSR */
/*     (Workspace: none needed) */

    if (ilvsl) {
	cggbak_("P", "L", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsl[vsl_offset], ldvsl, &ierr);
    }
    if (ilvsr) {
	cggbak_("P", "R", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &
		vsr[vsr_offset], ldvsr, &ierr);
    }

/*     Undo scaling */

    if (ilascl) {
	clascl_("U", &c__0, &c__0, &anrmto, &anrm, n, n, &a[a_offset], lda, &
		ierr);
	clascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alpha[1], n, &
		ierr);
    }

    if (ilbscl) {
	clascl_("U", &c__0, &c__0, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &
		ierr);
	clascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr);
    }

    if (wantst) {

/*        Check if reordering is correct */

	lastsl = TRUE;
	*sdim = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    cursl = (*selctg)(&alpha[i__], &beta[i__]);
	    if (cursl) {
		++(*sdim);
	    }
	    if (cursl && ! lastsl) {
		*info = *n + 2;
	    }
	    lastsl = cursl;
/* L20: */
	}

    }

L30:

    work[1].r = (float) lwkopt, work[1].i = 0.f;

    return 0;

/*     End of CGGES */

} /* cgges_ */
