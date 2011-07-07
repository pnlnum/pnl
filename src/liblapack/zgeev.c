/* zgeev.f -- translated by f2c (version 20061008).
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
static int c__0 = 0;
static int c_n1 = -1;

 int zgeev_(char *jobvl, char *jobvr, int *n, 
	doublecomplex *a, int *lda, doublecomplex *w, doublecomplex *vl, 
	int *ldvl, doublecomplex *vr, int *ldvr, doublecomplex *work, 
	int *lwork, double *rwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, 
	    i__2, i__3;
    double d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double sqrt(double), d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    int i__, k, ihi;
    double scl;
    int ilo;
    double dum[1], eps;
    doublecomplex tmp;
    int ibal;
    char side[1];
    double anrm;
    int ierr, itau, iwrk, nout;
    extern int lsame_(char *, char *);
    extern  int zscal_(int *, doublecomplex *, 
	    doublecomplex *, int *), dlabad_(double *, double *);
    extern double dznrm2_(int *, doublecomplex *, int *);
    int scalea;
    extern double dlamch_(char *);
    double cscale;
    extern  int zgebak_(char *, char *, int *, int *, 
	    int *, double *, int *, doublecomplex *, int *, 
	    int *), zgebal_(char *, int *, 
	    doublecomplex *, int *, int *, int *, double *, 
	    int *);
    extern int idamax_(int *, double *, int *);
    extern  int xerbla_(char *, int *);
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *);
    int select[1];
    extern  int zdscal_(int *, double *, 
	    doublecomplex *, int *);
    double bignum;
    extern double zlange_(char *, int *, int *, doublecomplex *, 
	    int *, double *);
    extern  int zgehrd_(int *, int *, int *, 
	    doublecomplex *, int *, doublecomplex *, doublecomplex *, 
	    int *, int *), zlascl_(char *, int *, int *, 
	    double *, double *, int *, int *, doublecomplex *, 
	     int *, int *), zlacpy_(char *, int *, 
	    int *, doublecomplex *, int *, doublecomplex *, int *);
    int minwrk, maxwrk;
    int wantvl;
    double smlnum;
    int hswork, irwork;
    extern  int zhseqr_(char *, char *, int *, int *, 
	    int *, doublecomplex *, int *, doublecomplex *, 
	    doublecomplex *, int *, doublecomplex *, int *, int *), ztrevc_(char *, char *, int *, int *, 
	    doublecomplex *, int *, doublecomplex *, int *, 
	    doublecomplex *, int *, int *, int *, doublecomplex *, 
	     double *, int *);
    int lquery, wantvr;
    extern  int zunghr_(int *, int *, int *, 
	    doublecomplex *, int *, doublecomplex *, doublecomplex *, 
	    int *, int *);


/*  -- LAPACK driver routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  ZGEEV computes for an N-by-N complex nonsymmetric matrix A, the */
/*  eigenvalues and, optionally, the left and/or right eigenvectors. */

/*  The right eigenvector v(j) of A satisfies */
/*                   A * v(j) = lambda(j) * v(j) */
/*  where lambda(j) is its eigenvalue. */
/*  The left eigenvector u(j) of A satisfies */
/*                u(j)**H * A = lambda(j) * u(j)**H */
/*  where u(j)**H denotes the conjugate transpose of u(j). */

/*  The computed eigenvectors are normalized to have Euclidean norm */
/*  equal to 1 and largest component float. */

/*  Arguments */
/*  ========= */

/*  JOBVL   (input) CHARACTER*1 */
/*          = 'N': left eigenvectors of A are not computed; */
/*          = 'V': left eigenvectors of are computed. */

/*  JOBVR   (input) CHARACTER*1 */
/*          = 'N': right eigenvectors of A are not computed; */
/*          = 'V': right eigenvectors of A are computed. */

/*  N       (input) INTEGER */
/*          The order of the matrix A. N >= 0. */

/*  A       (input/output) COMPLEX*16 array, dimension (LDA,N) */
/*          On entry, the N-by-N matrix A. */
/*          On exit, A has been overwritten. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= MAX(1,N). */

/*  W       (output) COMPLEX*16 array, dimension (N) */
/*          W contains the computed eigenvalues. */

/*  VL      (output) COMPLEX*16 array, dimension (LDVL,N) */
/*          If JOBVL = 'V', the left eigenvectors u(j) are stored one */
/*          after another in the columns of VL, in the same order */
/*          as their eigenvalues. */
/*          If JOBVL = 'N', VL is not referenced. */
/*          u(j) = VL(:,j), the j-th column of VL. */

/*  LDVL    (input) INTEGER */
/*          The leading dimension of the array VL.  LDVL >= 1; if */
/*          JOBVL = 'V', LDVL >= N. */

/*  VR      (output) COMPLEX*16 array, dimension (LDVR,N) */
/*          If JOBVR = 'V', the right eigenvectors v(j) are stored one */
/*          after another in the columns of VR, in the same order */
/*          as their eigenvalues. */
/*          If JOBVR = 'N', VR is not referenced. */
/*          v(j) = VR(:,j), the j-th column of VR. */

/*  LDVR    (input) INTEGER */
/*          The leading dimension of the array VR.  LDVR >= 1; if */
/*          JOBVR = 'V', LDVR >= N. */

/*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

/*  LWORK   (input) INTEGER */
/*          The dimension of the array WORK.  LWORK >= MAX(1,2*N). */
/*          For good performance, LWORK must generally be larger. */

/*          If LWORK = -1, then a workspace query is assumed; the routine */
/*          only calculates the optimal size of the WORK array, returns */
/*          this value as the first entry of the WORK array, and no error */
/*          message related to LWORK is issued by XERBLA. */

/*  RWORK   (workspace) DOUBLE PRECISION array, dimension (2*N) */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/*          > 0:  if INFO = i, the QR algorithm failed to compute all the */
/*                eigenvalues, and no eigenvectors have been computed; */
/*                elements and i+1:N of W contain eigenvalues which have */
/*                converged. */

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

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    vl_dim1 = *ldvl;
    vl_offset = 1 + vl_dim1;
    vl -= vl_offset;
    vr_dim1 = *ldvr;
    vr_offset = 1 + vr_dim1;
    vr -= vr_offset;
    --work;
    --rwork;

    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    wantvl = lsame_(jobvl, "V");
    wantvr = lsame_(jobvr, "V");
    if (! wantvl && ! lsame_(jobvl, "N")) {
	*info = -1;
    } else if (! wantvr && ! lsame_(jobvr, "N")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < MAX(1,*n)) {
	*info = -5;
    } else if (*ldvl < 1 || wantvl && *ldvl < *n) {
	*info = -8;
    } else if (*ldvr < 1 || wantvr && *ldvr < *n) {
	*info = -10;
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       CWorkspace refers to complex workspace, and RWorkspace to float */
/*       workspace. NB refers to the optimal block size for the */
/*       immediately following subroutine, as returned by ILAENV. */
/*       HSWORK refers to the workspace preferred by ZHSEQR, as */
/*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N, */
/*       the worst case.) */

    if (*info == 0) {
	if (*n == 0) {
	    minwrk = 1;
	    maxwrk = 1;
	} else {
	    maxwrk = *n + *n * ilaenv_(&c__1, "ZGEHRD", " ", n, &c__1, n, &
		    c__0);
	    minwrk = *n << 1;
	    if (wantvl) {
/* Computing MAX */
		i__1 = maxwrk, i__2 = *n + (*n - 1) * ilaenv_(&c__1, "ZUNGHR", 
			 " ", n, &c__1, n, &c_n1);
		maxwrk = MAX(i__1,i__2);
		zhseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &w[1], &vl[
			vl_offset], ldvl, &work[1], &c_n1, info);
	    } else if (wantvr) {
/* Computing MAX */
		i__1 = maxwrk, i__2 = *n + (*n - 1) * ilaenv_(&c__1, "ZUNGHR", 
			 " ", n, &c__1, n, &c_n1);
		maxwrk = MAX(i__1,i__2);
		zhseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &w[1], &vr[
			vr_offset], ldvr, &work[1], &c_n1, info);
	    } else {
		zhseqr_("E", "N", n, &c__1, n, &a[a_offset], lda, &w[1], &vr[
			vr_offset], ldvr, &work[1], &c_n1, info);
	    }
	    hswork = (int) work[1].r;
/* Computing MAX */
	    i__1 = MAX(maxwrk,hswork);
	    maxwrk = MAX(i__1,minwrk);
	}
	work[1].r = (double) maxwrk, work[1].i = 0.;

	if (*lwork < minwrk && ! lquery) {
	    *info = -12;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZGEEV ", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Get machine constants */

    eps = dlamch_("P");
    smlnum = dlamch_("S");
    bignum = 1. / smlnum;
    dlabad_(&smlnum, &bignum);
    smlnum = sqrt(smlnum) / eps;
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

    anrm = zlange_("M", n, n, &a[a_offset], lda, dum);
    scalea = FALSE;
    if (anrm > 0. && anrm < smlnum) {
	scalea = TRUE;
	cscale = smlnum;
    } else if (anrm > bignum) {
	scalea = TRUE;
	cscale = bignum;
    }
    if (scalea) {
	zlascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &
		ierr);
    }

/*     Balance the matrix */
/*     (CWorkspace: none) */
/*     (RWorkspace: need N) */

    ibal = 1;
    zgebal_("B", n, &a[a_offset], lda, &ilo, &ihi, &rwork[ibal], &ierr);

/*     Reduce to upper Hessenberg form */
/*     (CWorkspace: need 2*N, prefer N+N*NB) */
/*     (RWorkspace: none) */

    itau = 1;
    iwrk = itau + *n;
    i__1 = *lwork - iwrk + 1;
    zgehrd_(n, &ilo, &ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1, 
	     &ierr);

    if (wantvl) {

/*        Want left eigenvectors */
/*        Copy Householder vectors to VL */

	*(unsigned char *)side = 'L';
	zlacpy_("L", n, n, &a[a_offset], lda, &vl[vl_offset], ldvl)
		;

/*        Generate unitary matrix in VL */
/*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
/*        (RWorkspace: none) */

	i__1 = *lwork - iwrk + 1;
	zunghr_(n, &ilo, &ihi, &vl[vl_offset], ldvl, &work[itau], &work[iwrk], 
		 &i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VL */
/*        (CWorkspace: need 1, prefer HSWORK (see comments) ) */
/*        (RWorkspace: none) */

	iwrk = itau;
	i__1 = *lwork - iwrk + 1;
	zhseqr_("S", "V", n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vl[
		vl_offset], ldvl, &work[iwrk], &i__1, info);

	if (wantvr) {

/*           Want left and right eigenvectors */
/*           Copy Schur vectors to VR */

	    *(unsigned char *)side = 'B';
	    zlacpy_("F", n, n, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr);
	}

    } else if (wantvr) {

/*        Want right eigenvectors */
/*        Copy Householder vectors to VR */

	*(unsigned char *)side = 'R';
	zlacpy_("L", n, n, &a[a_offset], lda, &vr[vr_offset], ldvr)
		;

/*        Generate unitary matrix in VR */
/*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
/*        (RWorkspace: none) */

	i__1 = *lwork - iwrk + 1;
	zunghr_(n, &ilo, &ihi, &vr[vr_offset], ldvr, &work[itau], &work[iwrk], 
		 &i__1, &ierr);

/*        Perform QR iteration, accumulating Schur vectors in VR */
/*        (CWorkspace: need 1, prefer HSWORK (see comments) ) */
/*        (RWorkspace: none) */

	iwrk = itau;
	i__1 = *lwork - iwrk + 1;
	zhseqr_("S", "V", n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vr[
		vr_offset], ldvr, &work[iwrk], &i__1, info);

    } else {

/*        Compute eigenvalues only */
/*        (CWorkspace: need 1, prefer HSWORK (see comments) ) */
/*        (RWorkspace: none) */

	iwrk = itau;
	i__1 = *lwork - iwrk + 1;
	zhseqr_("E", "N", n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vr[
		vr_offset], ldvr, &work[iwrk], &i__1, info);
    }

/*     If INFO > 0 from ZHSEQR, then quit */

    if (*info > 0) {
	goto L50;
    }

    if (wantvl || wantvr) {

/*        Compute left and/or right eigenvectors */
/*        (CWorkspace: need 2*N) */
/*        (RWorkspace: need 2*N) */

	irwork = ibal + *n;
	ztrevc_(side, "B", select, n, &a[a_offset], lda, &vl[vl_offset], ldvl, 
		 &vr[vr_offset], ldvr, n, &nout, &work[iwrk], &rwork[irwork], 
		&ierr);
    }

    if (wantvl) {

/*        Undo balancing of left eigenvectors */
/*        (CWorkspace: none) */
/*        (RWorkspace: need N) */

	zgebak_("B", "L", n, &ilo, &ihi, &rwork[ibal], n, &vl[vl_offset], 
		ldvl, &ierr);

/*        Normalize left eigenvectors and make largest component float */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    scl = 1. / dznrm2_(n, &vl[i__ * vl_dim1 + 1], &c__1);
	    zdscal_(n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
	    i__2 = *n;
	    for (k = 1; k <= i__2; ++k) {
		i__3 = k + i__ * vl_dim1;
/* Computing 2nd power */
		d__1 = vl[i__3].r;
/* Computing 2nd power */
		d__2 = d_imag(&vl[k + i__ * vl_dim1]);
		rwork[irwork + k - 1] = d__1 * d__1 + d__2 * d__2;
/* L10: */
	    }
	    k = idamax_(n, &rwork[irwork], &c__1);
	    d_cnjg(&z__2, &vl[k + i__ * vl_dim1]);
	    d__1 = sqrt(rwork[irwork + k - 1]);
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
	    tmp.r = z__1.r, tmp.i = z__1.i;
	    zscal_(n, &tmp, &vl[i__ * vl_dim1 + 1], &c__1);
	    i__2 = k + i__ * vl_dim1;
	    i__3 = k + i__ * vl_dim1;
	    d__1 = vl[i__3].r;
	    z__1.r = d__1, z__1.i = 0.;
	    vl[i__2].r = z__1.r, vl[i__2].i = z__1.i;
/* L20: */
	}
    }

    if (wantvr) {

/*        Undo balancing of right eigenvectors */
/*        (CWorkspace: none) */
/*        (RWorkspace: need N) */

	zgebak_("B", "R", n, &ilo, &ihi, &rwork[ibal], n, &vr[vr_offset], 
		ldvr, &ierr);

/*        Normalize right eigenvectors and make largest component float */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    scl = 1. / dznrm2_(n, &vr[i__ * vr_dim1 + 1], &c__1);
	    zdscal_(n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
	    i__2 = *n;
	    for (k = 1; k <= i__2; ++k) {
		i__3 = k + i__ * vr_dim1;
/* Computing 2nd power */
		d__1 = vr[i__3].r;
/* Computing 2nd power */
		d__2 = d_imag(&vr[k + i__ * vr_dim1]);
		rwork[irwork + k - 1] = d__1 * d__1 + d__2 * d__2;
/* L30: */
	    }
	    k = idamax_(n, &rwork[irwork], &c__1);
	    d_cnjg(&z__2, &vr[k + i__ * vr_dim1]);
	    d__1 = sqrt(rwork[irwork + k - 1]);
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
	    tmp.r = z__1.r, tmp.i = z__1.i;
	    zscal_(n, &tmp, &vr[i__ * vr_dim1 + 1], &c__1);
	    i__2 = k + i__ * vr_dim1;
	    i__3 = k + i__ * vr_dim1;
	    d__1 = vr[i__3].r;
	    z__1.r = d__1, z__1.i = 0.;
	    vr[i__2].r = z__1.r, vr[i__2].i = z__1.i;
/* L40: */
	}
    }

/*     Undo scaling if necessary */

L50:
    if (scalea) {
	i__1 = *n - *info;
/* Computing MAX */
	i__3 = *n - *info;
	i__2 = MAX(i__3,1);
	zlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &w[*info + 1]
, &i__2, &ierr);
	if (*info > 0) {
	    i__1 = ilo - 1;
	    zlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &w[1], n, 
		     &ierr);
	}
    }

    work[1].r = (double) maxwrk, work[1].i = 0.;
    return 0;

/*     End of ZGEEV */

} /* zgeev_ */
