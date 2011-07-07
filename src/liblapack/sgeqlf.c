/* sgeqlf.f -- translated by f2c (version 20061008).
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
static int c__3 = 3;
static int c__2 = 2;

 int sgeqlf_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *lwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    int i__, k, ib, nb, ki, kk, mu, nu, nx, iws, nbmin, iinfo;
    extern  int sgeql2_(int *, int *, float *, int 
	    *, float *, float *, int *), slarfb_(char *, char *, char *, 
	    char *, int *, int *, int *, float *, int *, float *
, int *, float *, int *, float *, int *), xerbla_(char *, int *);
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *);
    extern  int slarft_(char *, char *, int *, int *, 
	    float *, int *, float *, float *, int *);
    int ldwork, lwkopt;
    int lquery;


/*  -- LAPACK routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SGEQLF computes a QL factorization of a float M-by-N matrix A: */
/*  A = Q * L. */

/*  Arguments */
/*  ========= */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix A.  M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A.  N >= 0. */

/*  A       (input/output) REAL array, dimension (LDA,N) */
/*          On entry, the M-by-N matrix A. */
/*          On exit, */
/*          if m >= n, the lower triangle of the subarray */
/*          A(m-n+1:m,1:n) contains the N-by-N lower triangular matrix L; */
/*          if m <= n, the elements on and below the (n-m)-th */
/*          superdiagonal contain the M-by-N lower trapezoidal matrix L; */
/*          the remaining elements, with the array TAU, represent the */
/*          orthogonal matrix Q as a product of elementary reflectors */
/*          (see Further Details). */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= MAX(1,M). */

/*  TAU     (output) REAL array, dimension (MIN(M,N)) */
/*          The scalar factors of the elementary reflectors (see Further */
/*          Details). */

/*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK)) */
/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

/*  LWORK   (input) INTEGER */
/*          The dimension of the array WORK.  LWORK >= MAX(1,N). */
/*          For optimum performance LWORK >= N*NB, where NB is the */
/*          optimal blocksize. */

/*          If LWORK = -1, then a workspace query is assumed; the routine */
/*          only calculates the optimal size of the WORK array, returns */
/*          this value as the first entry of the WORK array, and no error */
/*          message related to LWORK is issued by XERBLA. */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */

/*  Further Details */
/*  =============== */

/*  The matrix Q is represented as a product of elementary reflectors */

/*     Q = H(k) . . . H(2) H(1), where k = MIN(m,n). */

/*  Each H(i) has the form */

/*     H(i) = I - tau * v * v' */

/*  where tau is a float scalar, and v is a float vector with */
/*  v(m-k+i+1:m) = 0 and v(m-k+i) = 1; v(1:m-k+i-1) is stored on exit in */
/*  A(1:m-k+i-1,n-k+i), and tau in TAU(i). */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < MAX(1,*m)) {
	*info = -4;
    }

    if (*info == 0) {
	k = MIN(*m,*n);
	if (k == 0) {
	    lwkopt = 1;
	} else {
	    nb = ilaenv_(&c__1, "SGEQLF", " ", m, n, &c_n1, &c_n1);
	    lwkopt = *n * nb;
	}
	work[1] = (float) lwkopt;

	if (*lwork < MAX(1,*n) && ! lquery) {
	    *info = -7;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SGEQLF", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (k == 0) {
	return 0;
    }

    nbmin = 2;
    nx = 1;
    iws = *n;
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
	i__1 = 0, i__2 = ilaenv_(&c__3, "SGEQLF", " ", m, n, &c_n1, &c_n1);
	nx = MAX(i__1,i__2);
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = ilaenv_(&c__2, "SGEQLF", " ", m, n, &c_n1, &
			c_n1);
		nbmin = MAX(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially. */
/*        The last kk columns are handled by the block method. */

	ki = (k - nx - 1) / nb * nb;
/* Computing MIN */
	i__1 = k, i__2 = ki + nb;
	kk = MIN(i__1,i__2);

	i__1 = k - kk + 1;
	i__2 = -nb;
	for (i__ = k - kk + ki + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ 
		+= i__2) {
/* Computing MIN */
	    i__3 = k - i__ + 1;
	    ib = MIN(i__3,nb);

/*           Compute the QL factorization of the current block */
/*           A(1:m-k+i+ib-1,n-k+i:n-k+i+ib-1) */

	    i__3 = *m - k + i__ + ib - 1;
	    sgeql2_(&i__3, &ib, &a[(*n - k + i__) * a_dim1 + 1], lda, &tau[
		    i__], &work[1], &iinfo);
	    if (*n - k + i__ > 1) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i+ib-1) . . . H(i+1) H(i) */

		i__3 = *m - k + i__ + ib - 1;
		slarft_("Backward", "Columnwise", &i__3, &ib, &a[(*n - k + 
			i__) * a_dim1 + 1], lda, &tau[i__], &work[1], &ldwork);

/*              Apply H' to A(1:m-k+i+ib-1,1:n-k+i-1) from the left */

		i__3 = *m - k + i__ + ib - 1;
		i__4 = *n - k + i__ - 1;
		slarfb_("Left", "Transpose", "Backward", "Columnwise", &i__3, 
			&i__4, &ib, &a[(*n - k + i__) * a_dim1 + 1], lda, &
			work[1], &ldwork, &a[a_offset], lda, &work[ib + 1], &
			ldwork);
	    }
/* L10: */
	}
	mu = *m - k + i__ + nb - 1;
	nu = *n - k + i__ + nb - 1;
    } else {
	mu = *m;
	nu = *n;
    }

/*     Use unblocked code to factor the last or only block */

    if (mu > 0 && nu > 0) {
	sgeql2_(&mu, &nu, &a[a_offset], lda, &tau[1], &work[1], &iinfo);
    }

    work[1] = (float) iws;
    return 0;

/*     End of SGEQLF */

} /* sgeqlf_ */
