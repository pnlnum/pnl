/* dggsvd.f -- translated by f2c (version 20061008).
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

 int dggsvd_(char *jobu, char *jobv, char *jobq, int *m, 
	int *n, int *p, int *k, int *l, double *a, 
	int *lda, double *b, int *ldb, double *alpha, 
	double *beta, double *u, int *ldu, double *v, int 
	*ldv, double *q, int *ldq, double *work, int *iwork, 
	int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, u_dim1, 
	    u_offset, v_dim1, v_offset, i__1, i__2;

    /* Local variables */
    int i__, j;
    double ulp;
    int ibnd;
    double tola;
    int isub;
    double tolb, unfl, temp, smax;
    extern int lsame_(char *, char *);
    double anorm, bnorm;
    extern  int dcopy_(int *, double *, int *, 
	    double *, int *);
    int wantq, wantu, wantv;
    extern double dlamch_(char *), dlange_(char *, int *, 
	    int *, double *, int *, double *);
    extern  int dtgsja_(char *, char *, char *, int *, 
	    int *, int *, int *, int *, double *, int 
	    *, double *, int *, double *, double *, 
	    double *, double *, double *, int *, double *, 
	     int *, double *, int *, double *, int *, 
	    int *);
    int ncycle;
    extern  int xerbla_(char *, int *), dggsvp_(
	    char *, char *, char *, int *, int *, int *, 
	    double *, int *, double *, int *, double *, 
	    double *, int *, int *, double *, int *, 
	    double *, int *, double *, int *, int *, 
	    double *, double *, int *);


/*  -- LAPACK driver routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGGSVD computes the generalized singular value decomposition (GSVD) */
/*  of an M-by-N float matrix A and P-by-N float matrix B: */

/*      U'*A*Q = D1*( 0 R ),    V'*B*Q = D2*( 0 R ) */

/*  where U, V and Q are orthogonal matrices, and Z' is the transpose */
/*  of Z.  Let K+L = the effective numerical rank of the matrix (A',B')', */
/*  then R is a K+L-by-K+L nonsingular upper triangular matrix, D1 and */
/*  D2 are M-by-(K+L) and P-by-(K+L) "diagonal" matrices and of the */
/*  following structures, respectively: */

/*  If M-K-L >= 0, */

/*                      K  L */
/*         D1 =     K ( I  0 ) */
/*                  L ( 0  C ) */
/*              M-K-L ( 0  0 ) */

/*                    K  L */
/*         D2 =   L ( 0  S ) */
/*              P-L ( 0  0 ) */

/*                  N-K-L  K    L */
/*    ( 0 R ) = K (  0   R11  R12 ) */
/*              L (  0    0   R22 ) */

/*  where */

/*    C = diag( ALPHA(K+1), ... , ALPHA(K+L) ), */
/*    S = diag( BETA(K+1),  ... , BETA(K+L) ), */
/*    C**2 + S**2 = I. */

/*    R is stored in A(1:K+L,N-K-L+1:N) on exit. */

/*  If M-K-L < 0, */

/*                    K M-K K+L-M */
/*         D1 =   K ( I  0    0   ) */
/*              M-K ( 0  C    0   ) */

/*                      K M-K K+L-M */
/*         D2 =   M-K ( 0  S    0  ) */
/*              K+L-M ( 0  0    I  ) */
/*                P-L ( 0  0    0  ) */

/*                     N-K-L  K   M-K  K+L-M */
/*    ( 0 R ) =     K ( 0    R11  R12  R13  ) */
/*                M-K ( 0     0   R22  R23  ) */
/*              K+L-M ( 0     0    0   R33  ) */

/*  where */

/*    C = diag( ALPHA(K+1), ... , ALPHA(M) ), */
/*    S = diag( BETA(K+1),  ... , BETA(M) ), */
/*    C**2 + S**2 = I. */

/*    (R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N), and R33 is stored */
/*    ( 0  R22 R23 ) */
/*    in B(M-K+1:L,N+M-K-L+1:N) on exit. */

/*  The routine computes C, S, R, and optionally the orthogonal */
/*  transformation matrices U, V and Q. */

/*  In particular, if B is an N-by-N nonsingular matrix, then the GSVD of */
/*  A and B implicitly gives the SVD of A*inv(B): */
/*                       A*inv(B) = U*(D1*inv(D2))*V'. */
/*  If ( A',B')' has orthonormal columns, then the GSVD of A and B is */
/*  also equal to the CS decomposition of A and B. Furthermore, the GSVD */
/*  can be used to derive the solution of the eigenvalue problem: */
/*                       A'*A x = lambda* B'*B x. */
/*  In some literature, the GSVD of A and B is presented in the form */
/*                   U'*A*X = ( 0 D1 ),   V'*B*X = ( 0 D2 ) */
/*  where U and V are orthogonal and X is nonsingular, D1 and D2 are */
/*  ``diagonal''.  The former GSVD form can be converted to the latter */
/*  form by taking the nonsingular matrix X as */

/*                       X = Q*( I   0    ) */
/*                             ( 0 inv(R) ). */

/*  Arguments */
/*  ========= */

/*  JOBU    (input) CHARACTER*1 */
/*          = 'U':  Orthogonal matrix U is computed; */
/*          = 'N':  U is not computed. */

/*  JOBV    (input) CHARACTER*1 */
/*          = 'V':  Orthogonal matrix V is computed; */
/*          = 'N':  V is not computed. */

/*  JOBQ    (input) CHARACTER*1 */
/*          = 'Q':  Orthogonal matrix Q is computed; */
/*          = 'N':  Q is not computed. */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix A.  M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrices A and B.  N >= 0. */

/*  P       (input) INTEGER */
/*          The number of rows of the matrix B.  P >= 0. */

/*  K       (output) INTEGER */
/*  L       (output) INTEGER */
/*          On exit, K and L specify the dimension of the subblocks */
/*          described in the Purpose section. */
/*          K + L = effective numerical rank of (A',B')'. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the M-by-N matrix A. */
/*          On exit, A contains the triangular matrix R, or part of R. */
/*          See Purpose for details. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A. LDA >= MAX(1,M). */

/*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*          On entry, the P-by-N matrix B. */
/*          On exit, B contains the triangular matrix R if M-K-L < 0. */
/*          See Purpose for details. */

/*  LDB     (input) INTEGER */
/*          The leading dimension of the array B. LDB >= MAX(1,P). */

/*  ALPHA   (output) DOUBLE PRECISION array, dimension (N) */
/*  BETA    (output) DOUBLE PRECISION array, dimension (N) */
/*          On exit, ALPHA and BETA contain the generalized singular */
/*          value pairs of A and B; */
/*            ALPHA(1:K) = 1, */
/*            BETA(1:K)  = 0, */
/*          and if M-K-L >= 0, */
/*            ALPHA(K+1:K+L) = C, */
/*            BETA(K+1:K+L)  = S, */
/*          or if M-K-L < 0, */
/*            ALPHA(K+1:M)=C, ALPHA(M+1:K+L)=0 */
/*            BETA(K+1:M) =S, BETA(M+1:K+L) =1 */
/*          and */
/*            ALPHA(K+L+1:N) = 0 */
/*            BETA(K+L+1:N)  = 0 */

/*  U       (output) DOUBLE PRECISION array, dimension (LDU,M) */
/*          If JOBU = 'U', U contains the M-by-M orthogonal matrix U. */
/*          If JOBU = 'N', U is not referenced. */

/*  LDU     (input) INTEGER */
/*          The leading dimension of the array U. LDU >= MAX(1,M) if */
/*          JOBU = 'U'; LDU >= 1 otherwise. */

/*  V       (output) DOUBLE PRECISION array, dimension (LDV,P) */
/*          If JOBV = 'V', V contains the P-by-P orthogonal matrix V. */
/*          If JOBV = 'N', V is not referenced. */

/*  LDV     (input) INTEGER */
/*          The leading dimension of the array V. LDV >= MAX(1,P) if */
/*          JOBV = 'V'; LDV >= 1 otherwise. */

/*  Q       (output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*          If JOBQ = 'Q', Q contains the N-by-N orthogonal matrix Q. */
/*          If JOBQ = 'N', Q is not referenced. */

/*  LDQ     (input) INTEGER */
/*          The leading dimension of the array Q. LDQ >= MAX(1,N) if */
/*          JOBQ = 'Q'; LDQ >= 1 otherwise. */

/*  WORK    (workspace) DOUBLE PRECISION array, */
/*                      dimension (MAX(3*N,M,P)+N) */

/*  IWORK   (workspace/output) INTEGER array, dimension (N) */
/*          On exit, IWORK stores the sorting information. More */
/*          precisely, the following loop will sort ALPHA */
/*             for I = K+1, MIN(M,K+L) */
/*                 swap ALPHA(I) and ALPHA(IWORK(I)) */
/*             endfor */
/*          such that ALPHA(1) >= ALPHA(2) >= ... >= ALPHA(N). */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/*          > 0:  if INFO = 1, the Jacobi-type procedure failed to */
/*                converge.  For further details, see subroutine DTGSJA. */

/*  Internal Parameters */
/*  =================== */

/*  TOLA    DOUBLE PRECISION */
/*  TOLB    DOUBLE PRECISION */
/*          TOLA and TOLB are the thresholds to determine the effective */
/*          rank of (A',B')'. Generally, they are set to */
/*                   TOLA = MAX(M,N)*norm(A)*MAZHEPS, */
/*                   TOLB = MAX(P,N)*norm(B)*MAZHEPS. */
/*          The size of TOLA and TOLB may affect the size of backward */
/*          errors of the decomposition. */

/*  Further Details */
/*  =============== */

/*  2-96 Based on modifications by */
/*     Ming Gu and Huan Ren, Computer Science Division, University of */
/*     California at Berkeley, USA */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --alpha;
    --beta;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --work;
    --iwork;

    /* Function Body */
    wantu = lsame_(jobu, "U");
    wantv = lsame_(jobv, "V");
    wantq = lsame_(jobq, "Q");

    *info = 0;
    if (! (wantu || lsame_(jobu, "N"))) {
	*info = -1;
    } else if (! (wantv || lsame_(jobv, "N"))) {
	*info = -2;
    } else if (! (wantq || lsame_(jobq, "N"))) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*p < 0) {
	*info = -6;
    } else if (*lda < MAX(1,*m)) {
	*info = -10;
    } else if (*ldb < MAX(1,*p)) {
	*info = -12;
    } else if (*ldu < 1 || wantu && *ldu < *m) {
	*info = -16;
    } else if (*ldv < 1 || wantv && *ldv < *p) {
	*info = -18;
    } else if (*ldq < 1 || wantq && *ldq < *n) {
	*info = -20;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGGSVD", &i__1);
	return 0;
    }

/*     Compute the Frobenius norm of matrices A and B */

    anorm = dlange_("1", m, n, &a[a_offset], lda, &work[1]);
    bnorm = dlange_("1", p, n, &b[b_offset], ldb, &work[1]);

/*     Get machine precision and set up threshold for determining */
/*     the effective numerical rank of the matrices A and B. */

    ulp = dlamch_("Precision");
    unfl = dlamch_("Safe Minimum");
    tola = MAX(*m,*n) * MAX(anorm,unfl) * ulp;
    tolb = MAX(*p,*n) * MAX(bnorm,unfl) * ulp;

/*     Preprocessing */

    dggsvp_(jobu, jobv, jobq, m, p, n, &a[a_offset], lda, &b[b_offset], ldb, &
	    tola, &tolb, k, l, &u[u_offset], ldu, &v[v_offset], ldv, &q[
	    q_offset], ldq, &iwork[1], &work[1], &work[*n + 1], info);

/*     Compute the GSVD of two upper "triangular" matrices */

    dtgsja_(jobu, jobv, jobq, m, p, n, k, l, &a[a_offset], lda, &b[b_offset], 
	    ldb, &tola, &tolb, &alpha[1], &beta[1], &u[u_offset], ldu, &v[
	    v_offset], ldv, &q[q_offset], ldq, &work[1], &ncycle, info);

/*     Sort the singular values and store the pivot indices in IWORK */
/*     Copy ALPHA to WORK, then sort ALPHA in WORK */

    dcopy_(n, &alpha[1], &c__1, &work[1], &c__1);
/* Computing MIN */
    i__1 = *l, i__2 = *m - *k;
    ibnd = MIN(i__1,i__2);
    i__1 = ibnd;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Scan for largest ALPHA(K+I) */

	isub = i__;
	smax = work[*k + i__];
	i__2 = ibnd;
	for (j = i__ + 1; j <= i__2; ++j) {
	    temp = work[*k + j];
	    if (temp > smax) {
		isub = j;
		smax = temp;
	    }
/* L10: */
	}
	if (isub != i__) {
	    work[*k + isub] = work[*k + i__];
	    work[*k + i__] = smax;
	    iwork[*k + i__] = *k + isub;
	} else {
	    iwork[*k + i__] = *k + i__;
	}
/* L20: */
    }

    return 0;

/*     End of DGGSVD */

} /* dggsvd_ */
