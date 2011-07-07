/* dggrqf.f -- translated by f2c (version 20061008).
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

 int dggrqf_(int *m, int *p, int *n, double *
	a, int *lda, double *taua, double *b, int *ldb, 
	double *taub, double *work, int *lwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    int nb, nb1, nb2, nb3, lopt;
    extern  int dgeqrf_(int *, int *, double *, 
	    int *, double *, double *, int *, int *), 
	    dgerqf_(int *, int *, double *, int *, double 
	    *, double *, int *, int *), xerbla_(char *, int *);
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *);
    extern  int dormrq_(char *, char *, int *, int *, 
	    int *, double *, int *, double *, double *, 
	    int *, double *, int *, int *);
    int lwkopt;
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

/*  DGGRQF computes a generalized RQ factorization of an M-by-N matrix A */
/*  and a P-by-N matrix B: */

/*              A = R*Q,        B = Z*T*Q, */

/*  where Q is an N-by-N orthogonal matrix, Z is a P-by-P orthogonal */
/*  matrix, and R and T assume one of the forms: */

/*  if M <= N,  R = ( 0  R12 ) M,   or if M > N,  R = ( R11 ) M-N, */
/*                   N-M  M                           ( R21 ) N */
/*                                                       N */

/*  where R12 or R21 is upper triangular, and */

/*  if P >= N,  T = ( T11 ) N  ,   or if P < N,  T = ( T11  T12 ) P, */
/*                  (  0  ) P-N                         P   N-P */
/*                     N */

/*  where T11 is upper triangular. */

/*  In particular, if B is square and nonsingular, the GRQ factorization */
/*  of A and B implicitly gives the RQ factorization of A*inv(B): */

/*               A*inv(B) = (R*inv(T))*Z' */

/*  where inv(B) denotes the inverse of the matrix B, and Z' denotes the */
/*  transpose of the matrix Z. */

/*  Arguments */
/*  ========= */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix A.  M >= 0. */

/*  P       (input) INTEGER */
/*          The number of rows of the matrix B.  P >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrices A and B. N >= 0. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the M-by-N matrix A. */
/*          On exit, if M <= N, the upper triangle of the subarray */
/*          A(1:M,N-M+1:N) contains the M-by-M upper triangular matrix R; */
/*          if M > N, the elements on and above the (M-N)-th subdiagonal */
/*          contain the M-by-N upper trapezoidal matrix R; the remaining */
/*          elements, with the array TAUA, represent the orthogonal */
/*          matrix Q as a product of elementary reflectors (see Further */
/*          Details). */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A. LDA >= MAX(1,M). */

/*  TAUA    (output) DOUBLE PRECISION array, dimension (MIN(M,N)) */
/*          The scalar factors of the elementary reflectors which */
/*          represent the orthogonal matrix Q (see Further Details). */

/*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*          On entry, the P-by-N matrix B. */
/*          On exit, the elements on and above the diagonal of the array */
/*          contain the MIN(P,N)-by-N upper trapezoidal matrix T (T is */
/*          upper triangular if P >= N); the elements below the diagonal, */
/*          with the array TAUB, represent the orthogonal matrix Z as a */
/*          product of elementary reflectors (see Further Details). */

/*  LDB     (input) INTEGER */
/*          The leading dimension of the array B. LDB >= MAX(1,P). */

/*  TAUB    (output) DOUBLE PRECISION array, dimension (MIN(P,N)) */
/*          The scalar factors of the elementary reflectors which */
/*          represent the orthogonal matrix Z (see Further Details). */

/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

/*  LWORK   (input) INTEGER */
/*          The dimension of the array WORK. LWORK >= MAX(1,N,M,P). */
/*          For optimum performance LWORK >= MAX(N,M,P)*MAX(NB1,NB2,NB3), */
/*          where NB1 is the optimal blocksize for the RQ factorization */
/*          of an M-by-N matrix, NB2 is the optimal blocksize for the */
/*          QR factorization of a P-by-N matrix, and NB3 is the optimal */
/*          blocksize for a call of DORMRQ. */

/*          If LWORK = -1, then a workspace query is assumed; the routine */
/*          only calculates the optimal size of the WORK array, returns */
/*          this value as the first entry of the WORK array, and no error */
/*          message related to LWORK is issued by XERBLA. */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INF0= -i, the i-th argument had an illegal value. */

/*  Further Details */
/*  =============== */

/*  The matrix Q is represented as a product of elementary reflectors */

/*     Q = H(1) H(2) . . . H(k), where k = MIN(m,n). */

/*  Each H(i) has the form */

/*     H(i) = I - taua * v * v' */

/*  where taua is a float scalar, and v is a float vector with */
/*  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit in */
/*  A(m-k+i,1:n-k+i-1), and taua in TAUA(i). */
/*  To form Q explicitly, use LAPACK subroutine DORGRQ. */
/*  To use Q to update another matrix, use LAPACK subroutine DORMRQ. */

/*  The matrix Z is represented as a product of elementary reflectors */

/*     Z = H(1) H(2) . . . H(k), where k = MIN(p,n). */

/*  Each H(i) has the form */

/*     H(i) = I - taub * v * v' */

/*  where taub is a float scalar, and v is a float vector with */
/*  v(1:i-1) = 0 and v(i) = 1; v(i+1:p) is stored on exit in B(i+1:p,i), */
/*  and taub in TAUB(i). */
/*  To form Z explicitly, use LAPACK subroutine DORGQR. */
/*  To use Z to update another matrix, use LAPACK subroutine DORMQR. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --taua;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --taub;
    --work;

    /* Function Body */
    *info = 0;
    nb1 = ilaenv_(&c__1, "DGERQF", " ", m, n, &c_n1, &c_n1);
    nb2 = ilaenv_(&c__1, "DGEQRF", " ", p, n, &c_n1, &c_n1);
    nb3 = ilaenv_(&c__1, "DORMRQ", " ", m, n, p, &c_n1);
/* Computing MAX */
    i__1 = MAX(nb1,nb2);
    nb = MAX(i__1,nb3);
/* Computing MAX */
    i__1 = MAX(*n,*m);
    lwkopt = MAX(i__1,*p) * nb;
    work[1] = (double) lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*p < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < MAX(1,*m)) {
	*info = -5;
    } else if (*ldb < MAX(1,*p)) {
	*info = -8;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = MAX(1,*m), i__1 = MAX(i__1,*p);
	if (*lwork < MAX(i__1,*n) && ! lquery) {
	    *info = -11;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGGRQF", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     RQ factorization of M-by-N matrix A: A = R*Q */

    dgerqf_(m, n, &a[a_offset], lda, &taua[1], &work[1], lwork, info);
    lopt = (int) work[1];

/*     Update B := B*Q' */

    i__1 = MIN(*m,*n);
/* Computing MAX */
    i__2 = 1, i__3 = *m - *n + 1;
    dormrq_("Right", "Transpose", p, n, &i__1, &a[MAX(i__2, i__3)+ a_dim1], 
	    lda, &taua[1], &b[b_offset], ldb, &work[1], lwork, info);
/* Computing MAX */
    i__1 = lopt, i__2 = (int) work[1];
    lopt = MAX(i__1,i__2);

/*     QR factorization of P-by-N matrix B: B = Z*T */

    dgeqrf_(p, n, &b[b_offset], ldb, &taub[1], &work[1], lwork, info);
/* Computing MAX */
    i__1 = lopt, i__2 = (int) work[1];
    work[1] = (double) MAX(i__1,i__2);

    return 0;

/*     End of DGGRQF */

} /* dggrqf_ */
