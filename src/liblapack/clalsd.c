/* clalsd.f -- translated by f2c (version 20061008).
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
static int c__1 = 1;
static int c__0 = 0;
static float c_b10 = 1.f;
static float c_b35 = 0.f;

 int clalsd_(char *uplo, int *smlsiz, int *n, int 
	*nrhs, float *d__, float *e, complex *b, int *ldb, float *rcond, 
	int *rank, complex *work, float *rwork, int *iwork, int *
	info)
{
    /* System generated locals */
    int b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    float r__1;
    complex q__1;

    /* Builtin functions */
    double r_imag(complex *), log(double), r_sign(float *, float *);

    /* Local variables */
    int c__, i__, j, k;
    float r__;
    int s, u, z__;
    float cs;
    int bx;
    float sn;
    int st, vt, nm1, st1;
    float eps;
    int iwk;
    float tol;
    int difl, difr;
    float rcnd;
    int jcol, irwb, perm, nsub, nlvl, sqre, bxst, jrow, irwu, jimag, 
	    jfloat;
    extern  int sgemm_(char *, char *, int *, int *, 
	    int *, float *, float *, int *, float *, int *, float *, 
	    float *, int *);
    int irwib;
    extern  int ccopy_(int *, complex *, int *, 
	    complex *, int *);
    int poles, sizei, irwrb, nsize;
    extern  int csrot_(int *, complex *, int *, 
	    complex *, int *, float *, float *);
    int irwvt, icmpq1, icmpq2;
    extern  int clalsa_(int *, int *, int *, 
	    int *, complex *, int *, complex *, int *, float *, 
	    int *, float *, int *, float *, float *, float *, float *, 
	    int *, int *, int *, int *, float *, float *, float *
, float *, int *, int *), clascl_(char *, int *, 
	    int *, float *, float *, int *, int *, complex *, 
	    int *, int *);
    extern double slamch_(char *);
    extern  int slasda_(int *, int *, int *, 
	    int *, float *, float *, float *, int *, float *, int *, 
	    float *, float *, float *, float *, int *, int *, int *, 
	    int *, float *, float *, float *, float *, int *, int *), 
	    clacpy_(char *, int *, int *, complex *, int *, 
	    complex *, int *), claset_(char *, int *, int 
	    *, complex *, complex *, complex *, int *), xerbla_(
	    char *, int *), slascl_(char *, int *, int *, 
	    float *, float *, int *, int *, float *, int *, int *
);
    extern int isamax_(int *, float *, int *);
    int givcol;
    extern  int slasdq_(char *, int *, int *, int 
	    *, int *, int *, float *, float *, float *, int *, float *
, int *, float *, int *, float *, int *), 
	    slaset_(char *, int *, int *, float *, float *, float *, 
	    int *), slartg_(float *, float *, float *, float *, float *
);
    float orgnrm;
    int givnum;
    extern double slanst_(char *, int *, float *, float *);
    extern  int slasrt_(char *, int *, float *, int *);
    int givptr, nrwork, irwwrk, smlszp;


/*  -- LAPACK routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  CLALSD uses the singular value decomposition of A to solve the least */
/*  squares problem of finding X to minimize the Euclidean norm of each */
/*  column of A*X-B, where A is N-by-N upper bidiagonal, and X and B */
/*  are N-by-NRHS. The solution X overwrites B. */

/*  The singular values of A smaller than RCOND times the largest */
/*  singular value are treated as zero in solving the least squares */
/*  problem; in this case a minimum norm solution is returned. */
/*  The actual singular values are returned in D in ascending order. */

/*  This code makes very mild assumptions about floating point */
/*  arithmetic. It will work on machines with a guard digit in */
/*  add/subtract, or on those binary machines without guard digits */
/*  which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2. */
/*  It could conceivably fail on hexadecimal or decimal machines */
/*  without guard digits, but we know of none. */

/*  Arguments */
/*  ========= */

/*  UPLO   (input) CHARACTER*1 */
/*         = 'U': D and E define an upper bidiagonal matrix. */
/*         = 'L': D and E define a  lower bidiagonal matrix. */

/*  SMLSIZ (input) INTEGER */
/*         The maximum size of the subproblems at the bottom of the */
/*         computation tree. */

/*  N      (input) INTEGER */
/*         The dimension of the  bidiagonal matrix.  N >= 0. */

/*  NRHS   (input) INTEGER */
/*         The number of columns of B. NRHS must be at least 1. */

/*  D      (input/output) REAL array, dimension (N) */
/*         On entry D contains the main diagonal of the bidiagonal */
/*         matrix. On exit, if INFO = 0, D contains its singular values. */

/*  E      (input/output) REAL array, dimension (N-1) */
/*         Contains the super-diagonal entries of the bidiagonal matrix. */
/*         On exit, E has been destroyed. */

/*  B      (input/output) COMPLEX array, dimension (LDB,NRHS) */
/*         On input, B contains the right hand sides of the least */
/*         squares problem. On output, B contains the solution X. */

/*  LDB    (input) INTEGER */
/*         The leading dimension of B in the calling subprogram. */
/*         LDB must be at least MAX(1,N). */

/*  RCOND  (input) REAL */
/*         The singular values of A less than or equal to RCOND times */
/*         the largest singular value are treated as zero in solving */
/*         the least squares problem. If RCOND is negative, */
/*         machine precision is used instead. */
/*         For example, if diag(S)*X=B were the least squares problem, */
/*         where diag(S) is a diagonal matrix of singular values, the */
/*         solution would be X(i) = B(i) / S(i) if S(i) is greater than */
/*         RCOND*MAX(S), and X(i) = 0 if S(i) is less than or equal to */
/*         RCOND*MAX(S). */

/*  RANK   (output) INTEGER */
/*         The number of singular values of A greater than RCOND times */
/*         the largest singular value. */

/*  WORK   (workspace) COMPLEX array, dimension (N * NRHS). */

/*  RWORK  (workspace) REAL array, dimension at least */
/*         (9*N + 2*N*SMLSIZ + 8*N*NLVL + 3*SMLSIZ*NRHS + (SMLSIZ+1)**2), */
/*         where */
/*         NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 ) */

/*  IWORK  (workspace) INTEGER array, dimension (3*N*NLVL + 11*N). */

/*  INFO   (output) INTEGER */
/*         = 0:  successful exit. */
/*         < 0:  if INFO = -i, the i-th argument had an illegal value. */
/*         > 0:  The algorithm failed to compute an singular value while */
/*               working on the submatrix lying in rows and columns */
/*               INFO/(N+1) through MOD(INFO,N+1). */

/*  Further Details */
/*  =============== */

/*  Based on contributions by */
/*     Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/*       California at Berkeley, USA */
/*     Osni Marques, LBNL/NERSC, USA */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --d__;
    --e;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;
    --rwork;
    --iwork;

    /* Function Body */
    *info = 0;

    if (*n < 0) {
	*info = -3;
    } else if (*nrhs < 1) {
	*info = -4;
    } else if (*ldb < 1 || *ldb < *n) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CLALSD", &i__1);
	return 0;
    }

    eps = slamch_("Epsilon");

/*     Set up the tolerance. */

    if (*rcond <= 0.f || *rcond >= 1.f) {
	rcnd = eps;
    } else {
	rcnd = *rcond;
    }

    *rank = 0;

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    } else if (*n == 1) {
	if (d__[1] == 0.f) {
	    claset_("A", &c__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb);
	} else {
	    *rank = 1;
	    clascl_("G", &c__0, &c__0, &d__[1], &c_b10, &c__1, nrhs, &b[
		    b_offset], ldb, info);
	    d__[1] = ABS(d__[1]);
	}
	return 0;
    }

/*     Rotate the matrix if it is lower bidiagonal. */

    if (*(unsigned char *)uplo == 'L') {
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    slartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
	    d__[i__] = r__;
	    e[i__] = sn * d__[i__ + 1];
	    d__[i__ + 1] = cs * d__[i__ + 1];
	    if (*nrhs == 1) {
		csrot_(&c__1, &b[i__ + b_dim1], &c__1, &b[i__ + 1 + b_dim1], &
			c__1, &cs, &sn);
	    } else {
		rwork[(i__ << 1) - 1] = cs;
		rwork[i__ * 2] = sn;
	    }
/* L10: */
	}
	if (*nrhs > 1) {
	    i__1 = *nrhs;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = *n - 1;
		for (j = 1; j <= i__2; ++j) {
		    cs = rwork[(j << 1) - 1];
		    sn = rwork[j * 2];
		    csrot_(&c__1, &b[j + i__ * b_dim1], &c__1, &b[j + 1 + i__ 
			    * b_dim1], &c__1, &cs, &sn);
/* L20: */
		}
/* L30: */
	    }
	}
    }

/*     Scale. */

    nm1 = *n - 1;
    orgnrm = slanst_("M", n, &d__[1], &e[1]);
    if (orgnrm == 0.f) {
	claset_("A", n, nrhs, &c_b1, &c_b1, &b[b_offset], ldb);
	return 0;
    }

    slascl_("G", &c__0, &c__0, &orgnrm, &c_b10, n, &c__1, &d__[1], n, info);
    slascl_("G", &c__0, &c__0, &orgnrm, &c_b10, &nm1, &c__1, &e[1], &nm1, 
	    info);

/*     If N is smaller than the minimum divide size SMLSIZ, then solve */
/*     the problem with another solver. */

    if (*n <= *smlsiz) {
	irwu = 1;
	irwvt = irwu + *n * *n;
	irwwrk = irwvt + *n * *n;
	irwrb = irwwrk;
	irwib = irwrb + *n * *nrhs;
	irwb = irwib + *n * *nrhs;
	slaset_("A", n, n, &c_b35, &c_b10, &rwork[irwu], n);
	slaset_("A", n, n, &c_b35, &c_b10, &rwork[irwvt], n);
	slasdq_("U", &c__0, n, n, n, &c__0, &d__[1], &e[1], &rwork[irwvt], n, 
		&rwork[irwu], n, &rwork[irwwrk], &c__1, &rwork[irwwrk], info);
	if (*info != 0) {
	    return 0;
	}

/*        In the float version, B is passed to SLASDQ and multiplied */
/*        internally by Q'. Here B is complex and that product is */
/*        computed below in two steps (float and imaginary parts). */

	j = irwb - 1;
	i__1 = *nrhs;
	for (jcol = 1; jcol <= i__1; ++jcol) {
	    i__2 = *n;
	    for (jrow = 1; jrow <= i__2; ++jrow) {
		++j;
		i__3 = jrow + jcol * b_dim1;
		rwork[j] = b[i__3].r;
/* L40: */
	    }
/* L50: */
	}
	sgemm_("T", "N", n, nrhs, n, &c_b10, &rwork[irwu], n, &rwork[irwb], n, 
		 &c_b35, &rwork[irwrb], n);
	j = irwb - 1;
	i__1 = *nrhs;
	for (jcol = 1; jcol <= i__1; ++jcol) {
	    i__2 = *n;
	    for (jrow = 1; jrow <= i__2; ++jrow) {
		++j;
		rwork[j] = r_imag(&b[jrow + jcol * b_dim1]);
/* L60: */
	    }
/* L70: */
	}
	sgemm_("T", "N", n, nrhs, n, &c_b10, &rwork[irwu], n, &rwork[irwb], n, 
		 &c_b35, &rwork[irwib], n);
	jfloat = irwrb - 1;
	jimag = irwib - 1;
	i__1 = *nrhs;
	for (jcol = 1; jcol <= i__1; ++jcol) {
	    i__2 = *n;
	    for (jrow = 1; jrow <= i__2; ++jrow) {
		++jfloat;
		++jimag;
		i__3 = jrow + jcol * b_dim1;
		i__4 = jfloat;
		i__5 = jimag;
		q__1.r = rwork[i__4], q__1.i = rwork[i__5];
		b[i__3].r = q__1.r, b[i__3].i = q__1.i;
/* L80: */
	    }
/* L90: */
	}

	tol = rcnd * (r__1 = d__[isamax_(n, &d__[1], &c__1)], ABS(r__1));
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (d__[i__] <= tol) {
		claset_("A", &c__1, nrhs, &c_b1, &c_b1, &b[i__ + b_dim1], ldb);
	    } else {
		clascl_("G", &c__0, &c__0, &d__[i__], &c_b10, &c__1, nrhs, &b[
			i__ + b_dim1], ldb, info);
		++(*rank);
	    }
/* L100: */
	}

/*        Since B is complex, the following call to SGEMM is performed */
/*        in two steps (float and imaginary parts). That is for V * B */
/*        (in the float version of the code V' is stored in WORK). */

/*        CALL SGEMM( 'T', 'N', N, NRHS, N, ONE, WORK, N, B, LDB, ZERO, */
/*    $               WORK( NWORK ), N ) */

	j = irwb - 1;
	i__1 = *nrhs;
	for (jcol = 1; jcol <= i__1; ++jcol) {
	    i__2 = *n;
	    for (jrow = 1; jrow <= i__2; ++jrow) {
		++j;
		i__3 = jrow + jcol * b_dim1;
		rwork[j] = b[i__3].r;
/* L110: */
	    }
/* L120: */
	}
	sgemm_("T", "N", n, nrhs, n, &c_b10, &rwork[irwvt], n, &rwork[irwb], 
		n, &c_b35, &rwork[irwrb], n);
	j = irwb - 1;
	i__1 = *nrhs;
	for (jcol = 1; jcol <= i__1; ++jcol) {
	    i__2 = *n;
	    for (jrow = 1; jrow <= i__2; ++jrow) {
		++j;
		rwork[j] = r_imag(&b[jrow + jcol * b_dim1]);
/* L130: */
	    }
/* L140: */
	}
	sgemm_("T", "N", n, nrhs, n, &c_b10, &rwork[irwvt], n, &rwork[irwb], 
		n, &c_b35, &rwork[irwib], n);
	jfloat = irwrb - 1;
	jimag = irwib - 1;
	i__1 = *nrhs;
	for (jcol = 1; jcol <= i__1; ++jcol) {
	    i__2 = *n;
	    for (jrow = 1; jrow <= i__2; ++jrow) {
		++jfloat;
		++jimag;
		i__3 = jrow + jcol * b_dim1;
		i__4 = jfloat;
		i__5 = jimag;
		q__1.r = rwork[i__4], q__1.i = rwork[i__5];
		b[i__3].r = q__1.r, b[i__3].i = q__1.i;
/* L150: */
	    }
/* L160: */
	}

/*        Unscale. */

	slascl_("G", &c__0, &c__0, &c_b10, &orgnrm, n, &c__1, &d__[1], n, 
		info);
	slasrt_("D", n, &d__[1], info);
	clascl_("G", &c__0, &c__0, &orgnrm, &c_b10, n, nrhs, &b[b_offset], 
		ldb, info);

	return 0;
    }

/*     Book-keeping and setting up some constants. */

    nlvl = (int) (log((float) (*n) / (float) (*smlsiz + 1)) / log(2.f)) + 1;

    smlszp = *smlsiz + 1;

    u = 1;
    vt = *smlsiz * *n + 1;
    difl = vt + smlszp * *n;
    difr = difl + nlvl * *n;
    z__ = difr + (nlvl * *n << 1);
    c__ = z__ + nlvl * *n;
    s = c__ + *n;
    poles = s + *n;
    givnum = poles + (nlvl << 1) * *n;
    nrwork = givnum + (nlvl << 1) * *n;
    bx = 1;

    irwrb = nrwork;
    irwib = irwrb + *smlsiz * *nrhs;
    irwb = irwib + *smlsiz * *nrhs;

    sizei = *n + 1;
    k = sizei + *n;
    givptr = k + *n;
    perm = givptr + *n;
    givcol = perm + nlvl * *n;
    iwk = givcol + (nlvl * *n << 1);

    st = 1;
    sqre = 0;
    icmpq1 = 1;
    icmpq2 = 0;
    nsub = 0;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((r__1 = d__[i__], ABS(r__1)) < eps) {
	    d__[i__] = r_sign(&eps, &d__[i__]);
	}
/* L170: */
    }

    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((r__1 = e[i__], ABS(r__1)) < eps || i__ == nm1) {
	    ++nsub;
	    iwork[nsub] = st;

/*           Subproblem found. First determine its size and then */
/*           apply divide and conquer on it. */

	    if (i__ < nm1) {

/*              A subproblem with E(I) small for I < NM1. */

		nsize = i__ - st + 1;
		iwork[sizei + nsub - 1] = nsize;
	    } else if ((r__1 = e[i__], ABS(r__1)) >= eps) {

/*              A subproblem with E(NM1) not too small but I = NM1. */

		nsize = *n - st + 1;
		iwork[sizei + nsub - 1] = nsize;
	    } else {

/*              A subproblem with E(NM1) small. This implies an */
/*              1-by-1 subproblem at D(N), which is not solved */
/*              explicitly. */

		nsize = i__ - st + 1;
		iwork[sizei + nsub - 1] = nsize;
		++nsub;
		iwork[nsub] = *n;
		iwork[sizei + nsub - 1] = 1;
		ccopy_(nrhs, &b[*n + b_dim1], ldb, &work[bx + nm1], n);
	    }
	    st1 = st - 1;
	    if (nsize == 1) {

/*              This is a 1-by-1 subproblem and is not solved */
/*              explicitly. */

		ccopy_(nrhs, &b[st + b_dim1], ldb, &work[bx + st1], n);
	    } else if (nsize <= *smlsiz) {

/*              This is a small subproblem and is solved by SLASDQ. */

		slaset_("A", &nsize, &nsize, &c_b35, &c_b10, &rwork[vt + st1], 
			 n);
		slaset_("A", &nsize, &nsize, &c_b35, &c_b10, &rwork[u + st1], 
			n);
		slasdq_("U", &c__0, &nsize, &nsize, &nsize, &c__0, &d__[st], &
			e[st], &rwork[vt + st1], n, &rwork[u + st1], n, &
			rwork[nrwork], &c__1, &rwork[nrwork], info)
			;
		if (*info != 0) {
		    return 0;
		}

/*              In the float version, B is passed to SLASDQ and multiplied */
/*              internally by Q'. Here B is complex and that product is */
/*              computed below in two steps (float and imaginary parts). */

		j = irwb - 1;
		i__2 = *nrhs;
		for (jcol = 1; jcol <= i__2; ++jcol) {
		    i__3 = st + nsize - 1;
		    for (jrow = st; jrow <= i__3; ++jrow) {
			++j;
			i__4 = jrow + jcol * b_dim1;
			rwork[j] = b[i__4].r;
/* L180: */
		    }
/* L190: */
		}
		sgemm_("T", "N", &nsize, nrhs, &nsize, &c_b10, &rwork[u + st1]
, n, &rwork[irwb], &nsize, &c_b35, &rwork[irwrb], &
			nsize);
		j = irwb - 1;
		i__2 = *nrhs;
		for (jcol = 1; jcol <= i__2; ++jcol) {
		    i__3 = st + nsize - 1;
		    for (jrow = st; jrow <= i__3; ++jrow) {
			++j;
			rwork[j] = r_imag(&b[jrow + jcol * b_dim1]);
/* L200: */
		    }
/* L210: */
		}
		sgemm_("T", "N", &nsize, nrhs, &nsize, &c_b10, &rwork[u + st1]
, n, &rwork[irwb], &nsize, &c_b35, &rwork[irwib], &
			nsize);
		jfloat = irwrb - 1;
		jimag = irwib - 1;
		i__2 = *nrhs;
		for (jcol = 1; jcol <= i__2; ++jcol) {
		    i__3 = st + nsize - 1;
		    for (jrow = st; jrow <= i__3; ++jrow) {
			++jfloat;
			++jimag;
			i__4 = jrow + jcol * b_dim1;
			i__5 = jfloat;
			i__6 = jimag;
			q__1.r = rwork[i__5], q__1.i = rwork[i__6];
			b[i__4].r = q__1.r, b[i__4].i = q__1.i;
/* L220: */
		    }
/* L230: */
		}

		clacpy_("A", &nsize, nrhs, &b[st + b_dim1], ldb, &work[bx + 
			st1], n);
	    } else {

/*              A large problem. Solve it using divide and conquer. */

		slasda_(&icmpq1, smlsiz, &nsize, &sqre, &d__[st], &e[st], &
			rwork[u + st1], n, &rwork[vt + st1], &iwork[k + st1], 
			&rwork[difl + st1], &rwork[difr + st1], &rwork[z__ + 
			st1], &rwork[poles + st1], &iwork[givptr + st1], &
			iwork[givcol + st1], n, &iwork[perm + st1], &rwork[
			givnum + st1], &rwork[c__ + st1], &rwork[s + st1], &
			rwork[nrwork], &iwork[iwk], info);
		if (*info != 0) {
		    return 0;
		}
		bxst = bx + st1;
		clalsa_(&icmpq2, smlsiz, &nsize, nrhs, &b[st + b_dim1], ldb, &
			work[bxst], n, &rwork[u + st1], n, &rwork[vt + st1], &
			iwork[k + st1], &rwork[difl + st1], &rwork[difr + st1]
, &rwork[z__ + st1], &rwork[poles + st1], &iwork[
			givptr + st1], &iwork[givcol + st1], n, &iwork[perm + 
			st1], &rwork[givnum + st1], &rwork[c__ + st1], &rwork[
			s + st1], &rwork[nrwork], &iwork[iwk], info);
		if (*info != 0) {
		    return 0;
		}
	    }
	    st = i__ + 1;
	}
/* L240: */
    }

/*     Apply the singular values and treat the tiny ones as zero. */

    tol = rcnd * (r__1 = d__[isamax_(n, &d__[1], &c__1)], ABS(r__1));

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Some of the elements in D can be negative because 1-by-1 */
/*        subproblems were not solved explicitly. */

	if ((r__1 = d__[i__], ABS(r__1)) <= tol) {
	    claset_("A", &c__1, nrhs, &c_b1, &c_b1, &work[bx + i__ - 1], n);
	} else {
	    ++(*rank);
	    clascl_("G", &c__0, &c__0, &d__[i__], &c_b10, &c__1, nrhs, &work[
		    bx + i__ - 1], n, info);
	}
	d__[i__] = (r__1 = d__[i__], ABS(r__1));
/* L250: */
    }

/*     Now apply back the right singular vectors. */

    icmpq2 = 1;
    i__1 = nsub;
    for (i__ = 1; i__ <= i__1; ++i__) {
	st = iwork[i__];
	st1 = st - 1;
	nsize = iwork[sizei + i__ - 1];
	bxst = bx + st1;
	if (nsize == 1) {
	    ccopy_(nrhs, &work[bxst], n, &b[st + b_dim1], ldb);
	} else if (nsize <= *smlsiz) {

/*           Since B and BX are complex, the following call to SGEMM */
/*           is performed in two steps (float and imaginary parts). */

/*           CALL SGEMM( 'T', 'N', NSIZE, NRHS, NSIZE, ONE, */
/*    $                  RWORK( VT+ST1 ), N, RWORK( BXST ), N, ZERO, */
/*    $                  B( ST, 1 ), LDB ) */

	    j = bxst - *n - 1;
	    jfloat = irwb - 1;
	    i__2 = *nrhs;
	    for (jcol = 1; jcol <= i__2; ++jcol) {
		j += *n;
		i__3 = nsize;
		for (jrow = 1; jrow <= i__3; ++jrow) {
		    ++jfloat;
		    i__4 = j + jrow;
		    rwork[jfloat] = work[i__4].r;
/* L260: */
		}
/* L270: */
	    }
	    sgemm_("T", "N", &nsize, nrhs, &nsize, &c_b10, &rwork[vt + st1], 
		    n, &rwork[irwb], &nsize, &c_b35, &rwork[irwrb], &nsize);
	    j = bxst - *n - 1;
	    jimag = irwb - 1;
	    i__2 = *nrhs;
	    for (jcol = 1; jcol <= i__2; ++jcol) {
		j += *n;
		i__3 = nsize;
		for (jrow = 1; jrow <= i__3; ++jrow) {
		    ++jimag;
		    rwork[jimag] = r_imag(&work[j + jrow]);
/* L280: */
		}
/* L290: */
	    }
	    sgemm_("T", "N", &nsize, nrhs, &nsize, &c_b10, &rwork[vt + st1], 
		    n, &rwork[irwb], &nsize, &c_b35, &rwork[irwib], &nsize);
	    jfloat = irwrb - 1;
	    jimag = irwib - 1;
	    i__2 = *nrhs;
	    for (jcol = 1; jcol <= i__2; ++jcol) {
		i__3 = st + nsize - 1;
		for (jrow = st; jrow <= i__3; ++jrow) {
		    ++jfloat;
		    ++jimag;
		    i__4 = jrow + jcol * b_dim1;
		    i__5 = jfloat;
		    i__6 = jimag;
		    q__1.r = rwork[i__5], q__1.i = rwork[i__6];
		    b[i__4].r = q__1.r, b[i__4].i = q__1.i;
/* L300: */
		}
/* L310: */
	    }
	} else {
	    clalsa_(&icmpq2, smlsiz, &nsize, nrhs, &work[bxst], n, &b[st + 
		    b_dim1], ldb, &rwork[u + st1], n, &rwork[vt + st1], &
		    iwork[k + st1], &rwork[difl + st1], &rwork[difr + st1], &
		    rwork[z__ + st1], &rwork[poles + st1], &iwork[givptr + 
		    st1], &iwork[givcol + st1], n, &iwork[perm + st1], &rwork[
		    givnum + st1], &rwork[c__ + st1], &rwork[s + st1], &rwork[
		    nrwork], &iwork[iwk], info);
	    if (*info != 0) {
		return 0;
	    }
	}
/* L320: */
    }

/*     Unscale and sort the singular values. */

    slascl_("G", &c__0, &c__0, &c_b10, &orgnrm, n, &c__1, &d__[1], n, info);
    slasrt_("D", n, &d__[1], info);
    clascl_("G", &c__0, &c__0, &orgnrm, &c_b10, n, nrhs, &b[b_offset], ldb, 
	    info);

    return 0;

/*     End of CLALSD */

} /* clalsd_ */
