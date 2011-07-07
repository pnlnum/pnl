/* ctgevc.f -- translated by f2c (version 20061008).
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

 int ctgevc_(char *side, char *howmny, int *select, 
	int *n, complex *s, int *lds, complex *p, int *ldp, 
	complex *vl, int *ldvl, complex *vr, int *ldvr, int *mm, 
	int *m, complex *work, float *rwork, int *info)
{
    /* System generated locals */
    int p_dim1, p_offset, s_dim1, s_offset, vl_dim1, vl_offset, vr_dim1, 
	    vr_offset, i__1, i__2, i__3, i__4, i__5;
    float r__1, r__2, r__3, r__4, r__5, r__6;
    complex q__1, q__2, q__3, q__4;

    /* Builtin functions */
    double r_imag(complex *);
    void r_cnjg(complex *, complex *);

    /* Local variables */
    complex d__;
    int i__, j;
    complex ca, cb;
    int je, im, jr;
    float big;
    int lsa, lsb;
    float ulp;
    complex sum;
    int ibeg, ieig, iend;
    float dmin__;
    int isrc;
    float temp;
    complex suma, sumb;
    float xmax, scale;
    int ilall;
    int iside;
    float sbeta;
    extern int lsame_(char *, char *);
    extern  int cgemv_(char *, int *, int *, complex *
, complex *, int *, complex *, int *, complex *, complex *
, int *);
    float small;
    int compl;
    float anorm, bnorm;
    int compr, ilbbad;
    float acoefa, bcoefa, acoeff;
    complex bcoeff;
    int ilback;
    extern  int slabad_(float *, float *);
    float ascale, bscale;
    extern /* Complex */ VOID cladiv_(complex *, complex *, complex *);
    extern double slamch_(char *);
    complex salpha;
    float safmin;
    extern  int xerbla_(char *, int *);
    float bignum;
    int ilcomp;
    int ihwmny;


/*  -- LAPACK routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */


/*  Purpose */
/*  ======= */

/*  CTGEVC computes some or all of the right and/or left eigenvectors of */
/*  a pair of complex matrices (S,P), where S and P are upper triangular. */
/*  Matrix pairs of this type are produced by the generalized Schur */
/*  factorization of a complex matrix pair (A,B): */

/*     A = Q*S*Z**H,  B = Q*P*Z**H */

/*  as computed by CGGHRD + CHGEQZ. */

/*  The right eigenvector x and the left eigenvector y of (S,P) */
/*  corresponding to an eigenvalue w are defined by: */

/*     S*x = w*P*x,  (y**H)*S = w*(y**H)*P, */

/*  where y**H denotes the conjugate tranpose of y. */
/*  The eigenvalues are not input to this routine, but are computed */
/*  directly from the diagonal elements of S and P. */

/*  This routine returns the matrices X and/or Y of right and left */
/*  eigenvectors of (S,P), or the products Z*X and/or Q*Y, */
/*  where Z and Q are input matrices. */
/*  If Q and Z are the unitary factors from the generalized Schur */
/*  factorization of a matrix pair (A,B), then Z*X and Q*Y */
/*  are the matrices of right and left eigenvectors of (A,B). */

/*  Arguments */
/*  ========= */

/*  SIDE    (input) CHARACTER*1 */
/*          = 'R': compute right eigenvectors only; */
/*          = 'L': compute left eigenvectors only; */
/*          = 'B': compute both right and left eigenvectors. */

/*  HOWMNY  (input) CHARACTER*1 */
/*          = 'A': compute all right and/or left eigenvectors; */
/*          = 'B': compute all right and/or left eigenvectors, */
/*                 backtransformed by the matrices in VR and/or VL; */
/*          = 'S': compute selected right and/or left eigenvectors, */
/*                 specified by the int array SELECT. */

/*  SELECT  (input) LOGICAL array, dimension (N) */
/*          If HOWMNY='S', SELECT specifies the eigenvectors to be */
/*          computed.  The eigenvector corresponding to the j-th */
/*          eigenvalue is computed if SELECT(j) = .TRUE.. */
/*          Not referenced if HOWMNY = 'A' or 'B'. */

/*  N       (input) INTEGER */
/*          The order of the matrices S and P.  N >= 0. */

/*  S       (input) COMPLEX array, dimension (LDS,N) */
/*          The upper triangular matrix S from a generalized Schur */
/*          factorization, as computed by CHGEQZ. */

/*  LDS     (input) INTEGER */
/*          The leading dimension of array S.  LDS >= MAX(1,N). */

/*  P       (input) COMPLEX array, dimension (LDP,N) */
/*          The upper triangular matrix P from a generalized Schur */
/*          factorization, as computed by CHGEQZ.  P must have float */
/*          diagonal elements. */

/*  LDP     (input) INTEGER */
/*          The leading dimension of array P.  LDP >= MAX(1,N). */

/*  VL      (input/output) COMPLEX array, dimension (LDVL,MM) */
/*          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must */
/*          contain an N-by-N matrix Q (usually the unitary matrix Q */
/*          of left Schur vectors returned by CHGEQZ). */
/*          On exit, if SIDE = 'L' or 'B', VL contains: */
/*          if HOWMNY = 'A', the matrix Y of left eigenvectors of (S,P); */
/*          if HOWMNY = 'B', the matrix Q*Y; */
/*          if HOWMNY = 'S', the left eigenvectors of (S,P) specified by */
/*                      SELECT, stored consecutively in the columns of */
/*                      VL, in the same order as their eigenvalues. */
/*          Not referenced if SIDE = 'R'. */

/*  LDVL    (input) INTEGER */
/*          The leading dimension of array VL.  LDVL >= 1, and if */
/*          SIDE = 'L' or 'l' or 'B' or 'b', LDVL >= N. */

/*  VR      (input/output) COMPLEX array, dimension (LDVR,MM) */
/*          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must */
/*          contain an N-by-N matrix Q (usually the unitary matrix Z */
/*          of right Schur vectors returned by CHGEQZ). */
/*          On exit, if SIDE = 'R' or 'B', VR contains: */
/*          if HOWMNY = 'A', the matrix X of right eigenvectors of (S,P); */
/*          if HOWMNY = 'B', the matrix Z*X; */
/*          if HOWMNY = 'S', the right eigenvectors of (S,P) specified by */
/*                      SELECT, stored consecutively in the columns of */
/*                      VR, in the same order as their eigenvalues. */
/*          Not referenced if SIDE = 'L'. */

/*  LDVR    (input) INTEGER */
/*          The leading dimension of the array VR.  LDVR >= 1, and if */
/*          SIDE = 'R' or 'B', LDVR >= N. */

/*  MM      (input) INTEGER */
/*          The number of columns in the arrays VL and/or VR. MM >= M. */

/*  M       (output) INTEGER */
/*          The number of columns in the arrays VL and/or VR actually */
/*          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M */
/*          is set to N.  Each selected eigenvector occupies one column. */

/*  WORK    (workspace) COMPLEX array, dimension (2*N) */

/*  RWORK   (workspace) REAL array, dimension (2*N) */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit. */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */

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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test the input parameters */

    /* Parameter adjustments */
    --select;
    s_dim1 = *lds;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    p_dim1 = *ldp;
    p_offset = 1 + p_dim1;
    p -= p_offset;
    vl_dim1 = *ldvl;
    vl_offset = 1 + vl_dim1;
    vl -= vl_offset;
    vr_dim1 = *ldvr;
    vr_offset = 1 + vr_dim1;
    vr -= vr_offset;
    --work;
    --rwork;

    /* Function Body */
    if (lsame_(howmny, "A")) {
	ihwmny = 1;
	ilall = TRUE;
	ilback = FALSE;
    } else if (lsame_(howmny, "S")) {
	ihwmny = 2;
	ilall = FALSE;
	ilback = FALSE;
    } else if (lsame_(howmny, "B")) {
	ihwmny = 3;
	ilall = TRUE;
	ilback = TRUE;
    } else {
	ihwmny = -1;
    }

    if (lsame_(side, "R")) {
	iside = 1;
	compl = FALSE;
	compr = TRUE;
    } else if (lsame_(side, "L")) {
	iside = 2;
	compl = TRUE;
	compr = FALSE;
    } else if (lsame_(side, "B")) {
	iside = 3;
	compl = TRUE;
	compr = TRUE;
    } else {
	iside = -1;
    }

    *info = 0;
    if (iside < 0) {
	*info = -1;
    } else if (ihwmny < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -4;
    } else if (*lds < MAX(1,*n)) {
	*info = -6;
    } else if (*ldp < MAX(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CTGEVC", &i__1);
	return 0;
    }

/*     Count the number of eigenvectors */

    if (! ilall) {
	im = 0;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (select[j]) {
		++im;
	    }
/* L10: */
	}
    } else {
	im = *n;
    }

/*     Check diagonal of B */

    ilbbad = FALSE;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (r_imag(&p[j + j * p_dim1]) != 0.f) {
	    ilbbad = TRUE;
	}
/* L20: */
    }

    if (ilbbad) {
	*info = -7;
    } else if (compl && *ldvl < *n || *ldvl < 1) {
	*info = -10;
    } else if (compr && *ldvr < *n || *ldvr < 1) {
	*info = -12;
    } else if (*mm < im) {
	*info = -13;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CTGEVC", &i__1);
	return 0;
    }

/*     Quick return if possible */

    *m = im;
    if (*n == 0) {
	return 0;
    }

/*     Machine Constants */

    safmin = slamch_("Safe minimum");
    big = 1.f / safmin;
    slabad_(&safmin, &big);
    ulp = slamch_("Epsilon") * slamch_("Base");
    small = safmin * *n / ulp;
    big = 1.f / small;
    bignum = 1.f / (safmin * *n);

/*     Compute the 1-norm of each column of the strictly upper triangular */
/*     part of A and B to check for possible overflow in the triangular */
/*     solver. */

    i__1 = s_dim1 + 1;
    anorm = (r__1 = s[i__1].r, ABS(r__1)) + (r__2 = r_imag(&s[s_dim1 + 1]), 
	    ABS(r__2));
    i__1 = p_dim1 + 1;
    bnorm = (r__1 = p[i__1].r, ABS(r__1)) + (r__2 = r_imag(&p[p_dim1 + 1]), 
	    ABS(r__2));
    rwork[1] = 0.f;
    rwork[*n + 1] = 0.f;
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	rwork[j] = 0.f;
	rwork[*n + j] = 0.f;
	i__2 = j - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + j * s_dim1;
	    rwork[j] += (r__1 = s[i__3].r, ABS(r__1)) + (r__2 = r_imag(&s[
		    i__ + j * s_dim1]), ABS(r__2));
	    i__3 = i__ + j * p_dim1;
	    rwork[*n + j] += (r__1 = p[i__3].r, ABS(r__1)) + (r__2 = r_imag(&
		    p[i__ + j * p_dim1]), ABS(r__2));
/* L30: */
	}
/* Computing MAX */
	i__2 = j + j * s_dim1;
	r__3 = anorm, r__4 = rwork[j] + ((r__1 = s[i__2].r, ABS(r__1)) + (
		r__2 = r_imag(&s[j + j * s_dim1]), ABS(r__2)));
	anorm = MAX(r__3,r__4);
/* Computing MAX */
	i__2 = j + j * p_dim1;
	r__3 = bnorm, r__4 = rwork[*n + j] + ((r__1 = p[i__2].r, ABS(r__1)) 
		+ (r__2 = r_imag(&p[j + j * p_dim1]), ABS(r__2)));
	bnorm = MAX(r__3,r__4);
/* L40: */
    }

    ascale = 1.f / MAX(anorm,safmin);
    bscale = 1.f / MAX(bnorm,safmin);

/*     Left eigenvectors */

    if (compl) {
	ieig = 0;

/*        Main loop over eigenvalues */

	i__1 = *n;
	for (je = 1; je <= i__1; ++je) {
	    if (ilall) {
		ilcomp = TRUE;
	    } else {
		ilcomp = select[je];
	    }
	    if (ilcomp) {
		++ieig;

		i__2 = je + je * s_dim1;
		i__3 = je + je * p_dim1;
		if ((r__2 = s[i__2].r, ABS(r__2)) + (r__3 = r_imag(&s[je + 
			je * s_dim1]), ABS(r__3)) <= safmin && (r__1 = p[
			i__3].r, ABS(r__1)) <= safmin) {

/*                 Singular matrix pencil -- return unit eigenvector */

		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
			i__3 = jr + ieig * vl_dim1;
			vl[i__3].r = 0.f, vl[i__3].i = 0.f;
/* L50: */
		    }
		    i__2 = ieig + ieig * vl_dim1;
		    vl[i__2].r = 1.f, vl[i__2].i = 0.f;
		    goto L140;
		}

/*              Non-singular eigenvalue: */
/*              Compute coefficients  a  and  b  in */
/*                   H */
/*                 y  ( a A - b B ) = 0 */

/* Computing MAX */
		i__2 = je + je * s_dim1;
		i__3 = je + je * p_dim1;
		r__4 = ((r__2 = s[i__2].r, ABS(r__2)) + (r__3 = r_imag(&s[je 
			+ je * s_dim1]), ABS(r__3))) * ascale, r__5 = (r__1 =
			 p[i__3].r, ABS(r__1)) * bscale, r__4 = MAX(r__4,
			r__5);
		temp = 1.f / MAX(r__4,safmin);
		i__2 = je + je * s_dim1;
		q__2.r = temp * s[i__2].r, q__2.i = temp * s[i__2].i;
		q__1.r = ascale * q__2.r, q__1.i = ascale * q__2.i;
		salpha.r = q__1.r, salpha.i = q__1.i;
		i__2 = je + je * p_dim1;
		sbeta = temp * p[i__2].r * bscale;
		acoeff = sbeta * ascale;
		q__1.r = bscale * salpha.r, q__1.i = bscale * salpha.i;
		bcoeff.r = q__1.r, bcoeff.i = q__1.i;

/*              Scale to avoid underflow */

		lsa = ABS(sbeta) >= safmin && ABS(acoeff) < small;
		lsb = (r__1 = salpha.r, ABS(r__1)) + (r__2 = r_imag(&salpha),
			 ABS(r__2)) >= safmin && (r__3 = bcoeff.r, ABS(r__3)
			) + (r__4 = r_imag(&bcoeff), ABS(r__4)) < small;

		scale = 1.f;
		if (lsa) {
		    scale = small / ABS(sbeta) * MIN(anorm,big);
		}
		if (lsb) {
/* Computing MAX */
		    r__3 = scale, r__4 = small / ((r__1 = salpha.r, ABS(r__1)
			    ) + (r__2 = r_imag(&salpha), ABS(r__2))) * MIN(
			    bnorm,big);
		    scale = MAX(r__3,r__4);
		}
		if (lsa || lsb) {
/* Computing MIN */
/* Computing MAX */
		    r__5 = 1.f, r__6 = ABS(acoeff), r__5 = MAX(r__5,r__6), 
			    r__6 = (r__1 = bcoeff.r, ABS(r__1)) + (r__2 = 
			    r_imag(&bcoeff), ABS(r__2));
		    r__3 = scale, r__4 = 1.f / (safmin * MAX(r__5,r__6));
		    scale = MIN(r__3,r__4);
		    if (lsa) {
			acoeff = ascale * (scale * sbeta);
		    } else {
			acoeff = scale * acoeff;
		    }
		    if (lsb) {
			q__2.r = scale * salpha.r, q__2.i = scale * salpha.i;
			q__1.r = bscale * q__2.r, q__1.i = bscale * q__2.i;
			bcoeff.r = q__1.r, bcoeff.i = q__1.i;
		    } else {
			q__1.r = scale * bcoeff.r, q__1.i = scale * bcoeff.i;
			bcoeff.r = q__1.r, bcoeff.i = q__1.i;
		    }
		}

		acoefa = ABS(acoeff);
		bcoefa = (r__1 = bcoeff.r, ABS(r__1)) + (r__2 = r_imag(&
			bcoeff), ABS(r__2));
		xmax = 1.f;
		i__2 = *n;
		for (jr = 1; jr <= i__2; ++jr) {
		    i__3 = jr;
		    work[i__3].r = 0.f, work[i__3].i = 0.f;
/* L60: */
		}
		i__2 = je;
		work[i__2].r = 1.f, work[i__2].i = 0.f;
/* Computing MAX */
		r__1 = ulp * acoefa * anorm, r__2 = ulp * bcoefa * bnorm, 
			r__1 = MAX(r__1,r__2);
		dmin__ = MAX(r__1,safmin);

/*                                              H */
/*              Triangular solve of  (a A - b B)  y = 0 */

/*                                      H */
/*              (rowwise in  (a A - b B) , or columnwise in a A - b B) */

		i__2 = *n;
		for (j = je + 1; j <= i__2; ++j) {

/*                 Compute */
/*                       j-1 */
/*                 SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k) */
/*                       k=je */
/*                 (Scale if necessary) */

		    temp = 1.f / xmax;
		    if (acoefa * rwork[j] + bcoefa * rwork[*n + j] > bignum * 
			    temp) {
			i__3 = j - 1;
			for (jr = je; jr <= i__3; ++jr) {
			    i__4 = jr;
			    i__5 = jr;
			    q__1.r = temp * work[i__5].r, q__1.i = temp * 
				    work[i__5].i;
			    work[i__4].r = q__1.r, work[i__4].i = q__1.i;
/* L70: */
			}
			xmax = 1.f;
		    }
		    suma.r = 0.f, suma.i = 0.f;
		    sumb.r = 0.f, sumb.i = 0.f;

		    i__3 = j - 1;
		    for (jr = je; jr <= i__3; ++jr) {
			r_cnjg(&q__3, &s[jr + j * s_dim1]);
			i__4 = jr;
			q__2.r = q__3.r * work[i__4].r - q__3.i * work[i__4]
				.i, q__2.i = q__3.r * work[i__4].i + q__3.i * 
				work[i__4].r;
			q__1.r = suma.r + q__2.r, q__1.i = suma.i + q__2.i;
			suma.r = q__1.r, suma.i = q__1.i;
			r_cnjg(&q__3, &p[jr + j * p_dim1]);
			i__4 = jr;
			q__2.r = q__3.r * work[i__4].r - q__3.i * work[i__4]
				.i, q__2.i = q__3.r * work[i__4].i + q__3.i * 
				work[i__4].r;
			q__1.r = sumb.r + q__2.r, q__1.i = sumb.i + q__2.i;
			sumb.r = q__1.r, sumb.i = q__1.i;
/* L80: */
		    }
		    q__2.r = acoeff * suma.r, q__2.i = acoeff * suma.i;
		    r_cnjg(&q__4, &bcoeff);
		    q__3.r = q__4.r * sumb.r - q__4.i * sumb.i, q__3.i = 
			    q__4.r * sumb.i + q__4.i * sumb.r;
		    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
		    sum.r = q__1.r, sum.i = q__1.i;

/*                 Form x(j) = - SUM / conjg( a*S(j,j) - b*P(j,j) ) */

/*                 with scaling and perturbation of the denominator */

		    i__3 = j + j * s_dim1;
		    q__3.r = acoeff * s[i__3].r, q__3.i = acoeff * s[i__3].i;
		    i__4 = j + j * p_dim1;
		    q__4.r = bcoeff.r * p[i__4].r - bcoeff.i * p[i__4].i, 
			    q__4.i = bcoeff.r * p[i__4].i + bcoeff.i * p[i__4]
			    .r;
		    q__2.r = q__3.r - q__4.r, q__2.i = q__3.i - q__4.i;
		    r_cnjg(&q__1, &q__2);
		    d__.r = q__1.r, d__.i = q__1.i;
		    if ((r__1 = d__.r, ABS(r__1)) + (r__2 = r_imag(&d__), 
			    ABS(r__2)) <= dmin__) {
			q__1.r = dmin__, q__1.i = 0.f;
			d__.r = q__1.r, d__.i = q__1.i;
		    }

		    if ((r__1 = d__.r, ABS(r__1)) + (r__2 = r_imag(&d__), 
			    ABS(r__2)) < 1.f) {
			if ((r__1 = sum.r, ABS(r__1)) + (r__2 = r_imag(&sum),
				 ABS(r__2)) >= bignum * ((r__3 = d__.r, ABS(
				r__3)) + (r__4 = r_imag(&d__), ABS(r__4)))) {
			    temp = 1.f / ((r__1 = sum.r, ABS(r__1)) + (r__2 =
				     r_imag(&sum), ABS(r__2)));
			    i__3 = j - 1;
			    for (jr = je; jr <= i__3; ++jr) {
				i__4 = jr;
				i__5 = jr;
				q__1.r = temp * work[i__5].r, q__1.i = temp * 
					work[i__5].i;
				work[i__4].r = q__1.r, work[i__4].i = q__1.i;
/* L90: */
			    }
			    xmax = temp * xmax;
			    q__1.r = temp * sum.r, q__1.i = temp * sum.i;
			    sum.r = q__1.r, sum.i = q__1.i;
			}
		    }
		    i__3 = j;
		    q__2.r = -sum.r, q__2.i = -sum.i;
		    cladiv_(&q__1, &q__2, &d__);
		    work[i__3].r = q__1.r, work[i__3].i = q__1.i;
/* Computing MAX */
		    i__3 = j;
		    r__3 = xmax, r__4 = (r__1 = work[i__3].r, ABS(r__1)) + (
			    r__2 = r_imag(&work[j]), ABS(r__2));
		    xmax = MAX(r__3,r__4);
/* L100: */
		}

/*              Back transform eigenvector if HOWMNY='B'. */

		if (ilback) {
		    i__2 = *n + 1 - je;
		    cgemv_("N", n, &i__2, &c_b2, &vl[je * vl_dim1 + 1], ldvl, 
			    &work[je], &c__1, &c_b1, &work[*n + 1], &c__1);
		    isrc = 2;
		    ibeg = 1;
		} else {
		    isrc = 1;
		    ibeg = je;
		}

/*              Copy and scale eigenvector into column of VL */

		xmax = 0.f;
		i__2 = *n;
		for (jr = ibeg; jr <= i__2; ++jr) {
/* Computing MAX */
		    i__3 = (isrc - 1) * *n + jr;
		    r__3 = xmax, r__4 = (r__1 = work[i__3].r, ABS(r__1)) + (
			    r__2 = r_imag(&work[(isrc - 1) * *n + jr]), ABS(
			    r__2));
		    xmax = MAX(r__3,r__4);
/* L110: */
		}

		if (xmax > safmin) {
		    temp = 1.f / xmax;
		    i__2 = *n;
		    for (jr = ibeg; jr <= i__2; ++jr) {
			i__3 = jr + ieig * vl_dim1;
			i__4 = (isrc - 1) * *n + jr;
			q__1.r = temp * work[i__4].r, q__1.i = temp * work[
				i__4].i;
			vl[i__3].r = q__1.r, vl[i__3].i = q__1.i;
/* L120: */
		    }
		} else {
		    ibeg = *n + 1;
		}

		i__2 = ibeg - 1;
		for (jr = 1; jr <= i__2; ++jr) {
		    i__3 = jr + ieig * vl_dim1;
		    vl[i__3].r = 0.f, vl[i__3].i = 0.f;
/* L130: */
		}

	    }
L140:
	    ;
	}
    }

/*     Right eigenvectors */

    if (compr) {
	ieig = im + 1;

/*        Main loop over eigenvalues */

	for (je = *n; je >= 1; --je) {
	    if (ilall) {
		ilcomp = TRUE;
	    } else {
		ilcomp = select[je];
	    }
	    if (ilcomp) {
		--ieig;

		i__1 = je + je * s_dim1;
		i__2 = je + je * p_dim1;
		if ((r__2 = s[i__1].r, ABS(r__2)) + (r__3 = r_imag(&s[je + 
			je * s_dim1]), ABS(r__3)) <= safmin && (r__1 = p[
			i__2].r, ABS(r__1)) <= safmin) {

/*                 Singular matrix pencil -- return unit eigenvector */

		    i__1 = *n;
		    for (jr = 1; jr <= i__1; ++jr) {
			i__2 = jr + ieig * vr_dim1;
			vr[i__2].r = 0.f, vr[i__2].i = 0.f;
/* L150: */
		    }
		    i__1 = ieig + ieig * vr_dim1;
		    vr[i__1].r = 1.f, vr[i__1].i = 0.f;
		    goto L250;
		}

/*              Non-singular eigenvalue: */
/*              Compute coefficients  a  and  b  in */

/*              ( a A - b B ) x  = 0 */

/* Computing MAX */
		i__1 = je + je * s_dim1;
		i__2 = je + je * p_dim1;
		r__4 = ((r__2 = s[i__1].r, ABS(r__2)) + (r__3 = r_imag(&s[je 
			+ je * s_dim1]), ABS(r__3))) * ascale, r__5 = (r__1 =
			 p[i__2].r, ABS(r__1)) * bscale, r__4 = MAX(r__4,
			r__5);
		temp = 1.f / MAX(r__4,safmin);
		i__1 = je + je * s_dim1;
		q__2.r = temp * s[i__1].r, q__2.i = temp * s[i__1].i;
		q__1.r = ascale * q__2.r, q__1.i = ascale * q__2.i;
		salpha.r = q__1.r, salpha.i = q__1.i;
		i__1 = je + je * p_dim1;
		sbeta = temp * p[i__1].r * bscale;
		acoeff = sbeta * ascale;
		q__1.r = bscale * salpha.r, q__1.i = bscale * salpha.i;
		bcoeff.r = q__1.r, bcoeff.i = q__1.i;

/*              Scale to avoid underflow */

		lsa = ABS(sbeta) >= safmin && ABS(acoeff) < small;
		lsb = (r__1 = salpha.r, ABS(r__1)) + (r__2 = r_imag(&salpha),
			 ABS(r__2)) >= safmin && (r__3 = bcoeff.r, ABS(r__3)
			) + (r__4 = r_imag(&bcoeff), ABS(r__4)) < small;

		scale = 1.f;
		if (lsa) {
		    scale = small / ABS(sbeta) * MIN(anorm,big);
		}
		if (lsb) {
/* Computing MAX */
		    r__3 = scale, r__4 = small / ((r__1 = salpha.r, ABS(r__1)
			    ) + (r__2 = r_imag(&salpha), ABS(r__2))) * MIN(
			    bnorm,big);
		    scale = MAX(r__3,r__4);
		}
		if (lsa || lsb) {
/* Computing MIN */
/* Computing MAX */
		    r__5 = 1.f, r__6 = ABS(acoeff), r__5 = MAX(r__5,r__6), 
			    r__6 = (r__1 = bcoeff.r, ABS(r__1)) + (r__2 = 
			    r_imag(&bcoeff), ABS(r__2));
		    r__3 = scale, r__4 = 1.f / (safmin * MAX(r__5,r__6));
		    scale = MIN(r__3,r__4);
		    if (lsa) {
			acoeff = ascale * (scale * sbeta);
		    } else {
			acoeff = scale * acoeff;
		    }
		    if (lsb) {
			q__2.r = scale * salpha.r, q__2.i = scale * salpha.i;
			q__1.r = bscale * q__2.r, q__1.i = bscale * q__2.i;
			bcoeff.r = q__1.r, bcoeff.i = q__1.i;
		    } else {
			q__1.r = scale * bcoeff.r, q__1.i = scale * bcoeff.i;
			bcoeff.r = q__1.r, bcoeff.i = q__1.i;
		    }
		}

		acoefa = ABS(acoeff);
		bcoefa = (r__1 = bcoeff.r, ABS(r__1)) + (r__2 = r_imag(&
			bcoeff), ABS(r__2));
		xmax = 1.f;
		i__1 = *n;
		for (jr = 1; jr <= i__1; ++jr) {
		    i__2 = jr;
		    work[i__2].r = 0.f, work[i__2].i = 0.f;
/* L160: */
		}
		i__1 = je;
		work[i__1].r = 1.f, work[i__1].i = 0.f;
/* Computing MAX */
		r__1 = ulp * acoefa * anorm, r__2 = ulp * bcoefa * bnorm, 
			r__1 = MAX(r__1,r__2);
		dmin__ = MAX(r__1,safmin);

/*              Triangular solve of  (a A - b B) x = 0  (columnwise) */

/*              WORK(1:j-1) contains sums w, */
/*              WORK(j+1:JE) contains x */

		i__1 = je - 1;
		for (jr = 1; jr <= i__1; ++jr) {
		    i__2 = jr;
		    i__3 = jr + je * s_dim1;
		    q__2.r = acoeff * s[i__3].r, q__2.i = acoeff * s[i__3].i;
		    i__4 = jr + je * p_dim1;
		    q__3.r = bcoeff.r * p[i__4].r - bcoeff.i * p[i__4].i, 
			    q__3.i = bcoeff.r * p[i__4].i + bcoeff.i * p[i__4]
			    .r;
		    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
		    work[i__2].r = q__1.r, work[i__2].i = q__1.i;
/* L170: */
		}
		i__1 = je;
		work[i__1].r = 1.f, work[i__1].i = 0.f;

		for (j = je - 1; j >= 1; --j) {

/*                 Form x(j) := - w(j) / d */
/*                 with scaling and perturbation of the denominator */

		    i__1 = j + j * s_dim1;
		    q__2.r = acoeff * s[i__1].r, q__2.i = acoeff * s[i__1].i;
		    i__2 = j + j * p_dim1;
		    q__3.r = bcoeff.r * p[i__2].r - bcoeff.i * p[i__2].i, 
			    q__3.i = bcoeff.r * p[i__2].i + bcoeff.i * p[i__2]
			    .r;
		    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
		    d__.r = q__1.r, d__.i = q__1.i;
		    if ((r__1 = d__.r, ABS(r__1)) + (r__2 = r_imag(&d__), 
			    ABS(r__2)) <= dmin__) {
			q__1.r = dmin__, q__1.i = 0.f;
			d__.r = q__1.r, d__.i = q__1.i;
		    }

		    if ((r__1 = d__.r, ABS(r__1)) + (r__2 = r_imag(&d__), 
			    ABS(r__2)) < 1.f) {
			i__1 = j;
			if ((r__1 = work[i__1].r, ABS(r__1)) + (r__2 = 
				r_imag(&work[j]), ABS(r__2)) >= bignum * ((
				r__3 = d__.r, ABS(r__3)) + (r__4 = r_imag(&
				d__), ABS(r__4)))) {
			    i__1 = j;
			    temp = 1.f / ((r__1 = work[i__1].r, ABS(r__1)) + 
				    (r__2 = r_imag(&work[j]), ABS(r__2)));
			    i__1 = je;
			    for (jr = 1; jr <= i__1; ++jr) {
				i__2 = jr;
				i__3 = jr;
				q__1.r = temp * work[i__3].r, q__1.i = temp * 
					work[i__3].i;
				work[i__2].r = q__1.r, work[i__2].i = q__1.i;
/* L180: */
			    }
			}
		    }

		    i__1 = j;
		    i__2 = j;
		    q__2.r = -work[i__2].r, q__2.i = -work[i__2].i;
		    cladiv_(&q__1, &q__2, &d__);
		    work[i__1].r = q__1.r, work[i__1].i = q__1.i;

		    if (j > 1) {

/*                    w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling */

			i__1 = j;
			if ((r__1 = work[i__1].r, ABS(r__1)) + (r__2 = 
				r_imag(&work[j]), ABS(r__2)) > 1.f) {
			    i__1 = j;
			    temp = 1.f / ((r__1 = work[i__1].r, ABS(r__1)) + 
				    (r__2 = r_imag(&work[j]), ABS(r__2)));
			    if (acoefa * rwork[j] + bcoefa * rwork[*n + j] >= 
				    bignum * temp) {
				i__1 = je;
				for (jr = 1; jr <= i__1; ++jr) {
				    i__2 = jr;
				    i__3 = jr;
				    q__1.r = temp * work[i__3].r, q__1.i = 
					    temp * work[i__3].i;
				    work[i__2].r = q__1.r, work[i__2].i = 
					    q__1.i;
/* L190: */
				}
			    }
			}

			i__1 = j;
			q__1.r = acoeff * work[i__1].r, q__1.i = acoeff * 
				work[i__1].i;
			ca.r = q__1.r, ca.i = q__1.i;
			i__1 = j;
			q__1.r = bcoeff.r * work[i__1].r - bcoeff.i * work[
				i__1].i, q__1.i = bcoeff.r * work[i__1].i + 
				bcoeff.i * work[i__1].r;
			cb.r = q__1.r, cb.i = q__1.i;
			i__1 = j - 1;
			for (jr = 1; jr <= i__1; ++jr) {
			    i__2 = jr;
			    i__3 = jr;
			    i__4 = jr + j * s_dim1;
			    q__3.r = ca.r * s[i__4].r - ca.i * s[i__4].i, 
				    q__3.i = ca.r * s[i__4].i + ca.i * s[i__4]
				    .r;
			    q__2.r = work[i__3].r + q__3.r, q__2.i = work[
				    i__3].i + q__3.i;
			    i__5 = jr + j * p_dim1;
			    q__4.r = cb.r * p[i__5].r - cb.i * p[i__5].i, 
				    q__4.i = cb.r * p[i__5].i + cb.i * p[i__5]
				    .r;
			    q__1.r = q__2.r - q__4.r, q__1.i = q__2.i - 
				    q__4.i;
			    work[i__2].r = q__1.r, work[i__2].i = q__1.i;
/* L200: */
			}
		    }
/* L210: */
		}

/*              Back transform eigenvector if HOWMNY='B'. */

		if (ilback) {
		    cgemv_("N", n, &je, &c_b2, &vr[vr_offset], ldvr, &work[1], 
			     &c__1, &c_b1, &work[*n + 1], &c__1);
		    isrc = 2;
		    iend = *n;
		} else {
		    isrc = 1;
		    iend = je;
		}

/*              Copy and scale eigenvector into column of VR */

		xmax = 0.f;
		i__1 = iend;
		for (jr = 1; jr <= i__1; ++jr) {
/* Computing MAX */
		    i__2 = (isrc - 1) * *n + jr;
		    r__3 = xmax, r__4 = (r__1 = work[i__2].r, ABS(r__1)) + (
			    r__2 = r_imag(&work[(isrc - 1) * *n + jr]), ABS(
			    r__2));
		    xmax = MAX(r__3,r__4);
/* L220: */
		}

		if (xmax > safmin) {
		    temp = 1.f / xmax;
		    i__1 = iend;
		    for (jr = 1; jr <= i__1; ++jr) {
			i__2 = jr + ieig * vr_dim1;
			i__3 = (isrc - 1) * *n + jr;
			q__1.r = temp * work[i__3].r, q__1.i = temp * work[
				i__3].i;
			vr[i__2].r = q__1.r, vr[i__2].i = q__1.i;
/* L230: */
		    }
		} else {
		    iend = 0;
		}

		i__1 = *n;
		for (jr = iend + 1; jr <= i__1; ++jr) {
		    i__2 = jr + ieig * vr_dim1;
		    vr[i__2].r = 0.f, vr[i__2].i = 0.f;
/* L240: */
		}

	    }
L250:
	    ;
	}
    }

    return 0;

/*     End of CTGEVC */

} /* ctgevc_ */
