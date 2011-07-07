/* clarfg.f -- translated by f2c (version 20061008).
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

static complex c_b5 = {1.f,0.f};

 int clarfg_(int *n, complex *alpha, complex *x, int *
	incx, complex *tau)
{
    /* System generated locals */
    int i__1;
    float r__1, r__2;
    complex q__1, q__2;

    /* Builtin functions */
    double r_imag(complex *), r_sign(float *, float *);

    /* Local variables */
    int j, knt;
    float beta;
    extern  int cscal_(int *, complex *, complex *, 
	    int *);
    float alphi, alphr, xnorm;
    extern double scnrm2_(int *, complex *, int *), slapy3_(float *
, float *, float *);
    extern /* Complex */ VOID cladiv_(complex *, complex *, complex *);
    extern double slamch_(char *);
    extern  int csscal_(int *, float *, complex *, int 
	    *);
    float safmin, rsafmn;


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  CLARFG generates a complex elementary reflector H of order n, such */
/*  that */

/*        H' * ( alpha ) = ( beta ),   H' * H = I. */
/*             (   x   )   (   0  ) */

/*  where alpha and beta are scalars, with beta float, and x is an */
/*  (n-1)-element complex vector. H is represented in the form */

/*        H = I - tau * ( 1 ) * ( 1 v' ) , */
/*                      ( v ) */

/*  where tau is a complex scalar and v is a complex (n-1)-element */
/*  vector. Note that H is not hermitian. */

/*  If the elements of x are all zero and alpha is float, then tau = 0 */
/*  and H is taken to be the unit matrix. */

/*  Otherwise  1 <= float(tau) <= 2  and  ABS(tau-1) <= 1 . */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The order of the elementary reflector. */

/*  ALPHA   (input/output) COMPLEX */
/*          On entry, the value alpha. */
/*          On exit, it is overwritten with the value beta. */

/*  X       (input/output) COMPLEX array, dimension */
/*                         (1+(N-2)*ABS(INCX)) */
/*          On entry, the vector x. */
/*          On exit, it is overwritten with the vector v. */

/*  INCX    (input) INTEGER */
/*          The increment between elements of X. INCX > 0. */

/*  TAU     (output) COMPLEX */
/*          The value tau. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*n <= 0) {
	tau->r = 0.f, tau->i = 0.f;
	return 0;
    }

    i__1 = *n - 1;
    xnorm = scnrm2_(&i__1, &x[1], incx);
    alphr = alpha->r;
    alphi = r_imag(alpha);

    if (xnorm == 0.f && alphi == 0.f) {

/*        H  =  I */

	tau->r = 0.f, tau->i = 0.f;
    } else {

/*        general case */

	r__1 = slapy3_(&alphr, &alphi, &xnorm);
	beta = -r_sign(&r__1, &alphr);
	safmin = slamch_("S") / slamch_("E");
	rsafmn = 1.f / safmin;

	knt = 0;
	if (ABS(beta) < safmin) {

/*           XNORM, BETA may be inaccurate; scale X and recompute them */

L10:
	    ++knt;
	    i__1 = *n - 1;
	    csscal_(&i__1, &rsafmn, &x[1], incx);
	    beta *= rsafmn;
	    alphi *= rsafmn;
	    alphr *= rsafmn;
	    if (ABS(beta) < safmin) {
		goto L10;
	    }

/*           New BETA is at most 1, at least SAFMIN */

	    i__1 = *n - 1;
	    xnorm = scnrm2_(&i__1, &x[1], incx);
	    q__1.r = alphr, q__1.i = alphi;
	    alpha->r = q__1.r, alpha->i = q__1.i;
	    r__1 = slapy3_(&alphr, &alphi, &xnorm);
	    beta = -r_sign(&r__1, &alphr);
	}
	r__1 = (beta - alphr) / beta;
	r__2 = -alphi / beta;
	q__1.r = r__1, q__1.i = r__2;
	tau->r = q__1.r, tau->i = q__1.i;
	q__2.r = alpha->r - beta, q__2.i = alpha->i;
	cladiv_(&q__1, &c_b5, &q__2);
	alpha->r = q__1.r, alpha->i = q__1.i;
	i__1 = *n - 1;
	cscal_(&i__1, alpha, &x[1], incx);

/*        If ALPHA is subnormal, it may lose relative accuracy */

	i__1 = knt;
	for (j = 1; j <= i__1; ++j) {
	    beta *= safmin;
/* L20: */
	}
	alpha->r = beta, alpha->i = 0.f;
    }

    return 0;

/*     End of CLARFG */

} /* clarfg_ */
