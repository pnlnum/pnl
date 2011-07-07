/* clarfp.f -- translated by f2c (version 20061008).
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

 int clarfp_(int *n, complex *alpha, complex *x, int *
	incx, complex *tau)
{
    /* System generated locals */
    int i__1, i__2;
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
    extern double scnrm2_(int *, complex *, int *), slapy2_(float *
, float *), slapy3_(float *, float *, float *);
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

/*  CLARFP generates a complex elementary reflector H of order n, such */
/*  that */

/*        H' * ( alpha ) = ( beta ),   H' * H = I. */
/*             (   x   )   (   0  ) */

/*  where alpha and beta are scalars, beta is float and non-negative, and */
/*  x is an (n-1)-element complex vector.  H is represented in the form */

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

/*        H  =  [1-alpha/ABS(alpha) 0; 0 I], sign chosen so ALPHA >= 0. */

	if (alphi == 0.f) {
	    if (alphr >= 0.f) {
/*              When TAU.eq.ZERO, the vector is special-cased to be */
/*              all zeros in the application routines.  We do not need */
/*              to clear it. */
		tau->r = 0.f, tau->i = 0.f;
	    } else {
/*              However, the application routines rely on explicit */
/*              zero checks when TAU.ne.ZERO, and we must clear X. */
		tau->r = 2.f, tau->i = 0.f;
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = (j - 1) * *incx + 1;
		    x[i__2].r = 0.f, x[i__2].i = 0.f;
		}
		q__1.r = -alpha->r, q__1.i = -alpha->i;
		alpha->r = q__1.r, alpha->i = q__1.i;
	    }
	} else {
/*           Only "reflecting" the diagonal entry to be float and non-negative. */
	    xnorm = slapy2_(&alphr, &alphi);
	    r__1 = 1.f - alphr / xnorm;
	    r__2 = -alphi / xnorm;
	    q__1.r = r__1, q__1.i = r__2;
	    tau->r = q__1.r, tau->i = q__1.i;
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = (j - 1) * *incx + 1;
		x[i__2].r = 0.f, x[i__2].i = 0.f;
	    }
	    alpha->r = xnorm, alpha->i = 0.f;
	}
    } else {

/*        general case */

	r__1 = slapy3_(&alphr, &alphi, &xnorm);
	beta = r_sign(&r__1, &alphr);
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
	    beta = r_sign(&r__1, &alphr);
	}
	q__1.r = alpha->r + beta, q__1.i = alpha->i;
	alpha->r = q__1.r, alpha->i = q__1.i;
	if (beta < 0.f) {
	    beta = -beta;
	    q__2.r = -alpha->r, q__2.i = -alpha->i;
	    q__1.r = q__2.r / beta, q__1.i = q__2.i / beta;
	    tau->r = q__1.r, tau->i = q__1.i;
	} else {
	    alphr = alphi * (alphi / alpha->r);
	    alphr += xnorm * (xnorm / alpha->r);
	    r__1 = alphr / beta;
	    r__2 = -alphi / beta;
	    q__1.r = r__1, q__1.i = r__2;
	    tau->r = q__1.r, tau->i = q__1.i;
	    r__1 = -alphr;
	    q__1.r = r__1, q__1.i = alphi;
	    alpha->r = q__1.r, alpha->i = q__1.i;
	}
	cladiv_(&q__1, &c_b5, alpha);
	alpha->r = q__1.r, alpha->i = q__1.i;
	i__1 = *n - 1;
	cscal_(&i__1, alpha, &x[1], incx);

/*        If BETA is subnormal, it may lose relative accuracy */

	i__1 = knt;
	for (j = 1; j <= i__1; ++j) {
	    beta *= safmin;
/* L20: */
	}
	alpha->r = beta, alpha->i = 0.f;
    }

    return 0;

/*     End of CLARFP */

} /* clarfp_ */
