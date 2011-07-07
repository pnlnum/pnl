/* slanv2.f -- translated by f2c (version 20061008).
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

static float c_b4 = 1.f;

 int slanv2_(float *a, float *b, float *c__, float *d__, float *
	rt1r, float *rt1i, float *rt2r, float *rt2i, float *cs, float *sn)
{
    /* System generated locals */
    float r__1, r__2;

    /* Builtin functions */
    double r_sign(float *, float *), sqrt(double);

    /* Local variables */
    float p, z__, aa, bb, cc, dd, cs1, sn1, sab, sac, eps, tau, temp, scale, 
	    bcmax, bcmis, sigma;
    extern double slapy2_(float *, float *), slamch_(char *);


/*  -- LAPACK driver routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SLANV2 computes the Schur factorization of a float 2-by-2 nonsymmetric */
/*  matrix in standard form: */

/*       [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ] */
/*       [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ] */

/*  where either */
/*  1) CC = 0 so that AA and DD are float eigenvalues of the matrix, or */
/*  2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex */
/*  conjugate eigenvalues. */

/*  Arguments */
/*  ========= */

/*  A       (input/output) REAL */
/*  B       (input/output) REAL */
/*  C       (input/output) REAL */
/*  D       (input/output) REAL */
/*          On entry, the elements of the input matrix. */
/*          On exit, they are overwritten by the elements of the */
/*          standardised Schur form. */

/*  RT1R    (output) REAL */
/*  RT1I    (output) REAL */
/*  RT2R    (output) REAL */
/*  RT2I    (output) REAL */
/*          The float and imaginary parts of the eigenvalues. If the */
/*          eigenvalues are a complex conjugate pair, RT1I > 0. */

/*  CS      (output) REAL */
/*  SN      (output) REAL */
/*          Parameters of the rotation matrix. */

/*  Further Details */
/*  =============== */

/*  Modified by V. Sima, Research Institute for Informatics, Bucharest, */
/*  Romania, to reduce the risk of cancellation errors, */
/*  when computing float eigenvalues, and to ensure, if possible, that */
/*  ABS(RT1R) >= ABS(RT2R). */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    eps = slamch_("P");
    if (*c__ == 0.f) {
	*cs = 1.f;
	*sn = 0.f;
	goto L10;

    } else if (*b == 0.f) {

/*        Swap rows and columns */

	*cs = 0.f;
	*sn = 1.f;
	temp = *d__;
	*d__ = *a;
	*a = temp;
	*b = -(*c__);
	*c__ = 0.f;
	goto L10;
    } else if (*a - *d__ == 0.f && r_sign(&c_b4, b) != r_sign(&c_b4, c__)) {
	*cs = 1.f;
	*sn = 0.f;
	goto L10;
    } else {

	temp = *a - *d__;
	p = temp * .5f;
/* Computing MAX */
	r__1 = ABS(*b), r__2 = ABS(*c__);
	bcmax = MAX(r__1,r__2);
/* Computing MIN */
	r__1 = ABS(*b), r__2 = ABS(*c__);
	bcmis = MIN(r__1,r__2) * r_sign(&c_b4, b) * r_sign(&c_b4, c__);
/* Computing MAX */
	r__1 = ABS(p);
	scale = MAX(r__1,bcmax);
	z__ = p / scale * p + bcmax / scale * bcmis;

/*        If Z is of the order of the machine accuracy, postpone the */
/*        decision on the nature of eigenvalues */

	if (z__ >= eps * 4.f) {

/*           Real eigenvalues. Compute A and D. */

	    r__1 = sqrt(scale) * sqrt(z__);
	    z__ = p + r_sign(&r__1, &p);
	    *a = *d__ + z__;
	    *d__ -= bcmax / z__ * bcmis;

/*           Compute B and the rotation matrix */

	    tau = slapy2_(c__, &z__);
	    *cs = z__ / tau;
	    *sn = *c__ / tau;
	    *b -= *c__;
	    *c__ = 0.f;
	} else {

/*           Complex eigenvalues, or float (almost) equal eigenvalues. */
/*           Make diagonal elements equal. */

	    sigma = *b + *c__;
	    tau = slapy2_(&sigma, &temp);
	    *cs = sqrt((ABS(sigma) / tau + 1.f) * .5f);
	    *sn = -(p / (tau * *cs)) * r_sign(&c_b4, &sigma);

/*           Compute [ AA  BB ] = [ A  B ] [ CS -SN ] */
/*                   [ CC  DD ]   [ C  D ] [ SN  CS ] */

	    aa = *a * *cs + *b * *sn;
	    bb = -(*a) * *sn + *b * *cs;
	    cc = *c__ * *cs + *d__ * *sn;
	    dd = -(*c__) * *sn + *d__ * *cs;

/*           Compute [ A  B ] = [ CS  SN ] [ AA  BB ] */
/*                   [ C  D ]   [-SN  CS ] [ CC  DD ] */

	    *a = aa * *cs + cc * *sn;
	    *b = bb * *cs + dd * *sn;
	    *c__ = -aa * *sn + cc * *cs;
	    *d__ = -bb * *sn + dd * *cs;

	    temp = (*a + *d__) * .5f;
	    *a = temp;
	    *d__ = temp;

	    if (*c__ != 0.f) {
		if (*b != 0.f) {
		    if (r_sign(&c_b4, b) == r_sign(&c_b4, c__)) {

/*                    Real eigenvalues: reduce to upper triangular form */

			sab = sqrt((ABS(*b)));
			sac = sqrt((ABS(*c__)));
			r__1 = sab * sac;
			p = r_sign(&r__1, c__);
			tau = 1.f / sqrt((r__1 = *b + *c__, ABS(r__1)));
			*a = temp + p;
			*d__ = temp - p;
			*b -= *c__;
			*c__ = 0.f;
			cs1 = sab * tau;
			sn1 = sac * tau;
			temp = *cs * cs1 - *sn * sn1;
			*sn = *cs * sn1 + *sn * cs1;
			*cs = temp;
		    }
		} else {
		    *b = -(*c__);
		    *c__ = 0.f;
		    temp = *cs;
		    *cs = -(*sn);
		    *sn = temp;
		}
	    }
	}

    }

L10:

/*     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I). */

    *rt1r = *a;
    *rt2r = *d__;
    if (*c__ == 0.f) {
	*rt1i = 0.f;
	*rt2i = 0.f;
    } else {
	*rt1i = sqrt((ABS(*b))) * sqrt((ABS(*c__)));
	*rt2i = -(*rt1i);
    }
    return 0;

/*     End of SLANV2 */

} /* slanv2_ */
