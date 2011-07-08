/* zlags2.f -- translated by f2c (version 20061008).
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

 int zlags2_(int *upper, double *a1, doublecomplex *
	a2, double *a3, double *b1, doublecomplex *b2, double *b3, 
	 double *csu, doublecomplex *snu, double *csv, doublecomplex *
	snv, double *csq, doublecomplex *snq)
{
    /* System generated locals */
    double d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    double z_abs(doublecomplex *), d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    double a;
    doublecomplex b, c__;
    double d__;
    doublecomplex r__, d1;
    double s1, s2, fb, fc;
    doublecomplex ua11, ua12, ua21, ua22, vb11, vb12, vb21, vb22;
    double csl, csr, snl, snr, aua11, aua12, aua21, aua22, avb12, avb11, 
	    avb21, avb22, ua11r, ua22r, vb11r, vb22r;
    extern  int dlasv2_(double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *), zlartg_(doublecomplex *
, doublecomplex *, double *, doublecomplex *, doublecomplex *)
	    ;


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  ZLAGS2 computes 2-by-2 unitary matrices U, V and Q, such */
/*  that if ( UPPER ) then */

/*            U'*A*Q = U'*( A1 A2 )*Q = ( x  0  ) */
/*                        ( 0  A3 )     ( x  x  ) */
/*  and */
/*            V'*B*Q = V'*( B1 B2 )*Q = ( x  0  ) */
/*                        ( 0  B3 )     ( x  x  ) */

/*  or if ( .NOT.UPPER ) then */

/*            U'*A*Q = U'*( A1 0  )*Q = ( x  x  ) */
/*                        ( A2 A3 )     ( 0  x  ) */
/*  and */
/*            V'*B*Q = V'*( B1 0  )*Q = ( x  x  ) */
/*                        ( B2 B3 )     ( 0  x  ) */
/*  where */

/*    U = (     CSU      SNU ), V = (     CSV     SNV ), */
/*        ( -CONJG(SNU)  CSU )      ( -CONJG(SNV) CSV ) */

/*    Q = (     CSQ      SNQ ) */
/*        ( -CONJG(SNQ)  CSQ ) */

/*  Z' denotes the conjugate transpose of Z. */

/*  The rows of the transformed A and B are parallel. Moreover, if the */
/*  input 2-by-2 matrix A is not zero, then the transformed (1,1) entry */
/*  of A is not zero. If the input matrices A and B are both not zero, */
/*  then the transformed (2,2) element of B is not zero, except when the */
/*  first rows of input A and B are parallel and the second rows are */
/*  zero. */

/*  Arguments */
/*  ========= */

/*  UPPER   (input) LOGICAL */
/*          = .TRUE.: the input matrices A and B are upper triangular. */
/*          = .FALSE.: the input matrices A and B are lower triangular. */

/*  A1      (input) DOUBLE PRECISION */
/*  A2      (input) COMPLEX*16 */
/*  A3      (input) DOUBLE PRECISION */
/*          On entry, A1, A2 and A3 are elements of the input 2-by-2 */
/*          upper (lower) triangular matrix A. */

/*  B1      (input) DOUBLE PRECISION */
/*  B2      (input) COMPLEX*16 */
/*  B3      (input) DOUBLE PRECISION */
/*          On entry, B1, B2 and B3 are elements of the input 2-by-2 */
/*          upper (lower) triangular matrix B. */

/*  CSU     (output) DOUBLE PRECISION */
/*  SNU     (output) COMPLEX*16 */
/*          The desired unitary matrix U. */

/*  CSV     (output) DOUBLE PRECISION */
/*  SNV     (output) COMPLEX*16 */
/*          The desired unitary matrix V. */

/*  CSQ     (output) DOUBLE PRECISION */
/*  SNQ     (output) COMPLEX*16 */
/*          The desired unitary matrix Q. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
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

    if (*upper) {

/*        Input matrices A and B are upper triangular matrices */

/*        Form matrix C = A*adj(B) = ( a b ) */
/*                                   ( 0 d ) */

	a = *a1 * *b3;
	d__ = *a3 * *b1;
	z__2.r = *b1 * a2->r, z__2.i = *b1 * a2->i;
	z__3.r = *a1 * b2->r, z__3.i = *a1 * b2->i;
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	b.r = z__1.r, b.i = z__1.i;
	fb = z_abs(&b);

/*        Transform complex 2-by-2 matrix C to float matrix by unitary */
/*        diagonal matrix diag(1,D1). */

	d1.r = 1., d1.i = 0.;
	if (fb != 0.) {
	    z__1.r = b.r / fb, z__1.i = b.i / fb;
	    d1.r = z__1.r, d1.i = z__1.i;
	}

/*        The SVD of float 2 by 2 triangular C */

/*         ( CSL -SNL )*( A B )*(  CSR  SNR ) = ( R 0 ) */
/*         ( SNL  CSL ) ( 0 D ) ( -SNR  CSR )   ( 0 T ) */

	dlasv2_(&a, &fb, &d__, &s1, &s2, &snr, &csr, &snl, &csl);

	if (ABS(csl) >= ABS(snl) || ABS(csr) >= ABS(snr)) {

/*           Compute the (1,1) and (1,2) elements of U'*A and V'*B, */
/*           and (1,2) element of |U|'*|A| and |V|'*|B|. */

	    ua11r = csl * *a1;
	    z__2.r = csl * a2->r, z__2.i = csl * a2->i;
	    z__4.r = snl * d1.r, z__4.i = snl * d1.i;
	    z__3.r = *a3 * z__4.r, z__3.i = *a3 * z__4.i;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    ua12.r = z__1.r, ua12.i = z__1.i;

	    vb11r = csr * *b1;
	    z__2.r = csr * b2->r, z__2.i = csr * b2->i;
	    z__4.r = snr * d1.r, z__4.i = snr * d1.i;
	    z__3.r = *b3 * z__4.r, z__3.i = *b3 * z__4.i;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    vb12.r = z__1.r, vb12.i = z__1.i;

	    aua12 = ABS(csl) * ((d__1 = a2->r, ABS(d__1)) + (d__2 = d_imag(a2)
		    , ABS(d__2))) + ABS(snl) * ABS(*a3);
	    avb12 = ABS(csr) * ((d__1 = b2->r, ABS(d__1)) + (d__2 = d_imag(b2)
		    , ABS(d__2))) + ABS(snr) * ABS(*b3);

/*           zero (1,2) elements of U'*A and V'*B */

	    if (ABS(ua11r) + ((d__1 = ua12.r, ABS(d__1)) + (d__2 = d_imag(&
		    ua12), ABS(d__2))) == 0.) {
		z__2.r = vb11r, z__2.i = 0.;
		z__1.r = -z__2.r, z__1.i = -z__2.i;
		d_cnjg(&z__3, &vb12);
		zlartg_(&z__1, &z__3, csq, snq, &r__);
	    } else if (ABS(vb11r) + ((d__1 = vb12.r, ABS(d__1)) + (d__2 = 
		    d_imag(&vb12), ABS(d__2))) == 0.) {
		z__2.r = ua11r, z__2.i = 0.;
		z__1.r = -z__2.r, z__1.i = -z__2.i;
		d_cnjg(&z__3, &ua12);
		zlartg_(&z__1, &z__3, csq, snq, &r__);
	    } else if (aua12 / (ABS(ua11r) + ((d__1 = ua12.r, ABS(d__1)) + (
		    d__2 = d_imag(&ua12), ABS(d__2)))) <= avb12 / (ABS(vb11r) 
		    + ((d__3 = vb12.r, ABS(d__3)) + (d__4 = d_imag(&vb12), 
		    ABS(d__4))))) {
		z__2.r = ua11r, z__2.i = 0.;
		z__1.r = -z__2.r, z__1.i = -z__2.i;
		d_cnjg(&z__3, &ua12);
		zlartg_(&z__1, &z__3, csq, snq, &r__);
	    } else {
		z__2.r = vb11r, z__2.i = 0.;
		z__1.r = -z__2.r, z__1.i = -z__2.i;
		d_cnjg(&z__3, &vb12);
		zlartg_(&z__1, &z__3, csq, snq, &r__);
	    }

	    *csu = csl;
	    z__2.r = -d1.r, z__2.i = -d1.i;
	    z__1.r = snl * z__2.r, z__1.i = snl * z__2.i;
	    snu->r = z__1.r, snu->i = z__1.i;
	    *csv = csr;
	    z__2.r = -d1.r, z__2.i = -d1.i;
	    z__1.r = snr * z__2.r, z__1.i = snr * z__2.i;
	    snv->r = z__1.r, snv->i = z__1.i;

	} else {

/*           Compute the (2,1) and (2,2) elements of U'*A and V'*B, */
/*           and (2,2) element of |U|'*|A| and |V|'*|B|. */

	    d_cnjg(&z__4, &d1);
	    z__3.r = -z__4.r, z__3.i = -z__4.i;
	    z__2.r = snl * z__3.r, z__2.i = snl * z__3.i;
	    z__1.r = *a1 * z__2.r, z__1.i = *a1 * z__2.i;
	    ua21.r = z__1.r, ua21.i = z__1.i;
	    d_cnjg(&z__5, &d1);
	    z__4.r = -z__5.r, z__4.i = -z__5.i;
	    z__3.r = snl * z__4.r, z__3.i = snl * z__4.i;
	    z__2.r = z__3.r * a2->r - z__3.i * a2->i, z__2.i = z__3.r * a2->i 
		    + z__3.i * a2->r;
	    d__1 = csl * *a3;
	    z__1.r = z__2.r + d__1, z__1.i = z__2.i;
	    ua22.r = z__1.r, ua22.i = z__1.i;

	    d_cnjg(&z__4, &d1);
	    z__3.r = -z__4.r, z__3.i = -z__4.i;
	    z__2.r = snr * z__3.r, z__2.i = snr * z__3.i;
	    z__1.r = *b1 * z__2.r, z__1.i = *b1 * z__2.i;
	    vb21.r = z__1.r, vb21.i = z__1.i;
	    d_cnjg(&z__5, &d1);
	    z__4.r = -z__5.r, z__4.i = -z__5.i;
	    z__3.r = snr * z__4.r, z__3.i = snr * z__4.i;
	    z__2.r = z__3.r * b2->r - z__3.i * b2->i, z__2.i = z__3.r * b2->i 
		    + z__3.i * b2->r;
	    d__1 = csr * *b3;
	    z__1.r = z__2.r + d__1, z__1.i = z__2.i;
	    vb22.r = z__1.r, vb22.i = z__1.i;

	    aua22 = ABS(snl) * ((d__1 = a2->r, ABS(d__1)) + (d__2 = d_imag(a2)
		    , ABS(d__2))) + ABS(csl) * ABS(*a3);
	    avb22 = ABS(snr) * ((d__1 = b2->r, ABS(d__1)) + (d__2 = d_imag(b2)
		    , ABS(d__2))) + ABS(csr) * ABS(*b3);

/*           zero (2,2) elements of U'*A and V'*B, and then swap. */

	    if ((d__1 = ua21.r, ABS(d__1)) + (d__2 = d_imag(&ua21), ABS(d__2))
		     + ((d__3 = ua22.r, ABS(d__3)) + (d__4 = d_imag(&ua22), 
		    ABS(d__4))) == 0.) {
		d_cnjg(&z__2, &vb21);
		z__1.r = -z__2.r, z__1.i = -z__2.i;
		d_cnjg(&z__3, &vb22);
		zlartg_(&z__1, &z__3, csq, snq, &r__);
	    } else if ((d__1 = vb21.r, ABS(d__1)) + (d__2 = d_imag(&vb21), 
		    ABS(d__2)) + z_abs(&vb22) == 0.) {
		d_cnjg(&z__2, &ua21);
		z__1.r = -z__2.r, z__1.i = -z__2.i;
		d_cnjg(&z__3, &ua22);
		zlartg_(&z__1, &z__3, csq, snq, &r__);
	    } else if (aua22 / ((d__1 = ua21.r, ABS(d__1)) + (d__2 = d_imag(&
		    ua21), ABS(d__2)) + ((d__3 = ua22.r, ABS(d__3)) + (d__4 = 
		    d_imag(&ua22), ABS(d__4)))) <= avb22 / ((d__5 = vb21.r, 
		    ABS(d__5)) + (d__6 = d_imag(&vb21), ABS(d__6)) + ((d__7 = 
		    vb22.r, ABS(d__7)) + (d__8 = d_imag(&vb22), ABS(d__8))))) 
		    {
		d_cnjg(&z__2, &ua21);
		z__1.r = -z__2.r, z__1.i = -z__2.i;
		d_cnjg(&z__3, &ua22);
		zlartg_(&z__1, &z__3, csq, snq, &r__);
	    } else {
		d_cnjg(&z__2, &vb21);
		z__1.r = -z__2.r, z__1.i = -z__2.i;
		d_cnjg(&z__3, &vb22);
		zlartg_(&z__1, &z__3, csq, snq, &r__);
	    }

	    *csu = snl;
	    z__1.r = csl * d1.r, z__1.i = csl * d1.i;
	    snu->r = z__1.r, snu->i = z__1.i;
	    *csv = snr;
	    z__1.r = csr * d1.r, z__1.i = csr * d1.i;
	    snv->r = z__1.r, snv->i = z__1.i;

	}

    } else {

/*        Input matrices A and B are lower triangular matrices */

/*        Form matrix C = A*adj(B) = ( a 0 ) */
/*                                   ( c d ) */

	a = *a1 * *b3;
	d__ = *a3 * *b1;
	z__2.r = *b3 * a2->r, z__2.i = *b3 * a2->i;
	z__3.r = *a3 * b2->r, z__3.i = *a3 * b2->i;
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	c__.r = z__1.r, c__.i = z__1.i;
	fc = z_abs(&c__);

/*        Transform complex 2-by-2 matrix C to float matrix by unitary */
/*        diagonal matrix diag(d1,1). */

	d1.r = 1., d1.i = 0.;
	if (fc != 0.) {
	    z__1.r = c__.r / fc, z__1.i = c__.i / fc;
	    d1.r = z__1.r, d1.i = z__1.i;
	}

/*        The SVD of float 2 by 2 triangular C */

/*         ( CSL -SNL )*( A 0 )*(  CSR  SNR ) = ( R 0 ) */
/*         ( SNL  CSL ) ( C D ) ( -SNR  CSR )   ( 0 T ) */

	dlasv2_(&a, &fc, &d__, &s1, &s2, &snr, &csr, &snl, &csl);

	if (ABS(csr) >= ABS(snr) || ABS(csl) >= ABS(snl)) {

/*           Compute the (2,1) and (2,2) elements of U'*A and V'*B, */
/*           and (2,1) element of |U|'*|A| and |V|'*|B|. */

	    z__4.r = -d1.r, z__4.i = -d1.i;
	    z__3.r = snr * z__4.r, z__3.i = snr * z__4.i;
	    z__2.r = *a1 * z__3.r, z__2.i = *a1 * z__3.i;
	    z__5.r = csr * a2->r, z__5.i = csr * a2->i;
	    z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
	    ua21.r = z__1.r, ua21.i = z__1.i;
	    ua22r = csr * *a3;

	    z__4.r = -d1.r, z__4.i = -d1.i;
	    z__3.r = snl * z__4.r, z__3.i = snl * z__4.i;
	    z__2.r = *b1 * z__3.r, z__2.i = *b1 * z__3.i;
	    z__5.r = csl * b2->r, z__5.i = csl * b2->i;
	    z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
	    vb21.r = z__1.r, vb21.i = z__1.i;
	    vb22r = csl * *b3;

	    aua21 = ABS(snr) * ABS(*a1) + ABS(csr) * ((d__1 = a2->r, ABS(d__1)
		    ) + (d__2 = d_imag(a2), ABS(d__2)));
	    avb21 = ABS(snl) * ABS(*b1) + ABS(csl) * ((d__1 = b2->r, ABS(d__1)
		    ) + (d__2 = d_imag(b2), ABS(d__2)));

/*           zero (2,1) elements of U'*A and V'*B. */

	    if ((d__1 = ua21.r, ABS(d__1)) + (d__2 = d_imag(&ua21), ABS(d__2))
		     + ABS(ua22r) == 0.) {
		z__1.r = vb22r, z__1.i = 0.;
		zlartg_(&z__1, &vb21, csq, snq, &r__);
	    } else if ((d__1 = vb21.r, ABS(d__1)) + (d__2 = d_imag(&vb21), 
		    ABS(d__2)) + ABS(vb22r) == 0.) {
		z__1.r = ua22r, z__1.i = 0.;
		zlartg_(&z__1, &ua21, csq, snq, &r__);
	    } else if (aua21 / ((d__1 = ua21.r, ABS(d__1)) + (d__2 = d_imag(&
		    ua21), ABS(d__2)) + ABS(ua22r)) <= avb21 / ((d__3 = 
		    vb21.r, ABS(d__3)) + (d__4 = d_imag(&vb21), ABS(d__4)) + 
		    ABS(vb22r))) {
		z__1.r = ua22r, z__1.i = 0.;
		zlartg_(&z__1, &ua21, csq, snq, &r__);
	    } else {
		z__1.r = vb22r, z__1.i = 0.;
		zlartg_(&z__1, &vb21, csq, snq, &r__);
	    }

	    *csu = csr;
	    d_cnjg(&z__3, &d1);
	    z__2.r = -z__3.r, z__2.i = -z__3.i;
	    z__1.r = snr * z__2.r, z__1.i = snr * z__2.i;
	    snu->r = z__1.r, snu->i = z__1.i;
	    *csv = csl;
	    d_cnjg(&z__3, &d1);
	    z__2.r = -z__3.r, z__2.i = -z__3.i;
	    z__1.r = snl * z__2.r, z__1.i = snl * z__2.i;
	    snv->r = z__1.r, snv->i = z__1.i;

	} else {

/*           Compute the (1,1) and (1,2) elements of U'*A and V'*B, */
/*           and (1,1) element of |U|'*|A| and |V|'*|B|. */

	    d__1 = csr * *a1;
	    d_cnjg(&z__4, &d1);
	    z__3.r = snr * z__4.r, z__3.i = snr * z__4.i;
	    z__2.r = z__3.r * a2->r - z__3.i * a2->i, z__2.i = z__3.r * a2->i 
		    + z__3.i * a2->r;
	    z__1.r = d__1 + z__2.r, z__1.i = z__2.i;
	    ua11.r = z__1.r, ua11.i = z__1.i;
	    d_cnjg(&z__3, &d1);
	    z__2.r = snr * z__3.r, z__2.i = snr * z__3.i;
	    z__1.r = *a3 * z__2.r, z__1.i = *a3 * z__2.i;
	    ua12.r = z__1.r, ua12.i = z__1.i;

	    d__1 = csl * *b1;
	    d_cnjg(&z__4, &d1);
	    z__3.r = snl * z__4.r, z__3.i = snl * z__4.i;
	    z__2.r = z__3.r * b2->r - z__3.i * b2->i, z__2.i = z__3.r * b2->i 
		    + z__3.i * b2->r;
	    z__1.r = d__1 + z__2.r, z__1.i = z__2.i;
	    vb11.r = z__1.r, vb11.i = z__1.i;
	    d_cnjg(&z__3, &d1);
	    z__2.r = snl * z__3.r, z__2.i = snl * z__3.i;
	    z__1.r = *b3 * z__2.r, z__1.i = *b3 * z__2.i;
	    vb12.r = z__1.r, vb12.i = z__1.i;

	    aua11 = ABS(csr) * ABS(*a1) + ABS(snr) * ((d__1 = a2->r, ABS(d__1)
		    ) + (d__2 = d_imag(a2), ABS(d__2)));
	    avb11 = ABS(csl) * ABS(*b1) + ABS(snl) * ((d__1 = b2->r, ABS(d__1)
		    ) + (d__2 = d_imag(b2), ABS(d__2)));

/*           zero (1,1) elements of U'*A and V'*B, and then swap. */

	    if ((d__1 = ua11.r, ABS(d__1)) + (d__2 = d_imag(&ua11), ABS(d__2))
		     + ((d__3 = ua12.r, ABS(d__3)) + (d__4 = d_imag(&ua12), 
		    ABS(d__4))) == 0.) {
		zlartg_(&vb12, &vb11, csq, snq, &r__);
	    } else if ((d__1 = vb11.r, ABS(d__1)) + (d__2 = d_imag(&vb11), 
		    ABS(d__2)) + ((d__3 = vb12.r, ABS(d__3)) + (d__4 = d_imag(
		    &vb12), ABS(d__4))) == 0.) {
		zlartg_(&ua12, &ua11, csq, snq, &r__);
	    } else if (aua11 / ((d__1 = ua11.r, ABS(d__1)) + (d__2 = d_imag(&
		    ua11), ABS(d__2)) + ((d__3 = ua12.r, ABS(d__3)) + (d__4 = 
		    d_imag(&ua12), ABS(d__4)))) <= avb11 / ((d__5 = vb11.r, 
		    ABS(d__5)) + (d__6 = d_imag(&vb11), ABS(d__6)) + ((d__7 = 
		    vb12.r, ABS(d__7)) + (d__8 = d_imag(&vb12), ABS(d__8))))) 
		    {
		zlartg_(&ua12, &ua11, csq, snq, &r__);
	    } else {
		zlartg_(&vb12, &vb11, csq, snq, &r__);
	    }

	    *csu = snr;
	    d_cnjg(&z__2, &d1);
	    z__1.r = csr * z__2.r, z__1.i = csr * z__2.i;
	    snu->r = z__1.r, snu->i = z__1.i;
	    *csv = snl;
	    d_cnjg(&z__2, &d1);
	    z__1.r = csl * z__2.r, z__1.i = csl * z__2.i;
	    snv->r = z__1.r, snv->i = z__1.i;

	}

    }

    return 0;

/*     End of ZLAGS2 */

} /* zlags2_ */
