/* dlags2.f -- translated by f2c (version 20061008).
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

 int dlags2_(int *upper, double *a1, double *a2, 
	double *a3, double *b1, double *b2, double *b3, 
	double *csu, double *snu, double *csv, double *snv, 
	double *csq, double *snq)
{
    /* System generated locals */
    double d__1;

    /* Local variables */
    double a, b, c__, d__, r__, s1, s2, ua11, ua12, ua21, ua22, vb11, 
	    vb12, vb21, vb22, csl, csr, snl, snr, aua11, aua12, aua21, aua22, 
	    avb11, avb12, avb21, avb22, ua11r, ua22r, vb11r, vb22r;
    extern  int dlasv2_(double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *), dlartg_(double *, 
	    double *, double *, double *, double *);


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAGS2 computes 2-by-2 orthogonal matrices U, V and Q, such */
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

/*  The rows of the transformed A and B are parallel, where */

/*    U = (  CSU  SNU ), V = (  CSV SNV ), Q = (  CSQ   SNQ ) */
/*        ( -SNU  CSU )      ( -SNV CSV )      ( -SNQ   CSQ ) */

/*  Z' denotes the transpose of Z. */


/*  Arguments */
/*  ========= */

/*  UPPER   (input) LOGICAL */
/*          = .TRUE.: the input matrices A and B are upper triangular. */
/*          = .FALSE.: the input matrices A and B are lower triangular. */

/*  A1      (input) DOUBLE PRECISION */
/*  A2      (input) DOUBLE PRECISION */
/*  A3      (input) DOUBLE PRECISION */
/*          On entry, A1, A2 and A3 are elements of the input 2-by-2 */
/*          upper (lower) triangular matrix A. */

/*  B1      (input) DOUBLE PRECISION */
/*  B2      (input) DOUBLE PRECISION */
/*  B3      (input) DOUBLE PRECISION */
/*          On entry, B1, B2 and B3 are elements of the input 2-by-2 */
/*          upper (lower) triangular matrix B. */

/*  CSU     (output) DOUBLE PRECISION */
/*  SNU     (output) DOUBLE PRECISION */
/*          The desired orthogonal matrix U. */

/*  CSV     (output) DOUBLE PRECISION */
/*  SNV     (output) DOUBLE PRECISION */
/*          The desired orthogonal matrix V. */

/*  CSQ     (output) DOUBLE PRECISION */
/*  SNQ     (output) DOUBLE PRECISION */
/*          The desired orthogonal matrix Q. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    if (*upper) {

/*        Input matrices A and B are upper triangular matrices */

/*        Form matrix C = A*adj(B) = ( a b ) */
/*                                   ( 0 d ) */

	a = *a1 * *b3;
	d__ = *a3 * *b1;
	b = *a2 * *b1 - *a1 * *b2;

/*        The SVD of float 2-by-2 triangular C */

/*         ( CSL -SNL )*( A B )*(  CSR  SNR ) = ( R 0 ) */
/*         ( SNL  CSL ) ( 0 D ) ( -SNR  CSR )   ( 0 T ) */

	dlasv2_(&a, &b, &d__, &s1, &s2, &snr, &csr, &snl, &csl);

	if (ABS(csl) >= ABS(snl) || ABS(csr) >= ABS(snr)) {

/*           Compute the (1,1) and (1,2) elements of U'*A and V'*B, */
/*           and (1,2) element of |U|'*|A| and |V|'*|B|. */

	    ua11r = csl * *a1;
	    ua12 = csl * *a2 + snl * *a3;

	    vb11r = csr * *b1;
	    vb12 = csr * *b2 + snr * *b3;

	    aua12 = ABS(csl) * ABS(*a2) + ABS(snl) * ABS(*a3);
	    avb12 = ABS(csr) * ABS(*b2) + ABS(snr) * ABS(*b3);

/*           zero (1,2) elements of U'*A and V'*B */

	    if (ABS(ua11r) + ABS(ua12) != 0.) {
		if (aua12 / (ABS(ua11r) + ABS(ua12)) <= avb12 / (ABS(vb11r) + 
			ABS(vb12))) {
		    d__1 = -ua11r;
		    dlartg_(&d__1, &ua12, csq, snq, &r__);
		} else {
		    d__1 = -vb11r;
		    dlartg_(&d__1, &vb12, csq, snq, &r__);
		}
	    } else {
		d__1 = -vb11r;
		dlartg_(&d__1, &vb12, csq, snq, &r__);
	    }

	    *csu = csl;
	    *snu = -snl;
	    *csv = csr;
	    *snv = -snr;

	} else {

/*           Compute the (2,1) and (2,2) elements of U'*A and V'*B, */
/*           and (2,2) element of |U|'*|A| and |V|'*|B|. */

	    ua21 = -snl * *a1;
	    ua22 = -snl * *a2 + csl * *a3;

	    vb21 = -snr * *b1;
	    vb22 = -snr * *b2 + csr * *b3;

	    aua22 = ABS(snl) * ABS(*a2) + ABS(csl) * ABS(*a3);
	    avb22 = ABS(snr) * ABS(*b2) + ABS(csr) * ABS(*b3);

/*           zero (2,2) elements of U'*A and V'*B, and then swap. */

	    if (ABS(ua21) + ABS(ua22) != 0.) {
		if (aua22 / (ABS(ua21) + ABS(ua22)) <= avb22 / (ABS(vb21) + 
			ABS(vb22))) {
		    d__1 = -ua21;
		    dlartg_(&d__1, &ua22, csq, snq, &r__);
		} else {
		    d__1 = -vb21;
		    dlartg_(&d__1, &vb22, csq, snq, &r__);
		}
	    } else {
		d__1 = -vb21;
		dlartg_(&d__1, &vb22, csq, snq, &r__);
	    }

	    *csu = snl;
	    *snu = csl;
	    *csv = snr;
	    *snv = csr;

	}

    } else {

/*        Input matrices A and B are lower triangular matrices */

/*        Form matrix C = A*adj(B) = ( a 0 ) */
/*                                   ( c d ) */

	a = *a1 * *b3;
	d__ = *a3 * *b1;
	c__ = *a2 * *b3 - *a3 * *b2;

/*        The SVD of float 2-by-2 triangular C */

/*         ( CSL -SNL )*( A 0 )*(  CSR  SNR ) = ( R 0 ) */
/*         ( SNL  CSL ) ( C D ) ( -SNR  CSR )   ( 0 T ) */

	dlasv2_(&a, &c__, &d__, &s1, &s2, &snr, &csr, &snl, &csl);

	if (ABS(csr) >= ABS(snr) || ABS(csl) >= ABS(snl)) {

/*           Compute the (2,1) and (2,2) elements of U'*A and V'*B, */
/*           and (2,1) element of |U|'*|A| and |V|'*|B|. */

	    ua21 = -snr * *a1 + csr * *a2;
	    ua22r = csr * *a3;

	    vb21 = -snl * *b1 + csl * *b2;
	    vb22r = csl * *b3;

	    aua21 = ABS(snr) * ABS(*a1) + ABS(csr) * ABS(*a2);
	    avb21 = ABS(snl) * ABS(*b1) + ABS(csl) * ABS(*b2);

/*           zero (2,1) elements of U'*A and V'*B. */

	    if (ABS(ua21) + ABS(ua22r) != 0.) {
		if (aua21 / (ABS(ua21) + ABS(ua22r)) <= avb21 / (ABS(vb21) + 
			ABS(vb22r))) {
		    dlartg_(&ua22r, &ua21, csq, snq, &r__);
		} else {
		    dlartg_(&vb22r, &vb21, csq, snq, &r__);
		}
	    } else {
		dlartg_(&vb22r, &vb21, csq, snq, &r__);
	    }

	    *csu = csr;
	    *snu = -snr;
	    *csv = csl;
	    *snv = -snl;

	} else {

/*           Compute the (1,1) and (1,2) elements of U'*A and V'*B, */
/*           and (1,1) element of |U|'*|A| and |V|'*|B|. */

	    ua11 = csr * *a1 + snr * *a2;
	    ua12 = snr * *a3;

	    vb11 = csl * *b1 + snl * *b2;
	    vb12 = snl * *b3;

	    aua11 = ABS(csr) * ABS(*a1) + ABS(snr) * ABS(*a2);
	    avb11 = ABS(csl) * ABS(*b1) + ABS(snl) * ABS(*b2);

/*           zero (1,1) elements of U'*A and V'*B, and then swap. */

	    if (ABS(ua11) + ABS(ua12) != 0.) {
		if (aua11 / (ABS(ua11) + ABS(ua12)) <= avb11 / (ABS(vb11) + 
			ABS(vb12))) {
		    dlartg_(&ua12, &ua11, csq, snq, &r__);
		} else {
		    dlartg_(&vb12, &vb11, csq, snq, &r__);
		}
	    } else {
		dlartg_(&vb12, &vb11, csq, snq, &r__);
	    }

	    *csu = snr;
	    *snu = csr;
	    *csv = snl;
	    *snv = csl;

	}

    }

    return 0;

/*     End of DLAGS2 */

} /* dlags2_ */
