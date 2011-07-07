/* cla_geamv.f -- translated by f2c (version 20061008).
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

 int cla_geamv__(int *trans, int *m, int *n, float 
	*alpha, complex *a, int *lda, complex *x, int *incx, float *
	beta, float *y, int *incy)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;
    float r__1, r__2;

    /* Builtin functions */
    double r_imag(complex *), r_sign(float *, float *);

    /* Local variables */
    extern int ilatrans_(char *);
    int i__, j;
    int symb_zero__;
    int iy, jx, kx, ky, info;
    float temp;
    int lenx, leny;
    float safe1;
    extern double slamch_(char *);
    extern  int xerbla_(char *, int *);


/*     -- LAPACK routine (version 3.2)                                 -- */
/*     -- Contributed by James Demmel, Deaglan Halligan, Yozo Hida and -- */
/*     -- Jason Riedy of Univ. of California Berkeley.                 -- */
/*     -- November 2008                                                -- */

/*     -- LAPACK is a software package provided by Univ. of Tennessee, -- */
/*     -- Univ. of California Berkeley and NAG Ltd.                    -- */

/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  CLA_GEAMV  performs one of the matrix-vector operations */

/*          y := alpha*ABS(A)*ABS(x) + beta*ABS(y), */
/*     or   y := alpha*ABS(A)'*ABS(x) + beta*ABS(y), */

/*  where alpha and beta are scalars, x and y are vectors and A is an */
/*  m by n matrix. */

/*  This function is primarily used in calculating error bounds. */
/*  To protect against underflow during evaluation, components in */
/*  the resulting vector are perturbed away from zero by (N+1) */
/*  times the underflow threshold.  To prevent unnecessarily large */
/*  errors for block-structure embedded in general matrices, */
/*  "symbolically" zero components are not perturbed.  A zero */
/*  entry is considered "symbolic" if all multiplications involved */
/*  in computing that entry have at least one zero multiplicand. */

/*  Parameters */
/*  ========== */

/*  TRANS  - INTEGER */
/*           On entry, TRANS specifies the operation to be performed as */
/*           follows: */

/*             BLAS_NO_TRANS      y := alpha*ABS(A)*ABS(x) + beta*ABS(y) */
/*             BLAS_TRANS         y := alpha*ABS(A')*ABS(x) + beta*ABS(y) */
/*             BLAS_CONJ_TRANS    y := alpha*ABS(A')*ABS(x) + beta*ABS(y) */

/*           Unchanged on exit. */

/*  M      - INTEGER */
/*           On entry, M specifies the number of rows of the matrix A. */
/*           M must be at least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER */
/*           On entry, N specifies the number of columns of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - REAL */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - COMPLEX          array of DIMENSION ( LDA, n ) */
/*           Before entry, the leading m by n part of the array A must */
/*           contain the matrix of coefficients. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. LDA must be at least */
/*           MAX( 1, m ). */
/*           Unchanged on exit. */

/*  X      - COMPLEX          array of DIMENSION at least */
/*           ( 1 + ( n - 1 )*ABS( INCX ) ) when TRANS = 'N' or 'n' */
/*           and at least */
/*           ( 1 + ( m - 1 )*ABS( INCX ) ) otherwise. */
/*           Before entry, the incremented array X must contain the */
/*           vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  BETA   - REAL */
/*           On entry, BETA specifies the scalar beta. When BETA is */
/*           supplied as zero then Y need not be set on input. */
/*           Unchanged on exit. */

/*  Y      - REAL             array of DIMENSION at least */
/*           ( 1 + ( m - 1 )*ABS( INCY ) ) when TRANS = 'N' or 'n' */
/*           and at least */
/*           ( 1 + ( n - 1 )*ABS( INCY ) ) otherwise. */
/*           Before entry with BETA non-zero, the incremented array Y */
/*           must contain the vector y. On exit, Y is overwritten by the */
/*           updated vector y. */

/*  INCY   - INTEGER */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*     .. */
/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function Definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (! (*trans == ilatrans_("N") || *trans == ilatrans_("T") || *trans == ilatrans_("C"))) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*lda < MAX(1,*m)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("CLA_GEAMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == 0.f && *beta == 1.f) {
	return 0;
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

    if (*trans == ilatrans_("N")) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }

/*     Set SAFE1 essentially to be the underflow threshold times the */
/*     number of additions in each row. */

    safe1 = slamch_("Safe minimum");
    safe1 = (*n + 1) * safe1;

/*     Form  y := alpha*ABS(A)*ABS(x) + beta*ABS(y). */

/*     The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to */
/*     the inexact flag.  Still doesn't help change the iteration order */
/*     to per-column. */

    iy = ky;
    if (*incx == 1) {
	i__1 = leny;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (*beta == 0.f) {
		symb_zero__ = TRUE;
		y[iy] = 0.f;
	    } else if (y[iy] == 0.f) {
		symb_zero__ = TRUE;
	    } else {
		symb_zero__ = FALSE;
		y[iy] = *beta * (r__1 = y[iy], ABS(r__1));
	    }
	    if (*alpha != 0.f) {
		i__2 = lenx;
		for (j = 1; j <= i__2; ++j) {
		    if (*trans == ilatrans_("N")) {
			i__3 = i__ + j * a_dim1;
			temp = (r__1 = a[i__3].r, ABS(r__1)) + (r__2 = 
				r_imag(&a[i__ + j * a_dim1]), ABS(r__2));
		    } else {
			i__3 = j + i__ * a_dim1;
			temp = (r__1 = a[i__3].r, ABS(r__1)) + (r__2 = 
				r_imag(&a[j + i__ * a_dim1]), ABS(r__2));
		    }
		    i__3 = j;
		    symb_zero__ = symb_zero__ && (x[i__3].r == 0.f && x[i__3]
			    .i == 0.f || temp == 0.f);
		    i__3 = j;
		    y[iy] += *alpha * ((r__1 = x[i__3].r, ABS(r__1)) + (r__2 
			    = r_imag(&x[j]), ABS(r__2))) * temp;
		}
	    }
	    if (! symb_zero__) {
		y[iy] += r_sign(&safe1, &y[iy]);
	    }
	    iy += *incy;
	}
    } else {
	i__1 = leny;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (*beta == 0.f) {
		symb_zero__ = TRUE;
		y[iy] = 0.f;
	    } else if (y[iy] == 0.f) {
		symb_zero__ = TRUE;
	    } else {
		symb_zero__ = FALSE;
		y[iy] = *beta * (r__1 = y[iy], ABS(r__1));
	    }
	    if (*alpha != 0.f) {
		jx = kx;
		i__2 = lenx;
		for (j = 1; j <= i__2; ++j) {
		    if (*trans == ilatrans_("N")) {
			i__3 = i__ + j * a_dim1;
			temp = (r__1 = a[i__3].r, ABS(r__1)) + (r__2 = 
				r_imag(&a[i__ + j * a_dim1]), ABS(r__2));
		    } else {
			i__3 = j + i__ * a_dim1;
			temp = (r__1 = a[i__3].r, ABS(r__1)) + (r__2 = 
				r_imag(&a[j + i__ * a_dim1]), ABS(r__2));
		    }
		    i__3 = jx;
		    symb_zero__ = symb_zero__ && (x[i__3].r == 0.f && x[i__3]
			    .i == 0.f || temp == 0.f);
		    i__3 = jx;
		    y[iy] += *alpha * ((r__1 = x[i__3].r, ABS(r__1)) + (r__2 
			    = r_imag(&x[jx]), ABS(r__2))) * temp;
		    jx += *incx;
		}
	    }
	    if (! symb_zero__) {
		y[iy] += r_sign(&safe1, &y[iy]);
	    }
	    iy += *incy;
	}
    }

    return 0;

/*     End of CLA_GEAMV */

} /* cla_geamv__ */
