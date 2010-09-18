/*
 * This file contains the routine dqwgts from the quadpack distribution
 * available on www.netlib.org
 *
 * This routine has been translated from Fortran using f2c and then slightly
 * modified by Jérôme Lelong for integration into Pnl. All Quadpack names have
 * been prefixed by pnl to avoid possible clash in names. Calls to d1mach have
 * been replaced by calls to dlamch coming from Lapack
 */

#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"

double pnl_dqwgts(double *x, double *a, double *b, double *
	alfa, double *beta, int *integr)
{
    /* System generated locals */
    double ret_val;

    /* Local variables */
    double xma, bmx;

/* ***begin prologue  dqwgts */
/* ***refer to dqk15w */
/* ***routines called  (none) */
/* ***revision date  810101   (yymmdd) */
/* ***keywords  weight function, algebraico-logarithmic */
/*             end-point singularities */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  this function subprogram is used together with the */
/*            routine dqaws and defines the weight function. */
/* ***end prologue  dqwgts */

/* ***first executable statement  dqwgts */
    xma = *x - *a;
    bmx = *b - *x;
    ret_val = pow(xma, *alfa) * pow(bmx, *beta);
    switch (*integr) {
	case 1:  goto L40;
	case 2:  goto L10;
	case 3:  goto L20;
	case 4:  goto L30;
    }
L10:
    ret_val *= log(xma);
    goto L40;
L20:
    ret_val *= log(bmx);
    goto L40;
L30:
    ret_val = ret_val * log(xma) * log(bmx);
L40:
    return ret_val;
} /* pnl_dqwgts */

