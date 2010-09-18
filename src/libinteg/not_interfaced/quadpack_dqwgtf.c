/*
 * This file contains the routine dqwgtf from the quadpack distribution
 * available on www.netlib.org
 *
 * This routine has been translated from Fortran using f2c and then slightly
 * modified by Jérôme Lelong for integration into Pnl. All Quadpack names have
 * been prefixed by pnl to avoid possible clash in names. Calls to d1mach have
 * been replaced by calls to dlamch coming from Lapack
 */

#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"

double pnl_dqwgtf(double *x, double *omega, double *p2, 
	double *p3, double *p4, int *integr)
{
    /* System generated locals */
    double ret_val;

    /* Local variables */
    double omx;

/* ***begin prologue  dqwgtf */
/* ***refer to   dqk15w */
/* ***routines called  (none) */
/* ***revision date 810101   (yymmdd) */
/* ***keywords  cos or sin in weight function */
/* ***author  piessens,robert, appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. * progr. div. - k.u.leuven */
/* ***end prologue  dqwgtf */

/* ***first executable statement  dqwgtf */
    omx = *omega * *x;
    switch (*integr) {
	case 1:  goto L10;
	case 2:  goto L20;
    }
L10:
    ret_val = cos(omx);
    goto L30;
L20:
    ret_val = sin(omx);
L30:
    return ret_val;
} /* pnl_dqwgtf */

