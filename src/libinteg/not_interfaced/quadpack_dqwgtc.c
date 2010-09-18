/*
 * This file contains the routine dqwgtc from the quadpack distribution
 * available on www.netlib.org
 *
 * This routine has been translated from Fortran using f2c and then slightly
 * modified by Jérôme Lelong for integration into Pnl. All Quadpack names have
 * been prefixed by pnl to avoid possible clash in names. Calls to d1mach have
 * been replaced by calls to dlamch coming from Lapack
 */

#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"

double pnl_dqwgtc(double *x, double *c__, double *p2, double 
	*p3, double *p4, int *kp)
{
    /* System generated locals */
    double ret_val;

/* ***begin prologue  dqwgtc */
/* ***refer to dqk15w */
/* ***routines called  (none) */
/* ***revision date  810101   (yymmdd) */
/* ***keywords  weight function, cauchy principal value */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  this function subprogram is used together with the */
/*            routine qawc and defines the weight function. */
/* ***end prologue  dqwgtc */

/* ***first executable statement  dqwgtc */
    ret_val = 1. / (*x - *c__);
    return ret_val;
} /* pnl_dqwgtc */

