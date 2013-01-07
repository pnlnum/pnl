/*  
 *  This file is a collection of wrappers around the
 *  Amos Fortran library of functions that take complex
 *  variables (see www.netlib.org) so that they can be called from
 *  the cephes library of corresponding name but work with complex
 *  arguments.
 * 
 *  Author:  Travis E. Oliphant
 *            oliphant.travis@altavista.net
 *
 */

#ifndef _AMOS_WRAPPERS_H
#define _AMOS_WRAPPERS_H

#include "pnl/pnl_complex.h" /* just for dcomplex ! */

/* already defined in math.h */
/* #define DOMAIN		1	|+ argument domain error +| */
/* #define SING		2	|+ argument singularity +| */
/* #define OVERFLOW	3	|+ overflow range error +| */
/* #define UNDERFLOW	4	|+ underflow range error +| */
/* #define TLOSS		5	|+ total loss of precision +| */
/* #define PLOSS		6	|+ partial loss of precision +| */
#define TOOMANY         7       /* too many iterations */
#define MAXITER        500


static int ierr_to_mtherr( int nz, int ierr);
static int mtherr(char *name, int code); /* from libcephes */

#define DO_MTHERR(name) if (nz !=0 || ierr !=0) mtherr(name, ierr_to_mtherr(nz, ierr))
#define CADDR(z) (double *)&z.r, (double*)&z.i
#define F2C_CST(z) (double *)&z->r, (double *)&z->i

/* extern int cairy_wrap(dcomplex z, dcomplex *ai, dcomplex *aip, dcomplex *bi, dcomplex *bip);
 * extern int cairye_wrap(dcomplex z, dcomplex *ai, dcomplex *aip, dcomplex *bi, dcomplex *bip); */


#endif



  







