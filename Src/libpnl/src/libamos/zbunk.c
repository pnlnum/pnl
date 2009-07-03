#include "amos.h"

/* Copyright, Donald E. Amos: sandia national laboratories
 *            from slatec library or amos library.
 *
 * issued by sandia laboratories, a prime contractor to the 
 * united states department of energy 
 * notice:
 * this report was prepared as an account of work sponsored by the 
 * united states government.  neither the united states nor the 
 * united states department of energy, nor any of their 
 * employees, nor any of their contractors, subcontractors, or their 
 * employees, makes any warranty, express or implied, or assumes any 
 * legal liability or responsibility for the accuracy, completeness 
 * or usefulness of any information, apparatus, product or process 
 * disclosed, or represents that its use would not infringe 
 * privately owned rights. 
 *
 *
 * this code has been approved for unlimited release. 
 */



/*
 *    zbunk computes the k bessel function for fnu.gt.fnul. 
 *    according to the uniform asymptotic expansion for k(fnu,z) 
 *    in zunk1 and the expansion for h(2,fnu,z) in zunk2 
 * 
 */

/* Subroutine */ int
amos_zbunk (double *zr, double *zi, double *fnu, int *kode, int *mr,const int *n,
	    double *yr, double *yi, int *nz, double *tol, double *elim,
	    double *alim)
{
  double ax, ay;

  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  *nz = 0;
  ax = fabs (*zr) * 1.7321;
  ay = fabs (*zi);
  if (ay > ax)
    {
      goto L10;
    }
  /*----------------------------------------------------------------------- 
   *    ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN 
   *    -PI/3.LE.ARG(Z).LE.PI/3 
   *----------------------------------------------------------------------- 
   */
  amos_zunk1 (zr, zi, fnu, kode, mr, n, &yr[1], &yi[1], nz, tol, elim, alim);
  goto L20;
 L10:
  /*----------------------------------------------------------------------- 
   *    ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU 
   *    APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I 
   *    AND HPI=PI/2 
   *----------------------------------------------------------------------- 
   */
  amos_zunk2 (zr, zi, fnu, kode, mr, n, &yr[1], &yi[1], nz, tol, elim, alim);
 L20:
  return 0;
}				/* zbunk_ */
