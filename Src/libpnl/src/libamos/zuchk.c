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
 *     y enters as a scaled quantity whose magnitude is greater than 
 *     exp(-alim)=ascle=1.0e+3*d1mach(1)/tol. the test is made to see 
 *     if the magnitude of the real or imaginary part would underflow 
 *     when y is scaled (by tol) to its proper value. y is accepted 
 *     if the underflow is at least one precision below the magnitude 
 *     of the largest component; otherwise the phase angle does not have 
 *     absolute accuracy and an underflow is assumed. 
 *
 */

int amos_zuchk (double *yr, double *yi, int *nz, double *ascle, double *tol)
{
  double wi, ss, st, wr;
  *nz = 0;
  wr = fabs (*yr);
  wi = fabs (*yi);
  st = MIN (wr, wi);
  if (st > *ascle)
    {
      return 0;
    }
  ss = MAX (wr, wi);
  st /= *tol;
  if (ss < st)
    {
      *nz = 1;
    }
  return 0;
}	

