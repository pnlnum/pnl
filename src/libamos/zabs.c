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


/**
 * @param zr: a double pointer 
 * @param zi: a double pointer 
 * 
 * computes the absolute value or magnitude of a double 
 * precision complex variable cmplx(zr,zi) 
 * 
 * 
 * @return a double 
 **/

double amos_azabs (const double *zr,const double *zi)
{
  double q, s, u, v;
  u = fabs (*zr);
  v = fabs (*zi);
  s = u + v;
  /*
   *    s*1.0d0 makes an unnormalized underflow on cdc machines into a 
   *    true floating zero 
   */
  s *= 1.;
  if (s == 0.)
    {
      return 0;
    }
  if (u > v)
    {
      q = v / u;
      return  u * sqrt (q * q + 1.);
    }
  q = u / v;
  return v * sqrt (q * q + 1.);
}

