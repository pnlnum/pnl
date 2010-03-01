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
 *    double precision complex square root, b=csqrt(a) 
 *
 */

int amos_azsqrt (const double *ar,const double *ai, double *br, double *bi)
{
  static const double drt = .7071067811865475244008443621;
  static const double dpi = 3.141592653589793238462643383;
  
  /* Local variables */
  double zm, dtheta;

  zm = amos_azabs (ar, ai);
  zm = sqrt (zm);
  if (*ar == 0.)
    {
      goto L10;
    }
  if (*ai == 0.)
    {
      goto L20;
    }
  dtheta = atan (*ai / *ar);
  if (dtheta <= 0.)
    {
      goto L40;
    }
  if (*ar < 0.)
    {
      dtheta -= dpi;
    }
  goto L50;
 L10:
  if (*ai > 0.)
    {
      goto L60;
    }
  if (*ai < 0.)
    {
      goto L70;
    }
  *br = 0.;
  *bi = 0.;
  return 0;
 L20:
  if (*ar > 0.)
    {
      goto L30;
    }
  *br = 0.;
  *bi = sqrt ((fabs (*ar)));
  return 0;
 L30:
  *br = sqrt (*ar);
  *bi = 0.;
  return 0;
 L40:
  if (*ar < 0.)
    {
      dtheta += dpi;
    }
 L50:
  dtheta *= .5;
  *br = zm * cos (dtheta);
  *bi = zm * sin (dtheta);
  return 0;
 L60:
  *br = zm * drt;
  *bi = zm * drt;
  return 0;
 L70:
  *br = zm * drt;
  *bi = -zm * drt;
  return 0;
}				/* azsqrt_ */
