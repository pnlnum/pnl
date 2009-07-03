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
 *    double precision complex logarithm b=clog(a) 
 *    ierr=0,normal return      ierr=1, z=cmplx(0.0,0.0) 
 *
 */

int amos_azlog (const double *ar,const double *ai, double *br, double *bi, int *ierr)
{
  static const double dpi = 3.141592653589793238462643383;
  static const double dhpi = 1.570796326794896619231321696;
  double zm, dtheta;

  *ierr = 0;
  if (*ar == 0.)
    {
      if (*ai == 0.)
	{
	  *ierr = 1;
	  return 0;
	}
      *bi = dhpi;
      *br = log ((fabs (*ai)));
      if (*ai < 0.)
	{
	  *bi = -(*bi);
	}
      return 0;
    }
  if (*ai == 0.)
    {
      if (*ar > 0.)
	{
	  *br = log (*ar);
	  *bi = 0.;
	  return 0;
	}
      *br = log ((fabs (*ar)));
      *bi = dpi;
      return 0;
    }
  dtheta = atan (*ai / *ar);
  if (dtheta <= 0.)
    {
      if (*ar < 0.)
	{
	  dtheta += dpi;
	}
    }
  else 
    {
      if (*ar < 0.)
	{
	  dtheta -= dpi;
	}
    }
  zm = amos_azabs (ar, ai);
  *br = log (zm);
  *bi = dtheta;
  return 0;
}

