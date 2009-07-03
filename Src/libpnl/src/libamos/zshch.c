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
 *    zshch computes the complex hyperbolic functions csh=sinh(x+i*y) 
 *    and cch=cosh(x+i*y), where i**2=-1. 
 *
 * 
 */

int amos_zshch (double *zr, double *zi, double *cshr, double *cshi, double *cchr,
		double *cchi)
{

  /* Local variables */
  double ch, cn, sh, sn;

  sh = sinh (*zr);
  ch = cosh (*zr);
  sn = sin (*zi);
  cn = cos (*zi);
  *cshr = sh * cn;
  *cshi = ch * sn;
  *cchr = ch * cn;
  *cchi = sh * sn;
  return 0;
}		
