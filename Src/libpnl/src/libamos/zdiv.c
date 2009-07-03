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
 *    double precision complex divide c=a/b. 
 */

int amos_zdiv (const double *ar,const  double *ai, const double *br, const double *bi, double *cr, double *ci)
{
  double ca, cb, cc, cd, bm;

  bm = 1. / amos_azabs (br, bi);
  cc = *br * bm;
  cd = *bi * bm;
  ca = (*ar * cc + *ai * cd) * bm;
  cb = (*ai * cc - *ar * cd) * bm;
  *cr = ca;
  *ci = cb;
  return 0;
}

