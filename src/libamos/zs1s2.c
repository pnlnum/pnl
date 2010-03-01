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
 *    zs1s2 tests for a possible underflow resulting from the 
 *    addition of the i and k functions in the analytic con- 
 *    tinuation formula where s1=k function and s2=i function. 
 *    on kode=1 the i and k functions are different orders of 
 *    magnitude, but for kode=2 they can be of the same order 
 *    of magnitude and the maximum must be at least one 
 *    precision above the underflow limit. 
 *
 */

/* Subroutine */ int
amos_zs1s2 (double *zrr, double *zri, double *s1r, double *s1i, double *s2r,
	    double *s2i, int *nz, double *ascle, double *alim, int *iuf)
{
  /* Initialized data */

  static const double zeror = 0.;
  static const double zeroi = 0.;

  /* Local variables */
  int idum;
  double aa, c1i, as1, as2, c1r, aln, s1di, s1dr;

  *nz = 0;
  as1 = amos_azabs (s1r, s1i);
  as2 = amos_azabs (s2r, s2i);
  if (*s1r == 0. && *s1i == 0.)
    {
      goto L10;
    }
  if (as1 == 0.)
    {
      goto L10;
    }
  aln = -(*zrr) - *zrr + log (as1);
  s1dr = *s1r;
  s1di = *s1i;
  *s1r = zeror;
  *s1i = zeroi;
  as1 = zeror;
  if (aln < -(*alim))
    {
      goto L10;
    }
  amos_azlog (&s1dr, &s1di, &c1r, &c1i, &idum);
  c1r = c1r - *zrr - *zrr;
  c1i = c1i - *zri - *zri;
  amos_azexp (&c1r, &c1i, s1r, s1i);
  as1 = amos_azabs (s1r, s1i);
  ++(*iuf);
 L10:
  aa = MAX (as1, as2);
  if (aa > *ascle)
    {
      return 0;
    }
  *s1r = zeror;
  *s1i = zeroi;
  *s2r = zeror;
  *s2i = zeroi;
  *nz = 1;
  *iuf = 0;
  return 0;
}				/* zs1s2_ */
