/* Minpack Copyright Notice (1999) University of Chicago.  All rights reserved
 * 
 * Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the
 * following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above
 * copyright notice, this list of conditions and the following
 * disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above
 * copyright notice, this list of conditions and the following
 * disclaimer in the documentation and/or other materials
 * provided with the distribution.
 * 
 * 3. The end-user documentation included with the
 * redistribution, if any, must include the following
 * acknowledgment:
 * 
 *    "This product includes software developed by the
 *    University of Chicago, as Operator of Argonne National
 *    Laboratory.
 * 
 * Alternately, this acknowledgment may appear in the software
 * itself, if and wherever such third-party acknowledgments
 * normally appear.
 * 
 * 4. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
 * WITHOUT WARRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
 * UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
 * THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
 * OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
 * OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
 * USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
 * THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
 * DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
 * UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL
 * BE CORRECTED.
 * 
 * 5. LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT
 * HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
 * ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
 * INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
 * ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
 * PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
 * SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
 * (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE,
 * EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
 * POSSIBILITY OF SUCH LOSS OR DAMAGES.
 */

#include <math.h>
#include "cminpack.h"

/* 
 * function enorm
 * 
 * given an n-vector x, this function calculates the
 * euclidean norm of x.
 * 
 * the euclidean norm is computed by accumulating the sum of
 * squares in three different sums. the sums of squares for the
 * small and large components are scaled so that no overflows
 * occur. non-destructive underflows are permitted. underflows
 * and overflows do not occur in the computation of the unscaled
 * sum of squares for the intermediate components.
 * the definitions of small, intermediate and large components
 * depend on two constants, rdwarf and rgiant. the main
 * restrictions on these constants are that rdwarf**2 not
 * underflow and rgiant**2 not overflow. the constants
 * given here are suitable for every known computer.
 * the function statement is
 *   double precision function pnl_minpack_enorm(n,x)
 * where
 *   n is a positive integer input variable.
 *   x is an input array of length n.
 * argonne national laboratory. minpack project. march 1980.
 * burton s. garbow, kenneth e. hillstrom, jorge j. more
 */

double pnl_minpack_enorm(int n, const double *x)
{

  double rdwarf=3.834e-20;
  double rgiant=1.304e19;

  /* System generated locals */
  double ret_val, d__1;

  /* Local variables */
  int i;
  double s1, s2, s3, xabs, x1max, x3max, agiant, floatn;

  /* Function Body */
  s1 = 0.;
  s2 = 0.;
  s3 = 0.;
  x1max = 0.;
  x3max = 0.;
  floatn = (double) (n);
  agiant = rgiant / floatn;
  for ( i = 0 ; i < n; i++ ) 
    {
      xabs = fabs(x[i]);
      if (xabs > rdwarf && xabs < agiant)
        {
          /* sum for intermediate components. */
          d__1 = xabs;
          s2 += d__1 * d__1;
          continue;
        }
      if (xabs <= rdwarf) 
        {
          /* sum for small components. */
          if (xabs <= x3max)
            {
              if (xabs != 0.) 
                {
                  d__1 = xabs / x3max;
                  s3 += d__1 * d__1;
                }
            }
          else
            {
              d__1 = x3max / xabs;
              s3 = 1. + s3 * (d__1 * d__1);
              x3max = xabs;
            }
          continue;
        }

      /* sum for large components. */

      if (xabs > x1max)
        {
          d__1 = x1max / xabs;
          s1 = 1. + s1 * (d__1 * d__1);
          x1max = xabs;
        }
      else
        {
          d__1 = xabs / x1max;
          s1 += d__1 * d__1;
        }
    }

  /* calculation of norm. */
  if (s1 == 0.)
    {
      if (s2 == 0.)
        {
          ret_val = x3max * sqrt(s3);
        }
      else
        {
          if (s2 >= x3max) 
            {
              ret_val = sqrt(s2 * (1. + x3max / s2 * (x3max * s3)));
            }
          else 
            {
              ret_val = sqrt(x3max * (s2 / x3max + x3max * s3));
            }
        }
    }
  else
    {
      ret_val = x1max * sqrt(s1 + s2 / x1max / x1max);
    }
  return ret_val;
} 

