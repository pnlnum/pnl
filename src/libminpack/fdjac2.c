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

/* subroutine fdjac2
 * 
 * this subroutine computes a forward-difference approximation
 * to the m by n jacobian matrix associated with a specified
 * problem of m functions in n variables.
 * 
 * the subroutine statement is
 * 
 *   subroutine pnl_minpack_fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
 * 
 * where
 * 
 *   fcn is the name of the user-supplied subroutine which
 *     calculates the functions. fcn must be declared
 *     in an external statement in the user calling
 *     program, and should be written as follows.
 * 
 *     subroutine fcn(m,n,x,fvec,iflag)
 *     integer m,n,iflag
 *     double precision x(n),fvec(m)
 *     ----------
 *     calculate the functions at x and
 *     return this vector in fvec.
 *     ----------
 *     return
 *     end
 * 
 *     the value of iflag should not be changed by fcn unless
 *     the user wants to terminate execution of fdjac2.
 *     in this case set iflag to a negative integer.
 * 
 *   m is a positive integer input variable set to the number
 *     of functions.
 * 
 *   n is a positive integer input variable set to the number
 *     of variables. n must not exceed m.
 * 
 *   x is an input array of length n.
 * 
 *   fvec is an input array of length m which must contain the
 *     functions evaluated at x.
 * 
 *   fjac is an output m by n array which contains the
 *     approximation to the jacobian matrix evaluated at x.
 * 
 *   ldfjac is a positive integer input variable not less than m
 *     which specifies the leading dimension of the array fjac.
 * 
 *   iflag is an integer variable which can be used to terminate
 *     the execution of fdjac2. see description of fcn.
 * 
 *   epsfcn is an input variable used in determining a suitable
 *     step length for the forward-difference approximation. this
 *     approximation assumes that the relative errors in the
 *     functions are of the order of epsfcn. if epsfcn is less
 *     than the machine precision, it is assumed that the relative
 *     errors in the functions are of the order of the machine
 *     precision.
 * 
 *   wa is a work array of length m.
 * 
 * subprograms called
 * 
 *   user-supplied ...... fcn
 * 
 *   minpack-supplied ... dpmpar
 * 
 *   fortran-supplied ... dabs,dmax1,dsqrt
 * 
 * argonne national laboratory. minpack project. march 1980.
 * burton s. garbow, kenneth e. hillstrom, jorge j. more
 * 
 */

 int pnl_minpack_fdjac2(minpack_func_mn fcn, void *p, int m, int n, double
                        *x, const double *fvec, double *fjac, int ldfjac,
                        double epsfcn, double *wa)
{
  /* System generated locals */
  int fjac_dim1, fjac_offset;

  /* Local variables */
  double h__;
  int i, j;
  double eps, temp, epsmch;
  int iflag;

  /* Parameter adjustments */
  --wa;
  --fvec;
  --x;
  iflag = 0; /* to avoid a warning */
  fjac_dim1 = ldfjac;
  fjac_offset = 1 + fjac_dim1 * 1;
  fjac -= fjac_offset;

  /* Function Body */

  /*     epsmch is the machine precision. */

  epsmch = pnl_minpack_dpmpar(1);

  eps = sqrt((MAX(epsfcn,epsmch)));
  for (j = 1; j <= n; ++j) 
    {
      temp = x[j];
      h__ = eps * fabs(temp);
      if (h__ == 0.)
        {
          h__ = eps;
        }
      x[j] = temp + h__;
      iflag = (*fcn)(p, m, n, &x[1], &wa[1], 1);
      if (iflag < 0) {
        return iflag;
      }
      x[j] = temp;
      for (i = 1; i<= m; ++i)
        {
          fjac[i + j * fjac_dim1] = (wa[i] - fvec[i]) / h__;
        }
    }
  return iflag;
} 

