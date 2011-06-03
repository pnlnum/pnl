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

/*
 * Common source of lmder.c and lmdif.c
 *
 * This file must not be compiled. It is included in lmdif.c and lmder.c
 * Note that either the macro _LMDER_ or _LMDIF_ must be defined before
 * inclusion
 */


{
  /* Initialized data */

  double p1 = .1;
  double p5 = .5;
  double p25 = .25;
  double p75 = .75;
  double p0001 = 1e-4;

  /* System generated locals */
  int fjac_dim1, fjac_offset, i__1, i__2;
  double d__1, d__2, d__3;

  /* Local variables */
  int i__, j, l;
  double par, sum;
  int iter;
  double temp=0., temp1, temp2;
  int iflag;
  double delta=0.0;
  double ratio;
  double fnorm, gnorm;
  double pnorm, xnorm=0.0, fnorm1, actred, dirder, epsmch, prered;
  int info;



  /* Parameter adjustments */
  --wa4;
  --fvec;
  --wa3;
  --wa2;
  --wa1;
  --qtf;
  --ipvt;
  --diag;
  --x;
  fjac_dim1 = ldfjac;
  fjac_offset = 1 + fjac_dim1 * 1;
  fjac -= fjac_offset;

  /* Function Body */

  /*     epsmch is the machine precision. */

  epsmch = pnl_minpack_dpmpar(1);

  info = 0;
  iflag = 0;
  *nfev = 0;
#ifdef _LMDER_
  *njev = 0;
#endif

  /*     check the input parameters for errors. */

  if (n <= 0 || m < n || ldfjac < m || ftol < 0. || xtol < 0. || 
      gtol < 0. || maxfev <= 0 || factor <= 0.) {
    goto L300;
  }
  if (mode != 2) {
    goto L20;
  }
  i__1 = n;
  for (j = 1; j <= i__1; ++j) {
    if (diag[j] <= 0.) {
      goto L300;
    }
    /* L10: */
  }
L20:

  /*     evaluate the function at the starting point */
  /*     and calculate its norm. */

#ifdef _LMDER_ 
  iflag = (*fcn)(p, m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, 1);
#elif defined(_LMDIF_)
  iflag = (*fcn)(p, m, n, &x[1], &fvec[1], 1);
#endif
  *nfev = 1;
  if (iflag < 0) {
    goto L300;
  }
  fnorm = pnl_minpack_enorm(m, &fvec[1]);

  /*     initialize levenberg-marquardt parameter and iteration counter. */

  par = 0.;
  iter = 1;

  /*     beginning of the outer loop. */

L30:

  /*        calculate the jacobian matrix. */

#ifdef _LMDER_
  iflag = (*fcn)(p, m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, 2);
  ++(*njev);
#elif defined(_LMDIF_)
  iflag = pnl_minpack_fdjac2(fcn, p, m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac,
                             epsfcn, &wa4[1]);
  *nfev += n;
#endif

  if (iflag < 0) {
    goto L300;
  }

  /*        if requested, call fcn to enable printing of iterates. */

  if (nprint <= 0) {
    goto L40;
  }
  iflag = 0;
  if ((iter - 1) % nprint == 0) {

#ifdef _LMDER_
    iflag = (*fcn)(p, m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, 0);
#elif defined(_LMDIF_)
    iflag = (*fcn)(p, m, n, &x[1], &fvec[1], 0);
#endif
  }
  if (iflag < 0) {
    goto L300;
  }
L40:

  /*        compute the qr factorization of the jacobian. */

  pnl_minpack_qrfac(m, n, &fjac[fjac_offset], ldfjac, TRUE, &ipvt[1], n, &wa1[1], &
                    wa2[1], &wa3[1]);

  /*        on the first iteration and if mode is 1, scale according */
  /*        to the norms of the columns of the initial jacobian. */

  if (iter != 1) {
    goto L80;
  }
  if (mode == 2) {
    goto L60;
  }
  i__1 = n;
  for (j = 1; j <= i__1; ++j) {
    diag[j] = wa2[j];
    if (wa2[j] == 0.) {
      diag[j] = 1.;
    }
    /* L50: */
  }
L60:

  /*        on the first iteration, calculate the norm of the scaled x */
  /*        and initialize the step bound delta. */

  i__1 = n;
  for (j = 1; j <= i__1; ++j) {
    wa3[j] = diag[j] * x[j];
    /* L70: */
  }
  xnorm = pnl_minpack_enorm(n, &wa3[1]);
  delta = factor * xnorm;
  if (delta == 0.) {
    delta = factor;
  }
L80:

  /*        form (q transpose)*fvec and store the first n components in */
  /*        qtf. */

  i__1 = m;
  for (i__ = 1; i__ <= i__1; ++i__) {
    wa4[i__] = fvec[i__];
    /* L90: */
  }
  i__1 = n;
  for (j = 1; j <= i__1; ++j) {
    if (fjac[j + j * fjac_dim1] == 0.) {
      goto L120;
    }
    sum = 0.;
    i__2 = m;
    for (i__ = j; i__ <= i__2; ++i__) {
      sum += fjac[i__ + j * fjac_dim1] * wa4[i__];
      /* L100: */
    }
    temp = -sum / fjac[j + j * fjac_dim1];
    i__2 = m;
    for (i__ = j; i__ <= i__2; ++i__) {
      wa4[i__] += fjac[i__ + j * fjac_dim1] * temp;
      /* L110: */
    }
L120:
    fjac[j + j * fjac_dim1] = wa1[j];
    qtf[j] = wa4[j];
    /* L130: */
  }

  /*        compute the norm of the scaled gradient. */

  gnorm = 0.;
  if (fnorm == 0.) {
    goto L170;
  }
  i__1 = n;
  for (j = 1; j <= i__1; ++j) {
    l = ipvt[j];
    if (wa2[l] == 0.) {
      goto L150;
    }
    sum = 0.;
    i__2 = j;
    for (i__ = 1; i__ <= i__2; ++i__) {
      sum += fjac[i__ + j * fjac_dim1] * (qtf[i__] / fnorm);
      /* L140: */
    }
    /* Computing MAX */
    d__2 = gnorm, d__3 = fabs(sum / wa2[l]);
    gnorm = MAX(d__2,d__3);
L150:
    /* L160: */
    ;
  }
L170:

  /*        test for convergence of the gradient norm. */

  if (gnorm <= gtol) {
    info = 4;
  }
  if (info != 0) {
    goto L300;
  }

  /*        rescale if necessary. */

  if (mode == 2) {
    goto L190;
  }
  i__1 = n;
  for (j = 1; j <= i__1; ++j) {
    /* Computing MAX */
    d__1 = diag[j], d__2 = wa2[j];
    diag[j] = MAX(d__1,d__2);
    /* L180: */
  }
L190:

  /*        beginning of the inner loop. */

L200:

  /*           determine the levenberg-marquardt parameter. */

  pnl_minpack_lmpar(n, &fjac[fjac_offset], ldfjac, &ipvt[1], &diag[1], &qtf[1], delta,
                    &par, &wa1[1], &wa2[1], &wa3[1], &wa4[1]);

  /*           store the direction p and x + p. calculate the norm of p. */

  i__1 = n;
  for (j = 1; j <= i__1; ++j) {
    wa1[j] = -wa1[j];
    wa2[j] = x[j] + wa1[j];
    wa3[j] = diag[j] * wa1[j];
    /* L210: */
  }
  pnorm = pnl_minpack_enorm(n, &wa3[1]);

  /*           on the first iteration, adjust the initial step bound. */

  if (iter == 1) {
    delta = MIN(delta,pnorm);
  }

  /*           evaluate the function at x + p and calculate its norm. */

#ifdef _LMDER_
  iflag = (*fcn)(p, m, n, &wa2[1], &wa4[1], &fjac[fjac_offset], ldfjac, 1);
#elif defined(_LMDIF_)
  iflag = (*fcn)(p, m, n, &wa2[1], &wa4[1], 1);
#endif
  ++(*nfev);
  if (iflag < 0) {
    goto L300;
  }
  fnorm1 = pnl_minpack_enorm(m, &wa4[1]);

  /*           compute the scaled actual reduction. */

  actred = -1.;
  if (p1 * fnorm1 < fnorm) {
    /* Computing 2nd power */
    d__1 = fnorm1 / fnorm;
    actred = 1. - d__1 * d__1;
  }

  /*           compute the scaled predicted reduction and */
  /*           the scaled directional derivative. */

  i__1 = n;
  for (j = 1; j <= i__1; ++j) {
    wa3[j] = 0.;
    l = ipvt[j];
    temp = wa1[l];
    i__2 = j;
    for (i__ = 1; i__ <= i__2; ++i__) {
      wa3[i__] += fjac[i__ + j * fjac_dim1] * temp;
      /* L220: */
    }
    /* L230: */
  }
  temp1 = pnl_minpack_enorm(n, &wa3[1]) / fnorm;
  temp2 = sqrt(par) * pnorm / fnorm;
  /* Computing 2nd power */
  d__1 = temp1;
  /* Computing 2nd power */
  d__2 = temp2;
  prered = d__1 * d__1 + d__2 * d__2 / p5;
  /* Computing 2nd power */
  d__1 = temp1;
  /* Computing 2nd power */
  d__2 = temp2;
  dirder = -(d__1 * d__1 + d__2 * d__2);

  /*           compute the ratio of the actual to the predicted */
  /*           reduction. */

  ratio = 0.;
  if (prered != 0.) {
    ratio = actred / prered;
  }

  /*           update the step bound. */

  if (ratio > p25) {
    goto L240;
  }
  if (actred >= 0.) {
    temp = p5;
  }
  if (actred < 0.) {
    temp = p5 * dirder / (dirder + p5 * actred);
  }
  if (p1 * fnorm1 >= fnorm || temp < p1) {
    temp = p1;
  }
  /* Computing MIN */
  d__1 = delta, d__2 = pnorm / p1;
  delta = temp * MIN(d__1,d__2);
  par /= temp;
  goto L260;
L240:
  if (par != 0. && ratio < p75) {
    goto L250;
  }
  delta = pnorm / p5;
  par = p5 * par;
L250:
L260:

  /*           test for successful iteration. */

  if (ratio < p0001) {
    goto L290;
  }

  /*           successful iteration. update x, fvec, and their norms. */

  i__1 = n;
  for (j = 1; j <= i__1; ++j) {
    x[j] = wa2[j];
    wa2[j] = diag[j] * x[j];
    /* L270: */
  }
  i__1 = m;
  for (i__ = 1; i__ <= i__1; ++i__) {
    fvec[i__] = wa4[i__];
    /* L280: */
  }
  xnorm = pnl_minpack_enorm(n, &wa2[1]);
  fnorm = fnorm1;
  ++iter;
L290:

  /*           tests for convergence. */

  if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1.) {
    info = 1;
  }
  if (delta <= xtol * xnorm) {
    info = 2;
  }
  if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1. && info 
      == 2) {
    info = 3;
  }
  if (info != 0) {
    goto L300;
  }

  /*           tests for termination and stringent tolerances. */

  if (*nfev >= maxfev) {
    info = 5;
  }
  if (fabs(actred) <= epsmch && prered <= epsmch && p5 * ratio <= 1.) {
    info = 6;
  }
  if (delta <= epsmch * xnorm) {
    info = 7;
  }
  if (gnorm <= epsmch) {
    info = 8;
  }
  if (info != 0) {
    goto L300;
  }

  /*           end of the inner loop. repeat if iteration unsuccessful. */

  if (ratio < p0001) {
    goto L200;
  }

  /*        end of the outer loop. */

  goto L30;
L300:

  /*     termination, either normal or user imposed. */

  if (iflag < 0) {
    info = iflag;
  }
  iflag = 0;
  if (nprint > 0) {
#ifdef _LMDER_
    iflag = (*fcn)(p, m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, 0);
#elif defined(_LMDIF_)
    iflag = (*fcn)(p, m, n, &x[1], &fvec[1], 0);
#endif
  }
  return info;
}

