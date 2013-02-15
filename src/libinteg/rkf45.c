
/*
 * Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as 
 * published by  the Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License  along with this program.  If not, see
 * <http://www.gnu.org/licenses/>. 
 *
 * This code is based on the original Fortran code 
 *
 * fehlberg fourth-fifth order runge-kutta method                     
 * 
 * written by h.a.watts and l.f.shampine                              
 * sandia laboratories                                                
 * albuquerque,new mexico                                             
 * 
 * It has been translated using f2c and many parts have been hand rewritten
 * by Jérôme Lelong to obtain a more easily readable and maintenable code.
 *
 */


#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_integration.h"

static int rkfs (PnlODEFunc *f, double *y, double *t, double *tout,
          double *relerr, double *abserr, int *iflag, double *yp,
          double *h, double *f1, double *f2, double *f3, double *f4,
          double *f5, double *savre, double *savae, int *nfe, int *kop,
          int *init, int *jflag, int *kflag);
static int fehl (PnlODEFunc *f,  double *y, double *t, double *h,
                 double *yp, double *f1, double *f2, double *f3, double
                 *f4, double *f5, double *s);



/*
 * r4_epsilon returns the roundoff unit for the double type.
 * 
 * Discussion:
 * 
 *     The roundoff unit is a number R which is a power of 2 with the
 *     property that, to the precision of the computer's arithmetic,
 *       1 < 1 + R
 *     but
 *       1 = ( 1 + R / 2 )
 * 
 * Licensing: This code is distributed under the GNU LGPL license.
 * 
 * Modified: 01 July 2004
 * 
 * Author: John Burkardt
 * 
 * Parameters:
 * 
 *     Output, the double round-off unit.
 */
static double r4_epsilon ( void )
{
  double value = 1.0;
  while ( 1.0 < 1.0 + value )
  {
    value /= 2.0;
  }
  value *= 2.0;
  return value;
}


static void init_err (double *err)
{
  if ( *err == 0. )
    {
      *err = sqrt (r4_epsilon ());
    }
}

/* 
 *      fehlberg fourth-fifth order runge-kutta method
 * 
 *      written by h.a.watts and l.f.shampine
 *                    sandia laboratories
 *                   albuquerque,new mexico
 * 
 *     rkf45 is primarily designed to solve non-stiff and mildly stiff
 *     differential equations when derivative evaluations are inexpensive.
 *     rkf45 should generally not be used when the user is demanding
 *     high accuracy.
 * 
 *  abstract
 * 
 *     subroutine  rkf45  integrates a system of neqn first order
 *     ordinary differential equations of the form
 *              dy(i)/dt = f(t,y(1),y(2),...,y(neqn))
 *               where the y(i) are given at t .
 *     typically the subroutine is used to integrate from t to tout but it
 *     can be used as a one-step integrator to advance the solution a
 *     single step in the direction of tout.  on return the parameters in
 *     the call list are set for continuing the integration. the user has
 *     only to call rkf45 again (and perhaps define a new value for tout).
 *     actually, rkf45 is an interfacing routine which calls subroutine
 *     rkfs for the solution.  rkfs in turn calls subroutine  fehl which
 *     computes an approximate solution over one step.
 * 
 *     rkf45  uses the runge-kutta-fehlberg (4,5)  method described
 *     in the reference
 *     e.fehlberg , low-order classical runge-kutta formulas with stepsize
 *                  control , nasa tr r-315
 * 
 *     the performance of rkf45 is illustrated in the reference
 *     l.f.shampine,h.a.watts,s.davenport, solving non-stiff ordinary
 *                  differential equations-the state of the art ,
 *                  sandia laboratories report sand75-0182 ,
 *                  to appear in siam review.
 * 
 * 
 *     the parameters represent-
 *       f -- subroutine f(t,y,yp) to evaluate derivatives yp(i)=dy(i)/dt
 *       neqn -- number of equations to be integrated
 *       y(*) -- solution vector at t
 *       t -- independent variable
 *       tout -- output point at which solution is desired
 *       relerr,abserr -- relative and absolute error tolerances for local
 *             error test. at each step the code requires that
 *                  abs(local error) .le. relerr*abs(y) + abserr
 *             for each component of the local error and solution vectors
 *       iflag -- indicator for status of integration
 *       work(*) -- array to hold information internal to rkf45 which is
 *             necessary for subsequent calls. must be dimensioned
 *             at least  3+6*neqn
 *       iwork(*) -- integer array used to hold information internal to
 *             rkf45 which is necessary for subsequent calls. must be
 *             dimensioned at least  5
 * 
 * 
 *   first call to rkf45
 * 
 *     the user must provide storage in his calling program for the arrays
 *     in the call list  -      y(neqn) , work(3+6*neqn) , iwork(5)  ,
 *     declare f in an external statement, supply subroutine f(t,y,yp) and
 *     initialize the following parameters-
 * 
 *       neqn -- number of equations to be integrated.  (neqn .ge. 1)
 *       y(*) -- vector of initial conditions
 *       t -- starting point of integration , must be a variable
 *       tout -- output point at which solution is desired.
 *             t=tout is allowed on the first call only, in which case
 *             rkf45 returns with iflag=2 if continuation is possible.
 *       relerr,abserr -- relative and absolute local error tolerances
 *             which must be non-negative. relerr must be a variable while
 *             abserr may be a constant. the code should normally not be
 *             used with relative error control smaller than about 1.e-8 .
 *             to avoid limiting precision difficulties the code requires
 *             relerr to be larger than an internally computed relative
 *             error parameter which is machine dependent. in particular,
 *             pure absolute error is not permitted. if a smaller than
 *             allowable value of relerr is attempted, rkf45 increases
 *             relerr appropriately and returns control to the user before
 *             continuing the integration.
 *       iflag -- +1,-1  indicator to initialize the code for each new
 *             problem. normal input is +1. the user should set iflag=-1
 *             only when one-step integrator control is essential. in this
 *             case, rkf45 attempts to advance the solution a single step
 *             in the direction of tout each time it is called. since this
 *             mode of operation results in extra computing overhead, it
 *             should be avoided unless needed.
 * 
 * 
 *   output from rkf45
 * 
 *       y(*) -- solution at t
 *       t -- last point reached in integration.
 *       iflag = 2 -- integration reached tout. indicates successful return
 *                    and is the normal mode for continuing integration.
 *             =-2 -- a single successful step in the direction of tout
 *                    has been taken. normal mode for continuing
 *                    integration one step at a time.
 *             = 3 -- integration was not completed because relative error
 *                    tolerance was too small. relerr has been increased
 *                    appropriately for continuing.
 *             = 4 -- integration was not completed because more than
 *                    3000 derivative evaluations were needed. this
 *                    is approximately 500 steps.
 *             = 5 -- integration was not completed because solution
 *                    vanished making a pure relative error test
 *                    impossible. must use non-zero abserr to continue.
 *                    using the one-step integration mode for one step
 *                    is a good way to proceed.
 *             = 6 -- integration was not completed because requested
 *                    accuracy could not be achieved using smallest
 *                    allowable stepsize. user must increase the error
 *                    tolerance before continued integration can be
 *                    attempted.
 *             = 7 -- it is likely that rkf45 is inefficient for solving
 *                    this problem. too much output is restricting the
 *                    natural stepsize choice. use the one-step integrator
 *                    mode.
 *             = 8 -- invalid input parameters
 *                    this indicator occurs if any of the following is
 *                    satisfied -   neqn .le. 0
 *                                  t=tout  and  iflag .ne. +1 or -1
 *                                  relerr or abserr .lt. 0.
 *                                  iflag .eq. 0  or  .lt. -2  or  .gt. 8
 *       work(*),iwork(*) -- information which is usually of no interest
 *                    to the user but necessary for subsequent calls.
 *                    work(1),...,work(neqn) contain the first derivatives
 *                    of the solution vector y at t. work(neqn+1) contains
 *                    the stepsize h to be attempted on the next step.
 *                    iwork(1) contains the derivative evaluation counter.
 * 
 * 
 *   subsequent calls to rkf45
 * 
 *     subroutine rkf45 returns with all information needed to continue
 *     the integration. if the integration reached tout, the user need onl
 *     define a new tout and call rkf45 again. in the one-step integrator
 *     mode (iflag=-2) the user must keep in mind that each step taken is
 *     in the direction of the current tout. upon reaching tout (indicated
 *     by changing iflag to 2),the user must then define a new tout and
 *     reset iflag to -2 to continue in the one-step integrator mode.
 * 
 *     if the integration was not completed but the user still wants to
 *     continue (iflag=3,4 cases), he just calls rkf45 again. with iflag=3
 *     the relerr parameter has been adjusted appropriately for continuing
 *     the integration. in the case of iflag=4 the function counter will
 *     be reset to 0 and another 3000 function evaluations are allowed.
 * 
 *     however,in the case iflag=5, the user must first alter the error
 *     criterion to use a positive value of abserr before integration can
 *     proceed. if he does not,execution is terminated.
 * 
 *     also,in the case iflag=6, it is necessary for the user to reset
 *     iflag to 2 (or -2 when the one-step integration mode is being used)
 *     as well as increasing either abserr,relerr or both before the
 *     integration can be continued. if this is not done, execution will
 *     be terminated. the occurrence of iflag=6 indicates a trouble spot
 *     (solution is changing rapidly,singularity may be present) and it
 *     often is inadvisable to continue.
 * 
 *     if iflag=7 is encountered, the user should use the one-step
 *     integration mode with the stepsize determined by the code or
 *     consider switching to the adams codes de/step,intrp. if the user
 *     insists upon continuing the integration with rkf45, he must reset
 *     iflag to 2 before calling rkf45 again. otherwise,execution will be
 *     terminated.
 * 
 *     if iflag=8 is obtained, integration can not be continued unless
 *     the invalid input parameters are corrected.
 * 
 *     it should be noted that the arrays work,iwork contain information
 *     required for subsequent integration. accordingly, work and iwork
 *     should not be altered.
 */

int pnl_ode_rkf45 (PnlODEFunc *f, double *y, double t, double tout, 
                   double relerr, double abserr, int *flag)
{
  int k1, k2, k3, k4, k5, k6, k1m, neqn;
  double *work; int *iwork;
  neqn = f->neqn;
  *flag = 1;

  /* init errors */
  init_err (&relerr);
  init_err (&abserr);

  if ( (work = malloc ((3 + 6 * neqn) * sizeof(double))) == NULL ) return FAIL; 
  if ( (iwork = malloc (5 * sizeof(int))) == NULL ) return FAIL; 

  /* Compute indices for the splitting of the work array */
  k1m = neqn;
  k1 = k1m + 1;
  k2 = k1 + neqn;
  k3 = k2 + neqn;
  k4 = k3 + neqn;
  k5 = k4 + neqn;
  k6 = k5 + neqn;

  rkfs (f, y, &t, &tout, &relerr, &abserr, flag, &work[0],
        &work[k1m], &work[k1], &work[k2], &work[k3], &work[k4],
        &work[k5], &work[k6], &work[k6 + 1], &iwork[0], &iwork[1],
        &iwork[2], &iwork[3], &iwork[4]);

  free (work); work = NULL;
  free (iwork); iwork = NULL;
  if ( *flag == 2 ) return OK; else return FAIL;
}

int pnl_ode_rkf45_step (PnlODEFunc *f, double *y, double *t,
                        double tout, double *relerr, double abserr, double
                        *work, int *iwork, int *iflag)
{
  int k1, k2, k3, k4, k5, k6, k1m, neqn;
  neqn = f->neqn;
  *iflag = - 2;
  /* Compute indices for the splitting of the work array */
  k1m = neqn;
  k1 = k1m + 1;
  k2 = k1 + neqn;
  k3 = k2 + neqn;
  k4 = k3 + neqn;
  k5 = k4 + neqn;
  k6 = k5 + neqn;

  /* init errors */
  init_err (relerr);
  init_err (&abserr);

  rkfs (f, y, t, &tout, relerr, &abserr, iflag, &work[0],
        &work[k1m], &work[k1], &work[k2], &work[k3], &work[k4],
        &work[k5], &work[k6], &work[k6 + 1], &iwork[0], &iwork[1],
        &iwork[2], &iwork[3], &iwork[4]);
  if ( *iflag == -2 ) return OK; else return FALSE;
}


/*
 * rkfs integrates a system of first order ordinary differential    
 * equations as described in the comments for rkf45 .               
 * the arrays yp,f1,f2,f3,f4,and f5 (of dimension at least neqn) and
 * the variables h,savre,savae,nfe,kop,init,jflag,and kflag are used
 * internally by the code and appear in the call list to eliminate  
 * local retention of variables between calls. accordingly, they    
 * should not be altered. items of possible interest are            
 * yp - derivative of solution vector at t                          
 * h  - an appropriate stepsize to be used for the next ste p
 * nfe- counter on the number of derivative function evaluations    
 */
static int rkfs (PnlODEFunc *f,  double *y, double *t, double *tout,
          double *relerr, double *abserr, int *iflag, double *yp,
          double *h, double *f1, double *f2, double *f3, double *f4,
          double *f5, double *savre, double *savae, int *nfe, int *kop,
          int *init, int *jflag, int *kflag)
{
  /* Initialized data */

  int neqn = f->neqn;
  static double remin = 1e-12;
  static int maxnfe = 3000;

  /* System generated locals */
  double d__1, d__2, d__3, d__4;

  /* Builtin functions */
  /* Subroutine */ int s_stop (char *, unsigned long);

  /* Local variables */
  static double a;
  static int k;
  static double s, ae, ee, dt, et, u26, rer, tol, ypk;
  static double hmin, toln;
  static int mflag;
  static double scale, eeoet;

  static int hfaild;
  static double esttol, twoeps;
  static int output;



  /*  remin is the minimum acceptable value of relerr.  attempts */
  /*  to obtain higher accuracy with this subroutine are usually */
  /*  very expensive and often unsuccessful. */


  /*
   * the expense is controlled by restricting the number
   * of function evaluations to be approximately maxnfe.
   * as set, this corresponds to about 500 steps.
   *
   * here two constants emboding the machine epsilon is present
   * twoesp is set to twice the machine epsilon while u26 is set
   * to 26 times the machine epsilon
   */
  twoeps = 2. * pnl_dlamch ("p");
  u26 = twoeps * 13.;


  /*     check input parameters */
  if (neqn < 1)
    {
      *iflag = 8;
      return FAIL;
    }
  if (*relerr < 0. || *abserr < 0.)
    {
      *iflag = 8;
      return FAIL;
    }
  mflag = ABS (*iflag);

  /*     is this the first call */
  if (mflag == 1)
    {
      goto L50;
    }

  /*     check continuation possibilities */

  if (*t == *tout && *kflag != 3)
    {
      *iflag = 8;
      return FAIL;
    }
  if (mflag != 2)
    {
      goto L25;
    }

  /*     iflag = +2 or -2 */
  if (*kflag == 3)
    {
      goto L45;
    }
  if (*init == 0)
    {
      goto L45;
    }
  if (*kflag == 4)
    {
      goto L40;
    }
  if (*kflag == 5 && *abserr == 0.)
    {
      goto L30;
    }
  if (*kflag == 6 && *relerr <= *savre && *abserr <= *savae)
    {
      goto L30;
    }
  goto L50;

  /*     iflag = 3,4,5,6,7 or 8 */
L25:
  if (*iflag == 3)
    {
      goto L45;
    }
  if (*iflag == 4)
    {
      goto L40;
    }
  if (*iflag == 5 && *abserr > 0.)
    {
      goto L45;
    }

  /*     integration cannot be continued since user did not respond to */
  /*     the instructions pertaining to iflag=5,6,7 or 8 */
L30:
  return 0;

  /*     reset function evaluation counter */
L40:
  *nfe = 0;
  if (mflag == 2)
    {
      goto L50;
    }

  /*     reset flag value from previous call */
L45:
  *iflag = *jflag;
  if (*kflag == 3)
    {
      mflag = ABS (*iflag);
    }

  /*     save input iflag and set continuation flag value for subsequent */
  /*     input checking */
L50:
  *jflag = *iflag;
  *kflag = 0;

  /*     save relerr and abserr for checking input on subsequent calls */
  *savre = *relerr;
  *savae = *abserr;

  /*     restrict relative error tolerance to be at least as large as */
  /*     2*eps+remin to avoid limiting precision difficulties arising */
  /*     from impossible accuracy requests */

  rer = twoeps + remin;
  if (*relerr < rer)
    {
      /*     relative error tolerance too small */
      *relerr = rer;
      *iflag = 3;
      *kflag = 3;
      return FAIL;
    }

  dt = *tout - *t;

  if (mflag == 1)
    {
      *init = 0;
      *kop = 0;

      a = *t;
      PNL_EVAL_ODEFUNC (f, a, y, yp);
      *nfe = 1;
      if (*t == *tout)
        {
          *iflag = 2;
          return OK;
        }
    }
  /* initialization -- */
  /* set initialization completion indicator,init */
  /* set indicator for too many output points,kop */
  /* evaluate initial derivatives */
  /* set counter for function evaluations,nfe */
  /* evaluate initial derivatives */
  /* set counter for function evaluations,nfe */
  /* estimate starting stepsize */
  if (*init == 0)
    {
      *init = 1;
      *h = ABS (dt);
      toln = 0.;
      for (k = 0; k < neqn; ++k)
        {
          tol = *relerr * (d__1 = y[k], ABS (d__1)) + *abserr;
          if (tol > 0.)
            {
              toln = tol;
              ypk = ABS (yp[k]);
              if (ypk * pnl_pow_i (*h, 5)  > tol)
                {
                  *h = pow (tol / ypk, 0.2);
                }
            }
        }

      if (toln <= 0.)
        {
          *h = 0.;
        }
      d__3 = ABS (*t), d__4 = ABS (dt);
      d__1 = *h, d__2 = u26 * MAX (d__3, d__4);
      *h = MAX (d__1, d__2);
      *jflag = 2 * PNL_SIGN (*iflag);

    }


  /*     set stepsize for integration in the direction from t to tout */
  *h = ABS (*h) * PNL_SIGN (dt);

  /*     test to see if rkf45 is being severely impacted by too many */
  /*     output points */

  if (ABS (*h) >= ABS (dt) * 2.)
    {
      ++(*kop);
    }
  /*     unnecessary frequency of output */
  if (*kop == 100)
    {
      *kop = 0;
      *iflag = 7;
      return FAIL;
    }

  /*     if too close to output point,extrapolate and return */
  if (ABS (dt) <= u26 * ABS (*t))
    {
      for (k = 0; k < neqn; ++k)
        {
          /* L90: */
          y[k] += dt * yp[k];
        }
      a = *tout;
      PNL_EVAL_ODEFUNC (f, a, y, yp);
      ++(*nfe);
      *t = *tout;
      *iflag = 2;
      return OK;
    }
  /*     initialize output point indicator */

  output = FALSE;

  /*     to avoid premature underflow in the error tolerance function, */
  /*     scale the error tolerances */

  scale = 2. / *relerr;
  ae = scale * *abserr;

  /*     step by step integration */

L100:
  hfaild = FALSE;

  /*     set smallest allowable stepsize */

  hmin = u26 * ABS (*t);

  /*     adjust stepsize if necessary to hit the output point. */
  /*     look ahead two steps to avoid drastic changes in the stepsize and */
  /*     thus lessen the impact of output points on the code. */

  dt = *tout - *t;
  if (ABS (dt) >= ABS (*h) * 2.)
    {
      goto L200;
    }
  if (ABS (dt) > ABS (*h))
    {
      goto L150;
    }

  /*     the next successful step will complete the integration to the */
  /*     output point */

  output = TRUE;
  *h = dt;
  goto L200;

L150:
  *h = dt * .5;



  /*     core integrator for taking a single step */

  /*     the tolerances have been scaled to avoid premature underflow in */
  /*     computing the error tolerance function et. */
  /*     to avoid problems with zero crossings,relative error is measured */
  /*     using the average of the magnitudes of the solution at the */
  /*     beginning and end of a step. */
  /*     the error estimate formula has been grouped to control loss of */
  /*     significance. */
  /*     to distinguish the various arguments, h is not permitted */
  /*     to become smaller than 26 units of roundoff in t. */
  /*     practical limits on the change in the stepsize are enforced to */
  /*     smooth the stepsize selection process and to avoid excessive */
  /*     chattering on problems having discontinuities. */
  /*     to prevent unnecessary failures, the code uses 9/10 the stepsize */
  /*     it estimates will succeed. */
  /*     after a step failure, the stepsize is not allowed to increase for */
  /*     the next attempted step. this makes the code more efficient on */
  /*     problems having discontinuities and more effective in general */
  /*     since local extrapolation is being used and extra caution seems */
  /*     warranted. */


  /*     test number of derivative function evaluations. */
  /*     if okay,try to advance the integration from t to t+h */

L200:
  if (*nfe > maxnfe)
    {
      /*     too much work */
      *iflag = 4;
      *kflag = 4;
      return FAIL;
    }
  /*     advance an approximate solution over one step of length h */

  fehl (f, y, t, h, yp, f1, f2, f3, f4, f5, f1);
  *nfe += 5;

  /*     compute and test allowable tolerances versus local error estimates */
  /*     and remove scaling of tolerances. note that relative error is */
  /*     measured with respect to the average of the magnitudes of the */
  /*     solution at the beginning and end of the step. */

  eeoet = 0.;
  for (k = 0; k < neqn; ++k)
    {
      et = (d__1 = y[k], ABS (d__1)) + (d__2 = f1[k], ABS (d__2)) + ae;
      if (et <= 0.)
        {
          *iflag = 5;
          return FAIL;
        }
      else
        {
          ee = (d__1 = yp[k] * -2090. + (f3[k] * 21970. - f4[k] * 15048.) +
                (f2[k] * 22528. - f5[k] * 27360.), ABS (d__1));
          d__1 = eeoet, d__2 = ee / et;
          eeoet = MAX (d__1, d__2);
        }
    }

  esttol = ABS (*h) * eeoet * scale / 752400.;

  if (esttol <= 1.)
    {
      goto L260;
    }


  /*  unsuccessful step */
  /*  reduce the stepsize , try again */
  /*  the decrease is limited to a factor of 1/10 */

  hfaild = TRUE;
  output = FALSE;
  if (esttol < 59049.)
    {
      s = .9 / pow (esttol, 0.2);
    }
  else
    {
      s = 0.1;
    }
  *h = s * *h;
  if (ABS (*h) > hmin)
    {
      goto L200;
    }

  /*     requested error unattainable at smallest allowable stepsize */
  *iflag = 6;
  *kflag = 6;
  return FAIL;


  /* successful step */
  /* store solution at t+h */
  /* and evaluate derivatives there */

L260:
  *t += *h;
  for (k = 0; k < neqn; ++k)
    {
      /* L270: */
      y[k] = f1[k];
    }
  a = *t;
  PNL_EVAL_ODEFUNC (f, a, y, yp);
  ++(*nfe);


  /* choose next stepsize */
  /* the increase is limited to a factor of 5 */
  /* if step failure has just occurred, next */
  /* stepsize is not allowed to increase */

  s = 5.;
  if (esttol > 1.889568e-4)
    {
      s = .9 / pow (esttol, 0.2);
    }
  if (hfaild)
    {
      s = MIN (s, 1.);
    }
  /* Computing MAX */
  d__2 = s * ABS (*h);
  d__1 = MAX (d__2, hmin);
  *h = ABS (d__1) * PNL_SIGN (*h);

  /* end of core integrator */


  /* should we take another step */

  if (output)
    {
      *t = *tout;
      *iflag = 2;
      return OK;
    }
  if (*iflag > 0)
    {
      goto L100;
    }


  /* integration successfully completed */

  /* one-step mode */
  *iflag = -2;
  return OK;

}	

/*
 *    fehl integrates a system of neqn first order
 *    ordinary differential equations of the form
 *             dy(i)/dt=f(t,y(1),---,y(neqn))
 *    where the initial values y(i) and the initial derivatives
 *    yp(i) are specified at the starting point t. fehl advances
 *    the solution over the fixed step h and returns
 *    the fifth order (sixth order accurate locally) solution
 *    approximation at t+h in array s(i).
 *    f1,---,f5 are arrays of dimension neqn which are needed
 *    for internal storage.
 *    the formulas have been grouped to control loss of significance.
 *    fehl should be called with an h not smaller than 13 units of
 *    roundoff in t so that the various independent arguments can be
 *    distinguished.
 */
static int fehl (PnlODEFunc *f, double *y, double *t, double *h,
                 double *yp, double *f1, double *f2, double *f3, double
                 *f4, double *f5, double *s)
{
  int neqn, k;
  double x, ch;

  neqn = f->neqn;

  /* Function Body */
  ch = *h / 4.;
  for (k = 0; k < neqn; k++)
    {
      f5[k] = y[k] + ch * yp[k];
    }
  x = *t + ch;
  PNL_EVAL_ODEFUNC (f, x, f5, f1);

  ch = *h * 3. / 32.;
  for (k = 0; k < neqn; k++)
    {
      f5[k] = y[k] + ch * (yp[k] + f1[k] * 3.);
    }
  x = *t + *h * 3. / 8.;
  PNL_EVAL_ODEFUNC (f, x, f5, f2);

  ch = *h / 2197.;
  for (k = 0; k < neqn; k++)
    {
      f5[k] = y[k] + ch * (yp[k] * 1932. + (f2[k] * 7296. - f1[k] * 7200.));
    }
  x = *t + *h * 12. / 13.;
  PNL_EVAL_ODEFUNC (f, x, f5, f3);

  ch = *h / 4104.;
  for (k = 0; k < neqn; k++)
    {
      f5[k] = y[k] + ch * (yp[k] * 8341. - f3[k] * 845. + (f2[k] * 29440. -
                                                           f1[k] * 32832.));
    }
  x = *t + *h;
  PNL_EVAL_ODEFUNC (f, x, f5, f4);

  ch = *h / 20520.;
  for (k = 0; k < neqn; k++)
    {
      f1[k] = y[k] + ch * (yp[k] * -6080. + (f3[k] * 9295. - f4[k] * 5643.)
                           + (f1[k] * 41040. - f2[k] * 28352.));
    }
  x = *t + *h / 2.;
  PNL_EVAL_ODEFUNC (f, x, f1, f5);

  /*     compute approximate solution at t+h */

  ch = *h / 7618050.;
  for (k = 0; k < neqn; k++)
    {
      s[k] = y[k] + ch * (yp[k] * 902880. + (f3[k] * 3855735. - f4[k] *
                                             1371249.) + (f2[k] * 3953664. +
                                                          f5[k] * 277020.));
    }

  return OK;
}
