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
 */


#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_internals.h"
#include "cminpack.h"

/**
 * This is wrapper for using PnlFuncRnFuncRn pointer in pnl_root_fsolve instead
 * of minpack_func_nn
 *
 * @param pnl_func this will actually be a pointer to a PnlRnFuncRnDFunc object
 * @param x the point at which the function is evaluated
 * @param n size of the vector x
 * @param fvec an array of size n, containing on output f(x)
 * @param fjac an array of size nxn, containing on output Jac(f)(x) using a
 * column major order
 * @param ldfjac is the number of rows of the Jacobian matrix
 * @param iflag internal minpack variable. if iflag=0, the function does
 * nothing. If iflag=1, it computes fvec. If iflag=2, it computes fjac
 * @return an integer. If the returned value is negative the function hybrd
 * terminates
 */
static int hybrj_fcn (void *pnl_func, int n, const double *x, double *fvec,
                      double *fjac, int ldfjac, int iflag)
{
  PnlVect X, Fvec;
  PnlMat Fjac;
  PnlRnFuncRnDFunc *f;
  f = (PnlRnFuncRnDFunc *) pnl_func;
  X = pnl_vect_wrap_array (x, n);
  if ( iflag == 1) 
    {
      Fvec = pnl_vect_wrap_array (fvec, n);
      PNL_EVAL_RNFUNCRN(f, &X, &Fvec);
    }
  else if ( iflag == 2 )
    {
      Fjac = pnl_mat_wrap_array (fjac, n, n);
      PNL_EVAL_RNFUNCRN_DF(f, &X, &Fjac);
      /*
       * Because Minpack uses the Fortran column-wise storage, 
       * we need to transpose the Jacobian matrix to convert from 
       * row-wise to column-wise storage
       */
      pnl_mat_sq_transpose (&Fjac);
    }
  return 0;
}


/**
 * This is wrapper for using PnlFuncRnFuncRn pointer in pnl_root_fsolve instead
 * of minpack_func_nn
 *
 * @param pnl_func this will actually be a pointer to a PnlRnFuncRn object
 * @param x the point at which the function is evaluated
 * @param n size of the vector x
 * @param fvec a vector of size n, containing on output f(x)
 * @param iflag internal minpack variable. if iflag=0, the function does
 * nothing
 * @return an integer. If the returned value is negative the function hybrd
 * terminates
 */
static int hybrd_fcn (void *pnl_func, int n, const double *x, double *fvec, int iflag)
{
  PnlVect X, Fvec;
  PnlRnFuncRn *f;
  f = (PnlRnFuncRn*) pnl_func;
  if ( iflag == 0 ) return 0;
  X = pnl_vect_wrap_array (x, n);
  Fvec = pnl_vect_wrap_array (fvec, n);
  PNL_EVAL_RNFUNCRN(f, &X, &Fvec);
  return 0;
}

/**
 * This is wrapper for using a PnRnFuncRm pointer in pnl_root_fsolve_lsq instead
 * of minpack_func_mn
 *
 * @param pnl_func this will actually be a pointer to a PnlRnFuncRm object
 * @param x the point at which the function is evaluated
 * @param m number of componenets of f (may be different from n)
 * @param n size of the vector x
 * @param fvec a vector of size m, containing on output f(x)
 * @param iflag internal minpack variable. if iflag=0, the function does
 * nothing
 * @return an integer. If the returned value is negative the function hybrd
 * terminates
 */
static int lmdif_fcn (void *pnl_func, int m, int n, const double *x, double *fvec, int iflag)
{
  PnlVect X, Fvec;
  PnlRnFuncRm *f;
  f = (PnlRnFuncRm*) pnl_func;
  if ( iflag == 0 ) return 0;
  X = pnl_vect_wrap_array (x, n);
  Fvec = pnl_vect_wrap_array (fvec, m);
  PNL_EVAL_RNFUNCRM(f, &X, &Fvec);
  return 0;
}

/**
 * This is wrapper for using a PnRnFuncRm pointer in pnl_root_fsolve_lsq instead
 * of minpack_func_mn
 *
 * @param pnl_func this will actually be a pointer to a PnlRnFuncRm object
 * @param x the point at which the function is evaluated
 * @param m number of componenets of f (may be different from n)
 * @param n size of the vector x
 * @param fvec a vector of size m, containing on output f(x)
 * @param fjac contains the Jacobian matrix on exit stored in a column major
 * order
 * @param ldfjac is the number of rows of the Jacobian matrix
 * @param iflag internal minpack variable. if iflag=0, the function does
 * nothing
 * @return an integer. If the returned value is negative the function hybrd
 * terminates
 */
static int lmder_fcn (void *pnl_func, int m, int n, const double *x, double *fvec, 
                      double *fjac, int ldfjac, int iflag)
{
  PnlVect X, Fvec;
  PnlRnFuncRmDFunc *f;
  f = (PnlRnFuncRmDFunc *) pnl_func;
  X = pnl_vect_wrap_array (x, n);
  if ( iflag == 1) 
    {
      Fvec = pnl_vect_wrap_array (fvec, m);
      PNL_EVAL_RNFUNCRM(f, &X, &Fvec);
    }
  else if ( iflag == 2 )
    {
      int i, j;
      double *jac;
      PnlMat Jac;
      jac = MALLOC_DOUBLE(m*n);
      Jac = pnl_mat_wrap_array (jac, m, n);
      PNL_EVAL_RNFUNCRN_DF(f, &X, &Jac);
      /*
       * Because Minpack uses the Fortran column-wise storage, 
       * we need to transpose the Jacobian matrix to convert from 
       * row-wise to column-wise storage
       */
      for ( i=0 ; i<m ; i++ )
        {
          for ( j=0 ; j<n ; j++ )
            {
              fjac[i*n+j] = jac[j*m+i];
            }
        }
      FREE(jac);
    }
  return 0;
}

/**
 * Compute the root of a function 
 * 
 * @param f a pointer to a PnlRnFuncRnDFunc. This an object for storing a
 * fonction f:R^n -> R^n and its Jacobian.
 * @param x on input, contains the initial value of the algorithm. On
 * output, contains the solution
 * @param fx on output, contains the value of f at x, ideally it should be
 * a vector of zeros
 * @param xtol the minimum relative error between two consecutive iterates.
 * If the lower bound if reached, the function returns 1
 * @param maxfev maximum number of evaluation of f. If the number is
 * reached, the function returns with the value 2
 * @param nfev number of evaluations of f in the algorithm
 * @param scale a vector of scale factors to make the components of the
 * solution roughly of the same order, once they have been scaled
 * @param error_msg a boolean TRUE or FALSE. If TRUE, a message is printed
 * if the hybrd function did not return properly
 *
 * @return OK or FAIL (if FAIL, use error_msg=TRUE to know what happened)
 */
int pnl_root_fsolve (PnlRnFuncRnDFunc *f, PnlVect *x, PnlVect *fx,  double xtol, 
                int maxfev, int *nfev, PnlVect *scale, int error_msg)
{

  double *wa, *fjac;
  double eps, epsfcn, factor;
  int with_jac, info, mode, nprint;
  int i, n, ml, lr, lwa, mu, njev;

  eps = pnl_minpack_dpmpar (1);
  epsfcn = 0.;
  n = x->size;
  factor = 100.;

  /*
   * Some pre-treatment on the parameters
   */
  wa = fjac = NULL;
  pnl_vect_resize (fx, n);
  with_jac = ( f->DF == NULL ) ? 0 : 1;
  if ( xtol <= 0 )
    {
      xtol = sqrt (2 * eps);
    }
  if ( maxfev <=0 )
    {
      maxfev = with_jac ? 100 * (n + 1) : 200 * (n + 1);
    }

  /*
   * Beginning of the routine
   */
  if ( with_jac )
    {
      lr = n * (n + 1) / 2;
      lwa = (n * (n + 13)) / 2;
      mode = 1; /* By default, internal scaling is used */
      if ((fjac = MALLOC_DOUBLE(n*n)) == NULL) return -1;
      if ((wa = MALLOC_DOUBLE(lwa)) == NULL) return -1;

      if ( scale != NULL )
        {
          /* We are using scaling */
          for ( i=0 ; i<n ; i++ ) wa[i] = PNL_GET(scale, i); 
          mode = 2;
        }
      nprint = 0;
      info = pnl_minpack_hybrj(hybrj_fcn, f, n, x->array, fx->array, fjac, n, xtol, maxfev, 
                               wa, mode, factor, nprint, nfev, &njev, &wa[6*n], lr,
                               &wa[n], &wa[2*n], &wa[3*n], &wa[4*n], &wa[5*n]);
      /*
       * tweak info
       */
      if ( info == 5 ) info=4;
      FREE(wa); FREE(fjac);
    }
  else
    {
      ml = mu = n - 1;
      lr = n * (n + 1) / 2;
      lwa = (n * (3 * n + 13)) / 2;
      mode = 1; /* By default, internal scaling is used */
      if ((wa = MALLOC_DOUBLE(lwa)) == NULL) return -1;

      if ( scale != NULL )
        {
          /* We are using scaling */
          for ( i=0 ; i<n ; i++ ) wa[i] = PNL_GET(scale, i); 
          mode = 2;
        }
      nprint = 0;
      info = pnl_minpack_hybrd(hybrd_fcn, f, n, x->array, fx->array, xtol, maxfev, 
                               ml, mu, epsfcn, wa, mode, factor,
                               nprint, nfev, &wa[6*n+lr], n, &wa[6*n], lr,
                               &wa[n], &wa[2*n], &wa[3*n], &wa[4*n], &wa[5*n]);
      FREE(wa);
    }

  /*
   * Error treatments
   */
  if ( error_msg == FALSE )
    {
      if ( info == 1 ) return OK; else return FAIL;
    }
  switch ( info )
    {
      case 0:
        printf ("Improper input parameters.\n");
        break;
      case 2: 
        printf("Relative error between two consecutive iterates is at most xtol.\n");
        break;
      case 3:
        printf("Number of calls to fcn has reached or exceeded maxfev=%d.\n", maxfev);
        break;
      case 4:
        printf("Iteration is not making good progress, as measured by the\nimprovement from the last five jacobian evaluations.\n");
        break;
      case 5:
        printf("Iteration is not making good progress, as measured by the\nimprovement from the last ten iterations.\n");
        break;
      default:
        if ( info < 0 ) 
          {
            printf ("Execution was aborted.\n");
          }
    }
  if ( info == 1 ) return OK; else return FAIL;

}

/**
 * Compute the root of  a sum of squares
 * 
 * @param f a pointer to a PnlRnFuncRmDFunc. This an object for storing a
 * fonction f:R^n -> R^m and its Jacobian.
 * @param x on input, contains the initial value of the algorithm. On
 * output, contains the solution.
 * @param m is the number of equations.
 * @param fx on output, contains the value of f at x, ideally it should be
 * a vector of zeros
 * @param xtol is the minimum relative error between two consecutive iterates.
 * If the lower bound if reached, the function returns 1
 * @param ftol minimum relative error in the computation of the sum
 * @param gtol minimum relative error in the computation of the gradient of the sum
 * @param maxfev maximum number of evaluation of f. If the number is
 * reached, the function returns with the value 2
 * @param nfev number of evaluations of f in the algorithm
 * @param scale a vector of scale factors to make the components of the
 * solution roughly of the same order, once they have been scaled
 * @param error_msg a boolean TRUE or FALSE. If TRUE, a message is printed
 * if the lmdif or lmder function did not return properly
 *
 * @return OK or FAIL (if FAIL, use error_msg=TRUE to know what happened)
 */
int pnl_root_fsolve_lsq (PnlRnFuncRmDFunc *f, PnlVect *x, int m, PnlVect *fx,  double xtol, 
                    double ftol, double gtol, int maxfev, int *nfev, 
                    PnlVect *scale, int error_msg)
{

  double *wa, *fjac, eps, epsfcn, factor;
  int    with_jac, info, mode, nprint, i, n, lwa, mp5n, njev;
  int    msg;

  eps = pnl_minpack_dpmpar (1);
  epsfcn = 0.;
  n = x->size;
  factor = 100.;
  msg = (error_msg == TRUE);

  /*
   * Some pre-treatment on the parameters
   */
  wa = fjac = NULL;
  pnl_vect_resize (fx, n);
  with_jac = ( f->DF == NULL ) ? 0 : 1;
  if ( xtol <= 0 ) { xtol = sqrt (2 * eps); }
  if ( gtol <= 0 ) { gtol = 0.; }
  if ( ftol <= 0 ) { ftol = sqrt (2 * eps); }
  if ( maxfev <=0 ) { maxfev = with_jac ? 100 * (n + 1) : 200 * (n + 1); }

  /*
   * Beginning of the routine
   */
  if ( with_jac )
    {
      int *ipvt;
      lwa = 5 * n + m;
      mode = 1; /* By default, internal scaling is used */
      if ((fjac = MALLOC_DOUBLE(m*n)) == NULL) return -1;
      if ((wa = MALLOC_DOUBLE(lwa)) == NULL) return -1;
      if ((ipvt = MALLOC_INT(n)) == NULL) return -1;

      if ( scale != NULL )
        {
          /* We are using scaling */
          for ( i=0 ; i<n ; i++ ) wa[i] = PNL_GET(scale, i); 
          mode = 2;
        }
      nprint = 0;
      info = pnl_minpack_lmder(lmder_fcn, f, m, n, x->array, fx->array,
                               fjac, m, ftol, xtol, gtol, maxfev, wa, mode,
                               factor, nprint, nfev, &njev, ipvt, &wa[n],
                               &wa[2*n], &wa[3*n], &wa[4*n], &wa[5*n]);
      FREE(fjac); FREE(wa); FREE(ipvt);

    }
  else
    {
      int *iwa;
      lwa = m* n + 5 * n + m;
      mp5n = m + n * 5;
      mode = 1; /* By default, internal scaling is used */
      if ((wa = MALLOC_DOUBLE(lwa)) == NULL) return -1;
      if ((iwa = MALLOC_INT(n)) == NULL) return -1;

      if ( scale != NULL )
        {
          /* We are using scaling */
          for ( i=0 ; i<n ; i++ ) wa[i] = PNL_GET(scale, i); 
          mode = 2;
        }
      nprint = 0;
      info = pnl_minpack_lmdif(lmdif_fcn, f, m, n, x->array, fx->array, ftol, 
                               xtol, gtol, maxfev, epsfcn, wa, mode, factor,
                               nprint, nfev, &wa[mp5n], m, iwa, &wa[n], 
                               &wa[2*n], &wa[3*n], &wa[4*n], &wa[5*n]);
      FREE(wa); FREE(iwa);
    }

  /*
   * Error treatments
   */
  switch ( info )
    {
    case 0:
      if (msg) printf ("Improper input parameters.\n");
      return FAIL;
      break;
    case 1:
      if (msg) printf("Both actual and predicted relative reductions in the sum\nof squares are at most ftol.\n");
      return OK;
      break;
    case 2:  
      if (msg) printf ("Relative error between two consecutive iterates is at most xtol.\n");
      return OK;
      break;
    case 3: 
      if (msg) printf("Both actual and predicted relative reductions in the sum\nof squares are at most ftol.\n");
      if (msg) printf ("Relative error between two consecutive iterates is at most xtol.\n");
      return OK;
      break;
    case 4: 
      if (msg) printf("The cosine of the angle between fvec and any\ncolumn of the jacobian is at most gtol in\nabsolute value.\n");
      return OK;
      break;
    case 5:
      if (msg) printf("Number of calls to fcn has reached or exceeded maxfev=%d.\n", maxfev);
      return FAIL;
      break;
    case 6: 
      if (msg) printf("ftol is too small. No further reduction in\nthe sum of squares is possible.\n");
      return OK;
      break;
    case 7: 
      if (msg) printf("xtol is too small. No further improvement in\nthe approximate solution x is possible.\n");
      return OK;
      break;
    case 8: 
      if (msg) printf ("gtol is too small. fvec is orthogonal to the\ncolumns of the jacobian to machine precision.\n");
      return OK;
      break;

    default:
      if ( info < 0 ) 
        {
          printf ("Execution was aborted.\n");
        }
      return FAIL;
    }

}

