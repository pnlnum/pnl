
/*************************************************************************/
/* Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>            */
/*                                                                       */
/* This program is free software; you can redistribute it and/or modify  */
/* it under the terms of the GNU Lesser General Public License as        */ 
/* published by  the Free Software Foundation; either version 3 of the   */
/* License, or (at your option) any later version.                       */
/*                                                                       */
/* This program is distributed in the hope that it will be useful,       */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/* GNU Lesser General Public License for more details.                   */
/*                                                                       */
/* You should have received a copy of the GNU Lesser General Public      */
/* License  along with this program.  If not, see                        */
/* <http://www.gnu.org/licenses/>.                                       */ 
/*************************************************************************/

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_internals.h"


/*
 * Some external declarations from QuadPack
 */

extern int pnl_dqng(PnlFunc *F, double *a, double *b, double *epsabs,
  double *epsrel, double *result, double *abserr, int *neval, int *ier);

extern int pnl_dqagie(PnlFunc *f, double *bound, int *inf, 
	double *epsabs, double *epsrel, int *limit, double *
	result, double *abserr, int *neval, int *ier, double *
	alist__, double *blist, double *rlist, double *elist, 
  int *iord, int *last);

extern int pnl_dqagse(PnlFunc *f, double *a, double *b, double 
	*epsabs, double *epsrel, int *limit, double *result, 
	double *abserr, int *neval, int *ier, double *alist__,
	 double *blist, double *rlist, double *elist, int *
   iord, int *last);

extern int pnl_dqagpe(PnlFunc * f, double *a, double *b, int *
	npts2, double *points, double *epsabs, double *epsrel, 
	int *limit, double *result, double *abserr, int *
	neval, int *ier, double *alist__, double *blist, 
	double *rlist, double *elist, double *pts, int *iord, 
  int *level, int *ndin, int *last);


/*
 * Some static variables used for 2D integration
 */
static double xsav,xepsabs,xepsrel,y0sav,y1sav;
static const PnlFunc2D* globalfunc;
static int status, xneval;

/**
 * Integration over a finite interval using a non-adaptive Gauss Konrod
 * procedure
 *
 * @param *F a PnlFunc to be integrated 
 * @param a lower bound for integration (may be -Inf)
 * @param b upper bound for integration (may be -Inf)
 * @param epsabs maximum absolute error accepted
 * @param epsrel maximum relative error accepted
 * @param result the result of the integration
 * @param abserr the absolute error of the computation
 * @param neval number of function evaluations
 * @return  OK or FAIL if the required precision cannot be attained
 */
int pnl_integration_qng (PnlFunc *F,
                        double a, double b,
                        double epsabs, double epsrel,
                        double *result, double *abserr,
                        int *neval)
{
  int ier;
  pnl_dqng(F, &a, &b, &epsabs, &epsrel, result, abserr, neval, &ier);

  if ( ier != 0 ) return FAIL;
  return OK;
}

/**
 * Integration over a finite interval using a non-adaptive Gauss Konrod
 * procedure
 *
 * This function is a synonymous of pnl_integration_qng for backward
 * compatibility
 *
 * @param *F a PnlFunc to be integrated 
 * @param a lower bound for integration (may be -Inf)
 * @param b upper bound for integration (may be -Inf)
 * @param epsabs maximum absolute error accepted
 * @param epsrel maximum relative error accepted
 * @param result the result of the integration
 * @param abserr the absolute error of the computation
 * @param neval number of function evaluations
 * @return  OK or FAIL if the required precision cannot be attained
 */
int pnl_integration_GK (PnlFunc *F,
                        double a, double b,
                        double epsabs, double epsrel,
                        double *result, double *abserr,
                        int *neval)
{
  int ier;
  pnl_dqng(F, &a, &b, &epsabs, &epsrel, result, abserr, neval, &ier);

  if ( ier != 0 ) return FAIL;
  return OK;
}




/**
 * Integration over a finite or non finite interval
 *
 * @param *f a PnlFunc to be integrated 
 * @param a lower bound for integration (may be -Inf)
 * @param b upper bound for integration (may be -Inf)
 * @param epsabs maximum absolute error accepted
 * @param epsrel maximum relative error accepted
 * @param limit number of subintervals of [a, b] used during the computation
 * @param result the result of the integration
 * @param abserr the absolute error of the computation
 * @param neval number of function evaluations
 * @return  OK or FAIL if the required precision cannot be attained
 */
int pnl_integration_qag (PnlFunc *f, double a, double b, double epsabs,
                          double epsrel, int limit, double *result, 
                          double *abserr, int *neval)
{
  double *alist, *blist, *rlist, *elist;
  int *iord;
  int last, ier, infbounds = 0, sign = 1;
  double bound = 0.0;

  if ( a == b )
    {
      *result = 0.;
      *abserr = 0.;
      *neval  = 0;
      return OK;
    }
  if ( a > b )
    {
      double tmp = a;
      a = b;
      b = tmp;
      sign = -1;
    }
  /*
   * From now, a < b
   */

  if ( pnl_isfinite(b) && pnl_isfinite(a) )
    {
      infbounds = 0;
    }
  else 
    {
      if ( pnl_isfinite(a) ) /* b is +Inf */
        {
          infbounds = 1; 
          bound = a;
        }
      else if ( pnl_isfinite(b) ) /* a is -Inf */
        {
          infbounds = -1;
          bound = b;
        }
      else /* a is -Inf and b is +Inf */
        {
          infbounds = 2;
        }
    }
  
  if ( limit == 0 ) limit = 750;

  /* allocate some arrays internally used by dqagse */
  alist = MALLOC_DOUBLE(limit);
  blist = MALLOC_DOUBLE(limit);
  rlist = MALLOC_DOUBLE(limit);
  elist = MALLOC_DOUBLE(limit);
  iord  = MALLOC_INT(limit); 

  if ( alist == NULL || blist == NULL || rlist == NULL ||
      elist == NULL || iord == NULL )
    {
      FREE (alist); FREE (blist); FREE (rlist); FREE (elist);
      FREE (iord);
    }

  if (infbounds == 0)
    {
      pnl_dqagse (f, &a, &b, &epsabs, &epsrel, &limit, result, abserr, neval, &ier, 
                  alist, blist, rlist, elist, iord, &last);
    }
  else
    {
      pnl_dqagie (f, &bound, &infbounds, &epsabs, &epsrel, &limit, result, 
                  abserr, neval, &ier, alist, blist, rlist, elist, iord, &last);
    }
  if ( sign == -1 ) *result = - (*result);
  FREE (alist); FREE (blist); FREE (rlist); FREE (elist);
  FREE (iord);

  if ( ier != 0 ) return FAIL;
  return OK;
}


/**
 * Integration over a finite or non finite interval of a function with known
 * singular points
 *
 * @param *f a PnlFunc to be integrated 
 * @param a lower bound for integration (may be -Inf)
 * @param b upper bound for integration (may be -Inf)
 * @param singularities a vector containing the singular points of f over the
 * interval (a,b). Note that a and b should not part of this vector
 * @param epsabs maximum absolute error accepted
 * @param epsrel maximum relative error accepted
 * @param limit number of subintervals of [a, b] used during the computation
 * @param result the result of the integration
 * @param abserr the absolute error of the computation
 * @param neval number of function evaluations
 * @return  OK or FAIL if the required precision cannot be attained
 */
int pnl_integration_qagp (PnlFunc *f, double a, double b, const PnlVect *singularities,
                          double epsabs, double epsrel, int limit, double *result, 
                          double *abserr, int *neval)
{
  double *alist, *blist, *rlist, *elist, *points, *pts;
  int *iord, *level, *ndim;
  int i, last, ier, npts2, sign = 1;

  if ( a == b )
    {
      *result = 0.;
      *abserr = 0.;
      *neval  = 0;
      return OK;
    }
  if ( a > b )
    {
      double tmp = a;
      a = b;
      b = tmp;
      sign = -1;
    }
  /*
   * From now, a < b
   */
  npts2 = singularities->size + 2;
  if ( limit == 0 ) limit = 750;

  /* allocate some arrays internally used by dqagse */
  alist  = MALLOC_DOUBLE(limit);
  blist  = MALLOC_DOUBLE(limit);
  rlist  = MALLOC_DOUBLE(limit);
  elist  = MALLOC_DOUBLE(limit);
  level  = MALLOC_INT(limit);
  iord   = MALLOC_INT(limit);
  points = MALLOC_DOUBLE(npts2);
  pts    = MALLOC_DOUBLE(npts2);
  ndim   = MALLOC_INT(npts2);

  if ( alist == NULL || blist == NULL || rlist == NULL || elist == NULL ||
       level == NULL || iord == NULL || pts == NULL || ndim == NULL ||
       points == NULL )
    {
      FREE(alist); FREE(blist); FREE(rlist); FREE(elist);
      FREE(level); FREE(iord); FREE(pts); FREE(ndim);
      FREE(points);
    }

  for ( i=0 ; i<npts2-2 ; i++ ) { points[i] = GET(singularities, i); }

  pnl_dqagpe (f, &a, &b, &npts2, points, &epsabs, &epsrel, &limit, 
              result, abserr, neval, &ier, alist, blist, rlist, 
              elist, pts, iord, level, ndim, &last);

  if ( sign == -1 ) *result = - (*result);
  FREE (alist); FREE (blist); FREE (rlist); FREE (elist);
  FREE(level); FREE(iord); FREE(pts); FREE(ndim);
  FREE (points);

  if ( ier != 0 ) return FAIL;
  return OK;
}

/*
 * 2D integration
 */

static double func1D(double y, void *params) 
{
  return globalfunc->F(xsav,y, params);
}

/*
 * A 1D wrapper for 2D integration
 */
static double int_1d(double x, void *params) 
{
  PnlFunc wrap_1d;
  double res,err;
  int it;
  wrap_1d.F= func1D;
  wrap_1d.params = params;
  xsav=x;
  if (pnl_integration_qng (&wrap_1d,y0sav,y1sav,xepsabs,xepsrel,&res,&err,&it)!=OK)
    {
      status = FAIL;
      return 0.;
    }
  xneval += it;
  return res;
}

/**
 * Integration over a two dimensional finite rectangle
 *
 * @param *F a PnlFunc2D to be integrated 
 * @param x0 lower left hand side of the domain
 * @param x1 lower right hand side of the domain
 * @param y0 upper left hand side of the domain
 * @param y1 upper right hand side of the domain
 * @param epsabs tolerance on the absolute error
 * @param epsrel tolerance on the relative error
 * @param result results of integration
 * @param abserr error of integration
 * @param neval number of function evaluations
 * @return  OK or FAIL if the required precision cannot be achieved
 */
int pnl_integration_qng_2d (PnlFunc2D *F,
                          double x0, double x1,
                          double y0,double y1,
                          double epsabs, double epsrel,
                          double * result, double * abserr,
                          int * neval)
{
  PnlFunc func_1d;
  func_1d.F = &int_1d;
  func_1d.params = NULL;
  globalfunc=F;
  xepsabs=epsabs;
  xepsrel=epsrel;
  y0sav=y0;
  y1sav=y1;
  status = OK;
  xneval = 0;
  if (pnl_integration_GK (&func_1d,x0,x1,epsabs,epsrel,result,abserr,neval) == FAIL)
    return FAIL;
  else 
    {
      *neval = xneval;
      return status;
    }
}


/**
 * Integration over a two dimensional finite rectangle
 *
 * This function is a synonymous of pnl_integration_qng_2d for backward
 * compatibility
 *
 * @param *F a PnlFunc2D to be integrated 
 * @param x0 lower left hand side of the domain
 * @param x1 lower right hand side of the domain
 * @param y0 upper left hand side of the domain
 * @param y1 upper right hand side of the domain
 * @param epsabs tolerance on the absolute error
 * @param epsrel tolerance on the relative error
 * @param result results of integration
 * @param abserr error of integration
 * @param neval number of function evaluations
 * @return  OK or FAIL if the required precision cannot be achieved
 */
int pnl_integration_GK2D (PnlFunc2D *F,
                          double x0, double x1,
                          double y0,double y1,
                          double epsabs, double epsrel,
                          double * result, double * abserr,
                          int * neval)
{
  return pnl_integration_qng_2d (F, x0, x1, y0, y1, epsabs, epsrel, 
                                 result, abserr, neval);
}

