/*                                                      beta.c
 *
 *      Beta function
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, b, y, beta();
 *
 * y = beta( a, b );
 *
 *
 *
 * DESCRIPTION:
 *
 *                   -     -
 *                  | (a) | (b)
 * beta( a, b )  =  -----------.
 *                     -
 *                    | (a+b)
 *
 * For large arguments the logarithm of the function is
 * evaluated using lgam(), then exponentiated.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC        0,30        1700       7.7e-15     1.5e-15
 *    IEEE       0,30       30000       8.1e-14     1.1e-14
 *
 * ERROR MESSAGES:
 *
 *   message         condition          value returned
 * beta overflow    log(beta) > MAXLOG       0.0
 *                  a or b <0 integer        0.0
 *
 */

/*                                                      beta.c  */


/*
  Cephes Math Library Release 2.0:  April, 1987
  Copyright 1984, 1987 by Stephen L. Moshier
  Direct inquiries to 30 Frost Street, Cambridge, MA 02140

  Modified by Jérôme Lelong <jerome.lelong@gmail.com> April, 2011
  to make the code thread-safe (global variable sgngam removed)
*/

#include "mconf.h"

#ifdef UNK
#define MAXGAM 34.84425627277176174
#endif
#ifdef DEC
#define MAXGAM 34.84425627277176174
#endif
#ifdef IBMPC
#define MAXGAM 171.624376956302725
#endif
#ifdef MIEEE
#define MAXGAM 171.624376956302725
#endif

extern int pnl_sf_log_gamma_sgn(double, double *, int *);

double beta( double a, double b )
{
  double y;
  int sign;

  sign = 1;

  if( a <= 0.0 )
    {
      if( a == floor(a) )
        goto over;
    }
  if( b <= 0.0 )
    {
      if( b == floor(b) )
        goto over;
    }


  y = a + b;
  if( fabs(y) > MAXGAM )
    {
      double lg;
      int sgngam;
      pnl_sf_log_gamma_sgn(y, &lg, &sgngam);
      y = lg;
      sign *= sgngam; /* keep track of the sign */
      pnl_sf_log_gamma_sgn(b, &lg, &sgngam);
      y = lg - y;
      sign *= sgngam;
      pnl_sf_log_gamma_sgn(a, &lg, &sgngam);
      y = lg + y;
      sign *= sgngam;
      if( y > MAXLOG )
        {
        over:
          mtherr( "beta", OVERFLOW );
          return( sign * MAXNUM );
        }
      return( sign * exp(y) );
    }

  y = pnl_sf_gamma(y);
  if( y == 0.0 )
    goto over;

  if( a > b )
    {
      y = pnl_sf_gamma(a)/y;
      y *= pnl_sf_gamma(b);
    }
  else
    {
      y = pnl_sf_gamma(b)/y;
      y *= pnl_sf_gamma(a);
    }

  return(y);
}



/* Natural log of |beta|.  Return the sign of beta in sgngam.  */

int pnl_sf_lbeta (double a, double b, double *res, int *sgn)
{
  double y;
  int sign;

  sign = 1;

  if( a <= 0.0 )
    {
      if( a == floor(a) )
        goto over;
    }
  if( b <= 0.0 )
    {
      if( b == floor(b) )
        goto over;
    }


  y = a + b;
  if( fabs(y) > MAXGAM )
    {
      double lg;
      pnl_sf_log_gamma_sgn (y, &lg, sgn);
      y = lg;
      sign *= *sgn; /* keep track of the sign */
      pnl_sf_log_gamma_sgn (b, &lg, sgn);
      y = lg - y;
      sign *= *sgn;
      pnl_sf_log_gamma_sgn (a, &lg, sgn);
      y = lg + y;
      sign *= *sgn;
      *sgn = sign;
      *res = y;
      return OK;
    }

  y = pnl_sf_gamma(y);
  if( y == 0.0 )
    {
    over:
      mtherr( "lbeta", OVERFLOW );
      return( sign * MAXNUM );
    }

  if( a > b )
    {
      y = pnl_sf_gamma(a)/y;
      y *= pnl_sf_gamma(b);
    }
  else
    {
      y = pnl_sf_gamma(b)/y;
      y *= pnl_sf_gamma(a);
    }

  if( y < 0 )
    {
      *sgn = -1;
      *res = log (-y);
    }
  else
    {
      *sgn = 1;
      *res = log (y);
    }
  return OK;
}

double lbeta( double a, double b )
{
  double res;
  int sgn;
  pnl_sf_lbeta (a, b, &res, &sgn);
  return res;
}


