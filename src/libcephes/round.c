/*                                                      round.c
 *
 *      Round double to nearest or even integer valued double
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, round();
 *
 * y = round(x);
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the nearest integer to x as a double precision
 * floating point result.  If x ends in 0.5 exactly, the
 * nearest even integer is chosen.
 * 
 *
 *
 * ACCURACY:
 *
 * If x is greater than 1/(2*MACHEP), its closest machine
 * representation is already an integer, so rounding does
 * not change it.
 */

/*
  Cephes Math Library Release 2.1:  January, 1989
  Copyright 1984, 1987, 1989 by Stephen L. Moshier
  Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include "mconf.h"

double pnl_round(double x)
{
  double y, r;

  /* Largest integer <= x */
  y = floor(x);

  /* Fractional part */
  r = x - y;

  if( r > 0.5 )
    {
      /* Round up to nearest. */
      y += 1.;
    }
  else if( r == 0.5 )
    {
      /* Round to even */
      r = y - 2.0 * floor( 0.5 * y );
      if( r == 1.0 )
        {
          y += 1.0;
        }
    }
  /* Else round down. */

  return(y);
}
