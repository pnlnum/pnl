/*                                                      hyp2f1.c
 *
 *      Gauss hypergeometric function   F
 *                                     2 1
 *
 *
 * SYNOPSIS:
 *
 * double a, b, c, x, y, hyp2f1();
 *
 * y = hyp2f1( a, b, c, x );
 *
 *
 * DESCRIPTION:
 *
 *
 *  hyp2f1( a, b, c, x )  =   F ( a, b; c; x )
 *                           2 1
 *
 *           inf.
 *            -   a(a+1)...(a+k) b(b+1)...(b+k)   k+1
 *   =  1 +   >   -----------------------------  x   .
 *            -         c(c+1)...(c+k) (k+1)!
 *          k = 0
 *
 *  Cases addressed are
 *      Tests and escapes for negative integer a, b, or c
 *      Linear transformation if c - a or c - b negative integer
 *      Special case c = a or c = b
 *      Linear transformation for  x near +1
 *      Transformation for x < -0.5
 *      Psi function expansion if x > 0.5 and c - a - b integer
 *      Conditionally, a recurrence on c to make c-a-b > 0
 *
 *      x < -1  AMS 15.3.7 transformation applied (Travis Oliphant)
 *         valid for b,a,c,(b-a) != integer and (c-a),(c-b) != negative integer
 *
 * x >= 1 is rejected (unless special cases are present)
 *
 * The parameters a, b, c are considered to be integer
 * valued if they are within 1.0e-14 of the nearest integer
 * (1.0e-13 for IEEE arithmetic).
 *
 * ACCURACY:
 *
 *
 *               Relative error (-1 < x < 1):
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -1,7        230000      1.2e-11     5.2e-14
 *
 * Several special cases also tested with a, b, c in
 * the range -7 to 7.
 *
 * ERROR MESSAGES:
 *
 * A "partial loss of precision" message is printed if
 * the internally estimated relative error exceeds 1^-12.
 * A "singularity" message is printed on overflow or
 * in cases not addressed (such as x < -1).
 */

/*                                                      hyp2f1  */


/*
  Cephes Math Library Release 2.8:  June, 2000
  Copyright 1984, 1987, 1992, 2000 by Stephen L. Moshier
*/


#include "mconf.h"

#ifdef DEC
#define EPS 1.0e-14
#define EPS2 1.0e-11
#endif

#ifdef IBMPC
#define EPS 1.0e-13
#define EPS2 1.0e-10
#endif

#ifdef MIEEE
#define EPS 1.0e-13
#define EPS2 1.0e-10
#endif

#ifdef UNK
#define EPS 1.0e-13
#define EPS2 1.0e-10
#endif

#define ETHRESH 1.0e-12

extern double psi ( double );
static double hyt2f1(double, double, double, double, double *);
static double hys2f1(double, double, double, double, double *);
double hyp2f1(double, double, double, double);

double hyp2f1(double a,double b,double c, double x)
{
  double d, d1, d2, e;
  double p, q, r, s, y, ax;
  double ia, ib, ic, id, err;
  int i, aid;
  int neg_int_a = 0, neg_int_b = 0;
  int neg_int_ca_or_cb = 0;

  err = 0.0;
  ax = fabs(x);
  s = 1.0 - x;
  ia = round(a); /* nearest integer to a */
  ib = round(b);

  if (x == 0.0) {
    return 1.0;
  }

  d = c - a - b;
  if (d <= -1) {
    return pow(s, d) * hyp2f1(c-a, c-b, c, x);
  }
  if (d <= 0 && x == 1)
    goto hypdiv;

  if (a <= 0 && fabs(a-ia) < EPS ) {      /* a is a negative integer */
    neg_int_a = 1;
  }

  if (b <= 0 && fabs(b-ib) < EPS ) {      /* b is a negative integer */
    neg_int_b = 1;
  }

  if (ax < 1.0 || x == -1.0) {
    /* 2F1(a,b;b;x) = (1-x)**(-a) */
    if( fabs(b-c) < EPS ) {         /* b = c */
      y = pow( s, -a );       /* s to the -a power */
      goto hypdon;
    }
    if( fabs(a-c) < EPS ) {         /* a = c */
      y = pow( s, -b );       /* s to the -b power */
      goto hypdon;
    }
  }



  if( c <= 0.0 )
    {
      ic = round(c);  /* nearest integer to c */
      if( fabs(c-ic) < EPS )          /* c is a negative integer */
        {
          /* check if termination before explosion */
          if( neg_int_a && (ia > ic) )
            goto hypok;
          if( neg_int_b && (ib > ic) )
            goto hypok;
          goto hypdiv;
        }
    }

  if (neg_int_a || neg_int_b)             /* function is a polynomial */
    goto hypok;

  if (x < -1.0) {
    double t1;
    t1 = fabs(b - a);
    if (fabs(t1 - round(t1)) < EPS) {
      /* this transformation has a pole for b-a= +-integer,
         so we average around it.
      */
      return 0.5*(hyp2f1(a, b*(1+1e-9), c, x) + hyp2f1(a, b*(1-1e-9), c, x));
    }
    p = hyp2f1(a, 1-c+a, 1-b+a, 1.0/x);
    q = hyp2f1(b, 1-c+b, 1-a+b, 1.0/x);
    p *= pow(-x, -a);
    q *= pow(-x, -b);
    t1 = pnl_sf_gamma(c);
    s = t1*pnl_sf_gamma(b-a)/(pnl_sf_gamma(b)*pnl_sf_gamma(c-a));
    y = t1*pnl_sf_gamma(a-b)/(pnl_sf_gamma(a)*pnl_sf_gamma(c-b));
    return s*p + y*q;
  }

  if( ax > 1.0 )                  /* series diverges      */
    goto hypdiv;

  p = c - a;
  ia = round(p); /* nearest integer to c-a */
  if( (ia <= 0.0) && (fabs(p-ia) < EPS) ) /* negative int c - a */
    neg_int_ca_or_cb = 1;

  r = c - b;
  ib = round(r); /* nearest integer to c-b */
  if( (ib <= 0.0) && (fabs(r-ib) < EPS) ) /* negative int c - b */
    neg_int_ca_or_cb = 1;

  id = round(d); /* nearest integer to d */
  q = fabs(d-id);

  /* Thanks to Christian Burger <BURGER@DMRHRZ11.HRZ.Uni-Marburg.DE>
   * for reporting a bug here.  */
  if( fabs(ax-1.0) < EPS ) {                      /* |x| == 1.0   */
    if( x > 0.0 ) {
      if (neg_int_ca_or_cb) {
        if( d >= 0.0 )
          goto hypf;
        else
          goto hypdiv;
      }
      if( d <= 0.0 )
        goto hypdiv;
      y = pnl_sf_gamma(c)*pnl_sf_gamma(d)/(pnl_sf_gamma(p)*pnl_sf_gamma(r));
      goto hypdon;
    }
    if( d <= -1.0 )
      goto hypdiv;
  }

  /* Conditionally make d > 0 by recurrence on c
   * AMS55 #15.2.27
   */
  if( d < 0.0 )
    {
      /* Try the power series first */
      y = hyt2f1( a, b, c, x, &err );
      if( err < ETHRESH )
        goto hypdon;
      /* Apply the recurrence if power series fails */
      err = 0.0;
      aid = 2 - id;
      e = c + aid;
      d2 = hyp2f1(a,b,e,x);
      d1 = hyp2f1(a,b,e+1.0,x);
      q = a + b + 1.0;
      for( i=0; i<aid; i++ )
        {
          r = e - 1.0;
          y = (e*(r-(2.0*e-q)*x)*d2 + (e-a)*(e-b)*x*d1)/(e*r*s);
          e = r;
          d1 = d2;
          d2 = y;
        }
      goto hypdon;
    }


  if (neg_int_ca_or_cb)
    goto hypf; /* negative integer c-a or c-b */

 hypok:
  y = hyt2f1( a, b, c, x, &err );


 hypdon:
  if( err > ETHRESH )
    {
      mtherr( "hyp2f1", PLOSS );
      /*      printf( "Estimated err = %.2e\n", err ); */
    }
  return(y);

  /* The transformation for c-a or c-b negative integer
   * AMS55 #15.3.3
   */
 hypf:
  y = pow( s, d ) * hys2f1( c-a, c-b, c, x, &err );
  goto hypdon;

  /* The alarm exit */
 hypdiv:
  mtherr( "hyp2f1", OVERFLOW );
  return( MAXNUM );
}






/* Apply transformations for |x| near 1
 * then call the power series
 */
static double hyt2f1(double a,double b,double c,double x, double *loss)
{
  double p, q, r, s, t, y, d, err, err1;
  double ax, id, d1, d2, e, y1;
  int i, aid;

  err = 0.0;
  s = 1.0 - x;
  if( x < -0.5 )
    {
      if( b > a )
        y = pow( s, -a ) * hys2f1( a, c-b, c, -x/s, &err );

      else
        y = pow( s, -b ) * hys2f1( c-a, b, c, -x/s, &err );

      goto done;
    }

  d = c - a - b;
  id = round(d);  /* nearest integer to d */

  if( x > 0.9 )
    {
      if( fabs(d-id) > EPS ) /* test for integer c-a-b */
        {
          /* Try the power series first */
          y = hys2f1( a, b, c, x, &err );
          if( err < ETHRESH )
            goto done;
          /* If power series fails, then apply AMS55 #15.3.6 */
          q = hys2f1( a, b, 1.0-d, s, &err );     
          q *= pnl_sf_gamma(d) /(pnl_sf_gamma(c-a) * pnl_sf_gamma(c-b));
          r = pow(s,d) * hys2f1( c-a, c-b, d+1.0, s, &err1 );
          r *= pnl_sf_gamma(-d)/(pnl_sf_gamma(a) * pnl_sf_gamma(b));
          y = q + r;

          q = fabs(q); /* estimate cancellation error */
          r = fabs(r);
          if( q > r )
            r = q;
          err += err1 + (MACHEP*r)/y;

          y *= pnl_sf_gamma(c);
          goto done;
        }
      else
        {
          /* Psi function expansion, AMS55 #15.3.10, #15.3.11, #15.3.12 */
          if( id >= 0.0 )
            {
              e = d;
              d1 = d;
              d2 = 0.0;
              aid = id;
            }
          else
            {
              e = -d;
              d1 = 0.0;
              d2 = d;
              aid = -id;
            }

          ax = log(s);

          /* sum for t = 0 */
          y = pnl_sf_psi(1.0) + pnl_sf_psi(1.0+e) - pnl_sf_psi(a+d1) - pnl_sf_psi(b+d1) - ax;
          y /= pnl_sf_gamma(e+1.0);

          p = (a+d1) * (b+d1) * s / pnl_sf_gamma(e+2.0); /* Poch for t=1 */
          t = 1.0;
          do
            {
              r = pnl_sf_psi(1.0+t) + pnl_sf_psi(1.0+t+e) - pnl_sf_psi(a+t+d1)
                - pnl_sf_psi(b+t+d1) - ax;
              q = p * r;
              y += q;
              p *= s * (a+t+d1) / (t+1.0);
              p *= (b+t+d1) / (t+1.0+e);
              t += 1.0;
            }
          while( fabs(q/y) > EPS );


          if( id == 0.0 )
            {
              y *= pnl_sf_gamma(c)/(pnl_sf_gamma(a)*pnl_sf_gamma(b));
              goto psidon;
            }

          y1 = 1.0;

          if( aid == 1 )
            goto nosum;

          t = 0.0;
          p = 1.0;
          for( i=1; i<aid; i++ )
            {
              r = 1.0-e+t;
              p *= s * (a+t+d2) * (b+t+d2) / r;
              t += 1.0;
              p /= t;
              y1 += p;
            }
        nosum:
          p = pnl_sf_gamma(c);
          y1 *= pnl_sf_gamma(e) * p / (pnl_sf_gamma(a+d1) * pnl_sf_gamma(b+d1));

          y *= p / (pnl_sf_gamma(a+d2) * pnl_sf_gamma(b+d2));
          if( (aid & 1) != 0 )
            y = -y;

          q = pow( s, id );       /* s to the id power */
          if( id > 0.0 )
            y *= q;
          else
            y1 *= q;

          y += y1;
        psidon:
          goto done;
        }

    }

  /* Use defining power series if no special cases */
  y = hys2f1( a, b, c, x, &err );

 done:
  *loss = err;
  return(y);
}





/* Defining power series expansion of Gauss hypergeometric function */

static double hys2f1(double a,double b,double c,double x,double *loss )
{
  double f, g, h, k, m, s, u, umax;
  int i;

  i = 0;
  umax = 0.0;
  f = a;
  g = b;
  h = c;
  s = 1.0;
  u = 1.0;
  k = 0.0;
  do
    {
      if( fabs(h) < EPS )
        {
          *loss = 1.0;
          return( MAXNUM );
        }
      m = k + 1.0;
      u = u * ((f+k) * (g+k) * x / ((h+k) * m));
      s += u;
      k = fabs(u);  /* remember largest term summed */
      if( k > umax )
        umax = k;
      k = m;
      if( ++i > 10000 ) /* should never happen */
        {
          *loss = 1.0;
          return(s);
        }
    }
  while( fabs(u/s) > MACHEP );

  /* return estimated relative error */
  *loss = (MACHEP*umax)/fabs(s) + (MACHEP*i);

  return(s);
}
