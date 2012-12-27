/*  
 *  This file is a collection of wrappers around the
 *  Amos Fortran library of functions that take complex
 *  variables (see www.netlib.org) so that they can be called from
 *  the cephes library of corresponding name but work with complex
 *  arguments.
 * 
 *  Author:  Travis E. Oliphant
 *            oliphant.travis@altavista.net
 *
 * This can be freely distributed provided this notification remain and it is
 * understood that no warranty is expressed or implied.
 *
 *
 * These wrappers have been adapted to Nsp by Jean-Philippe Chancelier (Dec
 * 2007) and later modified by Jerome Lelong (Jan 2009) to deal with negative
 * orders.
 *
 */

#include <stdlib.h>
#include <math.h>
#include "amos.h"
#include "amos_wrappers.h"
#include "pnl/pnl_specfun.h"

static int print_error_messages = 1;

/* Notice: the order of appearance of the following
 * messages is bound to the error codes defined
 * in amos_wrappers.h.
 */
static char *ermsg[8] = {
  "unknown",      /* error code 0 */
  "domain",       /* error code 1 */
  "singularity",  /* et seq.      */
  "overflow",
  "underflow",
  "total loss of precision",
  "partial loss of precision",
  "too many iterations"
};


static int mtherr(char *name, int code)
{

  /* Display string passed by calling program,
   * which is supposed to be the name of the
   * function in which the error occurred:
   */

  /* Display error message defined
   * by the code argument.
   */
  if( (code <= 0) || (code >= 8) )
    code = 0;
  if (print_error_messages) {
    printf( "\n%s ", name );
    printf( "%s error\n", ermsg[code] );
  }

  /* Return to calling
   * program
   */
  return( 0 );
}

static int ierr_to_mtherr( int nz, int ierr) 
{
  /* Return mtherr equivalents for ierr values */
  if (nz != 0) return UNDERFLOW;
  switch (ierr) {
  case 1:
    return DOMAIN;
  case 2:
    return OVERFLOW;
  case 3:
    return PLOSS;
  case 4:
    return TLOSS;
  case 5:   /* Algorithm termination condition not met */
    return TLOSS;    
  }
  return 0;
}

/* int cairy_wrap(dcomplex z, dcomplex *ai, dcomplex *aip, dcomplex *bi, dcomplex *bip) 
 * {
 *   int id = 0;
 *   int ierr = 0;
 *   int kode = 1;
 *   int nz;
 * 
 *   pnl_zairy(CADDR(z), &id, &kode, F2C_CST(ai), &nz, &ierr);
 *   DO_MTHERR("airy:");
 *   pnl_zbiry(CADDR(z), &id, &kode, F2C_CST(bi),  &ierr);
 *   DO_MTHERR("airy:");
 *   id = 1;
 *   pnl_zairy(CADDR(z), &id, &kode, F2C_CST(aip), &nz, &ierr);
 *   DO_MTHERR("airy:");
 *   pnl_zbiry(CADDR(z), &id, &kode, F2C_CST(bip),  &ierr);
 *   DO_MTHERR("airy:");
 *   return 0;
 * }
 * 
 * int cairye_wrap(dcomplex z, dcomplex *ai, dcomplex *aip, dcomplex *bi, dcomplex *bip) 
 * {
 *   int id = 0;
 *   int kode = 2;        /\* Exponential scaling *\/
 *   int nz, ierr;
 * 
 *   pnl_zairy(CADDR(z), &id, &kode, F2C_CST(ai), &nz, &ierr);
 *   DO_MTHERR("airye:");
 *   pnl_zbiry(CADDR(z), &id, &kode, F2C_CST(bi),  &ierr);
 *   DO_MTHERR("airye:");
 *   
 *   id = 1;
 *   pnl_zairy(CADDR(z), &id, &kode, F2C_CST(aip), &nz, &ierr);
 *   DO_MTHERR("airye:");
 *   pnl_zbiry(CADDR(z), &id, &kode, F2C_CST(bip),  &ierr);
 *   DO_MTHERR("airye:");
 *   return 0;
 * } */

/**
 * Complex Modified Bessel function of the first kind
 *
 * @param z a complex number
 * @param v a real number, the order of the Bessel function
 *
 */
dcomplex pnl_complex_bessel_i ( double v, dcomplex z ) 
{
  int nz, ierr;
  dcomplex cy;
  int n = 1;
  int kode = 1;
  if (v >= 0)
    {
      pnl_zbesi(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
      DO_MTHERR("iv:");
    }
  else
    {
      dcomplex aux1, aux2;
      aux1 = pnl_complex_bessel_i (-v, z);
      aux2 = pnl_complex_bessel_k (-v, z);
      aux2 = CRmul (aux2, M_2_PI * sin(M_PI * -v));
      cy = Cadd (aux1, aux2);
    }
  return cy;
}

/**
 * Complex Modified Bessel function of the first kind
 * divided by exp(|Creal(z)|)
 *
 * @param z a complex number
 * @param v a real number, the order of the Bessel function
 *
 */
dcomplex pnl_complex_bessel_i_scaled( double v, dcomplex z ) 
{
  int nz, ierr;
  dcomplex cy;
  int n = 1;
  int kode = 2;
  if (v >= 0)
    {
      pnl_zbesi(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
      DO_MTHERR("ive:");
    }
  else
    {
      dcomplex aux1, aux2;
      aux1 = pnl_complex_bessel_i_scaled (-v, z);
      aux2 = pnl_complex_bessel_k (-v, z);
      aux2 = CRmul (aux2, M_2_PI * sin(M_PI * -v) * exp ( -fabs ( Creal(z) ) ));
      cy = Cadd (aux1, aux2);
    }
  return cy;
}


/**
 * Complex  Bessel function of the first kind
 *
 * @param z a complex number
 * @param v a real number, the order of the Bessel function
 *
 */
dcomplex pnl_complex_bessel_j( double v, dcomplex z ) 
{
  int nz, ierr;
  dcomplex cy;
  int n = 1;
  int kode = 1;
  if (v >= 0)
    {
      pnl_zbesj(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
      DO_MTHERR("jv:");
    }
  else
    {
      dcomplex aux1, aux2;
      aux1 = pnl_complex_bessel_j (-v, z);
      aux1 = CRmul (aux1, cos(M_PI * v));
      aux2 = pnl_complex_bessel_y (-v, z);
      aux2 = CRmul (aux2, sin(M_PI * v));
      cy = Cadd (aux1, aux2);
    }
  return cy;
}

/**
 * Complex  Bessel function of the first kind
 * divided by  exp(|Imag(z)|)
 *
 * @param z a complex number
 * @param v a real number, the order of the Bessel function
 *
 */
dcomplex pnl_complex_bessel_j_scaled ( double v, dcomplex z ) 
{
  int nz, ierr;
  dcomplex cy;
  int n = 1;
  int kode = 2;
  if (v >= 0)
    {
      pnl_zbesj(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
      DO_MTHERR("jve:");
    }
  else
    {
      dcomplex aux1, aux2;
      aux1 = pnl_complex_bessel_j_scaled (-v, z);
      aux1 = CRmul (aux1, cos(M_PI * v));
      aux2 = pnl_complex_bessel_y_scaled (-v, z);
      aux2 = CRmul (aux2, sin(M_PI * v));
      cy = Cadd (aux1, aux2);
    }
  return cy;
}

/**
 * Complex  Bessel function of the second kind
 *
 * @param z a complex number
 * @param v a real number, the order of the Bessel function
 *
 */  
dcomplex pnl_complex_bessel_y( double v, dcomplex z ) 
{
  int nz, ierr;
  dcomplex cy, cwork;
  int n = 1;
  int kode = 1;
  if (v >= 0)
    {
      pnl_zbesy(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, CADDR(cwork), &ierr);
      DO_MTHERR("yv:");
    }
  else
    {
      dcomplex aux1, aux2;
      aux1 = pnl_complex_bessel_y (-v, z);
      aux1 = CRmul (aux1, cos(M_PI * v));
      aux2 = pnl_complex_bessel_j (-v, z);
      aux2 = CRmul (aux2, sin(M_PI * -v));
      cy = Cadd (aux1, aux2);
    }
  return cy;
}

/**
 * Complex  Bessel function of the second kind
 * divided by  exp(|Imag(z)|)
 *
 * @param z a complex number
 * @param v a real number, the order of the Bessel function
 *
 */
dcomplex pnl_complex_bessel_y_scaled( double v, dcomplex z ) 
{
  int nz, ierr;
  dcomplex cy, cwork;
  int n = 1;
  int kode = 2;
  if (v >= 0)
    {
      pnl_zbesy(CADDR(z), &v, &kode, &n, CADDR(cy), &nz, CADDR(cwork), &ierr);
      DO_MTHERR("yve:");
    }
  else
    {
      dcomplex aux1, aux2;
      aux1 = pnl_complex_bessel_y_scaled (-v, z);
      aux1 = CRmul (aux1, cos(M_PI * v));
      aux2 = pnl_complex_bessel_j_scaled (-v, z);
      aux2 = CRmul (aux2, sin(M_PI * -v));
      cy = Cadd (aux1, aux2);
    }
  return cy;
}


/**
 * Complex  Modified Bessel function of the third kind
 *
 * @param z a complex number
 * @param v a real number, the order of the Bessel function
 *
 */
dcomplex pnl_complex_bessel_k( double v, dcomplex z ) 
{
  int nz, ierr;
  dcomplex cy;
  int n = 1;
  int kode = 1;
  double nu = fabs (v);
  pnl_zbesk(CADDR(z), &nu,  &kode, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("kv:");
  return cy;
}

/**
 * Complex Modified Bessel function of the third kind
 * multiplied by exp(z)
 *
 * @param z a complex number
 * @param v a real number, the order of the Bessel function
 *
 */
dcomplex pnl_complex_bessel_k_scaled( double v, dcomplex z ) 
{
  int nz, ierr;
  dcomplex cy;
  int n = 1;
  int kode = 2;
  double nu = fabs (v);
  pnl_zbesk(CADDR(z), &nu, &kode, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("kve:");
  return cy;
}
  
/**
 * Complex   Hankel function of the first kind
 *
 * @param z a complex number
 * @param v a real number, the order of the Bessel function
 *
 */
dcomplex pnl_complex_bessel_h1( double v, dcomplex z ) 
{
  int nz, ierr;
  dcomplex cy;
  int n = 1, kode = 1, m = 1;
  if (v >= 0)
    {
      pnl_zbesh(CADDR(z), &v,  &kode, &m, &n, CADDR(cy), &nz, &ierr);
      DO_MTHERR("hankel1:");
    }
  else
    {
      cy = pnl_complex_bessel_h1 (-v, z );
      cy = Cmul (cy, CIexp (M_PI * -v) );
    }
  return cy;
}

/**
 * Complex Hankel function of the first kind
 * divided by Cexp( I * z)
 *
 * @param z a complex number
 * @param v a real number, the order of the Bessel function
 *
 */
dcomplex pnl_complex_bessel_h1_scaled( double v, dcomplex z ) 
{
  int nz, ierr;
  dcomplex cy;
  int n = 1, kode = 2, m = 1;
  if (v >= 0)
    {
      pnl_zbesh(CADDR(z), &v, &kode, &m, &n, CADDR(cy), &nz, &ierr);
      DO_MTHERR("hankel1e:");
    }
  else
    {
      cy = pnl_complex_bessel_h1_scaled (-v, z );
      cy = Cmul (cy, CIexp (M_PI * -v) );
    }
  return cy;
}
  
/**
 * Complex Hankel function of the first kind
 *
 * @param z a complex number
 * @param v a real number, the order of the Bessel function
 *
 */
dcomplex pnl_complex_bessel_h2( double v, dcomplex z ) 
{
  int nz, ierr;
  dcomplex cy;
  int n = 1, kode = 1, m = 2;
  if (v >= 0)
    {
      pnl_zbesh(CADDR(z), &v,  &kode, &m, &n, CADDR(cy), &nz, &ierr);
      DO_MTHERR("hankel2:");
    }
  else
    {
      cy = pnl_complex_bessel_h2 (-v, z );
      cy = Cmul (cy, CIexp (M_PI * v) );
    }
  return cy;
}

/**
 * Complex Hankel function of the second kind
 * multiplied by Cexp( I * z)
 *
 * @param z a complex number
 * @param v a real number, the order of the Bessel function
 *
 */
dcomplex pnl_complex_bessel_h2_scaled( double v, dcomplex z ) 
{
  int nz, ierr;
  dcomplex cy;
  int n = 1, kode = 2, m = 2;
  if (v >= 0)
    {
      pnl_zbesh(CADDR(z), &v, &kode, &m, &n, CADDR(cy), &nz, &ierr);
      DO_MTHERR("hankel2e:");
    }
  else
    {
      cy = pnl_complex_bessel_h2_scaled (-v, z );
      cy = Cmul (cy, CIexp (M_PI * v) );
    }
  return cy;
}

/**
 * Compute the ratio of modified Bessel functions of the first kind
 * I_{v+1} / I_v
 * 
 * @param v a real number, the order of the Bessel function
 * @param x a complex number
 * 
 * @return I_{v+1}(x) / I_v(x)
 */
dcomplex pnl_complex_bessel_rati (double v, dcomplex x)
{
  int n;
  double d__1, tol;
  dcomplex cy;
  n = 1;
  d__1 = pnl_d1mach (4);
  tol = MAX (d__1, 1e-18);

  pnl_zrati (CADDR(x), &v, &n, CADDR(cy), &tol);
  return cy;
}


/* real Bessel functions */
/**
 *  Modified Bessel function of the first kind
 *
 * @param x a real number
 * @param v a real number, the order of the Bessel function
 *
 */
double pnl_bessel_i ( double v, double x ) 
{
  int nz, ierr;
  dcomplex cy, z = Complex (x, 0.);
  int n = 1, kode = 1;
  if (v >= 0)
    {
      pnl_zbesi(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
      DO_MTHERR("iv:");
      return cy.r;
    }
  else
    {
      double  aux1, aux2;
      aux1 = pnl_bessel_i (-v, x);
      aux2 = M_2_PI * sin (M_PI * -v) * pnl_bessel_k (-v, x);
      return aux1 + aux2;
    }
}

/**
 * Modified Bessel function of the first kind
 * divided by exp(|x|)
 *
 * @param x a real number
 * @param v a real number, the order of the Bessel function
 *
 */
double pnl_bessel_i_scaled( double v, double x ) 
{
  int nz, ierr;
  dcomplex cy, z = Complex (x, 0.);
  int n = 1, kode = 2;
  if (v >= 0)
    {
      pnl_zbesi(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
      DO_MTHERR("ive:");
      return cy.r;
    }
  else
    {
      double  aux1, aux2;
      aux1 = pnl_bessel_i_scaled (-v, x);
      aux2 = M_2_PI * sin (M_PI * -v) * pnl_bessel_k (-v, x) * exp(-fabs(x));
      return aux1 + aux2;
    }
}


/**
 *   Bessel function of the first kind
 *
 * @param x a real number
 * @param v a real number, the order of the Bessel function
 *
 */
double pnl_bessel_j( double v, double x ) 
{
  int nz, ierr;
  dcomplex cy, z = Complex (x, 0.);
  int n = 1, kode = 1;
  if (v >= 0)
    {
      pnl_zbesj(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
      DO_MTHERR("jv:");
      return cy.r;
    }
  else
    {
      double  aux1, aux2;
      aux1 = cos (M_PI * v) * pnl_bessel_j (-v, x);
      aux2 = sin (M_PI * v) * pnl_bessel_y (-v, x);
      return aux1 + aux2;
    }
}

/**
 * Bessel function of the first kind
 * same as pnl_bessel_j
 *
 * @param x a real number
 * @param v a real number, the order of the Bessel function
 *
 */
double pnl_bessel_j_scaled ( double v, double x ) 
{
  int nz, ierr;
  dcomplex cy, z = Complex (x, 0.);
  int n = 1, kode = 2;
  if (v >= 0)
    {
      pnl_zbesj(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
      DO_MTHERR("jv:");
      return cy.r;
    }
  else
    {
      double  aux1, aux2;
      aux1 = cos (M_PI * v) * pnl_bessel_j_scaled (-v, x);
      aux2 = sin (M_PI * v) * pnl_bessel_y_scaled (-v, x);
      return aux1 +aux2;
    }
}

/**
 *   Bessel function of the second kind
 *
 * @param x a real number
 * @param v a real number, the order of the Bessel function
 *
 */  
double pnl_bessel_y( double v, double x ) 
{
  int nz, ierr;
  dcomplex cy, z, cwork;
  int n = 1, kode = 1;
  z = Complex (x, 0.);
  if (v >= 0)
    {
      pnl_zbesy(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, CADDR(cwork), &ierr);
      DO_MTHERR("yv:");
      return cy.r;
    }
  else
    {
      double  aux1, aux2;
      aux1 = cos (M_PI * v) * pnl_bessel_y (-v, x);
      aux2 = sin (M_PI * -v) * pnl_bessel_j (-v, x);
      return aux1 +aux2;
    }
}

/**
 * Scaled   Bessel function of the second kind
 * same as pnl_bessel_y
 *
 * @param x a real number
 * @param v a real number, the order of the Bessel function
 *
 */
double pnl_bessel_y_scaled( double v, double x ) 
{
  int nz, ierr;
  dcomplex cy, z, cwork;
  int n = 1;
  int kode = 2;
  z = Complex (x, 0.);
  if (v >= 0)
    {
      pnl_zbesy(CADDR(z), &v, &kode, &n, CADDR(cy), &nz, CADDR(cwork), &ierr);
      DO_MTHERR("yve:");
      return cy.r;
    }
  else
    {
      double  aux1, aux2;
      aux1 = cos (M_PI * v) * pnl_bessel_y_scaled (-v, x);
      aux2 = sin (M_PI * -v) * pnl_bessel_j_scaled (-v, x);
      return aux1 +aux2;
    }
}


/**
 *   Modified Bessel function of the third kind
 *
 * @param x a real number
 * @param v a real number, the order of the Bessel function
 *
 */
double pnl_bessel_k( double v, double x ) 
{
  int nz, ierr;
  dcomplex cy, z = Complex (x, 0.);
  int n = 1;
  int kode = 1;
  double nu = fabs(v);
  pnl_zbesk(CADDR(z), &nu,  &kode, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("kv:");
  return cy.r;
}

/**
 * Modified Bessel function of the third kind
 * multiplied by exp(x)
 *
 * @param x a real number
 * @param v a real number, the order of the Bessel function
 *
 */
double pnl_bessel_k_scaled( double v, double x ) 
{
  int nz, ierr;
  dcomplex cy, z = Complex (x, 0.);
  int n = 1;
  int kode = 2;
  double nu = fabs (v);
  pnl_zbesk(CADDR(z), &nu, &kode, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("kve:");
  return cy.r;
}
  
/**
 * Hankel function of the first kind
 *
 * @param x a real number
 * @param v a real number, the order of the Bessel function
 *
 */
dcomplex pnl_bessel_h1( double v, double x ) 
{
  int nz, ierr;
  dcomplex cy, z = Complex (x, 0.);
  int n = 1;
  int kode = 1;
  int m = 1;
  if (v >= 0)
    {
      pnl_zbesh(CADDR(z), &v,  &kode, &m, &n, CADDR(cy), &nz, &ierr);
      DO_MTHERR("hankel1:");
    }
  else
    {
      cy = pnl_bessel_h1 (-v, x );
      cy = Cmul (cy, CIexp (M_PI * -v) );
    }
  return cy;
}

/**
 * Hankel function of the first kind
 * divided by Cexp( I * x)
 *
 * @param x a real number
 * @param v a real number, the order of the Bessel function
 *
 */
dcomplex pnl_bessel_h1_scaled( double v, double x ) 
{
  int nz, ierr;
  dcomplex cy, z = Complex (x, 0.);
  int n = 1;
  int kode = 2;
  int m = 1;
  if (v >= 0)
    {
      pnl_zbesh(CADDR(z), &v,  &kode, &m, &n, CADDR(cy), &nz, &ierr);
      DO_MTHERR("hankel1:");
    }
  else
    {
      cy = pnl_bessel_h1_scaled (-v, x );
      cy = Cmul (cy, CIexp (M_PI * -v) );
    }
  return cy;
}
  
/**
 *  Hankel function of the first kind
 *
 * @param x a real number
 * @param v a real number, the order of the Bessel function
 *
 */
dcomplex pnl_bessel_h2( double v, double x ) 
{
  int nz, ierr;
  dcomplex cy, z = Complex (x, 0.);
  int n = 1;
  int kode = 1;
  int m = 2;
  if (v >= 0)
    {
      pnl_zbesh(CADDR(z), &v,  &kode, &m, &n, CADDR(cy), &nz, &ierr);
      DO_MTHERR("hankel2:");
    }
  else
    {
      cy = pnl_bessel_h2 (-v, x );
      cy = Cmul (cy, CIexp (M_PI * v) );
    }
  return cy;
}

/**
 * Hankel function of the second kind
 * multiplied by Cexp( I * x)
 *
 * @param x a real number
 * @param v a real number, the order of the Bessel function
 *
 */
dcomplex pnl_bessel_h2_scaled( double v, double x ) 
{
  int nz, ierr;
  dcomplex cy, z = Complex (x, 0.);
  int n = 1;
  int kode = 2;
  int m = 2;
  if (v >= 0)
    {
      pnl_zbesh(CADDR(z), &v,  &kode, &m, &n, CADDR(cy), &nz, &ierr);
      DO_MTHERR("hankel2:");
    }
  else
    {
      cy = pnl_bessel_h2_scaled (-v, x );
      cy = Cmul (cy, CIexp (M_PI * v) );
    }
  return cy;
}

/**
 * Compute the ratio of modified Bessel functions of the first kind
 * I_{v+1} / I_v
 * 
 * @param v a real number, the order of the Bessel function
 * @param x a real number
 * 
 * @return I_{v+1}(x) / I_v(x)
 */
double pnl_bessel_rati (double v, double x)
{
  double d__1, tol;
  int  n;
  dcomplex cy, z;
  n = 1;
  z = Complex (x, 0.);
  d__1 = pnl_d1mach (4);
  tol = MAX (d__1, 1e-18);

  pnl_zrati (CADDR(z), &v, &n, CADDR(cy), &tol);
  return cy.r;
}
