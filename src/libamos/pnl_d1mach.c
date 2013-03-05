#include <stdio.h>
#include "amos.h" 

extern double pnl_dlamch (char *cmach);

/* 
 * emulation of the needed values for d1mach and i1mach 
 * using nsp functions (jpc Dec 2007)
 */

/*
 *  D1MACH can be used to obtain machine-dependent parameters 
 *  for the local machine environment.  It is a function 
 *  subprogram with one (input) argument, and can be called 
 *  as follows, for example 
 * 
 *       D = D1MACH(I) 
 * 
 *  where I=1,...,5.  The (output) value of D above is 
 *  determined by the (input) value of I.  The results for 
 *  various values of I are discussed below. 
 * 
 *  D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude. 
 *  D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude. 
 *  D1MACH( 3) = B**(-T), the smallest relative spacing. 
 *  D1MACH( 4) = B**(1-T), the largest relative spacing. 
 *  D1MACH( 5) = LOG10(B) 
 * 
 *  Assume double precision numbers are represented in the T-digit, 
 *  base-B form 
 * 
 *             sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) 
 * 
 *  where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and 
 *  EMIN .LE. E .LE. EMAX. 
 * 
 *  The values of B, T, EMIN and EMAX are provided in I1MACH as 
 *  follows: 
 *  I1MACH(10) = B, the base. 
 *  I1MACH(14) = T, the number of base-B digits. 
 *  I1MACH(15) = EMIN, the smallest exponent E. 
 *  I1MACH(16) = EMAX, the largest exponent E. 
 */

double pnl_d1mach (int i)
{
  switch (i) 
    {
    case 1: return pnl_dlamch("u");
    case 2: return pnl_dlamch("o");
    case 3: return pnl_dlamch("e");
    case 4: return pnl_dlamch("p");
    case 5: return log(pnl_dlamch("b"))/log(10.0);
    default :
            printf("d1mach called with wrong argument\n");
            return 0;
    }
}

/*
 *    I1MACH can be used to obtain machine-dependent parameters 
 *    for the local machine environment.  It is a function 
 *    subroutine with one (input) argument, and can be called 
 *    as follows, for example 
 * 
 *         K = I1MACH(I) 
 * 
 *    where I=1,...,16.  The (output) value of K above is 
 *    determined by the (input) value of I.  The results for 
 *    various values of I are discussed below. 
 * 
 * I/O unit numbers. 
 *   I1MACH( 1) = the standard input unit. 
 *   I1MACH( 2) = the standard output unit. 
 *   I1MACH( 3) = the standard punch unit. 
 *   I1MACH( 4) = the standard error message unit. 
 * 
 * Words. 
 *   I1MACH( 5) = the number of bits per int storage unit. 
 *   I1MACH( 6) = the number of characters per int storage unit. 
 * 
 * Ints. 
 *   assume ints are represented in the S-digit, base-A form 
 * 
 *              sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) ) 
 * 
 *              where 0 .LE. X(I) .LT. A for I=0,...,S-1. 
 *   I1MACH( 7) = A, the base. 
 *   I1MACH( 8) = S, the number of base-A digits. 
 *   I1MACH( 9) = A**S - 1, the largest magnitude. 
 * 
 * Floating-Point Numbers. 
 *   Assume floating-point numbers are represented in the T-digit, 
 *   base-B form 
 *              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) 
 * 
 *              where 0 .LE. X(I) .LT. B for I=1,...,T, 
 *              0 .LT. X(1), and EMIN .LE. E .LE. EMAX. 
 *   I1MACH(10) = B, the base. 
 * 
 * Single-Precision 
 *   I1MACH(11) = T, the number of base-B digits. 
 *   I1MACH(12) = EMIN, the smallest exponent E. 
 *   I1MACH(13) = EMAX, the largest exponent E. 
 * 
 * Double-Precision 
 *   I1MACH(14) = T, the number of base-B digits. 
 *   I1MACH(15) = EMIN, the smallest exponent E. 
 *   I1MACH(16) = EMAX, the largest exponent E. 
 * 
 */

static int largestint(void)
{
  static int first=0,large;
  if ( first == 0) 
    {
      if (sizeof(int)==sizeof(long)) large = (int) LONG_MAX ;
      else if (sizeof(int)==sizeof(short)) large = SHRT_MAX;
      else large = 2147483647 ; /** using default value **/
      first++;
    }
  return large;
}


int amos_i1mach (int i)
{
  switch (i) 
    {
    case 9: return largestint(); 
    case 14: return pnl_dlamch("n");
    case 15: return pnl_dlamch("m");
    case 16: return pnl_dlamch("l");
    default :
             printf("1imach called with wrong argument\n");
             return 0;
    }
  return 0;
}

int pnl_ipmpar(int i)
{
  switch (i) 
    {
    case 3: return largestint();
    case 4: return (int) pnl_dlamch("b");
    case 9: return (int) pnl_dlamch("m");
    case 10: return (int) pnl_dlamch("l");
    default :
             printf("ipmpar called with wrong argument\n");
             return 0;
    }
  return 0;
}

