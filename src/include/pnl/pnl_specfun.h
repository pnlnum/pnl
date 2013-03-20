#ifndef _PNL_SPECFUN_H 
#define _PNL_SPECFUN_H 

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "pnl/pnl_complex.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**
 * \defgroup PnlSpecFun  Special Functions
 * \brief Approximations of some special functions.
 * We provide approximations for incomplete Gamma functions, error functions,
 * exponential integrals and hypergeomtric functions
 */
/*@{*/

extern void pnl_deactivate_mtherr (void);
extern void pnl_activate_mtherr (void);


/**
 * \defgroup PnlBessel  Bessel functions
 * \brief Approximations of the Bessel functions
 */

/**
 * \ingroup PnlBessel
 */
/*@{*/

/* complex Bessel functions from libamos */

extern dcomplex pnl_complex_bessel_i ( double v, dcomplex z );
extern dcomplex pnl_complex_bessel_i_scaled( double v, dcomplex z );
extern dcomplex pnl_complex_bessel_j( double v, dcomplex z );
extern dcomplex pnl_complex_bessel_j_scaled ( double v, dcomplex z );
extern dcomplex pnl_complex_bessel_y( double v, dcomplex z );
extern dcomplex pnl_complex_bessel_y_scaled( double v, dcomplex z );
extern dcomplex pnl_complex_bessel_k( double v, dcomplex z );
extern dcomplex pnl_complex_bessel_k_scaled( double v, dcomplex z );
extern dcomplex pnl_complex_bessel_h1( double v, dcomplex z );
extern dcomplex pnl_complex_bessel_h1_scaled( double v, dcomplex z );
extern dcomplex pnl_complex_bessel_h2( double v, dcomplex z );
extern dcomplex pnl_complex_bessel_h2_scaled( double v, dcomplex z );
extern dcomplex pnl_complex_bessel_rati (double v, dcomplex x);


/* real Bessel functions from libamos */

extern double   pnl_bessel_i ( double v, double x );
extern double   pnl_bessel_i_scaled( double v, double x );
extern double   pnl_bessel_j( double v, double x );
extern double   pnl_bessel_j_scaled ( double v, double x );
extern double   pnl_bessel_y( double v, double x );
extern double   pnl_bessel_y_scaled( double v, double x );
extern double   pnl_bessel_k( double v, double x );
extern double   pnl_bessel_k_scaled( double v, double x );
extern dcomplex pnl_bessel_h1( double v, double x );
extern dcomplex pnl_bessel_h1_scaled( double v, double x );
extern dcomplex pnl_bessel_h2( double v, double x );
extern dcomplex pnl_bessel_h2_scaled( double v, double x );
extern double pnl_bessel_rati (double v, double x);

/* @} */

extern double pnl_sf_gamma_inc(double a,double x);
extern double pnl_sf_gamma_inc_P(double a,double x);
extern double pnl_sf_gamma_inc_Q(double a,double x);
extern double pnl_sf_expint_En(int n,double x);
/* extern double pnl_sf_expint_Ei(double x); */


extern double pnl_sf_fact (int n);
extern double pnl_sf_choose (int n , int k);
extern double pnl_sf_gamma (double);
extern double pnl_sf_log_gamma (double);
extern int pnl_sf_log_gamma_sgn(double x, double *res, int *sgn);
extern double pnl_sf_log_erf (double x);
extern double pnl_sf_log_erfc (double x);
extern double pnl_sf_erf (double x);
extern double pnl_sf_erfc (double x);
extern double pnl_sf_psi (double x);

extern double pnl_sf_hyperg_2F1 (double a, double b, double c, double x);
extern double pnl_sf_hyperg_1F1 (double a, double b, double x);
extern double pnl_sf_hyperg_2F0 (double a, double b, double x);
extern double pnl_sf_hyperg_0F1 (double c, double x);
extern double pnl_sf_hyperg_U (double a, double b, double x);

/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_SPECFUN_H */ 
