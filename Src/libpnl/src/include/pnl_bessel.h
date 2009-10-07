#ifndef _PNL_BESSEL
#define _PNL_BESSEL

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


#include "pnl_complex.h"

/**
 * \defgroup PnlBessel  Bessel functions
 * \brief Approximations of the Bessel functions
 */

/*@{*/

/* complex Bessel functions from libamos */

extern fcomplex pnl_complex_bessel_i ( double v, fcomplex z );
extern fcomplex pnl_complex_bessel_i_scaled( double v, fcomplex z );
extern fcomplex pnl_complex_bessel_j( double v, fcomplex z );
extern fcomplex pnl_complex_bessel_j_scaled ( double v, fcomplex z );
extern fcomplex pnl_complex_bessel_y( double v, fcomplex z );
extern fcomplex pnl_complex_bessel_y_scaled( double v, fcomplex z );
extern fcomplex pnl_complex_bessel_k( double v, fcomplex z );
extern fcomplex pnl_complex_bessel_k_scaled( double v, fcomplex z );
extern fcomplex pnl_complex_bessel_h1( double v, fcomplex z );
extern fcomplex pnl_complex_bessel_h1_scaled( double v, fcomplex z );
extern fcomplex pnl_complex_bessel_h2( double v, fcomplex z );
extern fcomplex pnl_complex_bessel_h2_scaled( double v, fcomplex z );


/* real Bessel functions from libamos */

extern double   pnl_bessel_i ( double v, double x );
extern double   pnl_bessel_i_scaled( double v, double x );
extern double   pnl_bessel_j( double v, double x );
extern double   pnl_bessel_j_scaled ( double v, double x );
extern double   pnl_bessel_y( double v, double x );
extern double   pnl_bessel_y_scaled( double v, double x );
extern double   pnl_bessel_k( double v, double x );
extern double   pnl_bessel_k_scaled( double v, double x );
extern fcomplex pnl_bessel_h1( double v, double x );
extern fcomplex pnl_bessel_h1_scaled( double v, double x );
extern fcomplex pnl_bessel_h2( double v, double x );
extern fcomplex pnl_bessel_h2_scaled( double v, double x );

/* @} */
#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif
