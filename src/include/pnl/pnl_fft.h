#ifndef _PNL_FFT_H
#define _PNL_FFT_H

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_complex.h"
#include "pnl/pnl_vector.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**
 * \defgroup FFT Discrete Fourier Transform
 */

/*@{*/

extern int pnl_fft_inplace (PnlVectComplex * data);
extern int pnl_ifft_inplace (PnlVectComplex * data);
extern int pnl_fft (const PnlVectComplex * in, PnlVectComplex * out);
extern int pnl_ifft (const PnlVectComplex * in, PnlVectComplex * out);
extern int pnl_fft2(double *re, double *im, int n);
extern int pnl_ifft2(double *re, double *im, int n);

extern int pnl_real_fft_inplace(double *data, int n);
extern int pnl_real_ifft_inplace(double *data, int n);
extern int pnl_real_fft(const PnlVect *in, PnlVectComplex *out);
extern int pnl_real_ifft(const PnlVectComplex *in, PnlVect *out);
extern int pnl_real_fft2(double *re, double *im, int n);
extern int pnl_real_ifft2(double *re, double *im, int n);

extern int pnl_ifft2d_inplace (PnlMatComplex *data);
extern int pnl_fft2d_inplace (PnlMatComplex *data);
extern int pnl_real_fft2d(const PnlMat *in, PnlMatComplex *out);
extern int pnl_real_ifft2d(const PnlMatComplex *in, PnlMat *out);

/*@}*/
#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_FFT_H */
