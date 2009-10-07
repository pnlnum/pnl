/*************************************************************************/
/* Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>            */
/*                                                                       */
/* This program is free software; you can redistribute it and/or modify  */
/* it under the terms of the GNU General Public License as published by  */
/* the Free Software Foundation; either version 3 of the License, or     */
/* (at your option) any later version.                                   */
/*                                                                       */
/* This program is distributed in the hope that it will be useful,       */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/* GNU General Public License for more details.                          */
/*                                                                       */
/* You should have received a copy of the GNU General Public License     */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/*************************************************************************/

#ifndef PNL_FFT_H
#define PNL_FFT_H


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl_mathtools.h"
#include "pnl_complex.h"
#include "pnl_vector.h"

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
/*@}*/
#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* PNL_FFT_H */
