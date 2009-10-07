
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

#ifndef PNL_LAPLACE
#define PNL_LAPLACE


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl_mathtools.h"
#include "pnl_complex.h"
#include "pnl_vector.h"

extern double pnl_ilap_cdf_euler(PnlCmplxFunc *f, double t, double h, int N, int M);
extern double pnl_ilap_euler(PnlCmplxFunc *f, double t, int N, int M);
extern void pnl_ilap_fft(PnlVect *res, PnlCmplxFunc *f, double T, double eps);
extern double pnl_ilap_gs_basic (PnlFunc *fhat, double t, int n);
extern double pnl_ilap_gs (PnlFunc *fhat, double t, int n);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* PNL_LAPLACE */
