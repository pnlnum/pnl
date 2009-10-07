#ifndef _PNL_INTEGRATIO
#define _PNL_INTEGRATION


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl_mathtools.h"

/* Gauss-Kronrod-Patterson quadrature coefficients for use in
   quadpack routine qng. These coefficients were calculated with
   101 decimal digit arithmetic by L. W. Fullerton, Bell Labs, Nov
   1981. */

extern int pnl_integration_GK (const PnlFunc *F,
                               double x0, double x1,
                               double epsabs, double epsrel,
                               double * result, double * abserr,
                               int * neval);

extern int pnl_integration_GK2D (const PnlFunc2D *F,
                                 double x0, double x1,
                                 double y0,double y1,
                                 double epsabs, double epsrel,
                                 double * result, double * abserr,
                                 int * neval);


extern double pnl_integration (const PnlFunc *F, double x0, double x1, int n, char *meth);
extern double pnl_integration_2D (const PnlFunc2D *F, double x0, double x1,
                                  double y0,double y1, int nx, int ny, char *meth);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_INTEGRATION */
