#ifndef _PNL_INTEGRATION_H
#define _PNL_INTEGRATION_H


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl/pnl_mathtools.h"

extern double pnl_integration (PnlFunc *F, double x0, double x1, 
                  int n, char *meth);

extern double pnl_integration_2d (PnlFunc2D *F, double x0, double x1,
                  double y0,double y1, int nx, int ny, char *meth);

extern int pnl_integration_qng (PnlFunc *F, double x0, double x1,
                  double epsabs, double epsrel, double * result,
                  double * abserr, int * neval);

extern int pnl_integration_GK (PnlFunc *F, double x0, double x1,
                  double epsabs, double epsrel, double * result,
                  double * abserr, int * neval);

extern int pnl_integration_qng_2d (PnlFunc2D *F, double x0, double x1,
                  double y0,double y1, double epsabs, double epsrel,
                  double * result, double * abserr, int * neval);

extern int pnl_integration_GK2D (PnlFunc2D *F, double x0, double x1,
                  double y0,double y1, double epsabs, double epsrel,
                  double * result, double * abserr, int * neval);

extern int pnl_integration_qag (PnlFunc *f, double a, double b,
                   double epsabs, double epsrel, int limit, double *result, 
                   double *abserr, int *neval);

extern int pnl_integration_qagp (PnlFunc *f, double a, double b, const PnlVect *singularities,
                   double epsabs, double epsrel, int limit, double *result, 
                   double *abserr, int *neval);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _PNL_INTEGRATION_H */
