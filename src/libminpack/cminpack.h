#ifndef __CMINPACK_H__
#define __CMINPACK_H__

#include "pnl/pnl_mathtools.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* Declarations for minpack */

/* Function types: */
/* The first argument can be used to store extra function parameters, thus */
/* avoiding the use of global variables. */
/* the iflag parameter is input-only (with respect to the FORTRAN */
/*  version), the output iflag value is the return value of the function. */
/* If iflag=0, the function shoulkd just print the current values (see */
/* the nprint parameters below). */
  
/* for hybrd1 and hybrd: */
/*         calculate the functions at x and */
/*         return this vector in fvec. */
/* return a negative value to terminate hybrd1/hybrd */
typedef int (*minpack_func_nn)(void *p, int n, const double *x, double
                               *fvec, int iflag );

/* for hybrj1 and hybrj */
/*         if iflag = 1 calculate the functions at x and */
/*         return this vector in fvec. do not alter fjac. */
/*         if iflag = 2 calculate the jacobian at x and */
/*         return this matrix in fjac. do not alter fvec. */
/* return a negative value to terminate hybrj1/hybrj */
typedef int (*minpack_funcder_nn)(void *p, int n, const double *x, double
                                  *fvec, double *fjac, int ldfjac, int
                                  iflag );

/* for lmdif1 and lmdif */
/*         calculate the functions at x and */
/*         return this vector in fvec. */
/* return a negative value to terminate lmdif1/lmdif */
typedef int (*minpack_func_mn)(void *p, int m, int n, const double *x,
                               double *fvec, int iflag );

/* for lmder1 and lmder */
/*         if iflag = 1 calculate the functions at x and */
/*         return this vector in fvec. do not alter fjac. */
/*         if iflag = 2 calculate the jacobian at x and */
/*         return this matrix in fjac. do not alter fvec. */
/* return a negative value to terminate lmder1/lmder */
typedef int (*minpack_funcder_mn)(void *p, int m, int n, const double *x,
                                  double *fvec, double *fjac, int ldfjac,
                                  int iflag );

/* for lmstr1 and lmstr */
/*         if iflag = 1 calculate the functions at x and */
/*         return this vector in fvec. */
/*         if iflag = i calculate the (i-1)-st row of the */
/*         jacobian at x and return this vector in fjrow. */
/* return a negative value to terminate lmstr1/lmstr */
typedef int (*minpack_funcderstr_mn)(void *p, int m, int n, const double
                                     *x, double *fvec, double *fjrow, int
                                     iflag );


/* MINPACK functions: */
/* the info parameter was removed from most functions: the return */
/* value of the function is used instead. */
/* The argument 'p' can be used to store extra function parameters, thus */
/* avoiding the use of global variables. You can also think of it as a */
/* 'this' pointer a la C++. */

/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (Jacobian calculated by
   a forward-difference approximation, more general). */
int pnl_minpack_hybrd ( minpack_func_nn fcn, void *p, int n, double *x,
                        double *fvec, double xtol, int maxfev, int ml, int
                        mu, double epsfcn, double *diag, int mode, double
                        factor, int nprint, int *nfev, double *fjac, int
                        ldfjac, double *r, int lr, double *qtf, double
                        *wa1, double *wa2, double *wa3, double *wa4);
          
/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (user-supplied Jacobian,
   more general) */
int pnl_minpack_hybrj ( minpack_funcder_nn fcn, void *p, int n, double *x,
                        double *fvec, double *fjac, int ldfjac, double
                        xtol, int maxfev, double *diag, int mode, double
                        factor, int nprint, int *nfev, int *njev, double
                        *r, int lr, double *qtf, double *wa1, double *wa2,
                        double *wa3, double *wa4 );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (Jacobian calculated by a forward-difference approximation, more
   general) */
int pnl_minpack_lmdif ( minpack_func_mn fcn, void *p, int m, int n, double
                        *x, double *fvec, double ftol, double xtol, double
                        gtol, int maxfev, double epsfcn, double *diag, int
                        mode, double factor, int nprint, int *nfev, double
                        *fjac, int ldfjac, int *ipvt, double *qtf, double
                        *wa1, double *wa2, double *wa3, double *wa4 );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian, more general) */
int pnl_minpack_lmder ( minpack_funcder_mn fcn, void *p, int m, int n,
                        double *x, double *fvec, double *fjac, int ldfjac,
                        double ftol, double xtol, double gtol, int maxfev,
                        double *diag, int mode, double factor, int nprint,
                        int *nfev, int *njev, int *ipvt, double *qtf,
                        double *wa1, double *wa2, double *wa3, double *wa4
                        );

void pnl_minpack_chkder ( int m, int n, const double *x, double *fvec,
                          double *fjac, int ldfjac, double *xp, double
                          *fvecp, int mode, double *err  );

double pnl_minpack_dpmpar ( int i );

double pnl_minpack_enorm ( int n, const double *x );

/* compute a forward-difference approximation to the m by n jacobian
   matrix associated with a specified problem of m functions in n
   variables. */
int pnl_minpack_fdjac2(minpack_func_mn fcn, void *p, int m, int n, double
                       *x, const double *fvec, double *fjac, int ldfjac,
                       double epsfcn, double *wa);

/* compute a forward-difference approximation to the n by n jacobian
   matrix associated with a specified problem of n functions in n
   variables. if the jacobian has a banded form, then function
   evaluations are saved by only approximating the nonzero terms. */
int pnl_minpack_fdjac1(minpack_func_nn fcn, void *p, int n, double *x,
                       const double *fvec, double *fjac, int ldfjac, int
                       ml, int mu, double epsfcn, double *wa1, double
                       *wa2);

/* internal MINPACK subroutines */
void pnl_minpack_dogleg(int n, const double *r__, int lr, const double
                        *diag, const double *qtb, double delta, double *x,
                        double *wa1, double *wa2);
void pnl_minpack_qrfac(int m, int n, double *a, int lda, int pivot, int
                       *ipvt, int lipvt, double *rdiag, double *acnorm,
                       double *wa);
void pnl_minpack_qrsolv(int n, double *r__, int ldr, const int *ipvt, const
                        double *diag, const double *qtb, double *x, double
                        *sdiag, double *wa);
void pnl_minpack_qform(int m, int n, double *q, int
                       ldq, double *wa);
void pnl_minpack_r1updt(int m, int n, double *s, int ls, const double *u,
                        double *v, double *w, int *sing);
void pnl_minpack_r1mpyq(int m, int n, double *a, int lda, const double *v,
                        const double *w);
void pnl_minpack_lmpar(int n, double *r__, int ldr, const int *ipvt, const
                       double *diag, const double *qtb, double delta,
                       double *par, double *x, double *sdiag, double *wa1,
                       double *wa2);
void pnl_minpack_rwupdt(int n, double *r__, int ldr, const double *w,
                        double *b, double *alpha, double *cos__, double
                        *sin__);
void pnl_minpack_covar(int n, double *r__, int ldr, const int *ipvt, double
                       tol, double *wa);
#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* __CMINPACK_H__ */
