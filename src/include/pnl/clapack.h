#ifndef _LAPACK_H_
#define _LAPACK_H_

/*
 * Header files for calling Lapack routines from C code
 * This file comes from clapack and has been slighty modified
 * for a better integration into Pnl
 */

#include "pnl/pnl_complex.h"
#include "pnl/pnl_internals.h"

#ifdef __cplusplus 	
extern "C" {	
#endif		


typedef int (*L_fp)(void);

extern int C2F(caxpy)(int *n, dcomplex *ca, dcomplex *cx, int * incx, dcomplex *cy, int *incy);

extern int C2F(ccopy)(int *n, dcomplex *cx, int *incx, dcomplex * cy, int *incy);

void C2F(cdotc)(dcomplex * ret_val, int *n, dcomplex *cx, int *incx, dcomplex *cy, int *incy);

void C2F(cdotu)(dcomplex * ret_val, int *n, dcomplex *cx, int *incx, dcomplex *cy, int *incy);

extern int C2F(cgbmv)(char *trans, int *m, int *n, int *kl, 
                      int *ku, dcomplex *alpha, dcomplex *a, int *lda, dcomplex *x, 
                      int *incx, dcomplex *beta, dcomplex *y, int *incy);

extern int C2F(cgemm)(char *transa, char *transb, int *m, int *
                      n, int *k, dcomplex *alpha, dcomplex *a, int *lda, dcomplex *b, 
                      int *ldb, dcomplex *beta, dcomplex *c__, int *ldc);

extern int C2F(cgemv)(char *trans, int *m, int *n, dcomplex *
                      alpha, dcomplex *a, int *lda, dcomplex *x, int *incx, dcomplex *
                      beta, dcomplex *y, int *incy);

extern int C2F(cgerc)(int *m, int *n, dcomplex *alpha, dcomplex *
                      x, int *incx, dcomplex *y, int *incy, dcomplex *a, int *lda);

extern int C2F(cgeru)(int *m, int *n, dcomplex *alpha, dcomplex *
                      x, int *incx, dcomplex *y, int *incy, dcomplex *a, int *lda);

extern int C2F(chbmv)(char *uplo, int *n, int *k, dcomplex *
                      alpha, dcomplex *a, int *lda, dcomplex *x, int *incx, dcomplex *
                      beta, dcomplex *y, int *incy);

extern int C2F(chemm)(char *side, char *uplo, int *m, int *n, 
                      dcomplex *alpha, dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                      dcomplex *beta, dcomplex *c__, int *ldc);

extern int C2F(chemv)(char *uplo, int *n, dcomplex *alpha, dcomplex *
                      a, int *lda, dcomplex *x, int *incx, dcomplex *beta, dcomplex *y, 
                      int *incy);

extern int C2F(cher)(char *uplo, int *n, double *alpha, dcomplex *x, 
                     int *incx, dcomplex *a, int *lda);

extern int C2F(cher2)(char *uplo, int *n, dcomplex *alpha, dcomplex *
                      x, int *incx, dcomplex *y, int *incy, dcomplex *a, int *lda);

extern int C2F(cher2k)(char *uplo, char *trans, int *n, int *k, 
                       dcomplex *alpha, dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       double *beta, dcomplex *c__, int *ldc);

extern int C2F(cherk)(char *uplo, char *trans, int *n, int *k, 
                      double *alpha, dcomplex *a, int *lda, double *beta, dcomplex *c__, 
                      int *ldc);

extern int C2F(chpmv)(char *uplo, int *n, dcomplex *alpha, dcomplex *
                      ap, dcomplex *x, int *incx, dcomplex *beta, dcomplex *y, int *
                      incy);

extern int C2F(chpr)(char *uplo, int *n, double *alpha, dcomplex *x, 
                     int *incx, dcomplex *ap);

extern int C2F(chpr2)(char *uplo, int *n, dcomplex *alpha, dcomplex *
                      x, int *incx, dcomplex *y, int *incy, dcomplex *ap);

extern int C2F(crotg)(dcomplex *ca, dcomplex *cb, double *c__, dcomplex *s);

extern int C2F(cscal)(int *n, dcomplex *ca, dcomplex *cx, int *
                      incx);

extern int C2F(csrot)(int *n, dcomplex *cx, int *incx, dcomplex *
                      cy, int *incy, double *c__, double *s);

extern int C2F(csscal)(int *n, double *sa, dcomplex *cx, int *incx);

extern int C2F(cswap)(int *n, dcomplex *cx, int *incx, dcomplex *
                      cy, int *incy);

extern int C2F(csymm)(char *side, char *uplo, int *m, int *n, 
                      dcomplex *alpha, dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                      dcomplex *beta, dcomplex *c__, int *ldc);

extern int C2F(csyr2k)(char *uplo, char *trans, int *n, int *k, 
                       dcomplex *alpha, dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       dcomplex *beta, dcomplex *c__, int *ldc);

extern int C2F(csyrk)(char *uplo, char *trans, int *n, int *k, 
                      dcomplex *alpha, dcomplex *a, int *lda, dcomplex *beta, dcomplex *c__, 
                      int *ldc);

extern int C2F(ctbmv)(char *uplo, char *trans, char *diag, int *n, 
                      int *k, dcomplex *a, int *lda, dcomplex *x, int *incx);

extern int C2F(ctbsv)(char *uplo, char *trans, char *diag, int *n, 
                      int *k, dcomplex *a, int *lda, dcomplex *x, int *incx);

extern int C2F(ctpmv)(char *uplo, char *trans, char *diag, int *n, 
                      dcomplex *ap, dcomplex *x, int *incx);

extern int C2F(ctpsv)(char *uplo, char *trans, char *diag, int *n, 
                      dcomplex *ap, dcomplex *x, int *incx);

extern int C2F(ctrmm)(char *side, char *uplo, char *transa, char *diag, 
                      int *m, int *n, dcomplex *alpha, dcomplex *a, int *lda, 
                      dcomplex *b, int *ldb);

extern int C2F(ctrmv)(char *uplo, char *trans, char *diag, int *n, 
                      dcomplex *a, int *lda, dcomplex *x, int *incx);

extern int C2F(ctrsm)(char *side, char *uplo, char *transa, char *diag, 
                      int *m, int *n, dcomplex *alpha, dcomplex *a, int *lda, 
                      dcomplex *b, int *ldb);

extern int C2F(ctrsv)(char *uplo, char *trans, char *diag, int *n, 
                      dcomplex *a, int *lda, dcomplex *x, int *incx);

double C2F(dasum)(int *n, double *dx, int *incx);

extern int C2F(daxpy)(int *n, double *da, double *dx, 
                      int *incx, double *dy, int *incy);

double C2F(dcabs1)(dcomplex *z__);

extern int C2F(dcopy)(int *n, double *dx, int *incx, 
                      double *dy, int *incy);

double C2F(ddot)(int *n, double *dx, int *incx, double *dy, 
                 int *incy);

extern int C2F(dgbmv)(char *trans, int *m, int *n, int *kl, 
                      int *ku, double *alpha, double *a, int *lda, 
                      double *x, int *incx, double *beta, double *y, 
                      int *incy);

extern int C2F(dgemm)(char *transa, char *transb, int *m, int *
                      n, int *k, double *alpha, double *a, int *lda, 
                      double *b, int *ldb, double *beta, double *c__, 
                      int *ldc);

extern int C2F(dgemv)(char *trans, int *m, int *n, double *
                      alpha, double *a, int *lda, double *x, int *incx, 
                      double *beta, double *y, int *incy);

extern int C2F(dger)(int *m, int *n, double *alpha, 
                     double *x, int *incx, double *y, int *incy, 
                     double *a, int *lda);

double C2F(dnrm2)(int *n, double *x, int *incx);

extern int C2F(drot)(int *n, double *dx, int *incx, 
                     double *dy, int *incy, double *c__, double *s);

extern int C2F(drotg)(double *da, double *db, double *c__, 
                      double *s);

extern int C2F(drotm)(int *n, double *dx, int *incx, 
                      double *dy, int *incy, double *dparam);

extern int C2F(drotmg)(double *dd1, double *dd2, double *
                       dx1, double *dy1, double *dparam);

extern int C2F(dsbmv)(char *uplo, int *n, int *k, double *
                      alpha, double *a, int *lda, double *x, int *incx, 
                      double *beta, double *y, int *incy);

extern int C2F(dscal)(int *n, double *da, double *dx, 
                      int *incx);

double C2F(dsdot)(int *n, double *sx, int *incx, double *sy, int *
                  incy);

extern int C2F(dspmv)(char *uplo, int *n, double *alpha, 
                      double *ap, double *x, int *incx, double *beta, 
                      double *y, int *incy);

extern int C2F(dspr)(char *uplo, int *n, double *alpha, 
                     double *x, int *incx, double *ap);

extern int C2F(dspr2)(char *uplo, int *n, double *alpha, 
                      double *x, int *incx, double *y, int *incy, 
                      double *ap);

extern int C2F(dswap)(int *n, double *dx, int *incx, 
                      double *dy, int *incy);

extern int C2F(dsymm)(char *side, char *uplo, int *m, int *n, 
                      double *alpha, double *a, int *lda, double *b, 
                      int *ldb, double *beta, double *c__, int *ldc);

extern int C2F(dsymv)(char *uplo, int *n, double *alpha, 
                      double *a, int *lda, double *x, int *incx, double 
                      *beta, double *y, int *incy);

extern int C2F(dsyr)(char *uplo, int *n, double *alpha, 
                     double *x, int *incx, double *a, int *lda);

extern int C2F(dsyr2)(char *uplo, int *n, double *alpha, 
                      double *x, int *incx, double *y, int *incy, 
                      double *a, int *lda);

extern int C2F(dsyr2k)(char *uplo, char *trans, int *n, int *k, 
                       double *alpha, double *a, int *lda, double *b, 
                       int *ldb, double *beta, double *c__, int *ldc);

extern int C2F(dsyrk)(char *uplo, char *trans, int *n, int *k, 
                      double *alpha, double *a, int *lda, double *beta, 
                      double *c__, int *ldc);

extern int C2F(dtbmv)(char *uplo, char *trans, char *diag, int *n, 
                      int *k, double *a, int *lda, double *x, int *incx);

extern int C2F(dtbsv)(char *uplo, char *trans, char *diag, int *n, 
                      int *k, double *a, int *lda, double *x, int *incx);

extern int C2F(dtpmv)(char *uplo, char *trans, char *diag, int *n, 
                      double *ap, double *x, int *incx);

extern int C2F(dtpsv)(char *uplo, char *trans, char *diag, int *n, 
                      double *ap, double *x, int *incx);

extern int C2F(dtrmm)(char *side, char *uplo, char *transa, char *diag, 
                      int *m, int *n, double *alpha, double *a, int *
                      lda, double *b, int *ldb);

extern int C2F(dtrmv)(char *uplo, char *trans, char *diag, int *n, 
                      double *a, int *lda, double *x, int *incx);

extern int C2F(dtrsm)(char *side, char *uplo, char *transa, char *diag, 
                      int *m, int *n, double *alpha, double *a, int *
                      lda, double *b, int *ldb);

extern int C2F(dtrsv)(char *uplo, char *trans, char *diag, int *n, 
                      double *a, int *lda, double *x, int *incx);

double C2F(dzasum)(int *n, dcomplex *zx, int *incx);

double C2F(dznrm2)(int *n, dcomplex *x, int *incx);

int C2F(icamax)(int *n, dcomplex *cx, int *incx);

int C2F(idamax)(int *n, double *dx, int *incx);

int C2F(isamax)(int *n, double *sx, int *incx);

int C2F(izamax)(int *n, dcomplex *zx, int *incx);

int C2F(lsame)(char *ca, char *cb);

double C2F(sasum)(int *n, double *sx, int *incx);

extern int C2F(saxpy)(int *n, double *sa, double *sx, int *incx, 
                      double *sy, int *incy);

double C2F(scabs1)(dcomplex *z__);

double C2F(scasum)(int *n, dcomplex *cx, int *incx);

double C2F(scnrm2)(int *n, dcomplex *x, int *incx);

extern int C2F(scopy)(int *n, double *sx, int *incx, double *sy, 
                      int *incy);

double C2F(sdot)(int *n, double *sx, int *incx, double *sy, int *incy);

double C2F(sdsdot)(int *n, double *sb, double *sx, int *incx, double *sy, 
                   int *incy);

extern int C2F(sgbmv)(char *trans, int *m, int *n, int *kl, 
                      int *ku, double *alpha, double *a, int *lda, double *x, int *
                      incx, double *beta, double *y, int *incy);

extern int C2F(sgemm)(char *transa, char *transb, int *m, int *
                      n, int *k, double *alpha, double *a, int *lda, double *b, int *
                      ldb, double *beta, double *c__, int *ldc);

extern int C2F(sgemv)(char *trans, int *m, int *n, double *alpha, 
                      double *a, int *lda, double *x, int *incx, double *beta, double *y, 
                      int *incy);

extern int C2F(sger)(int *m, int *n, double *alpha, double *x, 
                     int *incx, double *y, int *incy, double *a, int *lda);

double C2F(snrm2)(int *n, double *x, int *incx);

extern int C2F(srot)(int *n, double *sx, int *incx, double *sy, 
                     int *incy, double *c__, double *s);

extern int C2F(srotg)(double *sa, double *sb, double *c__, double *s);

extern int C2F(srotm)(int *n, double *sx, int *incx, double *sy, 
                      int *incy, double *sparam);

extern int C2F(srotmg)(double *sd1, double *sd2, double *sx1, double *sy1, double 
                       *sparam);

extern int C2F(ssbmv)(char *uplo, int *n, int *k, double *alpha, 
                      double *a, int *lda, double *x, int *incx, double *beta, double *y, 
                      int *incy);

extern int C2F(sscal)(int *n, double *sa, double *sx, int *incx);

extern int C2F(sspmv)(char *uplo, int *n, double *alpha, double *ap, 
                      double *x, int *incx, double *beta, double *y, int *incy);

extern int C2F(sspr)(char *uplo, int *n, double *alpha, double *x, 
                     int *incx, double *ap);

extern int C2F(sspr2)(char *uplo, int *n, double *alpha, double *x, 
                      int *incx, double *y, int *incy, double *ap);

extern int C2F(sswap)(int *n, double *sx, int *incx, double *sy, 
                      int *incy);

extern int C2F(ssymm)(char *side, char *uplo, int *m, int *n, 
                      double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, 
                      double *c__, int *ldc);

extern int C2F(ssymv)(char *uplo, int *n, double *alpha, double *a, 
                      int *lda, double *x, int *incx, double *beta, double *y, int *
                      incy);

extern int C2F(ssyr)(char *uplo, int *n, double *alpha, double *x, 
                     int *incx, double *a, int *lda);

extern int C2F(ssyr2)(char *uplo, int *n, double *alpha, double *x, 
                      int *incx, double *y, int *incy, double *a, int *lda);

extern int C2F(ssyr2k)(char *uplo, char *trans, int *n, int *k, 
                       double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, 
                       double *c__, int *ldc);

extern int C2F(ssyrk)(char *uplo, char *trans, int *n, int *k, 
                      double *alpha, double *a, int *lda, double *beta, double *c__, int *
                      ldc);

extern int C2F(stbmv)(char *uplo, char *trans, char *diag, int *n, 
                      int *k, double *a, int *lda, double *x, int *incx);

extern int C2F(stbsv)(char *uplo, char *trans, char *diag, int *n, 
                      int *k, double *a, int *lda, double *x, int *incx);

extern int C2F(stpmv)(char *uplo, char *trans, char *diag, int *n, 
                      double *ap, double *x, int *incx);

extern int C2F(stpsv)(char *uplo, char *trans, char *diag, int *n, 
                      double *ap, double *x, int *incx);

extern int C2F(strmm)(char *side, char *uplo, char *transa, char *diag, 
                      int *m, int *n, double *alpha, double *a, int *lda, double *b, 
                      int *ldb);

extern int C2F(strmv)(char *uplo, char *trans, char *diag, int *n, 
                      double *a, int *lda, double *x, int *incx);

extern int C2F(strsm)(char *side, char *uplo, char *transa, char *diag, 
                      int *m, int *n, double *alpha, double *a, int *lda, double *b, 
                      int *ldb);

extern int C2F(strsv)(char *uplo, char *trans, char *diag, int *n, 
                      double *a, int *lda, double *x, int *incx);

extern int C2F(xerbla)(char *srname, int *info);

extern int C2F(xerbla_array_)(char *srname_array__, int *srname_len__, int *info);

extern int C2F(zaxpy)(int *n, dcomplex *za, dcomplex *zx, 
                      int *incx, dcomplex *zy, int *incy);

extern int C2F(zcopy)(int *n, dcomplex *zx, int *incx, 
                      dcomplex *zy, int *incy);

void C2F(zdotc)(dcomplex * ret_val, int *n, 
                dcomplex *zx, int *incx, dcomplex *zy, int *incy);

void C2F(zdotu)(dcomplex * ret_val, int *n, 
                dcomplex *zx, int *incx, dcomplex *zy, int *incy);

extern int C2F(zdrot)(int *n, dcomplex *cx, int *incx, 
                      dcomplex *cy, int *incy, double *c__, double *s);

extern int C2F(zdscal)(int *n, double *da, dcomplex *zx, 
                       int *incx);

extern int C2F(zgbmv)(char *trans, int *m, int *n, int *kl, 
                      int *ku, dcomplex *alpha, dcomplex *a, int *lda, 
                      dcomplex *x, int *incx, dcomplex *beta, dcomplex *
                      y, int *incy);

extern int C2F(zgemm)(char *transa, char *transb, int *m, int *
                      n, int *k, dcomplex *alpha, dcomplex *a, int *lda, 
                      dcomplex *b, int *ldb, dcomplex *beta, dcomplex *
                      c__, int *ldc);

extern int C2F(zgemv)(char *trans, int *m, int *n, 
                      dcomplex *alpha, dcomplex *a, int *lda, dcomplex *
                      x, int *incx, dcomplex *beta, dcomplex *y, int *
                      incy);

extern int C2F(zgerc)(int *m, int *n, dcomplex *alpha, 
                      dcomplex *x, int *incx, dcomplex *y, int *incy, 
                      dcomplex *a, int *lda);

extern int C2F(zgeru)(int *m, int *n, dcomplex *alpha, 
                      dcomplex *x, int *incx, dcomplex *y, int *incy, 
                      dcomplex *a, int *lda);

extern int C2F(zhbmv)(char *uplo, int *n, int *k, dcomplex 
                      *alpha, dcomplex *a, int *lda, dcomplex *x, int *
                      incx, dcomplex *beta, dcomplex *y, int *incy);

extern int C2F(zhemm)(char *side, char *uplo, int *m, int *n, 
                      dcomplex *alpha, dcomplex *a, int *lda, dcomplex *
                      b, int *ldb, dcomplex *beta, dcomplex *c__, int *
                      ldc);

extern int C2F(zhemv)(char *uplo, int *n, dcomplex *alpha, 
                      dcomplex *a, int *lda, dcomplex *x, int *incx, 
                      dcomplex *beta, dcomplex *y, int *incy);

extern int C2F(zher)(char *uplo, int *n, double *alpha, 
                     dcomplex *x, int *incx, dcomplex *a, int *lda);

extern int C2F(zher2)(char *uplo, int *n, dcomplex *alpha, 
                      dcomplex *x, int *incx, dcomplex *y, int *incy, 
                      dcomplex *a, int *lda);

extern int C2F(zher2k)(char *uplo, char *trans, int *n, int *k, 
                       dcomplex *alpha, dcomplex *a, int *lda, dcomplex *
                       b, int *ldb, double *beta, dcomplex *c__, int *ldc);

extern int C2F(zherk)(char *uplo, char *trans, int *n, int *k, 
                      double *alpha, dcomplex *a, int *lda, double *beta, 
                      dcomplex *c__, int *ldc);

extern int C2F(zhpmv)(char *uplo, int *n, dcomplex *alpha, 
                      dcomplex *ap, dcomplex *x, int *incx, dcomplex *
                      beta, dcomplex *y, int *incy);

extern int C2F(zhpr)(char *uplo, int *n, double *alpha, 
                     dcomplex *x, int *incx, dcomplex *ap);

extern int C2F(zhpr2)(char *uplo, int *n, dcomplex *alpha, 
                      dcomplex *x, int *incx, dcomplex *y, int *incy, 
                      dcomplex *ap);

extern int C2F(zrotg)(dcomplex *ca, dcomplex *cb, double *
                      c__, dcomplex *s);

extern int C2F(zscal)(int *n, dcomplex *za, dcomplex *zx, 
                      int *incx);

extern int C2F(zswap)(int *n, dcomplex *zx, int *incx, 
                      dcomplex *zy, int *incy);

extern int C2F(zsymm)(char *side, char *uplo, int *m, int *n, 
                      dcomplex *alpha, dcomplex *a, int *lda, dcomplex *
                      b, int *ldb, dcomplex *beta, dcomplex *c__, int *
                      ldc);

extern int C2F(zsyr2k)(char *uplo, char *trans, int *n, int *k, 
                       dcomplex *alpha, dcomplex *a, int *lda, dcomplex *
                       b, int *ldb, dcomplex *beta, dcomplex *c__, int *
                       ldc);

extern int C2F(zsyrk)(char *uplo, char *trans, int *n, int *k, 
                      dcomplex *alpha, dcomplex *a, int *lda, dcomplex *
                      beta, dcomplex *c__, int *ldc);

extern int C2F(ztbmv)(char *uplo, char *trans, char *diag, int *n, 
                      int *k, dcomplex *a, int *lda, dcomplex *x, int 
                      *incx);

extern int C2F(ztbsv)(char *uplo, char *trans, char *diag, int *n, 
                      int *k, dcomplex *a, int *lda, dcomplex *x, int 
                      *incx);

extern int C2F(ztpmv)(char *uplo, char *trans, char *diag, int *n, 
                      dcomplex *ap, dcomplex *x, int *incx);

extern int C2F(ztpsv)(char *uplo, char *trans, char *diag, int *n, 
                      dcomplex *ap, dcomplex *x, int *incx);

extern int C2F(ztrmm)(char *side, char *uplo, char *transa, char *diag, 
                      int *m, int *n, dcomplex *alpha, dcomplex *a, 
                      int *lda, dcomplex *b, int *ldb);

extern int C2F(ztrmv)(char *uplo, char *trans, char *diag, int *n, 
                      dcomplex *a, int *lda, dcomplex *x, int *incx);

extern int C2F(ztrsm)(char *side, char *uplo, char *transa, char *diag, 
                      int *m, int *n, dcomplex *alpha, dcomplex *a, 
                      int *lda, dcomplex *b, int *ldb);

extern int C2F(ztrsv)(char *uplo, char *trans, char *diag, int *n, 
                      dcomplex *a, int *lda, dcomplex *x, int *incx);

extern int C2F(cbdsqr)(char *uplo, int *n, int *ncvt, int *
                       nru, int *ncc, double *d__, double *e, dcomplex *vt, int *ldvt, 
                       dcomplex *u, int *ldu, dcomplex *c__, int *ldc, double *rwork, 
                       int *info);

extern int C2F(cgbbrd)(char *vect, int *m, int *n, int *ncc, 
                       int *kl, int *ku, dcomplex *ab, int *ldab, double *d__, 
                       double *e, dcomplex *q, int *ldq, dcomplex *pt, int *ldpt, 
                       dcomplex *c__, int *ldc, dcomplex *work, double *rwork, int *info);

extern int C2F(cgbcon)(char *norm, int *n, int *kl, int *ku, 
                       dcomplex *ab, int *ldab, int *ipiv, double *anorm, double *rcond, 
                       dcomplex *work, double *rwork, int *info);

extern int C2F(cgbequ)(int *m, int *n, int *kl, int *ku, 
                       dcomplex *ab, int *ldab, double *r__, double *c__, double *rowcnd, double 
                       *colcnd, double *amax, int *info);

extern int C2F(cgbequb)(int *m, int *n, int *kl, int *
                        ku, dcomplex *ab, int *ldab, double *r__, double *c__, double *rowcnd, 
                        double *colcnd, double *amax, int *info);

extern int C2F(cgbrfs)(char *trans, int *n, int *kl, int *
                       ku, int *nrhs, dcomplex *ab, int *ldab, dcomplex *afb, int *
                       ldafb, int *ipiv, dcomplex *b, int *ldb, dcomplex *x, int *
                       ldx, double *ferr, double *berr, dcomplex *work, double *rwork, int *
                       info);

extern int C2F(cgbrfsx)(char *trans, char *equed, int *n, int *
                        kl, int *ku, int *nrhs, dcomplex *ab, int *ldab, dcomplex *
                        afb, int *ldafb, int *ipiv, double *r__, double *c__, dcomplex *b, 
                        int *ldb, dcomplex *x, int *ldx, double *rcond, double *berr, 
                        int *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__, 
                        int *nparams, double *params, dcomplex *work, double *rwork, int *
                        info);

extern int C2F(cgbsv)(int *n, int *kl, int *ku, int *
                      nrhs, dcomplex *ab, int *ldab, int *ipiv, dcomplex *b, int *
                      ldb, int *info);

extern int C2F(cgbsvx)(char *fact, char *trans, int *n, int *kl, 
                       int *ku, int *nrhs, dcomplex *ab, int *ldab, dcomplex *afb, 
                       int *ldafb, int *ipiv, char *equed, double *r__, double *c__, 
                       dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double 
                       *ferr, double *berr, dcomplex *work, double *rwork, int *info);

extern int C2F(cgbsvxx)(char *fact, char *trans, int *n, int *
                        kl, int *ku, int *nrhs, dcomplex *ab, int *ldab, dcomplex *
                        afb, int *ldafb, int *ipiv, char *equed, double *r__, double *c__, 
                        dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, 
                        double *rpvgrw, double *berr, int *n_err_bnds__, double *
                        err_bnds_norm__, double *err_bnds_comp__, int *nparams, double *
                        params, dcomplex *work, double *rwork, int *info);

extern int C2F(cgbtf2)(int *m, int *n, int *kl, int *ku, 
                       dcomplex *ab, int *ldab, int *ipiv, int *info);

extern int C2F(cgbtrf)(int *m, int *n, int *kl, int *ku, 
                       dcomplex *ab, int *ldab, int *ipiv, int *info);

extern int C2F(cgbtrs)(char *trans, int *n, int *kl, int *
                       ku, int *nrhs, dcomplex *ab, int *ldab, int *ipiv, dcomplex 
                       *b, int *ldb, int *info);

extern int C2F(cgebak)(char *job, char *side, int *n, int *ilo, 
                       int *ihi, double *scale, int *m, dcomplex *v, int *ldv, 
                       int *info);

extern int C2F(cgebal)(char *job, int *n, dcomplex *a, int *lda, 
                       int *ilo, int *ihi, double *scale, int *info);

extern int C2F(cgebd2)(int *m, int *n, dcomplex *a, int *lda, 
                       double *d__, double *e, dcomplex *tauq, dcomplex *taup, dcomplex *work, 
                       int *info);

extern int C2F(cgebrd)(int *m, int *n, dcomplex *a, int *lda, 
                       double *d__, double *e, dcomplex *tauq, dcomplex *taup, dcomplex *work, 
                       int *lwork, int *info);

extern int C2F(cgecon)(char *norm, int *n, dcomplex *a, int *lda, 
                       double *anorm, double *rcond, dcomplex *work, double *rwork, int *info);

extern int C2F(cgeequ)(int *m, int *n, dcomplex *a, int *lda, 
                       double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, 
                       int *info);

extern int C2F(cgeequb)(int *m, int *n, dcomplex *a, int *
                        lda, double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, 
                        int *info);

extern int C2F(cgees)(char *jobvs, char *sort, L_fp select, int *n, 
                      dcomplex *a, int *lda, int *sdim, dcomplex *w, dcomplex *vs, 
                      int *ldvs, dcomplex *work, int *lwork, double *rwork, int *
                      bwork, int *info);

extern int C2F(cgeesx)(char *jobvs, char *sort, L_fp select, char *
                       sense, int *n, dcomplex *a, int *lda, int *sdim, dcomplex *
                       w, dcomplex *vs, int *ldvs, double *rconde, double *rcondv, dcomplex *
                       work, int *lwork, double *rwork, int *bwork, int *info);

extern int C2F(cgeev)(char *jobvl, char *jobvr, int *n, dcomplex *a, 
                      int *lda, dcomplex *w, dcomplex *vl, int *ldvl, dcomplex *vr, 
                      int *ldvr, dcomplex *work, int *lwork, double *rwork, int *
                      info);

extern int C2F(cgeevx)(char *balanc, char *jobvl, char *jobvr, char *
                       sense, int *n, dcomplex *a, int *lda, dcomplex *w, dcomplex *vl, 
                       int *ldvl, dcomplex *vr, int *ldvr, int *ilo, int *ihi, 
                       double *scale, double *abnrm, double *rconde, double *rcondv, dcomplex *work, 
                       int *lwork, double *rwork, int *info);

extern int C2F(cgegs)(char *jobvsl, char *jobvsr, int *n, dcomplex *
                      a, int *lda, dcomplex *b, int *ldb, dcomplex *alpha, dcomplex *
                      beta, dcomplex *vsl, int *ldvsl, dcomplex *vsr, int *ldvsr, 
                      dcomplex *work, int *lwork, double *rwork, int *info);

extern int C2F(cgegv)(char *jobvl, char *jobvr, int *n, dcomplex *a, 
                      int *lda, dcomplex *b, int *ldb, dcomplex *alpha, dcomplex *beta, 
                      dcomplex *vl, int *ldvl, dcomplex *vr, int *ldvr, dcomplex *
                      work, int *lwork, double *rwork, int *info);

extern int C2F(cgehd2)(int *n, int *ilo, int *ihi, dcomplex *
                       a, int *lda, dcomplex *tau, dcomplex *work, int *info);

extern int C2F(cgehrd)(int *n, int *ilo, int *ihi, dcomplex *
                       a, int *lda, dcomplex *tau, dcomplex *work, int *lwork, int 
                       *info);

extern int C2F(cgelq2)(int *m, int *n, dcomplex *a, int *lda, 
                       dcomplex *tau, dcomplex *work, int *info);

extern int C2F(cgelqf)(int *m, int *n, dcomplex *a, int *lda, 
                       dcomplex *tau, dcomplex *work, int *lwork, int *info);

extern int C2F(cgels)(char *trans, int *m, int *n, int *
                      nrhs, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *
                      work, int *lwork, int *info);

extern int C2F(cgelsd)(int *m, int *n, int *nrhs, dcomplex *
                       a, int *lda, dcomplex *b, int *ldb, double *s, double *rcond, 
                       int *rank, dcomplex *work, int *lwork, double *rwork, int *
                       iwork, int *info);

extern int C2F(cgelss)(int *m, int *n, int *nrhs, dcomplex *
                       a, int *lda, dcomplex *b, int *ldb, double *s, double *rcond, 
                       int *rank, dcomplex *work, int *lwork, double *rwork, int *
                       info);

extern int C2F(cgelsx)(int *m, int *n, int *nrhs, dcomplex *
                       a, int *lda, dcomplex *b, int *ldb, int *jpvt, double *rcond, 
                       int *rank, dcomplex *work, double *rwork, int *info);

extern int C2F(cgelsy)(int *m, int *n, int *nrhs, dcomplex *
                       a, int *lda, dcomplex *b, int *ldb, int *jpvt, double *rcond, 
                       int *rank, dcomplex *work, int *lwork, double *rwork, int *
                       info);

extern int C2F(cgeql2)(int *m, int *n, dcomplex *a, int *lda, 
                       dcomplex *tau, dcomplex *work, int *info);

extern int C2F(cgeqlf)(int *m, int *n, dcomplex *a, int *lda, 
                       dcomplex *tau, dcomplex *work, int *lwork, int *info);

extern int C2F(cgeqp3)(int *m, int *n, dcomplex *a, int *lda, 
                       int *jpvt, dcomplex *tau, dcomplex *work, int *lwork, double *
                       rwork, int *info);

extern int C2F(cgeqpf)(int *m, int *n, dcomplex *a, int *lda, 
                       int *jpvt, dcomplex *tau, dcomplex *work, double *rwork, int *
                       info);

extern int C2F(cgeqr2)(int *m, int *n, dcomplex *a, int *lda, 
                       dcomplex *tau, dcomplex *work, int *info);

extern int C2F(cgeqrf)(int *m, int *n, dcomplex *a, int *lda, 
                       dcomplex *tau, dcomplex *work, int *lwork, int *info);

extern int C2F(cgerfs)(char *trans, int *n, int *nrhs, dcomplex *
                       a, int *lda, dcomplex *af, int *ldaf, int *ipiv, dcomplex *
                       b, int *ldb, dcomplex *x, int *ldx, double *ferr, double *berr, 
                       dcomplex *work, double *rwork, int *info);

extern int C2F(cgerfsx)(char *trans, char *equed, int *n, int *
                        nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, int *
                        ipiv, double *r__, double *c__, dcomplex *b, int *ldb, dcomplex *x, 
                        int *ldx, double *rcond, double *berr, int *n_err_bnds__, double *
                        err_bnds_norm__, double *err_bnds_comp__, int *nparams, double *
                        params, dcomplex *work, double *rwork, int *info);

extern int C2F(cgerq2)(int *m, int *n, dcomplex *a, int *lda, 
                       dcomplex *tau, dcomplex *work, int *info);

extern int C2F(cgerqf)(int *m, int *n, dcomplex *a, int *lda, 
                       dcomplex *tau, dcomplex *work, int *lwork, int *info);

extern int C2F(cgesc2)(int *n, dcomplex *a, int *lda, dcomplex *
                       rhs, int *ipiv, int *jpiv, double *scale);

extern int C2F(cgesdd)(char *jobz, int *m, int *n, dcomplex *a, 
                       int *lda, double *s, dcomplex *u, int *ldu, dcomplex *vt, int 
                       *ldvt, dcomplex *work, int *lwork, double *rwork, int *iwork, 
                       int *info);

extern int C2F(cgesv)(int *n, int *nrhs, dcomplex *a, int *
                      lda, int *ipiv, dcomplex *b, int *ldb, int *info);

extern int C2F(cgesvd)(char *jobu, char *jobvt, int *m, int *n, 
                       dcomplex *a, int *lda, double *s, dcomplex *u, int *ldu, dcomplex *
                       vt, int *ldvt, dcomplex *work, int *lwork, double *rwork, 
                       int *info);

extern int C2F(cgesvx)(char *fact, char *trans, int *n, int *
                       nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, int *
                       ipiv, char *equed, double *r__, double *c__, dcomplex *b, int *ldb, 
                       dcomplex *x, int *ldx, double *rcond, double *ferr, double *berr, 
                       dcomplex *work, double *rwork, int *info);

extern int C2F(cgesvxx)(char *fact, char *trans, int *n, int *
                        nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, int *
                        ipiv, char *equed, double *r__, double *c__, dcomplex *b, int *ldb, 
                        dcomplex *x, int *ldx, double *rcond, double *rpvgrw, double *berr, 
                        int *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__, 
                        int *nparams, double *params, dcomplex *work, double *rwork, int *
                        info);

extern int C2F(cgetc2)(int *n, dcomplex *a, int *lda, int *
                       ipiv, int *jpiv, int *info);

extern int C2F(cgetf2)(int *m, int *n, dcomplex *a, int *lda, 
                       int *ipiv, int *info);

extern int C2F(cgetrf)(int *m, int *n, dcomplex *a, int *lda, 
                       int *ipiv, int *info);

extern int C2F(cgetri)(int *n, dcomplex *a, int *lda, int *
                       ipiv, dcomplex *work, int *lwork, int *info);

extern int C2F(cgetrs)(char *trans, int *n, int *nrhs, dcomplex *
                       a, int *lda, int *ipiv, dcomplex *b, int *ldb, int *
                       info);

extern int C2F(cggbak)(char *job, char *side, int *n, int *ilo, 
                       int *ihi, double *lscale, double *rscale, int *m, dcomplex *v, 
                       int *ldv, int *info);

extern int C2F(cggbal)(char *job, int *n, dcomplex *a, int *lda, 
                       dcomplex *b, int *ldb, int *ilo, int *ihi, double *lscale, 
                       double *rscale, double *work, int *info);

extern int C2F(cgges)(char *jobvsl, char *jobvsr, char *sort, L_fp 
                      selctg, int *n, dcomplex *a, int *lda, dcomplex *b, int *
                      ldb, int *sdim, dcomplex *alpha, dcomplex *beta, dcomplex *vsl, 
                      int *ldvsl, dcomplex *vsr, int *ldvsr, dcomplex *work, int *
                      lwork, double *rwork, int *bwork, int *info);

extern int C2F(cggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp 
                       selctg, char *sense, int *n, dcomplex *a, int *lda, dcomplex *b, 
                       int *ldb, int *sdim, dcomplex *alpha, dcomplex *beta, dcomplex *
                       vsl, int *ldvsl, dcomplex *vsr, int *ldvsr, double *rconde, double 
                       *rcondv, dcomplex *work, int *lwork, double *rwork, int *iwork, 
                       int *liwork, int *bwork, int *info);

extern int C2F(cggev)(char *jobvl, char *jobvr, int *n, dcomplex *a, 
                      int *lda, dcomplex *b, int *ldb, dcomplex *alpha, dcomplex *beta, 
                      dcomplex *vl, int *ldvl, dcomplex *vr, int *ldvr, dcomplex *
                      work, int *lwork, double *rwork, int *info);

extern int C2F(cggevx)(char *balanc, char *jobvl, char *jobvr, char *
                       sense, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       dcomplex *alpha, dcomplex *beta, dcomplex *vl, int *ldvl, dcomplex *
                       vr, int *ldvr, int *ilo, int *ihi, double *lscale, double *
                       rscale, double *abnrm, double *bbnrm, double *rconde, double *rcondv, dcomplex 
                       *work, int *lwork, double *rwork, int *iwork, int *bwork, 
                       int *info);

extern int C2F(cggglm)(int *n, int *m, int *p, dcomplex *a, 
                       int *lda, dcomplex *b, int *ldb, dcomplex *d__, dcomplex *x, 
                       dcomplex *y, dcomplex *work, int *lwork, int *info);

extern int C2F(cgghrd)(char *compq, char *compz, int *n, int *
                       ilo, int *ihi, dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       dcomplex *q, int *ldq, dcomplex *z__, int *ldz, int *info);

extern int C2F(cgglse)(int *m, int *n, int *p, dcomplex *a, 
                       int *lda, dcomplex *b, int *ldb, dcomplex *c__, dcomplex *d__, 
                       dcomplex *x, dcomplex *work, int *lwork, int *info);

extern int C2F(cggqrf)(int *n, int *m, int *p, dcomplex *a, 
                       int *lda, dcomplex *taua, dcomplex *b, int *ldb, dcomplex *taub, 
                       dcomplex *work, int *lwork, int *info);

extern int C2F(cggrqf)(int *m, int *p, int *n, dcomplex *a, 
                       int *lda, dcomplex *taua, dcomplex *b, int *ldb, dcomplex *taub, 
                       dcomplex *work, int *lwork, int *info);

extern int C2F(cggsvd)(char *jobu, char *jobv, char *jobq, int *m, 
                       int *n, int *p, int *k, int *l, dcomplex *a, int *
                       lda, dcomplex *b, int *ldb, double *alpha, double *beta, dcomplex *u, 
                       int *ldu, dcomplex *v, int *ldv, dcomplex *q, int *ldq, 
                       dcomplex *work, double *rwork, int *iwork, int *info);

extern int C2F(cggsvp)(char *jobu, char *jobv, char *jobq, int *m, 
                       int *p, int *n, dcomplex *a, int *lda, dcomplex *b, int 
                       *ldb, double *tola, double *tolb, int *k, int *l, dcomplex *u, 
                       int *ldu, dcomplex *v, int *ldv, dcomplex *q, int *ldq, 
                       int *iwork, double *rwork, dcomplex *tau, dcomplex *work, int *
                       info);

extern int C2F(cgtcon)(char *norm, int *n, dcomplex *dl, dcomplex *
                       d__, dcomplex *du, dcomplex *du2, int *ipiv, double *anorm, double *
                       rcond, dcomplex *work, int *info);

extern int C2F(cgtrfs)(char *trans, int *n, int *nrhs, dcomplex *
                       dl, dcomplex *d__, dcomplex *du, dcomplex *dlf, dcomplex *df, dcomplex *
                       duf, dcomplex *du2, int *ipiv, dcomplex *b, int *ldb, dcomplex *
                       x, int *ldx, double *ferr, double *berr, dcomplex *work, double *rwork, 
                       int *info);

extern int C2F(cgtsv)(int *n, int *nrhs, dcomplex *dl, dcomplex *
                      d__, dcomplex *du, dcomplex *b, int *ldb, int *info);

extern int C2F(cgtsvx)(char *fact, char *trans, int *n, int *
                       nrhs, dcomplex *dl, dcomplex *d__, dcomplex *du, dcomplex *dlf, dcomplex *
                       df, dcomplex *duf, dcomplex *du2, int *ipiv, dcomplex *b, int *
                       ldb, dcomplex *x, int *ldx, double *rcond, double *ferr, double *berr, 
                       dcomplex *work, double *rwork, int *info);

extern int C2F(cgttrf)(int *n, dcomplex *dl, dcomplex *d__, dcomplex *
                       du, dcomplex *du2, int *ipiv, int *info);

extern int C2F(cgttrs)(char *trans, int *n, int *nrhs, dcomplex *
                       dl, dcomplex *d__, dcomplex *du, dcomplex *du2, int *ipiv, dcomplex *
                       b, int *ldb, int *info);

extern int C2F(cgtts2)(int *itrans, int *n, int *nrhs, 
                       dcomplex *dl, dcomplex *d__, dcomplex *du, dcomplex *du2, int *ipiv, 
                       dcomplex *b, int *ldb);

extern int C2F(chbev)(char *jobz, char *uplo, int *n, int *kd, 
                      dcomplex *ab, int *ldab, double *w, dcomplex *z__, int *ldz, 
                      dcomplex *work, double *rwork, int *info);

extern int C2F(chbevd)(char *jobz, char *uplo, int *n, int *kd, 
                       dcomplex *ab, int *ldab, double *w, dcomplex *z__, int *ldz, 
                       dcomplex *work, int *lwork, double *rwork, int *lrwork, int *
                       iwork, int *liwork, int *info);

extern int C2F(chbevx)(char *jobz, char *range, char *uplo, int *n, 
                       int *kd, dcomplex *ab, int *ldab, dcomplex *q, int *ldq, 
                       double *vl, double *vu, int *il, int *iu, double *abstol, int *
                       m, double *w, dcomplex *z__, int *ldz, dcomplex *work, double *rwork, 
                       int *iwork, int *ifail, int *info);

extern int C2F(chbgst)(char *vect, char *uplo, int *n, int *ka, 
                       int *kb, dcomplex *ab, int *ldab, dcomplex *bb, int *ldbb, 
                       dcomplex *x, int *ldx, dcomplex *work, double *rwork, int *info);

extern int C2F(chbgv)(char *jobz, char *uplo, int *n, int *ka, 
                      int *kb, dcomplex *ab, int *ldab, dcomplex *bb, int *ldbb, 
                      double *w, dcomplex *z__, int *ldz, dcomplex *work, double *rwork, 
                      int *info);

extern int C2F(chbgvd)(char *jobz, char *uplo, int *n, int *ka, 
                       int *kb, dcomplex *ab, int *ldab, dcomplex *bb, int *ldbb, 
                       double *w, dcomplex *z__, int *ldz, dcomplex *work, int *lwork, 
                       double *rwork, int *lrwork, int *iwork, int *liwork, 
                       int *info);

extern int C2F(chbgvx)(char *jobz, char *range, char *uplo, int *n, 
                       int *ka, int *kb, dcomplex *ab, int *ldab, dcomplex *bb, 
                       int *ldbb, dcomplex *q, int *ldq, double *vl, double *vu, int *
                       il, int *iu, double *abstol, int *m, double *w, dcomplex *z__, 
                       int *ldz, dcomplex *work, double *rwork, int *iwork, int *
                       ifail, int *info);

extern int C2F(chbtrd)(char *vect, char *uplo, int *n, int *kd, 
                       dcomplex *ab, int *ldab, double *d__, double *e, dcomplex *q, int *
                       ldq, dcomplex *work, int *info);

extern int C2F(checon)(char *uplo, int *n, dcomplex *a, int *lda, 
                       int *ipiv, double *anorm, double *rcond, dcomplex *work, int *
                       info);

extern int C2F(cheequb)(char *uplo, int *n, dcomplex *a, int *
                        lda, double *s, double *scond, double *amax, dcomplex *work, int *info);

extern int C2F(cheev)(char *jobz, char *uplo, int *n, dcomplex *a, 
                      int *lda, double *w, dcomplex *work, int *lwork, double *rwork, 
                      int *info);

extern int C2F(cheevd)(char *jobz, char *uplo, int *n, dcomplex *a, 
                       int *lda, double *w, dcomplex *work, int *lwork, double *rwork, 
                       int *lrwork, int *iwork, int *liwork, int *info);

extern int C2F(cheevr)(char *jobz, char *range, char *uplo, int *n, 
                       dcomplex *a, int *lda, double *vl, double *vu, int *il, int *
                       iu, double *abstol, int *m, double *w, dcomplex *z__, int *ldz, 
                       int *isuppz, dcomplex *work, int *lwork, double *rwork, int *
                       lrwork, int *iwork, int *liwork, int *info);

extern int C2F(cheevx)(char *jobz, char *range, char *uplo, int *n, 
                       dcomplex *a, int *lda, double *vl, double *vu, int *il, int *
                       iu, double *abstol, int *m, double *w, dcomplex *z__, int *ldz, 
                       dcomplex *work, int *lwork, double *rwork, int *iwork, int *
                       ifail, int *info);

extern int C2F(chegs2)(int *itype, char *uplo, int *n, dcomplex *
                       a, int *lda, dcomplex *b, int *ldb, int *info);

extern int C2F(chegst)(int *itype, char *uplo, int *n, dcomplex *
                       a, int *lda, dcomplex *b, int *ldb, int *info);

extern int C2F(chegv)(int *itype, char *jobz, char *uplo, int *
                      n, dcomplex *a, int *lda, dcomplex *b, int *ldb, double *w, 
                      dcomplex *work, int *lwork, double *rwork, int *info);

extern int C2F(chegvd)(int *itype, char *jobz, char *uplo, int *
                       n, dcomplex *a, int *lda, dcomplex *b, int *ldb, double *w, 
                       dcomplex *work, int *lwork, double *rwork, int *lrwork, int *
                       iwork, int *liwork, int *info);

extern int C2F(chegvx)(int *itype, char *jobz, char *range, char *
                       uplo, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       double *vl, double *vu, int *il, int *iu, double *abstol, int *
                       m, double *w, dcomplex *z__, int *ldz, dcomplex *work, int *lwork, 
                       double *rwork, int *iwork, int *ifail, int *info);

extern int C2F(cherfs)(char *uplo, int *n, int *nrhs, dcomplex *
                       a, int *lda, dcomplex *af, int *ldaf, int *ipiv, dcomplex *
                       b, int *ldb, dcomplex *x, int *ldx, double *ferr, double *berr, 
                       dcomplex *work, double *rwork, int *info);

extern int C2F(cherfsx)(char *uplo, char *equed, int *n, int *
                        nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, int *
                        ipiv, double *s, dcomplex *b, int *ldb, dcomplex *x, int *ldx, 
                        double *rcond, double *berr, int *n_err_bnds__, double *err_bnds_norm__, 
                        double *err_bnds_comp__, int *nparams, double *params, dcomplex *work, 
                        double *rwork, int *info);

extern int C2F(chesv)(char *uplo, int *n, int *nrhs, dcomplex *a, 
                      int *lda, int *ipiv, dcomplex *b, int *ldb, dcomplex *work, 
                      int *lwork, int *info);

extern int C2F(chesvx)(char *fact, char *uplo, int *n, int *
                       nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, int *
                       ipiv, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, 
                       double *ferr, double *berr, dcomplex *work, int *lwork, double *rwork, 
                       int *info);

extern int C2F(chesvxx)(char *fact, char *uplo, int *n, int *
                        nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, int *
                        ipiv, char *equed, double *s, dcomplex *b, int *ldb, dcomplex *x, 
                        int *ldx, double *rcond, double *rpvgrw, double *berr, int *
                        n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__, int *
                        nparams, double *params, dcomplex *work, double *rwork, int *info);

extern int C2F(chetd2)(char *uplo, int *n, dcomplex *a, int *lda, 
                       double *d__, double *e, dcomplex *tau, int *info);

extern int C2F(chetf2)(char *uplo, int *n, dcomplex *a, int *lda, 
                       int *ipiv, int *info);

extern int C2F(chetrd)(char *uplo, int *n, dcomplex *a, int *lda, 
                       double *d__, double *e, dcomplex *tau, dcomplex *work, int *lwork, 
                       int *info);

extern int C2F(chetrf)(char *uplo, int *n, dcomplex *a, int *lda, 
                       int *ipiv, dcomplex *work, int *lwork, int *info);

extern int C2F(chetri)(char *uplo, int *n, dcomplex *a, int *lda, 
                       int *ipiv, dcomplex *work, int *info);

extern int C2F(chetrs)(char *uplo, int *n, int *nrhs, dcomplex *
                       a, int *lda, int *ipiv, dcomplex *b, int *ldb, int *
                       info);

extern int C2F(chfrk)(char *transr, char *uplo, char *trans, int *n, 
                      int *k, double *alpha, dcomplex *a, int *lda, double *beta, 
                      dcomplex *c__);

extern int C2F(chgeqz)(char *job, char *compq, char *compz, int *n, 
                       int *ilo, int *ihi, dcomplex *h__, int *ldh, dcomplex *t, 
                       int *ldt, dcomplex *alpha, dcomplex *beta, dcomplex *q, int *ldq, 
                       dcomplex *z__, int *ldz, dcomplex *work, int *lwork, double *
                       rwork, int *info);

/* Character */ void C2F(chla_transtype_)(char *ret_val, int ret_val_len, 
                                          int *trans);

extern int C2F(chpcon)(char *uplo, int *n, dcomplex *ap, int *
                       ipiv, double *anorm, double *rcond, dcomplex *work, int *info);

extern int C2F(chpev)(char *jobz, char *uplo, int *n, dcomplex *ap, 
                      double *w, dcomplex *z__, int *ldz, dcomplex *work, double *rwork, 
                      int *info);

extern int C2F(chpevd)(char *jobz, char *uplo, int *n, dcomplex *ap, 
                       double *w, dcomplex *z__, int *ldz, dcomplex *work, int *lwork, 
                       double *rwork, int *lrwork, int *iwork, int *liwork, 
                       int *info);

extern int C2F(chpevx)(char *jobz, char *range, char *uplo, int *n, 
                       dcomplex *ap, double *vl, double *vu, int *il, int *iu, double *
                       abstol, int *m, double *w, dcomplex *z__, int *ldz, dcomplex *
                       work, double *rwork, int *iwork, int *ifail, int *info);

extern int C2F(chpgst)(int *itype, char *uplo, int *n, dcomplex *
                       ap, dcomplex *bp, int *info);

extern int C2F(chpgv)(int *itype, char *jobz, char *uplo, int *
                      n, dcomplex *ap, dcomplex *bp, double *w, dcomplex *z__, int *ldz, 
                      dcomplex *work, double *rwork, int *info);

extern int C2F(chpgvd)(int *itype, char *jobz, char *uplo, int *
                       n, dcomplex *ap, dcomplex *bp, double *w, dcomplex *z__, int *ldz, 
                       dcomplex *work, int *lwork, double *rwork, int *lrwork, int *
                       iwork, int *liwork, int *info);

extern int C2F(chpgvx)(int *itype, char *jobz, char *range, char *
                       uplo, int *n, dcomplex *ap, dcomplex *bp, double *vl, double *vu, 
                       int *il, int *iu, double *abstol, int *m, double *w, dcomplex *
                       z__, int *ldz, dcomplex *work, double *rwork, int *iwork, 
                       int *ifail, int *info);

extern int C2F(chprfs)(char *uplo, int *n, int *nrhs, dcomplex *
                       ap, dcomplex *afp, int *ipiv, dcomplex *b, int *ldb, dcomplex *x, 
                       int *ldx, double *ferr, double *berr, dcomplex *work, double *rwork, 
                       int *info);

extern int C2F(chpsv)(char *uplo, int *n, int *nrhs, dcomplex *
                      ap, int *ipiv, dcomplex *b, int *ldb, int *info);

extern int C2F(chpsvx)(char *fact, char *uplo, int *n, int *
                       nrhs, dcomplex *ap, dcomplex *afp, int *ipiv, dcomplex *b, int *
                       ldb, dcomplex *x, int *ldx, double *rcond, double *ferr, double *berr, 
                       dcomplex *work, double *rwork, int *info);

extern int C2F(chptrd)(char *uplo, int *n, dcomplex *ap, double *d__, 
                       double *e, dcomplex *tau, int *info);

extern int C2F(chptrf)(char *uplo, int *n, dcomplex *ap, int *
                       ipiv, int *info);

extern int C2F(chptri)(char *uplo, int *n, dcomplex *ap, int *
                       ipiv, dcomplex *work, int *info);

extern int C2F(chptrs)(char *uplo, int *n, int *nrhs, dcomplex *
                       ap, int *ipiv, dcomplex *b, int *ldb, int *info);

extern int C2F(chsein)(char *side, char *eigsrc, char *initv, int *
                       select, int *n, dcomplex *h__, int *ldh, dcomplex *w, dcomplex *
                       vl, int *ldvl, dcomplex *vr, int *ldvr, int *mm, int *
                       m, dcomplex *work, double *rwork, int *ifaill, int *ifailr, 
                       int *info);

extern int C2F(chseqr)(char *job, char *compz, int *n, int *ilo, 
                       int *ihi, dcomplex *h__, int *ldh, dcomplex *w, dcomplex *z__, 
                       int *ldz, dcomplex *work, int *lwork, int *info);

extern int C2F(cla_gbamv_)(int *trans, int *m, int *n, 
                           int *kl, int *ku, double *alpha, dcomplex *ab, int *ldab, 
                           dcomplex *x, int *incx, double *beta, double *y, int *incy);

double C2F(cla_gbrcond_c_)(char *trans, int *n, int *kl, int *ku, 
                           dcomplex *ab, int *ldab, dcomplex *afb, int *ldafb, int *
                           ipiv, double *c__, int *capply, int *info, dcomplex *work, double *
                           rwork, int trans_len);

double C2F(cla_gbrcond_x_)(char *trans, int *n, int *kl, int *ku, 
                           dcomplex *ab, int *ldab, dcomplex *afb, int *ldafb, int *
                           ipiv, dcomplex *x, int *info, dcomplex *work, double *rwork, int 
                           trans_len);

extern int C2F(cla_gbrfsx_extended_)(int *prec_type__, int *
                                     trans_type__, int *n, int *kl, int *ku, int *nrhs, 
                                     dcomplex *ab, int *ldab, dcomplex *afb, int *ldafb, int *
                                     ipiv, int *colequ, double *c__, dcomplex *b, int *ldb, dcomplex *
                                     y, int *ldy, double *berr_out__, int *n_norms__, double *errs_n__,
                                     double *errs_c__, dcomplex *res, double *ayb, dcomplex *dy, dcomplex *
                                     y_tail__, double *rcond, int *ithresh, double *rthresh, double *dz_ub__,
                                     int *ignore_cwise__, int *info);

double C2F(cla_gbrpvgrw_)(int *n, int *kl, int *ku, int *
                          ncols, dcomplex *ab, int *ldab, dcomplex *afb, int *ldafb);

extern int C2F(cla_geamv_)(int *trans, int *m, int *n, double 
                           *alpha, dcomplex *a, int *lda, dcomplex *x, int *incx, double *
                           beta, double *y, int *incy);

double C2F(cla_gercond_c_)(char *trans, int *n, dcomplex *a, int *lda, 
                           dcomplex *af, int *ldaf, int *ipiv, double *c__, int *capply,
                           int *info, dcomplex *work, double *rwork, int trans_len);

double C2F(cla_gercond_x_)(char *trans, int *n, dcomplex *a, int *lda, 
                           dcomplex *af, int *ldaf, int *ipiv, dcomplex *x, int *info, 
                           dcomplex *work, double *rwork, int trans_len);

extern int C2F(cla_gerfsx_extended_)(int *prec_type__, int *
                                     trans_type__, int *n, int *nrhs, dcomplex *a, int *lda, 
                                     dcomplex *af, int *ldaf, int *ipiv, int *colequ, double *c__,
                                     dcomplex *b, int *ldb, dcomplex *y, int *ldy, double *berr_out__,
                                     int *n_norms__, double *errs_n__, double *errs_c__, dcomplex *res, 
                                     double *ayb, dcomplex *dy, dcomplex *y_tail__, double *rcond, int *
                                     ithresh, double *rthresh, double *dz_ub__, int *ignore_cwise__, 
                                     int *info);

extern int C2F(cla_heamv_)(int *uplo, int *n, double *alpha, 
                           dcomplex *a, int *lda, dcomplex *x, int *incx, double *beta, double 
                           *y, int *incy);

double C2F(cla_hercond_c_)(char *uplo, int *n, dcomplex *a, int *lda, 
                           dcomplex *af, int *ldaf, int *ipiv, double *c__, int *capply,
                           int *info, dcomplex *work, double *rwork, int uplo_len);

double C2F(cla_hercond_x_)(char *uplo, int *n, dcomplex *a, int *lda, 
                           dcomplex *af, int *ldaf, int *ipiv, dcomplex *x, int *info, 
                           dcomplex *work, double *rwork, int uplo_len);

extern int C2F(cla_herfsx_extended_)(int *prec_type__, char *uplo, 
                                     int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *af, 
                                     int *ldaf, int *ipiv, int *colequ, double *c__, dcomplex *b, 
                                     int *ldb, dcomplex *y, int *ldy, double *berr_out__, int *
                                     n_norms__, double *errs_n__, double *errs_c__, dcomplex *res, double *ayb, 
                                     dcomplex *dy, dcomplex *y_tail__, double *rcond, int *ithresh, double *
                                     rthresh, double *dz_ub__, int *ignore_cwise__, int *info, 
                                     int uplo_len);

double C2F(cla_herpvgrw_)(char *uplo, int *n, int *info, dcomplex *a, 
                          int *lda, dcomplex *af, int *ldaf, int *ipiv, double *work, 
                          int uplo_len);

extern int C2F(cla_lin_berr_)(int *n, int *nz, int *nrhs, 
                              dcomplex *res, double *ayb, double *berr);

double C2F(cla_porcond_c_)(char *uplo, int *n, dcomplex *a, int *lda, 
                           dcomplex *af, int *ldaf, double *c__, int *capply, int *info,
                           dcomplex *work, double *rwork, int uplo_len);

double C2F(cla_porcond_x_)(char *uplo, int *n, dcomplex *a, int *lda, 
                           dcomplex *af, int *ldaf, dcomplex *x, int *info, dcomplex *work, 
                           double *rwork, int uplo_len);

extern int C2F(cla_porfsx_extended_)(int *prec_type__, char *uplo, 
                                     int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *af, 
                                     int *ldaf, int *colequ, double *c__, dcomplex *b, int *ldb, 
                                     dcomplex *y, int *ldy, double *berr_out__, int *n_norms__, double *
                                     errs_n__, double *errs_c__, dcomplex *res, double *ayb, dcomplex *dy, 
                                     dcomplex *y_tail__, double *rcond, int *ithresh, double *rthresh, double 
                                     *dz_ub__, int *ignore_cwise__, int *info, int uplo_len);

double C2F(cla_porpvgrw_)(char *uplo, int *ncols, dcomplex *a, int *
                          lda, dcomplex *af, int *ldaf, double *work, int uplo_len);

double C2F(cla_rpvgrw_)(int *n, int *ncols, dcomplex *a, int *lda, 
                        dcomplex *af, int *ldaf);

extern int C2F(cla_syamv_)(int *uplo, int *n, double *alpha, 
                           dcomplex *a, int *lda, dcomplex *x, int *incx, double *beta, double 
                           *y, int *incy);

double C2F(cla_syrcond_c_)(char *uplo, int *n, dcomplex *a, int *lda, 
                           dcomplex *af, int *ldaf, int *ipiv, double *c__, int *capply,
                           int *info, dcomplex *work, double *rwork, int uplo_len);

double C2F(cla_syrcond_x_)(char *uplo, int *n, dcomplex *a, int *lda, 
                           dcomplex *af, int *ldaf, int *ipiv, dcomplex *x, int *info, 
                           dcomplex *work, double *rwork, int uplo_len);

extern int C2F(cla_syrfsx_extended_)(int *prec_type__, char *uplo, 
                                     int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *af, 
                                     int *ldaf, int *ipiv, int *colequ, double *c__, dcomplex *b, 
                                     int *ldb, dcomplex *y, int *ldy, double *berr_out__, int *
                                     n_norms__, double *errs_n__, double *errs_c__, dcomplex *res, double *ayb, 
                                     dcomplex *dy, dcomplex *y_tail__, double *rcond, int *ithresh, double *
                                     rthresh, double *dz_ub__, int *ignore_cwise__, int *info, 
                                     int uplo_len);

double C2F(cla_syrpvgrw_)(char *uplo, int *n, int *info, dcomplex *a, 
                          int *lda, dcomplex *af, int *ldaf, int *ipiv, double *work, 
                          int uplo_len);

extern int C2F(cla_wwaddw_)(int *n, dcomplex *x, dcomplex *y, dcomplex 
                            *w);

extern int C2F(clabrd)(int *m, int *n, int *nb, dcomplex *a, 
                       int *lda, double *d__, double *e, dcomplex *tauq, dcomplex *taup, 
                       dcomplex *x, int *ldx, dcomplex *y, int *ldy);

extern int C2F(clacgv)(int *n, dcomplex *x, int *incx);

extern int C2F(clacn2)(int *n, dcomplex *v, dcomplex *x, double *est, 
                       int *kase, int *isave);

extern int C2F(clacon)(int *n, dcomplex *v, dcomplex *x, double *est, 
                       int *kase);

extern int C2F(clacp2)(char *uplo, int *m, int *n, double *a, 
                       int *lda, dcomplex *b, int *ldb);

extern int C2F(clacpy)(char *uplo, int *m, int *n, dcomplex *a, 
                       int *lda, dcomplex *b, int *ldb);

extern int C2F(clacrm)(int *m, int *n, dcomplex *a, int *lda, 
                       double *b, int *ldb, dcomplex *c__, int *ldc, double *rwork);

extern int C2F(clacrt)(int *n, dcomplex *cx, int *incx, dcomplex *
                       cy, int *incy, dcomplex *c__, dcomplex *s);

/* dcomplex */ void C2F(cladiv)(dcomplex * ret_val, dcomplex *x, dcomplex *y);

extern int C2F(claed0)(int *qsiz, int *n, double *d__, double *e, 
                       dcomplex *q, int *ldq, dcomplex *qstore, int *ldqs, double *rwork, 
                       int *iwork, int *info);

extern int C2F(claed7)(int *n, int *cutpnt, int *qsiz, 
                       int *tlvls, int *curlvl, int *curpbm, double *d__, dcomplex *
                       q, int *ldq, double *rho, int *indxq, double *qstore, int *
                       qptr, int *prmptr, int *perm, int *givptr, int *
                       givcol, double *givnum, dcomplex *work, double *rwork, int *iwork, 
                       int *info);

extern int C2F(claed8)(int *k, int *n, int *qsiz, dcomplex *
                       q, int *ldq, double *d__, double *rho, int *cutpnt, double *z__, 
                       double *dlamda, dcomplex *q2, int *ldq2, double *w, int *indxp, 
                       int *indx, int *indxq, int *perm, int *givptr, 
                       int *givcol, double *givnum, int *info);

extern int C2F(claein)(int *rightv, int *noinit, int *n, 
                       dcomplex *h__, int *ldh, dcomplex *w, dcomplex *v, dcomplex *b, 
                       int *ldb, double *rwork, double *eps3, double *smlnum, int *info);

extern int C2F(claesy)(dcomplex *a, dcomplex *b, dcomplex *c__, dcomplex *
                       rt1, dcomplex *rt2, dcomplex *evscal, dcomplex *cs1, dcomplex *sn1);

extern int C2F(claev2)(dcomplex *a, dcomplex *b, dcomplex *c__, double *rt1, 
                       double *rt2, double *cs1, dcomplex *sn1);

extern int C2F(clag2z)(int *m, int *n, dcomplex *sa, int *
                       ldsa, dcomplex *a, int *lda, int *info);

extern int C2F(clags2)(int *upper, double *a1, dcomplex *a2, double *a3, 
                       double *b1, dcomplex *b2, double *b3, double *csu, dcomplex *snu, double *csv, 
                       dcomplex *snv, double *csq, dcomplex *snq);

extern int C2F(clagtm)(char *trans, int *n, int *nrhs, double *
                       alpha, dcomplex *dl, dcomplex *d__, dcomplex *du, dcomplex *x, int *
                       ldx, double *beta, dcomplex *b, int *ldb);

extern int C2F(clahef)(char *uplo, int *n, int *nb, int *kb, 
                       dcomplex *a, int *lda, int *ipiv, dcomplex *w, int *ldw, 
                       int *info);

extern int C2F(clahqr)(int *wantt, int *wantz, int *n, 
                       int *ilo, int *ihi, dcomplex *h__, int *ldh, dcomplex *w, 
                       int *iloz, int *ihiz, dcomplex *z__, int *ldz, int *
                       info);

extern int C2F(clahr2)(int *n, int *k, int *nb, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *t, int *ldt, dcomplex *y, 
                       int *ldy);

extern int C2F(clahrd)(int *n, int *k, int *nb, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *t, int *ldt, dcomplex *y, 
                       int *ldy);

extern int C2F(claic1)(int *job, int *j, dcomplex *x, double *sest, 
                       dcomplex *w, dcomplex *gamma, double *sestpr, dcomplex *s, dcomplex *c__);

extern int C2F(clals0)(int *icompq, int *nl, int *nr, 
                       int *sqre, int *nrhs, dcomplex *b, int *ldb, dcomplex *bx, 
                       int *ldbx, int *perm, int *givptr, int *givcol, 
                       int *ldgcol, double *givnum, int *ldgnum, double *poles, double *
                       difl, double *difr, double *z__, int *k, double *c__, double *s, double *
                       rwork, int *info);

extern int C2F(clalsa)(int *icompq, int *smlsiz, int *n, 
                       int *nrhs, dcomplex *b, int *ldb, dcomplex *bx, int *ldbx, 
                       double *u, int *ldu, double *vt, int *k, double *difl, double *difr, 
                       double *z__, double *poles, int *givptr, int *givcol, int *
                       ldgcol, int *perm, double *givnum, double *c__, double *s, double *rwork, 
                       int *iwork, int *info);

extern int C2F(clalsd)(char *uplo, int *smlsiz, int *n, int 
                       *nrhs, double *d__, double *e, dcomplex *b, int *ldb, double *rcond, 
                       int *rank, dcomplex *work, double *rwork, int *iwork, int *
                       info);

double C2F(clangb)(char *norm, int *n, int *kl, int *ku, dcomplex *
                   ab, int *ldab, double *work);

double C2F(clange)(char *norm, int *m, int *n, dcomplex *a, int *
                   lda, double *work);

double C2F(clangt)(char *norm, int *n, dcomplex *dl, dcomplex *d__, dcomplex 
                   *du);

double C2F(clanhb)(char *norm, char *uplo, int *n, int *k, dcomplex *
                   ab, int *ldab, double *work);

double C2F(clanhe)(char *norm, char *uplo, int *n, dcomplex *a, int *
                   lda, double *work);

double C2F(clanhf)(char *norm, char *transr, char *uplo, int *n, dcomplex *
                   a, double *work);

double C2F(clanhp)(char *norm, char *uplo, int *n, dcomplex *ap, double *
                   work);

double C2F(clanhs)(char *norm, int *n, dcomplex *a, int *lda, double *
                   work);

double C2F(clanht)(char *norm, int *n, double *d__, dcomplex *e);

double C2F(clansb)(char *norm, char *uplo, int *n, int *k, dcomplex *
                   ab, int *ldab, double *work);

double C2F(clansp)(char *norm, char *uplo, int *n, dcomplex *ap, double *
                   work);

double C2F(clansy)(char *norm, char *uplo, int *n, dcomplex *a, int *
                   lda, double *work);

double C2F(clantb)(char *norm, char *uplo, char *diag, int *n, int *k, 
                   dcomplex *ab, int *ldab, double *work);

double C2F(clantp)(char *norm, char *uplo, char *diag, int *n, dcomplex *
                   ap, double *work);

double C2F(clantr)(char *norm, char *uplo, char *diag, int *m, int *n, 
                   dcomplex *a, int *lda, double *work);

extern int C2F(clapll)(int *n, dcomplex *x, int *incx, dcomplex *
                       y, int *incy, double *ssmin);

extern int C2F(clapmt)(int *forwrd, int *m, int *n, dcomplex 
                       *x, int *ldx, int *k);

extern int C2F(claqgb)(int *m, int *n, int *kl, int *ku, 
                       dcomplex *ab, int *ldab, double *r__, double *c__, double *rowcnd, double 
                       *colcnd, double *amax, char *equed);

extern int C2F(claqge)(int *m, int *n, dcomplex *a, int *lda, 
                       double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, char *
                       equed);

extern int C2F(claqhb)(char *uplo, int *n, int *kd, dcomplex *ab, 
                       int *ldab, double *s, double *scond, double *amax, char *equed);

extern int C2F(claqhe)(char *uplo, int *n, dcomplex *a, int *lda, 
                       double *s, double *scond, double *amax, char *equed);

extern int C2F(claqhp)(char *uplo, int *n, dcomplex *ap, double *s, 
                       double *scond, double *amax, char *equed);

extern int C2F(claqp2)(int *m, int *n, int *offset, dcomplex 
                       *a, int *lda, int *jpvt, dcomplex *tau, double *vn1, double *vn2, 
                       dcomplex *work);

extern int C2F(claqps)(int *m, int *n, int *offset, int 
                       *nb, int *kb, dcomplex *a, int *lda, int *jpvt, dcomplex *
                       tau, double *vn1, double *vn2, dcomplex *auxv, dcomplex *f, int *ldf);

extern int C2F(claqr0)(int *wantt, int *wantz, int *n, 
                       int *ilo, int *ihi, dcomplex *h__, int *ldh, dcomplex *w, 
                       int *iloz, int *ihiz, dcomplex *z__, int *ldz, dcomplex *
                       work, int *lwork, int *info);

extern int C2F(claqr1)(int *n, dcomplex *h__, int *ldh, dcomplex *
                       s1, dcomplex *s2, dcomplex *v);

extern int C2F(claqr2)(int *wantt, int *wantz, int *n, 
                       int *ktop, int *kbot, int *nw, dcomplex *h__, int *ldh, 
                       int *iloz, int *ihiz, dcomplex *z__, int *ldz, int *
                       ns, int *nd, dcomplex *sh, dcomplex *v, int *ldv, int *nh, 
                       dcomplex *t, int *ldt, int *nv, dcomplex *wv, int *ldwv, 
                       dcomplex *work, int *lwork);

extern int C2F(claqr3)(int *wantt, int *wantz, int *n, 
                       int *ktop, int *kbot, int *nw, dcomplex *h__, int *ldh, 
                       int *iloz, int *ihiz, dcomplex *z__, int *ldz, int *
                       ns, int *nd, dcomplex *sh, dcomplex *v, int *ldv, int *nh, 
                       dcomplex *t, int *ldt, int *nv, dcomplex *wv, int *ldwv, 
                       dcomplex *work, int *lwork);

extern int C2F(claqr4)(int *wantt, int *wantz, int *n, 
                       int *ilo, int *ihi, dcomplex *h__, int *ldh, dcomplex *w, 
                       int *iloz, int *ihiz, dcomplex *z__, int *ldz, dcomplex *
                       work, int *lwork, int *info);

extern int C2F(claqr5)(int *wantt, int *wantz, int *kacc22, 
                       int *n, int *ktop, int *kbot, int *nshfts, dcomplex *s, 
                       dcomplex *h__, int *ldh, int *iloz, int *ihiz, dcomplex *
                       z__, int *ldz, dcomplex *v, int *ldv, dcomplex *u, int *ldu, 
                       int *nv, dcomplex *wv, int *ldwv, int *nh, dcomplex *wh, 
                       int *ldwh);

extern int C2F(claqsb)(char *uplo, int *n, int *kd, dcomplex *ab, 
                       int *ldab, double *s, double *scond, double *amax, char *equed);

extern int C2F(claqsp)(char *uplo, int *n, dcomplex *ap, double *s, 
                       double *scond, double *amax, char *equed);

extern int C2F(claqsy)(char *uplo, int *n, dcomplex *a, int *lda, 
                       double *s, double *scond, double *amax, char *equed);

extern int C2F(clar1v)(int *n, int *b1, int *bn, double *
                       lambda, double *d__, double *l, double *ld, double *lld, double *pivmin, double *
                       gaptol, dcomplex *z__, int *wantnc, int *negcnt, double *ztz, 
                       double *mingma, int *r__, int *isuppz, double *nrminv, double *
                       resid, double *rqcorr, double *work);

extern int C2F(clar2v)(int *n, dcomplex *x, dcomplex *y, dcomplex *z__, 
                       int *incx, double *c__, dcomplex *s, int *incc);

extern int C2F(clarcm)(int *m, int *n, double *a, int *lda, 
                       dcomplex *b, int *ldb, dcomplex *c__, int *ldc, double *rwork);

extern int C2F(clarf)(char *side, int *m, int *n, dcomplex *v, 
                      int *incv, dcomplex *tau, dcomplex *c__, int *ldc, dcomplex *
                      work);

extern int C2F(clarfb)(char *side, char *trans, char *direct, char *
                       storev, int *m, int *n, int *k, dcomplex *v, int *ldv, 
                       dcomplex *t, int *ldt, dcomplex *c__, int *ldc, dcomplex *work, 
                       int *ldwork);

extern int C2F(clarfg)(int *n, dcomplex *alpha, dcomplex *x, int *
                       incx, dcomplex *tau);

extern int C2F(clarfp)(int *n, dcomplex *alpha, dcomplex *x, int *
                       incx, dcomplex *tau);

extern int C2F(clarft)(char *direct, char *storev, int *n, int *
                       k, dcomplex *v, int *ldv, dcomplex *tau, dcomplex *t, int *ldt);

extern int C2F(clarfx)(char *side, int *m, int *n, dcomplex *v, 
                       dcomplex *tau, dcomplex *c__, int *ldc, dcomplex *work);

extern int C2F(clargv)(int *n, dcomplex *x, int *incx, dcomplex *
                       y, int *incy, double *c__, int *incc);

extern int C2F(clarnv)(int *idist, int *iseed, int *n, 
                       dcomplex *x);

extern int C2F(clarrv)(int *n, double *vl, double *vu, double *d__, double *
                       l, double *pivmin, int *isplit, int *m, int *dol, int *
                       dou, double *minrgp, double *rtol1, double *rtol2, double *w, double *werr, 
                       double *wgap, int *iblock, int *indexw, double *gers, dcomplex *
                       z__, int *ldz, int *isuppz, double *work, int *iwork, 
                       int *info);

extern int C2F(clarscl2)(int *m, int *n, double *d__, dcomplex *x, 
                         int *ldx);

extern int C2F(clartg)(dcomplex *f, dcomplex *g, double *cs, dcomplex *sn, 
                       dcomplex *r__);

extern int C2F(clartv)(int *n, dcomplex *x, int *incx, dcomplex *
                       y, int *incy, double *c__, dcomplex *s, int *incc);

extern int C2F(clarz)(char *side, int *m, int *n, int *l, 
                      dcomplex *v, int *incv, dcomplex *tau, dcomplex *c__, int *ldc, 
                      dcomplex *work);

extern int C2F(clarzb)(char *side, char *trans, char *direct, char *
                       storev, int *m, int *n, int *k, int *l, dcomplex *v, 
                       int *ldv, dcomplex *t, int *ldt, dcomplex *c__, int *ldc, 
                       dcomplex *work, int *ldwork);

extern int C2F(clarzt)(char *direct, char *storev, int *n, int *
                       k, dcomplex *v, int *ldv, dcomplex *tau, dcomplex *t, int *ldt);

extern int C2F(clascl)(char *type__, int *kl, int *ku, double *
                       cfrom, double *cto, int *m, int *n, dcomplex *a, int *lda, 
                       int *info);

extern int C2F(clascl2)(int *m, int *n, double *d__, dcomplex *x, 
                        int *ldx);

extern int C2F(claset)(char *uplo, int *m, int *n, dcomplex *
                       alpha, dcomplex *beta, dcomplex *a, int *lda);

extern int C2F(clasr)(char *side, char *pivot, char *direct, int *m, 
                      int *n, double *c__, double *s, dcomplex *a, int *lda);

extern int C2F(classq)(int *n, dcomplex *x, int *incx, double *
                       scale, double *sumsq);

extern int C2F(claswp)(int *n, dcomplex *a, int *lda, int *
                       k1, int *k2, int *ipiv, int *incx);

extern int C2F(clasyf)(char *uplo, int *n, int *nb, int *kb, 
                       dcomplex *a, int *lda, int *ipiv, dcomplex *w, int *ldw, 
                       int *info);

extern int C2F(clatbs)(char *uplo, char *trans, char *diag, char *
                       normin, int *n, int *kd, dcomplex *ab, int *ldab, dcomplex *
                       x, double *scale, double *cnorm, int *info);

extern int C2F(clatdf)(int *ijob, int *n, dcomplex *z__, int 
                       *ldz, dcomplex *rhs, double *rdsum, double *rdscal, int *ipiv, int 
                       *jpiv);

extern int C2F(clatps)(char *uplo, char *trans, char *diag, char *
                       normin, int *n, dcomplex *ap, dcomplex *x, double *scale, double *cnorm, 
                       int *info);

extern int C2F(clatrd)(char *uplo, int *n, int *nb, dcomplex *a, 
                       int *lda, double *e, dcomplex *tau, dcomplex *w, int *ldw);

extern int C2F(clatrs)(char *uplo, char *trans, char *diag, char *
                       normin, int *n, dcomplex *a, int *lda, dcomplex *x, double *scale, 
                       double *cnorm, int *info);

extern int C2F(clatrz)(int *m, int *n, int *l, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work);

extern int C2F(clatzm)(char *side, int *m, int *n, dcomplex *v, 
                       int *incv, dcomplex *tau, dcomplex *c1, dcomplex *c2, int *ldc, 
                       dcomplex *work);

extern int C2F(clauu2)(char *uplo, int *n, dcomplex *a, int *lda, 
                       int *info);

extern int C2F(clauum)(char *uplo, int *n, dcomplex *a, int *lda, 
                       int *info);

extern int C2F(cpbcon)(char *uplo, int *n, int *kd, dcomplex *ab, 
                       int *ldab, double *anorm, double *rcond, dcomplex *work, double *rwork, 
                       int *info);

extern int C2F(cpbequ)(char *uplo, int *n, int *kd, dcomplex *ab, 
                       int *ldab, double *s, double *scond, double *amax, int *info);

extern int C2F(cpbrfs)(char *uplo, int *n, int *kd, int *
                       nrhs, dcomplex *ab, int *ldab, dcomplex *afb, int *ldafb, 
                       dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *ferr, double *
                       berr, dcomplex *work, double *rwork, int *info);

extern int C2F(cpbstf)(char *uplo, int *n, int *kd, dcomplex *ab, 
                       int *ldab, int *info);

extern int C2F(cpbsv)(char *uplo, int *n, int *kd, int *
                      nrhs, dcomplex *ab, int *ldab, dcomplex *b, int *ldb, int *
                      info);

extern int C2F(cpbsvx)(char *fact, char *uplo, int *n, int *kd, 
                       int *nrhs, dcomplex *ab, int *ldab, dcomplex *afb, int *
                       ldafb, char *equed, double *s, dcomplex *b, int *ldb, dcomplex *x, 
                       int *ldx, double *rcond, double *ferr, double *berr, dcomplex *work, 
                       double *rwork, int *info);

extern int C2F(cpbtf2)(char *uplo, int *n, int *kd, dcomplex *ab, 
                       int *ldab, int *info);

extern int C2F(cpbtrf)(char *uplo, int *n, int *kd, dcomplex *ab, 
                       int *ldab, int *info);

extern int C2F(cpbtrs)(char *uplo, int *n, int *kd, int *
                       nrhs, dcomplex *ab, int *ldab, dcomplex *b, int *ldb, int *
                       info);

extern int C2F(cpftrf)(char *transr, char *uplo, int *n, dcomplex *a, 
                       int *info);

extern int C2F(cpftri)(char *transr, char *uplo, int *n, dcomplex *a, 
                       int *info);

extern int C2F(cpftrs)(char *transr, char *uplo, int *n, int *
                       nrhs, dcomplex *a, dcomplex *b, int *ldb, int *info);

extern int C2F(cpocon)(char *uplo, int *n, dcomplex *a, int *lda, 
                       double *anorm, double *rcond, dcomplex *work, double *rwork, int *info);

extern int C2F(cpoequ)(int *n, dcomplex *a, int *lda, double *s, 
                       double *scond, double *amax, int *info);

extern int C2F(cpoequb)(int *n, dcomplex *a, int *lda, double *s, 
                        double *scond, double *amax, int *info);

extern int C2F(cporfs)(char *uplo, int *n, int *nrhs, dcomplex *
                       a, int *lda, dcomplex *af, int *ldaf, dcomplex *b, int *ldb, 
                       dcomplex *x, int *ldx, double *ferr, double *berr, dcomplex *work, 
                       double *rwork, int *info);

extern int C2F(cporfsx)(char *uplo, char *equed, int *n, int *
                        nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, double *s, 
                        dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double 
                        *berr, int *n_err_bnds__, double *err_bnds_norm__, double *
                        err_bnds_comp__, int *nparams, double *params, dcomplex *work, double *
                        rwork, int *info);

extern int C2F(cposv)(char *uplo, int *n, int *nrhs, dcomplex *a, 
                      int *lda, dcomplex *b, int *ldb, int *info);

extern int C2F(cposvx)(char *fact, char *uplo, int *n, int *
                       nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, char *
                       equed, double *s, dcomplex *b, int *ldb, dcomplex *x, int *ldx, 
                       double *rcond, double *ferr, double *berr, dcomplex *work, double *rwork, 
                       int *info);

extern int C2F(cposvxx)(char *fact, char *uplo, int *n, int *
                        nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, char *
                        equed, double *s, dcomplex *b, int *ldb, dcomplex *x, int *ldx, 
                        double *rcond, double *rpvgrw, double *berr, int *n_err_bnds__, double *
                        err_bnds_norm__, double *err_bnds_comp__, int *nparams, double *
                        params, dcomplex *work, double *rwork, int *info);

extern int C2F(cpotf2)(char *uplo, int *n, dcomplex *a, int *lda, 
                       int *info);

extern int C2F(cpotrf)(char *uplo, int *n, dcomplex *a, int *lda, 
                       int *info);

extern int C2F(cpotri)(char *uplo, int *n, dcomplex *a, int *lda, 
                       int *info);

extern int C2F(cpotrs)(char *uplo, int *n, int *nrhs, dcomplex *
                       a, int *lda, dcomplex *b, int *ldb, int *info);

extern int C2F(cppcon)(char *uplo, int *n, dcomplex *ap, double *anorm, 
                       double *rcond, dcomplex *work, double *rwork, int *info);

extern int C2F(cppequ)(char *uplo, int *n, dcomplex *ap, double *s, 
                       double *scond, double *amax, int *info);

extern int C2F(cpprfs)(char *uplo, int *n, int *nrhs, dcomplex *
                       ap, dcomplex *afp, dcomplex *b, int *ldb, dcomplex *x, int *ldx, 
                       double *ferr, double *berr, dcomplex *work, double *rwork, int *info);

extern int C2F(cppsv)(char *uplo, int *n, int *nrhs, dcomplex *
                      ap, dcomplex *b, int *ldb, int *info);

extern int C2F(cppsvx)(char *fact, char *uplo, int *n, int *
                       nrhs, dcomplex *ap, dcomplex *afp, char *equed, double *s, dcomplex *b, 
                       int *ldb, dcomplex *x, int *ldx, double *rcond, double *ferr, double 
                       *berr, dcomplex *work, double *rwork, int *info);

extern int C2F(cpptrf)(char *uplo, int *n, dcomplex *ap, int *
                       info);

extern int C2F(cpptri)(char *uplo, int *n, dcomplex *ap, int *
                       info);

extern int C2F(cpptrs)(char *uplo, int *n, int *nrhs, dcomplex *
                       ap, dcomplex *b, int *ldb, int *info);

extern int C2F(cpstf2)(char *uplo, int *n, dcomplex *a, int *lda, 
                       int *piv, int *rank, double *tol, double *work, int *info);

extern int C2F(cpstrf)(char *uplo, int *n, dcomplex *a, int *lda, 
                       int *piv, int *rank, double *tol, double *work, int *info);

extern int C2F(cptcon)(int *n, double *d__, dcomplex *e, double *anorm, 
                       double *rcond, double *rwork, int *info);

extern int C2F(cpteqr)(char *compz, int *n, double *d__, double *e, 
                       dcomplex *z__, int *ldz, double *work, int *info);

extern int C2F(cptrfs)(char *uplo, int *n, int *nrhs, double *d__, 
                       dcomplex *e, double *df, dcomplex *ef, dcomplex *b, int *ldb, dcomplex 
                       *x, int *ldx, double *ferr, double *berr, dcomplex *work, double *rwork, 
                       int *info);

extern int C2F(cptsv)(int *n, int *nrhs, double *d__, dcomplex *e, 
                      dcomplex *b, int *ldb, int *info);

extern int C2F(cptsvx)(char *fact, int *n, int *nrhs, double *d__, 
                       dcomplex *e, double *df, dcomplex *ef, dcomplex *b, int *ldb, dcomplex 
                       *x, int *ldx, double *rcond, double *ferr, double *berr, dcomplex *work, 
                       double *rwork, int *info);

extern int C2F(cpttrf)(int *n, double *d__, dcomplex *e, int *info);

extern int C2F(cpttrs)(char *uplo, int *n, int *nrhs, double *d__, 
                       dcomplex *e, dcomplex *b, int *ldb, int *info);

extern int C2F(cptts2)(int *iuplo, int *n, int *nrhs, double *
                       d__, dcomplex *e, dcomplex *b, int *ldb);

extern int C2F(crot)(int *n, dcomplex *cx, int *incx, dcomplex *
                     cy, int *incy, double *c__, dcomplex *s);

extern int C2F(cspcon)(char *uplo, int *n, dcomplex *ap, int *
                       ipiv, double *anorm, double *rcond, dcomplex *work, int *info);

extern int C2F(cspmv)(char *uplo, int *n, dcomplex *alpha, dcomplex *
                      ap, dcomplex *x, int *incx, dcomplex *beta, dcomplex *y, int *
                      incy);

extern int C2F(cspr)(char *uplo, int *n, dcomplex *alpha, dcomplex *x, 
                     int *incx, dcomplex *ap);

extern int C2F(csprfs)(char *uplo, int *n, int *nrhs, dcomplex *
                       ap, dcomplex *afp, int *ipiv, dcomplex *b, int *ldb, dcomplex *x, 
                       int *ldx, double *ferr, double *berr, dcomplex *work, double *rwork, 
                       int *info);

extern int C2F(cspsv)(char *uplo, int *n, int *nrhs, dcomplex *
                      ap, int *ipiv, dcomplex *b, int *ldb, int *info);

extern int C2F(cspsvx)(char *fact, char *uplo, int *n, int *
                       nrhs, dcomplex *ap, dcomplex *afp, int *ipiv, dcomplex *b, int *
                       ldb, dcomplex *x, int *ldx, double *rcond, double *ferr, double *berr, 
                       dcomplex *work, double *rwork, int *info);

extern int C2F(csptrf)(char *uplo, int *n, dcomplex *ap, int *
                       ipiv, int *info);

extern int C2F(csptri)(char *uplo, int *n, dcomplex *ap, int *
                       ipiv, dcomplex *work, int *info);

extern int C2F(csptrs)(char *uplo, int *n, int *nrhs, dcomplex *
                       ap, int *ipiv, dcomplex *b, int *ldb, int *info);

extern int C2F(csrscl)(int *n, double *sa, dcomplex *sx, int *incx);

extern int C2F(cstedc)(char *compz, int *n, double *d__, double *e, 
                       dcomplex *z__, int *ldz, dcomplex *work, int *lwork, double *
                       rwork, int *lrwork, int *iwork, int *liwork, int *
                       info);

extern int C2F(cstegr)(char *jobz, char *range, int *n, double *d__, 
                       double *e, double *vl, double *vu, int *il, int *iu, double *abstol, 
                       int *m, double *w, dcomplex *z__, int *ldz, int *isuppz, 
                       double *work, int *lwork, int *iwork, int *liwork, int *
                       info);

extern int C2F(cstein)(int *n, double *d__, double *e, int *m, double 
                       *w, int *iblock, int *isplit, dcomplex *z__, int *ldz, 
                       double *work, int *iwork, int *ifail, int *info);

extern int C2F(cstemr)(char *jobz, char *range, int *n, double *d__, 
                       double *e, double *vl, double *vu, int *il, int *iu, int *m, 
                       double *w, dcomplex *z__, int *ldz, int *nzc, int *isuppz, 
                       int *tryrac, double *work, int *lwork, int *iwork, int *
                       liwork, int *info);

extern int C2F(csteqr)(char *compz, int *n, double *d__, double *e, 
                       dcomplex *z__, int *ldz, double *work, int *info);

extern int C2F(csycon)(char *uplo, int *n, dcomplex *a, int *lda, 
                       int *ipiv, double *anorm, double *rcond, dcomplex *work, int *
                       info);

extern int C2F(csyequb)(char *uplo, int *n, dcomplex *a, int *
                        lda, double *s, double *scond, double *amax, dcomplex *work, int *info);

extern int C2F(csymv)(char *uplo, int *n, dcomplex *alpha, dcomplex *
                      a, int *lda, dcomplex *x, int *incx, dcomplex *beta, dcomplex *y, 
                      int *incy);

extern int C2F(csyr)(char *uplo, int *n, dcomplex *alpha, dcomplex *x, 
                     int *incx, dcomplex *a, int *lda);

extern int C2F(csyrfs)(char *uplo, int *n, int *nrhs, dcomplex *
                       a, int *lda, dcomplex *af, int *ldaf, int *ipiv, dcomplex *
                       b, int *ldb, dcomplex *x, int *ldx, double *ferr, double *berr, 
                       dcomplex *work, double *rwork, int *info);

extern int C2F(csyrfsx)(char *uplo, char *equed, int *n, int *
                        nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, int *
                        ipiv, double *s, dcomplex *b, int *ldb, dcomplex *x, int *ldx, 
                        double *rcond, double *berr, int *n_err_bnds__, double *err_bnds_norm__, 
                        double *err_bnds_comp__, int *nparams, double *params, dcomplex *work, 
                        double *rwork, int *info);

extern int C2F(csysv)(char *uplo, int *n, int *nrhs, dcomplex *a, 
                      int *lda, int *ipiv, dcomplex *b, int *ldb, dcomplex *work, 
                      int *lwork, int *info);

extern int C2F(csysvx)(char *fact, char *uplo, int *n, int *
                       nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, int *
                       ipiv, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, 
                       double *ferr, double *berr, dcomplex *work, int *lwork, double *rwork, 
                       int *info);

extern int C2F(csysvxx)(char *fact, char *uplo, int *n, int *
                        nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, int *
                        ipiv, char *equed, double *s, dcomplex *b, int *ldb, dcomplex *x, 
                        int *ldx, double *rcond, double *rpvgrw, double *berr, int *
                        n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__, int *
                        nparams, double *params, dcomplex *work, double *rwork, int *info);

extern int C2F(csytf2)(char *uplo, int *n, dcomplex *a, int *lda, 
                       int *ipiv, int *info);

extern int C2F(csytrf)(char *uplo, int *n, dcomplex *a, int *lda, 
                       int *ipiv, dcomplex *work, int *lwork, int *info);

extern int C2F(csytri)(char *uplo, int *n, dcomplex *a, int *lda, 
                       int *ipiv, dcomplex *work, int *info);

extern int C2F(csytrs)(char *uplo, int *n, int *nrhs, dcomplex *
                       a, int *lda, int *ipiv, dcomplex *b, int *ldb, int *
                       info);

extern int C2F(ctbcon)(char *norm, char *uplo, char *diag, int *n, 
                       int *kd, dcomplex *ab, int *ldab, double *rcond, dcomplex *work, 
                       double *rwork, int *info);

extern int C2F(ctbrfs)(char *uplo, char *trans, char *diag, int *n, 
                       int *kd, int *nrhs, dcomplex *ab, int *ldab, dcomplex *b, 
                       int *ldb, dcomplex *x, int *ldx, double *ferr, double *berr, 
                       dcomplex *work, double *rwork, int *info);

extern int C2F(ctbtrs)(char *uplo, char *trans, char *diag, int *n, 
                       int *kd, int *nrhs, dcomplex *ab, int *ldab, dcomplex *b, 
                       int *ldb, int *info);

extern int C2F(ctfsm)(char *transr, char *side, char *uplo, char *trans, 
                      char *diag, int *m, int *n, dcomplex *alpha, dcomplex *a, 
                      dcomplex *b, int *ldb);

extern int C2F(ctftri)(char *transr, char *uplo, char *diag, int *n, 
                       dcomplex *a, int *info);

extern int C2F(ctfttp)(char *transr, char *uplo, int *n, dcomplex *
                       arf, dcomplex *ap, int *info);

extern int C2F(ctfttr)(char *transr, char *uplo, int *n, dcomplex *
                       arf, dcomplex *a, int *lda, int *info);

extern int C2F(ctgevc)(char *side, char *howmny, int *select, 
                       int *n, dcomplex *s, int *lds, dcomplex *p, int *ldp, 
                       dcomplex *vl, int *ldvl, dcomplex *vr, int *ldvr, int *mm, 
                       int *m, dcomplex *work, double *rwork, int *info);

extern int C2F(ctgex2)(int *wantq, int *wantz, int *n, 
                       dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *q, 
                       int *ldq, dcomplex *z__, int *ldz, int *j1, int *info);

extern int C2F(ctgexc)(int *wantq, int *wantz, int *n, 
                       dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *q, 
                       int *ldq, dcomplex *z__, int *ldz, int *ifst, int *
                       ilst, int *info);

extern int C2F(ctgsen)(int *ijob, int *wantq, int *wantz, 
                       int *select, int *n, dcomplex *a, int *lda, dcomplex *b, 
                       int *ldb, dcomplex *alpha, dcomplex *beta, dcomplex *q, int *ldq, 
                       dcomplex *z__, int *ldz, int *m, double *pl, double *pr, double *
                       dif, dcomplex *work, int *lwork, int *iwork, int *liwork, 
                       int *info);

extern int C2F(ctgsja)(char *jobu, char *jobv, char *jobq, int *m, 
                       int *p, int *n, int *k, int *l, dcomplex *a, int *
                       lda, dcomplex *b, int *ldb, double *tola, double *tolb, double *alpha, 
                       double *beta, dcomplex *u, int *ldu, dcomplex *v, int *ldv, 
                       dcomplex *q, int *ldq, dcomplex *work, int *ncycle, int *
                       info);

extern int C2F(ctgsna)(char *job, char *howmny, int *select, 
                       int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       dcomplex *vl, int *ldvl, dcomplex *vr, int *ldvr, double *s, double 
                       *dif, int *mm, int *m, dcomplex *work, int *lwork, int 
                       *iwork, int *info);

extern int C2F(ctgsy2)(char *trans, int *ijob, int *m, int *
                       n, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *c__, 
                       int *ldc, dcomplex *d__, int *ldd, dcomplex *e, int *lde, 
                       dcomplex *f, int *ldf, double *scale, double *rdsum, double *rdscal, 
                       int *info);

extern int C2F(ctgsyl)(char *trans, int *ijob, int *m, int *
                       n, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *c__, 
                       int *ldc, dcomplex *d__, int *ldd, dcomplex *e, int *lde, 
                       dcomplex *f, int *ldf, double *scale, double *dif, dcomplex *work, 
                       int *lwork, int *iwork, int *info);

extern int C2F(ctpcon)(char *norm, char *uplo, char *diag, int *n, 
                       dcomplex *ap, double *rcond, dcomplex *work, double *rwork, int *info);

extern int C2F(ctprfs)(char *uplo, char *trans, char *diag, int *n, 
                       int *nrhs, dcomplex *ap, dcomplex *b, int *ldb, dcomplex *x, 
                       int *ldx, double *ferr, double *berr, dcomplex *work, double *rwork, 
                       int *info);

extern int C2F(ctptri)(char *uplo, char *diag, int *n, dcomplex *ap, 
                       int *info);

extern int C2F(ctptrs)(char *uplo, char *trans, char *diag, int *n, 
                       int *nrhs, dcomplex *ap, dcomplex *b, int *ldb, int *info);

extern int C2F(ctpttf)(char *transr, char *uplo, int *n, dcomplex *
                       ap, dcomplex *arf, int *info);

extern int C2F(ctpttr)(char *uplo, int *n, dcomplex *ap, dcomplex *a, 
                       int *lda, int *info);

extern int C2F(ctrcon)(char *norm, char *uplo, char *diag, int *n, 
                       dcomplex *a, int *lda, double *rcond, dcomplex *work, double *rwork, 
                       int *info);

extern int C2F(ctrevc)(char *side, char *howmny, int *select, 
                       int *n, dcomplex *t, int *ldt, dcomplex *vl, int *ldvl, 
                       dcomplex *vr, int *ldvr, int *mm, int *m, dcomplex *work, 
                       double *rwork, int *info);

extern int C2F(ctrexc)(char *compq, int *n, dcomplex *t, int *
                       ldt, dcomplex *q, int *ldq, int *ifst, int *ilst, int *
                       info);

extern int C2F(ctrrfs)(char *uplo, char *trans, char *diag, int *n, 
                       int *nrhs, dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       dcomplex *x, int *ldx, double *ferr, double *berr, dcomplex *work, double 
                       *rwork, int *info);

extern int C2F(ctrsen)(char *job, char *compq, int *select, int 
                       *n, dcomplex *t, int *ldt, dcomplex *q, int *ldq, dcomplex *w, 
                       int *m, double *s, double *sep, dcomplex *work, int *lwork, 
                       int *info);

extern int C2F(ctrsna)(char *job, char *howmny, int *select, 
                       int *n, dcomplex *t, int *ldt, dcomplex *vl, int *ldvl, 
                       dcomplex *vr, int *ldvr, double *s, double *sep, int *mm, int *
                       m, dcomplex *work, int *ldwork, double *rwork, int *info);

extern int C2F(ctrsyl)(char *trana, char *tranb, int *isgn, int 
                       *m, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       dcomplex *c__, int *ldc, double *scale, int *info);

extern int C2F(ctrti2)(char *uplo, char *diag, int *n, dcomplex *a, 
                       int *lda, int *info);

extern int C2F(ctrtri)(char *uplo, char *diag, int *n, dcomplex *a, 
                       int *lda, int *info);

extern int C2F(ctrtrs)(char *uplo, char *trans, char *diag, int *n, 
                       int *nrhs, dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       int *info);

extern int C2F(ctrttf)(char *transr, char *uplo, int *n, dcomplex *a, 
                       int *lda, dcomplex *arf, int *info);

extern int C2F(ctrttp)(char *uplo, int *n, dcomplex *a, int *lda, 
                       dcomplex *ap, int *info);

extern int C2F(ctzrqf)(int *m, int *n, dcomplex *a, int *lda, 
                       dcomplex *tau, int *info);

extern int C2F(ctzrzf)(int *m, int *n, dcomplex *a, int *lda, 
                       dcomplex *tau, dcomplex *work, int *lwork, int *info);

extern int C2F(cung2l)(int *m, int *n, int *k, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work, int *info);

extern int C2F(cung2r)(int *m, int *n, int *k, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work, int *info);

extern int C2F(cungbr)(char *vect, int *m, int *n, int *k, 
                       dcomplex *a, int *lda, dcomplex *tau, dcomplex *work, int *lwork, 
                       int *info);

extern int C2F(cunghr)(int *n, int *ilo, int *ihi, dcomplex *
                       a, int *lda, dcomplex *tau, dcomplex *work, int *lwork, int 
                       *info);

extern int C2F(cungl2)(int *m, int *n, int *k, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work, int *info);

extern int C2F(cunglq)(int *m, int *n, int *k, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work, int *lwork, int *
                       info);

extern int C2F(cungql)(int *m, int *n, int *k, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work, int *lwork, int *
                       info);

extern int C2F(cungqr)(int *m, int *n, int *k, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work, int *lwork, int *
                       info);

extern int C2F(cungr2)(int *m, int *n, int *k, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work, int *info);

extern int C2F(cungrq)(int *m, int *n, int *k, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work, int *lwork, int *
                       info);

extern int C2F(cungtr)(char *uplo, int *n, dcomplex *a, int *lda, 
                       dcomplex *tau, dcomplex *work, int *lwork, int *info);

extern int C2F(cunm2l)(char *side, char *trans, int *m, int *n, 
                       int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, 
                       int *ldc, dcomplex *work, int *info);

extern int C2F(cunm2r)(char *side, char *trans, int *m, int *n, 
                       int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, 
                       int *ldc, dcomplex *work, int *info);

extern int C2F(cunmbr)(char *vect, char *side, char *trans, int *m, 
                       int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, 
                       dcomplex *c__, int *ldc, dcomplex *work, int *lwork, int *
                       info);

extern int C2F(cunmhr)(char *side, char *trans, int *m, int *n, 
                       int *ilo, int *ihi, dcomplex *a, int *lda, dcomplex *tau, 
                       dcomplex *c__, int *ldc, dcomplex *work, int *lwork, int *
                       info);

extern int C2F(cunml2)(char *side, char *trans, int *m, int *n, 
                       int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, 
                       int *ldc, dcomplex *work, int *info);

extern int C2F(cunmlq)(char *side, char *trans, int *m, int *n, 
                       int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, 
                       int *ldc, dcomplex *work, int *lwork, int *info);

extern int C2F(cunmql)(char *side, char *trans, int *m, int *n, 
                       int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, 
                       int *ldc, dcomplex *work, int *lwork, int *info);

extern int C2F(cunmqr)(char *side, char *trans, int *m, int *n, 
                       int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, 
                       int *ldc, dcomplex *work, int *lwork, int *info);

extern int C2F(cunmr2)(char *side, char *trans, int *m, int *n, 
                       int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, 
                       int *ldc, dcomplex *work, int *info);

extern int C2F(cunmr3)(char *side, char *trans, int *m, int *n, 
                       int *k, int *l, dcomplex *a, int *lda, dcomplex *tau, 
                       dcomplex *c__, int *ldc, dcomplex *work, int *info);

extern int C2F(cunmrq)(char *side, char *trans, int *m, int *n, 
                       int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, 
                       int *ldc, dcomplex *work, int *lwork, int *info);

extern int C2F(cunmrz)(char *side, char *trans, int *m, int *n, 
                       int *k, int *l, dcomplex *a, int *lda, dcomplex *tau, 
                       dcomplex *c__, int *ldc, dcomplex *work, int *lwork, int *
                       info);

extern int C2F(cunmtr)(char *side, char *uplo, char *trans, int *m, 
                       int *n, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, 
                       int *ldc, dcomplex *work, int *lwork, int *info);

extern int C2F(cupgtr)(char *uplo, int *n, dcomplex *ap, dcomplex *
                       tau, dcomplex *q, int *ldq, dcomplex *work, int *info);

extern int C2F(cupmtr)(char *side, char *uplo, char *trans, int *m, 
                       int *n, dcomplex *ap, dcomplex *tau, dcomplex *c__, int *ldc, 
                       dcomplex *work, int *info);

extern int C2F(dbdsdc)(char *uplo, char *compq, int *n, double *
                       d__, double *e, double *u, int *ldu, double *vt, 
                       int *ldvt, double *q, int *iq, double *work, int *
                       iwork, int *info);

extern int C2F(dbdsqr)(char *uplo, int *n, int *ncvt, int *
                       nru, int *ncc, double *d__, double *e, double *vt, 
                       int *ldvt, double *u, int *ldu, double *c__, int *
                       ldc, double *work, int *info);

extern int C2F(ddisna)(char *job, int *m, int *n, double *
                       d__, double *sep, int *info);

extern int C2F(dgbbrd)(char *vect, int *m, int *n, int *ncc, 
                       int *kl, int *ku, double *ab, int *ldab, double *
                       d__, double *e, double *q, int *ldq, double *pt, 
                       int *ldpt, double *c__, int *ldc, double *work, 
                       int *info);

extern int C2F(dgbcon)(char *norm, int *n, int *kl, int *ku, 
                       double *ab, int *ldab, int *ipiv, double *anorm, 
                       double *rcond, double *work, int *iwork, int *info);

extern int C2F(dgbequ)(int *m, int *n, int *kl, int *ku, 
                       double *ab, int *ldab, double *r__, double *c__, 
                       double *rowcnd, double *colcnd, double *amax, int *
                       info);

extern int C2F(dgbequb)(int *m, int *n, int *kl, int *
                        ku, double *ab, int *ldab, double *r__, double *c__, 
                        double *rowcnd, double *colcnd, double *amax, int *
                        info);

extern int C2F(dgbrfs)(char *trans, int *n, int *kl, int *
                       ku, int *nrhs, double *ab, int *ldab, double *afb, 
                       int *ldafb, int *ipiv, double *b, int *ldb, 
                       double *x, int *ldx, double *ferr, double *berr, 
                       double *work, int *iwork, int *info);

extern int C2F(dgbrfsx)(char *trans, char *equed, int *n, int *
                        kl, int *ku, int *nrhs, double *ab, int *ldab, 
                        double *afb, int *ldafb, int *ipiv, double *r__, 
                        double *c__, double *b, int *ldb, double *x, int *
                        ldx, double *rcond, double *berr, int *n_err_bnds__, 
                        double *err_bnds_norm__, double *err_bnds_comp__, int *
                        nparams, double *params, double *work, int *iwork, 
                        int *info);

extern int C2F(dgbsv)(int *n, int *kl, int *ku, int *
                      nrhs, double *ab, int *ldab, int *ipiv, double *b, 
                      int *ldb, int *info);

extern int C2F(dgbsvx)(char *fact, char *trans, int *n, int *kl, 
                       int *ku, int *nrhs, double *ab, int *ldab, 
                       double *afb, int *ldafb, int *ipiv, char *equed, 
                       double *r__, double *c__, double *b, int *ldb, 
                       double *x, int *ldx, double *rcond, double *ferr, 
                       double *berr, double *work, int *iwork, int *info);

extern int C2F(dgbsvxx)(char *fact, char *trans, int *n, int *
                        kl, int *ku, int *nrhs, double *ab, int *ldab, 
                        double *afb, int *ldafb, int *ipiv, char *equed, 
                        double *r__, double *c__, double *b, int *ldb, 
                        double *x, int *ldx, double *rcond, double *rpvgrw, 
                        double *berr, int *n_err_bnds__, double *err_bnds_norm__, 
                        double *err_bnds_comp__, int *nparams, double *params, 
                        double *work, int *iwork, int *info);

extern int C2F(dgbtf2)(int *m, int *n, int *kl, int *ku, 
                       double *ab, int *ldab, int *ipiv, int *info);

extern int C2F(dgbtrf)(int *m, int *n, int *kl, int *ku, 
                       double *ab, int *ldab, int *ipiv, int *info);

extern int C2F(dgbtrs)(char *trans, int *n, int *kl, int *
                       ku, int *nrhs, double *ab, int *ldab, int *ipiv, 
                       double *b, int *ldb, int *info);

extern int C2F(dgebak)(char *job, char *side, int *n, int *ilo, 
                       int *ihi, double *scale, int *m, double *v, int *
                       ldv, int *info);

extern int C2F(dgebal)(char *job, int *n, double *a, int *
                       lda, int *ilo, int *ihi, double *scale, int *info);

extern int C2F(dgebd2)(int *m, int *n, double *a, int *
                       lda, double *d__, double *e, double *tauq, double *
                       taup, double *work, int *info);

extern int C2F(dgebrd)(int *m, int *n, double *a, int *
                       lda, double *d__, double *e, double *tauq, double *
                       taup, double *work, int *lwork, int *info);

extern int C2F(dgecon)(char *norm, int *n, double *a, int *
                       lda, double *anorm, double *rcond, double *work, int *
                       iwork, int *info);

extern int C2F(dgeequ)(int *m, int *n, double *a, int *
                       lda, double *r__, double *c__, double *rowcnd, double 
                       *colcnd, double *amax, int *info);

extern int C2F(dgeequb)(int *m, int *n, double *a, int *
                        lda, double *r__, double *c__, double *rowcnd, double 
                        *colcnd, double *amax, int *info);

extern int C2F(dgees)(char *jobvs, char *sort, L_fp select, int *n, 
                      double *a, int *lda, int *sdim, double *wr, 
                      double *wi, double *vs, int *ldvs, double *work, 
                      int *lwork, int *bwork, int *info);

extern int C2F(dgeesx)(char *jobvs, char *sort, L_fp select, char *
                       sense, int *n, double *a, int *lda, int *sdim, 
                       double *wr, double *wi, double *vs, int *ldvs, 
                       double *rconde, double *rcondv, double *work, int *
                       lwork, int *iwork, int *liwork, int *bwork, int *info);

extern int C2F(dgeev)(char *jobvl, char *jobvr, int *n, double *
                      a, int *lda, double *wr, double *wi, double *vl, 
                      int *ldvl, double *vr, int *ldvr, double *work, 
                      int *lwork, int *info);

extern int C2F(dgeevx)(char *balanc, char *jobvl, char *jobvr, char *
                       sense, int *n, double *a, int *lda, double *wr, 
                       double *wi, double *vl, int *ldvl, double *vr, 
                       int *ldvr, int *ilo, int *ihi, double *scale, 
                       double *abnrm, double *rconde, double *rcondv, double 
                       *work, int *lwork, int *iwork, int *info);

extern int C2F(dgegs)(char *jobvsl, char *jobvsr, int *n, 
                      double *a, int *lda, double *b, int *ldb, double *
                      alphar, double *alphai, double *beta, double *vsl, 
                      int *ldvsl, double *vsr, int *ldvsr, double *work, 
                      int *lwork, int *info);

extern int C2F(dgegv)(char *jobvl, char *jobvr, int *n, double *
                      a, int *lda, double *b, int *ldb, double *alphar, 
                      double *alphai, double *beta, double *vl, int *ldvl, 
                      double *vr, int *ldvr, double *work, int *lwork, 
                      int *info);

extern int C2F(dgehd2)(int *n, int *ilo, int *ihi, 
                       double *a, int *lda, double *tau, double *work, 
                       int *info);

extern int C2F(dgehrd)(int *n, int *ilo, int *ihi, 
                       double *a, int *lda, double *tau, double *work, 
                       int *lwork, int *info);

extern int C2F(dgejsv)(char *joba, char *jobu, char *jobv, char *jobr, 
                       char *jobt, char *jobp, int *m, int *n, double *a, 
                       int *lda, double *sva, double *u, int *ldu, 
                       double *v, int *ldv, double *work, int *lwork, 
                       int *iwork, int *info);

extern int C2F(dgelq2)(int *m, int *n, double *a, int *
                       lda, double *tau, double *work, int *info);

extern int C2F(dgelqf)(int *m, int *n, double *a, int *
                       lda, double *tau, double *work, int *lwork, int *info);

extern int C2F(dgels)(char *trans, int *m, int *n, int *
                      nrhs, double *a, int *lda, double *b, int *ldb, 
                      double *work, int *lwork, int *info);

extern int C2F(dgelsd)(int *m, int *n, int *nrhs, 
                       double *a, int *lda, double *b, int *ldb, double *
                       s, double *rcond, int *rank, double *work, int *lwork, 
                       int *iwork, int *info);

extern int C2F(dgelss)(int *m, int *n, int *nrhs, 
                       double *a, int *lda, double *b, int *ldb, double *
                       s, double *rcond, int *rank, double *work, int *lwork, 
                       int *info);

extern int C2F(dgelsx)(int *m, int *n, int *nrhs, 
                       double *a, int *lda, double *b, int *ldb, int *
                       jpvt, double *rcond, int *rank, double *work, int *
                       info);

extern int C2F(dgelsy)(int *m, int *n, int *nrhs, 
                       double *a, int *lda, double *b, int *ldb, int *
                       jpvt, double *rcond, int *rank, double *work, int *
                       lwork, int *info);

extern int C2F(dgeql2)(int *m, int *n, double *a, int *
                       lda, double *tau, double *work, int *info);

extern int C2F(dgeqlf)(int *m, int *n, double *a, int *
                       lda, double *tau, double *work, int *lwork, int *info);

extern int C2F(dgeqp3)(int *m, int *n, double *a, int *
                       lda, int *jpvt, double *tau, double *work, int *lwork, 
                       int *info);

extern int C2F(dgeqpf)(int *m, int *n, double *a, int *
                       lda, int *jpvt, double *tau, double *work, int *info);

extern int C2F(dgeqr2)(int *m, int *n, double *a, int *
                       lda, double *tau, double *work, int *info);

extern int C2F(dgeqrf)(int *m, int *n, double *a, int *
                       lda, double *tau, double *work, int *lwork, int *info);

extern int C2F(dgerfs)(char *trans, int *n, int *nrhs, 
                       double *a, int *lda, double *af, int *ldaf, int *
                       ipiv, double *b, int *ldb, double *x, int *ldx, 
                       double *ferr, double *berr, double *work, int *iwork, 
                       int *info);

extern int C2F(dgerfsx)(char *trans, char *equed, int *n, int *
                        nrhs, double *a, int *lda, double *af, int *ldaf, 
                        int *ipiv, double *r__, double *c__, double *b, 
                        int *ldb, double *x, int *ldx, double *rcond, 
                        double *berr, int *n_err_bnds__, double *err_bnds_norm__, 
                        double *err_bnds_comp__, int *nparams, double *params, 
                        double *work, int *iwork, int *info);

extern int C2F(dgerq2)(int *m, int *n, double *a, int *
                       lda, double *tau, double *work, int *info);

extern int C2F(dgerqf)(int *m, int *n, double *a, int *
                       lda, double *tau, double *work, int *lwork, int *info);

extern int C2F(dgesc2)(int *n, double *a, int *lda, 
                       double *rhs, int *ipiv, int *jpiv, double *scale);

extern int C2F(dgesdd)(char *jobz, int *m, int *n, double *
                       a, int *lda, double *s, double *u, int *ldu, 
                       double *vt, int *ldvt, double *work, int *lwork, 
                       int *iwork, int *info);

extern int C2F(dgesv)(int *n, int *nrhs, double *a, int 
                      *lda, int *ipiv, double *b, int *ldb, int *info);

extern int C2F(dgesvd)(char *jobu, char *jobvt, int *m, int *n, 
                       double *a, int *lda, double *s, double *u, int *
                       ldu, double *vt, int *ldvt, double *work, int *lwork, 
                       int *info);

extern int C2F(dgesvj)(char *joba, char *jobu, char *jobv, int *m, 
                       int *n, double *a, int *lda, double *sva, int *mv, 
                       double *v, int *ldv, double *work, int *lwork, 
                       int *info);

extern int C2F(dgesvx)(char *fact, char *trans, int *n, int *
                       nrhs, double *a, int *lda, double *af, int *ldaf, 
                       int *ipiv, char *equed, double *r__, double *c__, 
                       double *b, int *ldb, double *x, int *ldx, double *
                       rcond, double *ferr, double *berr, double *work, int *
                       iwork, int *info);

extern int C2F(dgesvxx)(char *fact, char *trans, int *n, int *
                        nrhs, double *a, int *lda, double *af, int *ldaf, 
                        int *ipiv, char *equed, double *r__, double *c__, 
                        double *b, int *ldb, double *x, int *ldx, double *
                        rcond, double *rpvgrw, double *berr, int *n_err_bnds__, 
                        double *err_bnds_norm__, double *err_bnds_comp__, int *
                        nparams, double *params, double *work, int *iwork, 
                        int *info);

extern int C2F(dgetc2)(int *n, double *a, int *lda, int 
                       *ipiv, int *jpiv, int *info);

extern int C2F(dgetf2)(int *m, int *n, double *a, int *
                       lda, int *ipiv, int *info);

extern int C2F(dgetrf)(int *m, int *n, double *a, int *
                       lda, int *ipiv, int *info);

extern int C2F(dgetri)(int *n, double *a, int *lda, int 
                       *ipiv, double *work, int *lwork, int *info);

extern int C2F(dgetrs)(char *trans, int *n, int *nrhs, 
                       double *a, int *lda, int *ipiv, double *b, int *
                       ldb, int *info);

extern int C2F(dggbak)(char *job, char *side, int *n, int *ilo, 
                       int *ihi, double *lscale, double *rscale, int *m, 
                       double *v, int *ldv, int *info);

extern int C2F(dggbal)(char *job, int *n, double *a, int *
                       lda, double *b, int *ldb, int *ilo, int *ihi, 
                       double *lscale, double *rscale, double *work, int *
                       info);

extern int C2F(dgges)(char *jobvsl, char *jobvsr, char *sort, L_fp 
                      selctg, int *n, double *a, int *lda, double *b, 
                      int *ldb, int *sdim, double *alphar, double *alphai, 
                      double *beta, double *vsl, int *ldvsl, double *vsr, 
                      int *ldvsr, double *work, int *lwork, int *bwork, 
                      int *info);

extern int C2F(dggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp 
                       selctg, char *sense, int *n, double *a, int *lda, 
                       double *b, int *ldb, int *sdim, double *alphar, 
                       double *alphai, double *beta, double *vsl, int *ldvsl, 
                       double *vsr, int *ldvsr, double *rconde, double *
                       rcondv, double *work, int *lwork, int *iwork, int *
                       liwork, int *bwork, int *info);

extern int C2F(dggev)(char *jobvl, char *jobvr, int *n, double *
                      a, int *lda, double *b, int *ldb, double *alphar, 
                      double *alphai, double *beta, double *vl, int *ldvl, 
                      double *vr, int *ldvr, double *work, int *lwork, 
                      int *info);

extern int C2F(dggevx)(char *balanc, char *jobvl, char *jobvr, char *
                       sense, int *n, double *a, int *lda, double *b, 
                       int *ldb, double *alphar, double *alphai, double *
                       beta, double *vl, int *ldvl, double *vr, int *ldvr, 
                       int *ilo, int *ihi, double *lscale, double *rscale, 
                       double *abnrm, double *bbnrm, double *rconde, double *
                       rcondv, double *work, int *lwork, int *iwork, int *
                       bwork, int *info);

extern int C2F(dggglm)(int *n, int *m, int *p, double *
                       a, int *lda, double *b, int *ldb, double *d__, 
                       double *x, double *y, double *work, int *lwork, 
                       int *info);

extern int C2F(dgghrd)(char *compq, char *compz, int *n, int *
                       ilo, int *ihi, double *a, int *lda, double *b, 
                       int *ldb, double *q, int *ldq, double *z__, int *
                       ldz, int *info);

extern int C2F(dgglse)(int *m, int *n, int *p, double *
                       a, int *lda, double *b, int *ldb, double *c__, 
                       double *d__, double *x, double *work, int *lwork, 
                       int *info);

extern int C2F(dggqrf)(int *n, int *m, int *p, double *
                       a, int *lda, double *taua, double *b, int *ldb, 
                       double *taub, double *work, int *lwork, int *info);

extern int C2F(dggrqf)(int *m, int *p, int *n, double *
                       a, int *lda, double *taua, double *b, int *ldb, 
                       double *taub, double *work, int *lwork, int *info);

extern int C2F(dggsvd)(char *jobu, char *jobv, char *jobq, int *m, 
                       int *n, int *p, int *k, int *l, double *a, 
                       int *lda, double *b, int *ldb, double *alpha, 
                       double *beta, double *u, int *ldu, double *v, int 
                       *ldv, double *q, int *ldq, double *work, int *iwork, 
                       int *info);

extern int C2F(dggsvp)(char *jobu, char *jobv, char *jobq, int *m, 
                       int *p, int *n, double *a, int *lda, double *b, 
                       int *ldb, double *tola, double *tolb, int *k, int 
                       *l, double *u, int *ldu, double *v, int *ldv, 
                       double *q, int *ldq, int *iwork, double *tau, 
                       double *work, int *info);

extern int C2F(dgsvj0)(char *jobv, int *m, int *n, double *
                       a, int *lda, double *d__, double *sva, int *mv, 
                       double *v, int *ldv, double *eps, double *sfmin, 
                       double *tol, int *nsweep, double *work, int *lwork, 
                       int *info);

extern int C2F(dgsvj1)(char *jobv, int *m, int *n, int *n1, 
                       double *a, int *lda, double *d__, double *sva, 
                       int *mv, double *v, int *ldv, double *eps, double 
                       *sfmin, double *tol, int *nsweep, double *work, int *
                       lwork, int *info);

extern int C2F(dgtcon)(char *norm, int *n, double *dl, 
                       double *d__, double *du, double *du2, int *ipiv, 
                       double *anorm, double *rcond, double *work, int *
                       iwork, int *info);

extern int C2F(dgtrfs)(char *trans, int *n, int *nrhs, 
                       double *dl, double *d__, double *du, double *dlf, 
                       double *df, double *duf, double *du2, int *ipiv, 
                       double *b, int *ldb, double *x, int *ldx, double *
                       ferr, double *berr, double *work, int *iwork, int *
                       info);

extern int C2F(dgtsv)(int *n, int *nrhs, double *dl, 
                      double *d__, double *du, double *b, int *ldb, int 
                      *info);

extern int C2F(dgtsvx)(char *fact, char *trans, int *n, int *
                       nrhs, double *dl, double *d__, double *du, double *
                       dlf, double *df, double *duf, double *du2, int *ipiv, 
                       double *b, int *ldb, double *x, int *ldx, double *
                       rcond, double *ferr, double *berr, double *work, int *
                       iwork, int *info);

extern int C2F(dgttrf)(int *n, double *dl, double *d__, 
                       double *du, double *du2, int *ipiv, int *info);

extern int C2F(dgttrs)(char *trans, int *n, int *nrhs, 
                       double *dl, double *d__, double *du, double *du2, 
                       int *ipiv, double *b, int *ldb, int *info);

extern int C2F(dgtts2)(int *itrans, int *n, int *nrhs, 
                       double *dl, double *d__, double *du, double *du2, 
                       int *ipiv, double *b, int *ldb);

extern int C2F(dhgeqz)(char *job, char *compq, char *compz, int *n, 
                       int *ilo, int *ihi, double *h__, int *ldh, double 
                       *t, int *ldt, double *alphar, double *alphai, double *
                       beta, double *q, int *ldq, double *z__, int *ldz, 
                       double *work, int *lwork, int *info);

extern int C2F(dhsein)(char *side, char *eigsrc, char *initv, int *
                       select, int *n, double *h__, int *ldh, double *wr, 
                       double *wi, double *vl, int *ldvl, double *vr, 
                       int *ldvr, int *mm, int *m, double *work, int *
                       ifaill, int *ifailr, int *info);

extern int C2F(dhseqr)(char *job, char *compz, int *n, int *ilo, 
                       int *ihi, double *h__, int *ldh, double *wr, 
                       double *wi, double *z__, int *ldz, double *work, 
                       int *lwork, int *info);

int C2F(disnan)(double *din);

extern int C2F(dla_gbamv_)(int *trans, int *m, int *n, 
                           int *kl, int *ku, double *alpha, double *ab, int *
                           ldab, double *x, int *incx, double *beta, double *y, 
                           int *incy);

double C2F(dla_gbrcond_)(char *trans, int *n, int *kl, int *ku, 
                         double *ab, int *ldab, double *afb, int *ldafb, 
                         int *ipiv, int *cmode, double *c__, int *info, 
                         double *work, int *iwork, int trans_len);

extern int C2F(dla_gbrfsx_extended_)(int *prec_type__, int *
                                     trans_type__, int *n, int *kl, int *ku, int *nrhs, 
                                     double *ab, int *ldab, double *afb, int *ldafb, 
                                     int *ipiv, int *colequ, double *c__, double *b, 
                                     int *ldb, double *y, int *ldy, double *berr_out__, 
                                     int *n_norms__, double *errs_n__, double *errs_c__, 
                                     double *res, double *ayb, double *dy, double *
                                     y_tail__, double *rcond, int *ithresh, double *rthresh, 
                                     double *dz_ub__, int *ignore_cwise__, int *info);

double C2F(dla_gbrpvgrw_)(int *n, int *kl, int *ku, int *
                          ncols, double *ab, int *ldab, double *afb, int *ldafb);

extern int C2F(dla_geamv_)(int *trans, int *m, int *n, 
                           double *alpha, double *a, int *lda, double *x, 
                           int *incx, double *beta, double *y, int *incy);

double C2F(dla_gercond_)(char *trans, int *n, double *a, int *lda,
                         double *af, int *ldaf, int *ipiv, int *cmode, 
                         double *c__, int *info, double *work, int *iwork, 
                         int trans_len);

extern int C2F(dla_gerfsx_extended_)(int *prec_type__, int *
                                     trans_type__, int *n, int *nrhs, double *a, int *lda, 
                                     double *af, int *ldaf, int *ipiv, int *colequ, 
                                     double *c__, double *b, int *ldb, double *y, int *
                                     ldy, double *berr_out__, int *n_norms__, double *errs_n__,
                                     double *errs_c__, double *res, double *ayb, double *
                                     dy, double *y_tail__, double *rcond, int *ithresh, 
                                     double *rthresh, double *dz_ub__, int *ignore_cwise__, 
                                     int *info);

extern int C2F(dla_lin_berr_)(int *n, int *nz, int *nrhs, 
                              double *res, double *ayb, double *berr);

double C2F(dla_porcond_)(char *uplo, int *n, double *a, int *lda, 
                         double *af, int *ldaf, int *cmode, double *c__, 
                         int *info, double *work, int *iwork, int uplo_len);

extern int C2F(dla_porfsx_extended_)(int *prec_type__, char *uplo, 
                                     int *n, int *nrhs, double *a, int *lda, double *
                                     af, int *ldaf, int *colequ, double *c__, double *b, 
                                     int *ldb, double *y, int *ldy, double *berr_out__, 
                                     int *n_norms__, double *errs_n__, double *errs_c__, 
                                     double *res, double *ayb, double *dy, double *
                                     y_tail__, double *rcond, int *ithresh, double *rthresh, 
                                     double *dz_ub__, int *ignore_cwise__, int *info, int 
                                     uplo_len);

double C2F(dla_porpvgrw_)(char *uplo, int *ncols, double *a, int *
                          lda, double *af, int *ldaf, double *work, int uplo_len);

double C2F(dla_rpvgrw_)(int *n, int *ncols, double *a, int *
                        lda, double *af, int *ldaf);

extern int C2F(dla_syamv_)(int *uplo, int *n, double *alpha,
                           double *a, int *lda, double *x, int *incx, 
                           double *beta, double *y, int *incy);

double C2F(dla_syrcond_)(char *uplo, int *n, double *a, int *lda, 
                         double *af, int *ldaf, int *ipiv, int *cmode, 
                         double *c__, int *info, double *work, int *iwork, 
                         int uplo_len);

extern int C2F(dla_syrfsx_extended_)(int *prec_type__, char *uplo, 
                                     int *n, int *nrhs, double *a, int *lda, double *
                                     af, int *ldaf, int *ipiv, int *colequ, double *c__, 
                                     double *b, int *ldb, double *y, int *ldy, double *
                                     berr_out__, int *n_norms__, double *errs_n__, double *
                                     errs_c__, double *res, double *ayb, double *dy, 
                                     double *y_tail__, double *rcond, int *ithresh, double 
                                     *rthresh, double *dz_ub__, int *ignore_cwise__, int *info,
                                     int uplo_len);

double C2F(dla_syrpvgrw_)(char *uplo, int *n, int *info, double *
                          a, int *lda, double *af, int *ldaf, int *ipiv, 
                          double *work, int uplo_len);

extern int C2F(dla_wwaddw_)(int *n, double *x, double *y, 
                            double *w);

extern int C2F(dlabad)(double *small, double *large);

extern int C2F(dlabrd)(int *m, int *n, int *nb, double *
                       a, int *lda, double *d__, double *e, double *tauq, 
                       double *taup, double *x, int *ldx, double *y, int 
                       *ldy);

extern int C2F(dlacn2)(int *n, double *v, double *x, 
                       int *isgn, double *est, int *kase, int *isave);

extern int C2F(dlacon)(int *n, double *v, double *x, 
                       int *isgn, double *est, int *kase);

extern int C2F(dlacpy)(char *uplo, int *m, int *n, double *
                       a, int *lda, double *b, int *ldb);

extern int C2F(dladiv)(double *a, double *b, double *c__, 
                       double *d__, double *p, double *q);

extern int C2F(dlae2)(double *a, double *b, double *c__, 
                      double *rt1, double *rt2);

extern int C2F(dlaebz)(int *ijob, int *nitmax, int *n, 
                       int *mmax, int *minp, int *nbmin, double *abstol, 
                       double *reltol, double *pivmin, double *d__, double *
                       e, double *e2, int *nval, double *ab, double *c__, 
                       int *mout, int *nab, double *work, int *iwork, 
                       int *info);

extern int C2F(dlaed0)(int *icompq, int *qsiz, int *n, 
                       double *d__, double *e, double *q, int *ldq, 
                       double *qstore, int *ldqs, double *work, int *iwork, 
                       int *info);

extern int C2F(dlaed1)(int *n, double *d__, double *q, 
                       int *ldq, int *indxq, double *rho, int *cutpnt, 
                       double *work, int *iwork, int *info);

extern int C2F(dlaed2)(int *k, int *n, int *n1, double *
                       d__, double *q, int *ldq, int *indxq, double *rho, 
                       double *z__, double *dlamda, double *w, double *q2, 
                       int *indx, int *indxc, int *indxp, int *coltyp, 
                       int *info);

extern int C2F(dlaed3)(int *k, int *n, int *n1, double *
                       d__, double *q, int *ldq, double *rho, double *dlamda, 
                       double *q2, int *indx, int *ctot, double *w, 
                       double *s, int *info);

extern int C2F(dlaed4)(int *n, int *i__, double *d__, 
                       double *z__, double *delta, double *rho, double *dlam, 
                       int *info);

extern int C2F(dlaed5)(int *i__, double *d__, double *z__, 
                       double *delta, double *rho, double *dlam);

extern int C2F(dlaed6)(int *kniter, int *orgati, double *
                       rho, double *d__, double *z__, double *finit, double *
                       tau, int *info);

extern int C2F(dlaed7)(int *icompq, int *n, int *qsiz, 
                       int *tlvls, int *curlvl, int *curpbm, double *d__, 
                       double *q, int *ldq, int *indxq, double *rho, int 
                       *cutpnt, double *qstore, int *qptr, int *prmptr, int *
                       perm, int *givptr, int *givcol, double *givnum, 
                       double *work, int *iwork, int *info);

extern int C2F(dlaed8)(int *icompq, int *k, int *n, int 
                       *qsiz, double *d__, double *q, int *ldq, int *indxq, 
                       double *rho, int *cutpnt, double *z__, double *dlamda, 
                       double *q2, int *ldq2, double *w, int *perm, int 
                       *givptr, int *givcol, double *givnum, int *indxp, int 
                       *indx, int *info);

extern int C2F(dlaed9)(int *k, int *kstart, int *kstop, 
                       int *n, double *d__, double *q, int *ldq, double *
                       rho, double *dlamda, double *w, double *s, int *lds, 
                       int *info);

extern int C2F(dlaeda)(int *n, int *tlvls, int *curlvl, 
                       int *curpbm, int *prmptr, int *perm, int *givptr, 
                       int *givcol, double *givnum, double *q, int *qptr, 
                       double *z__, double *ztemp, int *info);

extern int C2F(dlaein)(int *rightv, int *noinit, int *n, 
                       double *h__, int *ldh, double *wr, double *wi, 
                       double *vr, double *vi, double *b, int *ldb, 
                       double *work, double *eps3, double *smlnum, double *
                       bignum, int *info);

extern int C2F(dlaev2)(double *a, double *b, double *c__, 
                       double *rt1, double *rt2, double *cs1, double *sn1);

extern int C2F(dlaexc)(int *wantq, int *n, double *t, 
                       int *ldt, double *q, int *ldq, int *j1, int *n1, 
                       int *n2, double *work, int *info);

extern int C2F(dlag2)(double *a, int *lda, double *b, 
                      int *ldb, double *safmin, double *scale1, double *
                      scale2, double *wr1, double *wr2, double *wi);

extern int C2F(dlag2s)(int *m, int *n, double *a, int *
                       lda, double *sa, int *ldsa, int *info);

extern int C2F(dlags2)(int *upper, double *a1, double *a2, 
                       double *a3, double *b1, double *b2, double *b3, 
                       double *csu, double *snu, double *csv, double *snv, 
                       double *csq, double *snq);

extern int C2F(dlagtf)(int *n, double *a, double *lambda, 
                       double *b, double *c__, double *tol, double *d__, 
                       int *in, int *info);

extern int C2F(dlagtm)(char *trans, int *n, int *nrhs, 
                       double *alpha, double *dl, double *d__, double *du, 
                       double *x, int *ldx, double *beta, double *b, int 
                       *ldb);

extern int C2F(dlagts)(int *job, int *n, double *a, 
                       double *b, double *c__, double *d__, int *in, 
                       double *y, double *tol, int *info);

extern int C2F(dlagv2)(double *a, int *lda, double *b, 
                       int *ldb, double *alphar, double *alphai, double *
                       beta, double *csl, double *snl, double *csr, double *
                       snr);

extern int C2F(dlahqr)(int *wantt, int *wantz, int *n, 
                       int *ilo, int *ihi, double *h__, int *ldh, double 
                       *wr, double *wi, int *iloz, int *ihiz, double *z__, 
                       int *ldz, int *info);

extern int C2F(dlahr2)(int *n, int *k, int *nb, double *
                       a, int *lda, double *tau, double *t, int *ldt, 
                       double *y, int *ldy);

extern int C2F(dlahrd)(int *n, int *k, int *nb, double *
                       a, int *lda, double *tau, double *t, int *ldt, 
                       double *y, int *ldy);

extern int C2F(dlaic1)(int *job, int *j, double *x, 
                       double *sest, double *w, double *gamma, double *
                       sestpr, double *s, double *c__);

int C2F(dlaisnan)(double *din1, double *din2);

extern int C2F(dlaln2)(int *ltrans, int *na, int *nw, 
                       double *smin, double *ca, double *a, int *lda, 
                       double *d1, double *d2, double *b, int *ldb, 
                       double *wr, double *wi, double *x, int *ldx, 
                       double *scale, double *xnorm, int *info);

extern int C2F(dlals0)(int *icompq, int *nl, int *nr, 
                       int *sqre, int *nrhs, double *b, int *ldb, double 
                       *bx, int *ldbx, int *perm, int *givptr, int *givcol, 
                       int *ldgcol, double *givnum, int *ldgnum, double *
                       poles, double *difl, double *difr, double *z__, int *
                       k, double *c__, double *s, double *work, int *info);

extern int C2F(dlalsa)(int *icompq, int *smlsiz, int *n, 
                       int *nrhs, double *b, int *ldb, double *bx, int *
                       ldbx, double *u, int *ldu, double *vt, int *k, 
                       double *difl, double *difr, double *z__, double *
                       poles, int *givptr, int *givcol, int *ldgcol, int *
                       perm, double *givnum, double *c__, double *s, double *
                       work, int *iwork, int *info);

extern int C2F(dlalsd)(char *uplo, int *smlsiz, int *n, int 
                       *nrhs, double *d__, double *e, double *b, int *ldb, 
                       double *rcond, int *rank, double *work, int *iwork, 
                       int *info);

extern int C2F(dlamrg)(int *n1, int *n2, double *a, int 
                       *dtrd1, int *dtrd2, int *index);

int C2F(dlaneg)(int *n, double *d__, double *lld, double *
                sigma, double *pivmin, int *r__);

double C2F(dlangb)(char *norm, int *n, int *kl, int *ku, 
                   double *ab, int *ldab, double *work);

double C2F(dlange)(char *norm, int *m, int *n, double *a, int 
                   *lda, double *work);

double C2F(dlangt)(char *norm, int *n, double *dl, double *d__, 
                   double *du);

double C2F(dlanhs)(char *norm, int *n, double *a, int *lda, 
                   double *work);

double C2F(dlansb)(char *norm, char *uplo, int *n, int *k, double 
                   *ab, int *ldab, double *work);

double C2F(dlansf)(char *norm, char *transr, char *uplo, int *n, 
                   double *a, double *work);

double C2F(dlansp)(char *norm, char *uplo, int *n, double *ap, 
                   double *work);

double C2F(dlanst)(char *norm, int *n, double *d__, double *e);

double C2F(dlansy)(char *norm, char *uplo, int *n, double *a, int 
                   *lda, double *work);

double C2F(dlantb)(char *norm, char *uplo, char *diag, int *n, int *k, 
                   double *ab, int *ldab, double *work);

double C2F(dlantp)(char *norm, char *uplo, char *diag, int *n, double 
                   *ap, double *work);

double C2F(dlantr)(char *norm, char *uplo, char *diag, int *m, int *n, 
                   double *a, int *lda, double *work);

extern int C2F(dlanv2)(double *a, double *b, double *c__, 
                       double *d__, double *rt1r, double *rt1i, double *rt2r, 
                       double *rt2i, double *cs, double *sn);

extern int C2F(dlapll)(int *n, double *x, int *incx, 
                       double *y, int *incy, double *ssmin);

extern int C2F(dlapmt)(int *forwrd, int *m, int *n, 
                       double *x, int *ldx, int *k);

double C2F(dlapy2)(double *x, double *y);

double C2F(dlapy3)(double *x, double *y, double *z__);

extern int C2F(dlaqgb)(int *m, int *n, int *kl, int *ku, 
                       double *ab, int *ldab, double *r__, double *c__, 
                       double *rowcnd, double *colcnd, double *amax, char *equed);

extern int C2F(dlaqge)(int *m, int *n, double *a, int *
                       lda, double *r__, double *c__, double *rowcnd, double 
                       *colcnd, double *amax, char *equed);

extern int C2F(dlaqp2)(int *m, int *n, int *offset, 
                       double *a, int *lda, int *jpvt, double *tau, 
                       double *vn1, double *vn2, double *work);

extern int C2F(dlaqps)(int *m, int *n, int *offset, int 
                       *nb, int *kb, double *a, int *lda, int *jpvt, 
                       double *tau, double *vn1, double *vn2, double *auxv, 
                       double *f, int *ldf);

extern int C2F(dlaqr0)(int *wantt, int *wantz, int *n, 
                       int *ilo, int *ihi, double *h__, int *ldh, double 
                       *wr, double *wi, int *iloz, int *ihiz, double *z__, 
                       int *ldz, double *work, int *lwork, int *info);

extern int C2F(dlaqr1)(int *n, double *h__, int *ldh, 
                       double *sr1, double *si1, double *sr2, double *si2, 
                       double *v);

extern int C2F(dlaqr2)(int *wantt, int *wantz, int *n, 
                       int *ktop, int *kbot, int *nw, double *h__, int *
                       ldh, int *iloz, int *ihiz, double *z__, int *ldz, 
                       int *ns, int *nd, double *sr, double *si, double *
                       v, int *ldv, int *nh, double *t, int *ldt, int *
                       nv, double *wv, int *ldwv, double *work, int *lwork);

extern int C2F(dlaqr3)(int *wantt, int *wantz, int *n, 
                       int *ktop, int *kbot, int *nw, double *h__, int *
                       ldh, int *iloz, int *ihiz, double *z__, int *ldz, 
                       int *ns, int *nd, double *sr, double *si, double *
                       v, int *ldv, int *nh, double *t, int *ldt, int *
                       nv, double *wv, int *ldwv, double *work, int *lwork);

extern int C2F(dlaqr4)(int *wantt, int *wantz, int *n, 
                       int *ilo, int *ihi, double *h__, int *ldh, double 
                       *wr, double *wi, int *iloz, int *ihiz, double *z__, 
                       int *ldz, double *work, int *lwork, int *info);

extern int C2F(dlaqr5)(int *wantt, int *wantz, int *kacc22, 
                       int *n, int *ktop, int *kbot, int *nshfts, double 
                       *sr, double *si, double *h__, int *ldh, int *iloz, 
                       int *ihiz, double *z__, int *ldz, double *v, int *
                       ldv, double *u, int *ldu, int *nv, double *wv, 
                       int *ldwv, int *nh, double *wh, int *ldwh);

extern int C2F(dlaqsb)(char *uplo, int *n, int *kd, double *
                       ab, int *ldab, double *s, double *scond, double *amax, 
                       char *equed);

extern int C2F(dlaqsp)(char *uplo, int *n, double *ap, 
                       double *s, double *scond, double *amax, char *equed);

extern int C2F(dlaqsy)(char *uplo, int *n, double *a, int *
                       lda, double *s, double *scond, double *amax, char *equed);

extern int C2F(dlaqtr)(int *ltran, int *ldouble, int *n, 
                       double *t, int *ldt, double *b, double *w, double 
                       *scale, double *x, double *work, int *info);

extern int C2F(dlar1v)(int *n, int *b1, int *bn, double 
                       *lambda, double *d__, double *l, double *ld, double *
                       lld, double *pivmin, double *gaptol, double *z__, int 
                       *wantnc, int *negcnt, double *ztz, double *mingma, 
                       int *r__, int *isuppz, double *nrminv, double *resid, 
                       double *rqcorr, double *work);

extern int C2F(dlar2v)(int *n, double *x, double *y, 
                       double *z__, int *incx, double *c__, double *s, 
                       int *incc);

extern int C2F(dlarf)(char *side, int *m, int *n, double *v, 
                      int *incv, double *tau, double *c__, int *ldc, 
                      double *work);

extern int C2F(dlarfb)(char *side, char *trans, char *direct, char *
                       storev, int *m, int *n, int *k, double *v, int *
                       ldv, double *t, int *ldt, double *c__, int *ldc, 
                       double *work, int *ldwork);

extern int C2F(dlarfg)(int *n, double *alpha, double *x, 
                       int *incx, double *tau);

extern int C2F(dlarfp)(int *n, double *alpha, double *x, 
                       int *incx, double *tau);

extern int C2F(dlarft)(char *direct, char *storev, int *n, int *
                       k, double *v, int *ldv, double *tau, double *t, 
                       int *ldt);

extern int C2F(dlarfx)(char *side, int *m, int *n, double *
                       v, double *tau, double *c__, int *ldc, double *work);

extern int C2F(dlargv)(int *n, double *x, int *incx, 
                       double *y, int *incy, double *c__, int *incc);

extern int C2F(dlarnv)(int *idist, int *iseed, int *n, 
                       double *x);

extern int C2F(dlarra)(int *n, double *d__, double *e, 
                       double *e2, double *spltol, double *tnrm, int *nsplit, 
                       int *isplit, int *info);

extern int C2F(dlarrb)(int *n, double *d__, double *lld, 
                       int *ifirst, int *ilast, double *rtol1, double *rtol2, 
                       int *offset, double *w, double *wgap, double *werr, 
                       double *work, int *iwork, double *pivmin, double *
                       spdiam, int *twist, int *info);

extern int C2F(dlarrc)(char *jobt, int *n, double *vl, 
                       double *vu, double *d__, double *e, double *pivmin, 
                       int *eigcnt, int *lcnt, int *rcnt, int *info);

extern int C2F(dlarrd)(char *range, char *order, int *n, double 
                       *vl, double *vu, int *il, int *iu, double *gers, 
                       double *reltol, double *d__, double *e, double *e2, 
                       double *pivmin, int *nsplit, int *isplit, int *m, 
                       double *w, double *werr, double *wl, double *wu, 
                       int *iblock, int *indexw, double *work, int *iwork, 
                       int *info);

extern int C2F(dlarre)(char *range, int *n, double *vl, 
                       double *vu, int *il, int *iu, double *d__, double 
                       *e, double *e2, double *rtol1, double *rtol2, double *
                       spltol, int *nsplit, int *isplit, int *m, double *w, 
                       double *werr, double *wgap, int *iblock, int *indexw, 
                       double *gers, double *pivmin, double *work, int *
                       iwork, int *info);

extern int C2F(dlarrf)(int *n, double *d__, double *l, 
                       double *ld, int *clstrt, int *clend, double *w, 
                       double *wgap, double *werr, double *spdiam, double *
                       clgapl, double *clgapr, double *pivmin, double *sigma, 
                       double *dplus, double *lplus, double *work, int *info);

extern int C2F(dlarrj)(int *n, double *d__, double *e2, 
                       int *ifirst, int *ilast, double *rtol, int *offset, 
                       double *w, double *werr, double *work, int *iwork, 
                       double *pivmin, double *spdiam, int *info);

extern int C2F(dlarrk)(int *n, int *iw, double *gl, 
                       double *gu, double *d__, double *e2, double *pivmin, 
                       double *reltol, double *w, double *werr, int *info);

extern int C2F(dlarrr)(int *n, double *d__, double *e, 
                       int *info);

extern int C2F(dlarrv)(int *n, double *vl, double *vu, 
                       double *d__, double *l, double *pivmin, int *isplit, 
                       int *m, int *dol, int *dou, double *minrgp, 
                       double *rtol1, double *rtol2, double *w, double *werr, 
                       double *wgap, int *iblock, int *indexw, double *gers, 
                       double *z__, int *ldz, int *isuppz, double *work, 
                       int *iwork, int *info);

extern int C2F(dlarscl2)(int *m, int *n, double *d__, 
                         double *x, int *ldx);

extern int C2F(dlartg)(double *f, double *g, double *cs, 
                       double *sn, double *r__);

extern int C2F(dlartv)(int *n, double *x, int *incx, 
                       double *y, int *incy, double *c__, double *s, int 
                       *incc);

extern int C2F(dlaruv)(int *iseed, int *n, double *x);

extern int C2F(dlarz)(char *side, int *m, int *n, int *l, 
                      double *v, int *incv, double *tau, double *c__, 
                      int *ldc, double *work);

extern int C2F(dlarzb)(char *side, char *trans, char *direct, char *
                       storev, int *m, int *n, int *k, int *l, double *v, 
                       int *ldv, double *t, int *ldt, double *c__, int *
                       ldc, double *work, int *ldwork);

extern int C2F(dlarzt)(char *direct, char *storev, int *n, int *
                       k, double *v, int *ldv, double *tau, double *t, 
                       int *ldt);

extern int C2F(dlas2)(double *f, double *g, double *h__, 
                      double *ssmin, double *ssmax);

extern int C2F(dlascl)(char *type__, int *kl, int *ku, 
                       double *cfrom, double *cto, int *m, int *n, 
                       double *a, int *lda, int *info);

extern int C2F(dlascl2)(int *m, int *n, double *d__, 
                        double *x, int *ldx);

extern int C2F(dlasd0)(int *n, int *sqre, double *d__, 
                       double *e, double *u, int *ldu, double *vt, int *
                       ldvt, int *smlsiz, int *iwork, double *work, int *
                       info);

extern int C2F(dlasd1)(int *nl, int *nr, int *sqre, 
                       double *d__, double *alpha, double *beta, double *u, 
                       int *ldu, double *vt, int *ldvt, int *idxq, int *
                       iwork, double *work, int *info);

extern int C2F(dlasd2)(int *nl, int *nr, int *sqre, int 
                       *k, double *d__, double *z__, double *alpha, double *
                       beta, double *u, int *ldu, double *vt, int *ldvt, 
                       double *dsigma, double *u2, int *ldu2, double *vt2, 
                       int *ldvt2, int *idxp, int *idx, int *idxc, int *
                       idxq, int *coltyp, int *info);

extern int C2F(dlasd3)(int *nl, int *nr, int *sqre, int 
                       *k, double *d__, double *q, int *ldq, double *dsigma, 
                       double *u, int *ldu, double *u2, int *ldu2, 
                       double *vt, int *ldvt, double *vt2, int *ldvt2, 
                       int *idxc, int *ctot, double *z__, int *info);

extern int C2F(dlasd4)(int *n, int *i__, double *d__, 
                       double *z__, double *delta, double *rho, double *
                       sigma, double *work, int *info);

extern int C2F(dlasd5)(int *i__, double *d__, double *z__, 
                       double *delta, double *rho, double *dsigma, double *
                       work);

extern int C2F(dlasd6)(int *icompq, int *nl, int *nr, 
                       int *sqre, double *d__, double *vf, double *vl, 
                       double *alpha, double *beta, int *idxq, int *perm, 
                       int *givptr, int *givcol, int *ldgcol, double *givnum, 
                       int *ldgnum, double *poles, double *difl, double *
                       difr, double *z__, int *k, double *c__, double *s, 
                       double *work, int *iwork, int *info);

extern int C2F(dlasd7)(int *icompq, int *nl, int *nr, 
                       int *sqre, int *k, double *d__, double *z__, 
                       double *zw, double *vf, double *vfw, double *vl, 
                       double *vlw, double *alpha, double *beta, double *
                       dsigma, int *idx, int *idxp, int *idxq, int *perm, 
                       int *givptr, int *givcol, int *ldgcol, double *givnum, 
                       int *ldgnum, double *c__, double *s, int *info);

extern int C2F(dlasd8)(int *icompq, int *k, double *d__, 
                       double *z__, double *vf, double *vl, double *difl, 
                       double *difr, int *lddifr, double *dsigma, double *
                       work, int *info);

extern int C2F(dlasda)(int *icompq, int *smlsiz, int *n, 
                       int *sqre, double *d__, double *e, double *u, int 
                       *ldu, double *vt, int *k, double *difl, double *difr, 
                       double *z__, double *poles, int *givptr, int *givcol, 
                       int *ldgcol, int *perm, double *givnum, double *c__, 
                       double *s, double *work, int *iwork, int *info);

extern int C2F(dlasdq)(char *uplo, int *sqre, int *n, int *
                       ncvt, int *nru, int *ncc, double *d__, double *e, 
                       double *vt, int *ldvt, double *u, int *ldu, 
                       double *c__, int *ldc, double *work, int *info);

extern int C2F(dlasdt)(int *n, int *lvl, int *nd, int *
                       inode, int *ndiml, int *ndimr, int *msub);

extern int C2F(dlaset)(char *uplo, int *m, int *n, double *
                       alpha, double *beta, double *a, int *lda);

extern int C2F(dlasq1)(int *n, double *d__, double *e, 
                       double *work, int *info);

extern int C2F(dlasq2)(int *n, double *z__, int *info);

extern int C2F(dlasq3)(int *i0, int *n0, double *z__, 
                       int *pp, double *dmin__, double *sigma, double *desig, 
                       double *qmax, int *nfail, int *iter, int *ndiv, 
                       int *ieee, int *ttype, double *dmin1, double *dmin2, 
                       double *dn, double *dn1, double *dn2, double *g, 
                       double *tau);

extern int C2F(dlasq4)(int *i0, int *n0, double *z__, 
                       int *pp, int *n0in, double *dmin__, double *dmin1, 
                       double *dmin2, double *dn, double *dn1, double *dn2, 
                       double *tau, int *ttype, double *g);

extern int C2F(dlasq5)(int *i0, int *n0, double *z__, 
                       int *pp, double *tau, double *dmin__, double *dmin1, 
                       double *dmin2, double *dn, double *dnm1, double *dnm2, 
                       int *ieee);

extern int C2F(dlasq6)(int *i0, int *n0, double *z__, 
                       int *pp, double *dmin__, double *dmin1, double *dmin2, 
                       double *dn, double *dnm1, double *dnm2);

extern int C2F(dlasr)(char *side, char *pivot, char *direct, int *m, 
                      int *n, double *c__, double *s, double *a, int *
                      lda);

extern int C2F(dlasrt)(char *id, int *n, double *d__, int *
                       info);

extern int C2F(dlassq)(int *n, double *x, int *incx, 
                       double *scale, double *sumsq);

extern int C2F(dlasv2)(double *f, double *g, double *h__, 
                       double *ssmin, double *ssmax, double *snr, double *
                       csr, double *snl, double *csl);

extern int C2F(dlaswp)(int *n, double *a, int *lda, int 
                       *k1, int *k2, int *ipiv, int *incx);

extern int C2F(dlasy2)(int *ltranl, int *ltranr, int *isgn, 
                       int *n1, int *n2, double *tl, int *ldtl, double *
                       tr, int *ldtr, double *b, int *ldb, double *scale, 
                       double *x, int *ldx, double *xnorm, int *info);

extern int C2F(dlasyf)(char *uplo, int *n, int *nb, int *kb, 
                       double *a, int *lda, int *ipiv, double *w, int *
                       ldw, int *info);

extern int C2F(dlat2s)(char *uplo, int *n, double *a, int *
                       lda, double *sa, int *ldsa, int *info);

extern int C2F(dlatbs)(char *uplo, char *trans, char *diag, char *
                       normin, int *n, int *kd, double *ab, int *ldab, 
                       double *x, double *scale, double *cnorm, int *info);

extern int C2F(dlatdf)(int *ijob, int *n, double *z__, 
                       int *ldz, double *rhs, double *rdsum, double *rdscal, 
                       int *ipiv, int *jpiv);

extern int C2F(dlatps)(char *uplo, char *trans, char *diag, char *
                       normin, int *n, double *ap, double *x, double *scale, 
                       double *cnorm, int *info);

extern int C2F(dlatrd)(char *uplo, int *n, int *nb, double *
                       a, int *lda, double *e, double *tau, double *w, 
                       int *ldw);

extern int C2F(dlatrs)(char *uplo, char *trans, char *diag, char *
                       normin, int *n, double *a, int *lda, double *x, 
                       double *scale, double *cnorm, int *info);

extern int C2F(dlatrz)(int *m, int *n, int *l, double *
                       a, int *lda, double *tau, double *work);

extern int C2F(dlatzm)(char *side, int *m, int *n, double *
                       v, int *incv, double *tau, double *c1, double *c2, 
                       int *ldc, double *work);

extern int C2F(dlauu2)(char *uplo, int *n, double *a, int *
                       lda, int *info);

extern int C2F(dlauum)(char *uplo, int *n, double *a, int *
                       lda, int *info);

extern int C2F(dopgtr)(char *uplo, int *n, double *ap, 
                       double *tau, double *q, int *ldq, double *work, 
                       int *info);

extern int C2F(dopmtr)(char *side, char *uplo, char *trans, int *m, 
                       int *n, double *ap, double *tau, double *c__, int 
                       *ldc, double *work, int *info);

extern int C2F(dorg2l)(int *m, int *n, int *k, double *
                       a, int *lda, double *tau, double *work, int *info);

extern int C2F(dorg2r)(int *m, int *n, int *k, double *
                       a, int *lda, double *tau, double *work, int *info);

extern int C2F(dorgbr)(char *vect, int *m, int *n, int *k, 
                       double *a, int *lda, double *tau, double *work, 
                       int *lwork, int *info);

extern int C2F(dorghr)(int *n, int *ilo, int *ihi, 
                       double *a, int *lda, double *tau, double *work, 
                       int *lwork, int *info);

extern int C2F(dorgl2)(int *m, int *n, int *k, double *
                       a, int *lda, double *tau, double *work, int *info);

extern int C2F(dorglq)(int *m, int *n, int *k, double *
                       a, int *lda, double *tau, double *work, int *lwork, 
                       int *info);

extern int C2F(dorgql)(int *m, int *n, int *k, double *
                       a, int *lda, double *tau, double *work, int *lwork, 
                       int *info);

extern int C2F(dorgqr)(int *m, int *n, int *k, double *
                       a, int *lda, double *tau, double *work, int *lwork, 
                       int *info);

extern int C2F(dorgr2)(int *m, int *n, int *k, double *
                       a, int *lda, double *tau, double *work, int *info);

extern int C2F(dorgrq)(int *m, int *n, int *k, double *
                       a, int *lda, double *tau, double *work, int *lwork, 
                       int *info);

extern int C2F(dorgtr)(char *uplo, int *n, double *a, int *
                       lda, double *tau, double *work, int *lwork, int *info);

extern int C2F(dorm2l)(char *side, char *trans, int *m, int *n, 
                       int *k, double *a, int *lda, double *tau, double *
                       c__, int *ldc, double *work, int *info);

extern int C2F(dorm2r)(char *side, char *trans, int *m, int *n, 
                       int *k, double *a, int *lda, double *tau, double *
                       c__, int *ldc, double *work, int *info);

extern int C2F(dormbr)(char *vect, char *side, char *trans, int *m, 
                       int *n, int *k, double *a, int *lda, double *tau, 
                       double *c__, int *ldc, double *work, int *lwork, 
                       int *info);

extern int C2F(dormhr)(char *side, char *trans, int *m, int *n, 
                       int *ilo, int *ihi, double *a, int *lda, double *
                       tau, double *c__, int *ldc, double *work, int *lwork, 
                       int *info);

extern int C2F(dorml2)(char *side, char *trans, int *m, int *n, 
                       int *k, double *a, int *lda, double *tau, double *
                       c__, int *ldc, double *work, int *info);

extern int C2F(dormlq)(char *side, char *trans, int *m, int *n, 
                       int *k, double *a, int *lda, double *tau, double *
                       c__, int *ldc, double *work, int *lwork, int *info);

extern int C2F(dormql)(char *side, char *trans, int *m, int *n, 
                       int *k, double *a, int *lda, double *tau, double *
                       c__, int *ldc, double *work, int *lwork, int *info);

extern int C2F(dormqr)(char *side, char *trans, int *m, int *n, 
                       int *k, double *a, int *lda, double *tau, double *
                       c__, int *ldc, double *work, int *lwork, int *info);

extern int C2F(dormr2)(char *side, char *trans, int *m, int *n, 
                       int *k, double *a, int *lda, double *tau, double *
                       c__, int *ldc, double *work, int *info);

extern int C2F(dormr3)(char *side, char *trans, int *m, int *n, 
                       int *k, int *l, double *a, int *lda, double *tau, 
                       double *c__, int *ldc, double *work, int *info);

extern int C2F(dormrq)(char *side, char *trans, int *m, int *n, 
                       int *k, double *a, int *lda, double *tau, double *
                       c__, int *ldc, double *work, int *lwork, int *info);

extern int C2F(dormrz)(char *side, char *trans, int *m, int *n, 
                       int *k, int *l, double *a, int *lda, double *tau, 
                       double *c__, int *ldc, double *work, int *lwork, 
                       int *info);

extern int C2F(dormtr)(char *side, char *uplo, char *trans, int *m, 
                       int *n, double *a, int *lda, double *tau, double *
                       c__, int *ldc, double *work, int *lwork, int *info);

extern int C2F(dpbcon)(char *uplo, int *n, int *kd, double *
                       ab, int *ldab, double *anorm, double *rcond, double *
                       work, int *iwork, int *info);

extern int C2F(dpbequ)(char *uplo, int *n, int *kd, double *
                       ab, int *ldab, double *s, double *scond, double *amax, 
                       int *info);

extern int C2F(dpbrfs)(char *uplo, int *n, int *kd, int *
                       nrhs, double *ab, int *ldab, double *afb, int *ldafb, 
                       double *b, int *ldb, double *x, int *ldx, double *
                       ferr, double *berr, double *work, int *iwork, int *
                       info);

extern int C2F(dpbstf)(char *uplo, int *n, int *kd, double *
                       ab, int *ldab, int *info);

extern int C2F(dpbsv)(char *uplo, int *n, int *kd, int *
                      nrhs, double *ab, int *ldab, double *b, int *ldb, 
                      int *info);

extern int C2F(dpbsvx)(char *fact, char *uplo, int *n, int *kd, 
                       int *nrhs, double *ab, int *ldab, double *afb, 
                       int *ldafb, char *equed, double *s, double *b, int *
                       ldb, double *x, int *ldx, double *rcond, double *ferr, 
                       double *berr, double *work, int *iwork, int *info);

extern int C2F(dpbtf2)(char *uplo, int *n, int *kd, double *
                       ab, int *ldab, int *info);

extern int C2F(dpbtrf)(char *uplo, int *n, int *kd, double *
                       ab, int *ldab, int *info);

extern int C2F(dpbtrs)(char *uplo, int *n, int *kd, int *
                       nrhs, double *ab, int *ldab, double *b, int *ldb, 
                       int *info);

extern int C2F(dpftrf)(char *transr, char *uplo, int *n, double 
                       *a, int *info);

extern int C2F(dpftri)(char *transr, char *uplo, int *n, double 
                       *a, int *info);

extern int C2F(dpftrs)(char *transr, char *uplo, int *n, int *
                       nrhs, double *a, double *b, int *ldb, int *info);

extern int C2F(dpocon)(char *uplo, int *n, double *a, int *
                       lda, double *anorm, double *rcond, double *work, int *
                       iwork, int *info);

extern int C2F(dpoequ)(int *n, double *a, int *lda, 
                       double *s, double *scond, double *amax, int *info);

extern int C2F(dpoequb)(int *n, double *a, int *lda, 
                        double *s, double *scond, double *amax, int *info);

extern int C2F(dporfs)(char *uplo, int *n, int *nrhs, 
                       double *a, int *lda, double *af, int *ldaf, 
                       double *b, int *ldb, double *x, int *ldx, double *
                       ferr, double *berr, double *work, int *iwork, int *
                       info);

extern int C2F(dporfsx)(char *uplo, char *equed, int *n, int *
                        nrhs, double *a, int *lda, double *af, int *ldaf, 
                        double *s, double *b, int *ldb, double *x, int *
                        ldx, double *rcond, double *berr, int *n_err_bnds__, 
                        double *err_bnds_norm__, double *err_bnds_comp__, int *
                        nparams, double *params, double *work, int *iwork, 
                        int *info);

extern int C2F(dposv)(char *uplo, int *n, int *nrhs, double 
                      *a, int *lda, double *b, int *ldb, int *info);

extern int C2F(dposvx)(char *fact, char *uplo, int *n, int *
                       nrhs, double *a, int *lda, double *af, int *ldaf, 
                       char *equed, double *s, double *b, int *ldb, double *
                       x, int *ldx, double *rcond, double *ferr, double *
                       berr, double *work, int *iwork, int *info);

extern int C2F(dposvxx)(char *fact, char *uplo, int *n, int *
                        nrhs, double *a, int *lda, double *af, int *ldaf, 
                        char *equed, double *s, double *b, int *ldb, double *
                        x, int *ldx, double *rcond, double *rpvgrw, double *
                        berr, int *n_err_bnds__, double *err_bnds_norm__, double *
                        err_bnds_comp__, int *nparams, double *params, double *
                        work, int *iwork, int *info);

extern int C2F(dpotf2)(char *uplo, int *n, double *a, int *
                       lda, int *info);

extern int C2F(dpotrf)(char *uplo, int *n, double *a, int *
                       lda, int *info);

extern int C2F(dpotri)(char *uplo, int *n, double *a, int *
                       lda, int *info);

extern int C2F(dpotrs)(char *uplo, int *n, int *nrhs, 
                       double *a, int *lda, double *b, int *ldb, int *
                       info);

extern int C2F(dppcon)(char *uplo, int *n, double *ap, 
                       double *anorm, double *rcond, double *work, int *
                       iwork, int *info);

extern int C2F(dppequ)(char *uplo, int *n, double *ap, 
                       double *s, double *scond, double *amax, int *info);

extern int C2F(dpprfs)(char *uplo, int *n, int *nrhs, 
                       double *ap, double *afp, double *b, int *ldb, 
                       double *x, int *ldx, double *ferr, double *berr, 
                       double *work, int *iwork, int *info);

extern int C2F(dppsv)(char *uplo, int *n, int *nrhs, double 
                      *ap, double *b, int *ldb, int *info);

extern int C2F(dppsvx)(char *fact, char *uplo, int *n, int *
                       nrhs, double *ap, double *afp, char *equed, double *s, 
                       double *b, int *ldb, double *x, int *ldx, double *
                       rcond, double *ferr, double *berr, double *work, int *
                       iwork, int *info);

extern int C2F(dpptrf)(char *uplo, int *n, double *ap, int *
                       info);

extern int C2F(dpptri)(char *uplo, int *n, double *ap, int *
                       info);

extern int C2F(dpptrs)(char *uplo, int *n, int *nrhs, 
                       double *ap, double *b, int *ldb, int *info);

extern int C2F(dpstf2)(char *uplo, int *n, double *a, int *
                       lda, int *piv, int *rank, double *tol, double *work, 
                       int *info);

extern int C2F(dpstrf)(char *uplo, int *n, double *a, int *
                       lda, int *piv, int *rank, double *tol, double *work, 
                       int *info);

extern int C2F(dptcon)(int *n, double *d__, double *e, 
                       double *anorm, double *rcond, double *work, int *info);

extern int C2F(dpteqr)(char *compz, int *n, double *d__, 
                       double *e, double *z__, int *ldz, double *work, 
                       int *info);

extern int C2F(dptrfs)(int *n, int *nrhs, double *d__, 
                       double *e, double *df, double *ef, double *b, int 
                       *ldb, double *x, int *ldx, double *ferr, double *berr, 
                       double *work, int *info);

extern int C2F(dptsv)(int *n, int *nrhs, double *d__, 
                      double *e, double *b, int *ldb, int *info);

extern int C2F(dptsvx)(char *fact, int *n, int *nrhs, 
                       double *d__, double *e, double *df, double *ef, 
                       double *b, int *ldb, double *x, int *ldx, double *
                       rcond, double *ferr, double *berr, double *work, int *
                       info);

extern int C2F(dpttrf)(int *n, double *d__, double *e, 
                       int *info);

extern int C2F(dpttrs)(int *n, int *nrhs, double *d__, 
                       double *e, double *b, int *ldb, int *info);

extern int C2F(dptts2)(int *n, int *nrhs, double *d__, 
                       double *e, double *b, int *ldb);

extern int C2F(drscl)(int *n, double *sa, double *sx, 
                      int *incx);

extern int C2F(dsbev)(char *jobz, char *uplo, int *n, int *kd, 
                      double *ab, int *ldab, double *w, double *z__, 
                      int *ldz, double *work, int *info);

extern int C2F(dsbevd)(char *jobz, char *uplo, int *n, int *kd, 
                       double *ab, int *ldab, double *w, double *z__, 
                       int *ldz, double *work, int *lwork, int *iwork, 
                       int *liwork, int *info);

extern int C2F(dsbevx)(char *jobz, char *range, char *uplo, int *n, 
                       int *kd, double *ab, int *ldab, double *q, int *
                       ldq, double *vl, double *vu, int *il, int *iu, 
                       double *abstol, int *m, double *w, double *z__, 
                       int *ldz, double *work, int *iwork, int *ifail, 
                       int *info);

extern int C2F(dsbgst)(char *vect, char *uplo, int *n, int *ka, 
                       int *kb, double *ab, int *ldab, double *bb, int *
                       ldbb, double *x, int *ldx, double *work, int *info);

extern int C2F(dsbgv)(char *jobz, char *uplo, int *n, int *ka, 
                      int *kb, double *ab, int *ldab, double *bb, int *
                      ldbb, double *w, double *z__, int *ldz, double *work, 
                      int *info);

extern int C2F(dsbgvd)(char *jobz, char *uplo, int *n, int *ka, 
                       int *kb, double *ab, int *ldab, double *bb, int *
                       ldbb, double *w, double *z__, int *ldz, double *work, 
                       int *lwork, int *iwork, int *liwork, int *info);

extern int C2F(dsbgvx)(char *jobz, char *range, char *uplo, int *n, 
                       int *ka, int *kb, double *ab, int *ldab, double *
                       bb, int *ldbb, double *q, int *ldq, double *vl, 
                       double *vu, int *il, int *iu, double *abstol, int 
                       *m, double *w, double *z__, int *ldz, double *work, 
                       int *iwork, int *ifail, int *info);

extern int C2F(dsbtrd)(char *vect, char *uplo, int *n, int *kd, 
                       double *ab, int *ldab, double *d__, double *e, 
                       double *q, int *ldq, double *work, int *info);

extern int C2F(dsfrk)(char *transr, char *uplo, char *trans, int *n, 
                      int *k, double *alpha, double *a, int *lda, 
                      double *beta, double *c__);

extern int C2F(dsgesv)(int *n, int *nrhs, double *a, 
                       int *lda, int *ipiv, double *b, int *ldb, double *
                       x, int *ldx, double *work, double *swork, int *iter, 
                       int *info);

extern int C2F(dspcon)(char *uplo, int *n, double *ap, int *
                       ipiv, double *anorm, double *rcond, double *work, int 
                       *iwork, int *info);

extern int C2F(dspev)(char *jobz, char *uplo, int *n, double *
                      ap, double *w, double *z__, int *ldz, double *work, 
                      int *info);

extern int C2F(dspevd)(char *jobz, char *uplo, int *n, double *
                       ap, double *w, double *z__, int *ldz, double *work, 
                       int *lwork, int *iwork, int *liwork, int *info);

extern int C2F(dspevx)(char *jobz, char *range, char *uplo, int *n, 
                       double *ap, double *vl, double *vu, int *il, int *
                       iu, double *abstol, int *m, double *w, double *z__, 
                       int *ldz, double *work, int *iwork, int *ifail, 
                       int *info);

extern int C2F(dspgst)(int *itype, char *uplo, int *n, 
                       double *ap, double *bp, int *info);

extern int C2F(dspgv)(int *itype, char *jobz, char *uplo, int *
                      n, double *ap, double *bp, double *w, double *z__, 
                      int *ldz, double *work, int *info);

extern int C2F(dspgvd)(int *itype, char *jobz, char *uplo, int *
                       n, double *ap, double *bp, double *w, double *z__, 
                       int *ldz, double *work, int *lwork, int *iwork, 
                       int *liwork, int *info);

extern int C2F(dspgvx)(int *itype, char *jobz, char *range, char *
                       uplo, int *n, double *ap, double *bp, double *vl, 
                       double *vu, int *il, int *iu, double *abstol, int 
                       *m, double *w, double *z__, int *ldz, double *work, 
                       int *iwork, int *ifail, int *info);

extern int C2F(dsposv)(char *uplo, int *n, int *nrhs, 
                       double *a, int *lda, double *b, int *ldb, double *
                       x, int *ldx, double *work, double *swork, int *iter, 
                       int *info);

extern int C2F(dsprfs)(char *uplo, int *n, int *nrhs, 
                       double *ap, double *afp, int *ipiv, double *b, 
                       int *ldb, double *x, int *ldx, double *ferr, 
                       double *berr, double *work, int *iwork, int *info);

extern int C2F(dspsv)(char *uplo, int *n, int *nrhs, double 
                      *ap, int *ipiv, double *b, int *ldb, int *info);

extern int C2F(dspsvx)(char *fact, char *uplo, int *n, int *
                       nrhs, double *ap, double *afp, int *ipiv, double *b, 
                       int *ldb, double *x, int *ldx, double *rcond, 
                       double *ferr, double *berr, double *work, int *iwork, 
                       int *info);

extern int C2F(dsptrd)(char *uplo, int *n, double *ap, 
                       double *d__, double *e, double *tau, int *info);

extern int C2F(dsptrf)(char *uplo, int *n, double *ap, int *
                       ipiv, int *info);

extern int C2F(dsptri)(char *uplo, int *n, double *ap, int *
                       ipiv, double *work, int *info);

extern int C2F(dsptrs)(char *uplo, int *n, int *nrhs, 
                       double *ap, int *ipiv, double *b, int *ldb, int *
                       info);

extern int C2F(dstebz)(char *range, char *order, int *n, double 
                       *vl, double *vu, int *il, int *iu, double *abstol, 
                       double *d__, double *e, int *m, int *nsplit, 
                       double *w, int *iblock, int *isplit, double *work, 
                       int *iwork, int *info);

extern int C2F(dstedc)(char *compz, int *n, double *d__, 
                       double *e, double *z__, int *ldz, double *work, 
                       int *lwork, int *iwork, int *liwork, int *info);

extern int C2F(dstegr)(char *jobz, char *range, int *n, double *
                       d__, double *e, double *vl, double *vu, int *il, 
                       int *iu, double *abstol, int *m, double *w, 
                       double *z__, int *ldz, int *isuppz, double *work, 
                       int *lwork, int *iwork, int *liwork, int *info);

extern int C2F(dstein)(int *n, double *d__, double *e, 
                       int *m, double *w, int *iblock, int *isplit, 
                       double *z__, int *ldz, double *work, int *iwork, 
                       int *ifail, int *info);

extern int C2F(dstemr)(char *jobz, char *range, int *n, double *
                       d__, double *e, double *vl, double *vu, int *il, 
                       int *iu, int *m, double *w, double *z__, int *ldz, 
                       int *nzc, int *isuppz, int *tryrac, double *work, 
                       int *lwork, int *iwork, int *liwork, int *info);

extern int C2F(dsteqr)(char *compz, int *n, double *d__, 
                       double *e, double *z__, int *ldz, double *work, 
                       int *info);

extern int C2F(dsterf)(int *n, double *d__, double *e, 
                       int *info);

extern int C2F(dstev)(char *jobz, int *n, double *d__, 
                      double *e, double *z__, int *ldz, double *work, 
                      int *info);

extern int C2F(dstevd)(char *jobz, int *n, double *d__, 
                       double *e, double *z__, int *ldz, double *work, 
                       int *lwork, int *iwork, int *liwork, int *info);

extern int C2F(dstevr)(char *jobz, char *range, int *n, double *
                       d__, double *e, double *vl, double *vu, int *il, 
                       int *iu, double *abstol, int *m, double *w, 
                       double *z__, int *ldz, int *isuppz, double *work, 
                       int *lwork, int *iwork, int *liwork, int *info);

extern int C2F(dstevx)(char *jobz, char *range, int *n, double *
                       d__, double *e, double *vl, double *vu, int *il, 
                       int *iu, double *abstol, int *m, double *w, 
                       double *z__, int *ldz, double *work, int *iwork, 
                       int *ifail, int *info);

extern int C2F(dsycon)(char *uplo, int *n, double *a, int *
                       lda, int *ipiv, double *anorm, double *rcond, double *
                       work, int *iwork, int *info);

extern int C2F(dsyequb)(char *uplo, int *n, double *a, int *
                        lda, double *s, double *scond, double *amax, double *
                        work, int *info);

extern int C2F(dsyev)(char *jobz, char *uplo, int *n, double *a, 
                      int *lda, double *w, double *work, int *lwork, 
                      int *info);

extern int C2F(dsyevd)(char *jobz, char *uplo, int *n, double *
                       a, int *lda, double *w, double *work, int *lwork, 
                       int *iwork, int *liwork, int *info);

extern int C2F(dsyevr)(char *jobz, char *range, char *uplo, int *n, 
                       double *a, int *lda, double *vl, double *vu, int *
                       il, int *iu, double *abstol, int *m, double *w, 
                       double *z__, int *ldz, int *isuppz, double *work, 
                       int *lwork, int *iwork, int *liwork, int *info);

extern int C2F(dsyevx)(char *jobz, char *range, char *uplo, int *n, 
                       double *a, int *lda, double *vl, double *vu, int *
                       il, int *iu, double *abstol, int *m, double *w, 
                       double *z__, int *ldz, double *work, int *lwork, 
                       int *iwork, int *ifail, int *info);

extern int C2F(dsygs2)(int *itype, char *uplo, int *n, 
                       double *a, int *lda, double *b, int *ldb, int *
                       info);

extern int C2F(dsygst)(int *itype, char *uplo, int *n, 
                       double *a, int *lda, double *b, int *ldb, int *
                       info);

extern int C2F(dsygv)(int *itype, char *jobz, char *uplo, int *
                      n, double *a, int *lda, double *b, int *ldb, 
                      double *w, double *work, int *lwork, int *info);

extern int C2F(dsygvd)(int *itype, char *jobz, char *uplo, int *
                       n, double *a, int *lda, double *b, int *ldb, 
                       double *w, double *work, int *lwork, int *iwork, 
                       int *liwork, int *info);

extern int C2F(dsygvx)(int *itype, char *jobz, char *range, char *
                       uplo, int *n, double *a, int *lda, double *b, int 
                       *ldb, double *vl, double *vu, int *il, int *iu, 
                       double *abstol, int *m, double *w, double *z__, 
                       int *ldz, double *work, int *lwork, int *iwork, 
                       int *ifail, int *info);

extern int C2F(dsyrfs)(char *uplo, int *n, int *nrhs, 
                       double *a, int *lda, double *af, int *ldaf, int *
                       ipiv, double *b, int *ldb, double *x, int *ldx, 
                       double *ferr, double *berr, double *work, int *iwork, 
                       int *info);

extern int C2F(dsyrfsx)(char *uplo, char *equed, int *n, int *
                        nrhs, double *a, int *lda, double *af, int *ldaf, 
                        int *ipiv, double *s, double *b, int *ldb, double 
                        *x, int *ldx, double *rcond, double *berr, int *
                        n_err_bnds__, double *err_bnds_norm__, double *
                        err_bnds_comp__, int *nparams, double *params, double *
                        work, int *iwork, int *info);

extern int C2F(dsysv)(char *uplo, int *n, int *nrhs, double 
                      *a, int *lda, int *ipiv, double *b, int *ldb, 
                      double *work, int *lwork, int *info);

extern int C2F(dsysvx)(char *fact, char *uplo, int *n, int *
                       nrhs, double *a, int *lda, double *af, int *ldaf, 
                       int *ipiv, double *b, int *ldb, double *x, int *
                       ldx, double *rcond, double *ferr, double *berr, 
                       double *work, int *lwork, int *iwork, int *info);

extern int C2F(dsysvxx)(char *fact, char *uplo, int *n, int *
                        nrhs, double *a, int *lda, double *af, int *ldaf, 
                        int *ipiv, char *equed, double *s, double *b, int *
                        ldb, double *x, int *ldx, double *rcond, double *
                        rpvgrw, double *berr, int *n_err_bnds__, double *
                        err_bnds_norm__, double *err_bnds_comp__, int *nparams, 
                        double *params, double *work, int *iwork, int *info);

extern int C2F(dsytd2)(char *uplo, int *n, double *a, int *
                       lda, double *d__, double *e, double *tau, int *info);

extern int C2F(dsytf2)(char *uplo, int *n, double *a, int *
                       lda, int *ipiv, int *info);

extern int C2F(dsytrd)(char *uplo, int *n, double *a, int *
                       lda, double *d__, double *e, double *tau, double *
                       work, int *lwork, int *info);

extern int C2F(dsytrf)(char *uplo, int *n, double *a, int *
                       lda, int *ipiv, double *work, int *lwork, int *info);

extern int C2F(dsytri)(char *uplo, int *n, double *a, int *
                       lda, int *ipiv, double *work, int *info);

extern int C2F(dsytrs)(char *uplo, int *n, int *nrhs, 
                       double *a, int *lda, int *ipiv, double *b, int *
                       ldb, int *info);

extern int C2F(dtbcon)(char *norm, char *uplo, char *diag, int *n, 
                       int *kd, double *ab, int *ldab, double *rcond, 
                       double *work, int *iwork, int *info);

extern int C2F(dtbrfs)(char *uplo, char *trans, char *diag, int *n, 
                       int *kd, int *nrhs, double *ab, int *ldab, double 
                       *b, int *ldb, double *x, int *ldx, double *ferr, 
                       double *berr, double *work, int *iwork, int *info);

extern int C2F(dtbtrs)(char *uplo, char *trans, char *diag, int *n, 
                       int *kd, int *nrhs, double *ab, int *ldab, double 
                       *b, int *ldb, int *info);

extern int C2F(dtfsm)(char *transr, char *side, char *uplo, char *trans, 
                      char *diag, int *m, int *n, double *alpha, double *a, 
                      double *b, int *ldb);

extern int C2F(dtftri)(char *transr, char *uplo, char *diag, int *n, 
                       double *a, int *info);

extern int C2F(dtfttp)(char *transr, char *uplo, int *n, double 
                       *arf, double *ap, int *info);

extern int C2F(dtfttr)(char *transr, char *uplo, int *n, double 
                       *arf, double *a, int *lda, int *info);

extern int C2F(dtgevc)(char *side, char *howmny, int *select, 
                       int *n, double *s, int *lds, double *p, int *ldp, 
                       double *vl, int *ldvl, double *vr, int *ldvr, int 
                       *mm, int *m, double *work, int *info);

extern int C2F(dtgex2)(int *wantq, int *wantz, int *n, 
                       double *a, int *lda, double *b, int *ldb, double *
                       q, int *ldq, double *z__, int *ldz, int *j1, int *
                       n1, int *n2, double *work, int *lwork, int *info);

extern int C2F(dtgexc)(int *wantq, int *wantz, int *n, 
                       double *a, int *lda, double *b, int *ldb, double *
                       q, int *ldq, double *z__, int *ldz, int *ifst, 
                       int *ilst, double *work, int *lwork, int *info);

extern int C2F(dtgsen)(int *ijob, int *wantq, int *wantz, 
                       int *select, int *n, double *a, int *lda, double *
                       b, int *ldb, double *alphar, double *alphai, double *
                       beta, double *q, int *ldq, double *z__, int *ldz, 
                       int *m, double *pl, double *pr, double *dif, 
                       double *work, int *lwork, int *iwork, int *liwork, 
                       int *info);

extern int C2F(dtgsja)(char *jobu, char *jobv, char *jobq, int *m, 
                       int *p, int *n, int *k, int *l, double *a, 
                       int *lda, double *b, int *ldb, double *tola, 
                       double *tolb, double *alpha, double *beta, double *u, 
                       int *ldu, double *v, int *ldv, double *q, int *
                       ldq, double *work, int *ncycle, int *info);

extern int C2F(dtgsna)(char *job, char *howmny, int *select, 
                       int *n, double *a, int *lda, double *b, int *ldb, 
                       double *vl, int *ldvl, double *vr, int *ldvr, 
                       double *s, double *dif, int *mm, int *m, double *
                       work, int *lwork, int *iwork, int *info);

extern int C2F(dtgsy2)(char *trans, int *ijob, int *m, int *
                       n, double *a, int *lda, double *b, int *ldb, 
                       double *c__, int *ldc, double *d__, int *ldd, 
                       double *e, int *lde, double *f, int *ldf, double *
                       scale, double *rdsum, double *rdscal, int *iwork, int 
                       *pq, int *info);

extern int C2F(dtgsyl)(char *trans, int *ijob, int *m, int *
                       n, double *a, int *lda, double *b, int *ldb, 
                       double *c__, int *ldc, double *d__, int *ldd, 
                       double *e, int *lde, double *f, int *ldf, double *
                       scale, double *dif, double *work, int *lwork, int *
                       iwork, int *info);

extern int C2F(dtpcon)(char *norm, char *uplo, char *diag, int *n, 
                       double *ap, double *rcond, double *work, int *iwork, 
                       int *info);

extern int C2F(dtprfs)(char *uplo, char *trans, char *diag, int *n, 
                       int *nrhs, double *ap, double *b, int *ldb, 
                       double *x, int *ldx, double *ferr, double *berr, 
                       double *work, int *iwork, int *info);

extern int C2F(dtptri)(char *uplo, char *diag, int *n, double *
                       ap, int *info);

extern int C2F(dtptrs)(char *uplo, char *trans, char *diag, int *n, 
                       int *nrhs, double *ap, double *b, int *ldb, int *
                       info);

extern int C2F(dtpttf)(char *transr, char *uplo, int *n, double 
                       *ap, double *arf, int *info);

extern int C2F(dtpttr)(char *uplo, int *n, double *ap, 
                       double *a, int *lda, int *info);

extern int C2F(dtrcon)(char *norm, char *uplo, char *diag, int *n, 
                       double *a, int *lda, double *rcond, double *work, 
                       int *iwork, int *info);

extern int C2F(dtrevc)(char *side, char *howmny, int *select, 
                       int *n, double *t, int *ldt, double *vl, int *
                       ldvl, double *vr, int *ldvr, int *mm, int *m, 
                       double *work, int *info);

extern int C2F(dtrexc)(char *compq, int *n, double *t, int *
                       ldt, double *q, int *ldq, int *ifst, int *ilst, 
                       double *work, int *info);

extern int C2F(dtrrfs)(char *uplo, char *trans, char *diag, int *n, 
                       int *nrhs, double *a, int *lda, double *b, int *
                       ldb, double *x, int *ldx, double *ferr, double *berr, 
                       double *work, int *iwork, int *info);

extern int C2F(dtrsen)(char *job, char *compq, int *select, int 
                       *n, double *t, int *ldt, double *q, int *ldq, 
                       double *wr, double *wi, int *m, double *s, double 
                       *sep, double *work, int *lwork, int *iwork, int *
                       liwork, int *info);

extern int C2F(dtrsna)(char *job, char *howmny, int *select, 
                       int *n, double *t, int *ldt, double *vl, int *
                       ldvl, double *vr, int *ldvr, double *s, double *sep, 
                       int *mm, int *m, double *work, int *ldwork, int *
                       iwork, int *info);

extern int C2F(dtrsyl)(char *trana, char *tranb, int *isgn, int 
                       *m, int *n, double *a, int *lda, double *b, int *
                       ldb, double *c__, int *ldc, double *scale, int *info);

extern int C2F(dtrti2)(char *uplo, char *diag, int *n, double *
                       a, int *lda, int *info);

extern int C2F(dtrtri)(char *uplo, char *diag, int *n, double *
                       a, int *lda, int *info);

extern int C2F(dtrtrs)(char *uplo, char *trans, char *diag, int *n, 
                       int *nrhs, double *a, int *lda, double *b, int *
                       ldb, int *info);

extern int C2F(dtrttf)(char *transr, char *uplo, int *n, double 
                       *a, int *lda, double *arf, int *info);

extern int C2F(dtrttp)(char *uplo, int *n, double *a, int *
                       lda, double *ap, int *info);

extern int C2F(dtzrqf)(int *m, int *n, double *a, int *
                       lda, double *tau, int *info);

extern int C2F(dtzrzf)(int *m, int *n, double *a, int *
                       lda, double *tau, double *work, int *lwork, int *info);

double C2F(dzsum1)(int *n, dcomplex *cx, int *incx);

int C2F(icmax1)(int *n, dcomplex *cx, int *incx);

int C2F(ieeeck)(int *ispec, double *zero, double *one);

int C2F(ilaclc)(int *m, int *n, dcomplex *a, int *lda);

int C2F(ilaclr)(int *m, int *n, dcomplex *a, int *lda);

int C2F(iladiag)(char *diag);

int C2F(iladlc)(int *m, int *n, double *a, int *lda);

int C2F(iladlr)(int *m, int *n, double *a, int *lda);

int C2F(ilaenv)(int *ispec, char *name__, char *opts, int *n1, 
                int *n2, int *n3, int *n4);

int C2F(ilaprec)(char *prec);

int C2F(ilaslc)(int *m, int *n, double *a, int *lda);

int C2F(ilaslr)(int *m, int *n, double *a, int *lda);

int C2F(ilatrans)(char *trans);

int C2F(ilauplo)(char *uplo);

extern int C2F(ilaver)(int *vers_major__, int *vers_minor__, 
                       int *vers_patch__);

int C2F(ilazlc)(int *m, int *n, dcomplex *a, int *lda);

int C2F(ilazlr)(int *m, int *n, dcomplex *a, int *lda);

int C2F(iparmq)(int *ispec, char *name__, char *opts, int *n, int 
                *ilo, int *ihi, int *lwork);

int C2F(izmax1)(int *n, dcomplex *cx, int *incx);

int C2F(lsamen)(int *n, char *ca, char *cb);

int C2F(smaxloc)(double *a, int *dimm);

extern int C2F(sbdsdc)(char *uplo, char *compq, int *n, double *d__, 
                       double *e, double *u, int *ldu, double *vt, int *ldvt, double *q, 
                       int *iq, double *work, int *iwork, int *info);

extern int C2F(sbdsqr)(char *uplo, int *n, int *ncvt, int *
                       nru, int *ncc, double *d__, double *e, double *vt, int *ldvt, double *
                       u, int *ldu, double *c__, int *ldc, double *work, int *info);

double C2F(scsum1)(int *n, dcomplex *cx, int *incx);

extern int C2F(sdisna)(char *job, int *m, int *n, double *d__, 
                       double *sep, int *info);

extern int C2F(sgbbrd)(char *vect, int *m, int *n, int *ncc, 
                       int *kl, int *ku, double *ab, int *ldab, double *d__, double *
                       e, double *q, int *ldq, double *pt, int *ldpt, double *c__, int 
                       *ldc, double *work, int *info);

extern int C2F(sgbcon)(char *norm, int *n, int *kl, int *ku, 
                       double *ab, int *ldab, int *ipiv, double *anorm, double *rcond, 
                       double *work, int *iwork, int *info);

extern int C2F(sgbequ)(int *m, int *n, int *kl, int *ku, 
                       double *ab, int *ldab, double *r__, double *c__, double *rowcnd, double *
                       colcnd, double *amax, int *info);

extern int C2F(sgbequb)(int *m, int *n, int *kl, int *
                        ku, double *ab, int *ldab, double *r__, double *c__, double *rowcnd, double 
                        *colcnd, double *amax, int *info);

extern int C2F(sgbrfs)(char *trans, int *n, int *kl, int *
                       ku, int *nrhs, double *ab, int *ldab, double *afb, int *ldafb, 
                       int *ipiv, double *b, int *ldb, double *x, int *ldx, double *
                       ferr, double *berr, double *work, int *iwork, int *info);

extern int C2F(sgbrfsx)(char *trans, char *equed, int *n, int *
                        kl, int *ku, int *nrhs, double *ab, int *ldab, double *afb, 
                        int *ldafb, int *ipiv, double *r__, double *c__, double *b, int 
                        *ldb, double *x, int *ldx, double *rcond, double *berr, int *
                        n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__, int *
                        nparams, double *params, double *work, int *iwork, int *info);

extern int C2F(sgbsv)(int *n, int *kl, int *ku, int *
                      nrhs, double *ab, int *ldab, int *ipiv, double *b, int *ldb, 
                      int *info);

extern int C2F(sgbsvx)(char *fact, char *trans, int *n, int *kl, 
                       int *ku, int *nrhs, double *ab, int *ldab, double *afb, 
                       int *ldafb, int *ipiv, char *equed, double *r__, double *c__, 
                       double *b, int *ldb, double *x, int *ldx, double *rcond, double *ferr, 
                       double *berr, double *work, int *iwork, int *info);

extern int C2F(sgbsvxx)(char *fact, char *trans, int *n, int *
                        kl, int *ku, int *nrhs, double *ab, int *ldab, double *afb, 
                        int *ldafb, int *ipiv, char *equed, double *r__, double *c__, 
                        double *b, int *ldb, double *x, int *ldx, double *rcond, double *
                        rpvgrw, double *berr, int *n_err_bnds__, double *err_bnds_norm__, 
                        double *err_bnds_comp__, int *nparams, double *params, double *work, 
                        int *iwork, int *info);

extern int C2F(sgbtf2)(int *m, int *n, int *kl, int *ku, 
                       double *ab, int *ldab, int *ipiv, int *info);

extern int C2F(sgbtrf)(int *m, int *n, int *kl, int *ku, 
                       double *ab, int *ldab, int *ipiv, int *info);

extern int C2F(sgbtrs)(char *trans, int *n, int *kl, int *
                       ku, int *nrhs, double *ab, int *ldab, int *ipiv, double *b, 
                       int *ldb, int *info);

extern int C2F(sgebak)(char *job, char *side, int *n, int *ilo, 
                       int *ihi, double *scale, int *m, double *v, int *ldv, int 
                       *info);

extern int C2F(sgebal)(char *job, int *n, double *a, int *lda, 
                       int *ilo, int *ihi, double *scale, int *info);

extern int C2F(sgebd2)(int *m, int *n, double *a, int *lda, 
                       double *d__, double *e, double *tauq, double *taup, double *work, int *info);

extern int C2F(sgebrd)(int *m, int *n, double *a, int *lda, 
                       double *d__, double *e, double *tauq, double *taup, double *work, int *
                       lwork, int *info);

extern int C2F(sgecon)(char *norm, int *n, double *a, int *lda, 
                       double *anorm, double *rcond, double *work, int *iwork, int *info);

extern int C2F(sgeequ)(int *m, int *n, double *a, int *lda, 
                       double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, int 
                       *info);

extern int C2F(sgeequb)(int *m, int *n, double *a, int *lda, 
                        double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, int 
                        *info);

extern int C2F(sgees)(char *jobvs, char *sort, L_fp select, int *n, 
                      double *a, int *lda, int *sdim, double *wr, double *wi, double *vs, 
                      int *ldvs, double *work, int *lwork, int *bwork, int *
                      info);

extern int C2F(sgeesx)(char *jobvs, char *sort, L_fp select, char *
                       sense, int *n, double *a, int *lda, int *sdim, double *wr, 
                       double *wi, double *vs, int *ldvs, double *rconde, double *rcondv, double *
                       work, int *lwork, int *iwork, int *liwork, int *bwork, 
                       int *info);

extern int C2F(sgeev)(char *jobvl, char *jobvr, int *n, double *a, 
                      int *lda, double *wr, double *wi, double *vl, int *ldvl, double *vr, 
                      int *ldvr, double *work, int *lwork, int *info);

extern int C2F(sgeevx)(char *balanc, char *jobvl, char *jobvr, char *
                       sense, int *n, double *a, int *lda, double *wr, double *wi, double *
                       vl, int *ldvl, double *vr, int *ldvr, int *ilo, int *
                       ihi, double *scale, double *abnrm, double *rconde, double *rcondv, double *work, 
                       int *lwork, int *iwork, int *info);

extern int C2F(sgegs)(char *jobvsl, char *jobvsr, int *n, double *a, 
                      int *lda, double *b, int *ldb, double *alphar, double *alphai, double 
                      *beta, double *vsl, int *ldvsl, double *vsr, int *ldvsr, double *
                      work, int *lwork, int *info);

extern int C2F(sgegv)(char *jobvl, char *jobvr, int *n, double *a, 
                      int *lda, double *b, int *ldb, double *alphar, double *alphai, double 
                      *beta, double *vl, int *ldvl, double *vr, int *ldvr, double *work, 
                      int *lwork, int *info);

extern int C2F(sgehd2)(int *n, int *ilo, int *ihi, double *a, 
                       int *lda, double *tau, double *work, int *info);

extern int C2F(sgehrd)(int *n, int *ilo, int *ihi, double *a, 
                       int *lda, double *tau, double *work, int *lwork, int *info);

extern int C2F(sgejsv)(char *joba, char *jobu, char *jobv, char *jobr, 
                       char *jobt, char *jobp, int *m, int *n, double *a, int *lda, 
                       double *sva, double *u, int *ldu, double *v, int *ldv, double *work, 
                       int *lwork, int *iwork, int *info);

extern int C2F(sgelq2)(int *m, int *n, double *a, int *lda, 
                       double *tau, double *work, int *info);

extern int C2F(sgelqf)(int *m, int *n, double *a, int *lda, 
                       double *tau, double *work, int *lwork, int *info);

extern int C2F(sgels)(char *trans, int *m, int *n, int *
                      nrhs, double *a, int *lda, double *b, int *ldb, double *work, 
                      int *lwork, int *info);

extern int C2F(sgelsd)(int *m, int *n, int *nrhs, double *a, 
                       int *lda, double *b, int *ldb, double *s, double *rcond, int *
                       rank, double *work, int *lwork, int *iwork, int *info);

extern int C2F(sgelss)(int *m, int *n, int *nrhs, double *a, 
                       int *lda, double *b, int *ldb, double *s, double *rcond, int *
                       rank, double *work, int *lwork, int *info);

extern int C2F(sgelsx)(int *m, int *n, int *nrhs, double *a, 
                       int *lda, double *b, int *ldb, int *jpvt, double *rcond, 
                       int *rank, double *work, int *info);

extern int C2F(sgelsy)(int *m, int *n, int *nrhs, double *a, 
                       int *lda, double *b, int *ldb, int *jpvt, double *rcond, 
                       int *rank, double *work, int *lwork, int *info);

extern int C2F(sgeql2)(int *m, int *n, double *a, int *lda, 
                       double *tau, double *work, int *info);

extern int C2F(sgeqlf)(int *m, int *n, double *a, int *lda, 
                       double *tau, double *work, int *lwork, int *info);

extern int C2F(sgeqp3)(int *m, int *n, double *a, int *lda, 
                       int *jpvt, double *tau, double *work, int *lwork, int *info);

extern int C2F(sgeqpf)(int *m, int *n, double *a, int *lda, 
                       int *jpvt, double *tau, double *work, int *info);

extern int C2F(sgeqr2)(int *m, int *n, double *a, int *lda, 
                       double *tau, double *work, int *info);

extern int C2F(sgeqrf)(int *m, int *n, double *a, int *lda, 
                       double *tau, double *work, int *lwork, int *info);

extern int C2F(sgerfs)(char *trans, int *n, int *nrhs, double *a, 
                       int *lda, double *af, int *ldaf, int *ipiv, double *b, 
                       int *ldb, double *x, int *ldx, double *ferr, double *berr, double *
                       work, int *iwork, int *info);

extern int C2F(sgerfsx)(char *trans, char *equed, int *n, int *
                        nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv, 
                        double *r__, double *c__, double *b, int *ldb, double *x, int *ldx, 
                        double *rcond, double *berr, int *n_err_bnds__, double *err_bnds_norm__, 
                        double *err_bnds_comp__, int *nparams, double *params, double *work, 
                        int *iwork, int *info);

extern int C2F(sgerq2)(int *m, int *n, double *a, int *lda, 
                       double *tau, double *work, int *info);

extern int C2F(sgerqf)(int *m, int *n, double *a, int *lda, 
                       double *tau, double *work, int *lwork, int *info);

extern int C2F(sgesc2)(int *n, double *a, int *lda, double *rhs, 
                       int *ipiv, int *jpiv, double *scale);

extern int C2F(sgesdd)(char *jobz, int *m, int *n, double *a, 
                       int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt, 
                       double *work, int *lwork, int *iwork, int *info);

extern int C2F(sgesv)(int *n, int *nrhs, double *a, int *lda, 
                      int *ipiv, double *b, int *ldb, int *info);

extern int C2F(sgesvd)(char *jobu, char *jobvt, int *m, int *n, 
                       double *a, int *lda, double *s, double *u, int *ldu, double *vt, 
                       int *ldvt, double *work, int *lwork, int *info);

extern int C2F(sgesvj)(char *joba, char *jobu, char *jobv, int *m, 
                       int *n, double *a, int *lda, double *sva, int *mv, double *v, 
                       int *ldv, double *work, int *lwork, int *info);

extern int C2F(sgesvx)(char *fact, char *trans, int *n, int *
                       nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv, 
                       char *equed, double *r__, double *c__, double *b, int *ldb, double *x, 
                       int *ldx, double *rcond, double *ferr, double *berr, double *work, 
                       int *iwork, int *info);

extern int C2F(sgesvxx)(char *fact, char *trans, int *n, int *
                        nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv, 
                        char *equed, double *r__, double *c__, double *b, int *ldb, double *x, 
                        int *ldx, double *rcond, double *rpvgrw, double *berr, int *
                        n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__, int *
                        nparams, double *params, double *work, int *iwork, int *info);

extern int C2F(sgetc2)(int *n, double *a, int *lda, int *ipiv, 
                       int *jpiv, int *info);

extern int C2F(sgetf2)(int *m, int *n, double *a, int *lda, 
                       int *ipiv, int *info);

extern int C2F(sgetrf)(int *m, int *n, double *a, int *lda, 
                       int *ipiv, int *info);

extern int C2F(sgetri)(int *n, double *a, int *lda, int *ipiv, 
                       double *work, int *lwork, int *info);

extern int C2F(sgetrs)(char *trans, int *n, int *nrhs, double *a, 
                       int *lda, int *ipiv, double *b, int *ldb, int *info);

extern int C2F(sggbak)(char *job, char *side, int *n, int *ilo, 
                       int *ihi, double *lscale, double *rscale, int *m, double *v, 
                       int *ldv, int *info);

extern int C2F(sggbal)(char *job, int *n, double *a, int *lda, 
                       double *b, int *ldb, int *ilo, int *ihi, double *lscale, double 
                       *rscale, double *work, int *info);

extern int C2F(sgges)(char *jobvsl, char *jobvsr, char *sort, L_fp 
                      selctg, int *n, double *a, int *lda, double *b, int *ldb, 
                      int *sdim, double *alphar, double *alphai, double *beta, double *vsl, 
                      int *ldvsl, double *vsr, int *ldvsr, double *work, int *lwork, 
                      int *bwork, int *info);

extern int C2F(sggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp 
                       selctg, char *sense, int *n, double *a, int *lda, double *b, 
                       int *ldb, int *sdim, double *alphar, double *alphai, double *beta, 
                       double *vsl, int *ldvsl, double *vsr, int *ldvsr, double *rconde, 
                       double *rcondv, double *work, int *lwork, int *iwork, int *
                       liwork, int *bwork, int *info);

extern int C2F(sggev)(char *jobvl, char *jobvr, int *n, double *a, 
                      int *lda, double *b, int *ldb, double *alphar, double *alphai, double 
                      *beta, double *vl, int *ldvl, double *vr, int *ldvr, double *work, 
                      int *lwork, int *info);

extern int C2F(sggevx)(char *balanc, char *jobvl, char *jobvr, char *
                       sense, int *n, double *a, int *lda, double *b, int *ldb, double 
                       *alphar, double *alphai, double *beta, double *vl, int *ldvl, double *vr, 
                       int *ldvr, int *ilo, int *ihi, double *lscale, double *rscale, 
                       double *abnrm, double *bbnrm, double *rconde, double *rcondv, double *work, 
                       int *lwork, int *iwork, int *bwork, int *info);

extern int C2F(sggglm)(int *n, int *m, int *p, double *a, 
                       int *lda, double *b, int *ldb, double *d__, double *x, double *y, 
                       double *work, int *lwork, int *info);

extern int C2F(sgghrd)(char *compq, char *compz, int *n, int *
                       ilo, int *ihi, double *a, int *lda, double *b, int *ldb, double 
                       *q, int *ldq, double *z__, int *ldz, int *info);

extern int C2F(sgglse)(int *m, int *n, int *p, double *a, 
                       int *lda, double *b, int *ldb, double *c__, double *d__, double *x, 
                       double *work, int *lwork, int *info);

extern int C2F(sggqrf)(int *n, int *m, int *p, double *a, 
                       int *lda, double *taua, double *b, int *ldb, double *taub, double *
                       work, int *lwork, int *info);

extern int C2F(sggrqf)(int *m, int *p, int *n, double *a, 
                       int *lda, double *taua, double *b, int *ldb, double *taub, double *
                       work, int *lwork, int *info);

extern int C2F(sggsvd)(char *jobu, char *jobv, char *jobq, int *m, 
                       int *n, int *p, int *k, int *l, double *a, int *lda, 
                       double *b, int *ldb, double *alpha, double *beta, double *u, int *
                       ldu, double *v, int *ldv, double *q, int *ldq, double *work, 
                       int *iwork, int *info);

extern int C2F(sggsvp)(char *jobu, char *jobv, char *jobq, int *m, 
                       int *p, int *n, double *a, int *lda, double *b, int *ldb, 
                       double *tola, double *tolb, int *k, int *l, double *u, int *ldu, 
                       double *v, int *ldv, double *q, int *ldq, int *iwork, double *
                       tau, double *work, int *info);

extern int C2F(sgsvj0)(char *jobv, int *m, int *n, double *a, 
                       int *lda, double *d__, double *sva, int *mv, double *v, int *
                       ldv, double *eps, double *sfmin, double *tol, int *nsweep, double *work, 
                       int *lwork, int *info);

extern int C2F(sgsvj1)(char *jobv, int *m, int *n, int *n1, 
                       double *a, int *lda, double *d__, double *sva, int *mv, double *v, 
                       int *ldv, double *eps, double *sfmin, double *tol, int *nsweep, 
                       double *work, int *lwork, int *info);

extern int C2F(sgtcon)(char *norm, int *n, double *dl, double *d__, 
                       double *du, double *du2, int *ipiv, double *anorm, double *rcond, double *
                       work, int *iwork, int *info);

extern int C2F(sgtrfs)(char *trans, int *n, int *nrhs, double *dl, 
                       double *d__, double *du, double *dlf, double *df, double *duf, double *du2, 
                       int *ipiv, double *b, int *ldb, double *x, int *ldx, double *
                       ferr, double *berr, double *work, int *iwork, int *info);

extern int C2F(sgtsv)(int *n, int *nrhs, double *dl, double *d__, 
                      double *du, double *b, int *ldb, int *info);

extern int C2F(sgtsvx)(char *fact, char *trans, int *n, int *
                       nrhs, double *dl, double *d__, double *du, double *dlf, double *df, double *duf, 
                       double *du2, int *ipiv, double *b, int *ldb, double *x, int *
                       ldx, double *rcond, double *ferr, double *berr, double *work, int *iwork, 
                       int *info);

extern int C2F(sgttrf)(int *n, double *dl, double *d__, double *du, double *
                       du2, int *ipiv, int *info);

extern int C2F(sgttrs)(char *trans, int *n, int *nrhs, double *dl, 
                       double *d__, double *du, double *du2, int *ipiv, double *b, int *ldb, 
                       int *info);

extern int C2F(sgtts2)(int *itrans, int *n, int *nrhs, double 
                       *dl, double *d__, double *du, double *du2, int *ipiv, double *b, int *
                       ldb);

extern int C2F(shgeqz)(char *job, char *compq, char *compz, int *n, 
                       int *ilo, int *ihi, double *h__, int *ldh, double *t, int 
                       *ldt, double *alphar, double *alphai, double *beta, double *q, int *ldq, 
                       double *z__, int *ldz, double *work, int *lwork, int *info);

extern int C2F(shsein)(char *side, char *eigsrc, char *initv, int *
                       select, int *n, double *h__, int *ldh, double *wr, double *wi, double 
                       *vl, int *ldvl, double *vr, int *ldvr, int *mm, int *m, 
                       double *work, int *ifaill, int *ifailr, int *info);

extern int C2F(shseqr)(char *job, char *compz, int *n, int *ilo, 
                       int *ihi, double *h__, int *ldh, double *wr, double *wi, double *z__, 
                       int *ldz, double *work, int *lwork, int *info);

int C2F(sisnan)(double *sin__);

extern int C2F(sla_gbamv_)(int *trans, int *m, int *n, 
                           int *kl, int *ku, double *alpha, double *ab, int *ldab, double *
                           x, int *incx, double *beta, double *y, int *incy);

double C2F(sla_gbrcond_)(char *trans, int *n, int *kl, int *ku, 
                         double *ab, int *ldab, double *afb, int *ldafb, int *ipiv, 
                         int *cmode, double *c__, int *info, double *work, int *iwork, 
                         int trans_len);

extern int C2F(sla_gbrfsx_extended_)(int *prec_type__, int *
                                     trans_type__, int *n, int *kl, int *ku, int *nrhs, 
                                     double *ab, int *ldab, double *afb, int *ldafb, int *ipiv, 
                                     int *colequ, double *c__, double *b, int *ldb, double *y, int *
                                     ldy, double *berr_out__, int *n_norms__, double *errs_n__, double *
                                     errs_c__, double *res, double *ayb, double *dy, double *y_tail__, double *rcond,
                                     int *ithresh, double *rthresh, double *dz_ub__, int *
                                     ignore_cwise__, int *info);

double C2F(sla_gbrpvgrw_)(int *n, int *kl, int *ku, int *
                          ncols, double *ab, int *ldab, double *afb, int *ldafb);

extern int C2F(sla_geamv_)(int *trans, int *m, int *n, double 
                           *alpha, double *a, int *lda, double *x, int *incx, double *beta, 
                           double *y, int *incy);

double C2F(sla_gercond_)(char *trans, int *n, double *a, int *lda, double 
                         *af, int *ldaf, int *ipiv, int *cmode, double *c__, int 
                         *info, double *work, int *iwork, int trans_len);

extern int C2F(sla_gerfsx_extended_)(int *prec_type__, int *
                                     trans_type__, int *n, int *nrhs, double *a, int *lda, double *
                                     af, int *ldaf, int *ipiv, int *colequ, double *c__, double *b,
                                     int *ldb, double *y, int *ldy, double *berr_out__, int *
                                     n_norms__, double *errs_n__, double *errs_c__, double *res, double *ayb, double 
                                     *dy, double *y_tail__, double *rcond, int *ithresh, double *rthresh, 
                                     double *dz_ub__, int *ignore_cwise__, int *info);

extern int C2F(sla_lin_berr_)(int *n, int *nz, int *nrhs, 
                              double *res, double *ayb, double *berr);

double C2F(sla_porcond_)(char *uplo, int *n, double *a, int *lda, double *
                         af, int *ldaf, int *cmode, double *c__, int *info, double *
                         work, int *iwork, int uplo_len);

extern int C2F(sla_porfsx_extended_)(int *prec_type__, char *uplo, 
                                     int *n, int *nrhs, double *a, int *lda, double *af, int *
                                     ldaf, int *colequ, double *c__, double *b, int *ldb, double *y, 
                                     int *ldy, double *berr_out__, int *n_norms__, double *errs_n__, 
                                     double *errs_c__, double *res, double *ayb, double *dy, double *y_tail__, double *
                                     rcond, int *ithresh, double *rthresh, double *dz_ub__, int *
                                     ignore_cwise__, int *info, int uplo_len);

double C2F(sla_porpvgrw_)(char *uplo, int *ncols, double *a, int *lda, 
                          double *af, int *ldaf, double *work, int uplo_len);

double C2F(sla_rpvgrw_)(int *n, int *ncols, double *a, int *lda, 
                        double *af, int *ldaf);

extern int C2F(sla_syamv_)(int *uplo, int *n, double *alpha, double 
                           *a, int *lda, double *x, int *incx, double *beta, double *y, 
                           int *incy);

double C2F(sla_syrcond_)(char *uplo, int *n, double *a, int *lda, double *
                         af, int *ldaf, int *ipiv, int *cmode, double *c__, int *
                         info, double *work, int *iwork, int uplo_len);

extern int C2F(sla_syrfsx_extended_)(int *prec_type__, char *uplo, 
                                     int *n, int *nrhs, double *a, int *lda, double *af, int *
                                     ldaf, int *ipiv, int *colequ, double *c__, double *b, int *
                                     ldb, double *y, int *ldy, double *berr_out__, int *n_norms__, 
                                     double *errs_n__, double *errs_c__, double *res, double *ayb, double *dy, double *
                                     y_tail__, double *rcond, int *ithresh, double *rthresh, double *dz_ub__,
                                     int *ignore_cwise__, int *info, int uplo_len);

double C2F(sla_syrpvgrw_)(char *uplo, int *n, int *info, double *a, 
                          int *lda, double *af, int *ldaf, int *ipiv, double *work, 
                          int uplo_len);

extern int C2F(sla_wwaddw_)(int *n, double *x, double *y, double *w);

extern int C2F(slabad)(double *small, double *large);

extern int C2F(slabrd)(int *m, int *n, int *nb, double *a, 
                       int *lda, double *d__, double *e, double *tauq, double *taup, double *x, 
                       int *ldx, double *y, int *ldy);

extern int C2F(slacn2)(int *n, double *v, double *x, int *isgn, 
                       double *est, int *kase, int *isave);

extern int C2F(slacon)(int *n, double *v, double *x, int *isgn, 
                       double *est, int *kase);

extern int C2F(slacpy)(char *uplo, int *m, int *n, double *a, 
                       int *lda, double *b, int *ldb);

extern int C2F(sladiv)(double *a, double *b, double *c__, double *d__, double *p, 
                       double *q);

extern int C2F(slae2)(double *a, double *b, double *c__, double *rt1, double *rt2);

extern int C2F(slaebz)(int *ijob, int *nitmax, int *n, 
                       int *mmax, int *minp, int *nbmin, double *abstol, double *
                       reltol, double *pivmin, double *d__, double *e, double *e2, int *nval, 
                       double *ab, double *c__, int *mout, int *nab, double *work, int 
                       *iwork, int *info);

extern int C2F(slaed0)(int *icompq, int *qsiz, int *n, double 
                       *d__, double *e, double *q, int *ldq, double *qstore, int *ldqs, 
                       double *work, int *iwork, int *info);

extern int C2F(slaed1)(int *n, double *d__, double *q, int *ldq, 
                       int *indxq, double *rho, int *cutpnt, double *work, int *
                       iwork, int *info);

extern int C2F(slaed2)(int *k, int *n, int *n1, double *d__, 
                       double *q, int *ldq, int *indxq, double *rho, double *z__, double *
                       dlamda, double *w, double *q2, int *indx, int *indxc, int *
                       indxp, int *coltyp, int *info);

extern int C2F(slaed3)(int *k, int *n, int *n1, double *d__, 
                       double *q, int *ldq, double *rho, double *dlamda, double *q2, int *
                       indx, int *ctot, double *w, double *s, int *info);

extern int C2F(slaed4)(int *n, int *i__, double *d__, double *z__, 
                       double *delta, double *rho, double *dlam, int *info);

extern int C2F(slaed5)(int *i__, double *d__, double *z__, double *delta, 
                       double *rho, double *dlam);

extern int C2F(slaed6)(int *kniter, int *orgati, double *rho, 
                       double *d__, double *z__, double *finit, double *tau, int *info);

extern int C2F(slaed7)(int *icompq, int *n, int *qsiz, 
                       int *tlvls, int *curlvl, int *curpbm, double *d__, double *q, 
                       int *ldq, int *indxq, double *rho, int *cutpnt, double *
                       qstore, int *qptr, int *prmptr, int *perm, int *
                       givptr, int *givcol, double *givnum, double *work, int *iwork, 
                       int *info);

extern int C2F(slaed8)(int *icompq, int *k, int *n, int 
                       *qsiz, double *d__, double *q, int *ldq, int *indxq, double *rho, 
                       int *cutpnt, double *z__, double *dlamda, double *q2, int *ldq2, 
                       double *w, int *perm, int *givptr, int *givcol, double *
                       givnum, int *indxp, int *indx, int *info);

extern int C2F(slaed9)(int *k, int *kstart, int *kstop, 
                       int *n, double *d__, double *q, int *ldq, double *rho, double *dlamda, 
                       double *w, double *s, int *lds, int *info);

extern int C2F(slaeda)(int *n, int *tlvls, int *curlvl, 
                       int *curpbm, int *prmptr, int *perm, int *givptr, 
                       int *givcol, double *givnum, double *q, int *qptr, double *z__, 
                       double *ztemp, int *info);

extern int C2F(slaein)(int *rightv, int *noinit, int *n, 
                       double *h__, int *ldh, double *wr, double *wi, double *vr, double *vi, double 
                       *b, int *ldb, double *work, double *eps3, double *smlnum, double *bignum, 
                       int *info);

extern int C2F(slaev2)(double *a, double *b, double *c__, double *rt1, double *
                       rt2, double *cs1, double *sn1);

extern int C2F(slaexc)(int *wantq, int *n, double *t, int *
                       ldt, double *q, int *ldq, int *j1, int *n1, int *n2, 
                       double *work, int *info);

extern int C2F(slag2)(double *a, int *lda, double *b, int *ldb, 
                      double *safmin, double *scale1, double *scale2, double *wr1, double *wr2, double *
                      wi);

extern int C2F(slag2d)(int *m, int *n, double *sa, int *ldsa, 
                       double *a, int *lda, int *info);

extern int C2F(slags2)(int *upper, double *a1, double *a2, double *a3, 
                       double *b1, double *b2, double *b3, double *csu, double *snu, double *csv, double *
                       snv, double *csq, double *snq);

extern int C2F(slagtf)(int *n, double *a, double *lambda, double *b, double 
                       *c__, double *tol, double *d__, int *in, int *info);

extern int C2F(slagtm)(char *trans, int *n, int *nrhs, double *
                       alpha, double *dl, double *d__, double *du, double *x, int *ldx, double *
                       beta, double *b, int *ldb);

extern int C2F(slagts)(int *job, int *n, double *a, double *b, double 
                       *c__, double *d__, int *in, double *y, double *tol, int *info);

extern int C2F(slagv2)(double *a, int *lda, double *b, int *ldb, 
                       double *alphar, double *alphai, double *beta, double *csl, double *snl, double *
                       csr, double *snr);

extern int C2F(slahqr)(int *wantt, int *wantz, int *n, 
                       int *ilo, int *ihi, double *h__, int *ldh, double *wr, double *
                       wi, int *iloz, int *ihiz, double *z__, int *ldz, int *
                       info);

extern int C2F(slahr2)(int *n, int *k, int *nb, double *a, 
                       int *lda, double *tau, double *t, int *ldt, double *y, int *ldy);

extern int C2F(slahrd)(int *n, int *k, int *nb, double *a, 
                       int *lda, double *tau, double *t, int *ldt, double *y, int *ldy);

extern int C2F(slaic1)(int *job, int *j, double *x, double *sest, 
                       double *w, double *gamma, double *sestpr, double *s, double *c__);

int C2F(slaisnan)(double *sin1, double *sin2);

extern int C2F(slaln2)(int *ltrans, int *na, int *nw, double *
                       smin, double *ca, double *a, int *lda, double *d1, double *d2, double *b, 
                       int *ldb, double *wr, double *wi, double *x, int *ldx, double *scale, 
                       double *xnorm, int *info);

extern int C2F(slals0)(int *icompq, int *nl, int *nr, 
                       int *sqre, int *nrhs, double *b, int *ldb, double *bx, 
                       int *ldbx, int *perm, int *givptr, int *givcol, 
                       int *ldgcol, double *givnum, int *ldgnum, double *poles, double *
                       difl, double *difr, double *z__, int *k, double *c__, double *s, double *
                       work, int *info);

extern int C2F(slalsa)(int *icompq, int *smlsiz, int *n, 
                       int *nrhs, double *b, int *ldb, double *bx, int *ldbx, double *
                       u, int *ldu, double *vt, int *k, double *difl, double *difr, double *
                       z__, double *poles, int *givptr, int *givcol, int *ldgcol, 
                       int *perm, double *givnum, double *c__, double *s, double *work, int *
                       iwork, int *info);

extern int C2F(slalsd)(char *uplo, int *smlsiz, int *n, int 
                       *nrhs, double *d__, double *e, double *b, int *ldb, double *rcond, 
                       int *rank, double *work, int *iwork, int *info);

extern int C2F(slamrg)(int *n1, int *n2, double *a, int *
                       strd1, int *strd2, int *index);

int C2F(slaneg)(int *n, double *d__, double *lld, double *sigma, double *pivmin, 
                int *r__);

double C2F(slangb)(char *norm, int *n, int *kl, int *ku, double *ab, 
                   int *ldab, double *work);

double C2F(slange)(char *norm, int *m, int *n, double *a, int *lda, 
                   double *work);

double C2F(slangt)(char *norm, int *n, double *dl, double *d__, double *du);

double C2F(slanhs)(char *norm, int *n, double *a, int *lda, double *work);

double C2F(slansb)(char *norm, char *uplo, int *n, int *k, double *ab, 
                   int *ldab, double *work);

double C2F(slansf)(char *norm, char *transr, char *uplo, int *n, double *a, 
                   double *work);

double C2F(slansp)(char *norm, char *uplo, int *n, double *ap, double *work);

double C2F(slanst)(char *norm, int *n, double *d__, double *e);

double C2F(slansy)(char *norm, char *uplo, int *n, double *a, int *lda, 
                   double *work);

double C2F(slantb)(char *norm, char *uplo, char *diag, int *n, int *k, 
                   double *ab, int *ldab, double *work);

double C2F(slantp)(char *norm, char *uplo, char *diag, int *n, double *ap, 
                   double *work);

double C2F(slantr)(char *norm, char *uplo, char *diag, int *m, int *n, 
                   double *a, int *lda, double *work);

extern int C2F(slanv2)(double *a, double *b, double *c__, double *d__, double *
                       rt1r, double *rt1i, double *rt2r, double *rt2i, double *cs, double *sn);

extern int C2F(slapll)(int *n, double *x, int *incx, double *y, 
                       int *incy, double *ssmin);

extern int C2F(slapmt)(int *forwrd, int *m, int *n, double *x, 
                       int *ldx, int *k);

double C2F(slapy2)(double *x, double *y);

double C2F(slapy3)(double *x, double *y, double *z__);

extern int C2F(slaqgb)(int *m, int *n, int *kl, int *ku, 
                       double *ab, int *ldab, double *r__, double *c__, double *rowcnd, double *
                       colcnd, double *amax, char *equed);

extern int C2F(slaqge)(int *m, int *n, double *a, int *lda, 
                       double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, char *
                       equed);

extern int C2F(slaqp2)(int *m, int *n, int *offset, double *a, 
                       int *lda, int *jpvt, double *tau, double *vn1, double *vn2, double *
                       work);

extern int C2F(slaqps)(int *m, int *n, int *offset, int 
                       *nb, int *kb, double *a, int *lda, int *jpvt, double *tau, 
                       double *vn1, double *vn2, double *auxv, double *f, int *ldf);

extern int C2F(slaqr0)(int *wantt, int *wantz, int *n, 
                       int *ilo, int *ihi, double *h__, int *ldh, double *wr, double *
                       wi, int *iloz, int *ihiz, double *z__, int *ldz, double *work, 
                       int *lwork, int *info);

extern int C2F(slaqr1)(int *n, double *h__, int *ldh, double *sr1, 
                       double *si1, double *sr2, double *si2, double *v);

extern int C2F(slaqr2)(int *wantt, int *wantz, int *n, 
                       int *ktop, int *kbot, int *nw, double *h__, int *ldh, 
                       int *iloz, int *ihiz, double *z__, int *ldz, int *ns, 
                       int *nd, double *sr, double *si, double *v, int *ldv, int *nh, 
                       double *t, int *ldt, int *nv, double *wv, int *ldwv, double *
                       work, int *lwork);

extern int C2F(slaqr3)(int *wantt, int *wantz, int *n, 
                       int *ktop, int *kbot, int *nw, double *h__, int *ldh, 
                       int *iloz, int *ihiz, double *z__, int *ldz, int *ns, 
                       int *nd, double *sr, double *si, double *v, int *ldv, int *nh, 
                       double *t, int *ldt, int *nv, double *wv, int *ldwv, double *
                       work, int *lwork);

extern int C2F(slaqr4)(int *wantt, int *wantz, int *n, 
                       int *ilo, int *ihi, double *h__, int *ldh, double *wr, double *
                       wi, int *iloz, int *ihiz, double *z__, int *ldz, double *work, 
                       int *lwork, int *info);

extern int C2F(slaqr5)(int *wantt, int *wantz, int *kacc22, 
                       int *n, int *ktop, int *kbot, int *nshfts, double *sr, 
                       double *si, double *h__, int *ldh, int *iloz, int *ihiz, double 
                       *z__, int *ldz, double *v, int *ldv, double *u, int *ldu, 
                       int *nv, double *wv, int *ldwv, int *nh, double *wh, int *
                       ldwh);

extern int C2F(slaqsb)(char *uplo, int *n, int *kd, double *ab, 
                       int *ldab, double *s, double *scond, double *amax, char *equed);

extern int C2F(slaqsp)(char *uplo, int *n, double *ap, double *s, double *
                       scond, double *amax, char *equed);

extern int C2F(slaqsy)(char *uplo, int *n, double *a, int *lda, 
                       double *s, double *scond, double *amax, char *equed);

extern int C2F(slaqtr)(int *ltran, int *ldouble, int *n, double 
                       *t, int *ldt, double *b, double *w, double *scale, double *x, double *work, 
                       int *info);

extern int C2F(slar1v)(int *n, int *b1, int *bn, double *
                       lambda, double *d__, double *l, double *ld, double *lld, double *pivmin, double *
                       gaptol, double *z__, int *wantnc, int *negcnt, double *ztz, double *
                       mingma, int *r__, int *isuppz, double *nrminv, double *resid, 
                       double *rqcorr, double *work);

extern int C2F(slar2v)(int *n, double *x, double *y, double *z__, int 
                       *incx, double *c__, double *s, int *incc);

extern int C2F(slarf)(char *side, int *m, int *n, double *v, 
                      int *incv, double *tau, double *c__, int *ldc, double *work);

extern int C2F(slarfb)(char *side, char *trans, char *direct, char *
                       storev, int *m, int *n, int *k, double *v, int *ldv, 
                       double *t, int *ldt, double *c__, int *ldc, double *work, int *
                       ldwork);

extern int C2F(slarfg)(int *n, double *alpha, double *x, int *incx, 
                       double *tau);

extern int C2F(slarfp)(int *n, double *alpha, double *x, int *incx, 
                       double *tau);

extern int C2F(slarft)(char *direct, char *storev, int *n, int *
                       k, double *v, int *ldv, double *tau, double *t, int *ldt);

extern int C2F(slarfx)(char *side, int *m, int *n, double *v, 
                       double *tau, double *c__, int *ldc, double *work);

extern int C2F(slargv)(int *n, double *x, int *incx, double *y, 
                       int *incy, double *c__, int *incc);

extern int C2F(slarnv)(int *idist, int *iseed, int *n, double 
                       *x);

extern int C2F(slarra)(int *n, double *d__, double *e, double *e2, double *
                       spltol, double *tnrm, int *nsplit, int *isplit, int *info);

extern int C2F(slarrb)(int *n, double *d__, double *lld, int *
                       ifirst, int *ilast, double *rtol1, double *rtol2, int *offset, 
                       double *w, double *wgap, double *werr, double *work, int *iwork, double *
                       pivmin, double *spdiam, int *twist, int *info);

extern int C2F(slarrc)(char *jobt, int *n, double *vl, double *vu, double 
                       *d__, double *e, double *pivmin, int *eigcnt, int *lcnt, int *
                       rcnt, int *info);

extern int C2F(slarrd)(char *range, char *order, int *n, double *vl, 
                       double *vu, int *il, int *iu, double *gers, double *reltol, double *
                       d__, double *e, double *e2, double *pivmin, int *nsplit, int *
                       isplit, int *m, double *w, double *werr, double *wl, double *wu, int *
                       iblock, int *indexw, double *work, int *iwork, int *info);

extern int C2F(slarre)(char *range, int *n, double *vl, double *vu, 
                       int *il, int *iu, double *d__, double *e, double *e2, double *rtol1, 
                       double *rtol2, double *spltol, int *nsplit, int *isplit, int *
                       m, double *w, double *werr, double *wgap, int *iblock, int *indexw, 
                       double *gers, double *pivmin, double *work, int *iwork, int *info);

extern int C2F(slarrf)(int *n, double *d__, double *l, double *ld, 
                       int *clstrt, int *clend, double *w, double *wgap, double *werr, 
                       double *spdiam, double *clgapl, double *clgapr, double *pivmin, double *sigma, 
                       double *dplus, double *lplus, double *work, int *info);

extern int C2F(slarrj)(int *n, double *d__, double *e2, int *ifirst, 
                       int *ilast, double *rtol, int *offset, double *w, double *werr, 
                       double *work, int *iwork, double *pivmin, double *spdiam, int *info);

extern int C2F(slarrk)(int *n, int *iw, double *gl, double *gu, 
                       double *d__, double *e2, double *pivmin, double *reltol, double *w, double *werr, 
                       int *info);

extern int C2F(slarrr)(int *n, double *d__, double *e, int *info);

extern int C2F(slarrv)(int *n, double *vl, double *vu, double *d__, double *
                       l, double *pivmin, int *isplit, int *m, int *dol, int *
                       dou, double *minrgp, double *rtol1, double *rtol2, double *w, double *werr, 
                       double *wgap, int *iblock, int *indexw, double *gers, double *z__, 
                       int *ldz, int *isuppz, double *work, int *iwork, int *
                       info);

extern int C2F(slarscl2)(int *m, int *n, double *d__, double *x, 
                         int *ldx);

extern int C2F(slartg)(double *f, double *g, double *cs, double *sn, double *r__);

extern int C2F(slartv)(int *n, double *x, int *incx, double *y, 
                       int *incy, double *c__, double *s, int *incc);

extern int C2F(slaruv)(int *iseed, int *n, double *x);

extern int C2F(slarz)(char *side, int *m, int *n, int *l, 
                      double *v, int *incv, double *tau, double *c__, int *ldc, double *
                      work);

extern int C2F(slarzb)(char *side, char *trans, char *direct, char *
                       storev, int *m, int *n, int *k, int *l, double *v, 
                       int *ldv, double *t, int *ldt, double *c__, int *ldc, double *
                       work, int *ldwork);

extern int C2F(slarzt)(char *direct, char *storev, int *n, int *
                       k, double *v, int *ldv, double *tau, double *t, int *ldt);

extern int C2F(slas2)(double *f, double *g, double *h__, double *ssmin, double *
                      ssmax);

extern int C2F(slascl)(char *type__, int *kl, int *ku, double *
                       cfrom, double *cto, int *m, int *n, double *a, int *lda, 
                       int *info);

extern int C2F(slascl2)(int *m, int *n, double *d__, double *x, 
                        int *ldx);

extern int C2F(slasd0)(int *n, int *sqre, double *d__, double *e, 
                       double *u, int *ldu, double *vt, int *ldvt, int *smlsiz, 
                       int *iwork, double *work, int *info);

extern int C2F(slasd1)(int *nl, int *nr, int *sqre, double *
                       d__, double *alpha, double *beta, double *u, int *ldu, double *vt, 
                       int *ldvt, int *idxq, int *iwork, double *work, int *
                       info);

extern int C2F(slasd2)(int *nl, int *nr, int *sqre, int 
                       *k, double *d__, double *z__, double *alpha, double *beta, double *u, int *
                       ldu, double *vt, int *ldvt, double *dsigma, double *u2, int *ldu2, 
                       double *vt2, int *ldvt2, int *idxp, int *idx, int *idxc, 
                       int *idxq, int *coltyp, int *info);

extern int C2F(slasd3)(int *nl, int *nr, int *sqre, int 
                       *k, double *d__, double *q, int *ldq, double *dsigma, double *u, int *
                       ldu, double *u2, int *ldu2, double *vt, int *ldvt, double *vt2, 
                       int *ldvt2, int *idxc, int *ctot, double *z__, int *
                       info);

extern int C2F(slasd4)(int *n, int *i__, double *d__, double *z__, 
                       double *delta, double *rho, double *sigma, double *work, int *info);

extern int C2F(slasd5)(int *i__, double *d__, double *z__, double *delta, 
                       double *rho, double *dsigma, double *work);

extern int C2F(slasd6)(int *icompq, int *nl, int *nr, 
                       int *sqre, double *d__, double *vf, double *vl, double *alpha, double *beta, 
                       int *idxq, int *perm, int *givptr, int *givcol, 
                       int *ldgcol, double *givnum, int *ldgnum, double *poles, double *
                       difl, double *difr, double *z__, int *k, double *c__, double *s, double *
                       work, int *iwork, int *info);

extern int C2F(slasd7)(int *icompq, int *nl, int *nr, 
                       int *sqre, int *k, double *d__, double *z__, double *zw, double *vf, 
                       double *vfw, double *vl, double *vlw, double *alpha, double *beta, double *dsigma, 
                       int *idx, int *idxp, int *idxq, int *perm, int *
                       givptr, int *givcol, int *ldgcol, double *givnum, int *
                       ldgnum, double *c__, double *s, int *info);

extern int C2F(slasd8)(int *icompq, int *k, double *d__, double *
                       z__, double *vf, double *vl, double *difl, double *difr, int *lddifr, 
                       double *dsigma, double *work, int *info);

extern int C2F(slasda)(int *icompq, int *smlsiz, int *n, 
                       int *sqre, double *d__, double *e, double *u, int *ldu, double *vt, 
                       int *k, double *difl, double *difr, double *z__, double *poles, int *
                       givptr, int *givcol, int *ldgcol, int *perm, double *givnum, 
                       double *c__, double *s, double *work, int *iwork, int *info);

extern int C2F(slasdq)(char *uplo, int *sqre, int *n, int *
                       ncvt, int *nru, int *ncc, double *d__, double *e, double *vt, 
                       int *ldvt, double *u, int *ldu, double *c__, int *ldc, double *
                       work, int *info);

extern int C2F(slasdt)(int *n, int *lvl, int *nd, int *
                       inode, int *ndiml, int *ndimr, int *msub);

extern int C2F(slaset)(char *uplo, int *m, int *n, double *alpha, 
                       double *beta, double *a, int *lda);

extern int C2F(slasq1)(int *n, double *d__, double *e, double *work, 
                       int *info);

extern int C2F(slasq2)(int *n, double *z__, int *info);

extern int C2F(slasq3)(int *i0, int *n0, double *z__, int *pp, 
                       double *dmin__, double *sigma, double *desig, double *qmax, int *nfail, 
                       int *iter, int *ndiv, int *ieee, int *ttype, double *
                       dmin1, double *dmin2, double *dn, double *dn1, double *dn2, double *g, double *
                       tau);

extern int C2F(slasq4)(int *i0, int *n0, double *z__, int *pp, 
                       int *n0in, double *dmin__, double *dmin1, double *dmin2, double *dn, 
                       double *dn1, double *dn2, double *tau, int *ttype, double *g);

extern int C2F(slasq5)(int *i0, int *n0, double *z__, int *pp, 
                       double *tau, double *dmin__, double *dmin1, double *dmin2, double *dn, double *
                       dnm1, double *dnm2, int *ieee);

extern int C2F(slasq6)(int *i0, int *n0, double *z__, int *pp, 
                       double *dmin__, double *dmin1, double *dmin2, double *dn, double *dnm1, double *
                       dnm2);

extern int C2F(slasr)(char *side, char *pivot, char *direct, int *m, 
                      int *n, double *c__, double *s, double *a, int *lda);

extern int C2F(slasrt)(char *id, int *n, double *d__, int *info);

extern int C2F(slassq)(int *n, double *x, int *incx, double *scale, 
                       double *sumsq);

extern int C2F(slasv2)(double *f, double *g, double *h__, double *ssmin, double *
                       ssmax, double *snr, double *csr, double *snl, double *csl);

extern int C2F(slaswp)(int *n, double *a, int *lda, int *k1, 
                       int *k2, int *ipiv, int *incx);

extern int C2F(slasy2)(int *ltranl, int *ltranr, int *isgn, 
                       int *n1, int *n2, double *tl, int *ldtl, double *tr, int *
                       ldtr, double *b, int *ldb, double *scale, double *x, int *ldx, double 
                       *xnorm, int *info);

extern int C2F(slasyf)(char *uplo, int *n, int *nb, int *kb, 
                       double *a, int *lda, int *ipiv, double *w, int *ldw, int 
                       *info);

extern int C2F(slatbs)(char *uplo, char *trans, char *diag, char *
                       normin, int *n, int *kd, double *ab, int *ldab, double *x, 
                       double *scale, double *cnorm, int *info);

extern int C2F(slatdf)(int *ijob, int *n, double *z__, int *
                       ldz, double *rhs, double *rdsum, double *rdscal, int *ipiv, int *
                       jpiv);

extern int C2F(slatps)(char *uplo, char *trans, char *diag, char *
                       normin, int *n, double *ap, double *x, double *scale, double *cnorm, 
                       int *info);

extern int C2F(slatrd)(char *uplo, int *n, int *nb, double *a, 
                       int *lda, double *e, double *tau, double *w, int *ldw);

extern int C2F(slatrs)(char *uplo, char *trans, char *diag, char *
                       normin, int *n, double *a, int *lda, double *x, double *scale, double 
                       *cnorm, int *info);

extern int C2F(slatrz)(int *m, int *n, int *l, double *a, 
                       int *lda, double *tau, double *work);

extern int C2F(slatzm)(char *side, int *m, int *n, double *v, 
                       int *incv, double *tau, double *c1, double *c2, int *ldc, double *
                       work);

extern int C2F(slauu2)(char *uplo, int *n, double *a, int *lda, 
                       int *info);

extern int C2F(slauum)(char *uplo, int *n, double *a, int *lda, 
                       int *info);

extern int C2F(sopgtr)(char *uplo, int *n, double *ap, double *tau, 
                       double *q, int *ldq, double *work, int *info);

extern int C2F(sopmtr)(char *side, char *uplo, char *trans, int *m, 
                       int *n, double *ap, double *tau, double *c__, int *ldc, double *work, 
                       int *info);

extern int C2F(sorg2l)(int *m, int *n, int *k, double *a, 
                       int *lda, double *tau, double *work, int *info);

extern int C2F(sorg2r)(int *m, int *n, int *k, double *a, 
                       int *lda, double *tau, double *work, int *info);

extern int C2F(sorgbr)(char *vect, int *m, int *n, int *k, 
                       double *a, int *lda, double *tau, double *work, int *lwork, int 
                       *info);

extern int C2F(sorghr)(int *n, int *ilo, int *ihi, double *a, 
                       int *lda, double *tau, double *work, int *lwork, int *info);

extern int C2F(sorgl2)(int *m, int *n, int *k, double *a, 
                       int *lda, double *tau, double *work, int *info);

extern int C2F(sorglq)(int *m, int *n, int *k, double *a, 
                       int *lda, double *tau, double *work, int *lwork, int *info);

extern int C2F(sorgql)(int *m, int *n, int *k, double *a, 
                       int *lda, double *tau, double *work, int *lwork, int *info);

extern int C2F(sorgqr)(int *m, int *n, int *k, double *a, 
                       int *lda, double *tau, double *work, int *lwork, int *info);

extern int C2F(sorgr2)(int *m, int *n, int *k, double *a, 
                       int *lda, double *tau, double *work, int *info);

extern int C2F(sorgrq)(int *m, int *n, int *k, double *a, 
                       int *lda, double *tau, double *work, int *lwork, int *info);

extern int C2F(sorgtr)(char *uplo, int *n, double *a, int *lda, 
                       double *tau, double *work, int *lwork, int *info);

extern int C2F(sorm2l)(char *side, char *trans, int *m, int *n, 
                       int *k, double *a, int *lda, double *tau, double *c__, int *ldc, 
                       double *work, int *info);

extern int C2F(sorm2r)(char *side, char *trans, int *m, int *n, 
                       int *k, double *a, int *lda, double *tau, double *c__, int *ldc, 
                       double *work, int *info);

extern int C2F(sormbr)(char *vect, char *side, char *trans, int *m, 
                       int *n, int *k, double *a, int *lda, double *tau, double *c__, 
                       int *ldc, double *work, int *lwork, int *info);

extern int C2F(sormhr)(char *side, char *trans, int *m, int *n, 
                       int *ilo, int *ihi, double *a, int *lda, double *tau, double *
                       c__, int *ldc, double *work, int *lwork, int *info);

extern int C2F(sorml2)(char *side, char *trans, int *m, int *n, 
                       int *k, double *a, int *lda, double *tau, double *c__, int *ldc, 
                       double *work, int *info);

extern int C2F(sormlq)(char *side, char *trans, int *m, int *n, 
                       int *k, double *a, int *lda, double *tau, double *c__, int *ldc, 
                       double *work, int *lwork, int *info);

extern int C2F(sormql)(char *side, char *trans, int *m, int *n, 
                       int *k, double *a, int *lda, double *tau, double *c__, int *ldc, 
                       double *work, int *lwork, int *info);

extern int C2F(sormqr)(char *side, char *trans, int *m, int *n, 
                       int *k, double *a, int *lda, double *tau, double *c__, int *ldc, 
                       double *work, int *lwork, int *info);

extern int C2F(sormr2)(char *side, char *trans, int *m, int *n, 
                       int *k, double *a, int *lda, double *tau, double *c__, int *ldc, 
                       double *work, int *info);

extern int C2F(sormr3)(char *side, char *trans, int *m, int *n, 
                       int *k, int *l, double *a, int *lda, double *tau, double *c__, 
                       int *ldc, double *work, int *info);

extern int C2F(sormrq)(char *side, char *trans, int *m, int *n, 
                       int *k, double *a, int *lda, double *tau, double *c__, int *ldc, 
                       double *work, int *lwork, int *info);

extern int C2F(sormrz)(char *side, char *trans, int *m, int *n, 
                       int *k, int *l, double *a, int *lda, double *tau, double *c__, 
                       int *ldc, double *work, int *lwork, int *info);

extern int C2F(sormtr)(char *side, char *uplo, char *trans, int *m, 
                       int *n, double *a, int *lda, double *tau, double *c__, int *ldc, 
                       double *work, int *lwork, int *info);

extern int C2F(spbcon)(char *uplo, int *n, int *kd, double *ab, 
                       int *ldab, double *anorm, double *rcond, double *work, int *iwork, 
                       int *info);

extern int C2F(spbequ)(char *uplo, int *n, int *kd, double *ab, 
                       int *ldab, double *s, double *scond, double *amax, int *info);

extern int C2F(spbrfs)(char *uplo, int *n, int *kd, int *
                       nrhs, double *ab, int *ldab, double *afb, int *ldafb, double *b, 
                       int *ldb, double *x, int *ldx, double *ferr, double *berr, double *
                       work, int *iwork, int *info);

extern int C2F(spbstf)(char *uplo, int *n, int *kd, double *ab, 
                       int *ldab, int *info);

extern int C2F(spbsv)(char *uplo, int *n, int *kd, int *
                      nrhs, double *ab, int *ldab, double *b, int *ldb, int *info);

extern int C2F(spbsvx)(char *fact, char *uplo, int *n, int *kd, 
                       int *nrhs, double *ab, int *ldab, double *afb, int *ldafb, 
                       char *equed, double *s, double *b, int *ldb, double *x, int *ldx, 
                       double *rcond, double *ferr, double *berr, double *work, int *iwork, 
                       int *info);

extern int C2F(spbtf2)(char *uplo, int *n, int *kd, double *ab, 
                       int *ldab, int *info);

extern int C2F(spbtrf)(char *uplo, int *n, int *kd, double *ab, 
                       int *ldab, int *info);

extern int C2F(spbtrs)(char *uplo, int *n, int *kd, int *
                       nrhs, double *ab, int *ldab, double *b, int *ldb, int *info);

extern int C2F(spftrf)(char *transr, char *uplo, int *n, double *a, 
                       int *info);

extern int C2F(spftri)(char *transr, char *uplo, int *n, double *a, 
                       int *info);

extern int C2F(spftrs)(char *transr, char *uplo, int *n, int *
                       nrhs, double *a, double *b, int *ldb, int *info);

extern int C2F(spocon)(char *uplo, int *n, double *a, int *lda, 
                       double *anorm, double *rcond, double *work, int *iwork, int *info);

extern int C2F(spoequ)(int *n, double *a, int *lda, double *s, double 
                       *scond, double *amax, int *info);

extern int C2F(spoequb)(int *n, double *a, int *lda, double *s, 
                        double *scond, double *amax, int *info);

extern int C2F(sporfs)(char *uplo, int *n, int *nrhs, double *a, 
                       int *lda, double *af, int *ldaf, double *b, int *ldb, double *x, 
                       int *ldx, double *ferr, double *berr, double *work, int *iwork, 
                       int *info);

extern int C2F(sporfsx)(char *uplo, char *equed, int *n, int *
                        nrhs, double *a, int *lda, double *af, int *ldaf, double *s, double *
                        b, int *ldb, double *x, int *ldx, double *rcond, double *berr, 
                        int *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__, 
                        int *nparams, double *params, double *work, int *iwork, int *
                        info);

extern int C2F(sposv)(char *uplo, int *n, int *nrhs, double *a, 
                      int *lda, double *b, int *ldb, int *info);

extern int C2F(sposvx)(char *fact, char *uplo, int *n, int *
                       nrhs, double *a, int *lda, double *af, int *ldaf, char *equed, 
                       double *s, double *b, int *ldb, double *x, int *ldx, double *rcond, 
                       double *ferr, double *berr, double *work, int *iwork, int *info);

extern int C2F(sposvxx)(char *fact, char *uplo, int *n, int *
                        nrhs, double *a, int *lda, double *af, int *ldaf, char *equed, 
                        double *s, double *b, int *ldb, double *x, int *ldx, double *rcond, 
                        double *rpvgrw, double *berr, int *n_err_bnds__, double *
                        err_bnds_norm__, double *err_bnds_comp__, int *nparams, double *
                        params, double *work, int *iwork, int *info);

extern int C2F(spotf2)(char *uplo, int *n, double *a, int *lda, 
                       int *info);

extern int C2F(spotrf)(char *uplo, int *n, double *a, int *lda, 
                       int *info);

extern int C2F(spotri)(char *uplo, int *n, double *a, int *lda, 
                       int *info);

extern int C2F(spotrs)(char *uplo, int *n, int *nrhs, double *a, 
                       int *lda, double *b, int *ldb, int *info);

extern int C2F(sppcon)(char *uplo, int *n, double *ap, double *anorm, 
                       double *rcond, double *work, int *iwork, int *info);

extern int C2F(sppequ)(char *uplo, int *n, double *ap, double *s, double *
                       scond, double *amax, int *info);

extern int C2F(spprfs)(char *uplo, int *n, int *nrhs, double *ap, 
                       double *afp, double *b, int *ldb, double *x, int *ldx, double *ferr, 
                       double *berr, double *work, int *iwork, int *info);

extern int C2F(sppsv)(char *uplo, int *n, int *nrhs, double *ap, 
                      double *b, int *ldb, int *info);

extern int C2F(sppsvx)(char *fact, char *uplo, int *n, int *
                       nrhs, double *ap, double *afp, char *equed, double *s, double *b, int *
                       ldb, double *x, int *ldx, double *rcond, double *ferr, double *berr, double 
                       *work, int *iwork, int *info);

extern int C2F(spptrf)(char *uplo, int *n, double *ap, int *info);

extern int C2F(spptri)(char *uplo, int *n, double *ap, int *info);

extern int C2F(spptrs)(char *uplo, int *n, int *nrhs, double *ap, 
                       double *b, int *ldb, int *info);

extern int C2F(spstf2)(char *uplo, int *n, double *a, int *lda, 
                       int *piv, int *rank, double *tol, double *work, int *info);

extern int C2F(spstrf)(char *uplo, int *n, double *a, int *lda, 
                       int *piv, int *rank, double *tol, double *work, int *info);

extern int C2F(sptcon)(int *n, double *d__, double *e, double *anorm, 
                       double *rcond, double *work, int *info);

extern int C2F(spteqr)(char *compz, int *n, double *d__, double *e, 
                       double *z__, int *ldz, double *work, int *info);

extern int C2F(sptrfs)(int *n, int *nrhs, double *d__, double *e, 
                       double *df, double *ef, double *b, int *ldb, double *x, int *ldx, 
                       double *ferr, double *berr, double *work, int *info);

extern int C2F(sptsv)(int *n, int *nrhs, double *d__, double *e, 
                      double *b, int *ldb, int *info);

extern int C2F(sptsvx)(char *fact, int *n, int *nrhs, double *d__, 
                       double *e, double *df, double *ef, double *b, int *ldb, double *x, int 
                       *ldx, double *rcond, double *ferr, double *berr, double *work, int *info);

extern int C2F(spttrf)(int *n, double *d__, double *e, int *info);

extern int C2F(spttrs)(int *n, int *nrhs, double *d__, double *e, 
                       double *b, int *ldb, int *info);

extern int C2F(sptts2)(int *n, int *nrhs, double *d__, double *e, 
                       double *b, int *ldb);

extern int C2F(srscl)(int *n, double *sa, double *sx, int *incx);

extern int C2F(ssbev)(char *jobz, char *uplo, int *n, int *kd, 
                      double *ab, int *ldab, double *w, double *z__, int *ldz, double *work, 
                      int *info);

extern int C2F(ssbevd)(char *jobz, char *uplo, int *n, int *kd, 
                       double *ab, int *ldab, double *w, double *z__, int *ldz, double *work, 
                       int *lwork, int *iwork, int *liwork, int *info);

extern int C2F(ssbevx)(char *jobz, char *range, char *uplo, int *n, 
                       int *kd, double *ab, int *ldab, double *q, int *ldq, double *vl, 
                       double *vu, int *il, int *iu, double *abstol, int *m, double *
                       w, double *z__, int *ldz, double *work, int *iwork, int *
                       ifail, int *info);

extern int C2F(ssbgst)(char *vect, char *uplo, int *n, int *ka, 
                       int *kb, double *ab, int *ldab, double *bb, int *ldbb, double *
                       x, int *ldx, double *work, int *info);

extern int C2F(ssbgv)(char *jobz, char *uplo, int *n, int *ka, 
                      int *kb, double *ab, int *ldab, double *bb, int *ldbb, double *
                      w, double *z__, int *ldz, double *work, int *info);

extern int C2F(ssbgvd)(char *jobz, char *uplo, int *n, int *ka, 
                       int *kb, double *ab, int *ldab, double *bb, int *ldbb, double *
                       w, double *z__, int *ldz, double *work, int *lwork, int *
                       iwork, int *liwork, int *info);

extern int C2F(ssbgvx)(char *jobz, char *range, char *uplo, int *n, 
                       int *ka, int *kb, double *ab, int *ldab, double *bb, int *
                       ldbb, double *q, int *ldq, double *vl, double *vu, int *il, int 
                       *iu, double *abstol, int *m, double *w, double *z__, int *ldz, double 
                       *work, int *iwork, int *ifail, int *info);

extern int C2F(ssbtrd)(char *vect, char *uplo, int *n, int *kd, 
                       double *ab, int *ldab, double *d__, double *e, double *q, int *ldq, 
                       double *work, int *info);

extern int C2F(ssfrk)(char *transr, char *uplo, char *trans, int *n, 
                      int *k, double *alpha, double *a, int *lda, double *beta, double *
                      c__);

extern int C2F(sspcon)(char *uplo, int *n, double *ap, int *ipiv, 
                       double *anorm, double *rcond, double *work, int *iwork, int *info);

extern int C2F(sspev)(char *jobz, char *uplo, int *n, double *ap, 
                      double *w, double *z__, int *ldz, double *work, int *info);

extern int C2F(sspevd)(char *jobz, char *uplo, int *n, double *ap, 
                       double *w, double *z__, int *ldz, double *work, int *lwork, int 
                       *iwork, int *liwork, int *info);

extern int C2F(sspevx)(char *jobz, char *range, char *uplo, int *n, 
                       double *ap, double *vl, double *vu, int *il, int *iu, double *abstol, 
                       int *m, double *w, double *z__, int *ldz, double *work, int *
                       iwork, int *ifail, int *info);

extern int C2F(sspgst)(int *itype, char *uplo, int *n, double *ap, 
                       double *bp, int *info);

extern int C2F(sspgv)(int *itype, char *jobz, char *uplo, int *
                      n, double *ap, double *bp, double *w, double *z__, int *ldz, double *work, 
                      int *info);

extern int C2F(sspgvd)(int *itype, char *jobz, char *uplo, int *
                       n, double *ap, double *bp, double *w, double *z__, int *ldz, double *work, 
                       int *lwork, int *iwork, int *liwork, int *info);

extern int C2F(sspgvx)(int *itype, char *jobz, char *range, char *
                       uplo, int *n, double *ap, double *bp, double *vl, double *vu, int *il, 
                       int *iu, double *abstol, int *m, double *w, double *z__, int *
                       ldz, double *work, int *iwork, int *ifail, int *info);

extern int C2F(ssprfs)(char *uplo, int *n, int *nrhs, double *ap, 
                       double *afp, int *ipiv, double *b, int *ldb, double *x, int *
                       ldx, double *ferr, double *berr, double *work, int *iwork, int *
                       info);

extern int C2F(sspsv)(char *uplo, int *n, int *nrhs, double *ap, 
                      int *ipiv, double *b, int *ldb, int *info);

extern int C2F(sspsvx)(char *fact, char *uplo, int *n, int *
                       nrhs, double *ap, double *afp, int *ipiv, double *b, int *ldb, double 
                       *x, int *ldx, double *rcond, double *ferr, double *berr, double *work, 
                       int *iwork, int *info);

extern int C2F(ssptrd)(char *uplo, int *n, double *ap, double *d__, 
                       double *e, double *tau, int *info);

extern int C2F(ssptrf)(char *uplo, int *n, double *ap, int *ipiv, 
                       int *info);

extern int C2F(ssptri)(char *uplo, int *n, double *ap, int *ipiv, 
                       double *work, int *info);

extern int C2F(ssptrs)(char *uplo, int *n, int *nrhs, double *ap, 
                       int *ipiv, double *b, int *ldb, int *info);

extern int C2F(sstebz)(char *range, char *order, int *n, double *vl, 
                       double *vu, int *il, int *iu, double *abstol, double *d__, double *e, 
                       int *m, int *nsplit, double *w, int *iblock, int *
                       isplit, double *work, int *iwork, int *info);

extern int C2F(sstedc)(char *compz, int *n, double *d__, double *e, 
                       double *z__, int *ldz, double *work, int *lwork, int *iwork, 
                       int *liwork, int *info);

extern int C2F(sstegr)(char *jobz, char *range, int *n, double *d__, 
                       double *e, double *vl, double *vu, int *il, int *iu, double *abstol, 
                       int *m, double *w, double *z__, int *ldz, int *isuppz, double *
                       work, int *lwork, int *iwork, int *liwork, int *info);

extern int C2F(sstein)(int *n, double *d__, double *e, int *m, double 
                       *w, int *iblock, int *isplit, double *z__, int *ldz, double *
                       work, int *iwork, int *ifail, int *info);

extern int C2F(sstemr)(char *jobz, char *range, int *n, double *d__, 
                       double *e, double *vl, double *vu, int *il, int *iu, int *m, 
                       double *w, double *z__, int *ldz, int *nzc, int *isuppz, 
                       int *tryrac, double *work, int *lwork, int *iwork, int *
                       liwork, int *info);

extern int C2F(ssteqr)(char *compz, int *n, double *d__, double *e, 
                       double *z__, int *ldz, double *work, int *info);

extern int C2F(ssterf)(int *n, double *d__, double *e, int *info);

extern int C2F(sstev)(char *jobz, int *n, double *d__, double *e, double *
                      z__, int *ldz, double *work, int *info);

extern int C2F(sstevd)(char *jobz, int *n, double *d__, double *e, double 
                       *z__, int *ldz, double *work, int *lwork, int *iwork, 
                       int *liwork, int *info);

extern int C2F(sstevr)(char *jobz, char *range, int *n, double *d__, 
                       double *e, double *vl, double *vu, int *il, int *iu, double *abstol, 
                       int *m, double *w, double *z__, int *ldz, int *isuppz, double *
                       work, int *lwork, int *iwork, int *liwork, int *info);

extern int C2F(sstevx)(char *jobz, char *range, int *n, double *d__, 
                       double *e, double *vl, double *vu, int *il, int *iu, double *abstol, 
                       int *m, double *w, double *z__, int *ldz, double *work, int *
                       iwork, int *ifail, int *info);

extern int C2F(ssycon)(char *uplo, int *n, double *a, int *lda, 
                       int *ipiv, double *anorm, double *rcond, double *work, int *iwork, 
                       int *info);

extern int C2F(ssyequb)(char *uplo, int *n, double *a, int *lda, 
                        double *s, double *scond, double *amax, double *work, int *info);

extern int C2F(ssyev)(char *jobz, char *uplo, int *n, double *a, 
                      int *lda, double *w, double *work, int *lwork, int *info);

extern int C2F(ssyevd)(char *jobz, char *uplo, int *n, double *a, 
                       int *lda, double *w, double *work, int *lwork, int *iwork, 
                       int *liwork, int *info);

extern int C2F(ssyevr)(char *jobz, char *range, char *uplo, int *n, 
                       double *a, int *lda, double *vl, double *vu, int *il, int *iu, 
                       double *abstol, int *m, double *w, double *z__, int *ldz, int *
                       isuppz, double *work, int *lwork, int *iwork, int *liwork, 
                       int *info);

extern int C2F(ssyevx)(char *jobz, char *range, char *uplo, int *n, 
                       double *a, int *lda, double *vl, double *vu, int *il, int *iu, 
                       double *abstol, int *m, double *w, double *z__, int *ldz, double *
                       work, int *lwork, int *iwork, int *ifail, int *info);

extern int C2F(ssygs2)(int *itype, char *uplo, int *n, double *a, 
                       int *lda, double *b, int *ldb, int *info);

extern int C2F(ssygst)(int *itype, char *uplo, int *n, double *a, 
                       int *lda, double *b, int *ldb, int *info);

extern int C2F(ssygv)(int *itype, char *jobz, char *uplo, int *
                      n, double *a, int *lda, double *b, int *ldb, double *w, double *work, 
                      int *lwork, int *info);

extern int C2F(ssygvd)(int *itype, char *jobz, char *uplo, int *
                       n, double *a, int *lda, double *b, int *ldb, double *w, double *work, 
                       int *lwork, int *iwork, int *liwork, int *info);

extern int C2F(ssygvx)(int *itype, char *jobz, char *range, char *
                       uplo, int *n, double *a, int *lda, double *b, int *ldb, double *
                       vl, double *vu, int *il, int *iu, double *abstol, int *m, 
                       double *w, double *z__, int *ldz, double *work, int *lwork, int 
                       *iwork, int *ifail, int *info);

extern int C2F(ssyrfs)(char *uplo, int *n, int *nrhs, double *a, 
                       int *lda, double *af, int *ldaf, int *ipiv, double *b, 
                       int *ldb, double *x, int *ldx, double *ferr, double *berr, double *
                       work, int *iwork, int *info);

extern int C2F(ssyrfsx)(char *uplo, char *equed, int *n, int *
                        nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv, 
                        double *s, double *b, int *ldb, double *x, int *ldx, double *rcond, 
                        double *berr, int *n_err_bnds__, double *err_bnds_norm__, double *
                        err_bnds_comp__, int *nparams, double *params, double *work, int *
                        iwork, int *info);

extern int C2F(ssysv)(char *uplo, int *n, int *nrhs, double *a, 
                      int *lda, int *ipiv, double *b, int *ldb, double *work, 
                      int *lwork, int *info);

extern int C2F(ssysvx)(char *fact, char *uplo, int *n, int *
                       nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv, 
                       double *b, int *ldb, double *x, int *ldx, double *rcond, double *ferr, 
                       double *berr, double *work, int *lwork, int *iwork, int *
                       info);

extern int C2F(ssysvxx)(char *fact, char *uplo, int *n, int *
                        nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv, 
                        char *equed, double *s, double *b, int *ldb, double *x, int *ldx, 
                        double *rcond, double *rpvgrw, double *berr, int *n_err_bnds__, double *
                        err_bnds_norm__, double *err_bnds_comp__, int *nparams, double *
                        params, double *work, int *iwork, int *info);

extern int C2F(ssytd2)(char *uplo, int *n, double *a, int *lda, 
                       double *d__, double *e, double *tau, int *info);

extern int C2F(ssytf2)(char *uplo, int *n, double *a, int *lda, 
                       int *ipiv, int *info);

extern int C2F(ssytrd)(char *uplo, int *n, double *a, int *lda, 
                       double *d__, double *e, double *tau, double *work, int *lwork, int *
                       info);

extern int C2F(ssytrf)(char *uplo, int *n, double *a, int *lda, 
                       int *ipiv, double *work, int *lwork, int *info);

extern int C2F(ssytri)(char *uplo, int *n, double *a, int *lda, 
                       int *ipiv, double *work, int *info);

extern int C2F(ssytrs)(char *uplo, int *n, int *nrhs, double *a, 
                       int *lda, int *ipiv, double *b, int *ldb, int *info);

extern int C2F(stbcon)(char *norm, char *uplo, char *diag, int *n, 
                       int *kd, double *ab, int *ldab, double *rcond, double *work, 
                       int *iwork, int *info);

extern int C2F(stbrfs)(char *uplo, char *trans, char *diag, int *n, 
                       int *kd, int *nrhs, double *ab, int *ldab, double *b, int 
                       *ldb, double *x, int *ldx, double *ferr, double *berr, double *work, 
                       int *iwork, int *info);

extern int C2F(stbtrs)(char *uplo, char *trans, char *diag, int *n, 
                       int *kd, int *nrhs, double *ab, int *ldab, double *b, int 
                       *ldb, int *info);

extern int C2F(stfsm)(char *transr, char *side, char *uplo, char *trans, 
                      char *diag, int *m, int *n, double *alpha, double *a, double *b, 
                      int *ldb);

extern int C2F(stftri)(char *transr, char *uplo, char *diag, int *n, 
                       double *a, int *info);

extern int C2F(stfttp)(char *transr, char *uplo, int *n, double *arf, 
                       double *ap, int *info);

extern int C2F(stfttr)(char *transr, char *uplo, int *n, double *arf, 
                       double *a, int *lda, int *info);

extern int C2F(stgevc)(char *side, char *howmny, int *select, 
                       int *n, double *s, int *lds, double *p, int *ldp, double *vl, 
                       int *ldvl, double *vr, int *ldvr, int *mm, int *m, double 
                       *work, int *info);

extern int C2F(stgex2)(int *wantq, int *wantz, int *n, double 
                       *a, int *lda, double *b, int *ldb, double *q, int *ldq, double *
                       z__, int *ldz, int *j1, int *n1, int *n2, double *work, 
                       int *lwork, int *info);

extern int C2F(stgexc)(int *wantq, int *wantz, int *n, double 
                       *a, int *lda, double *b, int *ldb, double *q, int *ldq, double *
                       z__, int *ldz, int *ifst, int *ilst, double *work, int *
                       lwork, int *info);

extern int C2F(stgsen)(int *ijob, int *wantq, int *wantz, 
                       int *select, int *n, double *a, int *lda, double *b, int *
                       ldb, double *alphar, double *alphai, double *beta, double *q, int *ldq, 
                       double *z__, int *ldz, int *m, double *pl, double *pr, double *dif, 
                       double *work, int *lwork, int *iwork, int *liwork, int *
                       info);

extern int C2F(stgsja)(char *jobu, char *jobv, char *jobq, int *m, 
                       int *p, int *n, int *k, int *l, double *a, int *lda, 
                       double *b, int *ldb, double *tola, double *tolb, double *alpha, double *
                       beta, double *u, int *ldu, double *v, int *ldv, double *q, int *
                       ldq, double *work, int *ncycle, int *info);

extern int C2F(stgsna)(char *job, char *howmny, int *select, 
                       int *n, double *a, int *lda, double *b, int *ldb, double *vl, 
                       int *ldvl, double *vr, int *ldvr, double *s, double *dif, int *
                       mm, int *m, double *work, int *lwork, int *iwork, int *
                       info);

extern int C2F(stgsy2)(char *trans, int *ijob, int *m, int *
                       n, double *a, int *lda, double *b, int *ldb, double *c__, int *
                       ldc, double *d__, int *ldd, double *e, int *lde, double *f, int 
                       *ldf, double *scale, double *rdsum, double *rdscal, int *iwork, int 
                       *pq, int *info);

extern int C2F(stgsyl)(char *trans, int *ijob, int *m, int *
                       n, double *a, int *lda, double *b, int *ldb, double *c__, int *
                       ldc, double *d__, int *ldd, double *e, int *lde, double *f, int 
                       *ldf, double *scale, double *dif, double *work, int *lwork, int *
                       iwork, int *info);

extern int C2F(stpcon)(char *norm, char *uplo, char *diag, int *n, 
                       double *ap, double *rcond, double *work, int *iwork, int *info);

extern int C2F(stprfs)(char *uplo, char *trans, char *diag, int *n, 
                       int *nrhs, double *ap, double *b, int *ldb, double *x, int *ldx, 
                       double *ferr, double *berr, double *work, int *iwork, int *info);

extern int C2F(stptri)(char *uplo, char *diag, int *n, double *ap, 
                       int *info);

extern int C2F(stptrs)(char *uplo, char *trans, char *diag, int *n, 
                       int *nrhs, double *ap, double *b, int *ldb, int *info);

extern int C2F(stpttf)(char *transr, char *uplo, int *n, double *ap, 
                       double *arf, int *info);

extern int C2F(stpttr)(char *uplo, int *n, double *ap, double *a, 
                       int *lda, int *info);

extern int C2F(strcon)(char *norm, char *uplo, char *diag, int *n, 
                       double *a, int *lda, double *rcond, double *work, int *iwork, 
                       int *info);

extern int C2F(strevc)(char *side, char *howmny, int *select, 
                       int *n, double *t, int *ldt, double *vl, int *ldvl, double *vr, 
                       int *ldvr, int *mm, int *m, double *work, int *info);

extern int C2F(strexc)(char *compq, int *n, double *t, int *ldt, 
                       double *q, int *ldq, int *ifst, int *ilst, double *work, 
                       int *info);

extern int C2F(strrfs)(char *uplo, char *trans, char *diag, int *n, 
                       int *nrhs, double *a, int *lda, double *b, int *ldb, double *x, 
                       int *ldx, double *ferr, double *berr, double *work, int *iwork, 
                       int *info);

extern int C2F(strsen)(char *job, char *compq, int *select, int 
                       *n, double *t, int *ldt, double *q, int *ldq, double *wr, double *wi, 
                       int *m, double *s, double *sep, double *work, int *lwork, int *
                       iwork, int *liwork, int *info);

extern int C2F(strsna)(char *job, char *howmny, int *select, 
                       int *n, double *t, int *ldt, double *vl, int *ldvl, double *vr, 
                       int *ldvr, double *s, double *sep, int *mm, int *m, double *
                       work, int *ldwork, int *iwork, int *info);

extern int C2F(strsyl)(char *trana, char *tranb, int *isgn, int 
                       *m, int *n, double *a, int *lda, double *b, int *ldb, double *
                       c__, int *ldc, double *scale, int *info);

extern int C2F(strti2)(char *uplo, char *diag, int *n, double *a, 
                       int *lda, int *info);

extern int C2F(strtri)(char *uplo, char *diag, int *n, double *a, 
                       int *lda, int *info);

extern int C2F(strtrs)(char *uplo, char *trans, char *diag, int *n, 
                       int *nrhs, double *a, int *lda, double *b, int *ldb, int *
                       info);

extern int C2F(strttf)(char *transr, char *uplo, int *n, double *a, 
                       int *lda, double *arf, int *info);

extern int C2F(strttp)(char *uplo, int *n, double *a, int *lda, 
                       double *ap, int *info);

extern int C2F(stzrqf)(int *m, int *n, double *a, int *lda, 
                       double *tau, int *info);

extern int C2F(stzrzf)(int *m, int *n, double *a, int *lda, 
                       double *tau, double *work, int *lwork, int *info);

extern int C2F(xerbla)(char *srname, int *info);

extern int C2F(zbdsqr)(char *uplo, int *n, int *ncvt, int *
                       nru, int *ncc, double *d__, double *e, dcomplex *vt, 
                       int *ldvt, dcomplex *u, int *ldu, dcomplex *c__, 
                       int *ldc, double *rwork, int *info);

extern int C2F(zcgesv)(int *n, int *nrhs, dcomplex *a, 
                       int *lda, int *ipiv, dcomplex *b, int *ldb, 
                       dcomplex *x, int *ldx, dcomplex *work, dcomplex *swork, 
                       double *rwork, int *iter, int *info);

extern int C2F(zcposv)(char *uplo, int *n, int *nrhs, 
                       dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       dcomplex *x, int *ldx, dcomplex *work, dcomplex *swork, 
                       double *rwork, int *iter, int *info);

extern int C2F(zdrscl)(int *n, double *sa, dcomplex *sx, 
                       int *incx);

extern int C2F(zgbbrd)(char *vect, int *m, int *n, int *ncc, 
                       int *kl, int *ku, dcomplex *ab, int *ldab, 
                       double *d__, double *e, dcomplex *q, int *ldq, 
                       dcomplex *pt, int *ldpt, dcomplex *c__, int *ldc, 
                       dcomplex *work, double *rwork, int *info);

extern int C2F(zgbcon)(char *norm, int *n, int *kl, int *ku, 
                       dcomplex *ab, int *ldab, int *ipiv, double *anorm, 
                       double *rcond, dcomplex *work, double *rwork, int *
                       info);

extern int C2F(zgbequ)(int *m, int *n, int *kl, int *ku, 
                       dcomplex *ab, int *ldab, double *r__, double *c__, 
                       double *rowcnd, double *colcnd, double *amax, int *
                       info);

extern int C2F(zgbequb)(int *m, int *n, int *kl, int *
                        ku, dcomplex *ab, int *ldab, double *r__, double *
                        c__, double *rowcnd, double *colcnd, double *amax, 
                        int *info);

extern int C2F(zgbrfs)(char *trans, int *n, int *kl, int *
                       ku, int *nrhs, dcomplex *ab, int *ldab, dcomplex *
                       afb, int *ldafb, int *ipiv, dcomplex *b, int *ldb, 
                       dcomplex *x, int *ldx, double *ferr, double *berr, 
                       dcomplex *work, double *rwork, int *info);

extern int C2F(zgbrfsx)(char *trans, char *equed, int *n, int *
                        kl, int *ku, int *nrhs, dcomplex *ab, int *ldab, 
                        dcomplex *afb, int *ldafb, int *ipiv, double *r__, 
                        double *c__, dcomplex *b, int *ldb, dcomplex *x, 
                        int *ldx, double *rcond, double *berr, int *
                        n_err_bnds__, double *err_bnds_norm__, double *
                        err_bnds_comp__, int *nparams, double *params, dcomplex *
                        work, double *rwork, int *info);

extern int C2F(zgbsv)(int *n, int *kl, int *ku, int *
                      nrhs, dcomplex *ab, int *ldab, int *ipiv, dcomplex *
                      b, int *ldb, int *info);

extern int C2F(zgbsvx)(char *fact, char *trans, int *n, int *kl, 
                       int *ku, int *nrhs, dcomplex *ab, int *ldab, 
                       dcomplex *afb, int *ldafb, int *ipiv, char *equed, 
                       double *r__, double *c__, dcomplex *b, int *ldb, 
                       dcomplex *x, int *ldx, double *rcond, double *ferr, 
                       double *berr, dcomplex *work, double *rwork, int *
                       info);

extern int C2F(zgbsvxx)(char *fact, char *trans, int *n, int *
                        kl, int *ku, int *nrhs, dcomplex *ab, int *ldab, 
                        dcomplex *afb, int *ldafb, int *ipiv, char *equed, 
                        double *r__, double *c__, dcomplex *b, int *ldb, 
                        dcomplex *x, int *ldx, double *rcond, double *rpvgrw, 
                        double *berr, int *n_err_bnds__, double *err_bnds_norm__, 
                        double *err_bnds_comp__, int *nparams, double *params, 
                        dcomplex *work, double *rwork, int *info);

extern int C2F(zgbtf2)(int *m, int *n, int *kl, int *ku, 
                       dcomplex *ab, int *ldab, int *ipiv, int *info);

extern int C2F(zgbtrf)(int *m, int *n, int *kl, int *ku, 
                       dcomplex *ab, int *ldab, int *ipiv, int *info);

extern int C2F(zgbtrs)(char *trans, int *n, int *kl, int *
                       ku, int *nrhs, dcomplex *ab, int *ldab, int *ipiv, 
                       dcomplex *b, int *ldb, int *info);

extern int C2F(zgebak)(char *job, char *side, int *n, int *ilo, 
                       int *ihi, double *scale, int *m, dcomplex *v, 
                       int *ldv, int *info);

extern int C2F(zgebal)(char *job, int *n, dcomplex *a, int 
                       *lda, int *ilo, int *ihi, double *scale, int *info);

extern int C2F(zgebd2)(int *m, int *n, dcomplex *a, 
                       int *lda, double *d__, double *e, dcomplex *tauq, 
                       dcomplex *taup, dcomplex *work, int *info);

extern int C2F(zgebrd)(int *m, int *n, dcomplex *a, 
                       int *lda, double *d__, double *e, dcomplex *tauq, 
                       dcomplex *taup, dcomplex *work, int *lwork, int *
                       info);

extern int C2F(zgecon)(char *norm, int *n, dcomplex *a, 
                       int *lda, double *anorm, double *rcond, dcomplex *
                       work, double *rwork, int *info);

extern int C2F(zgeequ)(int *m, int *n, dcomplex *a, 
                       int *lda, double *r__, double *c__, double *rowcnd, 
                       double *colcnd, double *amax, int *info);

extern int C2F(zgeequb)(int *m, int *n, dcomplex *a, 
                        int *lda, double *r__, double *c__, double *rowcnd, 
                        double *colcnd, double *amax, int *info);

extern int C2F(zgees)(char *jobvs, char *sort, L_fp select, int *n, 
                      dcomplex *a, int *lda, int *sdim, dcomplex *w, 
                      dcomplex *vs, int *ldvs, dcomplex *work, int *lwork, 
                      double *rwork, int *bwork, int *info);

extern int C2F(zgeesx)(char *jobvs, char *sort, L_fp select, char *
                       sense, int *n, dcomplex *a, int *lda, int *sdim, 
                       dcomplex *w, dcomplex *vs, int *ldvs, double *
                       rconde, double *rcondv, dcomplex *work, int *lwork, 
                       double *rwork, int *bwork, int *info);

extern int C2F(zgeev)(char *jobvl, char *jobvr, int *n, 
                      dcomplex *a, int *lda, dcomplex *w, dcomplex *vl, 
                      int *ldvl, dcomplex *vr, int *ldvr, dcomplex *work, 
                      int *lwork, double *rwork, int *info);

extern int C2F(zgeevx)(char *balanc, char *jobvl, char *jobvr, char *
                       sense, int *n, dcomplex *a, int *lda, dcomplex *w, 
                       dcomplex *vl, int *ldvl, dcomplex *vr, int *ldvr, 
                       int *ilo, int *ihi, double *scale, double *abnrm, 
                       double *rconde, double *rcondv, dcomplex *work, int *
                       lwork, double *rwork, int *info);

extern int C2F(zgegs)(char *jobvsl, char *jobvsr, int *n, 
                      dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                      dcomplex *alpha, dcomplex *beta, dcomplex *vsl, 
                      int *ldvsl, dcomplex *vsr, int *ldvsr, dcomplex *
                      work, int *lwork, double *rwork, int *info);

extern int C2F(zgegv)(char *jobvl, char *jobvr, int *n, 
                      dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                      dcomplex *alpha, dcomplex *beta, dcomplex *vl, int 
                      *ldvl, dcomplex *vr, int *ldvr, dcomplex *work, int 
                      *lwork, double *rwork, int *info);

extern int C2F(zgehd2)(int *n, int *ilo, int *ihi, 
                       dcomplex *a, int *lda, dcomplex *tau, dcomplex *
                       work, int *info);

extern int C2F(zgehrd)(int *n, int *ilo, int *ihi, 
                       dcomplex *a, int *lda, dcomplex *tau, dcomplex *
                       work, int *lwork, int *info);

extern int C2F(zgelq2)(int *m, int *n, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work, int *info);

extern int C2F(zgelqf)(int *m, int *n, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work, int *lwork, 
                       int *info);

extern int C2F(zgels)(char *trans, int *m, int *n, int *
                      nrhs, dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                      dcomplex *work, int *lwork, int *info);

extern int C2F(zgelsd)(int *m, int *n, int *nrhs, 
                       dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       double *s, double *rcond, int *rank, dcomplex *work, 
                       int *lwork, double *rwork, int *iwork, int *info);

extern int C2F(zgelss)(int *m, int *n, int *nrhs, 
                       dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       double *s, double *rcond, int *rank, dcomplex *work, 
                       int *lwork, double *rwork, int *info);

extern int C2F(zgelsx)(int *m, int *n, int *nrhs, 
                       dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       int *jpvt, double *rcond, int *rank, dcomplex *work, 
                       double *rwork, int *info);

extern int C2F(zgelsy)(int *m, int *n, int *nrhs, 
                       dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       int *jpvt, double *rcond, int *rank, dcomplex *work, 
                       int *lwork, double *rwork, int *info);

extern int C2F(zgeql2)(int *m, int *n, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work, int *info);

extern int C2F(zgeqlf)(int *m, int *n, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work, int *lwork, 
                       int *info);

extern int C2F(zgeqp3)(int *m, int *n, dcomplex *a, 
                       int *lda, int *jpvt, dcomplex *tau, dcomplex *work, 
                       int *lwork, double *rwork, int *info);

extern int C2F(zgeqpf)(int *m, int *n, dcomplex *a, 
                       int *lda, int *jpvt, dcomplex *tau, dcomplex *work, 
                       double *rwork, int *info);

extern int C2F(zgeqr2)(int *m, int *n, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work, int *info);

extern int C2F(zgeqrf)(int *m, int *n, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work, int *lwork, 
                       int *info);

extern int C2F(zgerfs)(char *trans, int *n, int *nrhs, 
                       dcomplex *a, int *lda, dcomplex *af, int *ldaf, 
                       int *ipiv, dcomplex *b, int *ldb, dcomplex *x, 
                       int *ldx, double *ferr, double *berr, dcomplex *work, 
                       double *rwork, int *info);

extern int C2F(zgerfsx)(char *trans, char *equed, int *n, int *
                        nrhs, dcomplex *a, int *lda, dcomplex *af, int *
                        ldaf, int *ipiv, double *r__, double *c__, dcomplex *
                        b, int *ldb, dcomplex *x, int *ldx, double *rcond, 
                        double *berr, int *n_err_bnds__, double *err_bnds_norm__, 
                        double *err_bnds_comp__, int *nparams, double *params, 
                        dcomplex *work, double *rwork, int *info);

extern int C2F(zgerq2)(int *m, int *n, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work, int *info);

extern int C2F(zgerqf)(int *m, int *n, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work, int *lwork, 
                       int *info);

extern int C2F(zgesc2)(int *n, dcomplex *a, int *lda, 
                       dcomplex *rhs, int *ipiv, int *jpiv, double *scale);

extern int C2F(zgesdd)(char *jobz, int *m, int *n, 
                       dcomplex *a, int *lda, double *s, dcomplex *u, 
                       int *ldu, dcomplex *vt, int *ldvt, dcomplex *work, 
                       int *lwork, double *rwork, int *iwork, int *info);

extern int C2F(zgesv)(int *n, int *nrhs, dcomplex *a, 
                      int *lda, int *ipiv, dcomplex *b, int *ldb, int *
                      info);

extern int C2F(zgesvd)(char *jobu, char *jobvt, int *m, int *n, 
                       dcomplex *a, int *lda, double *s, dcomplex *u, 
                       int *ldu, dcomplex *vt, int *ldvt, dcomplex *work, 
                       int *lwork, double *rwork, int *info);

extern int C2F(zgesvx)(char *fact, char *trans, int *n, int *
                       nrhs, dcomplex *a, int *lda, dcomplex *af, int *
                       ldaf, int *ipiv, char *equed, double *r__, double *c__, 
                       dcomplex *b, int *ldb, dcomplex *x, int *ldx, 
                       double *rcond, double *ferr, double *berr, dcomplex *
                       work, double *rwork, int *info);

extern int C2F(zgesvxx)(char *fact, char *trans, int *n, int *
                        nrhs, dcomplex *a, int *lda, dcomplex *af, int *
                        ldaf, int *ipiv, char *equed, double *r__, double *c__, 
                        dcomplex *b, int *ldb, dcomplex *x, int *ldx, 
                        double *rcond, double *rpvgrw, double *berr, int *
                        n_err_bnds__, double *err_bnds_norm__, double *
                        err_bnds_comp__, int *nparams, double *params, dcomplex *
                        work, double *rwork, int *info);

extern int C2F(zgetc2)(int *n, dcomplex *a, int *lda, 
                       int *ipiv, int *jpiv, int *info);

extern int C2F(zgetf2)(int *m, int *n, dcomplex *a, 
                       int *lda, int *ipiv, int *info);

extern int C2F(zgetrf)(int *m, int *n, dcomplex *a, 
                       int *lda, int *ipiv, int *info);

extern int C2F(zgetri)(int *n, dcomplex *a, int *lda, 
                       int *ipiv, dcomplex *work, int *lwork, int *info);

extern int C2F(zgetrs)(char *trans, int *n, int *nrhs, 
                       dcomplex *a, int *lda, int *ipiv, dcomplex *b, 
                       int *ldb, int *info);

extern int C2F(zggbak)(char *job, char *side, int *n, int *ilo, 
                       int *ihi, double *lscale, double *rscale, int *m, 
                       dcomplex *v, int *ldv, int *info);

extern int C2F(zggbal)(char *job, int *n, dcomplex *a, int 
                       *lda, dcomplex *b, int *ldb, int *ilo, int *ihi, 
                       double *lscale, double *rscale, double *work, int *
                       info);

extern int C2F(zgges)(char *jobvsl, char *jobvsr, char *sort, L_fp 
                      selctg, int *n, dcomplex *a, int *lda, dcomplex *b, 
                      int *ldb, int *sdim, dcomplex *alpha, dcomplex *
                      beta, dcomplex *vsl, int *ldvsl, dcomplex *vsr, int 
                      *ldvsr, dcomplex *work, int *lwork, double *rwork, 
                      int *bwork, int *info);

extern int C2F(zggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp 
                       selctg, char *sense, int *n, dcomplex *a, int *lda, 
                       dcomplex *b, int *ldb, int *sdim, dcomplex *alpha, 
                       dcomplex *beta, dcomplex *vsl, int *ldvsl, 
                       dcomplex *vsr, int *ldvsr, double *rconde, double *
                       rcondv, dcomplex *work, int *lwork, double *rwork, 
                       int *iwork, int *liwork, int *bwork, int *info);

extern int C2F(zggev)(char *jobvl, char *jobvr, int *n, 
                      dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                      dcomplex *alpha, dcomplex *beta, dcomplex *vl, int 
                      *ldvl, dcomplex *vr, int *ldvr, dcomplex *work, int 
                      *lwork, double *rwork, int *info);

extern int C2F(zggevx)(char *balanc, char *jobvl, char *jobvr, char *
                       sense, int *n, dcomplex *a, int *lda, dcomplex *b, 
                       int *ldb, dcomplex *alpha, dcomplex *beta, 
                       dcomplex *vl, int *ldvl, dcomplex *vr, int *ldvr, 
                       int *ilo, int *ihi, double *lscale, double *rscale, 
                       double *abnrm, double *bbnrm, double *rconde, double *
                       rcondv, dcomplex *work, int *lwork, double *rwork, 
                       int *iwork, int *bwork, int *info);

extern int C2F(zggglm)(int *n, int *m, int *p, 
                       dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       dcomplex *d__, dcomplex *x, dcomplex *y, dcomplex 
                       *work, int *lwork, int *info);

extern int C2F(zgghrd)(char *compq, char *compz, int *n, int *
                       ilo, int *ihi, dcomplex *a, int *lda, dcomplex *b, 
                       int *ldb, dcomplex *q, int *ldq, dcomplex *z__, 
                       int *ldz, int *info);

extern int C2F(zgglse)(int *m, int *n, int *p, 
                       dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       dcomplex *c__, dcomplex *d__, dcomplex *x, 
                       dcomplex *work, int *lwork, int *info);

extern int C2F(zggqrf)(int *n, int *m, int *p, 
                       dcomplex *a, int *lda, dcomplex *taua, dcomplex *b, 
                       int *ldb, dcomplex *taub, dcomplex *work, int *
                       lwork, int *info);

extern int C2F(zggrqf)(int *m, int *p, int *n, 
                       dcomplex *a, int *lda, dcomplex *taua, dcomplex *b, 
                       int *ldb, dcomplex *taub, dcomplex *work, int *
                       lwork, int *info);

extern int C2F(zggsvd)(char *jobu, char *jobv, char *jobq, int *m, 
                       int *n, int *p, int *k, int *l, dcomplex *a, 
                       int *lda, dcomplex *b, int *ldb, double *alpha, 
                       double *beta, dcomplex *u, int *ldu, dcomplex *v, 
                       int *ldv, dcomplex *q, int *ldq, dcomplex *work, 
                       double *rwork, int *iwork, int *info);

extern int C2F(zggsvp)(char *jobu, char *jobv, char *jobq, int *m, 
                       int *p, int *n, dcomplex *a, int *lda, dcomplex 
                       *b, int *ldb, double *tola, double *tolb, int *k, 
                       int *l, dcomplex *u, int *ldu, dcomplex *v, int 
                       *ldv, dcomplex *q, int *ldq, int *iwork, double *
                       rwork, dcomplex *tau, dcomplex *work, int *info);

extern int C2F(zgtcon)(char *norm, int *n, dcomplex *dl, 
                       dcomplex *d__, dcomplex *du, dcomplex *du2, int *
                       ipiv, double *anorm, double *rcond, dcomplex *work, 
                       int *info);

extern int C2F(zgtrfs)(char *trans, int *n, int *nrhs, 
                       dcomplex *dl, dcomplex *d__, dcomplex *du, 
                       dcomplex *dlf, dcomplex *df, dcomplex *duf, 
                       dcomplex *du2, int *ipiv, dcomplex *b, int *ldb, 
                       dcomplex *x, int *ldx, double *ferr, double *berr, 
                       dcomplex *work, double *rwork, int *info);

extern int C2F(zgtsv)(int *n, int *nrhs, dcomplex *dl, 
                      dcomplex *d__, dcomplex *du, dcomplex *b, int *ldb, 
                      int *info);

extern int C2F(zgtsvx)(char *fact, char *trans, int *n, int *
                       nrhs, dcomplex *dl, dcomplex *d__, dcomplex *du, 
                       dcomplex *dlf, dcomplex *df, dcomplex *duf, 
                       dcomplex *du2, int *ipiv, dcomplex *b, int *ldb, 
                       dcomplex *x, int *ldx, double *rcond, double *ferr, 
                       double *berr, dcomplex *work, double *rwork, int *
                       info);

extern int C2F(zgttrf)(int *n, dcomplex *dl, dcomplex *
                       d__, dcomplex *du, dcomplex *du2, int *ipiv, int *
                       info);

extern int C2F(zgttrs)(char *trans, int *n, int *nrhs, 
                       dcomplex *dl, dcomplex *d__, dcomplex *du, 
                       dcomplex *du2, int *ipiv, dcomplex *b, int *ldb, 
                       int *info);

extern int C2F(zgtts2)(int *itrans, int *n, int *nrhs, 
                       dcomplex *dl, dcomplex *d__, dcomplex *du, 
                       dcomplex *du2, int *ipiv, dcomplex *b, int *ldb);

extern int C2F(zhbev)(char *jobz, char *uplo, int *n, int *kd, 
                      dcomplex *ab, int *ldab, double *w, dcomplex *z__, 
                      int *ldz, dcomplex *work, double *rwork, int *info);

extern int C2F(zhbevd)(char *jobz, char *uplo, int *n, int *kd, 
                       dcomplex *ab, int *ldab, double *w, dcomplex *z__, 
                       int *ldz, dcomplex *work, int *lwork, double *rwork, 
                       int *lrwork, int *iwork, int *liwork, int *info);

extern int C2F(zhbevx)(char *jobz, char *range, char *uplo, int *n, 
                       int *kd, dcomplex *ab, int *ldab, dcomplex *q, 
                       int *ldq, double *vl, double *vu, int *il, int *
                       iu, double *abstol, int *m, double *w, dcomplex *z__, 
                       int *ldz, dcomplex *work, double *rwork, int *iwork, 
                       int *ifail, int *info);

extern int C2F(zhbgst)(char *vect, char *uplo, int *n, int *ka, 
                       int *kb, dcomplex *ab, int *ldab, dcomplex *bb, 
                       int *ldbb, dcomplex *x, int *ldx, dcomplex *work, 
                       double *rwork, int *info);

extern int C2F(zhbgv)(char *jobz, char *uplo, int *n, int *ka, 
                      int *kb, dcomplex *ab, int *ldab, dcomplex *bb, 
                      int *ldbb, double *w, dcomplex *z__, int *ldz, 
                      dcomplex *work, double *rwork, int *info);

extern int C2F(zhbgvd)(char *jobz, char *uplo, int *n, int *ka, 
                       int *kb, dcomplex *ab, int *ldab, dcomplex *bb, 
                       int *ldbb, double *w, dcomplex *z__, int *ldz, 
                       dcomplex *work, int *lwork, double *rwork, int *
                       lrwork, int *iwork, int *liwork, int *info);

extern int C2F(zhbgvx)(char *jobz, char *range, char *uplo, int *n, 
                       int *ka, int *kb, dcomplex *ab, int *ldab, 
                       dcomplex *bb, int *ldbb, dcomplex *q, int *ldq, 
                       double *vl, double *vu, int *il, int *iu, double *
                       abstol, int *m, double *w, dcomplex *z__, int *ldz, 
                       dcomplex *work, double *rwork, int *iwork, int *
                       ifail, int *info);

extern int C2F(zhbtrd)(char *vect, char *uplo, int *n, int *kd, 
                       dcomplex *ab, int *ldab, double *d__, double *e, 
                       dcomplex *q, int *ldq, dcomplex *work, int *info);

extern int C2F(zhecon)(char *uplo, int *n, dcomplex *a, 
                       int *lda, int *ipiv, double *anorm, double *rcond, 
                       dcomplex *work, int *info);

extern int C2F(zheequb)(char *uplo, int *n, dcomplex *a, 
                        int *lda, double *s, double *scond, double *amax, 
                        dcomplex *work, int *info);

extern int C2F(zheev)(char *jobz, char *uplo, int *n, dcomplex 
                      *a, int *lda, double *w, dcomplex *work, int *lwork, 
                      double *rwork, int *info);

extern int C2F(zheevd)(char *jobz, char *uplo, int *n, 
                       dcomplex *a, int *lda, double *w, dcomplex *work, 
                       int *lwork, double *rwork, int *lrwork, int *iwork, 
                       int *liwork, int *info);

extern int C2F(zheevr)(char *jobz, char *range, char *uplo, int *n, 
                       dcomplex *a, int *lda, double *vl, double *vu, 
                       int *il, int *iu, double *abstol, int *m, double *
                       w, dcomplex *z__, int *ldz, int *isuppz, dcomplex *
                       work, int *lwork, double *rwork, int *lrwork, int *
                       iwork, int *liwork, int *info);

extern int C2F(zheevx)(char *jobz, char *range, char *uplo, int *n, 
                       dcomplex *a, int *lda, double *vl, double *vu, 
                       int *il, int *iu, double *abstol, int *m, double *
                       w, dcomplex *z__, int *ldz, dcomplex *work, int *
                       lwork, double *rwork, int *iwork, int *ifail, int *
                       info);

extern int C2F(zhegs2)(int *itype, char *uplo, int *n, 
                       dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       int *info);

extern int C2F(zhegst)(int *itype, char *uplo, int *n, 
                       dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       int *info);

extern int C2F(zhegv)(int *itype, char *jobz, char *uplo, int *
                      n, dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                      double *w, dcomplex *work, int *lwork, double *rwork, 
                      int *info);

extern int C2F(zhegvd)(int *itype, char *jobz, char *uplo, int *
                       n, dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       double *w, dcomplex *work, int *lwork, double *rwork, 
                       int *lrwork, int *iwork, int *liwork, int *info);

extern int C2F(zhegvx)(int *itype, char *jobz, char *range, char *
                       uplo, int *n, dcomplex *a, int *lda, dcomplex *b, 
                       int *ldb, double *vl, double *vu, int *il, int *
                       iu, double *abstol, int *m, double *w, dcomplex *z__, 
                       int *ldz, dcomplex *work, int *lwork, double *rwork, 
                       int *iwork, int *ifail, int *info);

extern int C2F(zherfs)(char *uplo, int *n, int *nrhs, 
                       dcomplex *a, int *lda, dcomplex *af, int *ldaf, 
                       int *ipiv, dcomplex *b, int *ldb, dcomplex *x, 
                       int *ldx, double *ferr, double *berr, dcomplex *work, 
                       double *rwork, int *info);

extern int C2F(zherfsx)(char *uplo, char *equed, int *n, int *
                        nrhs, dcomplex *a, int *lda, dcomplex *af, int *
                        ldaf, int *ipiv, double *s, dcomplex *b, int *ldb, 
                        dcomplex *x, int *ldx, double *rcond, double *berr, 
                        int *n_err_bnds__, double *err_bnds_norm__, double *
                        err_bnds_comp__, int *nparams, double *params, dcomplex *
                        work, double *rwork, int *info);

extern int C2F(zhesv)(char *uplo, int *n, int *nrhs, 
                      dcomplex *a, int *lda, int *ipiv, dcomplex *b, 
                      int *ldb, dcomplex *work, int *lwork, int *info);

extern int C2F(zhesvx)(char *fact, char *uplo, int *n, int *
                       nrhs, dcomplex *a, int *lda, dcomplex *af, int *
                       ldaf, int *ipiv, dcomplex *b, int *ldb, dcomplex *x, 
                       int *ldx, double *rcond, double *ferr, double *berr, 
                       dcomplex *work, int *lwork, double *rwork, int *info);

extern int C2F(zhesvxx)(char *fact, char *uplo, int *n, int *
                        nrhs, dcomplex *a, int *lda, dcomplex *af, int *
                        ldaf, int *ipiv, char *equed, double *s, dcomplex *b, 
                        int *ldb, dcomplex *x, int *ldx, double *rcond, 
                        double *rpvgrw, double *berr, int *n_err_bnds__, 
                        double *err_bnds_norm__, double *err_bnds_comp__, int *
                        nparams, double *params, dcomplex *work, double *rwork, 
                        int *info);

extern int C2F(zhetd2)(char *uplo, int *n, dcomplex *a, 
                       int *lda, double *d__, double *e, dcomplex *tau, 
                       int *info);

extern int C2F(zhetf2)(char *uplo, int *n, dcomplex *a, 
                       int *lda, int *ipiv, int *info);

extern int C2F(zhetrd)(char *uplo, int *n, dcomplex *a, 
                       int *lda, double *d__, double *e, dcomplex *tau, 
                       dcomplex *work, int *lwork, int *info);

extern int C2F(zhetrf)(char *uplo, int *n, dcomplex *a, 
                       int *lda, int *ipiv, dcomplex *work, int *lwork, 
                       int *info);

extern int C2F(zhetri)(char *uplo, int *n, dcomplex *a, 
                       int *lda, int *ipiv, dcomplex *work, int *info);

extern int C2F(zhetrs)(char *uplo, int *n, int *nrhs, 
                       dcomplex *a, int *lda, int *ipiv, dcomplex *b, 
                       int *ldb, int *info);

extern int C2F(zhfrk)(char *transr, char *uplo, char *trans, int *n, 
                      int *k, double *alpha, dcomplex *a, int *lda, 
                      double *beta, dcomplex *c__);

extern int C2F(zhgeqz)(char *job, char *compq, char *compz, int *n, 
                       int *ilo, int *ihi, dcomplex *h__, int *ldh, 
                       dcomplex *t, int *ldt, dcomplex *alpha, dcomplex *
                       beta, dcomplex *q, int *ldq, dcomplex *z__, int *
                       ldz, dcomplex *work, int *lwork, double *rwork, int *
                       info);

extern int C2F(zhpcon)(char *uplo, int *n, dcomplex *ap, 
                       int *ipiv, double *anorm, double *rcond, dcomplex *
                       work, int *info);

extern int C2F(zhpev)(char *jobz, char *uplo, int *n, dcomplex 
                      *ap, double *w, dcomplex *z__, int *ldz, dcomplex *
                      work, double *rwork, int *info);

extern int C2F(zhpevd)(char *jobz, char *uplo, int *n, 
                       dcomplex *ap, double *w, dcomplex *z__, int *ldz, 
                       dcomplex *work, int *lwork, double *rwork, int *
                       lrwork, int *iwork, int *liwork, int *info);

extern int C2F(zhpevx)(char *jobz, char *range, char *uplo, int *n, 
                       dcomplex *ap, double *vl, double *vu, int *il, 
                       int *iu, double *abstol, int *m, double *w, 
                       dcomplex *z__, int *ldz, dcomplex *work, double *
                       rwork, int *iwork, int *ifail, int *info);

extern int C2F(zhpgst)(int *itype, char *uplo, int *n, 
                       dcomplex *ap, dcomplex *bp, int *info);

extern int C2F(zhpgv)(int *itype, char *jobz, char *uplo, int *
                      n, dcomplex *ap, dcomplex *bp, double *w, dcomplex 
                      *z__, int *ldz, dcomplex *work, double *rwork, int *
                      info);

extern int C2F(zhpgvd)(int *itype, char *jobz, char *uplo, int *
                       n, dcomplex *ap, dcomplex *bp, double *w, dcomplex 
                       *z__, int *ldz, dcomplex *work, int *lwork, double *
                       rwork, int *lrwork, int *iwork, int *liwork, int *
                       info);

extern int C2F(zhpgvx)(int *itype, char *jobz, char *range, char *
                       uplo, int *n, dcomplex *ap, dcomplex *bp, double *
                       vl, double *vu, int *il, int *iu, double *abstol, 
                       int *m, double *w, dcomplex *z__, int *ldz, 
                       dcomplex *work, double *rwork, int *iwork, int *
                       ifail, int *info);

extern int C2F(zhprfs)(char *uplo, int *n, int *nrhs, 
                       dcomplex *ap, dcomplex *afp, int *ipiv, dcomplex *
                       b, int *ldb, dcomplex *x, int *ldx, double *ferr, 
                       double *berr, dcomplex *work, double *rwork, int *
                       info);

extern int C2F(zhpsv)(char *uplo, int *n, int *nrhs, 
                      dcomplex *ap, int *ipiv, dcomplex *b, int *ldb, 
                      int *info);

extern int C2F(zhpsvx)(char *fact, char *uplo, int *n, int *
                       nrhs, dcomplex *ap, dcomplex *afp, int *ipiv, 
                       dcomplex *b, int *ldb, dcomplex *x, int *ldx, 
                       double *rcond, double *ferr, double *berr, dcomplex *
                       work, double *rwork, int *info);

extern int C2F(zhptrd)(char *uplo, int *n, dcomplex *ap, 
                       double *d__, double *e, dcomplex *tau, int *info);

extern int C2F(zhptrf)(char *uplo, int *n, dcomplex *ap, 
                       int *ipiv, int *info);

extern int C2F(zhptri)(char *uplo, int *n, dcomplex *ap, 
                       int *ipiv, dcomplex *work, int *info);

extern int C2F(zhptrs)(char *uplo, int *n, int *nrhs, 
                       dcomplex *ap, int *ipiv, dcomplex *b, int *ldb, 
                       int *info);

extern int C2F(zhsein)(char *side, char *eigsrc, char *initv, int *
                       select, int *n, dcomplex *h__, int *ldh, dcomplex *
                       w, dcomplex *vl, int *ldvl, dcomplex *vr, int *ldvr, 
                       int *mm, int *m, dcomplex *work, double *rwork, 
                       int *ifaill, int *ifailr, int *info);

extern int C2F(zhseqr)(char *job, char *compz, int *n, int *ilo, 
                       int *ihi, dcomplex *h__, int *ldh, dcomplex *w, 
                       dcomplex *z__, int *ldz, dcomplex *work, int *lwork, 
                       int *info);

extern int C2F(zla_gbamv_)(int *trans, int *m, int *n, 
                           int *kl, int *ku, double *alpha, dcomplex *ab, 
                           int *ldab, dcomplex *x, int *incx, double *beta, 
                           double *y, int *incy);

double C2F(zla_gbrcond_c_)(char *trans, int *n, int *kl, int *ku, 
                           dcomplex *ab, int *ldab, dcomplex *afb, int *ldafb, 
                           int *ipiv, double *c__, int *capply, int *info, 
                           dcomplex *work, double *rwork, int trans_len);

double C2F(zla_gbrcond_x_)(char *trans, int *n, int *kl, int *ku, 
                           dcomplex *ab, int *ldab, dcomplex *afb, int *ldafb, 
                           int *ipiv, dcomplex *x, int *info, dcomplex *work, 
                           double *rwork, int trans_len);

extern int C2F(zla_gbrfsx_extended_)(int *prec_type__, int *
                                     trans_type__, int *n, int *kl, int *ku, int *nrhs, 
                                     dcomplex *ab, int *ldab, dcomplex *afb, int *ldafb, 
                                     int *ipiv, int *colequ, double *c__, dcomplex *b, 
                                     int *ldb, dcomplex *y, int *ldy, double *berr_out__, 
                                     int *n_norms__, double *errs_n__, double *errs_c__, 
                                     dcomplex *res, double *ayb, dcomplex *dy, dcomplex 
                                     *y_tail__, double *rcond, int *ithresh, double *rthresh, 
                                     double *dz_ub__, int *ignore_cwise__, int *info);

double C2F(zla_gbrpvgrw_)(int *n, int *kl, int *ku, int *
                          ncols, dcomplex *ab, int *ldab, dcomplex *afb, int *
                          ldafb);

extern int C2F(zla_geamv_)(int *trans, int *m, int *n, 
                           double *alpha, dcomplex *a, int *lda, dcomplex *x, 
                           int *incx, double *beta, double *y, int *incy);

double C2F(zla_gercond_c_)(char *trans, int *n, dcomplex *a, int 
                           *lda, dcomplex *af, int *ldaf, int *ipiv, double *
                           c__, int *capply, int *info, dcomplex *work, double *
                           rwork, int trans_len);

double C2F(zla_gercond_x_)(char *trans, int *n, dcomplex *a, int 
                           *lda, dcomplex *af, int *ldaf, int *ipiv, dcomplex *
                           x, int *info, dcomplex *work, double *rwork, int 
                           trans_len);

extern int C2F(zla_gerfsx_extended_)(int *prec_type__, int *
                                     trans_type__, int *n, int *nrhs, dcomplex *a, int *
                                     lda, dcomplex *af, int *ldaf, int *ipiv, int *colequ,
                                     double *c__, dcomplex *b, int *ldb, dcomplex *y, 
                                     int *ldy, double *berr_out__, int *n_norms__, double *
                                     errs_n__, double *errs_c__, dcomplex *res, double *ayb, 
                                     dcomplex *dy, dcomplex *y_tail__, double *rcond, 
                                     int *ithresh, double *rthresh, double *dz_ub__, int *
                                     ignore_cwise__, int *info);

extern int C2F(zla_heamv_)(int *uplo, int *n, double *alpha,
                           dcomplex *a, int *lda, dcomplex *x, int *incx, 
                           double *beta, double *y, int *incy);

double C2F(zla_hercond_c_)(char *uplo, int *n, dcomplex *a, int *
                           lda, dcomplex *af, int *ldaf, int *ipiv, double *c__,
                           int *capply, int *info, dcomplex *work, double *
                           rwork, int uplo_len);

double C2F(zla_hercond_x_)(char *uplo, int *n, dcomplex *a, int *
                           lda, dcomplex *af, int *ldaf, int *ipiv, dcomplex *
                           x, int *info, dcomplex *work, double *rwork, int 
                           uplo_len);

extern int C2F(zla_herfsx_extended_)(int *prec_type__, char *uplo, 
                                     int *n, int *nrhs, dcomplex *a, int *lda, 
                                     dcomplex *af, int *ldaf, int *ipiv, int *colequ, 
                                     double *c__, dcomplex *b, int *ldb, dcomplex *y, 
                                     int *ldy, double *berr_out__, int *n_norms__, double *
                                     errs_n__, double *errs_c__, dcomplex *res, double *ayb, 
                                     dcomplex *dy, dcomplex *y_tail__, double *rcond, 
                                     int *ithresh, double *rthresh, double *dz_ub__, int *
                                     ignore_cwise__, int *info, int uplo_len);

double C2F(zla_herpvgrw_)(char *uplo, int *n, int *info, 
                          dcomplex *a, int *lda, dcomplex *af, int *ldaf, 
                          int *ipiv, double *work, int uplo_len);

extern int C2F(zla_lin_berr_)(int *n, int *nz, int *nrhs, 
                              dcomplex *res, double *ayb, double *berr);

double C2F(zla_porcond_c_)(char *uplo, int *n, dcomplex *a, int *
                           lda, dcomplex *af, int *ldaf, double *c__, int *
                           capply, int *info, dcomplex *work, double *rwork, int 
                           uplo_len);

double C2F(zla_porcond_x_)(char *uplo, int *n, dcomplex *a, int *
                           lda, dcomplex *af, int *ldaf, dcomplex *x, int *
                           info, dcomplex *work, double *rwork, int uplo_len);

extern int C2F(zla_porfsx_extended_)(int *prec_type__, char *uplo, 
                                     int *n, int *nrhs, dcomplex *a, int *lda, 
                                     dcomplex *af, int *ldaf, int *colequ, double *c__, 
                                     dcomplex *b, int *ldb, dcomplex *y, int *ldy, 
                                     double *berr_out__, int *n_norms__, double *errs_n__, 
                                     double *errs_c__, dcomplex *res, double *ayb, 
                                     dcomplex *dy, dcomplex *y_tail__, double *rcond, 
                                     int *ithresh, double *rthresh, double *dz_ub__, int *
                                     ignore_cwise__, int *info, int uplo_len);

double C2F(zla_porpvgrw_)(char *uplo, int *ncols, dcomplex *a, 
                          int *lda, dcomplex *af, int *ldaf, double *work, 
                          int uplo_len);

double C2F(zla_rpvgrw_)(int *n, int *ncols, dcomplex *a, int 
                        *lda, dcomplex *af, int *ldaf);

extern int C2F(zla_syamv_)(int *uplo, int *n, double *alpha,
                           dcomplex *a, int *lda, dcomplex *x, int *incx, 
                           double *beta, double *y, int *incy);

double C2F(zla_syrcond_c_)(char *uplo, int *n, dcomplex *a, int *
                           lda, dcomplex *af, int *ldaf, int *ipiv, double *c__,
                           int *capply, int *info, dcomplex *work, double *
                           rwork, int uplo_len);

double C2F(zla_syrcond_x_)(char *uplo, int *n, dcomplex *a, int *
                           lda, dcomplex *af, int *ldaf, int *ipiv, dcomplex *
                           x, int *info, dcomplex *work, double *rwork, int 
                           uplo_len);

extern int C2F(zla_syrfsx_extended_)(int *prec_type__, char *uplo, 
                                     int *n, int *nrhs, dcomplex *a, int *lda, 
                                     dcomplex *af, int *ldaf, int *ipiv, int *colequ, 
                                     double *c__, dcomplex *b, int *ldb, dcomplex *y, 
                                     int *ldy, double *berr_out__, int *n_norms__, double *
                                     errs_n__, double *errs_c__, dcomplex *res, double *ayb, 
                                     dcomplex *dy, dcomplex *y_tail__, double *rcond, 
                                     int *ithresh, double *rthresh, double *dz_ub__, int *
                                     ignore_cwise__, int *info, int uplo_len);

double C2F(zla_syrpvgrw_)(char *uplo, int *n, int *info, 
                          dcomplex *a, int *lda, dcomplex *af, int *ldaf, 
                          int *ipiv, double *work, int uplo_len);

extern int C2F(zla_wwaddw_)(int *n, dcomplex *x, dcomplex 
                            *y, dcomplex *w);

extern int C2F(zlabrd)(int *m, int *n, int *nb, 
                       dcomplex *a, int *lda, double *d__, double *e, 
                       dcomplex *tauq, dcomplex *taup, dcomplex *x, int *
                       ldx, dcomplex *y, int *ldy);

extern int C2F(zlacgv)(int *n, dcomplex *x, int *incx);

extern int C2F(zlacn2)(int *n, dcomplex *v, dcomplex *x, 
                       double *est, int *kase, int *isave);

extern int C2F(zlacon)(int *n, dcomplex *v, dcomplex *x, 
                       double *est, int *kase);

extern int C2F(zlacp2)(char *uplo, int *m, int *n, double *
                       a, int *lda, dcomplex *b, int *ldb);

extern int C2F(zlacpy)(char *uplo, int *m, int *n, 
                       dcomplex *a, int *lda, dcomplex *b, int *ldb);

extern int C2F(zlacrm)(int *m, int *n, dcomplex *a, 
                       int *lda, double *b, int *ldb, dcomplex *c__, 
                       int *ldc, double *rwork);

extern int C2F(zlacrt)(int *n, dcomplex *cx, int *incx, 
                       dcomplex *cy, int *incy, dcomplex *c__, dcomplex *
                       s);

/* Double dcomplex */ void C2F(zladiv)(dcomplex * ret_val, dcomplex *x, 
                                       dcomplex *y);

extern int C2F(zlaed0)(int *qsiz, int *n, double *d__, 
                       double *e, dcomplex *q, int *ldq, dcomplex *qstore, 
                       int *ldqs, double *rwork, int *iwork, int *info);

extern int C2F(zlaed7)(int *n, int *cutpnt, int *qsiz, 
                       int *tlvls, int *curlvl, int *curpbm, double *d__, 
                       dcomplex *q, int *ldq, double *rho, int *indxq, 
                       double *qstore, int *qptr, int *prmptr, int *perm, 
                       int *givptr, int *givcol, double *givnum, dcomplex *
                       work, double *rwork, int *iwork, int *info);

extern int C2F(zlaed8)(int *k, int *n, int *qsiz, 
                       dcomplex *q, int *ldq, double *d__, double *rho, 
                       int *cutpnt, double *z__, double *dlamda, dcomplex *
                       q2, int *ldq2, double *w, int *indxp, int *indx, 
                       int *indxq, int *perm, int *givptr, int *givcol, 
                       double *givnum, int *info);

extern int C2F(zlaein)(int *rightv, int *noinit, int *n, 
                       dcomplex *h__, int *ldh, dcomplex *w, dcomplex *v, 
                       dcomplex *b, int *ldb, double *rwork, double *eps3, 
                       double *smlnum, int *info);

extern int C2F(zlaesy)(dcomplex *a, dcomplex *b, 
                       dcomplex *c__, dcomplex *rt1, dcomplex *rt2, 
                       dcomplex *evscal, dcomplex *cs1, dcomplex *sn1);

extern int C2F(zlaev2)(dcomplex *a, dcomplex *b, 
                       dcomplex *c__, double *rt1, double *rt2, double *cs1, 
                       dcomplex *sn1);

extern int C2F(zlag2c)(int *m, int *n, dcomplex *a, 
                       int *lda, dcomplex *sa, int *ldsa, int *info);

extern int C2F(zlags2)(int *upper, double *a1, dcomplex *
                       a2, double *a3, double *b1, dcomplex *b2, double *b3, 
                       double *csu, dcomplex *snu, double *csv, dcomplex *
                       snv, double *csq, dcomplex *snq);

extern int C2F(zlagtm)(char *trans, int *n, int *nrhs, 
                       double *alpha, dcomplex *dl, dcomplex *d__, 
                       dcomplex *du, dcomplex *x, int *ldx, double *beta, 
                       dcomplex *b, int *ldb);

extern int C2F(zlahef)(char *uplo, int *n, int *nb, int *kb, 
                       dcomplex *a, int *lda, int *ipiv, dcomplex *w, 
                       int *ldw, int *info);

extern int C2F(zlahqr)(int *wantt, int *wantz, int *n, 
                       int *ilo, int *ihi, dcomplex *h__, int *ldh, 
                       dcomplex *w, int *iloz, int *ihiz, dcomplex *z__, 
                       int *ldz, int *info);

extern int C2F(zlahr2)(int *n, int *k, int *nb, 
                       dcomplex *a, int *lda, dcomplex *tau, dcomplex *t, 
                       int *ldt, dcomplex *y, int *ldy);

extern int C2F(zlahrd)(int *n, int *k, int *nb, 
                       dcomplex *a, int *lda, dcomplex *tau, dcomplex *t, 
                       int *ldt, dcomplex *y, int *ldy);

extern int C2F(zlaic1)(int *job, int *j, dcomplex *x, 
                       double *sest, dcomplex *w, dcomplex *gamma, double *
                       sestpr, dcomplex *s, dcomplex *c__);

extern int C2F(zlals0)(int *icompq, int *nl, int *nr, 
                       int *sqre, int *nrhs, dcomplex *b, int *ldb, 
                       dcomplex *bx, int *ldbx, int *perm, int *givptr, 
                       int *givcol, int *ldgcol, double *givnum, int *ldgnum, 
                       double *poles, double *difl, double *difr, double *
                       z__, int *k, double *c__, double *s, double *rwork, 
                       int *info);

extern int C2F(zlalsa)(int *icompq, int *smlsiz, int *n, 
                       int *nrhs, dcomplex *b, int *ldb, dcomplex *bx, 
                       int *ldbx, double *u, int *ldu, double *vt, int *
                       k, double *difl, double *difr, double *z__, double *
                       poles, int *givptr, int *givcol, int *ldgcol, int *
                       perm, double *givnum, double *c__, double *s, double *
                       rwork, int *iwork, int *info);

extern int C2F(zlalsd)(char *uplo, int *smlsiz, int *n, int 
                       *nrhs, double *d__, double *e, dcomplex *b, int *ldb, 
                       double *rcond, int *rank, dcomplex *work, double *
                       rwork, int *iwork, int *info);

double C2F(zlangb)(char *norm, int *n, int *kl, int *ku, 
                   dcomplex *ab, int *ldab, double *work);

double C2F(zlange)(char *norm, int *m, int *n, dcomplex *a, 
                   int *lda, double *work);

double C2F(zlangt)(char *norm, int *n, dcomplex *dl, dcomplex *
                   d__, dcomplex *du);

double C2F(zlanhb)(char *norm, char *uplo, int *n, int *k, 
                   dcomplex *ab, int *ldab, double *work);

double C2F(zlanhe)(char *norm, char *uplo, int *n, dcomplex *a, 
                   int *lda, double *work);

double C2F(zlanhf)(char *norm, char *transr, char *uplo, int *n, 
                   dcomplex *a, double *work);

double C2F(zlanhp)(char *norm, char *uplo, int *n, dcomplex *ap, 
                   double *work);

double C2F(zlanhs)(char *norm, int *n, dcomplex *a, int *lda, 
                   double *work);

double C2F(zlanht)(char *norm, int *n, double *d__, dcomplex *e);

double C2F(zlansb)(char *norm, char *uplo, int *n, int *k, 
                   dcomplex *ab, int *ldab, double *work);

double C2F(zlansp)(char *norm, char *uplo, int *n, dcomplex *ap, 
                   double *work);

double C2F(zlansy)(char *norm, char *uplo, int *n, dcomplex *a, 
                   int *lda, double *work);

double C2F(zlantb)(char *norm, char *uplo, char *diag, int *n, int *k, 
                   dcomplex *ab, int *ldab, double *work);

double C2F(zlantp)(char *norm, char *uplo, char *diag, int *n, 
                   dcomplex *ap, double *work);

double C2F(zlantr)(char *norm, char *uplo, char *diag, int *m, int *n, 
                   dcomplex *a, int *lda, double *work);

extern int C2F(zlapll)(int *n, dcomplex *x, int *incx, 
                       dcomplex *y, int *incy, double *ssmin);

extern int C2F(zlapmt)(int *forwrd, int *m, int *n, 
                       dcomplex *x, int *ldx, int *k);

extern int C2F(zlaqgb)(int *m, int *n, int *kl, int *ku, 
                       dcomplex *ab, int *ldab, double *r__, double *c__, 
                       double *rowcnd, double *colcnd, double *amax, char *equed);

extern int C2F(zlaqge)(int *m, int *n, dcomplex *a, 
                       int *lda, double *r__, double *c__, double *rowcnd, 
                       double *colcnd, double *amax, char *equed);

extern int C2F(zlaqhb)(char *uplo, int *n, int *kd, 
                       dcomplex *ab, int *ldab, double *s, double *scond, 
                       double *amax, char *equed);

extern int C2F(zlaqhe)(char *uplo, int *n, dcomplex *a, 
                       int *lda, double *s, double *scond, double *amax, 
                       char *equed);

extern int C2F(zlaqhp)(char *uplo, int *n, dcomplex *ap, 
                       double *s, double *scond, double *amax, char *equed);

extern int C2F(zlaqp2)(int *m, int *n, int *offset, 
                       dcomplex *a, int *lda, int *jpvt, dcomplex *tau, 
                       double *vn1, double *vn2, dcomplex *work);

extern int C2F(zlaqps)(int *m, int *n, int *offset, int 
                       *nb, int *kb, dcomplex *a, int *lda, int *jpvt, 
                       dcomplex *tau, double *vn1, double *vn2, dcomplex *
                       auxv, dcomplex *f, int *ldf);

extern int C2F(zlaqr0)(int *wantt, int *wantz, int *n, 
                       int *ilo, int *ihi, dcomplex *h__, int *ldh, 
                       dcomplex *w, int *iloz, int *ihiz, dcomplex *z__, 
                       int *ldz, dcomplex *work, int *lwork, int *info);

extern int C2F(zlaqr1)(int *n, dcomplex *h__, int *ldh, 
                       dcomplex *s1, dcomplex *s2, dcomplex *v);

extern int C2F(zlaqr2)(int *wantt, int *wantz, int *n, 
                       int *ktop, int *kbot, int *nw, dcomplex *h__, 
                       int *ldh, int *iloz, int *ihiz, dcomplex *z__, 
                       int *ldz, int *ns, int *nd, dcomplex *sh, 
                       dcomplex *v, int *ldv, int *nh, dcomplex *t, 
                       int *ldt, int *nv, dcomplex *wv, int *ldwv, 
                       dcomplex *work, int *lwork);

extern int C2F(zlaqr3)(int *wantt, int *wantz, int *n, 
                       int *ktop, int *kbot, int *nw, dcomplex *h__, 
                       int *ldh, int *iloz, int *ihiz, dcomplex *z__, 
                       int *ldz, int *ns, int *nd, dcomplex *sh, 
                       dcomplex *v, int *ldv, int *nh, dcomplex *t, 
                       int *ldt, int *nv, dcomplex *wv, int *ldwv, 
                       dcomplex *work, int *lwork);

extern int C2F(zlaqr4)(int *wantt, int *wantz, int *n, 
                       int *ilo, int *ihi, dcomplex *h__, int *ldh, 
                       dcomplex *w, int *iloz, int *ihiz, dcomplex *z__, 
                       int *ldz, dcomplex *work, int *lwork, int *info);

extern int C2F(zlaqr5)(int *wantt, int *wantz, int *kacc22, 
                       int *n, int *ktop, int *kbot, int *nshfts, 
                       dcomplex *s, dcomplex *h__, int *ldh, int *iloz, 
                       int *ihiz, dcomplex *z__, int *ldz, dcomplex *v, 
                       int *ldv, dcomplex *u, int *ldu, int *nv, 
                       dcomplex *wv, int *ldwv, int *nh, dcomplex *wh, 
                       int *ldwh);

extern int C2F(zlaqsb)(char *uplo, int *n, int *kd, 
                       dcomplex *ab, int *ldab, double *s, double *scond, 
                       double *amax, char *equed);

extern int C2F(zlaqsp)(char *uplo, int *n, dcomplex *ap, 
                       double *s, double *scond, double *amax, char *equed);

extern int C2F(zlaqsy)(char *uplo, int *n, dcomplex *a, 
                       int *lda, double *s, double *scond, double *amax, 
                       char *equed);

extern int C2F(zlar1v)(int *n, int *b1, int *bn, double 
                       *lambda, double *d__, double *l, double *ld, double *
                       lld, double *pivmin, double *gaptol, dcomplex *z__, 
                       int *wantnc, int *negcnt, double *ztz, double *mingma, 
                       int *r__, int *isuppz, double *nrminv, double *resid, 
                       double *rqcorr, double *work);

extern int C2F(zlar2v)(int *n, dcomplex *x, dcomplex *y, 
                       dcomplex *z__, int *incx, double *c__, dcomplex *s, 
                       int *incc);

extern int C2F(zlarcm)(int *m, int *n, double *a, int *
                       lda, dcomplex *b, int *ldb, dcomplex *c__, int *ldc, 
                       double *rwork);

extern int C2F(zlarf)(char *side, int *m, int *n, dcomplex 
                      *v, int *incv, dcomplex *tau, dcomplex *c__, int *
                      ldc, dcomplex *work);

extern int C2F(zlarfb)(char *side, char *trans, char *direct, char *
                       storev, int *m, int *n, int *k, dcomplex *v, int 
                       *ldv, dcomplex *t, int *ldt, dcomplex *c__, int *
                       ldc, dcomplex *work, int *ldwork);

extern int C2F(zlarfg)(int *n, dcomplex *alpha, dcomplex *
                       x, int *incx, dcomplex *tau);

extern int C2F(zlarfp)(int *n, dcomplex *alpha, dcomplex *
                       x, int *incx, dcomplex *tau);

extern int C2F(zlarft)(char *direct, char *storev, int *n, int *
                       k, dcomplex *v, int *ldv, dcomplex *tau, dcomplex *
                       t, int *ldt);

extern int C2F(zlarfx)(char *side, int *m, int *n, 
                       dcomplex *v, dcomplex *tau, dcomplex *c__, int *
                       ldc, dcomplex *work);

extern int C2F(zlargv)(int *n, dcomplex *x, int *incx, 
                       dcomplex *y, int *incy, double *c__, int *incc);

extern int C2F(zlarnv)(int *idist, int *iseed, int *n, 
                       dcomplex *x);

extern int C2F(zlarrv)(int *n, double *vl, double *vu, 
                       double *d__, double *l, double *pivmin, int *isplit, 
                       int *m, int *dol, int *dou, double *minrgp, 
                       double *rtol1, double *rtol2, double *w, double *werr, 
                       double *wgap, int *iblock, int *indexw, double *gers, 
                       dcomplex *z__, int *ldz, int *isuppz, double *work, 
                       int *iwork, int *info);

extern int C2F(zlarscl2)(int *m, int *n, double *d__, 
                         dcomplex *x, int *ldx);

extern int C2F(zlartg)(dcomplex *f, dcomplex *g, double *
                       cs, dcomplex *sn, dcomplex *r__);

extern int C2F(zlartv)(int *n, dcomplex *x, int *incx, 
                       dcomplex *y, int *incy, double *c__, dcomplex *s, 
                       int *incc);

extern int C2F(zlarz)(char *side, int *m, int *n, int *l, 
                      dcomplex *v, int *incv, dcomplex *tau, dcomplex *
                      c__, int *ldc, dcomplex *work);

extern int C2F(zlarzb)(char *side, char *trans, char *direct, char *
                       storev, int *m, int *n, int *k, int *l, dcomplex 
                       *v, int *ldv, dcomplex *t, int *ldt, dcomplex *c__, 
                       int *ldc, dcomplex *work, int *ldwork);

extern int C2F(zlarzt)(char *direct, char *storev, int *n, int *
                       k, dcomplex *v, int *ldv, dcomplex *tau, dcomplex *
                       t, int *ldt);

extern int C2F(zlascl)(char *type__, int *kl, int *ku, 
                       double *cfrom, double *cto, int *m, int *n, 
                       dcomplex *a, int *lda, int *info);

extern int C2F(zlascl2)(int *m, int *n, double *d__, 
                        dcomplex *x, int *ldx);

extern int C2F(zlaset)(char *uplo, int *m, int *n, 
                       dcomplex *alpha, dcomplex *beta, dcomplex *a, int *
                       lda);

extern int C2F(zlasr)(char *side, char *pivot, char *direct, int *m, 
                      int *n, double *c__, double *s, dcomplex *a, 
                      int *lda);

extern int C2F(zlassq)(int *n, dcomplex *x, int *incx, 
                       double *scale, double *sumsq);

extern int C2F(zlaswp)(int *n, dcomplex *a, int *lda, 
                       int *k1, int *k2, int *ipiv, int *incx);

extern int C2F(zlasyf)(char *uplo, int *n, int *nb, int *kb, 
                       dcomplex *a, int *lda, int *ipiv, dcomplex *w, 
                       int *ldw, int *info);

extern int C2F(zlat2c)(char *uplo, int *n, dcomplex *a, 
                       int *lda, dcomplex *sa, int *ldsa, int *info);

extern int C2F(zlatbs)(char *uplo, char *trans, char *diag, char *
                       normin, int *n, int *kd, dcomplex *ab, int *ldab, 
                       dcomplex *x, double *scale, double *cnorm, int *info);

extern int C2F(zlatdf)(int *ijob, int *n, dcomplex *z__, 
                       int *ldz, dcomplex *rhs, double *rdsum, double *
                       rdscal, int *ipiv, int *jpiv);

extern int C2F(zlatps)(char *uplo, char *trans, char *diag, char *
                       normin, int *n, dcomplex *ap, dcomplex *x, double *
                       scale, double *cnorm, int *info);

extern int C2F(zlatrd)(char *uplo, int *n, int *nb, 
                       dcomplex *a, int *lda, double *e, dcomplex *tau, 
                       dcomplex *w, int *ldw);

extern int C2F(zlatrs)(char *uplo, char *trans, char *diag, char *
                       normin, int *n, dcomplex *a, int *lda, dcomplex *x, 
                       double *scale, double *cnorm, int *info);

extern int C2F(zlatrz)(int *m, int *n, int *l, 
                       dcomplex *a, int *lda, dcomplex *tau, dcomplex *
                       work);

extern int C2F(zlatzm)(char *side, int *m, int *n, 
                       dcomplex *v, int *incv, dcomplex *tau, dcomplex *
                       c1, dcomplex *c2, int *ldc, dcomplex *work);

extern int C2F(zlauu2)(char *uplo, int *n, dcomplex *a, 
                       int *lda, int *info);

extern int C2F(zlauum)(char *uplo, int *n, dcomplex *a, 
                       int *lda, int *info);

extern int C2F(zpbcon)(char *uplo, int *n, int *kd, 
                       dcomplex *ab, int *ldab, double *anorm, double *
                       rcond, dcomplex *work, double *rwork, int *info);

extern int C2F(zpbequ)(char *uplo, int *n, int *kd, 
                       dcomplex *ab, int *ldab, double *s, double *scond, 
                       double *amax, int *info);

extern int C2F(zpbrfs)(char *uplo, int *n, int *kd, int *
                       nrhs, dcomplex *ab, int *ldab, dcomplex *afb, int *
                       ldafb, dcomplex *b, int *ldb, dcomplex *x, int *ldx, 
                       double *ferr, double *berr, dcomplex *work, double *
                       rwork, int *info);

extern int C2F(zpbstf)(char *uplo, int *n, int *kd, 
                       dcomplex *ab, int *ldab, int *info);

extern int C2F(zpbsv)(char *uplo, int *n, int *kd, int *
                      nrhs, dcomplex *ab, int *ldab, dcomplex *b, int *
                      ldb, int *info);

extern int C2F(zpbsvx)(char *fact, char *uplo, int *n, int *kd, 
                       int *nrhs, dcomplex *ab, int *ldab, dcomplex *afb, 
                       int *ldafb, char *equed, double *s, dcomplex *b, int 
                       *ldb, dcomplex *x, int *ldx, double *rcond, double *
                       ferr, double *berr, dcomplex *work, double *rwork, 
                       int *info);

extern int C2F(zpbtf2)(char *uplo, int *n, int *kd, 
                       dcomplex *ab, int *ldab, int *info);

extern int C2F(zpbtrf)(char *uplo, int *n, int *kd, 
                       dcomplex *ab, int *ldab, int *info);

extern int C2F(zpbtrs)(char *uplo, int *n, int *kd, int *
                       nrhs, dcomplex *ab, int *ldab, dcomplex *b, int *
                       ldb, int *info);

extern int C2F(zpftrf)(char *transr, char *uplo, int *n, 
                       dcomplex *a, int *info);

extern int C2F(zpftri)(char *transr, char *uplo, int *n, 
                       dcomplex *a, int *info);

extern int C2F(zpftrs)(char *transr, char *uplo, int *n, int *
                       nrhs, dcomplex *a, dcomplex *b, int *ldb, int *info);

extern int C2F(zpocon)(char *uplo, int *n, dcomplex *a, 
                       int *lda, double *anorm, double *rcond, dcomplex *
                       work, double *rwork, int *info);

extern int C2F(zpoequ)(int *n, dcomplex *a, int *lda, 
                       double *s, double *scond, double *amax, int *info);

extern int C2F(zpoequb)(int *n, dcomplex *a, int *lda, 
                        double *s, double *scond, double *amax, int *info);

extern int C2F(zporfs)(char *uplo, int *n, int *nrhs, 
                       dcomplex *a, int *lda, dcomplex *af, int *ldaf, 
                       dcomplex *b, int *ldb, dcomplex *x, int *ldx, 
                       double *ferr, double *berr, dcomplex *work, double *
                       rwork, int *info);

extern int C2F(zporfsx)(char *uplo, char *equed, int *n, int *
                        nrhs, dcomplex *a, int *lda, dcomplex *af, int *
                        ldaf, double *s, dcomplex *b, int *ldb, dcomplex *x, 
                        int *ldx, double *rcond, double *berr, int *
                        n_err_bnds__, double *err_bnds_norm__, double *
                        err_bnds_comp__, int *nparams, double *params, dcomplex *
                        work, double *rwork, int *info);

extern int C2F(zposv)(char *uplo, int *n, int *nrhs, 
                      dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                      int *info);

extern int C2F(zposvx)(char *fact, char *uplo, int *n, int *
                       nrhs, dcomplex *a, int *lda, dcomplex *af, int *
                       ldaf, char *equed, double *s, dcomplex *b, int *ldb, 
                       dcomplex *x, int *ldx, double *rcond, double *ferr, 
                       double *berr, dcomplex *work, double *rwork, int *
                       info);

extern int C2F(zposvxx)(char *fact, char *uplo, int *n, int *
                        nrhs, dcomplex *a, int *lda, dcomplex *af, int *
                        ldaf, char *equed, double *s, dcomplex *b, int *ldb, 
                        dcomplex *x, int *ldx, double *rcond, double *rpvgrw, 
                        double *berr, int *n_err_bnds__, double *err_bnds_norm__, 
                        double *err_bnds_comp__, int *nparams, double *params, 
                        dcomplex *work, double *rwork, int *info);

extern int C2F(zpotf2)(char *uplo, int *n, dcomplex *a, 
                       int *lda, int *info);

extern int C2F(zpotrf)(char *uplo, int *n, dcomplex *a, 
                       int *lda, int *info);

extern int C2F(zpotri)(char *uplo, int *n, dcomplex *a, 
                       int *lda, int *info);

extern int C2F(zpotrs)(char *uplo, int *n, int *nrhs, 
                       dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       int *info);

extern int C2F(zppcon)(char *uplo, int *n, dcomplex *ap, 
                       double *anorm, double *rcond, dcomplex *work, double 
                       *rwork, int *info);

extern int C2F(zppequ)(char *uplo, int *n, dcomplex *ap, 
                       double *s, double *scond, double *amax, int *info);

extern int C2F(zpprfs)(char *uplo, int *n, int *nrhs, 
                       dcomplex *ap, dcomplex *afp, dcomplex *b, int *ldb, 
                       dcomplex *x, int *ldx, double *ferr, double *berr, 
                       dcomplex *work, double *rwork, int *info);

extern int C2F(zppsv)(char *uplo, int *n, int *nrhs, 
                      dcomplex *ap, dcomplex *b, int *ldb, int *info);

extern int C2F(zppsvx)(char *fact, char *uplo, int *n, int *
                       nrhs, dcomplex *ap, dcomplex *afp, char *equed, double *
                       s, dcomplex *b, int *ldb, dcomplex *x, int *ldx, 
                       double *rcond, double *ferr, double *berr, dcomplex *
                       work, double *rwork, int *info);

extern int C2F(zpptrf)(char *uplo, int *n, dcomplex *ap, 
                       int *info);

extern int C2F(zpptri)(char *uplo, int *n, dcomplex *ap, 
                       int *info);

extern int C2F(zpptrs)(char *uplo, int *n, int *nrhs, 
                       dcomplex *ap, dcomplex *b, int *ldb, int *info);

extern int C2F(zpstf2)(char *uplo, int *n, dcomplex *a, 
                       int *lda, int *piv, int *rank, double *tol, 
                       double *work, int *info);

extern int C2F(zpstrf)(char *uplo, int *n, dcomplex *a, 
                       int *lda, int *piv, int *rank, double *tol, 
                       double *work, int *info);

extern int C2F(zptcon)(int *n, double *d__, dcomplex *e, 
                       double *anorm, double *rcond, double *rwork, int *
                       info);

extern int C2F(zpteqr)(char *compz, int *n, double *d__, 
                       double *e, dcomplex *z__, int *ldz, double *work, 
                       int *info);

extern int C2F(zptrfs)(char *uplo, int *n, int *nrhs, 
                       double *d__, dcomplex *e, double *df, dcomplex *ef, 
                       dcomplex *b, int *ldb, dcomplex *x, int *ldx, 
                       double *ferr, double *berr, dcomplex *work, double *
                       rwork, int *info);

extern int C2F(zptsv)(int *n, int *nrhs, double *d__, 
                      dcomplex *e, dcomplex *b, int *ldb, int *info);

extern int C2F(zptsvx)(char *fact, int *n, int *nrhs, 
                       double *d__, dcomplex *e, double *df, dcomplex *ef, 
                       dcomplex *b, int *ldb, dcomplex *x, int *ldx, 
                       double *rcond, double *ferr, double *berr, dcomplex *
                       work, double *rwork, int *info);

extern int C2F(zpttrf)(int *n, double *d__, dcomplex *e, 
                       int *info);

extern int C2F(zpttrs)(char *uplo, int *n, int *nrhs, 
                       double *d__, dcomplex *e, dcomplex *b, int *ldb, 
                       int *info);

extern int C2F(zptts2)(int *iuplo, int *n, int *nrhs, 
                       double *d__, dcomplex *e, dcomplex *b, int *ldb);

extern int C2F(zrot)(int *n, dcomplex *cx, int *incx, 
                     dcomplex *cy, int *incy, double *c__, dcomplex *s);

extern int C2F(zspcon)(char *uplo, int *n, dcomplex *ap, 
                       int *ipiv, double *anorm, double *rcond, dcomplex *
                       work, int *info);

extern int C2F(zspmv)(char *uplo, int *n, dcomplex *alpha, 
                      dcomplex *ap, dcomplex *x, int *incx, dcomplex *
                      beta, dcomplex *y, int *incy);

extern int C2F(zspr)(char *uplo, int *n, dcomplex *alpha, 
                     dcomplex *x, int *incx, dcomplex *ap);

extern int C2F(zsprfs)(char *uplo, int *n, int *nrhs, 
                       dcomplex *ap, dcomplex *afp, int *ipiv, dcomplex *
                       b, int *ldb, dcomplex *x, int *ldx, double *ferr, 
                       double *berr, dcomplex *work, double *rwork, int *
                       info);

extern int C2F(zspsv)(char *uplo, int *n, int *nrhs, 
                      dcomplex *ap, int *ipiv, dcomplex *b, int *ldb, 
                      int *info);

extern int C2F(zspsvx)(char *fact, char *uplo, int *n, int *
                       nrhs, dcomplex *ap, dcomplex *afp, int *ipiv, 
                       dcomplex *b, int *ldb, dcomplex *x, int *ldx, 
                       double *rcond, double *ferr, double *berr, dcomplex *
                       work, double *rwork, int *info);

extern int C2F(zsptrf)(char *uplo, int *n, dcomplex *ap, 
                       int *ipiv, int *info);

extern int C2F(zsptri)(char *uplo, int *n, dcomplex *ap, 
                       int *ipiv, dcomplex *work, int *info);

extern int C2F(zsptrs)(char *uplo, int *n, int *nrhs, 
                       dcomplex *ap, int *ipiv, dcomplex *b, int *ldb, 
                       int *info);

extern int C2F(zstedc)(char *compz, int *n, double *d__, 
                       double *e, dcomplex *z__, int *ldz, dcomplex *work, 
                       int *lwork, double *rwork, int *lrwork, int *iwork, 
                       int *liwork, int *info);

extern int C2F(zstegr)(char *jobz, char *range, int *n, double *
                       d__, double *e, double *vl, double *vu, int *il, 
                       int *iu, double *abstol, int *m, double *w, 
                       dcomplex *z__, int *ldz, int *isuppz, double *work, 
                       int *lwork, int *iwork, int *liwork, int *info);

extern int C2F(zstein)(int *n, double *d__, double *e, 
                       int *m, double *w, int *iblock, int *isplit, 
                       dcomplex *z__, int *ldz, double *work, int *iwork, 
                       int *ifail, int *info);

extern int C2F(zstemr)(char *jobz, char *range, int *n, double *
                       d__, double *e, double *vl, double *vu, int *il, 
                       int *iu, int *m, double *w, dcomplex *z__, int *
                       ldz, int *nzc, int *isuppz, int *tryrac, double *work, 
                       int *lwork, int *iwork, int *liwork, int *info);

extern int C2F(zsteqr)(char *compz, int *n, double *d__, 
                       double *e, dcomplex *z__, int *ldz, double *work, 
                       int *info);

extern int C2F(zsycon)(char *uplo, int *n, dcomplex *a, 
                       int *lda, int *ipiv, double *anorm, double *rcond, 
                       dcomplex *work, int *info);

extern int C2F(zsyequb)(char *uplo, int *n, dcomplex *a, 
                        int *lda, double *s, double *scond, double *amax, 
                        dcomplex *work, int *info);

extern int C2F(zsymv)(char *uplo, int *n, dcomplex *alpha, 
                      dcomplex *a, int *lda, dcomplex *x, int *incx, 
                      dcomplex *beta, dcomplex *y, int *incy);

extern int C2F(zsyr)(char *uplo, int *n, dcomplex *alpha, 
                     dcomplex *x, int *incx, dcomplex *a, int *lda);

extern int C2F(zsyrfs)(char *uplo, int *n, int *nrhs, 
                       dcomplex *a, int *lda, dcomplex *af, int *ldaf, 
                       int *ipiv, dcomplex *b, int *ldb, dcomplex *x, 
                       int *ldx, double *ferr, double *berr, dcomplex *work, 
                       double *rwork, int *info);

extern int C2F(zsyrfsx)(char *uplo, char *equed, int *n, int *
                        nrhs, dcomplex *a, int *lda, dcomplex *af, int *
                        ldaf, int *ipiv, double *s, dcomplex *b, int *ldb, 
                        dcomplex *x, int *ldx, double *rcond, double *berr, 
                        int *n_err_bnds__, double *err_bnds_norm__, double *
                        err_bnds_comp__, int *nparams, double *params, dcomplex *
                        work, double *rwork, int *info);

extern int C2F(zsysv)(char *uplo, int *n, int *nrhs, 
                      dcomplex *a, int *lda, int *ipiv, dcomplex *b, 
                      int *ldb, dcomplex *work, int *lwork, int *info);

extern int C2F(zsysvx)(char *fact, char *uplo, int *n, int *
                       nrhs, dcomplex *a, int *lda, dcomplex *af, int *
                       ldaf, int *ipiv, dcomplex *b, int *ldb, dcomplex *x, 
                       int *ldx, double *rcond, double *ferr, double *berr, 
                       dcomplex *work, int *lwork, double *rwork, int *info);

extern int C2F(zsysvxx)(char *fact, char *uplo, int *n, int *
                        nrhs, dcomplex *a, int *lda, dcomplex *af, int *
                        ldaf, int *ipiv, char *equed, double *s, dcomplex *b, 
                        int *ldb, dcomplex *x, int *ldx, double *rcond, 
                        double *rpvgrw, double *berr, int *n_err_bnds__, 
                        double *err_bnds_norm__, double *err_bnds_comp__, int *
                        nparams, double *params, dcomplex *work, double *rwork, 
                        int *info);

extern int C2F(zsytf2)(char *uplo, int *n, dcomplex *a, 
                       int *lda, int *ipiv, int *info);

extern int C2F(zsytrf)(char *uplo, int *n, dcomplex *a, 
                       int *lda, int *ipiv, dcomplex *work, int *lwork, 
                       int *info);

extern int C2F(zsytri)(char *uplo, int *n, dcomplex *a, 
                       int *lda, int *ipiv, dcomplex *work, int *info);

extern int C2F(zsytrs)(char *uplo, int *n, int *nrhs, 
                       dcomplex *a, int *lda, int *ipiv, dcomplex *b, 
                       int *ldb, int *info);

extern int C2F(ztbcon)(char *norm, char *uplo, char *diag, int *n, 
                       int *kd, dcomplex *ab, int *ldab, double *rcond, 
                       dcomplex *work, double *rwork, int *info);

extern int C2F(ztbrfs)(char *uplo, char *trans, char *diag, int *n, 
                       int *kd, int *nrhs, dcomplex *ab, int *ldab, 
                       dcomplex *b, int *ldb, dcomplex *x, int *ldx, 
                       double *ferr, double *berr, dcomplex *work, double *
                       rwork, int *info);

extern int C2F(ztbtrs)(char *uplo, char *trans, char *diag, int *n, 
                       int *kd, int *nrhs, dcomplex *ab, int *ldab, 
                       dcomplex *b, int *ldb, int *info);

extern int C2F(ztfsm)(char *transr, char *side, char *uplo, char *trans, 
                      char *diag, int *m, int *n, dcomplex *alpha, 
                      dcomplex *a, dcomplex *b, int *ldb);

extern int C2F(ztftri)(char *transr, char *uplo, char *diag, int *n, 
                       dcomplex *a, int *info);

extern int C2F(ztfttp)(char *transr, char *uplo, int *n, 
                       dcomplex *arf, dcomplex *ap, int *info);

extern int C2F(ztfttr)(char *transr, char *uplo, int *n, 
                       dcomplex *arf, dcomplex *a, int *lda, int *info);

extern int C2F(ztgevc)(char *side, char *howmny, int *select, 
                       int *n, dcomplex *s, int *lds, dcomplex *p, int 
                       *ldp, dcomplex *vl, int *ldvl, dcomplex *vr, int *
                       ldvr, int *mm, int *m, dcomplex *work, double *rwork, 
                       int *info);

extern int C2F(ztgex2)(int *wantq, int *wantz, int *n, 
                       dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       dcomplex *q, int *ldq, dcomplex *z__, int *ldz, 
                       int *j1, int *info);

extern int C2F(ztgexc)(int *wantq, int *wantz, int *n, 
                       dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       dcomplex *q, int *ldq, dcomplex *z__, int *ldz, 
                       int *ifst, int *ilst, int *info);

extern int C2F(ztgsen)(int *ijob, int *wantq, int *wantz, 
                       int *select, int *n, dcomplex *a, int *lda, 
                       dcomplex *b, int *ldb, dcomplex *alpha, dcomplex *
                       beta, dcomplex *q, int *ldq, dcomplex *z__, int *
                       ldz, int *m, double *pl, double *pr, double *dif, 
                       dcomplex *work, int *lwork, int *iwork, int *liwork, 
                       int *info);

extern int C2F(ztgsja)(char *jobu, char *jobv, char *jobq, int *m, 
                       int *p, int *n, int *k, int *l, dcomplex *a, 
                       int *lda, dcomplex *b, int *ldb, double *tola, 
                       double *tolb, double *alpha, double *beta, dcomplex *
                       u, int *ldu, dcomplex *v, int *ldv, dcomplex *q, 
                       int *ldq, dcomplex *work, int *ncycle, int *info);

extern int C2F(ztgsna)(char *job, char *howmny, int *select, 
                       int *n, dcomplex *a, int *lda, dcomplex *b, int 
                       *ldb, dcomplex *vl, int *ldvl, dcomplex *vr, int *
                       ldvr, double *s, double *dif, int *mm, int *m, 
                       dcomplex *work, int *lwork, int *iwork, int *info);

extern int C2F(ztgsy2)(char *trans, int *ijob, int *m, int *
                       n, dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       dcomplex *c__, int *ldc, dcomplex *d__, int *ldd, 
                       dcomplex *e, int *lde, dcomplex *f, int *ldf, 
                       double *scale, double *rdsum, double *rdscal, int *
                       info);

extern int C2F(ztgsyl)(char *trans, int *ijob, int *m, int *
                       n, dcomplex *a, int *lda, dcomplex *b, int *ldb, 
                       dcomplex *c__, int *ldc, dcomplex *d__, int *ldd, 
                       dcomplex *e, int *lde, dcomplex *f, int *ldf, 
                       double *scale, double *dif, dcomplex *work, int *
                       lwork, int *iwork, int *info);

extern int C2F(ztpcon)(char *norm, char *uplo, char *diag, int *n, 
                       dcomplex *ap, double *rcond, dcomplex *work, double 
                       *rwork, int *info);

extern int C2F(ztprfs)(char *uplo, char *trans, char *diag, int *n, 
                       int *nrhs, dcomplex *ap, dcomplex *b, int *ldb, 
                       dcomplex *x, int *ldx, double *ferr, double *berr, 
                       dcomplex *work, double *rwork, int *info);

extern int C2F(ztptri)(char *uplo, char *diag, int *n, 
                       dcomplex *ap, int *info);

extern int C2F(ztptrs)(char *uplo, char *trans, char *diag, int *n, 
                       int *nrhs, dcomplex *ap, dcomplex *b, int *ldb, 
                       int *info);

extern int C2F(ztpttf)(char *transr, char *uplo, int *n, 
                       dcomplex *ap, dcomplex *arf, int *info);

extern int C2F(ztpttr)(char *uplo, int *n, dcomplex *ap, 
                       dcomplex *a, int *lda, int *info);

extern int C2F(ztrcon)(char *norm, char *uplo, char *diag, int *n, 
                       dcomplex *a, int *lda, double *rcond, dcomplex *
                       work, double *rwork, int *info);

extern int C2F(ztrevc)(char *side, char *howmny, int *select, 
                       int *n, dcomplex *t, int *ldt, dcomplex *vl, 
                       int *ldvl, dcomplex *vr, int *ldvr, int *mm, int 
                       *m, dcomplex *work, double *rwork, int *info);

extern int C2F(ztrexc)(char *compq, int *n, dcomplex *t, 
                       int *ldt, dcomplex *q, int *ldq, int *ifst, int *
                       ilst, int *info);

extern int C2F(ztrrfs)(char *uplo, char *trans, char *diag, int *n, 
                       int *nrhs, dcomplex *a, int *lda, dcomplex *b, 
                       int *ldb, dcomplex *x, int *ldx, double *ferr, 
                       double *berr, dcomplex *work, double *rwork, int *
                       info);

extern int C2F(ztrsen)(char *job, char *compq, int *select, int 
                       *n, dcomplex *t, int *ldt, dcomplex *q, int *ldq, 
                       dcomplex *w, int *m, double *s, double *sep, 
                       dcomplex *work, int *lwork, int *info);

extern int C2F(ztrsna)(char *job, char *howmny, int *select, 
                       int *n, dcomplex *t, int *ldt, dcomplex *vl, 
                       int *ldvl, dcomplex *vr, int *ldvr, double *s, 
                       double *sep, int *mm, int *m, dcomplex *work, 
                       int *ldwork, double *rwork, int *info);

extern int C2F(ztrsyl)(char *trana, char *tranb, int *isgn, int 
                       *m, int *n, dcomplex *a, int *lda, dcomplex *b, 
                       int *ldb, dcomplex *c__, int *ldc, double *scale, 
                       int *info);

extern int C2F(ztrti2)(char *uplo, char *diag, int *n, 
                       dcomplex *a, int *lda, int *info);

extern int C2F(ztrtri)(char *uplo, char *diag, int *n, 
                       dcomplex *a, int *lda, int *info);

extern int C2F(ztrtrs)(char *uplo, char *trans, char *diag, int *n, 
                       int *nrhs, dcomplex *a, int *lda, dcomplex *b, 
                       int *ldb, int *info);

extern int C2F(ztrttf)(char *transr, char *uplo, int *n, 
                       dcomplex *a, int *lda, dcomplex *arf, int *info);

extern int C2F(ztrttp)(char *uplo, int *n, dcomplex *a, 
                       int *lda, dcomplex *ap, int *info);

extern int C2F(ztzrqf)(int *m, int *n, dcomplex *a, 
                       int *lda, dcomplex *tau, int *info);

extern int C2F(ztzrzf)(int *m, int *n, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work, int *lwork, 
                       int *info);

extern int C2F(zung2l)(int *m, int *n, int *k, 
                       dcomplex *a, int *lda, dcomplex *tau, dcomplex *
                       work, int *info);

extern int C2F(zung2r)(int *m, int *n, int *k, 
                       dcomplex *a, int *lda, dcomplex *tau, dcomplex *
                       work, int *info);

extern int C2F(zungbr)(char *vect, int *m, int *n, int *k, 
                       dcomplex *a, int *lda, dcomplex *tau, dcomplex *
                       work, int *lwork, int *info);

extern int C2F(zunghr)(int *n, int *ilo, int *ihi, 
                       dcomplex *a, int *lda, dcomplex *tau, dcomplex *
                       work, int *lwork, int *info);

extern int C2F(zungl2)(int *m, int *n, int *k, 
                       dcomplex *a, int *lda, dcomplex *tau, dcomplex *
                       work, int *info);

extern int C2F(zunglq)(int *m, int *n, int *k, 
                       dcomplex *a, int *lda, dcomplex *tau, dcomplex *
                       work, int *lwork, int *info);

extern int C2F(zungql)(int *m, int *n, int *k, 
                       dcomplex *a, int *lda, dcomplex *tau, dcomplex *
                       work, int *lwork, int *info);

extern int C2F(zungqr)(int *m, int *n, int *k, 
                       dcomplex *a, int *lda, dcomplex *tau, dcomplex *
                       work, int *lwork, int *info);

extern int C2F(zungr2)(int *m, int *n, int *k, 
                       dcomplex *a, int *lda, dcomplex *tau, dcomplex *
                       work, int *info);

extern int C2F(zungrq)(int *m, int *n, int *k, 
                       dcomplex *a, int *lda, dcomplex *tau, dcomplex *
                       work, int *lwork, int *info);

extern int C2F(zungtr)(char *uplo, int *n, dcomplex *a, 
                       int *lda, dcomplex *tau, dcomplex *work, int *lwork, 
                       int *info);

extern int C2F(zunm2l)(char *side, char *trans, int *m, int *n, 
                       int *k, dcomplex *a, int *lda, dcomplex *tau, 
                       dcomplex *c__, int *ldc, dcomplex *work, int *info);

extern int C2F(zunm2r)(char *side, char *trans, int *m, int *n, 
                       int *k, dcomplex *a, int *lda, dcomplex *tau, 
                       dcomplex *c__, int *ldc, dcomplex *work, int *info);

extern int C2F(zunmbr)(char *vect, char *side, char *trans, int *m, 
                       int *n, int *k, dcomplex *a, int *lda, dcomplex 
                       *tau, dcomplex *c__, int *ldc, dcomplex *work, int *
                       lwork, int *info);

extern int C2F(zunmhr)(char *side, char *trans, int *m, int *n, 
                       int *ilo, int *ihi, dcomplex *a, int *lda, 
                       dcomplex *tau, dcomplex *c__, int *ldc, dcomplex *
                       work, int *lwork, int *info);

extern int C2F(zunml2)(char *side, char *trans, int *m, int *n, 
                       int *k, dcomplex *a, int *lda, dcomplex *tau, 
                       dcomplex *c__, int *ldc, dcomplex *work, int *info);

extern int C2F(zunmlq)(char *side, char *trans, int *m, int *n, 
                       int *k, dcomplex *a, int *lda, dcomplex *tau, 
                       dcomplex *c__, int *ldc, dcomplex *work, int *lwork, 
                       int *info);

extern int C2F(zunmql)(char *side, char *trans, int *m, int *n, 
                       int *k, dcomplex *a, int *lda, dcomplex *tau, 
                       dcomplex *c__, int *ldc, dcomplex *work, int *lwork, 
                       int *info);

extern int C2F(zunmqr)(char *side, char *trans, int *m, int *n, 
                       int *k, dcomplex *a, int *lda, dcomplex *tau, 
                       dcomplex *c__, int *ldc, dcomplex *work, int *lwork, 
                       int *info);

extern int C2F(zunmr2)(char *side, char *trans, int *m, int *n, 
                       int *k, dcomplex *a, int *lda, dcomplex *tau, 
                       dcomplex *c__, int *ldc, dcomplex *work, int *info);

extern int C2F(zunmr3)(char *side, char *trans, int *m, int *n, 
                       int *k, int *l, dcomplex *a, int *lda, dcomplex 
                       *tau, dcomplex *c__, int *ldc, dcomplex *work, int *
                       info);

extern int C2F(zunmrq)(char *side, char *trans, int *m, int *n, 
                       int *k, dcomplex *a, int *lda, dcomplex *tau, 
                       dcomplex *c__, int *ldc, dcomplex *work, int *lwork, 
                       int *info);

extern int C2F(zunmrz)(char *side, char *trans, int *m, int *n, 
                       int *k, int *l, dcomplex *a, int *lda, dcomplex 
                       *tau, dcomplex *c__, int *ldc, dcomplex *work, int *
                       lwork, int *info);

extern int C2F(zunmtr)(char *side, char *uplo, char *trans, int *m, 
                       int *n, dcomplex *a, int *lda, dcomplex *tau, 
                       dcomplex *c__, int *ldc, dcomplex *work, int *lwork, 
                       int *info);

extern int C2F(zupgtr)(char *uplo, int *n, dcomplex *ap, 
                       dcomplex *tau, dcomplex *q, int *ldq, dcomplex *
                       work, int *info);

extern int C2F(zupmtr)(char *side, char *uplo, char *trans, int *m, 
                       int *n, dcomplex *ap, dcomplex *tau, dcomplex *c__, 
                       int *ldc, dcomplex *work, int *info);

double C2F(dsecnd)(void);

extern int C2F(ilaver)(int *vers_major__, int *vers_minor__, int *vers_patch__);

int C2F(lsame)(char *ca, char *cb);

double C2F(second)(void);

double C2F(slamch)(char *cmach);

extern int C2F(slamc1)(int *beta, int *t, int *rnd, int *ieee1);

extern int C2F(slamc2)(int *beta, int *t, int *rnd, double *
                       eps, int *emin, double *rmin, int *emax, double *rmax);

extern double C2F(slamc3)(double *a, double *b);

extern int C2F(slamc4)(int *emin, double *start, int *base);

extern int C2F(slamc5)(int *beta, int *p, int *emin,
                       int *ieee, int *emax, double *rmax);


extern double C2F(dlamch)(char *cmach);

double pnl_dlamc3(double *a, double *b);

extern int C2F(ilaenv)(int *ispec, char *name__, char *opts, int *n1, 
                       int *n2, int *n3, int *n4);

#ifdef __cplusplus
}
#endif


#endif /* __CLAPACK_H */
