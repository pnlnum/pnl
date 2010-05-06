#ifndef _PNL_MATHTOOLS_H
#define _PNL_MATHTOOLS_H 


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "pnl_vector.h"

typedef unsigned char boolean;

extern int intapprox (double s);
extern double trunc(double x);
extern double Cnp(int n, int p);
extern double pnl_fact(int n);
extern double lgamma (double x);
extern double tgamma (double x);
extern double pnl_pow_i (double x, int n);
extern double pnl_nan (void);
extern double pnl_posinf (void);
extern double pnl_neginf (void);
extern int pnl_isnan (double x);
extern int pnl_isfinite (double x);
extern int pnl_isinf (double x);


#define false 0
#undef FALSE
#define FALSE 0
#define true 1
#undef TRUE
#define TRUE 1


/* The following constants are normally defined in /usr/include/{math.h,
 * values.h, limits.h}. In
 * case they are not
 */

#ifndef M_E
#define M_E            2.7182818284590452354   /* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E        1.4426950408889634074   /* log_2 e */
#endif

#ifndef M_LOG10E
#define M_LOG10E       0.43429448190325182765  /* log_10 e */
#endif

#ifndef M_LN2
#define M_LN2          0.69314718055994530942  /* log_e 2 */
#endif

#ifndef M_LN10
#define M_LN10         2.30258509299404568402  /* log_e 10 */
#endif

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2         1.57079632679489661923  /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4         0.78539816339744830962  /* pi/4 */
#endif

#ifndef M_1_PI
#define M_1_PI         0.31830988618379067154  /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI         0.63661977236758134308  /* 2/pi */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI     1.12837916709551257390  /* 2/sqrt(pi) */
#endif

#ifndef M_SQRT2
#define M_SQRT2        1.41421356237309504880  /* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2      0.70710678118654752440  /* 1/sqrt(2) */
#endif

#ifndef M_EULER
#define M_EULER        0.5772156649015328606
#endif

#ifndef M_1_SQRT2PI
#define M_1_SQRT2PI    3.9894228040143267791E-1 /* 1/sqrt(2 pi) */
#endif

#ifndef M_SQRT2_PI
#define M_SQRT2_PI     0.79788456080286535588   /* sqrt(2/pi) */
#endif

#ifndef M_2PI
#define M_2PI          6.28318530717958623199   /* 2 pi */
#endif

#ifndef M_SQRT2PI 
#define M_SQRT2PI      2.506628274631000502415  /* sqrt(2*pi) */
#endif

#ifndef INT_MAX
#define INT_MAX        2147483647
#endif
#define MAX_INT INT_MAX

#ifndef DBL_MAX
#define DBL_MAX		   1.79769313486231470e+308
#endif
#define DOUBLE_MAX DBL_MAX

#ifndef DBL_EPSILON
#define DBL_EPSILON        2.2204460492503131e-16
#endif

#define BIG_DOUBLE 1.0e6
#define PRECISION 1.0e-7 /*Precision for the localization of FD methods*/
#define INC 1.0e-5 /*Relative Increment for Delta-Hedging*/

#define MAXLOOPS 5000
#ifdef POW
#undef POW
#endif
#define POW(x,y) pow( (double) (x), (double) (y))

#ifdef SQR
#undef SQR
#endif
#define SQR(x) pnl_pow_i(x, 2)

#ifdef CUB
#undef CUB
#endif
#define CUB(x) pnl_pow_i(x, 3)

#ifdef ABS
#undef ABS
#endif
#define ABS(x) ( ((x) >= 0) ? (x) : -(x) )

/* MAX and MIN already defined in <sys/param.h>  */
/* if we are compiling for Windows (Visual or mingw32
 * cross-compiling), sys/param.h is not included */
#ifndef MAX
#define MAX(A,B) ( (A) > (B) ? (A):(B) )
#endif
#ifndef MIN
#define MIN(A,B) ( (A) < (B) ? (A):(B) )
#endif

#ifndef NAN
#define NAN (pnl_nan ())
#endif

#ifndef INFINITY
#define INFINITY (pnl_posinf ())
#endif

#define PNL_POSINF INFINITY
#define PNL_NEGINF (-INFINITY)

#ifndef isinf
#define isinf pnl_isinf
#endif
#ifndef isnan
#define isnan pnl_isnan
#endif
#ifndef isfinite 
#define isfinite pnl_isfinite
#endif


#define PNL_IS_ODD(n) ((n) & 1) /* last bit is 1 */
#define PNL_IS_EVEN(n) ((n) ^ 1) /* last bit is 0 */

#define PNL_ALTERNATE(n) (((n)&1) ? -1 : 1) /* (-1)^n */


#define PNL_ERROR(msg, func) {fprintf(stderr, "%s in function %s \n", msg, func); abort();}
#ifdef PNL_RANGE_CHECK_OFF
#define PNL_CHECK(eq, msg, func)
#else
#define PNL_CHECK(eq, msg, func)                                    \
  if (eq) {                                                         \
    fprintf(stderr, "%s in function %s \n", msg, func); abort();    \
  }
#endif

#define OK 0
#define WRONG 1
#define FAIL 1 /* synonym of WRONG (more menaningful) */

/*
 * f: R --> R
 * The function returns it value at x
 */
typedef struct {
  double (*function) (double x, void *params);
  void *params;
} PnlFunc ;
#define PNL_EVAL_FUNC(F, x) (*((F)->function))(x, (F)->params)

/*
 * f: R^2 --> R
 * The function returns it value at x
 */
typedef struct {
  double (*function) (double x, double y, void *params);
  void *params;
} PnlFunc2D ;
#define PNL_EVAL_FUNC2D(F, x, y) (*((F)->function))(x, y, (F)->params)

/*
 * f: R --> R
 * The function computes its value and derivative at x and stores them in f
 * and df respectively
 */
typedef struct {
  void (*function) (double x, double *f, double *df, void *params);
  void *params;
} PnlFuncDFunc ;
#define PNL_EVAL_FUNC_DFUNC(F, x, f, df) (*((F)->function))(x, f, df, (F)->params)

/*
 * f: R^n --> R
 * The function returns its value at x, which is a vector
 */
typedef struct {
  double (*function) (const PnlVect *x, void *params);
  void *params;
} PnlRnFuncR ;
#define PNL_EVAL_RNFUNCR(F, x) (*((F)->function))(x, (F)->params)

/*
 * f: R^n --> R^m
 * The function computes its value at x (vector of size n) and stores it in
 * f (vector of size m)
 */
typedef struct {
  void (*function) (const PnlVect *x, PnlVect *f, void *params);
  void *params;
} PnlRnFuncRm ;
#define PNL_EVAL_RNFUNCRM(F, x, f) (*((F)->function))(x, f, (F)->params)

/*
 * f: R^n --> R^n
 * The function computes its value at x (vector of size n) and stores it in
 * res (vector of size n)
 */
typedef struct {
  void (*function) (const PnlVect *x, PnlVect *f, void *params);
  void *params;
} PnlRnFuncRn ;
#define PNL_EVAL_RNFUNCRN(F, x, f) (*((F)->function))(x, f, (F)->params)


/*
 * f: R^n --> R^n
 * The function computes its value and gradient at x (vector of size n) 
 * and stores them in f (vector of size n) and df (matrix of size nxn) 
 * respectively
 */
typedef struct {
  void (*function) (const PnlVect *x, PnlVect *f, PnlMat *df, void *params);
  void *params;
} PnlRnFuncRnDFunc ;
#define PNL_EVAL_RNFUNCRN_DFUNC(F, x, f, df) (*((F)->function))(x, f, df, (F)->params)



extern void pnl_qsort (void *a, int n, int es, int lda, int *t, int ldt, int use_index, int (*cmp)(void const *, void const *));

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_MATHTOOLS_H */ 
