#ifndef _PNL_MATHTOOLS_H
#define _PNL_MATHTOOLS_H 


#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "pnl/pnl_vector.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

extern int pnl_itrunc(double x);
extern long int pnl_ltrunc(double x);
extern double pnl_trunc(double x);
extern int pnl_iround(double x);
extern long int pnl_lround(double x);
extern double pnl_round(double x);
extern double pnl_fact(int n);
extern double pnl_lgamma (double x);
extern double pnl_tgamma (double x);
extern double pnl_pow_i (double x, int n);
extern double pnl_nan (void);
extern double pnl_posinf (void);
extern double pnl_neginf (void);
extern int pnl_isnan (double x);
extern int pnl_isfinite (double x);
extern int pnl_isinf (double x);
extern double pnl_acosh (double x);
extern double pnl_asinh (double x);
extern double pnl_atanh (double x);
extern double pnl_log1p (double x);
extern double pnl_expm1 (double x);
extern double pnl_cosm1 (double x);
extern double pnl_round (double x);

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
#ifndef SQRT_DBL_MAX
#define SQRT_DBL_MAX   1.3407807929942596e+154
#endif
#define DOUBLE_MAX DBL_MAX

#ifndef DBL_MIN
#define DBL_MIN        2.2250738585072014e-308
#endif
#ifndef SQRT_DBL_MIN
#define SQRT_DBL_MIN  1.4916681462400413e-154 
#endif


#ifndef DBL_EPSILON
#define DBL_EPSILON        2.2204460492503131e-16
#endif

#ifndef SQRT_DBL_EPSILON
#define SQRT_DBL_EPSILON  1.4901161193847656e-08
#endif


/* #define BIG_DOUBLE 1.0e6 */
/* #define PRECISION 1.0e-7 |+Precision for the localization of FD methods+| */
/* #define INC 1.0e-5 |+Relative Increment for Delta-Hedging+| */

#define MAXLOOPS 5000
#ifdef POW
#undef POW
#endif
#define POW(x,y) pow( (double) (x), (double) (y))

#ifndef PNL_HAVE_EXP10
#undef exp10
#define exp10(x) exp(10.0, x)
#endif

#ifndef PNL_HAVE_ROUND
#undef round
#define round(x) pnl_round(x)
#endif

#ifndef PNL_HAVE_TRUNC
#undef trunc
#define trunc(x) pnl_trunc(x)
#endif

#ifdef SQR
#undef SQR
#endif
#define SQR(x) ((x) * (x))

#ifdef CUB
#undef CUB
#endif
#define CUB(x) ((x) * (x) * (x))

#ifdef ABS
#undef ABS
#endif
#define ABS(x) ( ((x) >= 0) ? (x) : -(x) )

#define PNL_SIGN(x) ( ((x) < 0) ? (-1) : (1))

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

#define PNL_NAN (pnl_nan ())
#define PNL_INF INFINITY
#define PNL_POSINF INFINITY
#define PNL_NEGINF (-INFINITY)

#define PNL_IS_ODD(n) ((n) & 1) /* last bit is 1 */
#define PNL_IS_EVEN(n) (!PNL_IS_ODD(n)) /* last bit is 0 */

#define PNL_ALTERNATE(n) (((n)&1) ? -1 : 1) /* (-1)^n */


#define PNL_ERROR(msg, func) {fprintf(stderr, "%s in function %s \n", msg, func); abort();}
#define PNL_MESSAGE_ERROR(msg, func) { if (pnl_message_is_on()) fprintf(stderr, "%s in function %s \n", msg, func); }
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

/**
 * f: R --> R
 * The function  pointer returns f(x)
 */
typedef struct
{
  double (*F) (double x, void *params);
  void *params;
} PnlFunc ;
#define PNL_EVAL_FUNC(Fstruct, x) (*((Fstruct)->F))(x, (Fstruct)->params)

/**
 * f: R^2 --> R
 * The function pointer returns f(x)
 */
typedef struct 
{
  double (*F) (double x, double y, void *params);
  void *params;
} PnlFunc2D ;
#define PNL_EVAL_FUNC2D(Fstruct, x, y) (*((Fstruct)->F))(x, y, (Fstruct)->params)

/**
 * f: R --> R
 * The function pointer computes f(x) and Df(x) and stores them in fx
 * and dfx respectively
 */
typedef struct 
{
  void (*F) (double x, double *fx, double *dfx, void *params);
  void *params;
} PnlFuncDFunc ;
#define PNL_EVAL_FUNC_FDF(Fstruct, x, fx, dfx) (*((Fstruct)->F))(x, fx, dfx, (Fstruct)->params)

/**
 * f: R^n --> R
 * The function pointer returns f(x)
 */
typedef struct 
{
  double (*F) (const PnlVect *x, void *params);
  void *params;
} PnlRnFuncR ;
#define PNL_EVAL_RNFUNCR(Fstruct, x) (*((Fstruct)->F))(x, (Fstruct)->params)

/**
 * f: R^n --> R^m
 * The function pointer computes the vector f(x) and stores it in
 * fx (vector of size m)
 */
typedef struct 
{
  void (*F) (const PnlVect *x, PnlVect *fx, void *params);
  void *params;
} PnlRnFuncRm ;
#define PNL_EVAL_RNFUNCRM(Fstruct, x, fx) (*((Fstruct)->F))(x, fx, (Fstruct)->params)

/**
 * Synonymous of PnlRnFuncRm for f:R^n --> R^n 
 */
typedef PnlRnFuncRm PnlRnFuncRn;
#define PNL_EVAL_RNFUNCRN  PNL_EVAL_RNFUNCRM


/**
 * f: R^n --> R^m
 * The function pointer computes the vector f(x) and stores it in fx
 * (vector of size m) 
 * The Dfunction pointer computes the matrix Df(x) and stores it in dfx
 * (matrix of size m x n) 
 */
typedef struct 
{
  void (*F) (const PnlVect *x, PnlVect *fx, void *params);
  void (*DF) (const PnlVect *x, PnlMat *dfx, void *params);
  void (*FDF) (const PnlVect *x, PnlVect *fx, PnlMat *dfx, void *params);
  void *params;
} PnlRnFuncRmDFunc ;
#define PNL_EVAL_RNFUNCRM_DF(Fstruct, x, dfx) (*((Fstruct)->DF))(x, dfx, (Fstruct)->params)
#define PNL_EVAL_RNFUNCRM_FDF(Fstruct, x, fx, dfx) (*((Fstruct)->FDF))(x, fx, dfx, (Fstruct)->params)

#define PNL_EVAL_RNFUNCRM_F_DF(Fstruct, x, fx, dfx) 					\
      if ( (Fstruct)->FDF != NULL ) 									\
        { 																\
          PNL_EVAL_RNFUNCRN_FDF (Fstruct, x, fx, dfx); 					\
        } 																\
      else 																\
        { 																\
          PNL_EVAL_RNFUNCRN (Fstruct, x, fx); 							\
          PNL_EVAL_RNFUNCRN_DF (Fstruct, x, dfx); 						\
        }




/**
 * Synonymous of PnlRnFuncRmDFunc for f:R^n --> R^m
 */
typedef PnlRnFuncRmDFunc PnlRnFuncRnDFunc;
#define PNL_EVAL_RNFUNCRN_DF PNL_EVAL_RNFUNCRM_DF
#define PNL_EVAL_RNFUNCRN_FDF PNL_EVAL_RNFUNCRM_FDF
#define PNL_EVAL_RNFUNCRN_F_DF PNL_EVAL_RNFUNCRM_F_DF

/**
 * ODE functions
 * yp_i (t) = dy_i(t)/dt
 *
 */
typedef struct
{
  /** 
   * yp_i (t) = dy_i(t)/dt
   * 
   * @param[in] neqn number of equations
   * @param[in] t the point at which to evaluate the function
   * @param[in] y a C array of size neqn
   * @param[out] yp a C array of size neqn (must be already allocated).
   * Contains on output dy(t)/dt
   * @param[in] params a generic pointer containing extra parameters involved
   * in the equation
   */
  void (*F) (int neqn, double t, const double *y, double *yp, void *params);
  int neqn; /*!< number of equations of the system */
  void *params;
} PnlODEFunc ;
#define PNL_EVAL_ODEFUNC(Fstruct, t, y, yp) (*((Fstruct)->F))((Fstruct)->neqn, t, y, yp, (Fstruct)->params)

extern void pnl_qsort (void *a, int n, int es, int lda, int *t, int ldt, int use_index, int (*cmp)(void const *, void const *));

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_MATHTOOLS_H */ 
