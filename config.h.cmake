#ifndef _PNL_CONFIG_H
#define _PNL_CONFIG_H
/* config.h.in.  Generated from cmake */

/* Define if a Fortran is detected */
#cmakedefine HAVE_FORTRAN_COMPILER

/* Define to dummy 'main' function (if any) required to link to the Fortran
   libraries. */
#cmakedefine F77_DUMMY_MAIN

/* Define to 1 if you have the 'exp10' function. */
#cmakedefine HAVE_EXP10

/* Define to 1 if you have the 'finite' function. */
#cmakedefine HAVE_FINITE

/* Define if you have inline */
#cmakedefine HAVE_INLINE

/* Define keyword for declaring inline functions */
#define PNL_INLINE_DECL @PNL_INLINE_DECL@

/* Define keyword for defining inline functions */
#define PNL_INLINE_FUNC @PNL_INLINE_FUNC@

/* Define to 1 if you have the 'isfinite' function. */
#cmakedefine HAVE_ISFINITE

/* Define to 1 if you have the 'isinf' function. */
#cmakedefine HAVE_ISINF

/* Define to 1 if you have the 'isnan' function. */
#cmakedefine HAVE_ISNAN

/* Define to 1 if you have the 'lgamma' function. */
#cmakedefine HAVE_LGAMMA

/* Define to 1 if you have the 'tgamma' function. */
#cmakedefine HAVE_TGAMMA

/* Define to 1 if you have the 'trunc' function. */
#cmakedefine HAVE_TRUNC

/* Define to 1 if you have 'dpstrf' function */
#cmakedefine HAVE_DPSTRF

/* turn off range checking by default internally */
#cmakedefine PNL_RANGE_CHECK_OFF

/* Runnning on Cygwin */
#cmakedefine _CYGWIN

/* Running on a Win32 system */
#cmakedefine _WIN32

/* Define to 1 if you use internal Blas */
#cmakedefine USE_INTERNAL_BLAS

#endif /* _PNL_CONFIG_H */

