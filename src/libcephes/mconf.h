/*							mconf.h
 *
 *	Common include file for math routines
 *
 *
 *
 * SYNOPSIS:
 *
 * #include "mconf.h"
 *
 *
 *
 * DESCRIPTION:
 *
 * This file contains definitions for error codes that are
 * passed to the common error handling routine mtherr()
 * (which see).
 *
 * The file also includes a conditional assembly definition
 * for the type of computer arithmetic (IEEE, DEC, Motorola
 * IEEE, or UNKnown).
 * 
 * For Digital Equipment PDP-11 and VAX computers, certain
 * IBM systems, and others that use numbers with a 56-bit
 * significand, the symbol DEC should be defined.  In this
 * mode, most floating point constants are given as arrays
 * of octal integers to eliminate decimal to binary conversion
 * errors that might be introduced by the compiler.
 *
 * For little-endian computers, such as IBM PC, that follow the
 * IEEE Standard for Binary Floating Point Arithmetic (ANSI/IEEE
 * Std 754-1985), the symbol IBMPC should be defined.  These
 * numbers have 53-bit significands.  In this mode, constants
 * are provided as arrays of hexadecimal 16 bit integers.
 *
 * Big-endian IEEE format is denoted MIEEE.  On some RISC
 * systems such as Sun SPARC, double precision constants
 * must be stored on 8-byte address boundaries.  Since integer
 * arrays may be aligned differently, the MIEEE configuration
 * may fail on such machines.
 *
 * To accommodate other types of computer arithmetic, all
 * constants are also provided in a normal decimal radix
 * which one can hope are correctly converted to a suitable
 * format by the available C language compiler.  To invoke
 * this mode, define the symbol UNK.
 *
 * An important difference among these modes is a predefined
 * set of machine arithmetic constants for each.  The numbers
 * MACHEP (the machine roundoff error), MAXNUM (largest number
 * represented), and several other parameters are preset by
 * the configuration symbol.  Check the file const.c to
 * ensure that these values are correct for your computer.
 *
 * Configurations NANS, INFINITIES, MINUSZERO, and DENORMAL
 * may fail on many systems.  Verify that they are supposed
 * to work on your computer.
 */

/*
Cephes Math Library Release 2.3:  June, 1995
Copyright 1984, 1987, 1989, 1995 by Stephen L. Moshier
*/

#ifndef CEPHES_MCONF_H
#define CEPHES_MCONF_H

#include "pnl/pnl_mathtools.h"
/* to prefix several functions with pnl_sp */
#include "pnl_names.h"

/* Constant definitions for math error conditions
 */

/* Already defined in math.h */
/* #define DOMAIN		1	|+ argument domain error +| */
/* #define SING		2	|+ argument singularity +| */
/* #define OVERFLOW	3	|+ overflow range error +| */
/* #define UNDERFLOW	4	|+ underflow range error +| */
/* #define TLOSS		5	|+ total loss of precision +| */
/* #define PLOSS		6	|+ partial loss of precision +| */
#define TOOMANY         7       /* too many iterations */
#define MAXITER        500

#define EDOM		33
#define ERANGE		34

/* SciPy note: by defining UNK, we prevent the compiler from
 * casting integers to floating point numbers.  If the Endianness
 * is detected incorrectly, this causes problems on some platforms.
 */
#define UNK 1

/* For 12-byte long doubles on an i386, pad a 16-bit short 0
 * to the end of real constants initialized by integer arrays.
 *
 * #define XPD 0,
 *
 * Otherwise, the type is 10 bytes long and XPD should be
 * defined blank (e.g., Microsoft C).
 *
 * #define XPD
 */
#define XPD 0,

/* Define to support tiny denormal numbers, else undefine. */
#define DENORMAL 1

/* Define to distinguish between -0.0 and +0.0.  */
#define MINUSZERO 1

#include "cephes_protos.h"

#define EULER  = M_EULER;        /* Euler constant */

#ifdef UNK
#if 1
#define MACHEP   1.11022302462515654042E-16   /* 2**-53 */
#else
#define MACHEP   1.38777878078144567553E-17   /* 2**-56 */
#endif
#define UFLOWTHRESH   2.22507385850720138309E-308 /* 2**-1022 */
#ifdef DENORMAL
#define MAXLOG   7.09782712893383996732E2     /* log(MAXNUM) */
/* #define MINLOG  -7.44440071921381262314E2 */     /* log(2**-1074) */
#define MINLOG  -7.451332191019412076235E2     /* log(2**-1075) */
#else
#define MAXLOG   7.08396418532264106224E2     /* log 2**1022 */
#define MINLOG  -7.08396418532264106224E2     /* log 2**-1022 */
#endif
#define MAXNUM   1.79769313486231570815E308    /* 2**1024*(1-MACHEP) */
#define PI       M_PI       /* pi */
#define PIO2     M_PI_2       /* pi/2 */
#define PIO4     M_PI_4    /* pi/4 */
#define SQRT2    M_SQRT2       /* sqrt(2) */
#define SQRTH    M_SQRT1_2    /* sqrt(2)/2 */
#define LOG2E    1.4426950408889634073599     /* 1/log(2) */
#define SQ2OPI   M_SQRT2_PI  /* sqrt( 2/pi ) */
#define LOGE2    M_LN2    /* log(2) */
#define LOGSQ2   3.46573590279972654709E-1    /* log(2)/2 */
#define THPIO4   2.35619449019234492885       /* 3*pi/4 */
#define TWOOPI   M_2_PI /* 2/pi */
#ifdef MINUSZERO
#define NEGZERO  -0.0
#else
#define NEGZERO  0.0
#endif
#endif

#endif /* CEPHES_MCONF_H */
