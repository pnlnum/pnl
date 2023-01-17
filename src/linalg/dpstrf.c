/*
 * Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by  the Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License  along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

/*
 * This file contains two functions from Lapack:
 *    dpstrf and dpstf2 (and dmaxloc)
 * to compute a Block Cholesky factorization for non
 * definite matrices. It is only used if PNL_HAVE_DPSTRF is
 * not defined by the configure script.
 * Note that this file is useless if the internal Blas & Lapack
 * libraries are used because these two functions are alredy
 * included there.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "pnl/pnl_config.h"
#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_machine.h"
#include "pnl/clapack.h"

#ifndef PNL_HAVE_DPSTRF


/**
 * Return PNL_TRUE if x is NaN and PNL_FALSE otherwise
 *
 * @param x the address of a real number
 *
 * @return  PNL_TRUE or PNL_FALSE
 */
int disnan_(double *x)
{
  return pnl_isnan(*x) ? PNL_TRUE : PNL_FALSE;
}

/* dpstrf.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
   on Microsoft Windows system, link with libf2c.lib;
   on Linux or Unix systems, link with .../path/to/libf2c.a -lm
   or, if you install libf2c.a in a standard place, with -lf2c -lm
   -- in that order, at the end of the command line, as in
   cc *.o -lf2c -lm
   Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

http://www.netlib.org/f2c/libf2c.zip
*/

/* Table of constant values */

static int c__1 = 1;
static int c_n1 = -1;
static double c_b22 = -1.;
static double c_b24 = 1.;

int C2F(dpstrf)(char *uplo, int *n, double *a, int *
                lda, int *piv, int *rank, double *tol, double *work,
                int *info)
{
  /* System generated locals */
  int a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
  double d__1;

  /* Builtin functions */
  double sqrt(double);

  /* Local variables */
  int i__, j, k, maxlocvar, jb, nb;
  double ajj;
  int pvt;
  extern  int dscal_(int *, double *, double *,
                     int *);
  extern int lsame_(char *, char *);
  extern  int dgemv_(char *, int *, int *,
                     double *, double *, int *, double *, int *,
                     double *, double *, int *);
  double dtemp;
  int itemp;
  extern  int dswap_(int *, double *, int *,
                     double *, int *);
  double dstop;
  int upper;
  extern  int dsyrk_(char *, char *, int *, int *,
                     double *, double *, int *, double *, double *,
                     int *);
  extern int dpstf2_(char *, int *, double *, int *, int *, int *,
                     double *, double *, int *);
  extern double dlamch_(char *);
  extern  int xerbla_(char *, int *);
  extern int ilaenv_(int *, char *, char *, int *, int *,
                     int *, int *);
  extern int dmaxloc_(double *, int *);


  /* -- LAPACK routine (version 3.2) -- */
  /*     Craig Lucas, University of Manchester / NAG Ltd. */
  /*     October, 2008 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DPSTRF computes the Cholesky factorization with complete */
  /*  pivoting of a float symmetric positive semidefinite matrix A. */

  /*  The factorization has the form */
  /*     P' * A * P = U' * U ,  if UPLO = 'U', */
  /*     P' * A * P = L  * L',  if UPLO = 'L', */
  /*  where U is an upper triangular matrix and L is lower triangular, and */
  /*  P is stored as vector PIV. */

  /*  This algorithm does not attempt to check that A is positive */
  /*  semidefinite. This version of the algorithm calls level 3 BLAS. */

  /*  Arguments */
  /*  ========= */

  /*  UPLO    (input) CHARACTER*1 */
  /*          Specifies whether the upper or lower triangular part of the */
  /*          symmetric matrix A is stored. */
  /*          = 'U':  Upper triangular */
  /*          = 'L':  Lower triangular */

  /*  N       (input) INTEGER */
  /*          The order of the matrix A.  N >= 0. */

  /*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
  /*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
  /*          n by n upper triangular part of A contains the upper */
  /*          triangular part of the matrix A, and the strictly lower */
  /*          triangular part of A is not referenced.  If UPLO = 'L', the */
  /*          leading n by n lower triangular part of A contains the lower */
  /*          triangular part of the matrix A, and the strictly upper */
  /*          triangular part of A is not referenced. */

  /*          On exit, if INFO = 0, the factor U or L from the Cholesky */
  /*          factorization as above. */

  /*  LDA     (input) INTEGER */
  /*          The leading dimension of the array A.  LDA >= MAX(1,N). */

  /*  PIV     (output) INTEGER array, dimension (N) */
  /*          PIV is such that the nonzero entries are P( PIV(K), K ) = 1. */

  /*  RANK    (output) INTEGER */
  /*          The rank of A given by the number of steps the algorithm */
  /*          completed. */

  /*  TOL     (input) DOUBLE PRECISION */
  /*          User defined tolerance. If TOL < 0, then N*U*MAX( A(K,K) ) */
  /*          will be used. The algorithm terminates at the (K-1)st step */
  /*          if the pivot <= TOL. */

  /*  WORK    DOUBLE PRECISION array, dimension (2*N) */
  /*          Work space. */

  /*  INFO    (output) INTEGER */
  /*          < 0: If INFO = -K, the K-th argument had an illegal value, */
  /*          = 0: algorithm completed successfully, and */
  /*          > 0: the matrix A is either rank deficient with computed rank */
  /*               as returned in RANK, or is indefinite.  See Section 7 of */
  /*               LAPACK Working Note #161 for further information. */

  /*  ===================================================================== */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input parameters. */

  /* Parameter adjustments */
  --work;
  --piv;
  a_dim1 = *lda;
  a_offset = 1 + a_dim1;
  a -= a_offset;

  /* Function Body */
  *info = 0;
  upper = lsame_(uplo, "U");
  if (! upper && ! lsame_(uplo, "L"))
    {
      *info = -1;
    }
  else if (*n < 0)
    {
      *info = -2;
    }
  else if (*lda < MAX(1, *n))
    {
      *info = -4;
    }
  if (*info != 0)
    {
      i__1 = -(*info);
      xerbla_("DPSTRF", &i__1);
      return 0;
    }

  /* Quick return if possible */

  if (*n == 0)
    {
      return 0;
    }

  /* Get block size */

  nb = ilaenv_(&c__1, "DPOTRF", uplo, n, &c_n1, &c_n1, &c_n1);
  if (nb <= 1 || nb >= *n)
    {

      /* Use unblocked code */

      dpstf2_(uplo, n, &a[a_dim1 + 1], lda, &piv[1], rank, tol, &work[1],
              info);
      goto L200;

    }
  else
    {

      /* Initialize PIV */

      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__)
        {
          piv[i__] = i__;
          /* L100: */
        }

      /* Compute stopping value */

      pvt = 1;
      ajj = a[pvt + pvt * a_dim1];
      i__1 = *n;
      for (i__ = 2; i__ <= i__1; ++i__)
        {
          if (a[i__ + i__ * a_dim1] > ajj)
            {
              pvt = i__;
              ajj = a[pvt + pvt * a_dim1];
            }
        }
      if (ajj == 0. || disnan_(&ajj))
        {
          *rank = 0;
          *info = 1;
          goto L200;
        }

      /* Compute stopping value if not supplied */

      if (*tol < 0.)
        {
          dstop = *n * dlamch_("Epsilon") * ajj;
        }
      else
        {
          dstop = *tol;
        }


      if (upper)
        {

          /* Compute the Cholesky factorization P' * A * P = U' * U */

          i__1 = *n;
          i__2 = nb;
          for (k = 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2)
            {

              /* Account for last block not being NB wide */

              /* Computing MIN */
              i__3 = nb, i__4 = *n - k + 1;
              jb = MIN(i__3, i__4);

              /* Set relevant part of first half of WORK to zero, */
              /* holds dot products */

              i__3 = *n;
              for (i__ = k; i__ <= i__3; ++i__)
                {
                  work[i__] = 0.;
                  /* L110: */
                }

              i__3 = k + jb - 1;
              for (j = k; j <= i__3; ++j)
                {

                  /* Find pivot, test for exit, else swap rows and columns */
                  /* Update dot products, compute possible pivots which are */
                  /* stored in the second half of WORK */

                  i__4 = *n;
                  for (i__ = j; i__ <= i__4; ++i__)
                    {

                      if (j > k)
                        {
                          /* Computing 2nd power */
                          d__1 = a[j - 1 + i__ * a_dim1];
                          work[i__] += d__1 * d__1;
                        }
                      work[*n + i__] = a[i__ + i__ * a_dim1] - work[i__];

                      /* L120: */
                    }

                  if (j > 1)
                    {
                      maxlocvar = (*n << 1) - (*n + j) + 1;
                      itemp = dmaxloc_(&work[*n + j], &maxlocvar);
                      pvt = itemp + j - 1;
                      ajj = work[*n + pvt];
                      if (ajj <= dstop || disnan_(&ajj))
                        {
                          a[j + j * a_dim1] = ajj;
                          goto L190;
                        }
                    }

                  if (j != pvt)
                    {

                      /* Pivot PNL_OK, so can now swap pivot rows and columns */

                      a[pvt + pvt * a_dim1] = a[j + j * a_dim1];
                      i__4 = j - 1;
                      dswap_(&i__4, &a[j * a_dim1 + 1], &c__1, &a[pvt *
                             a_dim1 + 1], &c__1);
                      if (pvt < *n)
                        {
                          i__4 = *n - pvt;
                          dswap_(&i__4, &a[j + (pvt + 1) * a_dim1], lda, &a[
                                   pvt + (pvt + 1) * a_dim1], lda);
                        }
                      i__4 = pvt - j - 1;
                      dswap_(&i__4, &a[j + (j + 1) * a_dim1], lda, &a[j + 1
                             + pvt * a_dim1], &c__1);

                      /* Swap dot products and PIV */

                      dtemp = work[j];
                      work[j] = work[pvt];
                      work[pvt] = dtemp;
                      itemp = piv[pvt];
                      piv[pvt] = piv[j];
                      piv[j] = itemp;
                    }

                  ajj = sqrt(ajj);
                  a[j + j * a_dim1] = ajj;

                  /* Compute elements J+1:N of row J. */

                  if (j < *n)
                    {
                      i__4 = j - k;
                      i__5 = *n - j;
                      dgemv_("Trans", &i__4, &i__5, &c_b22, &a[k + (j + 1) *
                             a_dim1], lda, &a[k + j * a_dim1], &c__1, &
                             c_b24, &a[j + (j + 1) * a_dim1], lda);
                      i__4 = *n - j;
                      d__1 = 1. / ajj;
                      dscal_(&i__4, &d__1, &a[j + (j + 1) * a_dim1], lda);
                    }

                  /* L130: */
                }

              /* Update trailing matrix, J already incremented */

              if (k + jb <= *n)
                {
                  i__3 = *n - j + 1;
                  dsyrk_("Upper", "Trans", &i__3, &jb, &c_b22, &a[k + j *
                         a_dim1], lda, &c_b24, &a[j + j * a_dim1], lda);
                }

              /* L140: */
            }

        }
      else
        {

          /* Compute the Cholesky factorization P' * A * P = L * L' */

          i__2 = *n;
          i__1 = nb;
          for (k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1)
            {

              /* Account for last block not being NB wide */

              /* Computing MIN */
              i__3 = nb, i__4 = *n - k + 1;
              jb = MIN(i__3, i__4);

              /* Set relevant part of first half of WORK to zero, */
              /* holds dot products */

              i__3 = *n;
              for (i__ = k; i__ <= i__3; ++i__)
                {
                  work[i__] = 0.;
                  /* L150: */
                }

              i__3 = k + jb - 1;
              for (j = k; j <= i__3; ++j)
                {

                  /* Find pivot, test for exit, else swap rows and columns */
                  /* Update dot products, compute possible pivots which are */
                  /* stored in the second half of WORK */

                  i__4 = *n;
                  for (i__ = j; i__ <= i__4; ++i__)
                    {

                      if (j > k)
                        {
                          /* Computing 2nd power */
                          d__1 = a[i__ + (j - 1) * a_dim1];
                          work[i__] += d__1 * d__1;
                        }
                      work[*n + i__] = a[i__ + i__ * a_dim1] - work[i__];

                      /* L160: */
                    }

                  if (j > 1)
                    {
                      maxlocvar = (*n << 1) - (*n + j) + 1;
                      itemp = dmaxloc_(&work[*n + j], &maxlocvar);
                      pvt = itemp + j - 1;
                      ajj = work[*n + pvt];
                      if (ajj <= dstop || disnan_(&ajj))
                        {
                          a[j + j * a_dim1] = ajj;
                          goto L190;
                        }
                    }

                  if (j != pvt)
                    {

                      /* Pivot PNL_OK, so can now swap pivot rows and columns */

                      a[pvt + pvt * a_dim1] = a[j + j * a_dim1];
                      i__4 = j - 1;
                      dswap_(&i__4, &a[j + a_dim1], lda, &a[pvt + a_dim1],
                             lda);
                      if (pvt < *n)
                        {
                          i__4 = *n - pvt;
                          dswap_(&i__4, &a[pvt + 1 + j * a_dim1], &c__1, &a[
                                   pvt + 1 + pvt * a_dim1], &c__1);
                        }
                      i__4 = pvt - j - 1;
                      dswap_(&i__4, &a[j + 1 + j * a_dim1], &c__1, &a[pvt +
                             (j + 1) * a_dim1], lda);

                      /* Swap dot products and PIV */

                      dtemp = work[j];
                      work[j] = work[pvt];
                      work[pvt] = dtemp;
                      itemp = piv[pvt];
                      piv[pvt] = piv[j];
                      piv[j] = itemp;
                    }

                  ajj = sqrt(ajj);
                  a[j + j * a_dim1] = ajj;

                  /* Compute elements J+1:N of column J. */

                  if (j < *n)
                    {
                      i__4 = *n - j;
                      i__5 = j - k;
                      dgemv_("No Trans", &i__4, &i__5, &c_b22, &a[j + 1 + k
                             * a_dim1], lda, &a[j + k * a_dim1], lda, &
                             c_b24, &a[j + 1 + j * a_dim1], &c__1);
                      i__4 = *n - j;
                      d__1 = 1. / ajj;
                      dscal_(&i__4, &d__1, &a[j + 1 + j * a_dim1], &c__1);
                    }

                  /* L170: */
                }

              /* Update trailing matrix, J already incremented */

              if (k + jb <= *n)
                {
                  i__3 = *n - j + 1;
                  dsyrk_("Lower", "No Trans", &i__3, &jb, &c_b22, &a[j + k *
                         a_dim1], lda, &c_b24, &a[j + j * a_dim1], lda);
                }

              /* L180: */
            }

        }
    }

  /* Ran to completion, A has full rank */

  *rank = *n;

  goto L200;
L190:

  /* Rank is the number of steps completed.  Set INFO = 1 to signal */
  /* that the factorization cannot be used to solve a system. */

  *rank = j - 1;
  *info = 1;

L200:
  return 0;

  /* End of DPSTRF */

} /* dpstrf_ */

/* dpstf2.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
   on Microsoft Windows system, link with libf2c.lib;
   on Linux or Unix systems, link with .../path/to/libf2c.a -lm
   or, if you install libf2c.a in a standard place, with -lf2c -lm
   -- in that order, at the end of the command line, as in
   cc *.o -lf2c -lm
   Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

http://www.netlib.org/f2c/libf2c.zip
*/


/* Table of constant values */

static double c_b16 = -1.;
static double c_b18 = 1.;

int C2F(dpstf2)(char *uplo, int *n, double *a, int *
                lda, int *piv, int *rank, double *tol, double *work,
                int *info)
{
  /* System generated locals */
  int a_dim1, a_offset, i__1, i__2, i__3;
  double d__1;

  /* Builtin functions */
  double sqrt(double);

  /* Local variables */
  int i__, j, maxlocval;
  double ajj;
  int pvt;
  extern  int dscal_(int *, double *, double *,
                     int *);
  extern int lsame_(char *, char *);
  extern  int dgemv_(char *, int *, int *,
                     double *, double *, int *, double *, int *,
                     double *, double *, int *);
  double dtemp;
  int itemp;
  extern  int dswap_(int *, double *, int *,
                     double *, int *);
  double dstop;
  int upper;
  extern double dlamch_(char *);
  extern int disnan_(double *);
  extern  int xerbla_(char *, int *);
  extern int dmaxloc_(double *, int *);


  /*  -- LAPACK PROTOTYPE routine (version 3.2) -- */
  /*     Craig Lucas, University of Manchester / NAG Ltd. */
  /*     October, 2008 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DPSTF2 computes the Cholesky factorization with complete */
  /*  pivoting of a float symmetric positive semidefinite matrix A. */

  /*  The factorization has the form */
  /*     P' * A * P = U' * U ,  if UPLO = 'U', */
  /*     P' * A * P = L  * L',  if UPLO = 'L', */
  /*  where U is an upper triangular matrix and L is lower triangular, and */
  /*  P is stored as vector PIV. */

  /*  This algorithm does not attempt to check that A is positive */
  /*  semidefinite. This version of the algorithm calls level 2 BLAS. */

  /*  Arguments */
  /*  ========= */

  /*  UPLO    (input) CHARACTER*1 */
  /*          Specifies whether the upper or lower triangular part of the */
  /*          symmetric matrix A is stored. */
  /*          = 'U':  Upper triangular */
  /*          = 'L':  Lower triangular */

  /*  N       (input) INTEGER */
  /*          The order of the matrix A.  N >= 0. */

  /*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
  /*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
  /*          n by n upper triangular part of A contains the upper */
  /*          triangular part of the matrix A, and the strictly lower */
  /*          triangular part of A is not referenced.  If UPLO = 'L', the */
  /*          leading n by n lower triangular part of A contains the lower */
  /*          triangular part of the matrix A, and the strictly upper */
  /*          triangular part of A is not referenced. */

  /*          On exit, if INFO = 0, the factor U or L from the Cholesky */
  /*          factorization as above. */

  /*  PIV     (output) INTEGER array, dimension (N) */
  /*          PIV is such that the nonzero entries are P( PIV(K), K ) = 1. */

  /*  RANK    (output) INTEGER */
  /*          The rank of A given by the number of steps the algorithm */
  /*          completed. */

  /*  TOL     (input) DOUBLE PRECISION */
  /*          User defined tolerance. If TOL < 0, then N*U*MAX( A( K,K ) ) */
  /*          will be used. The algorithm terminates at the (K-1)st step */
  /*          if the pivot <= TOL. */

  /*  LDA     (input) INTEGER */
  /*          The leading dimension of the array A.  LDA >= MAX(1,N). */

  /*  WORK    DOUBLE PRECISION array, dimension (2*N) */
  /*          Work space. */

  /*  INFO    (output) INTEGER */
  /*          < 0: If INFO = -K, the K-th argument had an illegal value, */
  /*          = 0: algorithm completed successfully, and */
  /*          > 0: the matrix A is either rank deficient with computed rank */
  /*               as returned in RANK, or is indefinite.  See Section 7 of */
  /*               LAPACK Working Note #161 for further information. */

  /*  ===================================================================== */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input parameters */

  /* Parameter adjustments */
  --work;
  --piv;
  a_dim1 = *lda;
  a_offset = 1 + a_dim1;
  a -= a_offset;

  /* Function Body */
  *info = 0;
  upper = lsame_(uplo, "U");
  if (! upper && ! lsame_(uplo, "L"))
    {
      *info = -1;
    }
  else if (*n < 0)
    {
      *info = -2;
    }
  else if (*lda < MAX(1, *n))
    {
      *info = -4;
    }
  if (*info != 0)
    {
      i__1 = -(*info);
      xerbla_("DPSTF2", &i__1);
      return 0;
    }

  /* Quick return if possible */

  if (*n == 0)
    {
      return 0;
    }

  /* Initialize PIV */

  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      piv[i__] = i__;
      /* L100: */
    }

  /* Compute stopping value */

  pvt = 1;
  ajj = a[pvt + pvt * a_dim1];
  i__1 = *n;
  for (i__ = 2; i__ <= i__1; ++i__)
    {
      if (a[i__ + i__ * a_dim1] > ajj)
        {
          pvt = i__;
          ajj = a[pvt + pvt * a_dim1];
        }
    }
  if (ajj == 0. || disnan_(&ajj))
    {
      *rank = 0;
      *info = 1;
      goto L170;
    }

  /* Compute stopping value if not supplied */

  if (*tol < 0.)
    {
      dstop = *n * dlamch_("Epsilon") * ajj;
    }
  else
    {
      dstop = *tol;
    }

  /* Set first half of WORK to zero, holds dot products */

  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      work[i__] = 0.;
      /* L110: */
    }

  if (upper)
    {

      /* Compute the Cholesky factorization P' * A * P = U' * U */

      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
        {

          /* Find pivot, test for exit, else swap rows and columns */
          /* Update dot products, compute possible pivots which are */
          /* stored in the second half of WORK */

          i__2 = *n;
          for (i__ = j; i__ <= i__2; ++i__)
            {

              if (j > 1)
                {
                  /* Computing 2nd power */
                  d__1 = a[j - 1 + i__ * a_dim1];
                  work[i__] += d__1 * d__1;
                }
              work[*n + i__] = a[i__ + i__ * a_dim1] - work[i__];

              /* L120: */
            }

          if (j > 1)
            {
              maxlocval = (*n << 1) - (*n + j) + 1;
              itemp = dmaxloc_(&work[*n + j], &maxlocval);
              pvt = itemp + j - 1;
              ajj = work[*n + pvt];
              if (ajj <= dstop || disnan_(&ajj))
                {
                  a[j + j * a_dim1] = ajj;
                  goto L160;
                }
            }

          if (j != pvt)
            {

              /* Pivot PNL_OK, so can now swap pivot rows and columns */

              a[pvt + pvt * a_dim1] = a[j + j * a_dim1];
              i__2 = j - 1;
              dswap_(&i__2, &a[j * a_dim1 + 1], &c__1, &a[pvt * a_dim1 + 1],
                     &c__1);
              if (pvt < *n)
                {
                  i__2 = *n - pvt;
                  dswap_(&i__2, &a[j + (pvt + 1) * a_dim1], lda, &a[pvt + (
                           pvt + 1) * a_dim1], lda);
                }
              i__2 = pvt - j - 1;
              dswap_(&i__2, &a[j + (j + 1) * a_dim1], lda, &a[j + 1 + pvt *
                     a_dim1], &c__1);

              /* Swap dot products and PIV */

              dtemp = work[j];
              work[j] = work[pvt];
              work[pvt] = dtemp;
              itemp = piv[pvt];
              piv[pvt] = piv[j];
              piv[j] = itemp;
            }

          ajj = sqrt(ajj);
          a[j + j * a_dim1] = ajj;

          /* Compute elements J+1:N of row J */

          if (j < *n)
            {
              i__2 = j - 1;
              i__3 = *n - j;
              dgemv_("Trans", &i__2, &i__3, &c_b16, &a[(j + 1) * a_dim1 + 1]
                     , lda, &a[j * a_dim1 + 1], &c__1, &c_b18, &a[j + (j +
                         1) * a_dim1], lda);
              i__2 = *n - j;
              d__1 = 1. / ajj;
              dscal_(&i__2, &d__1, &a[j + (j + 1) * a_dim1], lda);
            }

          /* L130: */
        }

    }
  else
    {

      /* Compute the Cholesky factorization P' * A * P = L * L' */

      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
        {

          /* Find pivot, test for exit, else swap rows and columns */
          /* Update dot products, compute possible pivots which are */
          /* stored in the second half of WORK */

          i__2 = *n;
          for (i__ = j; i__ <= i__2; ++i__)
            {

              if (j > 1)
                {
                  /* Computing 2nd power */
                  d__1 = a[i__ + (j - 1) * a_dim1];
                  work[i__] += d__1 * d__1;
                }
              work[*n + i__] = a[i__ + i__ * a_dim1] - work[i__];

              /* L140: */
            }

          if (j > 1)
            {
              maxlocval = (*n << 1) - (*n + j) + 1;
              itemp = dmaxloc_(&work[*n + j], &maxlocval);
              pvt = itemp + j - 1;
              ajj = work[*n + pvt];
              if (ajj <= dstop || disnan_(&ajj))
                {
                  a[j + j * a_dim1] = ajj;
                  goto L160;
                }
            }

          if (j != pvt)
            {

              /* Pivot PNL_OK, so can now swap pivot rows and columns */

              a[pvt + pvt * a_dim1] = a[j + j * a_dim1];
              i__2 = j - 1;
              dswap_(&i__2, &a[j + a_dim1], lda, &a[pvt + a_dim1], lda);
              if (pvt < *n)
                {
                  i__2 = *n - pvt;
                  dswap_(&i__2, &a[pvt + 1 + j * a_dim1], &c__1, &a[pvt + 1
                         + pvt * a_dim1], &c__1);
                }
              i__2 = pvt - j - 1;
              dswap_(&i__2, &a[j + 1 + j * a_dim1], &c__1, &a[pvt + (j + 1)
                     * a_dim1], lda);

              /* Swap dot products and PIV */

              dtemp = work[j];
              work[j] = work[pvt];
              work[pvt] = dtemp;
              itemp = piv[pvt];
              piv[pvt] = piv[j];
              piv[j] = itemp;
            }

          ajj = sqrt(ajj);
          a[j + j * a_dim1] = ajj;

          /* Compute elements J+1:N of column J */

          if (j < *n)
            {
              i__2 = *n - j;
              i__3 = j - 1;
              dgemv_("No Trans", &i__2, &i__3, &c_b16, &a[j + 1 + a_dim1],
                     lda, &a[j + a_dim1], lda, &c_b18, &a[j + 1 + j *
                         a_dim1], &c__1);
              i__2 = *n - j;
              d__1 = 1. / ajj;
              dscal_(&i__2, &d__1, &a[j + 1 + j * a_dim1], &c__1);
            }

          /* L150: */
        }

    }

  /* Ran to completion, A has full rank */

  *rank = *n;

  goto L170;
L160:

  /* Rank is number of steps completed.  Set INFO = 1 to signal */
  /* that the factorization cannot be used to solve a system. */

  *rank = j - 1;
  *info = 1;

L170:
  return 0;

  /* End of DPSTF2 */

} /* dpstf2_ */

int C2F(dmaxloc)(double *a, int *dimm)
{
  /* System generated locals */
  int ret_val, i__1;

  /* Local variables */
  int i__;
  double dmax__;



  /* Parameter adjustments */
  --a;

  /* Function Body */
  ret_val = 1;
  dmax__ = a[1];
  i__1 = *dimm;
  for (i__ = 2; i__ <= i__1; ++i__)
    {
      if (dmax__ < a[i__])
        {
          dmax__ = a[i__];
          ret_val = i__;
        }
      /* L20: */
    }
  return ret_val;
} /* dmaxloc_ */

#endif /* PNL_HAVE_DPSTRF */
