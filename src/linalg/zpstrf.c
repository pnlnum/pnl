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

static dcomplex c_b1 = {1., 0.};
static int c__1 = 1;
static int c_n1 = -1;
static double c_b29 = -1.;
static double c_b30 = 1.;

static void d_cnjg(dcomplex *r, dcomplex *z)
{
  double zi = z->i;
  r->r = z->r;
  r->i = -zi;
}

int zpstf2_(char *uplo, int *n, dcomplex *a,
            int *lda, int *piv, int *rank, double *tol,
            double *work, int *info)
{
  /* System generated locals */
  int a_dim1, a_offset, i__1, i__2, i__3;
  double d__1;
  dcomplex z__1, z__2;

  /* Builtin functions */
  void d_cnjg(dcomplex *, dcomplex *);
  double sqrt(double);

  /* Local variables */
  int i__, j, maxlocval;
  double ajj;
  int pvt;
  extern int lsame_(char *, char *);
  double dtemp;
  int itemp;
  extern  int zgemv_(char *, int *, int *,
                     dcomplex *, dcomplex *, int *, dcomplex *,
                     int *, dcomplex *, dcomplex *, int *);
  double dstop;
  int upper;
  dcomplex ztemp;
  extern  int zswap_(int *, dcomplex *, int *,
                     dcomplex *, int *);
  extern double dlamch_(char *);
  extern int disnan_(double *);
  extern  int xerbla_(char *, int *), zdscal_(
    int *, double *, dcomplex *, int *), zlacgv_(
      int *, dcomplex *, int *);
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

  /*  ZPSTF2 computes the Cholesky factorization with complete */
  /*  pivoting of a complex Hermitian positive semidefinite matrix A. */

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

  /*  A       (input/output) COMPLEX*16 array, dimension (LDA,N) */
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
      xerbla_("ZPSTF2", &i__1);
      return 0;
    }

  /*     Quick return if possible */

  if (*n == 0)
    {
      return 0;
    }

  /*     Initialize PIV */

  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      piv[i__] = i__;
      /* L100: */
    }

  /*     Compute stopping value */

  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      i__2 = i__ + i__ * a_dim1;
      work[i__] = a[i__2].r;
      /* L110: */
    }
  pvt = dmaxloc_(&work[1], n);
  i__1 = pvt + pvt * a_dim1;
  ajj = a[i__1].r;
  if (ajj == 0. || disnan_(&ajj))
    {
      *rank = 0;
      *info = 1;
      goto L200;
    }

  /*     Compute stopping value if not supplied */

  if (*tol < 0.)
    {
      dstop = *n * dlamch_("Epsilon") * ajj;
    }
  else
    {
      dstop = *tol;
    }

  /*     Set first half of WORK to zero, holds dot products */

  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      work[i__] = 0.;
      /* L120: */
    }

  if (upper)
    {

      /*        Compute the Cholesky factorization P' * A * P = U' * U */

      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
        {

          /*        Find pivot, test for exit, else swap rows and columns */
          /*        Update dot products, compute possible pivots which are */
          /*        stored in the second half of WORK */

          i__2 = *n;
          for (i__ = j; i__ <= i__2; ++i__)
            {

              if (j > 1)
                {
                  d_cnjg(&z__2, &a[j - 1 + i__ * a_dim1]);
                  i__3 = j - 1 + i__ * a_dim1;
                  z__1.r = z__2.r * a[i__3].r - z__2.i * a[i__3].i, z__1.i =
                             z__2.r * a[i__3].i + z__2.i * a[i__3].r;
                  work[i__] += z__1.r;
                }
              i__3 = i__ + i__ * a_dim1;
              work[*n + i__] = a[i__3].r - work[i__];

              /* L130: */
            }

          if (j > 1)
            {
              maxlocval = (*n << 1) - (*n + j) + 1;
              itemp = dmaxloc_(&work[*n + j], &maxlocval);
              pvt = itemp + j - 1;
              ajj = work[*n + pvt];
              if (ajj <= dstop || disnan_(&ajj))
                {
                  i__2 = j + j * a_dim1;
                  a[i__2].r = ajj, a[i__2].i = 0.;
                  goto L190;
                }
            }

          if (j != pvt)
            {

              /*              Pivot PNL_OK, so can now swap pivot rows and columns */

              i__2 = pvt + pvt * a_dim1;
              i__3 = j + j * a_dim1;
              a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
              i__2 = j - 1;
              zswap_(&i__2, &a[j * a_dim1 + 1], &c__1, &a[pvt * a_dim1 + 1],
                     &c__1);
              if (pvt < *n)
                {
                  i__2 = *n - pvt;
                  zswap_(&i__2, &a[j + (pvt + 1) * a_dim1], lda, &a[pvt + (
                           pvt + 1) * a_dim1], lda);
                }
              i__2 = pvt - 1;
              for (i__ = j + 1; i__ <= i__2; ++i__)
                {
                  d_cnjg(&z__1, &a[j + i__ * a_dim1]);
                  ztemp.r = z__1.r, ztemp.i = z__1.i;
                  i__3 = j + i__ * a_dim1;
                  d_cnjg(&z__1, &a[i__ + pvt * a_dim1]);
                  a[i__3].r = z__1.r, a[i__3].i = z__1.i;
                  i__3 = i__ + pvt * a_dim1;
                  a[i__3].r = ztemp.r, a[i__3].i = ztemp.i;
                  /* L140: */
                }
              i__2 = j + pvt * a_dim1;
              d_cnjg(&z__1, &a[j + pvt * a_dim1]);
              a[i__2].r = z__1.r, a[i__2].i = z__1.i;

              /*              Swap dot products and PIV */

              dtemp = work[j];
              work[j] = work[pvt];
              work[pvt] = dtemp;
              itemp = piv[pvt];
              piv[pvt] = piv[j];
              piv[j] = itemp;
            }

          ajj = sqrt(ajj);
          i__2 = j + j * a_dim1;
          a[i__2].r = ajj, a[i__2].i = 0.;

          /*           Compute elements J+1:N of row J */

          if (j < *n)
            {
              i__2 = j - 1;
              zlacgv_(&i__2, &a[j * a_dim1 + 1], &c__1);
              i__2 = j - 1;
              i__3 = *n - j;
              z__1.r = -1., z__1.i = -0.;
              zgemv_("Trans", &i__2, &i__3, &z__1, &a[(j + 1) * a_dim1 + 1],
                     lda, &a[j * a_dim1 + 1], &c__1, &c_b1, &a[j + (j + 1)
                         * a_dim1], lda);
              i__2 = j - 1;
              zlacgv_(&i__2, &a[j * a_dim1 + 1], &c__1);
              i__2 = *n - j;
              d__1 = 1. / ajj;
              zdscal_(&i__2, &d__1, &a[j + (j + 1) * a_dim1], lda);
            }

          /* L150: */
        }

    }
  else
    {

      /*        Compute the Cholesky factorization P' * A * P = L * L' */

      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
        {

          /*        Find pivot, test for exit, else swap rows and columns */
          /*        Update dot products, compute possible pivots which are */
          /*        stored in the second half of WORK */

          i__2 = *n;
          for (i__ = j; i__ <= i__2; ++i__)
            {

              if (j > 1)
                {
                  d_cnjg(&z__2, &a[i__ + (j - 1) * a_dim1]);
                  i__3 = i__ + (j - 1) * a_dim1;
                  z__1.r = z__2.r * a[i__3].r - z__2.i * a[i__3].i, z__1.i =
                             z__2.r * a[i__3].i + z__2.i * a[i__3].r;
                  work[i__] += z__1.r;
                }
              i__3 = i__ + i__ * a_dim1;
              work[*n + i__] = a[i__3].r - work[i__];

              /* L160: */
            }

          if (j > 1)
            {
              maxlocval = (*n << 1) - (*n + j) + 1;
              itemp = dmaxloc_(&work[*n + j], &maxlocval);
              pvt = itemp + j - 1;
              ajj = work[*n + pvt];
              if (ajj <= dstop || disnan_(&ajj))
                {
                  i__2 = j + j * a_dim1;
                  a[i__2].r = ajj, a[i__2].i = 0.;
                  goto L190;
                }
            }

          if (j != pvt)
            {

              /*              Pivot PNL_OK, so can now swap pivot rows and columns */

              i__2 = pvt + pvt * a_dim1;
              i__3 = j + j * a_dim1;
              a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
              i__2 = j - 1;
              zswap_(&i__2, &a[j + a_dim1], lda, &a[pvt + a_dim1], lda);
              if (pvt < *n)
                {
                  i__2 = *n - pvt;
                  zswap_(&i__2, &a[pvt + 1 + j * a_dim1], &c__1, &a[pvt + 1
                         + pvt * a_dim1], &c__1);
                }
              i__2 = pvt - 1;
              for (i__ = j + 1; i__ <= i__2; ++i__)
                {
                  d_cnjg(&z__1, &a[i__ + j * a_dim1]);
                  ztemp.r = z__1.r, ztemp.i = z__1.i;
                  i__3 = i__ + j * a_dim1;
                  d_cnjg(&z__1, &a[pvt + i__ * a_dim1]);
                  a[i__3].r = z__1.r, a[i__3].i = z__1.i;
                  i__3 = pvt + i__ * a_dim1;
                  a[i__3].r = ztemp.r, a[i__3].i = ztemp.i;
                  /* L170: */
                }
              i__2 = pvt + j * a_dim1;
              d_cnjg(&z__1, &a[pvt + j * a_dim1]);
              a[i__2].r = z__1.r, a[i__2].i = z__1.i;

              /*              Swap dot products and PIV */

              dtemp = work[j];
              work[j] = work[pvt];
              work[pvt] = dtemp;
              itemp = piv[pvt];
              piv[pvt] = piv[j];
              piv[j] = itemp;
            }

          ajj = sqrt(ajj);
          i__2 = j + j * a_dim1;
          a[i__2].r = ajj, a[i__2].i = 0.;

          /*           Compute elements J+1:N of column J */

          if (j < *n)
            {
              i__2 = j - 1;
              zlacgv_(&i__2, &a[j + a_dim1], lda);
              i__2 = *n - j;
              i__3 = j - 1;
              z__1.r = -1., z__1.i = -0.;
              zgemv_("No Trans", &i__2, &i__3, &z__1, &a[j + 1 + a_dim1],
                     lda, &a[j + a_dim1], lda, &c_b1, &a[j + 1 + j *
                         a_dim1], &c__1);
              i__2 = j - 1;
              zlacgv_(&i__2, &a[j + a_dim1], lda);
              i__2 = *n - j;
              d__1 = 1. / ajj;
              zdscal_(&i__2, &d__1, &a[j + 1 + j * a_dim1], &c__1);
            }

          /* L180: */
        }

    }

  /*     Ran to completion, A has full rank */

  *rank = *n;

  goto L200;
L190:

  /*     Rank is number of steps completed.  Set INFO = 1 to signal */
  /*     that the factorization cannot be used to solve a system. */

  *rank = j - 1;
  *info = 1;

L200:
  return 0;

  /*     End of ZPSTF2 */

} /* zpstf2_ */


int zpstrf_(char *uplo, int *n, dcomplex *a,
            int *lda, int *piv, int *rank, double *tol,
            double *work, int *info)
{
  /* System generated locals */
  int a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
  double d__1;
  dcomplex z__1, z__2;

  /* Builtin functions */
  void d_cnjg(dcomplex *, dcomplex *);
  double sqrt(double);

  /* Local variables */
  int i__, j, k, maxlocval, jb, nb;
  double ajj;
  int pvt;
  extern int lsame_(char *, char *);
  double dtemp;
  int itemp;
  extern  int zherk_(char *, char *, int *, int *,
                     double *, dcomplex *, int *, double *,
                     dcomplex *, int *), zgemv_(char *,
                         int *, int *, dcomplex *, dcomplex *, int *,
                         dcomplex *, int *, dcomplex *, dcomplex *,
                         int *);
  double dstop;
  int upper;
  dcomplex ztemp;
  extern  int zswap_(int *, dcomplex *, int *,
                     dcomplex *, int *);
  extern double dlamch_(char *);
  extern  int zpstf2_(char *, int *, dcomplex *,
                      int *, int *, int *, double *, double *,
                      int *);
  extern int disnan_(double *);
  extern  int xerbla_(char *, int *);
  extern int ilaenv_(int *, char *, char *, int *, int *,
                     int *, int *);
  extern  int zdscal_(int *, double *,
                      dcomplex *, int *), zlacgv_(int *, dcomplex *,
                          int *);
  extern int dmaxloc_(double *, int *);


  /*  -- LAPACK routine (version 3.2.1)                                  -- */

  /*  -- Contributed by Craig Lucas, University of Manchester / NAG Ltd. -- */
  /*  -- April 2009                                                      -- */

  /*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  ZPSTRF computes the Cholesky factorization with complete */
  /*  pivoting of a complex Hermitian positive semidefinite matrix A. */

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

  /*  A       (input/output) COMPLEX*16 array, dimension (LDA,N) */
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
      xerbla_("ZPSTRF", &i__1);
      return 0;
    }

  /*     Quick return if possible */

  if (*n == 0)
    {
      return 0;
    }

  /*     Get block size */

  nb = ilaenv_(&c__1, "ZPOTRF", uplo, n, &c_n1, &c_n1, &c_n1);
  if (nb <= 1 || nb >= *n)
    {

      /*        Use unblocked code */

      zpstf2_(uplo, n, &a[a_dim1 + 1], lda, &piv[1], rank, tol, &work[1],
              info);
      goto L230;

    }
  else
    {

      /*     Initialize PIV */

      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__)
        {
          piv[i__] = i__;
          /* L100: */
        }

      /*     Compute stopping value */

      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__)
        {
          i__2 = i__ + i__ * a_dim1;
          work[i__] = a[i__2].r;
          /* L110: */
        }
      pvt = dmaxloc_(&work[1], n);
      i__1 = pvt + pvt * a_dim1;
      ajj = a[i__1].r;
      if (ajj == 0. || disnan_(&ajj))
        {
          *rank = 0;
          *info = 1;
          goto L230;
        }

      /*     Compute stopping value if not supplied */

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

          /*           Compute the Cholesky factorization P' * A * P = U' * U */

          i__1 = *n;
          i__2 = nb;
          for (k = 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2)
            {

              /*              Account for last block not being NB wide */

              /* Computing MIN */
              i__3 = nb, i__4 = *n - k + 1;
              jb = MIN(i__3, i__4);

              /*              Set relevant part of first half of WORK to zero, */
              /*              holds dot products */

              i__3 = *n;
              for (i__ = k; i__ <= i__3; ++i__)
                {
                  work[i__] = 0.;
                  /* L120: */
                }

              i__3 = k + jb - 1;
              for (j = k; j <= i__3; ++j)
                {

                  /*              Find pivot, test for exit, else swap rows and columns */
                  /*              Update dot products, compute possible pivots which are */
                  /*              stored in the second half of WORK */

                  i__4 = *n;
                  for (i__ = j; i__ <= i__4; ++i__)
                    {

                      if (j > k)
                        {
                          d_cnjg(&z__2, &a[j - 1 + i__ * a_dim1]);
                          i__5 = j - 1 + i__ * a_dim1;
                          z__1.r = z__2.r * a[i__5].r - z__2.i * a[i__5].i,
                               z__1.i = z__2.r * a[i__5].i + z__2.i * a[
                                          i__5].r;
                          work[i__] += z__1.r;
                        }
                      i__5 = i__ + i__ * a_dim1;
                      work[*n + i__] = a[i__5].r - work[i__];

                      /* L130: */
                    }

                  if (j > 1)
                    {
                      maxlocval = (*n << 1) - (*n + j) + 1;
                      itemp = dmaxloc_(&work[*n + j], &maxlocval);
                      pvt = itemp + j - 1;
                      ajj = work[*n + pvt];
                      if (ajj <= dstop || disnan_(&ajj))
                        {
                          i__4 = j + j * a_dim1;
                          a[i__4].r = ajj, a[i__4].i = 0.;
                          goto L220;
                        }
                    }

                  if (j != pvt)
                    {

                      /*                    Pivot PNL_OK, so can now swap pivot rows and columns */

                      i__4 = pvt + pvt * a_dim1;
                      i__5 = j + j * a_dim1;
                      a[i__4].r = a[i__5].r, a[i__4].i = a[i__5].i;
                      i__4 = j - 1;
                      zswap_(&i__4, &a[j * a_dim1 + 1], &c__1, &a[pvt *
                             a_dim1 + 1], &c__1);
                      if (pvt < *n)
                        {
                          i__4 = *n - pvt;
                          zswap_(&i__4, &a[j + (pvt + 1) * a_dim1], lda, &a[
                                   pvt + (pvt + 1) * a_dim1], lda);
                        }
                      i__4 = pvt - 1;
                      for (i__ = j + 1; i__ <= i__4; ++i__)
                        {
                          d_cnjg(&z__1, &a[j + i__ * a_dim1]);
                          ztemp.r = z__1.r, ztemp.i = z__1.i;
                          i__5 = j + i__ * a_dim1;
                          d_cnjg(&z__1, &a[i__ + pvt * a_dim1]);
                          a[i__5].r = z__1.r, a[i__5].i = z__1.i;
                          i__5 = i__ + pvt * a_dim1;
                          a[i__5].r = ztemp.r, a[i__5].i = ztemp.i;
                          /* L140: */
                        }
                      i__4 = j + pvt * a_dim1;
                      d_cnjg(&z__1, &a[j + pvt * a_dim1]);
                      a[i__4].r = z__1.r, a[i__4].i = z__1.i;

                      /*                    Swap dot products and PIV */

                      dtemp = work[j];
                      work[j] = work[pvt];
                      work[pvt] = dtemp;
                      itemp = piv[pvt];
                      piv[pvt] = piv[j];
                      piv[j] = itemp;
                    }

                  ajj = sqrt(ajj);
                  i__4 = j + j * a_dim1;
                  a[i__4].r = ajj, a[i__4].i = 0.;

                  /*                 Compute elements J+1:N of row J. */

                  if (j < *n)
                    {
                      i__4 = j - 1;
                      zlacgv_(&i__4, &a[j * a_dim1 + 1], &c__1);
                      i__4 = j - k;
                      i__5 = *n - j;
                      z__1.r = -1., z__1.i = -0.;
                      zgemv_("Trans", &i__4, &i__5, &z__1, &a[k + (j + 1) *
                                                              a_dim1], lda, &a[k + j * a_dim1], &c__1, &
                             c_b1, &a[j + (j + 1) * a_dim1], lda);
                      i__4 = j - 1;
                      zlacgv_(&i__4, &a[j * a_dim1 + 1], &c__1);
                      i__4 = *n - j;
                      d__1 = 1. / ajj;
                      zdscal_(&i__4, &d__1, &a[j + (j + 1) * a_dim1], lda);
                    }

                  /* L150: */
                }

              /*              Update trailing matrix, J already incremented */

              if (k + jb <= *n)
                {
                  i__3 = *n - j + 1;
                  zherk_("Upper", "Conj Trans", &i__3, &jb, &c_b29, &a[k +
                         j * a_dim1], lda, &c_b30, &a[j + j * a_dim1], lda);
                }

              /* L160: */
            }

        }
      else
        {

          /*        Compute the Cholesky factorization P' * A * P = L * L' */

          i__2 = *n;
          i__1 = nb;
          for (k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1)
            {

              /*              Account for last block not being NB wide */

              /* Computing MIN */
              i__3 = nb, i__4 = *n - k + 1;
              jb = MIN(i__3, i__4);

              /*              Set relevant part of first half of WORK to zero, */
              /*              holds dot products */

              i__3 = *n;
              for (i__ = k; i__ <= i__3; ++i__)
                {
                  work[i__] = 0.;
                  /* L170: */
                }

              i__3 = k + jb - 1;
              for (j = k; j <= i__3; ++j)
                {

                  /*              Find pivot, test for exit, else swap rows and columns */
                  /*              Update dot products, compute possible pivots which are */
                  /*              stored in the second half of WORK */

                  i__4 = *n;
                  for (i__ = j; i__ <= i__4; ++i__)
                    {

                      if (j > k)
                        {
                          d_cnjg(&z__2, &a[i__ + (j - 1) * a_dim1]);
                          i__5 = i__ + (j - 1) * a_dim1;
                          z__1.r = z__2.r * a[i__5].r - z__2.i * a[i__5].i,
                               z__1.i = z__2.r * a[i__5].i + z__2.i * a[
                                          i__5].r;
                          work[i__] += z__1.r;
                        }
                      i__5 = i__ + i__ * a_dim1;
                      work[*n + i__] = a[i__5].r - work[i__];

                      /* L180: */
                    }

                  if (j > 1)
                    {
                      maxlocval = (*n << 1) - (*n + j) + 1;
                      itemp = dmaxloc_(&work[*n + j], &maxlocval);
                      pvt = itemp + j - 1;
                      ajj = work[*n + pvt];
                      if (ajj <= dstop || disnan_(&ajj))
                        {
                          i__4 = j + j * a_dim1;
                          a[i__4].r = ajj, a[i__4].i = 0.;
                          goto L220;
                        }
                    }

                  if (j != pvt)
                    {

                      /*                    Pivot PNL_OK, so can now swap pivot rows and columns */

                      i__4 = pvt + pvt * a_dim1;
                      i__5 = j + j * a_dim1;
                      a[i__4].r = a[i__5].r, a[i__4].i = a[i__5].i;
                      i__4 = j - 1;
                      zswap_(&i__4, &a[j + a_dim1], lda, &a[pvt + a_dim1],
                             lda);
                      if (pvt < *n)
                        {
                          i__4 = *n - pvt;
                          zswap_(&i__4, &a[pvt + 1 + j * a_dim1], &c__1, &a[
                                   pvt + 1 + pvt * a_dim1], &c__1);
                        }
                      i__4 = pvt - 1;
                      for (i__ = j + 1; i__ <= i__4; ++i__)
                        {
                          d_cnjg(&z__1, &a[i__ + j * a_dim1]);
                          ztemp.r = z__1.r, ztemp.i = z__1.i;
                          i__5 = i__ + j * a_dim1;
                          d_cnjg(&z__1, &a[pvt + i__ * a_dim1]);
                          a[i__5].r = z__1.r, a[i__5].i = z__1.i;
                          i__5 = pvt + i__ * a_dim1;
                          a[i__5].r = ztemp.r, a[i__5].i = ztemp.i;
                          /* L190: */
                        }
                      i__4 = pvt + j * a_dim1;
                      d_cnjg(&z__1, &a[pvt + j * a_dim1]);
                      a[i__4].r = z__1.r, a[i__4].i = z__1.i;


                      /*                    Swap dot products and PIV */

                      dtemp = work[j];
                      work[j] = work[pvt];
                      work[pvt] = dtemp;
                      itemp = piv[pvt];
                      piv[pvt] = piv[j];
                      piv[j] = itemp;
                    }

                  ajj = sqrt(ajj);
                  i__4 = j + j * a_dim1;
                  a[i__4].r = ajj, a[i__4].i = 0.;

                  /*                 Compute elements J+1:N of column J. */

                  if (j < *n)
                    {
                      i__4 = j - 1;
                      zlacgv_(&i__4, &a[j + a_dim1], lda);
                      i__4 = *n - j;
                      i__5 = j - k;
                      z__1.r = -1., z__1.i = -0.;
                      zgemv_("No Trans", &i__4, &i__5, &z__1, &a[j + 1 + k *
                             a_dim1], lda, &a[j + k * a_dim1], lda, &c_b1,
                             &a[j + 1 + j * a_dim1], &c__1);
                      i__4 = j - 1;
                      zlacgv_(&i__4, &a[j + a_dim1], lda);
                      i__4 = *n - j;
                      d__1 = 1. / ajj;
                      zdscal_(&i__4, &d__1, &a[j + 1 + j * a_dim1], &c__1);
                    }

                  /* L200: */
                }

              /*              Update trailing matrix, J already incremented */

              if (k + jb <= *n)
                {
                  i__3 = *n - j + 1;
                  zherk_("Lower", "No Trans", &i__3, &jb, &c_b29, &a[j + k *
                         a_dim1], lda, &c_b30, &a[j + j * a_dim1], lda);
                }

              /* L210: */
            }

        }
    }

  /*     Ran to completion, A has full rank */

  *rank = *n;

  goto L230;
L220:

  /*     Rank is the number of steps completed.  Set INFO = 1 to signal */
  /*     that the factorization cannot be used to solve a system. */

  *rank = j - 1;
  *info = 1;

L230:
  return 0;

  /*     End of ZPSTRF */

} /* zpstrf_ */

#endif
