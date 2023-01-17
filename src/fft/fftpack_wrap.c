
/************************************************************************/
/* Copyright Jérôme Lelong <jerome.lelong@gmail.com>                    */
/*                                                                      */
/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as       */
/* published by the Free Software Foundation, either version 3 of the   */
/* License, or (at your option) any later version.                      */
/*                                                                      */
/* This program is distributed in the hope that it will be useful, but  */
/* WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    */
/* Lesser General Public License for more details.                      */
/*                                                                      */
/* You should have received a copy of the GNU Lesser General Public     */
/* License along with this program.  If not, see                        */
/* <http://www.gnu.org/licenses/>.                                      */
/************************************************************************/

#include <string.h>
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_fft.h"
#include "fftpack.h"

/**
 * Allocate a workspace for FFT functions
 * 
 * @param func the function name
 * @param n the size of the input data
 * @return the allocated array or NULL
 */
double* pnl_fft_alloc_wspace(const char *func, int n)
{
  if (strcmp(func, "pnl_fft_inplace") == 0 || strcmp(func, "pnl_ifft_inplace") == 0
      || strcmp(func, "pnl_real_fft_inplace") == 0 || strcmp(func, "pnl_real_ifft_inplace") == 0
      || strcmp(func, "pnl_real_ifft") == 0)
  {
    return malloc((4 * n + 15) * sizeof(double));
  }
  else if (strcmp(func, "pnl_real_fft") == 0 )
  {
    return malloc((5 * n + 15) * sizeof(double));
  }
  else
  {
    perror("un");
  }
  return NULL;
}

/**
 * In-place Forward FFT
 *
 * @param data input complex vector. On output contains the FFT of the input vector
 * @return PNL_OK or PNL_FAIL
 */
int pnl_fft_inplace_with_wspace(PnlVectComplex *data, double *wspace)
{
  int n;
  n = data->size;

  cffti(n, wspace);
  cfftf(n, (double *)(data->array), wspace);

  return PNL_OK;
}

/**
 * In-place Inverse FFT
 *
 * @param data input complex vector. On output contains the inverse
 * FFT of the input vector
 * @return PNL_OK or PNL_FAIL
 */
int pnl_ifft_inplace_with_wspace(PnlVectComplex *data, double *wspace)
{
  int n;
  n = data->size;

  cffti(n, wspace);
  cfftb(n, (double *)(data->array), wspace);
  pnl_vect_complex_mult_double(data, 1.0 / (double) n);

  return PNL_OK;
}

/**
 *  Forward FFT of real valued sequence
 *
 * @param data an array of real numbers. On exit contains the FFT of the input sequence.
 * @param n size of re and im
 * @return PNL_OK or PNL_FAIL
 */
int pnl_real_fft_inplace_with_wspace(double *data, double *wspace, int n)
{
  rffti(n, wspace);
  rfftf(n, data, wspace);

  return PNL_OK;
}

/**
 * Inverse FFT of a sequence to be known to be the FFT a real valued sequence.
 *
 * @param data an array of real numbers. On exit contains the inverse FFT.
 * @param n size of re and im
 * @return PNL_OK or PNL_FAIL
 */
int pnl_real_ifft_inplace_with_wspace(double *data, double *wspace, int n)
{
  int i;

  rffti(n, wspace);
  rfftb(n, data, wspace);
  for (i = 0; i < n; i++)
    {
      data[i] /= (double) n;
    }
  return PNL_OK;
}

/**
 * Forward FFT of real valued sequence
 *
 * @param in an array of real numbers.
 * @param out a complex vector containing the FFT sequence
 * @return PNL_OK or PNL_FAIL
 */
int pnl_real_fft_with_wspace(const PnlVect *in, PnlVectComplex *out, double *wspace)
{
  int     n, i, l;
  double *data, *wsave;

  n = in->size;
  /* Split wspace in two parts */
  data = wspace;
  wsave = wspace + n;


  for (i = 0; i < n; i++)
    {
      data[i] = PNL_GET(in, i);
    }
  if (pnl_real_fft_inplace_with_wspace(data, wsave, n) == PNL_FAIL) return PNL_FAIL;

  pnl_vect_complex_resize(out, n);
  if (PNL_IS_ODD(n))
    {
      l = (n + 1) / 2;
    }
  else
    {
      l = n / 2;
    }

  LET_REAL(out, 0) = data[0];
  LET_IMAG(out, 0) = 0.;

  for (i = 1; i < l; i++)
    {
      LET_REAL(out, i) = data[2 * i - 1];
      LET_IMAG(out, i) = data[2 * i];
      LET_REAL(out, n - i) = data[2 * i - 1];
      LET_IMAG(out, n - i) = -data[2 * i];
    }
  if (PNL_IS_EVEN(n))
    {
      LET_REAL(out, l) = data[n - 1];
      LET_IMAG(out, l) = 0.;
    }
  return PNL_OK;
}

/**
 * Forward FFT of real valued sequence
 *
 * @param in an array of complex numbers representing the FFT of a real valued sequence
 * @param out a real vector containing the (real-) inverse FFT of in
 * @return PNL_OK or PNL_FAIL
 */
int pnl_real_ifft_with_wspace(const PnlVectComplex *in, PnlVect *out, double *wspace)
{
  int     n, i, l;

  n = in->size;
  pnl_vect_resize(out, n);

  LET(out, 0) = GET_REAL(in, 0);
  if (PNL_IS_ODD(n))
    {
      l = (n + 1) / 2;
    }
  else
    {
      l = n / 2;
    }
  for (i = 1; i < l; i++)
    {
      LET(out, 2 * i - 1) = GET_REAL(in, i);
      LET(out, 2 * i) = GET_IMAG(in, i);
    }
  if (PNL_IS_EVEN(n))
    {
      LET(out, n - 1) = GET_REAL(in, l);
    }

  if (pnl_real_ifft_inplace_with_wspace(out->array, wspace, n) == PNL_FAIL) return PNL_FAIL;

  return PNL_OK;
}


/**
 * In-place Forward FFT
 *
 * @param data input complex vector. On output contains the FFT of the input vector
 * @return PNL_OK or PNL_FAIL
 */
int pnl_fft_inplace(PnlVectComplex *data)
{
  int n = data->size;
  double *wspace = pnl_fft_alloc_wspace("pnl_fft_inplace", n);
  int res = pnl_fft_inplace_with_wspace(data, wspace);
  free(wspace);
  return res;
}

/**
 * In-place Inverse FFT
 *
 * @param data input complex vector. On output contains the inverse
 * FFT of the input vector
 * @return PNL_OK or PNL_FAIL
 */
int pnl_ifft_inplace(PnlVectComplex *data)
{
  int n = data->size;
  double *wspace = pnl_fft_alloc_wspace("pnl_ifft_inplace", n);
  int res = pnl_ifft_inplace_with_wspace(data, wspace);
  free(wspace);
  return res;
}

/**
 * Forward FFT
 *
 * @param in input complex vector
 * @param out on output contains the FFT of the input vector. This vector
 * must have already been allocated
 * @return PNL_OK or PNL_FAIL
 */
int pnl_fft(const PnlVectComplex *in, PnlVectComplex *out)
{
  pnl_vect_complex_clone(out, in);
  return pnl_fft_inplace(out);
}

/**
 * Inverse FFT
 *
 * @param in input complex vector
 * @param out on output contains the Inverse FFT of the input vector. This vector
 * must have already been allocated
 * @return PNL_OK or PNL_FAIL
 */
int pnl_ifft(const PnlVectComplex *in, PnlVectComplex *out)
{
  pnl_vect_complex_clone(out, in);
  return pnl_ifft_inplace(out);
}

/**
 * Forward FFT
 *
 * @param re real part of the data. On exit contains the real part of the FFT.
 * @param im imaginary part of the data. On exit contains the imaginary part of the FFT.
 * @param n size of re and im
 * @return PNL_OK or PNL_FAIL
 */
int pnl_fft2(double *re, double *im, int n)
{
  PnlVectComplex *in;
  int i;
  in = pnl_vect_complex_create(n);

  for (i = 0; i < n; i++)
    {
      in->array[i].r = re[i];
      in->array[i].i = im[i];
    }
  if (pnl_fft_inplace(in) == PNL_FAIL) return PNL_FAIL;

  for (i = 0; i < n; i++)
    {
      re[i] = in->array[i].r;
      im[i] = in->array[i].i;
    }
  pnl_vect_complex_free(&in);
  return PNL_OK;
}

/**
 * Inverse FFT
 *
 * @param re real part of the data. On exit contains the real part of the
 * inverse FFT.
 * @param im imaginary part of the data. On exit contains the imaginary part
 * of the inverse FFT.
 * @param n size of re and im
 * @return PNL_OK or PNL_FAIL
 */
int pnl_ifft2(double *re, double *im, int n)
{
  PnlVectComplex *in;
  int i;
  in = pnl_vect_complex_create(n);

  for (i = 0; i < n; i++)
    {
      in->array[i].r = re[i];
      in->array[i].i = im[i];
    }
  if (pnl_ifft_inplace(in) == PNL_FAIL) return PNL_FAIL;

  for (i = 0; i < n; i++)
    {
      re[i] = in->array[i].r;
      im[i] = in->array[i].i;
    }

  pnl_vect_complex_free(&in);
  return PNL_OK;
}

/**
 *  Forward FFT of real valued sequence
 *
 * @param data an array of real numbers. On exit contains the FFT of the input sequence.
 * @param n size of re and im
 * @return PNL_OK or PNL_FAIL
 */
int pnl_real_fft_inplace(double *data, int n)
{
  int res;
  double *wsave = pnl_fft_alloc_wspace("pnl_real_fft_inplace", n);
  if (!wsave) return PNL_FAIL;
  res = pnl_real_fft_inplace_with_wspace(data, wsave, n);
  free(wsave);
  return res;
}

/**
 * Inverse FFT of a sequence to be known to be the FFT a real valued sequence.
 *
 * @param data an array of real numbers. On exit contains the inverse FFT.
 * @param n size of re and im
 * @return PNL_OK or PNL_FAIL
 */
int pnl_real_ifft_inplace(double *data, int n)
{
  int res;
  double *wsave = pnl_fft_alloc_wspace("pnl_real_ifft_inplace", n);
  if (!wsave) return PNL_FAIL;
  res = pnl_real_ifft_inplace_with_wspace(data, wsave, n);
  free(wsave);
  return res;
}

/**
 * Forward FFT of real valued sequence
 *
 * @param in an array of real numbers.
 * @param out a complex vector containing the FFT sequence
 * @return PNL_OK or PNL_FAIL
 */
int pnl_real_fft(const PnlVect *in, PnlVectComplex *out)
{
  int res;
  double *wspace = pnl_fft_alloc_wspace("pnl_real_fft", in->size);
  if (!wspace) return PNL_FAIL;
  res = pnl_real_fft_with_wspace(in, out, wspace);
  free(wspace);
  return res;
}

/**
 * Forward FFT of real valued sequence
 *
 * @param in an array of complex numbers representing the FFT of a real valued sequence
 * @param out a real vector containing the (real-) inverse FFT of in
 * @return PNL_OK or PNL_FAIL
 */
int pnl_real_ifft(const PnlVectComplex *in, PnlVect *out)
{
  int res;
  double *wspace = pnl_fft_alloc_wspace("pnl_real_ifft", in->size);
  if (!wspace) return PNL_FAIL;
  res = pnl_real_ifft_with_wspace(in, out, wspace);
  free(wspace);
  return res;
}

/**
 * Forward Real FFT
 *
 * @param re real valued data data. On exit contains the real part of the FFT.
 * @param im unused on input. On exit contains the imaginary part of the FFT.
 * @param n size of re and im
 * @return PNL_OK or PNL_FAIL
 */
int pnl_real_fft2(double *re, double *im, int n)
{

  int i, l;
  if (pnl_real_fft_inplace(re, n) == PNL_FAIL) return PNL_FAIL;

  if (PNL_IS_ODD(n))
    {
      l = (n + 1) / 2;
    }
  else
    {
      l = n / 2;
    }
  im[0] = 0.;

  /* the two following loops must be decoupled, because they shuffle the
     elements around */
  for (i = 1; i < l; i++)
    {
      re[i] = re[2 * i - 1];
      im[i] = re[2 * i];
    }
  if (PNL_IS_EVEN(n))
    {
      re[l] = re[n - 1];
      im[l] = 0.;
    }
  for (i = 1; i < l; i++)
    {
      re[n - i] = re[i];
      im[n - i] = -im[i];
    }

  return PNL_OK;
}

/**
 * Inverse Real FFT
 *
 * @param re real part of the data. On exit contains the
 * inverse FFT which is real valued
 * @param im imaginary part of the data. Unusedon output
 * of the inverse FFT.
 * @param n size of re and im
 * @return PNL_OK or PNL_FAIL
 */
int pnl_real_ifft2(double *re, double *im, int n)
{
  int i, l;
  if (PNL_IS_ODD(n))
    {
      l = (n + 1) / 2;
    }
  else
    {
      l = n / 2;
    }
  if (PNL_IS_EVEN(n))
    {
      re[n - 1] = re[l];
    }
  for (i = l - 1; i >= 1; i--)
    {
      re[2 * i] = im[i];
      re[2 * i - 1] = re[i];
    }
  /* no imaginary in the data to be recovered */
  for (i = 0; i < n; i++)
    {
      im[i] = 0.;
    }

  if (pnl_real_ifft_inplace(re, n) == PNL_FAIL) return PNL_FAIL;
  return PNL_OK;
}

/**
 * Compute the 2D FFT a matrix
 *
 * @param data a complex matrix
 *
 * @return PNL_OK or PNL_FAIL
 */
int pnl_fft2d_inplace(PnlMatComplex *data)
{
  int i, j, n;
  double *wsave;
  dcomplex *row;


  n = MAX(data->m, data->n);
  if ((wsave = malloc((4 * n + 15) * sizeof(double))) == NULL) return PNL_FAIL;
  cffti(data->n, wsave);
  /*
   * First, compute the FFT of each row inplace
   */
  for (i = 0 ; i < data->m ; i++)
    {
      cfftf(data->n, (double *)(data->array + i * data->n), wsave);
    }

  /*
   * Second, compute the FFT of each column of the matrix resulting from
   * step 1
   */
  cffti(data->m, wsave);
  if ((row = malloc(sizeof(dcomplex) * data->m)) == NULL) return PNL_FAIL;
  for (j = 0 ; j < data->n ; j++)
    {
      for (i = 0 ; i < data->m ; i++) row[i] = PNL_MGET(data, i, j);
      cfftf(data->m, (double *)(row), wsave);
      for (i = 0 ; i < data->m ; i++)  PNL_MLET(data, i, j) = row[i];
    }

  free(wsave);
  free(row);
  return PNL_OK;
}

/**
 * Compute the 2D FFT a matrix
 *
 * @param in input complex matrix
 * @param out on output contains the FFT of the input matrix. This matrix
 * must have already been allocated
 * @return PNL_OK or PNL_FAIL
 */
int pnl_fft2d(const PnlMatComplex *in, PnlMatComplex *out)
{
  pnl_mat_complex_clone(out, in);
  return pnl_fft2d_inplace(out);
}

/**
 * Compute the 2D inverse (backward) FFT a matrix
 *
 * @param data a complex matrix
 *
 * @return PNL_OK or PNL_FAIL
 */
int pnl_ifft2d_inplace(PnlMatComplex *data)
{
  int i, j, n;
  double *wsave;
  dcomplex *col;


  n = MAX(data->m, data->n);
  if ((wsave = malloc((4 * n + 15) * sizeof(double))) == NULL) return PNL_FAIL;
  cffti(data->n, wsave);
  /*
   * First, compute the FFT of each col inplace
   */
  for (i = 0 ; i < data->m ; i++)
    {
      cfftb(data->n, (double *)(data->array + i * data->n), wsave);
    }

  /*
   * Second, compute the FFT of each column of the matrix resulting from
   * step 1
   */
  cffti(data->m, wsave);
  if ((col = malloc(sizeof(dcomplex) * data->m)) == NULL) return PNL_FAIL;
  for (j = 0 ; j < data->n ; j++)
    {
      /* Don't forget the renormalization from the previous step */
      for (i = 0 ; i < data->m ; i++)
        {
          col[i] = CRdiv(PNL_MGET(data, i, j), (double) data->m);
        }
      cfftb(data->m, (double *)(col), wsave);
      for (i = 0 ; i < data->m ; i++)
        {
          PNL_MLET(data, i, j) = CRdiv(col[i], (double) data->n);
        }
    }

  free(wsave);
  free(col);
  return PNL_OK;
}

/**
 * Compute the inverse 2D FFT a matrix
 *
 * @param in input complex matrix
 * @param out on output contains the FFT of the input matrix. This matrix
 * must have already been allocated
 * @return PNL_OK or PNL_FAIL
 */
int pnl_ifft2d(const PnlMatComplex *in, PnlMatComplex *out)
{
  pnl_mat_complex_clone(out, in);
  return pnl_ifft2d_inplace(out);
}

/**
 * Compute the 2D FFT a real matrix
 *
 * @param in input real matrix
 * @param out on output contains the FFT of the input matrix. This matrix
 * must have already been allocated
 * @return PNL_OK or PNL_FAIL
 */
int pnl_real_fft2d(const PnlMat *in, PnlMatComplex *out)
{

  int i, j, n , l;
  double *data, *wsave;
  dcomplex *col;
  pnl_mat_complex_resize(out, in->m, in->n);

  n = MAX(in->m, in->n);
  if ((wsave = malloc((4 * n + 15) * sizeof(double))) == NULL) return PNL_FAIL;
  if ((data = malloc(n * sizeof(double))) == NULL) return PNL_FAIL;

  /*
   * Compute the FFT of each row
   */
  rffti(in->n, wsave);
  if (PNL_IS_ODD(in->n))
    {
      l = (in->n + 1) / 2;
    }
  else
    {
      l = in->n / 2;
    }

  for (i = 0 ; i < in->m ; i++)
    {
      memcpy(data, in->array + i * in->n, in->n * sizeof(double));
      rfftf(in->n, data, wsave);

      /*
       * Copy output data
       */

      MLET_REAL(out, i, 0) = data[0];
      MLET_IMAG(out, i, 0) = 0.;

      for (j = 1; j < l; j++)
        {
          MLET_REAL(out, i, j) = data[2 * j - 1];
          MLET_IMAG(out, i, j) = data[2 * j];
          MLET_REAL(out, i, n - j) = data[2 * j - 1];
          MLET_IMAG(out, i, n - j) = -data[2 * j];
        }
      if (PNL_IS_EVEN(n))
        {
          MLET_REAL(out, i, l) = data[n - 1];
          MLET_IMAG(out, i, l) = 0.;
        }
    }

  /*
   * Second, compute the FFT of each column of the matrix resulting from
   * step 1
   */
  cffti(in->m, wsave);
  if ((col = malloc(sizeof(dcomplex) * in->m)) == NULL) return PNL_FAIL;
  for (j = 0 ; j < in->n ; j++)
    {
      for (i = 0 ; i < in->m ; i++) col[i] = PNL_MGET(out, i, j);
      cfftf(out->m, (double *)(col), wsave);
      for (i = 0 ; i < in->m ; i++)  PNL_MLET(out, i, j) = col[i];
    }

  free(wsave);
  free(col);
  free(data);
  return PNL_OK;
}

/**
 * Compute the inverse 2D FFT a complex matrix known to be the 2D FFT of a
 * real matrix
 *
 * @param in input complex matrix. This matrix is lost on output
 * @param out real matrix on output contains the 2D FFT of the input matrix. This matrix
 * must have already been allocated
 * @return PNL_OK or PNL_FAIL
 */
int pnl_real_ifft2d(const PnlMatComplex *in, PnlMat *out)
{

  int i, j, n , l;
  double *wsave, *col;
  PnlMatComplex *data;
  pnl_mat_resize(out, in->m, in->n);

  n = MAX(in->m, in->n);
  if ((wsave = malloc((4 * n + 15) * sizeof(double))) == NULL) return PNL_FAIL;

  /*
   * Compute the FFT of each row
   */
  data = pnl_mat_complex_copy(in);
  cffti(in->n, wsave);

  for (i = 0 ; i < in->m ; i++)
    {
      cfftb(in->n, (double *)(data->array + i * in->n), wsave);
    }

  /*
   * Second, compute the FFT of each column of the matrix resulting from
   * step 1
   */
  rffti(in->m, wsave);
  if (PNL_IS_ODD(in->m))
    {
      l = (in->m + 1) / 2;
    }
  else
    {
      l = in->m / 2;
    }
  if ((col = malloc(sizeof(double) * in->m)) == NULL) return PNL_FAIL;
  for (j = 0 ; j < in->n ; j++)
    {
      /* Don't forget the renormalization from the previous step */
      col[0] = MGET_REAL(data, 0, j) / in->n;
      for (i = 1; i < l; i++)
        {
          col[2 * i - 1] = MGET_REAL(data, i, j) / in->n;
          col[2 * i] = MGET_IMAG(data, i, j) / in->n;
        }
      if (PNL_IS_EVEN(in->m))
        {
          col[in->m - 1] = MGET_REAL(data, l, j) / in->n;
        }
      rfftb(out->m, col, wsave);
      for (i = 0 ; i < in->m ; i++)  PNL_MLET(out, i, j) = col[i] / in->m;
    }

  free(col);
  free(wsave);
  pnl_mat_complex_free(&data);
  return PNL_OK;

}

