
/************************************************************************/
/* Copyright Jérôme Lelong <jerome.lelong@gmail.com>                    */
/* Copyright David Pommier <david.pommier@gmail.com>                    */
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_complex.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_fft.h"
#include "tests_utils.h"


static void fftpack(void )
{  
  double abserr=2E-4;
  int kmax=2048;
  int Nmax=2*kmax;
  double L=4.0;
  double Delta=L/(double) kmax;
  PnlVect *error;
  PnlVectComplex * vv,*vv2, *vxi,*vxi2,*result,*result2,*inverse_fft,*inverse_fft2;
  int j,i;
  error=pnl_vect_create(Nmax);
  vv=pnl_vect_complex_create(Nmax);
  vv2=pnl_vect_complex_create(Nmax);
  vxi=pnl_vect_complex_create(Nmax);
  vxi2  = pnl_vect_complex_create(Nmax);
  result=pnl_vect_complex_create(Nmax);
  result2=pnl_vect_complex_create(Nmax);
  inverse_fft=pnl_vect_complex_create(Nmax);
  inverse_fft2=pnl_vect_complex_create(Nmax);
  /*function input for direct fourier transform */
  for ( j=0 ; j<Nmax; j++ )
    {
      double yj = (j-kmax)*Delta;
      pnl_vect_complex_set(vv,j,Complex(exp(-yj*yj/2.),0.0));
      pnl_vect_complex_set(vv2,j,Complex(-1.0*yj*Creal(pnl_vect_complex_get(vv,j)),0.0));
    }

  for( j=0;j<Nmax;j++)
    {
      double xij = 2.*M_PI*(j%kmax-kmax*(j/kmax))/(Delta*Nmax);
      pnl_vect_complex_set(vxi,j,Complex(exp(-xij*xij/2.),0.0));
      pnl_vect_complex_set(vxi2,j,Complex(0.0,xij*exp(-xij*xij/2.)));
    }
  pnl_fft(vv2,result2); 
  pnl_fft(vv,result);
  /*direct fourier tranform */
  for(i=0; i< Nmax; i++)
    {
      LET(error,i)=Cabs(Csub(RCmul(Delta*M_1_SQRT2PI*PNL_ALTERNATE(i),
                                   pnl_vect_complex_get(result,i))
                             ,pnl_vect_complex_get(vxi,i)));

    }
  pnl_test_eq_abs (pnl_vect_max(error), 0., abserr, "fft (Gauss)", ""); 
  for(i=0; i< Nmax; i++)
    {
      LET(error,i)=Cabs(Csub(RCmul(Delta*M_1_SQRT2PI*PNL_ALTERNATE(i),
                                   pnl_vect_complex_get(result2,i)),
                             pnl_vect_complex_get(vxi2,i)));

    }
  pnl_test_eq_abs (pnl_vect_max(error), 0., abserr, "fft (Gauss derivative)", ""); 



  pnl_ifft(result,inverse_fft); 
  pnl_ifft(result2,inverse_fft2); 
  for(i=0; i< Nmax; i++)
    {
      LET(error,i)=Cabs(Csub(pnl_vect_complex_get(inverse_fft,i),pnl_vect_complex_get(vv,i)));
    }
  pnl_test_eq_abs (pnl_vect_max(error), 0., abserr, "ifft (Gauss)", ""); 
  for(i=0; i< Nmax; i++)
    {
      LET(error,i)=Cabs(Csub(pnl_vect_complex_get(inverse_fft2,i),pnl_vect_complex_get(vv2,i)));
    }
  pnl_test_eq_abs (pnl_vect_max(error), 0., abserr, "ifft (Gauss derivative)", ""); 

  pnl_vect_complex_free(&inverse_fft);
  pnl_vect_complex_free(&inverse_fft2);
  pnl_vect_complex_free(&result);
  pnl_vect_complex_free(&result2);
  pnl_vect_complex_free(&vv);
  pnl_vect_complex_free(&vv2);
  pnl_vect_complex_free(&vxi);
  pnl_vect_complex_free(&vxi2);
  pnl_vect_free(&error);

}

static void pnl_fft_inplace_test(void)
{  

  int               kmax  = 32;
  int               Nmax  = 2*kmax;
  double            L     = 4.0;
  double            abserr= 1E-4;
  double            Delta = L/(double) kmax;
  PnlVect          *error, *result_r, *result_i, *inverse_fft_r, *inverse_fft_i;
  PnlVectComplex   *vv, *result,*inverse_fft;
  double           *re, *im;
  PnlVect           v_re, v_im;
  int               j,i;
  error=pnl_vect_create(Nmax);
  vv=pnl_vect_complex_create(Nmax);
  result=pnl_vect_complex_create(Nmax);
  inverse_fft=pnl_vect_complex_create(Nmax);
  re = malloc (sizeof(double) * Nmax);
  im = malloc (sizeof(double) * Nmax);
  /*function input for direct fourier transform */

  for (j=0; j<Nmax; j++)
    {
      double yj = (j-kmax)*Delta;
      pnl_vect_complex_set(vv,j,Complex(cos(yj), sin(yj)));
      re[j] = GET_REAL(vv, j);
      im[j] = GET_IMAG(vv, j);
    }
  pnl_vect_complex_clone (result, vv);
  pnl_fft_inplace(result);
  pnl_fft2(re, im , Nmax);
  result_r = pnl_vect_new ();
  result_i = pnl_vect_new ();
  pnl_vect_complex_split_in_vect (result, result_r, result_i);
  v_re = pnl_vect_wrap_array (re, Nmax);
  v_im = pnl_vect_wrap_array (im, Nmax);
  pnl_test_vect_eq_abs ( result_r, &v_re, abserr, "Comparison fft and fft2 (real part)", "");
  pnl_test_vect_eq_abs ( result_i, &v_im, abserr, "Comparison fft and fft2 (imag part)", "");

  pnl_vect_complex_clone (inverse_fft, result);
  pnl_ifft_inplace(inverse_fft);
  pnl_ifft2(re, im , Nmax);

  inverse_fft_r = pnl_vect_new ();
  inverse_fft_i = pnl_vect_new ();
  pnl_vect_complex_split_in_vect (inverse_fft, inverse_fft_r, inverse_fft_i); 
  for(i=0; i< Nmax; i++)
    {
      LET(error,i)=Cabs(Csub(pnl_vect_complex_get(inverse_fft,i),pnl_vect_complex_get(vv,i)));
    }
  pnl_test_eq_abs (pnl_vect_max(error), 0., abserr, "ifft (error)", "");

  pnl_vect_complex_free(&inverse_fft);
  pnl_vect_free(&inverse_fft_r);
  pnl_vect_free(&inverse_fft_i);
  pnl_vect_complex_free(&result);
  pnl_vect_free(&result_r);
  pnl_vect_free(&result_i);
  pnl_vect_complex_free(&vv);
  pnl_vect_free(&error);
  free (re); free (im);
}

static void pnl_fft_real_test(void )
{  
  int               kmax  = 32;
  int               Nmax  = 2*kmax;
  double            L     = 4.0;
  double            Delta = L/(double) kmax;
  PnlVect          *error, *vv, *inverse_fft;
  PnlVectComplex   *fft;
  int               j;


  error       = pnl_vect_create(Nmax);
  vv          = pnl_vect_create(Nmax);
  fft         = pnl_vect_complex_create(Nmax);
  inverse_fft = pnl_vect_create(Nmax);
  /*function input for direct fourier transform */

  for (j=0; j<Nmax; j++)
    {
      double yj = (j-kmax)*Delta;
      pnl_vect_set(vv,j,cos(yj)*sin(yj));
    }

  pnl_real_fft(vv, fft);

  pnl_real_ifft(fft, inverse_fft);
  pnl_test_vect_eq_abs (inverse_fft, vv, 1E-10, "Real ifft", "");

  pnl_vect_free(&inverse_fft);
  pnl_vect_complex_free(&fft);
  pnl_vect_free(&vv);
  pnl_vect_free(&error);
}

static void pnl_fft_real2_test(void )
{  
  int               kmax  = 32;
  int               Nmax  = 2*kmax;
  double            L     = 4.0;
  double            Delta = L/(double) kmax;
  PnlVect          *vv, *inverse_fft;
  PnlVectComplex   *fft;
  double           *re, *im, error;
  int               j;

  vv          = pnl_vect_create(Nmax);
  fft         = pnl_vect_complex_create(Nmax);
  inverse_fft = pnl_vect_create(Nmax);
  re          = malloc (Nmax * sizeof (double));
  im          = malloc (Nmax * sizeof (double));
  /*function input for direct fourier transform */

  for (j=0; j<Nmax; j++)
    {
      double yj = (j-kmax)*Delta;
      pnl_vect_set(vv,j,cos(yj)*sin(yj));
      re[j] = cos(yj)*sin(yj);
    }

  pnl_real_fft (vv, fft);
  pnl_real_fft2 (re, im, Nmax);
  error = 0.;
  for (j=0; j<Nmax; j++)
    {
      error = MAX (error, fabs(re[j] - GET_REAL(fft, j)) +  fabs(im[j] - GET_IMAG(fft, j)));
    }

  pnl_test_eq_abs (error, 0., 1E-10, "Comparison real_fft and real_fft2" ,"");

  pnl_real_ifft2 (re, im, Nmax);
  error = 0.;
  for (j=0; j<Nmax; j++)
    {
      error = MAX (error, fabs(re[j] - GET(vv, j)) +  fabs(im[j]));
    }
  pnl_test_eq_abs (error, 0., 1E-10, "real_ifft2" ,"");


  pnl_vect_free(&inverse_fft);
  pnl_vect_complex_free(&fft);
  pnl_vect_free(&vv);
  free (re);
  free (im);
}

static void fft_real (PnlVect * a,int fft_size, int sign)
{
  int i, opposite;
  double last, second;
  switch (sign)
    {
    case 1 : /* backward */
      second = pnl_vect_get (a, 1);
      opposite = 1;
      for ( i=1 ; i<fft_size - 1 ; i++ )
        {
          pnl_vect_set (a, i, opposite * pnl_vect_get (a, i + 1));
          opposite = - opposite;
        }
      pnl_vect_set (a , fft_size - 1, second);

      pnl_real_ifft_inplace (a->array, fft_size);
      break;
    case -1 : /* forward */
      pnl_real_fft_inplace (a->array, fft_size);
      last = pnl_vect_get (a, fft_size - 1);
      opposite = -1;
      for ( i=fft_size - 1 ; i>1 ; i-- )
        {
          pnl_vect_set (a, i, opposite * pnl_vect_get (a, i - 1));
          opposite = - opposite;
        }
      pnl_vect_set (a , 1, last);
      break;
    }
}

static void pnl_fft_real_inplace_test(void )
{  
  int               kmax  = 32;
  int               Nmax  = 2*kmax;
  double            L     = 4.0;
  double            Delta = L/(double) kmax;
  PnlVect          *vv, *fft;
  int               j;

  vv          = pnl_vect_create(Nmax);
  fft         = pnl_vect_create(Nmax);
  /*function input for direct fourier transform */

  for (j=0; j<Nmax; j++)
    {
      double yj = j*Delta;
      pnl_vect_set(vv,j, exp (-yj * yj / 2));
    }

  pnl_vect_clone (fft, vv);
  fft_real(vv, vv->size, -1);
  fft_real(vv, vv->size, 1);

  pnl_test_vect_eq_abs (vv, fft, 1E-5, "real_fft composed with real_ifft", "");

  pnl_vect_free(&fft);
  pnl_vect_free(&vv);
}

static void pnl_fft_2d_inplace_test(void)
{  
  PnlMatComplex *data, *data_fft, *data_ifft, *expected;
  dcomplex *in , *expected_fft; 
  double in_tmp[24] = {1, 0, 2, 0, 3, 0, 9, 0, 8, 0, 5, 0, 1, 0, 2, 0, 9, 0, 8, 0, 7, 0, 2, 0};
  double expected_fft_tmp[24]= { 
    57.0000000000, 0.0000000000,
    7.0000000000, -2.0000000000,
    1.0000000000, 0.0000000000,
    7.0000000000, 2.0000000000,
    -6.0000000000, 8.6602540378,
    -3.9019237886, 7.1698729811,
    -11.0000000000, 3.4641016151,
    -9.0980762114, -15.8301270189,
    -6.0000000000, - 8.6602540378,
    -9.0980762114, 15.8301270189,
    -11.0000000000, -3.4641016151,
    -3.9019237886, -7.1698729811
  };

  in = (dcomplex *) in_tmp;
  expected_fft = (dcomplex *) expected_fft_tmp;

  data = pnl_mat_complex_create_from_ptr (3, 4, in);
  expected = pnl_mat_complex_create_from_ptr (3, 4, expected_fft);
  data_fft = pnl_mat_complex_copy (data);
  pnl_fft2d_inplace (data_fft);
  data_ifft = pnl_mat_complex_copy (data_fft);
  pnl_ifft2d_inplace (data_ifft);

  pnl_test_mat_complex_eq_abs (data_fft, expected, 1E-5, "pnl_fft2d", "");
  pnl_test_mat_complex_eq_abs (data_ifft, data, 1E-5, "pnl_ifft2d", "");

  pnl_mat_complex_free (&data);
  pnl_mat_complex_free (&data_fft);
  pnl_mat_complex_free (&data_ifft);
  pnl_mat_complex_free (&expected);
}

static void pnl_fft_2d_real_test(void)
{  
  PnlMatComplex *data_fft, *expected;
  PnlMat *data, *data_ifft;
  dcomplex *expected_fft; 
  double in[12] = {1, 2, 3, 9, 8, 5, 1, 2, 9, 8, 7, 2};
  double expected_fft_tmp[24]= { 
    57.0000000000, 0.0000000000,
    7.0000000000, -2.0000000000,
    1.0000000000, 0.0000000000,
    7.0000000000, 2.0000000000,
    -6.0000000000, 8.6602540378,
    -3.9019237886, 7.1698729811,
    -11.0000000000, 3.4641016151,
    -9.0980762114, -15.8301270189,
    -6.0000000000, - 8.6602540378,
    -9.0980762114, 15.8301270189,
    -11.0000000000, -3.4641016151,
    -3.9019237886, -7.1698729811
  };

  expected_fft = (dcomplex *) expected_fft_tmp;
  data_fft= pnl_mat_complex_new ();
  data = pnl_mat_create_from_ptr (3, 4, in);
  expected = pnl_mat_complex_create_from_ptr (3, 4, expected_fft);
  pnl_real_fft2d (data, data_fft);
  pnl_test_mat_complex_eq_abs (data_fft, expected, 1E-5, "pnl_real_fft2d", "");

  data_ifft = pnl_mat_new ();
  pnl_real_ifft2d (data_fft, data_ifft);
  pnl_test_mat_eq_abs (data_ifft, data, 1E-5, "pnl_real_ifft2d", "");

  pnl_mat_free (&data);
  pnl_mat_complex_free (&data_fft);
  pnl_mat_free (&data_ifft);
  pnl_mat_complex_free (&expected);
}

int main (int argc, char **argv)
{
  /*  dft(); */
  pnl_test_init (argc, argv);
  fftpack();
  pnl_fft_inplace_test();
  pnl_fft_real_test();
  pnl_fft_real2_test();
  pnl_fft_real_inplace_test();
  /*  dft_real(); */
  pnl_fft_2d_inplace_test ();
  pnl_fft_2d_real_test ();
  exit (pnl_test_finalize ("DFT"));
}
