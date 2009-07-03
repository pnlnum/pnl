
/*************************************************************************/
/* Written and (C) by                                                    */
/* David Pommier <david.pommier@gmail.com>                               */
/* Jérôme Lelong <jerome.lelong@gmail.com>                               */
/*                                                                       */
/* This program is free software; you can redistribute it and/or modify  */
/* it under the terms of the GNU General Public License as published by  */
/* the Free Software Foundation; either version 3 of the License, or     */
/* (at your option) any later version.                                   */
/*                                                                       */
/* This program is distributed in the hope that it will be useful,       */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/* GNU General Public License for more details.                          */
/*                                                                       */
/* You should have received a copy of the GNU General Public License     */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl_complex.h"
#include "pnl_mathtools.h"
#include "pnl_matrix.h"
#include "pnl_vector.h"
#include "pnl_fft.h"
#include "tests.h"


/* static void dft(void )
 * {  
 *   /\*
 *     Short program to test discret fourier transform ...
 *   *\/
 *  
 *   int kmax=32;/\*16*1024; *\/
 *   int Nmax=2*kmax;
 *   double L=4.0;
 *   double Delta=L/(double) kmax;
 *   PnlVect * y,*xi,*error;
 *   PnlVectComplex * vv,*vv2, *vxi,*vxi2,*result,*result2,*inverse_fft,*inverse_fft2;
 *   int j,i;
 *   printf("\nTest discrete fourier transform - Basic algorithm\n");
 *   y=pnl_vect_create(Nmax);
 *   xi=pnl_vect_create(Nmax);
 *   error=pnl_vect_create(Nmax);
 *   vv=pnl_vect_complex_create(Nmax);
 *   vv2=pnl_vect_complex_create(Nmax);
 *   vxi=pnl_vect_complex_create(Nmax);
 *   vxi2  = pnl_vect_complex_create(Nmax);
 *   result=pnl_vect_complex_create(Nmax);
 *   result2=pnl_vect_complex_create(Nmax);
 *   inverse_fft=pnl_vect_complex_create(Nmax);
 *   inverse_fft2=pnl_vect_complex_create(Nmax);
 *   /\*function input for direct fourier transform *\/
 *   j=0;
 *   do
 *     {
 *       LET(y,j)=(j-kmax)*Delta;
 *       pnl_vect_complex_set(vv,j,Complex(exp(-GET(y,j)*GET(y,j)/2.),0.0)); 
 *       pnl_vect_complex_set(vv2,j,Complex(-1.0*LET(y,j)*Creal(pnl_vect_complex_get(vv,j)),0.0));
 *       j+=1;
 *     }while(j<Nmax);
 *   
 *   for( j=0;j<Nmax;j++)
 *     {
 *       {
 *         LET(xi,j)=2.*M_PI*(j%kmax-kmax*(j/kmax))/(Delta*Nmax);
 *         pnl_vect_complex_set(vxi,j,Complex(exp(-GET(xi,j)*GET(xi,j)/2.),0.0));
 *         /\* the result of direct fourier transform for compare with the fft calculate result *\/
 *         pnl_vect_complex_set(vxi2,j,Complex(0.0,-GET(xi,j)*exp(-GET(xi,j)*GET(xi,j)/2.)));
 *         /\* the result of direct fourier transform for compare with the fft calculate result *\/
 *       }
 *     }
 *   pnl_dft_complex(vv,result);
 *   pnl_dft_complex(vv2,result2); 
 *   /\*direct fourier tranform *\/
 *   for(i=0; i< Nmax; i++)
 *     {
 *       LET(error,i)=Cabs(Csub(RCmul(Delta*M_1_SQRT2PI*pow((-1.), i),
 *                                    pnl_vect_complex_get(result,i))
 *                              ,pnl_vect_complex_get(vxi,i)));
 *       
 *     }
 *   printf("  Test real part error = %f\n",pnl_vect_max(error));
 *   for(i=0; i< Nmax; i++)
 *     {
 *       LET(error,i)=Cabs(Csub(RCmul(Delta*M_1_SQRT2PI*pow((-1.), i),
 *                                    pnl_vect_complex_get(result2,i)),pnl_vect_complex_get(vxi2,i)));
 *       
 *     }
 *   printf("  imaginary part error = %f\n",pnl_vect_max(error));
 *   printf("\nTest inverse discrete fourier transform - Basic algorithm \n");
 *   pnl_dft_complex_inverse(result,inverse_fft); 
 *   pnl_dft_complex_inverse(result2,inverse_fft2);
 *   for(i=0; i< Nmax; i++)
 *     {
 *       LET(error,i)=Cabs(Csub(pnl_vect_complex_get(inverse_fft,i),pnl_vect_complex_get(vv,i)));
 *     }
 *   printf("  inverse fft of pure real function, error = %f \n",pnl_vect_max(error));
 *   for(i=0; i< Nmax; i++)
 *     {
 *       LET(error,i)=Cabs(Csub(pnl_vect_complex_get(inverse_fft2,i),pnl_vect_complex_get(vv2,i)));
 *     }
 *   printf("  inverse fft of pure imaginary function, error = %f \n",pnl_vect_max(error));
 *   
 *   pnl_vect_complex_free(&inverse_fft);
 *   pnl_vect_complex_free(&inverse_fft2);
 *   pnl_vect_complex_free(&result);
 *   pnl_vect_complex_free(&result2);
 *   pnl_vect_complex_free(&vv);
 *   pnl_vect_complex_free(&vv2);
 *   pnl_vect_complex_free(&vxi);
 *   pnl_vect_complex_free(&vxi2);
 *   pnl_vect_free(&xi);
 *   pnl_vect_free(&y);
 *   pnl_vect_free(&error);
 *   
 * } */

static void fftpack(void )
{  
  /*
    Short program to test discret fourier transform ...
  */
 
  int kmax=32;
  int Nmax=2*kmax;
  double L=4.0;
  double Delta=L/(double) kmax;
  PnlVect * y,*xi,*error;
  PnlVectComplex * vv,*vv2, *vxi,*vxi2,*result,*result2,*inverse_fft,*inverse_fft2;
  int j,i;
  printf("\nFFTPack\n");
  y=pnl_vect_create(Nmax);
  xi=pnl_vect_create(Nmax);
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
  j=0;
  do
    {
      double yj = (j-kmax)*Delta;
      LET(y,j)=yj;
      pnl_vect_complex_set(vv,j,Complex(exp(-GET(y,j)*GET(y,j)/2.),0.0));
      pnl_vect_complex_set(vv2,j,Complex(-1.0*LET(y,j)*Creal(pnl_vect_complex_get(vv,j)),0.0));
      j+=1;
    }while(j<Nmax);
  
  for( j=0;j<Nmax;j++)
    {
      
      LET(xi,j)=2.*M_PI*(j%kmax-kmax*(j/kmax))/(Delta*Nmax);
      pnl_vect_complex_set(vxi,j,Complex(exp(-GET(xi,j)*GET(xi,j)/2.),0.0));
      /* the result of direct fourier transform for compare with the fft
         calculate result */
      pnl_vect_complex_set(vxi2,j,Complex(0.0,-GET(xi,j)*exp(-GET(xi,j)*GET(xi,j)/2.)));
      /* the result of direct fourier transform for compare with the fft
         calculate result */
      
    }
  pnl_fft(vv2,result2); 
  pnl_fft(vv,result);
  /*direct fourier tranform */
  for(i=0; i< Nmax; i++)
    {
      LET(error,i)=Cabs(Csub(RCmul(Delta*M_1_SQRT2PI*pow((-1.), i),
                                   pnl_vect_complex_get(result,i))
                             ,pnl_vect_complex_get(vxi,i)));
      
    }
  printf("  real part error = %f \n",pnl_vect_max(error));
  for(i=0; i< Nmax; i++)
    {
      LET(error,i)=Cabs(Csub(RCmul(Delta*M_1_SQRT2PI*pow((-1.), i),
                                   pnl_vect_complex_get(result2,i)),
                             pnl_vect_complex_get(vxi2,i)));
      
    }
  printf("  imaginary part error = %f \n",pnl_vect_max(error));
  printf("inverse discrete fourier transform - FFTPack \n");
  pnl_ifft(result,inverse_fft); 
  pnl_ifft(result2,inverse_fft2); 
  for(i=0; i< Nmax; i++)
    {
      LET(error,i)=Cabs(Csub(pnl_vect_complex_get(inverse_fft,i),pnl_vect_complex_get(vv,i)));
    }
  printf("  inverse fft of pure real function, error = %f \n",pnl_vect_max(error));
  for(i=0; i< Nmax; i++)
    {
      LET(error,i)=Cabs(Csub(pnl_vect_complex_get(inverse_fft2,i),pnl_vect_complex_get(vv2,i)));
    }
  printf("  inverse fft of pure imaginary function, error = %f \n",pnl_vect_max(error));
  
  pnl_vect_complex_free(&inverse_fft);
  pnl_vect_complex_free(&inverse_fft2);
  pnl_vect_complex_free(&result);
  pnl_vect_complex_free(&result2);
  pnl_vect_complex_free(&vv);
  pnl_vect_complex_free(&vv2);
  pnl_vect_complex_free(&vxi);
  pnl_vect_complex_free(&vxi2);
  pnl_vect_free(&xi);
  pnl_vect_free(&y);
  pnl_vect_free(&error);
  
}

static void pnl_fft_inplace_test(void )
{  
 
  int               kmax  = 32;
  int               Nmax  = 2*kmax;
  double            L     = 4.0;
  double            Delta = L/(double) kmax;
  PnlVect          *error;
  PnlVectComplex   *vv, *result,*inverse_fft;
  double           *re, *im;
  int               j,i;
  printf("\nFFTPack -- in place\n");
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

  /* printf ("res = ");
     pnl_vect_complex_print_nsp (result); */
  for (j=0; j<Nmax; j++)
    {
      if ( GET_REAL(result, j) != re[j] || GET_IMAG(result, j) != im[j] )
        printf ("error in pnl_fft or pnl_fft2\n");
    }

  pnl_vect_complex_clone (inverse_fft, result);
  pnl_ifft_inplace(inverse_fft);
  pnl_ifft2(re, im , Nmax);

  for (j=0; j<Nmax; j++)
    {
      if ( GET_REAL(inverse_fft, j) != re[j] || GET_IMAG(inverse_fft, j) != im[j] )
        {
          printf ("error in pnl_fft or pnl_fft2\n");
          break;
        }
    }
  if (j == Nmax) {printf ("  orignal sequence perfectly recovered \n");}
    
  for(i=0; i< Nmax; i++)
    {
      LET(error,i)=Cabs(Csub(pnl_vect_complex_get(inverse_fft,i),pnl_vect_complex_get(vv,i)));
    }
  printf("  inverse fft, error = %f \n",pnl_vect_max(error));
  for (j=0; j<Nmax; j++)
    {
      if ( GET_REAL(inverse_fft, j) != re[j] || GET_IMAG(inverse_fft, j) != im[j] )
        {
          printf ("error in pnl_fft or pnl_fft2\n");
          break;
        }
    }
  if (j == Nmax) {printf ("  orignal sequence perfectly recovered \n");}
  
  pnl_vect_complex_free(&inverse_fft);
  pnl_vect_complex_free(&result);
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

  printf("\nFFTPack real sequences -- pnl_real_fft\n");

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

  /*  printf ("res = ");
      pnl_vect_complex_print_nsp (fft); */

  pnl_real_ifft(fft, inverse_fft);
  for (j=0; j<Nmax; j++)
    {
      if ( fabs(GET(inverse_fft, j) - GET(vv, j)) > 1E-10 )
        {
          printf ("error in pnl_real_fft\n");
          pnl_vect_minus_vect (inverse_fft, vv);
          pnl_vect_map_inplace (inverse_fft, fabs);
          printf ("recovering error : %f\n", pnl_vect_max(inverse_fft));
          break;
        }
    }
  if (j == Nmax) {printf ("  orignal sequence perfectly recovered \n");}

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

  printf("\nFFTPack real sequences -- pnl_real_fft2\n");

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

  if (error > 0.)
    {
      printf ("  difference between pnl_real_fft and pnl_real_fft2 : %f\n", error);
    }
  else {printf ("  perfect match between pnl_real_fft and pnl_real_fft2\n");}

  pnl_real_ifft2 (re, im, Nmax);
  error = 0.;
  for (j=0; j<Nmax; j++)
    {
      error = MAX (error, fabs(re[j] - GET(vv, j)) +  fabs(im[j]));
    }
  
  if (error > 1.E-10)
    {
      printf ("  recovering error in pnl_real_ifft2 : %f\n", error);
    }
  else {printf ("  initial data perfectly recovered\n");}

  pnl_vect_free(&inverse_fft);
  pnl_vect_complex_free(&fft);
  pnl_vect_free(&vv);
  free (re);
  free (im);
}


/* static void dft_real(void )
 * {  
 *   /\*
 *     Short program to test discret fourier transform on real functions  ...
 *   *\/
 *   int kmax=32;/\*16*1024; *\/
 *   int Nmax=2*kmax;
 *   double L=4.0;
 *   double Delta=L/(double) kmax;
 *   PnlVect * y,*xi,*error,*error2, * vv, * vv2, *vxi, *vxi2,*result,*result2,*inverse_fft,*inverse_fft2;
 *   int j,i;
 *   printf("\nTest discrete fourier transform - Creal FFT algorithm \n");
 *   y=pnl_vect_create(Nmax);
 *   xi=pnl_vect_create(kmax);
 *   error=pnl_vect_create(kmax);
 *   error2=pnl_vect_create(Nmax);
 *   vv=pnl_vect_create(Nmax);
 *   vxi=pnl_vect_create(kmax);
 *   vv2=pnl_vect_create(Nmax);
 *   vxi2=pnl_vect_create(kmax);
 *   /\*function input for direct fourier transform *\/
 *   j=0;
 *   do
 *     {
 *       LET(y,j)=(j-kmax)*Delta;
 *       LET(vv,j)=exp(-GET(y,j)*GET(y,j)/2.); 
 *       LET(vv2,j)=-GET(y,j)*exp(-GET(y,j)*GET(y,j)/2.); 
 *       j+=1;
 *     }while(j<Nmax);
 *   for( j=0;j<kmax;j++)
 *     {
 *       {
 *         LET(xi,j)=2.*M_PI*j/(Delta*Nmax);
 *         LET(vxi,j)=exp(-GET(xi,j)*GET(xi,j)/2.);    
 *         /\* the real result of direct fourier transform for compare with the realfft calculate result *\/
 *         LET(vxi2,j)=-GET(xi,j)*exp(-GET(xi,j)*GET(xi,j)/2.);
 *         /\* the imaginary result of direct fourier transform for compare with the realfft calculate result *\/
 *       }
 *     }
 *   result=pnl_vect_copy(vv);
 *   result2=pnl_vect_copy(vv2);
 *   pnl_fft_real(result, Nmax, -1); 
 *   pnl_fft_real(result2, Nmax, -1); 
 *   /\*direct fourier tranform *\/
 *   for(i=0; i< kmax; i++)
 *     {
 *       LET(error,i)=Cabs(Csub(RCmul(Delta*M_1_SQRT2PI*pow((-1.), i),
 *                                    Complex(LET(result,2*i),LET(result,2*i+1)))
 *                              ,Complex(LET(vxi,i),0.)));
 *     }
 *   printf("  real part error = %f \n",pnl_vect_max(error));
 *   for(i=0; i< kmax; i++)
 *     {
 *       LET(error,i)=Cabs(Csub(RCmul(Delta*M_1_SQRT2PI*pow((-1.), i),
 *                                    Complex(LET(result2,2*i),LET(result2,2*i+1)))
 *                              ,Complex(0.0,LET(vxi2,i))));
 *     }
 *   printf("  imaginary part error = %f \n",pnl_vect_max(error));
 *   printf("Test inverse discrete fourier transform - Creal FFT algorithm \n");
 *   inverse_fft=pnl_vect_copy(result);
 *   pnl_fft_real(inverse_fft, Nmax,1); 
 *   inverse_fft2=pnl_vect_copy(result2);
 *   pnl_fft_real(inverse_fft2, Nmax,1); 
 *   for(i=0; i< Nmax; i++)
 *     {
 *       LET(error2,i)=fabs(LET(inverse_fft,i)-LET(vv,i));
 *     }
 *   printf("  inverse fft of pure real function,error = %f \n",pnl_vect_max(error2));
 *   for(i=0; i< Nmax; i++)
 *     {
 *       LET(error2,i)=fabs(LET(inverse_fft2,i)-LET(vv2,i));
 *     }
 *   printf("  inverse fft of pure imaginary function,error = %f \n",pnl_vect_max(error2));
 *   pnl_vect_free(&result2);
 *   pnl_vect_free(&vv2);
 *   pnl_vect_free(&vxi2);
 *   pnl_vect_free(&result);
 *   pnl_vect_free(&vv);
 *   pnl_vect_free(&vxi);
 *   pnl_vect_free(&xi);
 *   pnl_vect_free(&y);
 *   pnl_vect_free(&error);
 *   pnl_vect_free(&error2);
 *   pnl_vect_free(&inverse_fft);
 *   pnl_vect_free(&inverse_fft2);
 * } */

void dft_test()
{
  /*  dft(); */
  fftpack();
  pnl_fft_inplace_test();
  pnl_fft_real_test();
  pnl_fft_real2_test();
  /*  dft_real(); */
}
