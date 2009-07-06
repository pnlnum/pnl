#include "pnl_fft.h"

/*
 * This file must not be compiled within PNL
 * It is only included in ../examples/dft_test.c for comparing old FFT
 * routines with FTPack
 */

extern void pnl_dft_complex_transform(const PnlVectComplex * data, 
                                      PnlVectComplex * result,
                                      const int sign);
extern void pnl_fft_real(PnlVect * a,int fft_size, int sign);

static void swap_double(double * x, double * y)
{
  double tmp = *x;
  *x = *y; *y = tmp;
}


/**
 * classical fast fourier tansform in \f$ O(N^2)\f$
 *
 * @param data   : function on grid
 * @param result : FFT(data)
 */
void pnl_dft_complex(const PnlVectComplex * data, 
                     PnlVectComplex * result)           
{ pnl_dft_complex_transform(data,result,1.0);}


/**
 * classical inverse fast fourier tansform without \f$ L_1\f$ normalisation in \f$ O(N^2)\f$
 *
 * @param data   : function on grid
 * @param result :  FFT^{-1}(data) 
 */
void pnl_dft_complex_backward(const PnlVectComplex * data, 
                              PnlVectComplex * result)
{ pnl_dft_complex_transform(data,result,-1.0);}

/**
 * classical inverse fast fourier tansform on complex in \f$ O(N^2)\f$
 *
 * @param data   : function on grid
 * @param result :  FFT^{-1}(data)
 */
void pnl_dft_complex_inverse(const PnlVectComplex * data, 
                             PnlVectComplex * result)
{
  pnl_dft_complex_transform(data,result,-1.0);
  pnl_vect_complex_mult_double(result,1.0/( double) (data->size));
}

/**
 * classical fast fourier tansform
 *
 * @param data   : function on grid
 * @param result :  FFT(data) or FFT^{-1}(data)
 * @param sign   : choice betwenn FFT (sign = -1) and FFT^{-1} (sign = 1)
 */
void pnl_dft_complex_transform(const PnlVectComplex * data, 
                               PnlVectComplex * result,
                               const int sign)
{
  int i, j, exponent;
  PnlVectComplex *twidle = pnl_vect_complex_create (data->size);
  const double d_theta = ((double) sign)*M_2PI/(double) (data->size);
  pnl_vect_complex_resize (result, data->size);

  for (i = 0; i < data->size; i++)
    {
      double theta = d_theta * (double) i;
      double w_real = cos (theta);
      double w_imag = sin (theta);
      pnl_vect_complex_set ( twidle, i, Complex(w_real, w_imag) );
    }

  
  for (i = 0; i < data->size; i++)
    {
      double sum_real = 0.;
      double sum_imag = 0.;

      exponent = 0;
      for (j = 0; j < (data->size); j++)
        {
          /* sum = exp(i theta) * data[j] */
          fcomplex w = pnl_vect_complex_get (twidle, (i*j)%data->size);
          double w_real = Creal(w);
          double w_imag = Cimag(w);
          fcomplex tmp = pnl_vect_complex_get(data,j);
          double data_real = Creal (tmp);
          double data_imag = Cimag (tmp);

          sum_real += w_real * data_real - w_imag * data_imag;
          sum_imag += w_real * data_imag + w_imag * data_real;
          
        }
      pnl_vect_complex_set(result,i,Complex(sum_real, sum_imag));
    }
  pnl_vect_complex_free (&twidle);
}


/**
 * in-place discret fourier transform for real functions
 * Extract for Numerical Recipes Paragraph 12.3 and adapted 
 *
 * Cooley-Tukey algorithm
 *
 * Calculates the Fourier transform of a set of n real-valuated data points
 * Replace these data stored in two Vector for rela and imaginaty part 
 * Routine also compute inverse Fourier Transform with same convention 
 *
 * @param a : real part of function 
 * @param b : image part of function
 * @param n : size of sub vector in which we do fft 
 * @param sign an integer to specify if FFT (sign=-1) or FFT^{-1} (sign=1) is
 * performed forward (-1) or backward (+1).
 * for lengths which are a power of 2
 */
void pnl_dft_real_expanded(double *a, double *b, int n, int sign)
{
  int i=0;
  PnlVect * data;
  data= pnl_vect_create(0);
  data->owner=0;
  data->size=n;
  if(sign==-1)
    {
          data->array=a;
          pnl_fft_real(data,n,sign);
          /* Distribution des valeurs: */
          /*-n/2+1 -n/2+2  0 1 2 n/2-2 n/2-1 */
          /* j=0=-n/2+1+i avec i = n/2-1 */
          b[0]=GET(data,0);
          b[n/2]=GET(data,1);
          for (i=1; i<n/2; i++)
            {
              b[i]=sign*GET(data,2*i+1);
              b[n-i]=GET(data,2*i);
            }
          a[0]=b[0];
          b[0]=0;
          a[n/2]=b[n/2];
          b[n/2]=0;
          for (i=1; i<n/2; i++)
            {
              a[i]=b[n-i];
              a[n-i]=a[i];
              b[n-i]= -b[i];
            }
        }
  else
    {
      data->array=b;
      for (i=n/2-1; i>0; i--)
        {
          LET(data,2*i+1)=sign*b[i];
          LET(data,2*i)=a[i];
        }
      b[0] = a[0];
      b[1] = a[n/2];
      pnl_fft_real(data,n,sign);
      a[0]=GET(data,0);
      for (i=1; i<n; i++)
        a[i]=GET(data,n-i);
     }
  pnl_vect_free(&data);
 
}

/**
 * in-place discret fourier transform for real functions
 * Extract for Numerical Recipes Paragraph 12.3 and adapted 
 *
 * Cooley-Tukey algorithm
 *
 * Calculates the Fourier transform of a set of n real-valuated data points
 * Replace these data stored in a Vector by the positive 
 * frequency half of its complex Fourier transform.
 * if i is even then \f$ a(2*i)+ \imath a(2i+1) \f$ is the ième    
 * coefficients of the Fourier transform 
 * if i is odd then \f$ a(N-2*i)+ \imath a(N-2i+1) \f$ is the ième    
 * coefficients of the Fourier transform
 * Routine also compute inverse Fourier Transform with same convention 
 *
 * @param a : real function on grid
 * @param fft_size : size of sub vector in which we do fft 
 * @param sign an integer to specify if FFT (sign=-1) or FFT^{-1} (sign=1) is
 * performed forward (-1) or backward (+1).
 * for lengths which are a power of 2
 */
void pnl_fft_real(PnlVect * a,int fft_size, int sign)
{
  double ttheta,theta,c1,c2,h1r;
  int half_fft_size;
  int i, i1,ii,jj, mmax, m,j,istep,isign;
  fcomplex tw,tw0,w,wp0,temp, h1,h2,*A_ptr,*B_ptr;
  if( fft_size==1 )
    {return;}
  ttheta = M_2PI/(double)(fft_size);
  c1 = 0.5;
  c2 = (sign==1)?0.5:-0.5;
  if (sign==1)
    {
      ttheta *= -1;
      tw0=Complex(-2.0*SQR(sin(0.5*ttheta)),sin(ttheta));
      tw=Cadd(CONE,tw0);
      for(i = 2; i <= fft_size/4+1; i++)
        {
          i1 = i+i-2;
          A_ptr=(fcomplex*) pnl_vect_lget(a,i1);
          B_ptr=(fcomplex*) pnl_vect_lget(a,fft_size-i1);
          h1=C_op_dapcb(c1,*A_ptr,*B_ptr);
          h2=C_op_idamcb(c2,*A_ptr,*B_ptr);
          temp=Cmul(tw,h2);
          *A_ptr = Cadd(h1,temp);
          *B_ptr= Conj(Csub(h1,temp));
          tw=Cadd(tw,Cmul(tw,tw0));/*Caddegal(tw,Cmul(tw,tw0)); */
        }
      h1r = GET(a,0);
      LET(a,0) = c1*(h1r+GET(a,1));
      LET(a,1) = c1*(h1r-GET(a,1));
    }
  isign = -sign; 
  half_fft_size = fft_size/2;
  j = 1;
  for(ii = 1; ii <= half_fft_size; ii++)
    {
      i = 2*ii-1;
      if( j>i )
        {
          swap_double(pnl_vect_lget(a,j-1),pnl_vect_lget(a,i-1));
          swap_double(pnl_vect_lget(a,j),pnl_vect_lget(a,i));
        }
      m = half_fft_size;
      while(m>=2&&j>m)
        {
          j -= m;
          m /= 2;
        }
      j = j+m;
    }
  mmax = 2;
  while(fft_size>mmax)
    {
      istep = 2*mmax;
      theta = 2*M_PI/(isign*mmax);
      w=CONE;
      wp0=Complex(-2.0*SQR(sin(0.5*theta)),sin(theta));
      for(ii = 1; ii <= mmax/2; ii++)
        {
          m = 2*ii-1;
          for(jj = 0; jj <= (fft_size-m)/istep; jj++)
            {
              i = m+jj*istep;
              A_ptr=(fcomplex*) pnl_vect_lget(a,i+mmax-1);
              B_ptr=(fcomplex*) pnl_vect_lget(a,i-1);
              temp=Cmul(w,*A_ptr);
              *A_ptr= Csub(*B_ptr,temp);
              *B_ptr= Cadd(*B_ptr,temp);
            }
          w=Cadd(w,Cmul(w,wp0));/*Caddegal(w,Cmul(w,wp0)); */
        }
      mmax = istep;
    }
  if(sign == 1)
    {
      for(i = 1; i <= 2*half_fft_size; i++)
        { LET(a,i-1) /= half_fft_size; }
    }
  if(sign == -1)
    {
      tw0=Complex(-2.0*SQR(sin(0.5*ttheta)),sin(ttheta));
      tw=Cadd(CONE,tw0);
      for(i = 2; i <= fft_size/4+1; i++)
        {
          
            i1 = i+i-2;
            A_ptr=(fcomplex*) pnl_vect_lget(a,i1);
            B_ptr=(fcomplex*) pnl_vect_lget(a,fft_size-i1);
            h1=C_op_dapcb(c1,*A_ptr,*B_ptr);
            h2=C_op_idamcb(c2,*A_ptr,*B_ptr);
            temp=Cmul(tw,h2);
            *A_ptr = Cadd(h1,temp);
            *B_ptr= Conj(Csub(h1,Cmul(tw,h2)));
            tw=Cadd(tw,Cmul(tw,tw0));/*Caddegal(tw,Cmul(tw,tw0)); */
        }
      h1r = GET(a,0);
      LET(a,0) = h1r+GET(a,1);
      LET(a,1) = h1r-GET(a,1);
    }
}


/****

     Blablab 

*/


/* static int pnl_fft_complex_pass_2 (PnlVectComplex * in,
 *                                    PnlVectComplex * out,
 *                                    const int sign,
 *                                    const uint product,
 *                                    const uint n,
 *                                    PnlVectComplex *twiddle)
 * {
 *   uint i = 0, j = 0;
 *   uint k, k1;
 * 
 *   const uint factor = 2;
 *   const uint m = n / factor;
 *   const uint q = n / product;
 *   const uint product_1 = product / factor;
 *   const uint jump = (factor - 1) * product_1;
 *   for (k = 0; k < q; k++)
 *     {
 *       fcomplex w;
 *       if (k == 0)
 *         {
 *           w =CONE;
 *         }
 *       else
 *         {
 *           if (sign == -1)
 *             {
 *               /\* forward tranform *\/
 *               w = pnl_vect_complex_get(twiddle,k - 1);
 *             }
 *           else
 *             {
 *               /\* backward tranform: w -> conjugate(w) *\/
 *               w = Conj(pnl_vect_complex_get(twiddle,k - 1));
 *             }
 *         }
 * 
 *       for (k1 = 0; k1 < product_1; k1++)
 *         {
 *           const fcomplex z0 = pnl_vect_complex_get(in,i);
 *           const fcomplex z1 = pnl_vect_complex_get(in,i+m);
 * 
 *           /\* compute x = W(2) z *\/
 * 
 *           /\* x0 = z0 + z1 *\/
 *           const fcomplex x0 =Cadd(z0,z1);
 *           /\* x1 = z0 - z1 *\/
 *           const fcomplex x1 = Csub(z0,z1);
 *           /\* apply twiddle factors *\/
 *           
 *           /\* out0 = 1 * x0 *\/
 *           pnl_vect_complex_set(out,j, x0);
 *                    
 *           /\* out1 = w * x1 *\/
 *           pnl_vect_complex_set(out,j+product_1,Cmul(w,x1));
 *           i++;
 *           j++;
 *         }
 *       j += jump;
 *     }
 *   return 0;
 * }
 * 
 * static int
 * pnl_fft_complex_pass_3(PnlVectComplex * in,
 *                        PnlVectComplex * out,
 *                        const int sign,
 *                        const uint product,
 *                        const uint n,
 *                        const PnlVectComplex *twiddle1,
 *                        const PnlVectComplex *twiddle2)
 * {
 *   uint i = 0, j = 0;
 *   uint k, k1;
 * 
 *   const uint factor = 3;
 *   const uint m = n / factor;
 *   const uint q = n / product;
 *   const uint product_1 = product / factor;
 *   const uint jump = (factor - 1) * product_1;
 * 
 *   const double tau = sqrt (3.0) / 2.0;
 * 
 *   for (k = 0; k < q; k++)
 *     {
 *       fcomplex w1,  w2;
 * 
 *       if (k == 0)
 *         {
 *           w1 = CONE;
 *           w2 = CONE;
 *         }
 *       else
 *         {
 *           if (sign == -1)
 *             {
 *               /\* forward tranform *\/
 *               w1 = pnl_vect_complex_get(twiddle1,k - 1);
 *               w2 = pnl_vect_complex_get(twiddle2,k - 1);
 *             }
 *           else
 *             {
 *               /\* backward tranform: w -> conjugate(w) *\/
 *               w1 = Conj(pnl_vect_complex_get(twiddle1,k - 1));
 *               w2 = Conj(pnl_vect_complex_get(twiddle2,k - 1));
 *             }
 *         }
 * 
 *       for (k1 = 0; k1 < product_1; k1++)
 *         {
 *           const fcomplex z0 = pnl_vect_complex_get(in,i);
 *           const fcomplex z1 = pnl_vect_complex_get(in,i+m);
 *           const fcomplex z2 = pnl_vect_complex_get(in,i+2*m);
 *           /\* compute x = W(3) z *\/
 *           /\* t1 = z1 + z2 *\/
 *           const fcomplex t1 = Cadd(z1,z2);
 *           /\* t2 = z0 - t1/2 *\/
 *           const fcomplex t2 = Csub(z0,RCmul(0.5,t1));
 *           /\* t3 = (+/-) sin(pi/3)*(z1 - z2) *\/
 *           const fcomplex t3 = RCmul(sign * tau,Csub(z1,z2));
 *           /\* x0 = z0 + t1 *\/
 *           const fcomplex x0 = Cadd(z0,t1);
 *           /\* x1 = t2 + i t3 *\/
 *           const fcomplex x1 = C_op_apib(t2,t3);
 *           /\* x2 = t2 - i t3 *\/
 *           const fcomplex x2 = C_op_amib(t2,t3);
 *           /\* apply twiddle factors *\/
 *           /\* to0 = 1 * x0 *\/
 *           pnl_vect_complex_set(out,j,x0);
 *           /\* to1 = w1 * x1 *\/
 *           pnl_vect_complex_set(out,j+product_1,Cmul(w1,x1));
 *           /\* to2 = w2 * x2 *\/
 *           pnl_vect_complex_set(out,j+2*product_1,Cmul(w2,x2));
 *           i++; j++;
 *         }
 *       j += jump;
 *     }
 *   return 0;
 * }
 * 
 * static int
 * pnl_fft_complex_pass_4 (PnlVectComplex * in,
 *                         PnlVectComplex * out,
 *                         const int sign,
 *                         const uint product,
 *                         const uint n,
 *                         const PnlVectComplex * twiddle1,
 *                         const PnlVectComplex * twiddle2,
 *                         const PnlVectComplex * twiddle3)
 * {
 *   uint i = 0, j = 0;
 *   uint k, k1;
 * 
 *   const uint factor = 4;
 *   const uint m = n / factor;
 *   const uint q = n / product;
 *   const uint p_1 = product / factor;
 *   const uint jump = (factor - 1) * p_1;
 * 
 *   for (k = 0; k < q; k++)
 *     {
 *       fcomplex w1, w2, w3;
 * 
 *       if (k == 0)
 *         {
 *           w1=CONE;
 *           w2=CONE;
 *           w3=CONE;
 *         }
 *       else
 *         {
 *           if (sign == -1)
 *             {
 *               /\* forward tranform *\/
 *               w1 =pnl_vect_complex_get(twiddle1,k - 1);
 *               w2 =pnl_vect_complex_get(twiddle2,k - 1);
 *               w3 =pnl_vect_complex_get(twiddle3,k - 1);
 *               
 *             }
 *           else
 *             {
 *               /\* backward tranform: w -> conjugate(w) *\/
 *               w1 =Conj(pnl_vect_complex_get(twiddle1,k - 1));
 *               w2 =Conj(pnl_vect_complex_get(twiddle2,k - 1));
 *               w3 =Conj(pnl_vect_complex_get(twiddle3,k - 1));
 *             }
 *         }
 * 
 *       for (k1 = 0; k1 < p_1; k1++)
 *         {
 *           const fcomplex z0 = pnl_vect_complex_get(in,i);
 *           const fcomplex z1 = pnl_vect_complex_get(in,i+m);
 *           const fcomplex z2 = pnl_vect_complex_get(in,i+2*m);
 *           const fcomplex z3 = pnl_vect_complex_get(in,i+3*m);
 * 
 *           /\* compute x = W(4) z *\/
 *           
 *           /\* t1 = z0 + z2 *\/
 *           const fcomplex t1 = Cadd(z0,z2);
 *           /\* t2 = z1 + z3 *\/
 *           const fcomplex t2 = Cadd(z1,z3);
 *           /\* t3 = z0 - z2 *\/
 *           const fcomplex t3 = Csub(z0,z2);
 *           /\* t4 = (+/-) (z1 - z3) *\/
 *           const fcomplex t4 = RCmul( sign,Csub(z1,z3));
 *           /\* x0 = t1 + t2 *\/
 *           const fcomplex x0 = Cadd(t1,t2);
 *           /\* x1 = t3 + i t4 *\/
 *           const fcomplex x1 = C_op_apib(t3,t4);
 *           /\* x2 = Csub(t1,t2) *\/
 *           const fcomplex x2 = Csub(t1,t2);
 *           /\* x3 = t3 - i t4 *\/
 *           const fcomplex x3 = C_op_amib(t3,t4);
 *           /\* apply twiddle factors *\/
 * 
 *           /\* to0 = 1 * x0 *\/
 *           pnl_vect_complex_set(out,j,x0);
 *           /\* to1 = w1 * x1 *\/
 *           pnl_vect_complex_set(out, j + p_1, Cmul(w1,x1));
 *           /\* to2 = w2 * x2 *\/
 *           pnl_vect_complex_set(out,j + 2 * p_1, Cmul(w2,x2));
 *           /\* to3 = w3 * x3 *\/
 *           pnl_vect_complex_set(out,j + 3 * p_1, Cmul(w3,x3));
 *           i++;
 *           j++;
 *         }
 *       j += jump;
 *     }
 *   return 0;
 * }
 * 
 * static int
 * pnl_fft_complex_pass_5 (PnlVectComplex * in,
 *                         PnlVectComplex * out,
 *                         const int sign,
 *                         const uint product,
 *                         const uint n,
 *                         const PnlVectComplex * twiddle1,
 *                         const PnlVectComplex * twiddle2,
 *                         const PnlVectComplex * twiddle3,
 *                         const PnlVectComplex * twiddle4)
 * {
 *   uint i = 0, j = 0;
 *   uint k, k1;
 * 
 *   const uint factor = 5;
 *   const uint m = n / factor;
 *   const uint q = n / product;
 *   const uint p_1 = product / factor;
 *   const uint jump = (factor - 1) * p_1;
 * 
 *   const double sin_2pi_by_5 = sin (2.0 * M_PI / 5.0);
 *   const double sin_2pi_by_10 = sin (2.0 * M_PI / 10.0);
 * 
 *   for (k = 0; k < q; k++)
 *     {
 * 
 *       fcomplex w1, w2, w3, w4;
 * 
 *       if (k == 0)
 *         {
 *           w1=CONE;
 *           w2=CONE;
 *           w3=CONE;
 *           w4=CONE;
 *         }
 *       else
 *         {
 *           if (sign == -1)
 *             {
 *               /\* forward tranform *\/
 *               w1 =pnl_vect_complex_get(twiddle1,k - 1);
 *               w2 =pnl_vect_complex_get(twiddle2,k - 1);
 *               w3 =pnl_vect_complex_get(twiddle3,k - 1);
 *               w4 =pnl_vect_complex_get(twiddle4,k - 1);
 *                           
 *             }
 *           else
 *             {
 *               /\* backward tranform: w -> conjugate(w) *\/
 *               w1 =Conj(pnl_vect_complex_get(twiddle1,k - 1));
 *               w2 =Conj(pnl_vect_complex_get(twiddle2,k - 1));
 *               w3 =Conj(pnl_vect_complex_get(twiddle3,k - 1));
 *               w4 =Conj(pnl_vect_complex_get(twiddle4,k - 1));
 *             }
 *         }
 * 
 *       for (k1 = 0; k1 < p_1; k1++)
 *         {
 *           const fcomplex z0 = pnl_vect_complex_get(in,i);
 *           const fcomplex z1 = pnl_vect_complex_get(in,i+m);
 *           const fcomplex z2 = pnl_vect_complex_get(in,i+2*m);
 *           const fcomplex z3 = pnl_vect_complex_get(in,i+3*m);
 *           const fcomplex z4 = pnl_vect_complex_get(in,i+4*m);
 *       
 *           /\* compute x = W(5) z *\/
 * 
 *           /\* t1 = z1 + z4 *\/
 *           const fcomplex t1 = Cadd(z1,z4);
 *           /\* t2 = z2 + z3 *\/
 *           const fcomplex t2 = Cadd(z2,z3);
 *           /\* t3 = z1 - z4 *\/
 *           const fcomplex t3 = Csub(z1,z4);
 *           /\* t4 = z2 - z3 *\/
 *           const fcomplex t4 = Csub(z2,z3);
 *           /\* t5 = t1 + t2 *\/
 *           const fcomplex t5 = Cadd(t1,t2);
 *           /\* t6 = (sqrt(5)/4)(t1 - t2) *\/
 *           const fcomplex t6 = RCmul((sqrt (5.0) / 4.0) , (Csub(t1,t2)));
 *           /\* t7 = z0 - ((t5)/4) *\/
 *           const fcomplex t7 = Csub(z0,RCmul(0.25,t5));
 *           /\* t8 = t7 + t6 *\/
 *           const fcomplex t8 = Cadd(t7,t6);
 *           /\* t9 = t7 - t6 *\/
 *           const fcomplex t9 = Csub(t7,t6);
 *           /\* t10 = sin(2 pi/5) t3 + sin(2 pi/10) t4 *\/
 *           const fcomplex t10 = RCmul(sign,Cadd(RCmul(sin_2pi_by_5, t3),RCmul(sin_2pi_by_10,t4)));
 *           /\* t11 = sin(2 pi/10) t3 - sin(2 pi/5) t4 *\/
 *           const fcomplex t11 = RCmul(sign,Csub(RCmul(sin_2pi_by_10, t3),RCmul(sin_2pi_by_5,t4)));
 *           /\* x0 = z0 + t5 *\/
 *           fcomplex x0 = Cadd(z0,t5);
 *           /\* x1 = t8 + i t10 *\/
 *           fcomplex x1 = C_op_apib(t8,t10);
 *           /\* x2 = t9 + i t11 *\/
 *           fcomplex x2 = C_op_apib(t9, t11);
 *           /\* x3 = t9 - i t11 *\/
 *           fcomplex x3 = C_op_amib(t9, t11);
 *           /\* x4 = t8 - i t10 *\/
 *           fcomplex x4 = C_op_amib(t8, t10);
 *           /\* apply twiddle factors *\/
 *           
 *           /\* to0 = 1 * x0 *\/
 *           pnl_vect_complex_set(out,j, x0);
 *           /\* to1 = w1 * x1 *\/
 *           pnl_vect_complex_set(out, j + p_1, Cmul(w1,x1));
 *           /\* to2 = w2 * x2 *\/
 *           pnl_vect_complex_set(out,j + 2 * p_1, Cmul(w2,x2));
 *           /\* to3 = w3 * x3 *\/
 *           pnl_vect_complex_set(out,j + 3 * p_1, Cmul(w3,x3));
 *           /\* to4 = w4 * x4 *\/
 *           pnl_vect_complex_set(out,j + 4 * p_1, Cmul(w4,x4));
 *           
 *           i++;
 *           j++;
 *         }
 *       j += jump;
 *     }
 *   return 0;
 * }
 * 
 * static int
 * pnl_fft_complex_pass_6 (PnlVectComplex * in,
 *                         PnlVectComplex * out,
 *                         const int sign,
 *                         const uint product,
 *                         const uint n,
 *                         const PnlVectComplex * twiddle1,
 *                         const PnlVectComplex * twiddle2,
 *                         const PnlVectComplex * twiddle3,
 *                         const PnlVectComplex * twiddle4,
 *                         const PnlVectComplex * twiddle5)
 * {
 * 
 *   uint i = 0, j = 0;
 *   uint k, k1;
 * 
 *   const uint factor = 6;
 *   const uint m = n / factor;
 *   const uint q = n / product;
 *   const uint p_1 = product / factor;
 *   const uint jump = (factor - 1) * p_1;
 * 
 *   const double tau = sqrt (3.0) / 2.0;
 * 
 *   for (k = 0; k < q; k++)
 *     {
 *       fcomplex w1, w2, w3, w4, w5;
 *       if (k == 0)
 *         {
 *           w1=CONE;
 *           w2=CONE;
 *           w3=CONE;
 *           w4=CONE;
 *           w5=CONE;
 *         }
 *       else
 *         {
 *           if (sign == -1)
 *             {
 *               /\* forward tranform *\/
 *               w1 =pnl_vect_complex_get(twiddle1,k - 1);
 *               w2 =pnl_vect_complex_get(twiddle2,k - 1);
 *               w3 =pnl_vect_complex_get(twiddle3,k - 1);
 *               w4 =pnl_vect_complex_get(twiddle4,k - 1);
 *               w5 =pnl_vect_complex_get(twiddle5,k - 1);
 *              
 *             }
 *           else
 *             {
 *               /\* backward tranform: w -> conjugate(w) *\/
 *               w1 =Conj(pnl_vect_complex_get(twiddle1,k - 1));
 *               w2 =Conj(pnl_vect_complex_get(twiddle2,k - 1));
 *               w3 =Conj(pnl_vect_complex_get(twiddle3,k - 1));
 *               w4 =Conj(pnl_vect_complex_get(twiddle4,k - 1));
 *               w5 =Conj(pnl_vect_complex_get(twiddle5,k - 1));
 *             }
 *         }
 * 
 *       for (k1 = 0; k1 < p_1; k1++)
 *         {
 *           const fcomplex z0 = pnl_vect_complex_get(in,i);
 *           const fcomplex z1 = pnl_vect_complex_get(in,i+m);
 *           const fcomplex z2 = pnl_vect_complex_get(in,i+2*m);
 *           const fcomplex z3 = pnl_vect_complex_get(in,i+3*m);
 *           const fcomplex z4 = pnl_vect_complex_get(in,i+4*m);
 *           const fcomplex z5 = pnl_vect_complex_get(in,i+5*m);
 *           
 *           /\* compute x = W(6) z *\/
 * 
 *           /\* W(6) is a combination of sums and differences of W(3) acting
 *              on the even and odd elements of z *\/
 *           
 *           /\* ta1 = z2 + z4 *\/
 *           const fcomplex ta1 = Cadd(z2,z4);
 *           /\* ta2 = z0 - ta1/2 *\/
 *           const fcomplex ta2 = Csub(z0,RCmul(0.5,ta1));
 *           /\* ta3 = (+/-) sin(pi/3)*(z2 - z4) *\/
 *           const fcomplex ta3 = RCmul((int) sign * tau, (Csub(z2,z4)));
 *           /\* a0 = z0 + ta1 *\/
 *           const fcomplex a0 = Cadd(z0,ta1);
 *           /\* a1 = ta2 + i ta3 *\/
 *           const fcomplex a1 = C_op_apib(ta2, ta3);
 *           /\* a2 = ta2 - i ta3 *\/
 *           const fcomplex a2 = C_op_amib(ta2, ta3);
 *           /\* tb1 = z5 + z1 *\/
 *           const fcomplex tb1 = Cadd(z5,z1);
 *           /\* tb2 = z3 - tb1/2 *\/
 *           const fcomplex tb2 = Csub(z3,RCmul(0.5,tb1));
 *           /\* tb3 = (+/-) sin(pi/3)*(z5 - z1) *\/
 *           const fcomplex tb3 = RCmul((int) sign * tau , (Csub(z5,z1)));
 *           /\* b0 = z3 + tb1 *\/
 *           const fcomplex b0 = Cadd(z3,tb1);
 *           /\* b1 = tb2 + i tb3 *\/
 *           const fcomplex b1 = C_op_apib(tb2, tb3);
 *           /\* b2 = tb2 - i tb3 *\/
 *           const fcomplex b2 = C_op_amib(tb2,tb3);
 *           /\* x0 = a0 + b0 *\/
 *           const fcomplex x0 = Cadd(a0,b0);
 *           /\* x4 = a1 + b1 *\/
 *           const fcomplex x4 = Cadd(a1,b1);
 *           /\* x2 = a2 + b2 *\/
 *           const fcomplex x2 = Cadd(a2,b2);
 *           /\* x3 = a0 - b0 *\/
 *           const fcomplex x3 = Csub(a0,b0);
 *           /\* x1 = a1 - b1 *\/
 *           const fcomplex x1 = Csub(a1,b1);
 *           /\* x5 = a2 - b2 *\/
 *           const fcomplex x5 = Csub(a2,b2);
 *           /\* apply twiddle factors *\/
 *           
 *           /\* to0 = 1 * x0 *\/
 *           pnl_vect_complex_set(out,j, x0);
 *           /\* to1 = w1 * x1 *\/
 *           pnl_vect_complex_set(out, j + p_1, Cmul(w1,x1));
 *           /\* to2 = w2 * x2 *\/
 *           pnl_vect_complex_set(out,j + 2 * p_1, Cmul(w2,x2));
 *           /\* to3 = w3 * x3 *\/
 *           pnl_vect_complex_set(out,j + 3 * p_1, Cmul(w3,x3));
 *           /\* to4 = w4 * x4 *\/
 *           pnl_vect_complex_set(out,j + 4 * p_1, Cmul(w4,x4));
 *           /\* to5 = w5 * x5 *\/
 *           pnl_vect_complex_set(out,j + 5 * p_1, Cmul(w5,x5));
 *           i++;
 *           j++;
 *         }
 *       j += jump;
 *     }
 *   return 0;
 * }
 * 
 * 
 * static int
 * pnl_fft_complex_pass_7 (PnlVectComplex * in,
 *                         PnlVectComplex * out,
 *                         const int sign,
 *                         const uint product,
 *                         const uint n,
 *                         const PnlVectComplex * twiddle1,
 *                         const PnlVectComplex * twiddle2,
 *                         const PnlVectComplex * twiddle3,
 *                         const PnlVectComplex * twiddle4,
 *                         const PnlVectComplex * twiddle5,
 *                         const PnlVectComplex * twiddle6)
 * {
 *   uint i = 0, j = 0;
 *   uint k, k1;
 * 
 *   const uint factor = 7;
 *   const uint m = n / factor;
 *   const uint q = n / product;
 *   const uint p_1 = product / factor;
 *   const uint jump = (factor - 1) * p_1;
 * 
 *   const double c1 = cos(1.0 * 2.0 * M_PI / 7.0) ;
 *   const double c2 = cos(2.0 * 2.0 * M_PI / 7.0) ;
 *   const double c3 = cos(3.0 * 2.0 * M_PI / 7.0) ;
 * 
 *   const double s1 = sin(1.0 * 2.0 * M_PI / 7.0) ;
 *   const double s2 = sin(2.0 * 2.0 * M_PI / 7.0) ;
 *   const double s3 = sin(3.0 * 2.0 * M_PI / 7.0) ;
 * 
 *   for (k = 0; k < q; k++)
 *     {
 *       fcomplex w1, w2, w3, w4, w5,w6;
 *       if (k == 0)
 *         {
 *           w1=CONE;
 *           w2=CONE;
 *           w3=CONE;
 *           w4=CONE;
 *           w5=CONE;
 *           w6=CONE;
 *         }
 *       else
 *         {
 *           if (sign == -1)
 *             {
 *               /\* forward tranform *\/
 *               w1 =pnl_vect_complex_get(twiddle1,k - 1);
 *               w2 =pnl_vect_complex_get(twiddle2,k - 1);
 *               w3 =pnl_vect_complex_get(twiddle3,k - 1);
 *               w4 =pnl_vect_complex_get(twiddle4,k - 1);
 *               w5 =pnl_vect_complex_get(twiddle5,k - 1);
 *               w6 =pnl_vect_complex_get(twiddle6,k - 1);
 *              
 *             }
 *           else
 *             {
 *               /\* backward tranform: w -> conjugate(w) *\/
 *               w1 =Conj(pnl_vect_complex_get(twiddle1,k - 1));
 *               w2 =Conj(pnl_vect_complex_get(twiddle2,k - 1));
 *               w3 =Conj(pnl_vect_complex_get(twiddle3,k - 1));
 *               w4 =Conj(pnl_vect_complex_get(twiddle4,k - 1));
 *               w5 =Conj(pnl_vect_complex_get(twiddle5,k - 1));
 *               w6 =Conj(pnl_vect_complex_get(twiddle6,k - 1));
 *             }
 *         }
 * 
 *       for (k1 = 0; k1 < p_1; k1++)
 *         {
 *           const fcomplex z0 = pnl_vect_complex_get(in,i);
 *           const fcomplex z1 = pnl_vect_complex_get(in,i+m);
 *           const fcomplex z2 = pnl_vect_complex_get(in,i+2*m);
 *           const fcomplex z3 = pnl_vect_complex_get(in,i+3*m);
 *           const fcomplex z4 = pnl_vect_complex_get(in,i+4*m);
 *           const fcomplex z5 = pnl_vect_complex_get(in,i+5*m);
 *           const fcomplex z6 = pnl_vect_complex_get(in,i+6*m);
 *           
 *           /\* compute x = W(7) z *\/
 *           /\* t0 = z1 + z6 *\/
 *           const fcomplex t0 = Cadd(z1,z6) ;
 *           /\* t1 = z1 - z6 *\/
 *           const fcomplex t1 = Csub(z1,z6) ;
 *           /\* t2 = z2 + z5 *\/
 *           const fcomplex t2 = Cadd(z2,z5) ;
 *           /\* t3 = z2 - z5 *\/
 *           const fcomplex t3 = Csub(z2,z5) ;
 *           /\* t4 = z4 + z3 *\/
 *           const fcomplex t4 =Cadd(z4,z3) ;
 *           /\* t5 = z4 - z3 *\/
 *           const fcomplex t5 = Csub(z4,z3) ;
 *           /\* t6 = t2 + t0 *\/
 *           const fcomplex t6 = Cadd(t2, t0) ;
 *           /\* t7 = t5 + t3 *\/
 *           const fcomplex t7 = Cadd(t5,t3) ;
 *           /\* b0 = z0 + t6 + t4 *\/
 *           const fcomplex b0 = Cadd(z0,Cadd(t6, t4)) ;
 *           /\* b1 = ((cos(2pi/7) + cos(4pi/7) + cos(6pi/7))/3-1) (t6 + t4) *\/
 *           const fcomplex b1 = RCmul(((c1 + c2 + c3)/3.0 - 1.0) , Cadd(t6 , t4));
 *           /\* b2 = ((2*cos(2pi/7) - cos(4pi/7) - cos(6pi/7))/3) (t0 - t4) *\/
 *           const fcomplex b2 = RCmul(((2.0 * c1 - c2 - c3)/3.0) ,(Csub(t0,t4)));
 *           /\* b3 = ((cos(2pi/7) - 2*cos(4pi/7) + cos(6pi/7))/3) (t4 - t2) *\/
 *           const fcomplex b3 = RCmul(((c1 - 2.0*c2 + c3)/3.0) , (Csub(t4,t2)));
 *           /\* b4 = ((cos(2pi/7) + cos(4pi/7) - 2*cos(6pi/7))/3) (t2 - t0) *\/
 *           const fcomplex b4 = RCmul(((c1 + c2 - 2.0 * c3)/3.0) , (Csub(t2,t0)));
 *           /\* b5 = sign * ((sin(2pi/7) + sin(4pi/7) - sin(6pi/7))/3) (t7 + t1) *\/
 *           const fcomplex b5 = RCmul(-(int)sign * ((s1 + s2 - s3)/3.0) , Cadd(t7 , t1)) ;
 *           /\* b6 = sign * ((2sin(2pi/7) - sin(4pi/7) + sin(6pi/7))/3) (t1 - t5) *\/
 *           const fcomplex b6 = RCmul(-(int)sign * ((2.0 * s1 - s2 + s3)/3.0) , (Csub(t1,t5))) ;
 *           /\* b7 = sign * ((sin(2pi/7) - 2sin(4pi/7) - sin(6pi/7))/3) (Csub(t5,t3)) *\/
 *           const fcomplex b7 = RCmul(-(int)sign * ((s1 - 2.0 * s2 - s3)/3.0) , (Csub(t5,t3))) ;
 *           /\* b8 = sign * ((sin(2pi/7) + sin(4pi/7) + 2sin(6pi/7))/3) (Csub(t3,t1)) *\/
 *           const fcomplex b8 = RCmul(-(int)sign * ((s1 + s2 + 2.0 * s3)/3.0) , (Csub(t3,t1))) ;
 *           
 *           /\* T0 = b0 + b1 *\/
 *           const fcomplex T0 = Cadd(b0 , b1) ;
 *           /\* T1 = b2 + b3 *\/
 *           const fcomplex T1 = Cadd(b2 , b3) ;
 *           /\* T2 = b4 - b3 *\/
 *           const fcomplex T2 = Csub(b4,b3) ;
 *           /\* T3 = -b2 - b4 *\/
 *           const fcomplex T3 = RCmul(-1.,Cadd(b2,b4));
 *           /\* T4 = b6 + b7 *\/
 *           const fcomplex T4 = Cadd(b6 , b7) ;
 *           /\* T5 = b8 - b7 *\/
 *           const fcomplex T5 = Csub(b8,b7) ;
 *           /\* T6 = -b8 - b6 *\/
 *           const fcomplex T6 = RCmul(-1,Cadd(b8,b6)) ;
 *           /\* T7 = T0 + T1 *\/
 *           const fcomplex T7 = Cadd(T0 , T1) ;
 *           /\* T8 = T0 + T2 *\/
 *           const fcomplex T8 = Cadd(T0 , T2) ;
 *           /\* T9 = T0 + T3 *\/
 *           const fcomplex T9 = Cadd(T0 , T3) ;
 *           /\* T10 = T4 + b5 *\/
 *           const fcomplex T10 = Cadd(T4 , b5) ;
 *           /\* T11 = T5 + b5 *\/
 *           const fcomplex T11 = Cadd(T5 , b5) ;
 *           /\* T12 = T6 + b5 *\/
 *           const fcomplex T12 = Cadd(T6 , b5) ;
 *           /\* x0 = b0 *\/
 *           const fcomplex x0 = b0 ;
 *           /\* x1 = T7 - i T10 *\/
 *           const fcomplex x1 = C_op_amib(T7,T10);
 *           /\* x2 = T9 - i T12 *\/
 *           const fcomplex x2 = C_op_amib(T9, T12) ;
 *           /\* x3 = T8 + i T11 *\/
 *           const fcomplex x3 = C_op_apib(T8,T11) ;
 *           /\* x4 = T8 - i T11 *\/
 *           const fcomplex x4 = C_op_amib(T8, T11) ;
 *           /\* x5 = T9 + i T12 *\/
 *           const fcomplex x5 = C_op_apib(T9,T12) ;
 *           /\* x6 = T7 + i T10 *\/
 *           const fcomplex x6 = C_op_apib(T7, T10) ;
 *          
 *           
 *           /\* apply twiddle factors *\/
 *           /\* to0 = 1 * x0 *\/
 *           pnl_vect_complex_set(out,j, x0);
 *           /\* to1 = w1 * x1 *\/
 *           pnl_vect_complex_set(out, j + p_1, Cmul(w1,x1));
 *           /\* to2 = w2 * x2 *\/
 *           pnl_vect_complex_set(out,j + 2 * p_1, Cmul(w2,x2));
 *           /\* to3 = w3 * x3 *\/
 *           pnl_vect_complex_set(out,j + 3 * p_1, Cmul(w3,x3));
 *           /\* to4 = w4 * x4 *\/
 *           pnl_vect_complex_set(out,j + 4 * p_1, Cmul(w4,x4));
 *           /\* to5 = w5 * x5 *\/
 *           pnl_vect_complex_set(out,j + 5 * p_1, Cmul(w5,x5));
 *           i++; j++;
 *           /\* to6 = w6 * x6 *\/
 *           pnl_vect_complex_set(out,j + 6 * p_1, Cmul(w6,x6));
 *           i++; j++;
 *         }
 *       j += jump;
 *     }
 *   return 0;
 * }
 * 
 * 
 * static int pnl_fft_complex_pass_n (PnlVectComplex * in,
 *                                    PnlVectComplex * out,
 *                                    const int sign,
 *                                    const uint factor,
 *                                    const uint product,
 *                                    const uint n,
 *                                    const PnlVectComplex * twiddle)
 * {
 *   uint i = 0, j = 0;
 *   uint k, k1;
 * 
 *   double bp, bm;
 * 
 *   const uint m = n / factor;
 *   const uint q = n / product;
 *   const uint p_1 = product / factor;
 *   const uint jump = (factor - 1) * p_1;
 * 
 *   uint e, e1;
 * 
 *   for (i = 0; i < m; i++)
 *     {
 *       pnl_vect_complex_set(out,i, pnl_vect_complex_get(in,i));
 *     }
 *   for (e = 1; e < (factor - 1) / 2 + 1; e++)
 *     {
 *       for (i = 0; i < m; i++)
 *         {
 *           const uint idx = i + e * m;
 *           const uint idxc = i + (factor - e) * m;
 *           pnl_vect_complex_set(out,idx, Cadd(pnl_vect_complex_get(in,idx),pnl_vect_complex_get(in,idxc)));
 *           pnl_vect_complex_set(out,idxc, Csub(pnl_vect_complex_get(in,idx),pnl_vect_complex_get(in,idxc)));
 *         }
 *     }
 * 
 *   /\* e = 0 *\/
 * 
 *   for (i=0 ; i<m; i++) 
 *     {
 *       pnl_vect_complex_set(in,i,pnl_vect_complex_get(out,i));
 *     }
 * 
 *   for (e1 = 1; e1 < (factor - 1) / 2 + 1; e1++)
 *     {
 *       for (i = 0; i < m; i++)
 *         {
 *           pnl_vect_complex_set(in,i, Cadd(pnl_vect_complex_get(in,i),pnl_vect_complex_get(out,i + e1*m))) ;
 *         }
 *     }
 * 
 *   for (e = 1; e < (factor-1)/2 + 1; e++)
 *     {
 *       uint idx = e*q ;
 *       const uint idx_step = e * q ;
 *       fcomplex w ;
 * 
 *       const uint em = e * m ;
 *       const uint ecm = (factor - e) * m ;
 * 
 *       for (i = 0; i < m; i++) 
 *         {
 *           pnl_vect_complex_set(in,i+em, pnl_vect_complex_get(out,i) );
 *           pnl_vect_complex_set(in,i+ecm, pnl_vect_complex_get(out,i)) ;
 *         }
 * 
 *       for (e1 = 1; e1 < (factor - 1) / 2 + 1; e1++)
 *         {
 *           if (idx == 0) {
 *             w = CONE ;
 *           }
 *           else {
 *             if (sign == -1) {
 *               w = pnl_vect_complex_get(twiddle,idx - 1) ;
 *             }
 *             else {
 *               w = Conj(pnl_vect_complex_get(twiddle,idx - 1)) ;
 *             }
 *           }
 * 
 *           for (i = 0; i < m; i++) 
 *             {
 *               const fcomplex xp = pnl_vect_complex_get(out,i + e1 * m);
 *               const fcomplex xm = pnl_vect_complex_get(out,i + (factor - e1) *m);
 *               const double ap = Creal(w) * Creal(xp) ;
 *               const double am = Cimag(w) * Cimag(xm) ; 
 * 
 *               fcomplex sum,sumc;
 *               sum.r= ap - am;
 *               sumc.r = ap + am;
 * 
 *               bp = Creal(w) * Cimag(xp) ;
 *               bm = Cimag(w) * Creal(xm) ;
 * 
 *               sum.i = bp + bm;
 *               sumc.i = bp - bm;
 * 
 *               pnl_vect_complex_set(in,i + em, Cadd(pnl_vect_complex_get(in,i + em),sum));
 *               pnl_vect_complex_set(in,i + ecm, Cadd(pnl_vect_complex_get(in,i + ecm),sumc));
 *             }
 *           idx += idx_step ;
 *           idx %= factor * q ;
 *         }
 *     }
 * 
 *   i = 0;
 *   j = 0;
 * 
 *   /\* k = 0 *\/
 *   for (k1 = 0; k1 < p_1; k1++)
 *     {
 *       pnl_vect_complex_set(out,k1, pnl_vect_complex_get(in,k1));
 *     }
 * 
 *   for (e1 = 1; e1 < factor; e1++)
 *     {
 *       for (k1 = 0; k1 < p_1; k1++)
 *         {
 *           pnl_vect_complex_set(out,k1 + e1 * p_1, pnl_vect_complex_get(in,k1 + e1 * m)) ;
 *         }
 *     }
 * 
 *   i = p_1 ;
 *   j = product ;
 * 
 *   for (k = 1; k < q; k++)
 *     {
 *       for (k1 = 0; k1 < p_1; k1++)
 *         {
 *           pnl_vect_complex_set(out,j, pnl_vect_complex_get(in,i));
 *           i++;
 *           j++;
 *         }
 *       j += jump;
 *     }
 * 
 *   i = p_1 ;
 *   j = product ;
 * 
 *   for (k = 1; k < q; k++)
 *     {
 *       for (k1 = 0; k1 < p_1; k1++)
 *         {
 *           for (e1 = 1; e1 < factor; e1++)
 *             {
 *               fcomplex x = pnl_vect_complex_get(in,i + e1 * m);
 *               fcomplex w;
 *               if (sign == -1) {
 *                 w = pnl_vect_complex_get(twiddle,(e1-1)*q + k-1) ;
 *                 
 *               }
 *               else {
 *                 w = Conj(pnl_vect_complex_get(twiddle,(e1-1)*q + k-1)) ;
 *               }
 * 
 *               pnl_vect_complex_set(out,j + e1 * p_1, Cmul(w, x));
 *             }
 *           i++;
 *           j++;
 *         }
 *       j += jump;
 *     }
 * 
 *   return 0;
 * } */
