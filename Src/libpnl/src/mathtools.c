#include <math.h>
#include "pnl_mathtools.h"
#include "pnl_cdf.h"


/** nearest integer round off function
 * @param s the value to be rounded
 * @return the rounded value casted to an integer. 1.5 is rounded to 1 and
 * -1.5 is rounded to -1
 */
int intapprox(double s)
{
  int r = (int) s;
  if(s>r+0.5){r=r+1;}
  if(s<r-0.5){r=r-1;}
  return r;
}


#ifndef HAVE_TRUC
double trunc(double x)
{
  if(x >= 0) return floor(x);
  return ceil(x);
}
#endif

/** Computes the binomial coefficients with an error smaller than 1 for large
 * values of the parameter n
 * @param n an integer
 * @param p an integer
 * @return C(n, p)
 */
double Cnp(int n, int p) 
{                          
  
  double z, iter;
  int i;
  z=0.0;

  if ((n-p<= -1) || (n<0) || (p<0)){
    return z;
  }
  else{
    if (p==0)
      z=1.0;
    else{
      z=1.0;
      iter=z;
      i=0;
      while(i<=p-1){
        iter=iter*(double)(n-i)/(p-i);
        i=i+1;
      }
      if ((iter-floor(iter))>=0.5)
        z=floor(iter)+1;
      else
        z=floor(iter);
    }
  }

  return z;
}

/** factorial function
 * @param n an integer
 * @return the factorial of n as a double
 */
double pnl_fact(int n)  
{
  int i;
  double z=1.0;
  if (n<=0) { z=1.0; }
  else { for ( i=1 ; i<=n ; i++) z=z*i; }
  return z;
}


#ifndef HAVE_LGAMMA
/* (C) Copr. 1986-92 Numerical Recipes Software A2.>$Y0%9j. */
double lgamma(double xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
                        24.01409824083091,-1.231739572450155,
                        0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
   
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
        
  return -tmp+log(2.5066282746310005*ser/x);
}
#endif

/* the tgamma function is part of C99, not available under
 * Windows */
#ifndef HAVE_TGAMMA
double tgamma(double x)
{
  return (x > 0) ? exp(lgamma(x)) : M_PI/tgamma(1-x)/sin(M_PI*(1-x));
} 

#endif


