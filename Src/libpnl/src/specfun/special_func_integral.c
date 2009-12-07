#include <math.h>
#include "pnl_mathtools.h"
#include "pnl_specfun.h"

/*************************************************************************/
/* Written and (C) by David Pommier <pommier.david@gmail.com>            */
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

#define SPECFUNC_EPS 3.0e-7
/*Relative accurancy for incomplete gamma function  */ 
#define SPECFUNC_ITMAX 100
/*Maximum number of iterations for incomplete gamma function  */ 
#define SPECFUNC_FPMIN 1.0e-30 
/*Number near the smallest representable floatting-point number  for incomplete gamma function*/

static void gamma_inc_compute_P(double * res,double a,double x,double lnga)
{
     
  int n;
  double ap,del,sum;
  if (x == 0.0)
    *res=0.0;
  else{
    ap=a;
    del=sum=1.0/a;
    for(n=1;n<=SPECFUNC_ITMAX;n++)
      {
        ++ap;
        del*=x/ap;
        sum+=del;
        if(fabs(del)<fabs(sum)*SPECFUNC_EPS)
          {
            *res=sum*exp(-x+a*log(x)-lnga);
            return;
          };
      }
    PNL_ERROR("Not comverged sum in ", "gamma_inc_compute");
  }
}

static void gamma_inc_compute_Q(double * res,double a,double x,double lnga)
{
     
  int n;
  double b,c,d,h,del,an;
  b=x+1.0-a;
  c=1/SPECFUNC_FPMIN;
  d=1.0/b;
  h=d;
  for(n=1;n<=SPECFUNC_ITMAX;n++)
    {
      an=-n*(n-a);
      b+=2;
      d=an*d+b;
      if(fabs(d)<SPECFUNC_FPMIN) d=SPECFUNC_FPMIN;
      c=b+an/c;
      if(fabs(c)<SPECFUNC_FPMIN) c=SPECFUNC_FPMIN;
      d=1/d;
      del=d*c;
      h*=del;
      if(fabs(del-1)<SPECFUNC_EPS) 
        {
          *res=exp(-x+a*log(x)-lnga)*h;
          return;
        }
    }
  PNL_ERROR("Not comverged sum in ", "gamma_inc_compute");
}


/** Incomplete Gamma Function factorial function
 *  Compute 
 *  Gamma(a,x) = int_x^{infty} e^{-u} u^{a-1} du 
 *
 *
 * @param a double
 * @param x double
 * @return Gamma(a,x)
 */
double pnl_gamma_inc_func(double a,double x)
{
  if (x==0) return tgamma (a);
  if (a==0) return pnl_expint_En(1,x);
  if (a<0)
    return (-pow(x,a)*exp(-x) + pnl_gamma_inc_func(a+1,x))/a;
  else
    {
      double res,lgamma_a=lgamma(a);
      if (x < 1.0+a){
        gamma_inc_compute_P(&res,a,x,lgamma_a);
        /* Use the series representation*/
        return (1-res)*exp(lgamma_a);
      } 
      else{
        gamma_inc_compute_Q(&res,a,x,lgamma_a);
        return res*exp(lgamma_a);
      }
    }
}


/** Incomplete Gamma Function factorial function
 *  Compute 
 *  G(a)   = Gamma(a)
 *  P(a,x) = 1/Gamma(a) int_0^x u^{a-1} e^{-u}  du 
 *  Q(a,x) = 1-P(a,x) = 1/Gamma(a) int_x^{infty} e^{-u} u^{a-1} du 
 *
 *
 * @param a double
 * @param x double
 * @param Ga a pointer on double
 * @param P a pointer on double
 * @param Q a pointer on double
 */
void pnl_gamma_inc(double a,double x,double * Ga,double *P,double *Q)
{
  double lgamma_a,res;
  res = 0.; /* to avoid a warning */
  if ((a <= 0.0) || (x < 0))
    PNL_ERROR("Invalid arguments", "gamma_inc");
 
  lgamma_a=lgamma(a);
  *Ga=exp(lgamma_a);
  if (x < 1.0+a){
    gamma_inc_compute_P(&res,a,x,lgamma_a);
    /* Use the series representation*/
    *P=res;
    *Q=1-res;
  } 
  else{
    gamma_inc_compute_Q(&res,a,x,lgamma_a);
    /* Use the continued fraction representation*/
    *P=1-res;
    *Q=res;
  }
}

/** Exponential integral is
 * = int_1^{infty} u^{-n} exp(-x u) du 
 * = x^{n-1} int_x^{infty} u^{-n} exp(-u) du 
 * = x^{n-1} Gamma(1-n,x)
 *
 *
 * @param n double
 * @param x double
 * @return int_1^{infty} u^{-n} exp(-x u) du
 */
double pnl_expint_En(int n,double x)
{
  int i,j,nm1;
  double a,b,c,d,del,fact,h,psi,ans;
  nm1=n-1;
  if(n< 0 || x <0.0 || (x==0.0 && (n==0.0 || n==1)))
    {
      PNL_ERROR("Invalid arguments", "expint_En");
      return 0.0;
    }
  if(n==0) return exp(-x)/x;
  if(x==0) return 1.0/nm1;
  if(x>1)/* Lentz's algotithm */
    {
      b=x+n;
      c=1/SPECFUNC_FPMIN;
      d=1.0/b;
      h=d;
      for(i=1;i<=SPECFUNC_ITMAX;i++)
        {
          a=-i*(nm1+i);
          b+=2;
          d=1.0/(a*d+b);
          c=b+a/c;
          del=d*c;
          h*=del;
          if(fabs(del-1)<SPECFUNC_EPS) 
            return h*exp(-x);
        }
      PNL_ERROR("Continued fraction failed ", "expint_En");
      return 0.0;
    }
  else/* Evaluated series */
    {
      ans=(nm1!=0)?1.0/nm1:-log(x)-M_EULER;
      fact=1;
      for(i=1;i<=SPECFUNC_ITMAX;i++)
        {
          fact*=-x/i;
          if(i!=nm1) del=-fact/(i-nm1);
          else
            {
              psi=-M_EULER;
              for(j=1;j<=nm1;j++) psi+=1.0/j;
              del=fact*(-log(x)+psi);
            }
          ans+=del;
          if(fabs(del)<fabs(ans)*SPECFUNC_EPS) 
            return ans;
        }
      PNL_ERROR("Serie expansion failed ", "expint_En");
      return 0.0;
    }
}


/** Principal value of the integral is
 * = - int_{-x}^{infty} exp(-u)/u du 
 * = int_{-infty}^x exp(u)/u du  x>0
 *
 * @param x double
 * @return - int_{-x}^{infty} exp(-u)/u du x>0
 */
double pnl_expint_Ei(double x)
{
  int i;
  double fact,prev,sum,term;
  if(x <= 0.0)
    {
      PNL_ERROR("Invalid arguents", "expint_Ei");
      return 0.0;
    }
  if(x<=SPECFUNC_FPMIN) return log(x)+M_EULER;
  if(x<=-log(SPECFUNC_EPS))
    {
      sum=0.0;
      fact=1.0;
      for(i=1;i<=SPECFUNC_ITMAX;i++)
        {
          fact*=x/i;
          term=fact/i;
          sum+=term;
          if(term<sum*SPECFUNC_EPS)
            return sum+log(x)+M_EULER;
        }
      PNL_ERROR("Series failed ", "expint_Ei");
      return 0.0;
    }

  else
    {
      sum=0.0;
      term=1;
      for(i=1;i<=SPECFUNC_ITMAX;i++)
        {
          prev=term;
          term*=i/x;
          if(term<SPECFUNC_EPS) break;
          if(term<prev) sum+=term;
          else
            {
              sum-=prev;
              break;
            }
        }
      return exp(x)*(1.0+sum)/x;
    }
}


#ifdef SPECFUNC_EPS
#undef SPECFUNC_EPS
#endif
#ifdef SPECFUNC_ITMAX
#undef SPECFUNC_ITMAX
#endif
#ifdef SPECFUNC_FPMIN
#undef SPECFUNC_FPMIN
#endif
