
/*************************************************************************/
/* Written and (C) by David Pommier <david.pommier@gmail.com>            */
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

#include "pnl_vector.h"
#include "pnl_matrix.h"
#include "pnl_matrix_uint.h"
#include "pnl_cdf.h"
#include "pnl_finance.h"

/* Old version of premia close formula  */
int CF_Call_BS(double s,double k,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta){
  double sigmasqrt,d1,d2,delta;
  
  sigmasqrt=sigma*sqrt(t);
  d1=(log(s/k)+(r-divid)*t)/sigmasqrt+sigmasqrt/2.;
  d2=d1-sigmasqrt;
  delta=exp(-divid*t)*cdf_nor(d1);
  /*Price*/
  *ptprice=s*delta-exp(-r*t)*k*cdf_nor(d2);
  /*Delta*/
  *ptdelta=delta;
  return OK;
}

int CF_Put_BS(double s,double k,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta)
{
  double sigmasqrt,d1,d2,delta;

  sigmasqrt=sigma*sqrt(t);
  d1=(log(s/k)+(r-divid)*t)/sigmasqrt+sigmasqrt/2.;
  d2=d1-sigmasqrt;
  delta=-exp(-divid*t)*cdf_nor(-d1);

  /*Price*/
  *ptprice=exp(-r*t)*k*cdf_nor(-d2)+delta*s;

  /*Delta*/
  *ptdelta=delta;
  return OK;
}

/* Compariason with Old version of premia close formula  */
int CF_Call_BS_Comp(double s,double k,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta)
{
  double sigmasqrt,d1,d2,delta;
  double bond=exp(-r*t);
  double forward=s*exp((r-divid)*t);
  sigmasqrt=sigma*sqrt(t);
  d1=(log(s/k)+(r-divid)*t)/sigmasqrt+sigmasqrt/2.;
  d2=d1-sigmasqrt;
  delta=exp(-divid*t)*cdf_nor(d1);
  /*Price*/
  *ptprice=s*delta-exp(-r*t)*k*cdf_nor(d2)-pnl_bs_call(sigma,bond,forward,k,t);
  /*Delta*/
  *ptdelta=delta-forward/s*pnl_bs_call_delta_forward(sigma,bond,forward,k,t);
  printf(" error in CF_Call_BS_Comp diff_price = %7.4f diff_delta = %7.4f \n",*ptprice,*ptdelta);
  return OK;
}

int CF_Put_BS_Comp(double s,double k,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta)
{
  double sigmasqrt,d1,d2,delta;
  double bond=exp(-r*t);
  double forward=s*exp((r-divid)*t);
  sigmasqrt=sigma*sqrt(t);
  d1=(log(s/k)+(r-divid)*t)/sigmasqrt+sigmasqrt/2.;
  d2=d1-sigmasqrt;
  delta=-exp(-divid*t)*cdf_nor(-d1);
  /*Price*/
  *ptprice=exp(-r*t)*k*cdf_nor(-d2)+delta*s-pnl_bs_put(sigma,bond,forward,k,t);
  /*Delta*/
  *ptdelta=delta-forward/s*pnl_bs_put_delta_forward(sigma,bond,forward,k,t);
  printf(" error in CF_Call_BS_Comp diff_price = %7.4f diff_delta = %7.4f \n",*ptprice,*ptdelta);
  return OK;
}



static void test_pnl_finance_function_call_put()
{
  double s=100;
  double k=98;
  double t=3;
  double r=0.05;
  double divid=0.01;
  double sigma=0.2;
  double sigma2=0;
  double ptprice;
  double ptdelta;
  double bond=exp(-r*t);
  double forward=s*exp((r-divid)*t);
  double Price;
  CF_Call_BS_Comp(s,k,t,r,divid,sigma,&ptprice,&ptdelta);
  CF_Put_BS_Comp(s,k,t,r,divid,sigma,&ptprice,&ptdelta);
  printf(" implied volatility expected = %7.4f ",sigma);
  Price= pnl_bs_call(sigma,bond,forward,k,t);
  sigma2=pnl_bs_implicit_vol (1,Price,bond, forward,k,t);
  printf(" computed = %7.4f \n",sigma2);

}

static void test_pnl_finance_function_vol_impli()
{
  int i,j,N=10;
  double s=100;
  double k0=100;
  double tmax=10;
  double r=0.05;
  double divid=0.01;
  double sigma0=0.2;
  PnlVect * Strike;
  PnlVect * Matu;
  PnlMatUint * IsCall;
  PnlMat * Price;
  PnlMat * VolImpli;
  PnlMat * VolImpli2;
  
  Strike=pnl_vect_create_from_double(N,k0);
  Matu=pnl_vect_create_from_double(N,tmax);
  IsCall=pnl_mat_uint_create(N,N);
  Price=pnl_mat_create_from_double(N,N,0.0);
  VolImpli=pnl_mat_create_from_double(N,N,0.0);
  VolImpli2=pnl_mat_create_from_double(N,N,0.0);
  for(i=0;i<Strike->size;i++)  
    LET(Strike,i)-=(-N/2+i)*k0/N;
  for(j=0;j<Matu->size;j++)
    LET(Matu,j)*=(double)(j+1)/(double)(N);
  for(i=0;i<Strike->size;i++)  
    for(j=0;j<Matu->size;j++)
      {
        double bond, forward;
        *(pnl_mat_uint_lget(IsCall,i,j))=(GET(Strike,i)<s)?1:0;
        MLET(VolImpli,i,j)=sigma0*(1+(GET(Strike,i)-s)/s) *exp(-0.1*GET(Matu,j));
        bond=exp(-r*GET(Matu,j));
        forward=s*exp((r-divid)*GET(Matu,j));
        MLET(Price,i,j)=pnl_bs_call_put (pnl_mat_uint_get(IsCall,i,j),MGET(VolImpli,i,j),bond,forward,
                                         GET(Strike,i),GET(Matu,j));
      }
  pnl_bs_matrix_implicit_vol(IsCall,Price,s,r,divid,Strike,Matu,VolImpli2);
  pnl_mat_print(VolImpli2);
  printf("\n");
  pnl_mat_print(VolImpli);
  
  pnl_mat_free(&VolImpli2);
  pnl_mat_free(&VolImpli);
  pnl_mat_free(&Price);
  pnl_mat_uint_free(&IsCall);
  pnl_vect_free(&Matu);
  pnl_vect_free(&Strike);
}


void finance_function_test(void)
{
  test_pnl_finance_function_call_put();
  test_pnl_finance_function_vol_impli();
}
