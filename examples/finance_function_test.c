
/************************************************************************/
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

#include "pnl/pnl_cdf.h"
#include "pnl/pnl_finance.h"
#include "tests_utils.h"

/* Old version of premia closed formula  */
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

static void test_pnl_finance_function_call_put()
{
  double s=100;
  double k=98;
  double t=3;
  double r=0.05;
  double divid=0.01;
  double sigma=0.2;
  double ptprice;
  double ptdelta;
  int error;
  CF_Put_BS(s,k,t,r,divid,sigma,&ptprice,&ptdelta);
  pnl_test_eq (pnl_bs_put(s,k,t,r,divid,sigma), ptprice, 1E-12, "pnl_bs_put", "");
  pnl_test_eq (pnl_bs_implicit_vol (0, ptprice, s, k, t, r, divid, &error), sigma,
               1E-12, "pnl_bs_implicit_vol", "");

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
  PnlVect *Strike;
  PnlVect *Matu;
  PnlMatInt *IsCall;
  PnlMat *Price;
  PnlMat *Vol;
  PnlMat *VolImpli;
  
  Strike = pnl_vect_create_from_scalar(N, k0);
  Matu = pnl_vect_create_from_scalar(N, tmax);
  IsCall = pnl_mat_int_create(N, N);
  Price = pnl_mat_create_from_scalar(N, N, 0.0);
  Vol = pnl_mat_create_from_scalar(N, N, 0.0);
  VolImpli = pnl_mat_create_from_scalar(N, N, 0.0);
  for(i=0;i<Strike->size;i++)  LET(Strike,i)-=(-N/2+i)*k0/N;
  for(j=0;j<Matu->size;j++) LET(Matu,j)*=(double)(j+1)/(double)(N);
  for(i=0;i<Strike->size;i++)  
    {
    for(j=0;j<Matu->size;j++)
      {
        *(pnl_mat_int_lget(IsCall, i, j)) = (GET(Strike, i)<s)?1:0;
        MLET(Vol, i, j) = sigma0*(1+(GET(Strike, i)-s)/s) *exp(-0.1*GET(Matu, j));
        MLET(Price, i, j) = pnl_bs_call_put (pnl_mat_int_get(IsCall, i, j), s, 
                                             GET(Strike,i),GET(Matu,j), 
                                             r, divid,MGET(Vol,i,j));
      }
    }
  pnl_bs_matrix_implicit_vol(IsCall,Price,s,r,divid,Strike,Matu,VolImpli);
  pnl_test_mat_eq (VolImpli,Vol, 1E-5, "pnl_bs_matrix_implicit_vol", "");
  
  pnl_mat_free(&VolImpli);
  pnl_mat_free(&Vol);
  pnl_mat_free(&Price);
  pnl_mat_int_free(&IsCall);
  pnl_vect_free(&Matu);
  pnl_vect_free(&Strike);
}


int main (int argc, char **argv)
{
  pnl_test_init (argc, argv);
  test_pnl_finance_function_call_put();
  test_pnl_finance_function_vol_impli();
  exit(pnl_test_finalize("Financial functions"));
}
