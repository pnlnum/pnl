
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "config.h"
#include "pnl_specfun.h"

#ifdef HAVE_GSL
extern double gsl_sf_expint_Ei (double x);
extern double gsl_sf_expint_En (int n,double x);
extern double gsl_sf_gamma_inc(double a,double x);
extern double gsl_sf_gamma_inc_Q (double a,double x);
extern double gsl_sf_gamma_inc_P (double a,double x);
#endif


static void exp_int_test ()
{
#ifdef HAVE_GSL
  int      n;
  double Gamma_gsl,IP_gsl,IQ_gsl,En_gsl,Ei_gsl,Gamma_neg_gsl;
  double a,x;
  double Gamma,IP,IQ,En,Ei,Gamma_neg;

  a=0.2;
  x=5.0;
  
  pnl_gamma_inc(a,x,&Gamma,&IP,&IQ);

  Gamma_neg=pnl_gamma_inc_func(-a,x);
  Ei=pnl_expint_Ei(x);

  printf( "error in Exponential Integrals\n");

  Gamma_gsl=gsl_sf_gamma_inc(a,0);
  Gamma_neg_gsl=gsl_sf_gamma_inc(-a,x);
  IQ_gsl=gsl_sf_gamma_inc_Q (a,x);
  IP_gsl=gsl_sf_gamma_inc_P (a,x);
  Ei_gsl=gsl_sf_expint_Ei (x);

  for(n=1;n<10;n++)
    {
      En=pnl_expint_En(n,x);
      En_gsl=gsl_sf_expint_En (n,x);
      printf( "  E_%d      = %7.4f - %7.4f = %7.4f \n",n,En,En_gsl,En-En_gsl);
    }
  printf( "Gamma     = %7.4f - %7.4f = %7.4f \n",Gamma,Gamma_gsl,Gamma-Gamma_gsl);
  printf( "Gamma_neg = %7.4f - %7.4f = %7.4f \n",Gamma_neg,Gamma_neg_gsl,Gamma_neg-Gamma_neg_gsl);
  printf( "  IQ      = %7.4f - %7.4f = %7.4f \n",IQ,IQ_gsl,IQ-IQ_gsl);
  printf( "  IP      = %7.4f - %7.4f = %7.4f \n",IP,IP_gsl,IP-IP_gsl);
  printf( "  Ei      = %7.4f - %7.4f = %7.4f \n",Ei,Ei_gsl,Ei-Ei_gsl);
  printf( "  En      = %7.4f - %7.4f = %7.4f \n",En,En_gsl,En-En_gsl);
#else
  printf("Tests for Exponential Integrals only available with GSL\n");
#endif
}

void special_func_test ()
{
  exp_int_test();
 
}
