/*************************************************************************/
/* Written and (C) by Céline Labart <labart@cmap.polytechnique.fr        */
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
#include "pnl_cdf.h"

/*test des fonctions pnl_cdf_ */


static void pnl_cdf_bet_test()
{
  int which;
  double p;
  double q;
  double x;
  double y;
  double a;
  double b;
  int status;
  double bound;
  printf("test de la fonction 'pnl_cdf_bet' : ");
  which=1;
  x=0.6;
  y=0.4;
  a=2.0;
  b=3.0;
  pnl_cdf_bet(&which,&p,&q,&x,&y,&a,&b,&status,&bound);
  printf("p=%f ",p);
  printf("q=%f \n",q);
  if (status != 0) printf("status=%d \n",status); 
}

static void pnl_cdf_bin_test()
{
  int which;
  double p;
  double q;
  double s;
  double xn;
  double pr;
  double ompr;
  int status;
  double bound;
  printf("test de la fonction 'pnl_cdf_bin' : ");
  which=1;
  s=2;
  xn=5;
  pr=0.7;
  ompr=0.3;
  pnl_cdf_bin(&which,&p,&q,&s,&xn,&pr,&ompr,&status,&bound);
  printf("p=%f ",p);
  printf("q=%f\n",q);
  if (status != 0) printf("status=%d \n",status); 
}

static void pnl_cdf_chi_test()
{
  int which;
  double p;
  double q;
  double x;
  double df;
  int status;
  double bound;
  printf("test de la fonction 'pnl_cdf_chi' : ");
  which=3;
  p=0.7;
  q=0.3;
  x=2;
  pnl_cdf_chi(&which,&p,&q,&x,&df,&status,&bound);
  printf("df=%f ",df); /* 1.6862067645 */
  if (status != 0) printf("status=%d \n",status);
  which = 2;
  df = 6;
  pnl_cdf_chi(&which,&p,&q,&x,&df,&status,&bound);
  printf("x=%f \n",x); /* 7.2311353317 */
  if (status != 0) printf("status=%d \n",status);
}

static void pnl_cdf_chn_test()
{
  int which;
  double p;
  double q;
  double x;
  double df;
  double pnonc;
  int status;
  double bound;
  printf("test de la fonction 'pnl_cdf_chn' : ");
  which=4;
  p=0.1;
  q=0.9;
  x=2;
  df=2;
  pnl_cdf_chn(&which,&p,&q,&x,&df,&pnonc,&status,&bound);
  printf("pnonc=%f \n",pnonc);
  if (status != 0) printf("status=%d \n",status);
}

static void pnl_cdf_f_test()
{
  int which;
  double p;
  double q;
  double f;
  double dfn;
  double dfd;
  int status;
  double bound;
  printf("test de la fonction 'pnl_cdf_f' : ");
  which=4;
  p=0.1;
  q=0.9;
  f=2;
  dfn=3;
  pnl_cdf_f(&which,&p,&q,&f,&dfn,&dfd,&status,&bound);
  printf("dfd=%f \n",dfd);  /* 05030479580 */
  if (status != 0) printf("status=%d \n",status);
}

static void pnl_cdf_fnc_test()
{
  int which;
  double p;
  double q;
  double f;
  double dfn;
  double dfd;
  double pnonc;
  int status;
  double bound;
  printf("test de la fonction 'pnl_cdf_fnc' : ");
  which=5;
  p=0.1;
  q=0.9;
  f=2;
  dfn=3;
  dfd=5;
  pnl_cdf_fnc(&which,&p,&q,&f,&dfn,&dfd,&pnonc,&status,&bound);
  printf("pnonc=%f \n",pnonc); /* 13.1787447557 */
  if (status != 0) printf("status=%d \n",status); 
}

static void pnl_cdf_gam_test()
{
  int which;
  double p;
  double q;
  double x;
  double shape;
  double scale;
  int status;
  double bound;
  printf("test de la fonction 'pnl_cdf_gam' : ");
  which=4;
  p=0.1;
  q=0.9;
  x=2;
  shape=3;
  pnl_cdf_gam(&which,&p,&q,&x,&shape,&scale,&status,&bound);
  printf("scale=%f \n",scale); /* 0.5510326641 */
  if (status != 0) printf("status=%d \n",status); 
}



static void pnl_cdf_nbn_test()
{
  int which;
  double p;
  double q;
  double s;
  double xn;
  double pr;
  double ompr;
  int status;
  double bound;
  printf("test de la fonction 'pnl_cdf_nbn' : ");
  which=1;
  s=2.0;
  xn=3.0;
  pr=0.7;
  ompr=0.3;
  pnl_cdf_nbn(&which,&p,&q,&s,&xn,&pr,&ompr,&status,&bound);
  printf("p=%f ",p);
  printf("q=%f \n",q);
  if (status != 0) printf("status=%d \n",status); 
}

static void pnl_cdf_nor_test()
{
  int which;
  double p;
  double q;
  double x;
  double mean;
  double sd;
  int status;
  double bound;
  printf("test de la fonction 'pnl_cdf_nor' : ");
  which=1;
  x=2.0;
  mean=0.0;
  sd=1.0;
  pnl_cdf_nor(&which,&p,&q,&x,&mean,&sd,&status,&bound);
  printf("p=%f ",p);
  printf("q=%f \n",q);
  if (status != 0) printf("status=%d \n",status); 
}

static void pnl_cdf_poi_test()
{
  int which;
  double p;
  double q;
  double s;
  double xlam;
  int status;
  double bound;
  printf("test de la fonction 'pnl_cdf_poi' : ");
  which=3;
  p=0.4;
  q=0.6;
  s=5.0;
  pnl_cdf_poi(&which,&p,&q,&s,&xlam,&status,&bound);
  printf("xlam=%f \n",xlam);
  if (status != 0) printf("status=%d \n",status); 
}

static void pnl_cdf_t_test()
{
  int which;
  double p;
  double q;
  double t;
  double df;
  int status;
  double bound;
  printf("test de la fonction 'pnl_cdf_t' : ");
  which=3;
  p=0.4;
  q=0.6;
  t=-5.0;
  pnl_cdf_t(&which,&p,&q,&t,&df,&status,&bound);
  printf("df=%f \n",df);
  if (status != 0) printf("status=%d \n",status); 
}

void cumulfunc_test()
{
  printf("\n");
  printf("TEST DES FONCTIONS DU FICHIER CUMULFUNC.C \n");
  printf("\n");
  pnl_cdf_bet_test();
  pnl_cdf_bin_test();
  pnl_cdf_chi_test();
  pnl_cdf_chn_test();
  pnl_cdf_f_test();
  pnl_cdf_fnc_test();
  pnl_cdf_gam_test();
  pnl_cdf_nbn_test();
  pnl_cdf_nor_test();
  pnl_cdf_poi_test();
  pnl_cdf_t_test();
}
