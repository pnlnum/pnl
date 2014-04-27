
/************************************************************************/
/* Copyright Céline Labart <labart@cmap.polytechnique.fr                */
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

/*
 * The expected values for the different tests have been computed using
 * Scilab with format (20)
 *
 * Nsp does not always provide exactly the same values, in particular for
 * computations with *which==5
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pnl/pnl_cdf.h"
#include "tests_utils.h"

static double abserr = 1E-8;

static void pnl_cdf_bet_test()
{
  int which;
  double p, pe;
  double q, qe;
  double x;
  double y;
  double a;
  double b;
  int status;
  double bound;
  which=1;
  x=0.6;
  y=0.4;
  a=2.0;
  b=3.0;
  pe = 0.820799999999999974;
  qe = 0.179200000000000081;
  pnl_cdf_bet(&which,&p,&q,&x,&y,&a,&b,&status,&bound);
  pnl_test_eq_abs (p, pe, abserr, "cdf_bet p", "");
  pnl_test_eq_abs (q, qe, abserr, "cdf_bet q", "");
}

static void pnl_cdf_bin_test()
{
  int which;
  double p, pe;
  double q, qe;
  double s;
  double xn;
  double pr;
  double ompr;
  int status;
  double bound;
  which=1;
  s=2;
  xn=5;
  pr=0.7;
  ompr=0.3;
  pe = 0.163079999999999919;
  qe = 0.836920000000000108;
  pnl_cdf_bin(&which,&p,&q,&s,&xn,&pr,&ompr,&status,&bound);
  pnl_test_eq_abs (p, pe, abserr, "cdf_bin p", "");
  pnl_test_eq_abs (q, qe, abserr, "cdf_bin q", "");
}

static void pnl_cdf_chi_test()
{
  int which;
  double p;
  double q;
  double x, xe;
  double df, dfe;
  int status;
  double bound;
  which=3;
  p=0.7;
  q=0.3;
  x=2;
  dfe = 1.686206763411698173;
  pnl_cdf_chi(&which,&p,&q,&x,&df,&status,&bound);
  pnl_test_eq_abs (df, dfe, abserr, "cdf_chi df", "");
  which = 2;
  df = 6;
  xe = 7.231135331731981530;
  pnl_cdf_chi(&which,&p,&q,&x,&df,&status,&bound);
  pnl_test_eq_abs (x, xe, abserr, "cdf_chi x", "");
}

static void pnl_cdf_chn_test()
{
  int which;
  double p;
  double q;
  double x;
  double df;
  double pnonc, pnonce;
  int status;
  double bound;
  which=4;
  p=0.1;
  q=0.9;
  x=2;
  df=2;
  pnonce = 5.812579835152692276;
  pnl_cdf_chn(&which,&p,&q,&x,&df,&pnonc,&status,&bound);
  pnl_test_eq_abs (pnonc, pnonce, abserr, "cdf_chn pnonc", "");
}

static void pnl_cdf_f_test()
{
  int which;
  double p;
  double q;
  double f;
  double dfn;
  double dfd, dfde;
  int status;
  double bound;
  which=4;
  p=0.1;
  q=0.9;
  f=2;
  dfn=3;
  dfde = 0.050304795796988773;
  pnl_cdf_f(&which,&p,&q,&f,&dfn,&dfd,&status,&bound);
  pnl_test_eq_abs (dfd, dfde, abserr, "cdf_f dfd", "");
}

static void pnl_cdf_fnc_test()
{
  int which;
  double p;
  double q;
  double f;
  double dfn;
  double dfd;
  double pnonc, pnonce;
  int status;
  double bound;
  which=5;
  p=0.1;
  q=0.9;
  f=2;
  dfn=3;
  dfd=5;
  pnonce = 13.17874475510091337;
  pnl_cdf_fnc(&which,&p,&q,&f,&dfn,&dfd,&pnonc,&status,&bound);
  pnl_test_eq_abs (pnonc, pnonce, abserr, "cdf_fnc pnonc", "");
}

static void pnl_cdf_gam_test()
{
  int which;
  double p;
  double q;
  double x;
  double shape;
  double scale, scalee;
  int status;
  double bound;
  which=4;
  p=0.1;
  q=0.9;
  x=2;
  shape=3;
  scalee = 0.551032664124660565;
  pnl_cdf_gam(&which,&p,&q,&x,&shape,&scale,&status,&bound);
  pnl_test_eq_abs (scale, scalee, abserr, "cdf_gam scale", "");
}

static void pnl_cdf_nbn_test()
{
  int which;
  double p, pe;
  double q, qe;
  double s;
  double xn;
  double pr;
  double ompr;
  int status;
  double bound;
  which=1;
  s=2.0;
  xn=3.0;
  pr=0.7;
  ompr=0.3;
  pe = 0.836920000000000108;
  qe = 0.163079999999999919;
  pnl_cdf_nbn(&which,&p,&q,&s,&xn,&pr,&ompr,&status,&bound);
  pnl_test_eq_abs (p, pe, abserr, "cdf_nbn p", "");
  pnl_test_eq_abs (q, qe, abserr, "cdf_nbn q", "");
}

static void pnl_cdf_nor_test()
{
  int which;
  double p, pe;
  double q, qe;
  double x;
  double mean;
  double sd;
  int status;
  double bound;
  which=1;
  x=2.0;
  mean=0.0;
  sd=1.0;
  pe = 0.977249868051820791;
  qe = 0.022750131948179212;
  pnl_cdf_nor(&which,&p,&q,&x,&mean,&sd,&status,&bound);
  pnl_test_eq_abs (p, pe, abserr, "cdf_nor p", "");
  pnl_test_eq_abs (q, qe, abserr, "cdf_nor q", "");
}

static void pnl_cdf_poi_test()
{
  int which;
  double p;
  double q;
  double s, se;
  double xlam, xlame;
  int status;
  double bound;
  which=3;
  p=0.4;
  q=0.6;
  se = s = 5.0;
  xlame = 6.291918983308751656;
  pnl_cdf_poi(&which,&p,&q,&s,&xlam,&status,&bound);
  pnl_test_eq_abs (xlam, xlame, abserr, "cdf_poi xlam", "");
  which = 2;
  pnl_cdf_poi(&which, &p, &q, &s, &xlam, &status, &bound);
  pnl_test_eq_abs (s, se, abserr, "cdf_poi inv", "");
}

static void pnl_cdf_t_test()
{
  int which;
  double p;
  double q;
  double t;
  double df, dfe;
  int status;
  double bound;
  which=3;
  p=0.4;
  q=0.6;
  t=-5.0;
  dfe = 0.060622362518168812;
  pnl_cdf_t(&which,&p,&q,&t,&df,&status,&bound);
  pnl_test_eq_abs (df, dfe, abserr, "cdf_t df", "");
}

int main (int argc, char **argv)
{
  pnl_test_init (argc, argv);
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
  exit (pnl_test_finalize ("CDF"));
}
