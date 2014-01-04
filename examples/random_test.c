
/************************************************************************/
/* Copyright Jérôme Lelong <jerome.lelong@gmail.com>                    */
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

#include "pnl/pnl_random.h"
#include  "tests_utils.h"

#define N_SAMPLES 100000

static double square (double x) { return x*x;}

static void reset_rng (PnlRng *rng, int s)
{
  if ( rng->rand_or_quasi == PNL_MC )
    {
      pnl_rng_sseed (rng, 123456);
    }
  else
    {
      if ( pnl_rng_sdim (rng, s) == FAIL )
        {
          printf ("Initialization error of generator %d\n", rng->type);
        }
    }
}

static void test_pnl_vect_rng(PnlRng *rng, const char *name)
{
  int samples = N_SAMPLES;
  double sum;
  int i;
  PnlVect *G;
  char str[256];
  G = pnl_vect_create(0);

  reset_rng (rng, 1);
  pnl_vect_rng_uni(G, samples, 0, 1, rng);
  pnl_vect_map_inplace(G, square);
  sum = pnl_vect_sum(G)/samples;
  sprintf (str, "E(U^2) with %s", name);
  pnl_test_eq_abs (sum, 1. / 3., 1E-2, str, "");

  /* Calling pnl_rng_init again ensures that the next samples will be drawn
     from the beginning of the sequence. Very important for QMC */
  reset_rng (rng, 1);
  sum=0.0;
  for (i=0; i<samples; i++)
    {
      sum += square(pnl_rng_normal(rng));
    }
  sprintf (str, "E(G^2) with %s", name);
  pnl_test_eq_abs (sum/samples, 1., 1E-2, str, "" );

  /* equivalently as above */
  reset_rng(rng, 1);
  pnl_vect_rng_normal(G, samples, rng);
  pnl_vect_map_inplace(G, square);
  sum = pnl_vect_sum(G)/samples;
  sprintf (str, "E(G^2) with %s (vectorized version)", name);
  pnl_test_eq_abs (sum, 1., 1E-2, str, "");

  pnl_vect_free(&G);
}

static void test_pnl_mat_rng(PnlRng *rng, const char *name)
{
  PnlMat *M;
  PnlVect *inf, *sup, *Vsum;
  int samples = N_SAMPLES;
  int dim = 10;
  int i;
  double sum;
  char str[256];
  inf = pnl_vect_create_from_scalar(dim, 0.0);
  sup = pnl_vect_create_from_scalar(dim, 1.0);
  M = pnl_mat_create(0,0);
  Vsum = pnl_vect_create (0);
  reset_rng (rng, dim);
  pnl_mat_rng_uni(M, samples, dim, inf, sup, rng);
  pnl_mat_sum_vect(Vsum, M, 'r');
  pnl_vect_div_scalar(Vsum, samples);
  /* we should ckeck that Vsum is full of 0.5 */
  pnl_vect_free(&Vsum);
  pnl_vect_free(&inf);
  pnl_vect_free(&sup);
    
    
  reset_rng (rng, dim);
  pnl_mat_rng_normal(M, samples, dim, rng);
  sum = 0.0;
  for (i=0; i<samples; i++)
    {
      sum += MGET(M,i,0)*MGET(M,i,8);
    }
  sprintf (str, "Covariance %s", name);
  pnl_test_eq_abs (sum/samples, 0., 1E-2, str, "");

  pnl_mat_free(&M);
}

static void test_pnl_rng_gauss(PnlRng *rng, const char *name)
{
  int dimension=10, samples=N_SAMPLES;
  int i;
  double g1,g2,sum=0.0;
  char str[256];
  reset_rng (rng, dimension);
  for (i=0; i<samples; i++)
    {
      g1=pnl_rng_gauss(dimension, CREATE, 0, rng);
      g2=pnl_rng_gauss(dimension, RETRIEVE, 8, rng);
      sum += g1*g2;
    }
  sprintf (str, "Covariance %s Gaussian case", name);
  pnl_test_eq_abs (sum/samples, 0., 1E-2, str, "");
}

static void test_rng ()
{
  int i;
  PnlRng *rng;
  i = 0;
  while (PnlRngArray[i] != NULL)
    {
      rng = PnlRngArray[i];
      test_pnl_vect_rng(rng, pnl_rng_get_name(PnlRngArray[i]->type));
      test_pnl_mat_rng(rng, pnl_rng_get_name(PnlRngArray[i]->type));
      test_pnl_rng_gauss(rng, pnl_rng_get_name(PnlRngArray[i]->type));
      i++;
    }
}

#define NGEN 3
static void std_call_dcmt ()
{
  int i, j, N=10000, count;
  dcmt_state **mts;
  double sum[NGEN], var[NGEN];

  
  mts = pnl_dcmt_create_array(NGEN,4172,&count);
  if (count != NGEN)
    {
      printf ("Error in creating dcmt\n"); abort ();
    }
  for ( j=0 ; j<count ; j++ )
    {
      sum[j] = var[j] = 0.;
      pnl_dcmt_sseed (mts[j], 1234);
    }

  for (i=0; i<N; i++) 
    {
      for ( j=0 ; j<count ; j++ )
        {
          double tmp;
          tmp = pnl_dcmt_genrand_double(mts[j]);
          sum[j] += tmp; var[j] += tmp * tmp;
        }
    }
  for ( j=0 ; j<count ; j++ )
    {
      pnl_test_eq_abs (sum[j]/N, 0.5, 1E-2, "DCMT E(U)", "");
      pnl_test_eq_abs (var[j]/N, 1. / 3., 1E-2,  "DCMT E(U^2)", "");
    }
  pnl_dcmt_free_array(mts, count);
}  

static void rand_call_dcmt ()
{
  int i, j, N=10000;
  PnlRng **rng;
  double sum[NGEN], var[NGEN];

  rng = malloc (NGEN * sizeof(PnlRng *));
  
  for ( j=0 ; j<NGEN ; j++ )
    {
      sum[j] = var[j] = 0.;
      if ( (rng[j] = pnl_rng_dcmt_create_id (j, 0)) == NULL)
        {
          perror ("Cannot create a DCMT.\n"); abort();
        }
      pnl_rng_sseed (rng[j], 1234);
    }

  for (i=0; i<N; i++) 
    {
      for ( j=0 ; j<NGEN ; j++ )
        {
          double tmp;
          tmp = pnl_rng_uni(rng[j]);
          sum[j] += tmp; var[j] += tmp * tmp;
        }
    }
  for ( j=0 ; j<NGEN ; j++ )
    {
      pnl_test_eq_abs (sum[j]/N, 0.5, 1E-2, "rng_dcmt E(U)", "");
      pnl_test_eq_abs (var[j]/N, 1. / 3., 1E-2,  "rng_dcmt E(U^2)", "");
      pnl_rng_free(&(rng[j]));
    }
  free (rng);
}  

int main (int argc, char **argv)
{
  pnl_test_init (argc, argv);
  test_rng();
  std_call_dcmt ();
  rand_call_dcmt ();
  exit (pnl_test_finalize("Random generators"));
}
