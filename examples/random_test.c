
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

static double square (double x) { return x*x;}

static void test_pnl_vect_rand(int type_generator)
{
  int samples = 10000;
  double sum;
  int i;
  PnlVect *G;
  G = pnl_vect_create(0);
  pnl_rand_init(type_generator, 1, samples);

  pnl_vect_rand_uni(G, samples, 0, 1, type_generator);
  pnl_vect_map_inplace(G, square);
  sum = pnl_vect_sum(G)/samples;
  printf("E (U^2) = %f (should be 1/3)\n", sum);

  /* Calling pnl_rand_init again ensures that the next samples will be drawn
     from the beginning of the sequence. Very important for QMC */
  pnl_rand_init(type_generator, 1, samples);    
  sum=0.0;
  for (i=0; i<samples; i++)
    {
      sum += square(pnl_rand_normal(type_generator));
    }
  printf("E (G^2) = %f (should be 1)\n", sum/samples);

  /* equivalently as above */
  pnl_rand_init(type_generator, 1, samples);    
  pnl_vect_rand_normal(G, samples, type_generator);
  pnl_vect_map_inplace(G, square);
  sum = pnl_vect_sum(G)/samples;
  printf("E (G^2) = %f (should be 1)\n", sum);

  pnl_vect_free(&G);
}

static void test_pnl_mat_rand(int type_generator)
{
  PnlMat *M;
  PnlVect *inf, *sup, *Vsum;
  int samples = 100000;
  int dim = 10;
  int i;
  double sum;
  inf = pnl_vect_create_from_double(dim, 0.0);
  sup = pnl_vect_create_from_double(dim, 1.0);
  M = pnl_mat_create(0,0);
  Vsum = pnl_vect_create (0);
  pnl_rand_init(type_generator, dim, samples);
  pnl_mat_rand_uni(M, samples, dim, inf, sup, type_generator);
  pnl_mat_sum_vect(Vsum, M, 'r');
  pnl_vect_div_double(Vsum, samples);
  pnl_vect_print(Vsum);
  printf("(should be 0.5)\n");
  pnl_vect_free(&Vsum);
  pnl_vect_free(&inf);
  pnl_vect_free(&sup);
    
    
  pnl_rand_init(type_generator, dim, samples);
  pnl_mat_rand_normal(M, samples, dim, type_generator);
  sum = 0.0;
  for (i=0; i<samples; i++)
    {
      sum += MGET(M,i,0)*MGET(M,i,8);
    }
  printf("Cov = %f (should be 0)\n", sum/samples);

  pnl_mat_free(&M);
}

static void test_pnl_rand_gauss(int type_generator)
{
  int dimension=10, samples=100000;
  int i;
  double g1,g2,sum=0.0;
  pnl_rand_init(type_generator, dimension, samples);
  for (i=0; i<samples; i++)
    {
      g1=pnl_rand_gauss(dimension, CREATE, 0, type_generator);
      g2=pnl_rand_gauss(dimension, RETRIEVE, 8, type_generator);
      sum += g1*g2;
    }
  printf("Cov = %f (should be 0)\n", sum/samples);
}

static void test_rng ()
{
  int i, type_gen;
  PnlRng *rng;
  i = 0;
  while (PnlRngArray[i].rng != NULL)
    {
      rng = PnlRngArray[i].rng;
      type_gen = rng->type;
      printf ("--> Generator %s\n", PnlRngArray[i].base.label);
      test_pnl_vect_rand(type_gen);
      test_pnl_mat_rand(type_gen);
      test_pnl_rand_gauss(type_gen);
      i++;
    }
}

static void rng_call ()
{
  int j, N=10000;
  double sum, var;
  PnlRng *rng;

  printf ("\n--> Test of the rng interface\n");
  rng = pnl_rng_create (PNL_RNG_MERSENNE);
  pnl_rng_sseed (rng, 4172);
  var = sum = 0.;
  
  for ( j=0 ; j<N ; j++ )
    {
      double tmp;
      tmp = pnl_rng_normal(rng);
      sum += tmp; var += tmp * tmp;
    }
  
    printf ("mean = %f (sould be 0.) \tvar = %f (souble be 1)\n", sum/N, var/N);
    pnl_rng_free (&rng);
}


#define NGEN 3
static void std_call_dcmt ()
{
  int i, j, N=10000, count;
  dcmt_state **mts;
  double sum[NGEN], var[NGEN];

  
  printf ("\n--> Test of direct call to DCMT\n");
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
      printf ("mean = %f (sould be 0.5) \tvar = %f (souble be 1/3)\n", sum[j]/N, var[j]/N);
    }
  pnl_dcmt_free_array(mts, count);
}  


static void rand_call_dcmt ()
{
  int i, j, N=10000;
  PnlRng **rng;
  double sum[NGEN], var[NGEN];

  printf ("\n--> Test of call to DCMT through pnl_rng_xxx\n");
  rng = malloc (NGEN * sizeof(PnlRng *));
  
  for ( j=0 ; j<NGEN ; j++ )
    {
      sum[j] = var[j] = 0.;
      /* if ( (rng[j] = pnl_rng_create (PNL_RNG_DCMT)) == NULL)
       *   {
       *     perror ("Cannot create a DCMT.\n"); abort();
       *   }
       * pnl_rng_sseed (rng[j], 1234); */
    }
  rng[0] = pnl_rng_create (PNL_RNG_DCMT);
  pnl_rng_sseed (rng[0], 1234);
  rng[1] = pnl_rng_create (PNL_RNG_DCMT);
  pnl_rng_sseed (rng[1], 1234);
  rng[2] = pnl_rng_create (PNL_RNG_DCMT);
  pnl_rng_sseed (rng[2], 1234);

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
      printf ("mean = %f (sould be 0.5) \tvar = %f (souble be 1/3)\n", sum[j]/N, var[j]/N);
      pnl_rng_free(&(rng[j]));
    }
  free (rng);
}  


void random_test(void)
{
  /* test_rng();
   * rng_call ();
   * std_call_dcmt (); */
  rand_call_dcmt ();
}
