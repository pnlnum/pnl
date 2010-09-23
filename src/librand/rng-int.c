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

#include <limits.h>
#include <math.h>
#include <stdio.h>


#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_random.h"

static char pnl_rng_label[] = "PnlRng";

void pnl_rng_free (PnlRng **rng)
{
  if ( *rng == NULL ) return;
  switch ( (*rng)->type )
    {
    case PNL_RNG_MERSENNE:
    case PNL_RNG_DCMT:
      if ( (*rng)->state != NULL )
        {
          free ((*rng)->state); (*rng)->state = NULL;
        }
      break;
    }
  free (*rng); *rng = NULL;
}

/**
 * Creates an empty rng
 */
PnlRng* pnl_rng_new ()
{
  PnlRng *rng;
  if ((rng = malloc (sizeof(PnlRng))) == NULL) return NULL;
  
  rng->object.type = PNL_TYPE_RNG;
  rng->object.label = pnl_rng_label;
  rng->object.parent_type = PNL_TYPE_RNG;
  rng->object.destroy = (destroy_func *) pnl_rng_free;

  rng->type = NULLINT;
  rng->id = NULLINT;
  rng->Compute = NULL;
  rng->rand_or_quasi = MC;
  rng->dimension = 0;
  rng->max_dim = INT_MAX;
  rng->counter = 0;
  rng->has_gauss = 0;
  rng->gauss = 0;
  rng->size_state = 0;
  rng->state = NULL;
  return rng;
}


static void MERSENNE(PnlRng *rng,double *sample)
{
  mt_state *state = (mt_state *) (rng->state);

  *sample = pnl_mt_genrand_double (state);
  rng->counter++;
}


/**
 * Creates a rng of the given type.
 * Note that the fields size_state and state are set to zero, which implies
 * that the created generator is unusable. It is only usefull to receive an
 * already workin generator.
 *
 * @param type the type of generator to create
 * @return a PnlRng or NULL if an error occurred
 */
PnlRng* pnl_rng_create (int type)
{
  PnlRng *rng;
  if ((rng = pnl_rng_new()) == NULL) return NULL;

  rng->type = type;
  rng->dimension = 0;
  rng->counter = 0;
  rng->has_gauss = 0;
  rng->gauss = 0;
  rng->size_state = 0;
  rng->state = NULL;
  switch (type)
    {
    case PNL_RNG_MERSENNE:
      rng->Compute = MERSENNE;
      rng->rand_or_quasi = MC;
      rng->size_state = sizeof(mt_state);
      rng->state = malloc(rng->size_state);
      break;
    default:
      free (rng); rng = NULL;
      printf("Unknown generator type\n");
    }
  return rng;
}


PnlRng** pnl_rng_dcmt_create_array (int n, ulong seed, int *count)
{
  PnlRng **rngtab;
  dcmt_state **mts;
  int i;

  mts = pnl_dcmt_create_array (n, seed, count);
  if (n != *count)
    {
      perror ("Not all generators could be created\n");
      abort ();
    }
  if ( (rngtab = malloc (n * sizeof(PnlRng *))) == NULL)
    {
      pnl_dcmt_free_array (mts, n); return NULL;
    }
  for ( i=0 ; i<n ; i++ )
    {
      rngtab[i] = pnl_rng_new ();
      pnl_rng_init (rngtab[i], PNL_RNG_DCMT);
      rngtab[i]->size_state = sizeof (dcmt_state);
      rngtab[i]->state = mts[i];
    }
  free (mts);
  return rngtab;
}

void pnl_rng_sseed (PnlRng *rng, ulong seed)
{
  switch (rng->type)
    {
    case PNL_RNG_MERSENNE :
    case PNL_RNG_MERSENNE_RANDOM_SEED :
      pnl_mt_sseed((mt_state *)(rng->state), seed);
      break;
    case PNL_RNG_DCMT :
      pnl_dcmt_sseed ((dcmt_state *)(rng->state), seed);
    }
  rng->counter=1;
  rng->has_gauss=0;
  rng->gauss=0.;
}
