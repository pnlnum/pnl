
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

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include "pnl/pnl_random.h"
#include "pnl/pnl_mpi.h"

#define NB_INT 5
#define NB_GEN 3
static char *binfile = "_tmp_poo.bin";

static void save_rng_array ()
{
  int i, j, n, count;
  PnlRng **rngtab;
  n = NB_GEN;
  rngtab = pnl_rng_dcmt_create_array (n, 4172, &count);
  if ( n != count )
    {
      perror ("Wrong number of generator created\n");
    }

  for ( i=0 ; i<n ; i++ )
    pnl_rng_sseed (rngtab[i], time(NULL));

  pnl_rng_save_to_file (rngtab, n, binfile);
  for ( i=0 ; i<n ; i++ )
    {
      for ( j=0 ; j<NB_INT ; j++ )
        printf ("%f ",pnl_rng_uni (rngtab[i]));
      printf("\n");
      pnl_rng_free (&(rngtab[i]));
    }
  free(rngtab);
}

static void load_rng_array ()
{
  int i, j;
  PnlRng **rngtab;
  rngtab = pnl_rng_create_from_file (binfile, NB_GEN);
  for ( i=0 ; i<NB_GEN ; i++ )
    {
      for ( j=0 ; j<NB_INT ; j++ )
        printf ("%f ",pnl_rng_uni (rngtab[i]));
      printf("\n");
      pnl_rng_free (&(rngtab[i]));
    }
  free(rngtab);
  unlink(binfile);
}

static int save_rng (PnlType t)
{
  PnlRng *rng1, *rng2;
  FILE *stream;
  int i;

  stream = fopen (binfile, "wb");

  rng1 = pnl_rng_create (t);
  pnl_rng_sseed (rng1, 1273);
  pnl_object_save (PNL_OBJECT(rng1), stream); 
  for ( i=0 ; i<NB_INT ; i++ )
    printf ("%f ",pnl_rng_uni (rng1));
  printf("\n");
  pnl_rng_free (&rng1);


  rng2 = pnl_rng_create (t);
  pnl_rng_sseed (rng2, 2273);
  pnl_object_save (PNL_OBJECT(rng2), stream); 
  for ( i=0 ; i<NB_INT ; i++ )
    printf ("%f ",pnl_rng_uni (rng2));
  printf("\n");
  pnl_rng_free (&rng2);

  fclose (stream);
  return MPI_SUCCESS;
}

static int load_rng ()
{

  PnlRng *rng1, *rng2;
  PnlList *L;
  FILE *stream;
  int i;

  stream = fopen (binfile, "rb");
  rng1 = PNL_RNG_OBJECT(pnl_object_load (stream)); 
  rng2 = PNL_RNG_OBJECT(pnl_object_load (stream)); 
  fclose (stream);

  for ( i=0 ; i<NB_INT ; i++ )
    printf ("%f ",pnl_rng_uni (rng1));
  printf("\n");
  pnl_rng_free (&rng1);


  for ( i=0 ; i<NB_INT ; i++ )
    printf ("%f ",pnl_rng_uni (rng2));
  printf("\n");
  pnl_rng_free (&rng2);

  stream = fopen (binfile, "rb");
  L = pnl_object_load_into_list (stream); 
  fclose (stream);

  rng1 = PNL_RNG_OBJECT(pnl_list_get (L, 0));
  rng2 = PNL_RNG_OBJECT(pnl_list_get (L, 1));

  for ( i=0 ; i<NB_INT ; i++ )
    printf ("%f ",pnl_rng_uni (rng1));
  printf("\n");

  for ( i=0 ; i<NB_INT ; i++ )
    printf ("%f ",pnl_rng_uni (rng2));
  printf("\n");

  pnl_list_free (&L);
  unlink(binfile);
  return MPI_SUCCESS;
}


int main (int argc, char **argv)
{
  PnlRngType t;
  MPI_Init (&argc, &argv);
  t = PNL_RNG_MERSENNE;
  printf ("--> save/load interface\n");
  save_rng(t);
  load_rng();
  printf ("--> save/load for array of rngs\n");
  save_rng_array();
  load_rng_array();
  MPI_Finalize();
  return OK;
}
