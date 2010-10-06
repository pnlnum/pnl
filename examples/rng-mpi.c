
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
#include <mpi.h>

#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_mpi.h"

#define SENDTAG 1
#define PNL_MPI_MESSAGE(info, msg)                             \
  if ( info != MPI_SUCCESS )                                   \
    {                                                          \
      if (pnl_message_is_on ()) printf (msg);                  \
      return info;                                             \
    }

#define NB_INT 10
#define NB_GEN 2
/*
 * Send / Recv examples
 */

static int send_dcmt ()
{
  PnlRng **rngtab;
  int i, j, count, *gen;

  /*
   * Create
   */
  rngtab = pnl_rng_dcmt_create_array (NB_GEN, 4172, &count);
  if ((gen = malloc (count * sizeof(int))) == NULL) return MPI_ERR_BUFFER;
  for ( i=0 ; i<count ; i++ )
    {
      pnl_rng_sseed (rngtab[i], 1273);
    }

  /*
   * Send
   */
  for ( i=0 ; i<count ; i++ )
    {
      printf ("sending gen %d\n", i);
      PNL_MPI_MESSAGE(pnl_object_mpi_send (PNL_OBJECT(rngtab[i]), 1, SENDTAG, MPI_COMM_WORLD), "error in sending rng");
    }

  /*
   * Use
   */
  for ( i=0 ; i<count ; i++ )
    {
      for ( j=0 ; j<NB_INT ; j++ )
        printf ("%f ",pnl_rng_uni (rngtab[i]));
      printf("\n");
    }

  /*
   * Destroy
   */
  for ( i=0 ; i<count ; i++ )
    {
      pnl_rng_free (&(rngtab[i]));
    }
  free (rngtab);
  return MPI_SUCCESS;
}


static int recv_dcmt ()
{
  PnlRng **rngtab;
  int i, j;
  MPI_Status status;

  /*
   * Receive and add
   */
  if ((rngtab = malloc (NB_GEN * sizeof(PnlRng *))) == NULL) return MPI_ERR_BUFFER;
  for ( i=0 ; i<NB_GEN ; i++ )
    {
      rngtab[i] = pnl_rng_new ();
      printf ("receiving gen %d\n", i);
      PNL_MPI_MESSAGE(pnl_object_mpi_recv (PNL_OBJECT(rngtab[i]), 0, SENDTAG, MPI_COMM_WORLD, &status),
                      "error in receving rng");
    }


  /*
   * Use
   */
  for ( i=0 ; i<NB_GEN ; i++ )
    {
      for ( j=0 ; j<NB_INT ; j++ )
        printf ("%f ",pnl_rng_uni (rngtab[i]));
      printf("\n");
    }
  
  /*
   * Destroy
   */
  for ( i=0 ; i<NB_GEN ; i++ )
    {
      pnl_rng_free (&(rngtab[i]));
    }
  free (rngtab);
  return MPI_SUCCESS;

}


static int send_rng ()
{
  PnlRng *rng;
  int i;

  rng = pnl_rng_create (PNL_RNG_MERSENNE);
  pnl_rng_sseed (rng, 1273);
  PNL_MPI_MESSAGE(pnl_object_mpi_send (PNL_OBJECT(rng), 1, SENDTAG, MPI_COMM_WORLD), "error in sending rng");

  for ( i=0 ; i<NB_INT ; i++ )
    printf ("%f ",pnl_rng_uni (rng));
  printf("\n");

  pnl_rng_free (&rng);
  return MPI_SUCCESS;
}

static int recv_rng ()
{
  PnlRng *rng;
  int j;
  MPI_Status status;

  rng = pnl_rng_new ();
  PNL_MPI_MESSAGE(pnl_object_mpi_recv (PNL_OBJECT(rng), 0, SENDTAG, MPI_COMM_WORLD, &status),
                  "error in receving rng");

  for ( j=0 ; j<NB_INT ; j++ )
    printf ("%f ",pnl_rng_uni (rng));
  
  printf("\n");
  pnl_rng_free (&rng);
  return MPI_SUCCESS;

}


int main(int argc, char *argv[])
{
  int rank, nproc;
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  pnl_message_on();
  if ( nproc != 2 )
    {
      if ( rank == 0 )
        {
          printf("Run the test with -np 2.\n");
        }
      MPI_Finalize ();
      exit (0);
    }

  if ( rank == 0 )
    {
      send_dcmt ();  MPI_Barrier(MPI_COMM_WORLD);
      send_rng ();  MPI_Barrier(MPI_COMM_WORLD);
    }
  else
    {
      recv_dcmt ();   MPI_Barrier(MPI_COMM_WORLD);
      recv_rng ();   MPI_Barrier(MPI_COMM_WORLD);
    }

  MPI_Finalize ();
  exit (0);
}
