
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

#include "pnl_random.h"
#include "pnl_vector.h"
#include "pnl_mpi.h"

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

static int send_rng ()
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
      gen[i] = pnl_rand_add (rngtab[i]);
      pnl_rand_init (gen[i], 1, 10);
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
        printf ("%f ",pnl_rand_uni (gen[i]));
      printf("\n");
    }

  /*
   * Destroy
   */
  for ( i=0 ; i<count ; i++ )
    {
      pnl_rand_remove (gen[i]);
    }
  free (gen);
  free (rngtab);
  return MPI_SUCCESS;
}


static int recv_rng ()
{
  PnlRng **rngtab;
  int i, j, *gen;
  MPI_Status status;

  /*
   * Receive and add
   */
  if ((rngtab = malloc (NB_GEN * sizeof(PnlRng *))) == NULL) return MPI_ERR_BUFFER;
  if ((gen = malloc (NB_GEN * sizeof(int))) == NULL) return MPI_ERR_BUFFER;
  for ( i=0 ; i<NB_GEN ; i++ )
    {
      rngtab[i] = pnl_rng_new ();
      printf ("receiving gen %d\n", i);
      PNL_MPI_MESSAGE(pnl_object_mpi_recv (PNL_OBJECT(rngtab[i]), 0, SENDTAG, MPI_COMM_WORLD, &status),
                      "error in receving rng");
      gen[i] = pnl_rand_add (rngtab[i]);
    }


  /*
   * Use
   */
  for ( i=0 ; i<NB_GEN ; i++ )
    {
      for ( j=0 ; j<NB_INT ; j++ )
        printf ("%f ",pnl_rand_uni (gen[i]));
      printf("\n");
    }
  
  /*
   * Destroy
   */
  for ( i=0 ; i<NB_GEN ; i++ )
    {
      pnl_rand_remove (gen[i]);
    }
  free (gen);
  free (rngtab);
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
      /* initializes some integers for the tests */
      send_rng (); 
    }
  else
    {
      recv_rng (); 
    }

  MPI_Finalize ();
  exit (0);
}
