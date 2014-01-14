
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
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_band_matrix.h"
#include "pnl/pnl_tridiag_matrix.h"
#include "pnl/pnl_basis.h"
#include "pnl/pnl_list.h"
#include "pnl/pnl_mpi.h"
#include "tests_utils.h"

#define SENDTAG 1
#define PNL_MPI_MESSAGE(info, msg)                             \
  if ( info != MPI_SUCCESS )                                   \
    {                                                          \
      if (pnl_message_is_on ()) printf (msg);                  \
      return info;                                             \
    }

#define NB_INT 100
static int IntArray[NB_INT];

static void init_int_array ()
{
  int i, gen;
  gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init (gen, NB_INT, 1);
  for ( i=0 ; i<NB_INT ; i++ ) { IntArray[i] = (int) (pnl_rand_uni (gen) * INT_MAX); }
}

/*
 * Send / Recv examples
 */

static int send_vector ()
{
  PnlVect *v;
  int info, n =5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  v = pnl_vect_new ();
  pnl_rand_init (gen, n, 1);
  pnl_vect_rand_uni (v, n, 0, 1, gen);
  printf ("Original vector \n"); pnl_vect_print (v); printf ("\n");
  info = pnl_object_mpi_send (PNL_OBJECT(v), 1, SENDTAG, MPI_COMM_WORLD);
  pnl_vect_free (&v);
  return info;
}

static int recv_vector ()
{
  MPI_Status status;
  PnlVect *v;
  int info;
  v = pnl_vect_new ();
  info = pnl_object_mpi_recv (PNL_OBJECT(v), 0, SENDTAG, MPI_COMM_WORLD, &status);
  printf ("Received vector \n"); pnl_vect_print (v); printf ("\n");
  pnl_vect_free (&v);
  return info;
}

static int send_int_vector ()
{
  PnlVectInt *v;
  int info, n =5;
  v = pnl_vect_int_create_from_ptr (n, IntArray);
  printf ("Original vector \n"); pnl_vect_int_print (v); printf ("\n");
  info = pnl_object_mpi_send (PNL_OBJECT(v), 1, SENDTAG, MPI_COMM_WORLD);
  pnl_vect_int_free (&v);
  return info;
}

static int recv_int_vector ()
{
  MPI_Status status;
  PnlVectInt *v;
  int info;
  v = pnl_vect_int_new ();
  info = pnl_object_mpi_recv (PNL_OBJECT(v), 0, SENDTAG, MPI_COMM_WORLD, &status);
  printf ("Received vector \n"); pnl_vect_int_print (v); printf ("\n");
  pnl_vect_int_free (&v);
  return info;
}

static int send_complex_vector ()
{
  PnlVect *v;
  PnlVectComplex *M;
  int info, m = 5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  v = pnl_vect_new ();
  pnl_rand_init (gen, m, 1);
  pnl_vect_rand_normal (v, 2*m, gen);
  M = pnl_vect_complex_create_from_ptr (m, (dcomplex *)v->array);
  printf ("Original vector \n"); pnl_vect_complex_print (M); printf ("\n");
  info = pnl_object_mpi_send (PNL_OBJECT(M), 1, SENDTAG, MPI_COMM_WORLD);
  pnl_vect_complex_free (&M);
  pnl_vect_free (&v);
  return info;
}

static int recv_complex_vector ()
{
  MPI_Status status;
  PnlVectComplex *M;
  int info;
  M = pnl_vect_complex_new ();
  info = pnl_object_mpi_recv (PNL_OBJECT(M), 0, SENDTAG, MPI_COMM_WORLD, &status);
  printf ("Received vector \n"); pnl_vect_complex_print (M); printf ("\n");
  pnl_vect_complex_free (&M);
  return info;
}

static int send_matrix ()
{
  PnlMat *M;
  int info, m = 4, n =5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  M = pnl_mat_new ();
  pnl_rand_init (gen, m, n);
  pnl_mat_rand_normal (M, m, n, gen);
  printf ("Original matrix \n"); pnl_mat_print (M); printf ("\n");
  info = pnl_object_mpi_send (PNL_OBJECT(M), 1, SENDTAG, MPI_COMM_WORLD);
  pnl_mat_free (&M);
  return info;
}

static int recv_matrix ()
{
  MPI_Status status;
  PnlMat *M;
  int info;
  M = pnl_mat_new ();
  info = pnl_object_mpi_recv (PNL_OBJECT(M), 0, SENDTAG, MPI_COMM_WORLD, &status);
  printf ("Received matrix \n"); pnl_mat_print (M); printf ("\n");
  pnl_mat_free (&M);
  return info;
}

static double set_zero_neg (double x) { return MAX(0., x); }

static int send_sp_matrix ()
{
  PnlMat *M;
  PnlSpMat *Sp;
  int info, m = 7, n =9;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  M = pnl_mat_new ();
  pnl_rand_init (gen, m, n);
  pnl_mat_rand_normal (M, m, n, gen);
  pnl_mat_map_inplace (M, set_zero_neg);
  Sp = pnl_sp_mat_create_from_mat (M);
  printf ("Original sparse matrix \n"); pnl_sp_mat_print (Sp); printf ("\n");
  info = pnl_object_mpi_send (PNL_OBJECT(Sp), 1, SENDTAG, MPI_COMM_WORLD);
  pnl_mat_free (&M);
  pnl_sp_mat_free (&Sp);
  return info;
}

static int recv_sp_matrix ()
{
  MPI_Status status;
  PnlSpMat *Sp;
  int info;
  Sp = pnl_sp_mat_new ();
  info = pnl_object_mpi_recv (PNL_OBJECT(Sp), 0, SENDTAG, MPI_COMM_WORLD, &status);
  printf ("Received sparse matrix \n"); pnl_sp_mat_print (Sp); printf ("\n");
  pnl_sp_mat_free (&Sp);
  return info;
}


static int send_int_matrix ()
{
  PnlMatInt *M;
  int info, n =5, m=4;
  M = pnl_mat_int_create_from_ptr (m, n, IntArray);
  printf ("Original matrix \n"); pnl_mat_int_print (M); printf ("\n");
  info = pnl_object_mpi_send (PNL_OBJECT(M), 1, SENDTAG, MPI_COMM_WORLD);
  pnl_mat_int_free (&M);
  return info;
}

static int recv_int_matrix ()
{
  MPI_Status status;
  PnlMatInt *M;
  int info;
  M = pnl_mat_int_new ();
  info = pnl_object_mpi_recv (PNL_OBJECT(M), 0, SENDTAG, MPI_COMM_WORLD, &status);
  printf ("Received matrix \n"); pnl_mat_int_print (M); printf ("\n");
  pnl_mat_int_free (&M);
  return info;
}

static int send_complex_matrix ()
{
  PnlVect *v;
  PnlMatComplex *M;
  int info, m = 4, n =5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  v = pnl_vect_new ();
  pnl_rand_init (gen, m, n);
  pnl_vect_rand_normal (v, 2*m*n, gen);
  M = pnl_mat_complex_create_from_ptr (m, n, (dcomplex *)v->array);
  printf ("Original matrix \n"); pnl_mat_complex_print (M); printf ("\n");
  info = pnl_object_mpi_send (PNL_OBJECT(M), 1, SENDTAG, MPI_COMM_WORLD);
  pnl_mat_complex_free (&M);
  pnl_vect_free (&v);
  return info;
}

static int recv_complex_matrix ()
{
  MPI_Status status;
  PnlMatComplex *M;
  int info;
  M = pnl_mat_complex_new ();
  info = pnl_object_mpi_recv (PNL_OBJECT(M), 0, SENDTAG, MPI_COMM_WORLD, &status);
  printf ("Received matrix \n"); pnl_mat_complex_print (M); printf ("\n");
  pnl_mat_complex_free (&M);
  return info;
}

static int send_bandmatrix ()
{
  PnlMat *M;
  PnlBandMat *BM;
  int info, m = 4, n =5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  M = pnl_mat_new ();
  pnl_rand_init (gen, m, n);
  pnl_mat_rand_normal (M, m, n, gen);
  BM = pnl_band_mat_create_from_mat (M, 3,1);
  printf ("Original Band matrix \n"); pnl_band_mat_print_as_full (BM); printf ("\n");
  info = pnl_object_mpi_send (PNL_OBJECT(BM), 1, SENDTAG, MPI_COMM_WORLD);
  pnl_mat_free (&M);
  pnl_band_mat_free (&BM);
  return info;
}

static int recv_bandmatrix ()
{
  MPI_Status status;
  PnlBandMat *M;
  int info;
  M = pnl_band_mat_new ();
  info = pnl_object_mpi_recv (PNL_OBJECT(M), 0, SENDTAG, MPI_COMM_WORLD, &status);
  printf ("Received Band matrix \n"); pnl_band_mat_print_as_full (M); printf ("\n");
  pnl_band_mat_free (&M);
  return info;
}

static int send_tridiagmatrix ()
{
  PnlMat *M;
  PnlTridiagMat *TM;
  int info, m = 5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  M = pnl_mat_new ();
  pnl_rand_init (gen, m, m);
  pnl_mat_rand_normal (M, m, m, gen);
  TM = pnl_tridiag_mat_create_from_mat (M);
  printf ("Original Tridiag matrix \n"); pnl_tridiag_mat_print (TM); printf ("\n");
  info = pnl_object_mpi_send (PNL_OBJECT(TM), 1, SENDTAG, MPI_COMM_WORLD);
  pnl_mat_free (&M);
  pnl_tridiag_mat_free (&TM);
  return info;
}

static int send_basis ()
{
  PnlBasis *B;
  int info, index, degree, spaced;
  index = PNL_BASIS_CANONICAL;
  degree = 4;
  spaced =3;
  B = pnl_basis_create_from_degree (index, degree, spaced);
  printf ("Original Basis \n"); pnl_basis_print (B); printf ("\n");
  info = pnl_object_mpi_send (PNL_OBJECT(B), 1, SENDTAG, MPI_COMM_WORLD);
  PNL_MPI_MESSAGE (info, "Basis not sent\n");
  pnl_basis_free (&B);
  return info;
}

static int ssend_basis ()
{
  PnlBasis *B;
  int info, index, degree, spaced;
  index = PNL_BASIS_CANONICAL;
  degree = 4;
  spaced =3;
  B = pnl_basis_create_from_degree (index, degree, spaced);
  printf ("Original Basis \n"); pnl_basis_print (B); printf ("\n");
  info = pnl_object_mpi_ssend (PNL_OBJECT(B), 1, SENDTAG, MPI_COMM_WORLD);
  PNL_MPI_MESSAGE (info, "Basis not sent\n");
  pnl_basis_free (&B);
  return info;
}

static int recv_tridiagmatrix ()
{
  MPI_Status status;
  PnlTridiagMat *M;
  int info;
  M = pnl_tridiag_mat_new ();
  info = pnl_object_mpi_recv (PNL_OBJECT(M), 0, SENDTAG, MPI_COMM_WORLD, &status);
  printf ("Received Tridiag matrix \n"); pnl_tridiag_mat_print (M); printf ("\n");
  pnl_tridiag_mat_free (&M);
  return info;
}

static int send_int_hmatrix ()
{
  PnlHmatInt *M;
  int info, m = 4, n =2, p=3;
  int dims[3]; dims[0] = m; dims[1] = n; dims[2] = p;
  M = pnl_hmat_int_create_from_ptr (3, dims, IntArray);
  printf ("Original Hmatrix \n"); pnl_hmat_int_print (M); printf ("\n");
  info = pnl_object_mpi_send (PNL_OBJECT(M), 1, SENDTAG, MPI_COMM_WORLD);
  printf ("Hmat sent\n");
  pnl_hmat_int_free (&M);
  return info;
}

static int recv_int_hmatrix ()
{
  MPI_Status status;
  PnlHmatInt *M;
  int info;
  M = pnl_hmat_int_new ();
  info = pnl_object_mpi_recv (PNL_OBJECT(M), 0, SENDTAG, MPI_COMM_WORLD, &status);
  printf ("Received Hmatrix \n"); pnl_hmat_int_print (M); printf ("\n");
  pnl_hmat_int_free (&M);
  return info;
}

static int recv_basis ()
{
  MPI_Status status;
  PnlBasis *B;
  int info;
  B = pnl_basis_new ();
  info = pnl_object_mpi_recv (PNL_OBJECT(B), 0, SENDTAG, MPI_COMM_WORLD, &status);
  printf ("Received basis\n"); pnl_basis_print (B); printf ("\n");
  pnl_basis_free (&B);
  return info;
}

static int send_list ()
{
  int info;
  PnlList *L;
  PnlMat *M;
  PnlTridiagMat *T;

  L = pnl_list_new ();
  M = pnl_mat_create_from_scalar (2, 3, 3.5);
  T = pnl_tridiag_mat_create_from_scalar (4, 0.5);
  pnl_list_insert_first (L, PNL_OBJECT(M));
  pnl_list_insert_first (L, PNL_OBJECT(T));
  pnl_list_print (L);
  info = pnl_object_mpi_send (PNL_OBJECT(L), 1, SENDTAG, MPI_COMM_WORLD);
  printf ("List sent\n");
  pnl_list_free (&L);
  return info;
}

static int recv_list ()
{
  MPI_Status status;
  PnlList *L;
  PnlMat *M;
  PnlTridiagMat *T;
  int info;
  L = pnl_list_new ();
  M = pnl_mat_new ();
  T = pnl_tridiag_mat_new ();
  pnl_list_insert_last (L, PNL_OBJECT(T));
  pnl_list_insert_last (L, PNL_OBJECT(M));
  info = pnl_object_mpi_recv (PNL_OBJECT(L), 0, SENDTAG, MPI_COMM_WORLD, &status);
  printf ("Received list\n"); pnl_list_print (L); printf ("\n");
  pnl_list_free (&L);
  return info;
}

/*
 * Isend / Irecv examples
 */

static int isend_matrix ()
{
  MPI_Request request;
  MPI_Status status;
  PnlMat *M;
  int info, m = 4, n =5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  M = pnl_mat_new ();
  pnl_rand_init (gen, m, n);
  pnl_mat_rand_normal (M, m, n, gen);
  printf ("Original matrix (Isend)\n"); pnl_mat_print (M); printf ("\n");
  info = pnl_object_mpi_isend (PNL_OBJECT(M), 1, SENDTAG, MPI_COMM_WORLD, &request);
  info = MPI_Wait (&request, &status);
  pnl_mat_free (&M);
  return info;
}

static int status_is_empty (MPI_Status *status)
{
  if ( status->MPI_TAG == MPI_ANY_TAG && status->MPI_SOURCE == MPI_ANY_SOURCE
       && status->MPI_ERROR == MPI_SUCCESS )
    return TRUE;
  else return FALSE;
}

static int irecv_matrix ()
{
  MPI_Request request;
  MPI_Status status;
  void *buf;
  PnlMat *M;
  int flag,info, size, pos=0;
  M = pnl_mat_new ();
  /* Wait until a message is available for reception */
  do
    {
      info = pnl_object_mpi_irecv (&buf, &size, 0, SENDTAG, MPI_COMM_WORLD, &flag, &request);
    }
  while (flag == FALSE);
  PNL_MPI_MESSAGE (info, "error in pnl_object_mpi_irecv\n");
  MPI_Test (&request, &flag, &status);
  if ( status_is_empty (&status) == TRUE )
    {
      printf ("Emtpy status\n");
      return MPI_ERR_REQUEST;
    }
  if ( flag == FALSE )
    {
      info = MPI_Wait (&request, &status);
    }
  pnl_object_mpi_unpack (PNL_OBJECT(M), buf, size, &pos, MPI_COMM_WORLD);
  printf ("Received matrix (Irecv)\n"); pnl_mat_print (M); printf ("\n");
  pnl_mat_free (&M);
  free (buf);
  return info;
}

/*
 * Bcast examples
 */

static int bcast_matrix (PnlMat *M)
{
  int rank, info, m = 4, n =5;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  if ( rank == 0 )
    {
      pnl_rand_init (gen, m, n);
      pnl_mat_rand_normal (M, m, n, gen);
      printf ("Original matrix (Bcast)\n"); pnl_mat_print (M); printf ("\n");
    }
  info = pnl_object_mpi_bcast (PNL_OBJECT(M), 0, MPI_COMM_WORLD);
  printf ("Received matrix (Bcast, rank = %d)\n", rank); pnl_mat_print (M); printf ("\n");
  return info;
}

/*
 * TridiagMatLU 
 */

static PnlTridiagMat* create_random_tridiag (n, gen)
{
  PnlVect *dl, *du, *d;
  PnlTridiagMat *M;
  d = pnl_vect_create (n);
  du = pnl_vect_create (n);
  dl = pnl_vect_create (n);
  pnl_vect_rand_uni (d, n, 0., 1., gen);
  pnl_vect_rand_uni (du, n-1, 0., 1., gen);
  pnl_vect_rand_uni (dl, n-1, 0., 1., gen);
  M = pnl_tridiag_mat_create_from_ptr (n, dl->array, d->array, du->array);
  pnl_vect_free (&d);
  pnl_vect_free (&dl);
  pnl_vect_free (&du);
  return M;
}

static int test_tridiag_mat_lu ()
{
  PnlVect *b, *x;
  PnlTridiagMat *M;
  int rank;
  int n = 5;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  if ( rank == 0 )
    {
      PnlTridiagMatLU *LU;
      int gen = PNL_RNG_MERSENNE_RANDOM_SEED;

      pnl_rand_init (gen, 1, 1);
      x = pnl_vect_create (n);
      b = pnl_vect_create (n);
      pnl_vect_rand_normal (x, n, gen);
      M = create_random_tridiag (n, gen);
      pnl_tridiag_mat_mult_vect_inplace (b, M, x);
      LU = pnl_tridiag_mat_lu_new ();
      pnl_tridiag_mat_lu_compute (LU, M);
      
      pnl_object_mpi_send (PNL_OBJECT(LU), 1, SENDTAG, MPI_COMM_WORLD);
      pnl_object_mpi_send (PNL_OBJECT(b), 1, SENDTAG, MPI_COMM_WORLD);
      pnl_object_mpi_send (PNL_OBJECT(x), 1, SENDTAG, MPI_COMM_WORLD);

      pnl_vect_free (&x);
      pnl_vect_free (&b);
      pnl_tridiag_mat_free (&M);
      pnl_tridiag_mat_lu_free (&LU);
    }
  else
    {
      PnlVect *x_save;
      PnlTridiagMatLU *LU;
      MPI_Status status;

      x_save = pnl_vect_create (n);
      x = pnl_vect_create (n);
      b = pnl_vect_create (n);
      LU = pnl_tridiag_mat_lu_new ();
      
      pnl_object_mpi_recv (PNL_OBJECT(LU), 0, SENDTAG, MPI_COMM_WORLD, &status);
      pnl_object_mpi_recv (PNL_OBJECT(b), 0, SENDTAG, MPI_COMM_WORLD, &status);
      pnl_object_mpi_recv (PNL_OBJECT(x_save), 0, SENDTAG, MPI_COMM_WORLD, &status);
      pnl_tridiag_mat_lu_syslin (x, LU, b);
      if ( pnl_test_vect_eq_abs (x, x_save, 1E-12, "tridiag_mat_lu (MIP tests)", "") == TRUE )
        {
          printf ("PnlTridiagMatLU Send/Receive OK\n");
        }
      else
        {
          printf ("PnlTridiagMatLU Send/Receive FAIL\n");
        }

      pnl_vect_free (&x);
      pnl_vect_free (&x_save);
      pnl_tridiag_mat_lu_free (&LU);
    }
  return OK;
}

static void test_reduce (int rank)
{
  int n = 9;
  if ( rank == 0 )
    {
      PnlVect *sum, *vect1, *vect2, *reduc;
      PnlRng *rng = pnl_rng_create (PNL_RNG_MERSENNE);
      pnl_rng_sseed (rng, rank + 1234);
      reduc = pnl_vect_new ();
      vect1 = pnl_vect_new ();
      vect2 = pnl_vect_new ();

      pnl_vect_rng_normal (vect1, n, rng);
      pnl_vect_rng_normal (vect2, n, rng);

      sum = pnl_vect_copy (vect1);
      pnl_vect_plus_vect (sum, vect2);
      pnl_object_mpi_send (PNL_OBJECT(vect2), 1, SENDTAG, MPI_COMM_WORLD);
      pnl_object_mpi_reduce (PNL_OBJECT(vect1), PNL_OBJECT(reduc), MPI_SUM, 0, MPI_COMM_WORLD);

      if ( pnl_test_vect_eq (sum, reduc, 1E-10, "mpi_reduce", "") == TRUE )
        {
          printf ("MPI_Reduce for PnlObject: OK\n");
        }
      else
        {
          printf ("MPI_Reduce for PnlObject: FAIL\n");
        }
      pnl_vect_free (&reduc);
      pnl_vect_free (&sum);
      pnl_vect_free (&vect1);
      pnl_vect_free (&vect2);
    }
  else
    {
      MPI_Status status;
      PnlVect *vect2 = pnl_vect_new ();
      pnl_object_mpi_recv (PNL_OBJECT(vect2), 0, SENDTAG, MPI_COMM_WORLD, &status);
      pnl_object_mpi_reduce (PNL_OBJECT(vect2), NULL, MPI_SUM, 0, MPI_COMM_WORLD);
      pnl_vect_free (&vect2);
    }
}

int main(int argc, char *argv[])
{
  PnlMat *M;
  int rank, nproc;
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  pnl_message_off();
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
      init_int_array ();
      send_vector (); MPI_Barrier (MPI_COMM_WORLD);
      send_int_vector (); MPI_Barrier (MPI_COMM_WORLD);
      send_complex_vector (); MPI_Barrier (MPI_COMM_WORLD);
      send_matrix (); MPI_Barrier (MPI_COMM_WORLD);
      send_sp_matrix (); MPI_Barrier (MPI_COMM_WORLD);
      send_int_matrix (); MPI_Barrier (MPI_COMM_WORLD);
      send_complex_matrix (); MPI_Barrier (MPI_COMM_WORLD);
      send_bandmatrix (); MPI_Barrier (MPI_COMM_WORLD);
      send_tridiagmatrix (); MPI_Barrier (MPI_COMM_WORLD);
      send_int_hmatrix (); MPI_Barrier (MPI_COMM_WORLD);
      isend_matrix (); MPI_Barrier (MPI_COMM_WORLD);
      send_matrix (); MPI_Barrier (MPI_COMM_WORLD);
      send_basis (); MPI_Barrier (MPI_COMM_WORLD);
      ssend_basis (); MPI_Barrier (MPI_COMM_WORLD);
      send_list (); MPI_Barrier (MPI_COMM_WORLD);
    }
  else
    {
      recv_vector (); MPI_Barrier (MPI_COMM_WORLD);
      recv_int_vector (); MPI_Barrier (MPI_COMM_WORLD);
      recv_complex_vector (); MPI_Barrier (MPI_COMM_WORLD);
      recv_matrix (); MPI_Barrier (MPI_COMM_WORLD);
      recv_sp_matrix (); MPI_Barrier (MPI_COMM_WORLD);
      recv_int_matrix (); MPI_Barrier (MPI_COMM_WORLD);
      recv_complex_matrix (); MPI_Barrier (MPI_COMM_WORLD);
      recv_bandmatrix (); MPI_Barrier (MPI_COMM_WORLD);
      recv_tridiagmatrix (); MPI_Barrier (MPI_COMM_WORLD);
      recv_int_hmatrix (); MPI_Barrier (MPI_COMM_WORLD);
      recv_matrix (); MPI_Barrier (MPI_COMM_WORLD);
      irecv_matrix (); MPI_Barrier (MPI_COMM_WORLD);
      recv_basis (); MPI_Barrier (MPI_COMM_WORLD);
      recv_basis (); MPI_Barrier (MPI_COMM_WORLD);
      recv_list (); MPI_Barrier (MPI_COMM_WORLD);
    }

  M = pnl_mat_new ();
  bcast_matrix (M);
  pnl_mat_free (&M);
  MPI_Barrier (MPI_COMM_WORLD);
  test_tridiag_mat_lu ();
  MPI_Barrier (MPI_COMM_WORLD);
  test_reduce (rank);

  MPI_Finalize ();
  exit (0);
}
