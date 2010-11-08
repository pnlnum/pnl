
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
#include <mpi.h>
#include "pnl/pnl_object.h"
#include "pnl/pnl_tridiag_matrix.h"
#include "pnl/pnl_band_matrix.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_internals.h"
#include "pnl/pnl_list.h"
#include "pnl/pnl_basis.h"
#include "pnl/pnl_random.h"
#include "pnl/pnl_mpi.h"

#define PNL_MPI_MESSAGE(info, msg)              \
  if ( info != MPI_SUCCESS )                    \
    {                                           \
      if (pnl_message_is_on ()) printf (msg);   \
      return info;                              \
    }


/*
 * MPI_Pack_size wrappers
 */
static int size_vector (const PnlObject *Obj, MPI_Comm comm, int *size);
static int size_matrix (const PnlObject *Obj, MPI_Comm comm, int *size);
static int size_tridiag_matrix (const PnlObject *Obj, MPI_Comm comm, int *size);
static int size_band_matrix (const PnlObject *Obj, MPI_Comm comm, int *size);
static int size_hmatrix (const PnlObject *Obj, MPI_Comm comm, int *size);
static int size_basis (const PnlObject *Obj, MPI_Comm comm, int *size);
static int size_rng (const PnlObject *Obj, MPI_Comm comm, int *size);
static int size_list (const PnlObject *Obj, MPI_Comm comm, int *size);

/*
 * MPI_Pack wrappers
 */
static int pack_vector (const PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_matrix (const PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_tridiag_matrix (const PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_band_matrix (const PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_hmatrix (const PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_basis (const PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_rng (const PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_list (const PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm);

/*
 * MPI_Unpack wrappers
 */
static int unpack_vector (PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_matrix (PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_tridiag_matrix (PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_band_matrix (PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_hmatrix (PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_basis (PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_rng (PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_list (PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm);

/**
 * Computes the length of the buffer needed to pack the PnlObject
 *
 * @param Obj a PnlObject actually containing a PnlVectorObject
 * @param comm an MPI Communicator
 * @param size the upper bound on the number of bytes needed to pack Obj
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int size_vector (const PnlObject *Obj, MPI_Comm comm, int *size)
{
  int info, count, mn;
  MPI_Datatype t = MPI_DATATYPE_NULL;
  PnlVectObject *V = PNL_VECT_OBJECT(Obj);
  *size = 0;
  mn = V->size;
  /* V->size */
  if((info=MPI_Pack_size(1,MPI_INT, comm,&count))) return(info);
  *size += count;
  switch (PNL_GET_TYPE (Obj))
    {
    case PNL_TYPE_VECTOR_DOUBLE : t = MPI_DOUBLE;
      break;
    case PNL_TYPE_VECTOR_COMPLEX : t = MPI_DOUBLE; mn *= 2;
      break;
    case PNL_TYPE_VECTOR_INT : t = MPI_INT;
      break;
    case PNL_TYPE_VECTOR_COMPACT :
      PNL_ERROR ("Not implemented", "size_vector");
      break;
    }
  info=MPI_Pack_size(mn,t,comm,&count);
  *size += count;
  return (info);
}

/**
 * Computes the length of the buffer needed to pack the PnlObject
 *
 * @param Obj a PnlObject actually containing a PnlMatObject
 * @param comm an MPI Communicator
 * @param size the upper bound on the number of bytes needed to pack Obj
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int size_matrix (const PnlObject *Obj, MPI_Comm comm, int *size)
{
  int info, count, mn;
  MPI_Datatype t = MPI_DATATYPE_NULL;
  PnlMatObject *M = PNL_MAT_OBJECT(Obj);
  *size = 0;
  /* M->m */
  if((info=MPI_Pack_size(1,MPI_INT, comm,&count))) return(info);
  *size += count;
  /* M->n */
  if((info=MPI_Pack_size(1,MPI_INT, comm,&count))) return(info);
  *size += count;
  /* mn is not needed because it is computed inside the resize function */
  mn = M->mn;
  switch (PNL_GET_TYPE (Obj))
    {
    case PNL_TYPE_MATRIX_DOUBLE : t=MPI_DOUBLE;
      break;
    case PNL_TYPE_MATRIX_COMPLEX : t=MPI_DOUBLE; mn *= 2;
      break;
    case PNL_TYPE_MATRIX_INT : t=MPI_INT;
      break;
    }
  info=MPI_Pack_size(mn,t,comm,&count);
  *size += count;
  return(info);
}

/**
 * Computes the length of the buffer needed to pack the PnlObject
 *
 * @param Obj a PnlObject actually containing a PnlMatObject
 * @param comm an MPI Communicator
 * @param size the upper bound on the number of bytes needed to pack Obj
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int size_tridiag_matrix (const PnlObject *Obj, MPI_Comm comm, int *size)
{
  int info, count;
  PnlTridiagMatObject *M = PNL_TRIDIAGMAT_OBJECT(Obj);
  *size = 0;
  /* M->size */
  if((info=MPI_Pack_size(1,MPI_INT,comm,&count))) return(info);
  *size += count;
  if (PNL_GET_TYPE(M)!=PNL_TYPE_TRIDIAG_MATRIX_DOUBLE)
    {
      PNL_ERROR ("Unknown type", "pack_tridiag_matrix");
    }
  /* M->{D,DU,DL} */
  if((info=MPI_Pack_size(M->size,MPI_DOUBLE,comm,&count))) return(info);
  *size += count;
  if((info=MPI_Pack_size(M->size-1,MPI_DOUBLE,comm,&count))) return(info);
  *size += count;
  info=MPI_Pack_size(M->size-1,MPI_DOUBLE,comm,&count);
  *size += count;
  return(info);
}

/**
 * Computes the length of the buffer needed to pack the PnlObject
 *
 * @param Obj a PnlObject actually containing a PnlBandMatObject
 * @param comm an MPI Communicator
 * @param size the upper bound on the number of bytes needed to pack Obj
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int size_band_matrix (const PnlObject *Obj, MPI_Comm comm, int *size)
{
  int info, count, n;
  PnlBandMatObject *M = PNL_BAND_MAT_OBJECT(Obj);
  *size = 0;
  /* M->{m,n,nu,nl} */
  if((info=MPI_Pack_size(1,MPI_INT, comm,&count))) return(info);
  *size += count;
  if((info=MPI_Pack_size(1,MPI_INT, comm,&count))) return(info);
  *size += count;
  if((info=MPI_Pack_size(1,MPI_INT, comm,&count))) return(info);
  *size += count;
  if((info=MPI_Pack_size(1,MPI_INT, comm,&count))) return(info);
  *size += count;
  /* m_band and n_band are not needed because they are computed inside the
     resize function */
  if (PNL_GET_TYPE(M)!=PNL_TYPE_BAND_MATRIX_DOUBLE)
    {
      PNL_ERROR ("Unknown type", "size_band_matrix");
    }
  n = M->m_band * M->n_band;
  info=MPI_Pack_size(n,MPI_DOUBLE, comm,&count);
  *size += count;
  return(info);
}

/**
 * Computes the length of the buffer needed to pack the PnlObject
 *
 * @param Obj a PnlObject actually containing a PnlHmatObject
 * @param comm an MPI Communicator
 * @param size the upper bound on the number of bytes needed to pack Obj
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int size_hmatrix (const PnlObject *Obj, MPI_Comm comm, int *size)
{
  int info, count, mn;
  MPI_Datatype t = MPI_DATATYPE_NULL;
  PnlHmatObject *M = PNL_HMAT_OBJECT(Obj);
  *size = 0;
  /* M->ndim */
  if((info=MPI_Pack_size(1, MPI_INT, comm, &count))) return(info);
  *size += count;
  /* M->dims */
  if((info=MPI_Pack_size(M->ndim, MPI_INT, comm, &count))) return(info);
  *size += count;
  /* mn is not needed because it is computed inside the resize function */
  mn = M->mn;
  switch (PNL_GET_TYPE (Obj))
    {
    case PNL_TYPE_HMATRIX_DOUBLE : t=MPI_DOUBLE;
      break;
    case PNL_TYPE_HMATRIX_COMPLEX : t=MPI_DOUBLE; mn *= 2;
      break;
    case PNL_TYPE_HMATRIX_INT : t=MPI_INT;
      break;
    }
  info=MPI_Pack_size(mn,t,comm,&count);
  *size += count;
  return(info);
}

/**
 * Computes the length of the buffer needed to pack the PnlObject
 *
 * @param Obj a PnlObject actually containing a PnlBasi
 * @param comm an MPI Communicator
 * @param size the upper bound on the number of bytes needed to pack Obj
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int size_basis (const PnlObject *Obj, MPI_Comm comm, int *size)
{
  int info, count;
  PnlBasis *B = PNL_BASIS_OBJECT(Obj);
  *size = 0;
  /* B->id */
  if((info=MPI_Pack_size(1, MPI_INT, comm, &count))) return(info);
  *size += count;
  /* B->T */
  if((info=pnl_object_mpi_pack_size(PNL_OBJECT(B->T),comm,&count))) return info;
  *size += count;
  /* B->isreduced */
  if((info=MPI_Pack_size(1, MPI_INT, comm, &count))) return(info);
  *size += count;
  if (B->isreduced == 1)
    {
      if((info=MPI_Pack_size(B->nb_variates, MPI_DOUBLE, comm, &count))) return(info);
      *size += 2*count;
    }
  return (info);
}

/**
 * Computes the length of the buffer needed to pack the PnlObject
 *
 * @param Obj a PnlObject actually containing a PnlRng
 * @param comm an MPI Communicator
 * @param size the upper bound on the number of bytes needed to pack Obj
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int size_rng (const PnlObject *Obj, MPI_Comm comm, int *size)
{
  int info, count;
  PnlRng *rng = PNL_RNG_OBJECT(Obj);
  *size = 0;

  /* rng->type */
  if((info=MPI_Pack_size(1, MPI_INT, comm, &count))) return(info);
  *size += count;
  /* rng->dimension */
  if((info=MPI_Pack_size(1, MPI_INT, comm, &count))) return(info);
  *size += count;
  /* rng->counter */
  if((info=MPI_Pack_size(1, MPI_INT, comm, &count))) return(info);
  *size += count;
  /* rng->has_gauss */
  if((info=MPI_Pack_size(1, MPI_INT, comm, &count))) return(info);
  *size += count;
  if (rng->has_gauss == 1)
    {
      /* rng->gauss */
      if((info=MPI_Pack_size(1, MPI_DOUBLE, comm, &count))) return(info);
      *size += count;
    }
  /* rng->state */
  if((info=pnl_rng_state_mpi_pack_size(rng, comm, &count))) return(info);
  *size += count;
  return (info);
}

/**
 * Computes the length of the buffer needed to pack the PnlObject
 *
 * @param Obj a PnlObject actually containing a PnlList
 * @param comm an MPI Communicator
 * @param size the upper bound on the number of bytes needed to pack Obj
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int size_list (const PnlObject *Obj, MPI_Comm comm, int *size)
{
  int info, count, i;
  PnlCell *C;
  PnlList *L = PNL_LIST_OBJECT(Obj);
  
  *size = 0;

  /* L->len */
  if((info=MPI_Pack_size(1, MPI_INT, comm, &count))) return(info);
  *size += count;

  /* Compute length for each element */
  C = L->first;
  for ( i=0 ; i<L->len ; i++ )
    {
      /* we store PNL_GET_TYPE(C->self) twice for unpacking easier */
      if((info=MPI_Pack_size(1, MPI_INT, comm, &count))) return(info);
      *size += count;
      if ((pnl_object_mpi_pack_size (C->self, comm, &count))) return(info);
      *size += count;
      C = C->next;
    }
  return (info);
}

/**
 * Packs a PnlVectObject
 *
 * @param Obj a PnlObject containing a PnlVectObject
 * @param buf an already allocated buffer of length bufsize used to store Obj
 * @param bufsize the size of the buffer. It must be at least equal to the value
 * computed by size_vector
 * @param comm an MPI Communicator
 * @param pos (in/out) the current position in buf
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int pack_vector (const PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int mn,info;
  MPI_Datatype t = MPI_DATATYPE_NULL;
  PnlVectObject *V = PNL_VECT_OBJECT(Obj);
  if ((info=MPI_Pack(&(V->size),1,MPI_INT,buf,bufsize,pos,comm))) return info;
  mn = V->size;
  switch (PNL_GET_TYPE (V))
    {
    case PNL_TYPE_VECTOR_DOUBLE : t = MPI_DOUBLE;
      break;
    case PNL_TYPE_VECTOR_COMPLEX : mn *= 2; t = MPI_DOUBLE;
      break;
    case PNL_TYPE_VECTOR_INT : t= MPI_INT;
      break;
    case PNL_TYPE_VECTOR_COMPACT :
      PNL_ERROR ("Not implemented", "pack_vector");
      break;
    }
  info=MPI_Pack(V->array,mn,t,buf,bufsize,pos,comm);
  return(info);
}

/**
 * Packs a PnlMatObject
 *
 * @param Obj a PnlObject containing a PnlMatObject
 * @param buf an already allocated buffer of length bufsize used to store Obj
 * @param bufsize the size of the buffer. It must be at least equal to the value
 * computed by size_matrix
 * @param comm an MPI Communicator
 * @param pos (in/out) the current position in buf
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int pack_matrix (const PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int mn,info;
  MPI_Datatype t = MPI_DATATYPE_NULL;
  PnlMatObject *M = PNL_MAT_OBJECT(Obj);
  if ((info=MPI_Pack(&(M->m), 1, MPI_INT, buf, bufsize, pos, comm))) return info;
  if ((info=MPI_Pack(&(M->n), 1, MPI_INT, buf, bufsize, pos, comm))) return info;
  /* mn is not needed because it is computed inside the resize function */
  mn = M->mn;
  switch (PNL_GET_TYPE (M))
    {
    case PNL_TYPE_MATRIX_DOUBLE : t = MPI_DOUBLE;
      break;
    case PNL_TYPE_MATRIX_COMPLEX : mn *= 2; t = MPI_DOUBLE;
      break;
    case PNL_TYPE_MATRIX_INT : t = MPI_INT;
      break;
    }
  info=MPI_Pack(M->array,mn,t,buf,bufsize,pos,comm);
  return(info);
}

/**
 * Packs a PnlTridiagMatObject
 *
 * @param Obj a PnlObject containing a PnlTridiagMatObject
 * @param buf an already allocated buffer of length bufsize used to store Obj
 * @param bufsize the size of the buffer. It must be at least equal to the value
 * computed by size_matrix
 * @param comm an MPI Communicator
 * @param pos (in/out) the current position in buf
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int pack_tridiag_matrix (const PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int n,info;
  PnlTridiagMatObject *M = PNL_TRIDIAGMAT_OBJECT(Obj);
  n = M->size;
  if ((info=MPI_Pack(&n, 1, MPI_INT, buf, bufsize, pos, comm))) return info;
  if (PNL_GET_TYPE(M)!=PNL_TYPE_TRIDIAG_MATRIX_DOUBLE)
    {
      PNL_ERROR ("Unknown type", "pack_tridiag_matrix");
    }
  if((info=MPI_Pack(M->D,n,MPI_DOUBLE,buf,bufsize,pos,comm)))return(info);
  if((info=MPI_Pack(M->DL,n-1,MPI_DOUBLE,buf,bufsize,pos,comm)))return(info);
  if((info=MPI_Pack(M->DU,n-1,MPI_DOUBLE,buf,bufsize,pos,comm)))return(info);
  return(info);
}

/**
 * Packs a PnlBandMatObject
 *
 * @param Obj a PnlObject containing a PnlBandMatObject
 * @param buf an already allocated buffer of length bufsize used to store Obj
 * @param bufsize the size of the buffer. It must be at least equal to the value
 * computed by size_matrix
 * @param comm an MPI Communicator
 * @param pos (in/out) the current position in buf
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int pack_band_matrix (const PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  PnlBandMatObject *M = PNL_BAND_MAT_OBJECT(Obj);
  if ((info=MPI_Pack(&M->m, 1, MPI_INT, buf, bufsize, pos, comm))) return info;
  if ((info=MPI_Pack(&M->n, 1, MPI_INT, buf, bufsize, pos, comm))) return info;
  if ((info=MPI_Pack(&M->nu, 1, MPI_INT, buf, bufsize, pos, comm))) return info;
  if ((info=MPI_Pack(&M->nl, 1, MPI_INT, buf, bufsize, pos, comm))) return info;
  /* m_band and n_band are not needed because they are computed inside the
     resize function */
  if (PNL_GET_TYPE(M)!=PNL_TYPE_BAND_MATRIX_DOUBLE)
    {
      PNL_ERROR ("Unknown type", "pack_band_matrix");
    }
  if((info=MPI_Pack(M->array,M->m_band*M->n_band,MPI_DOUBLE,buf,bufsize,pos,comm)))return(info);
  return(info);
}

/**
 * Packs a PnlHmatObject
 *
 * @param Obj a PnlObject containing a PnlHmatObject
 * @param buf an already allocated buffer of length bufsize used to store Obj
 * @param bufsize the size of the buffer. It must be at least equal to the value
 * computed by size_matrix
 * @param comm an MPI Communicator
 * @param pos (in/out) the current position in buf
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int pack_hmatrix (const PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  MPI_Datatype t = MPI_DATATYPE_NULL;
  int mn, info;
  PnlHmatObject *M = PNL_HMAT_OBJECT(Obj);
  if ((info=MPI_Pack(&M->ndim, 1, MPI_INT, buf, bufsize, pos, comm))) return info;
  if ((info=MPI_Pack(M->dims, M->ndim, MPI_INT, buf, bufsize, pos, comm))) return info;
  /* mn is not needed because it is computed inside the resize function */
  mn = M->mn;
  switch (PNL_GET_TYPE (M))
    {
    case PNL_TYPE_HMATRIX_DOUBLE :  t = MPI_DOUBLE;
      break;
    case PNL_TYPE_HMATRIX_COMPLEX : mn *= 2; t = MPI_DOUBLE;
      break;
    case PNL_TYPE_HMATRIX_INT : t = MPI_INT;
      break;
    }
  if((info=MPI_Pack(M->array,mn,t,buf,bufsize,pos,comm)))return(info);
  return(info);
}

/**
 * Packs a PnlBasis
 *
 * @param Obj a PnlObject containing a PnlBasis
 * @param buf an already allocated buffer of length bufsize used to store Obj
 * @param bufsize the size of the buffer. It must be at least equal to the value
 * computed by size_matrix
 * @param comm an MPI Communicator
 * @param pos (in/out) the current position in buf
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int pack_basis (const PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  PnlBasis *B = PNL_BASIS_OBJECT(Obj);
  if ((info=MPI_Pack(&B->id, 1, MPI_INT, buf, bufsize, pos, comm))) return info;
  if ((info=pnl_object_mpi_pack(PNL_OBJECT(B->T), buf, bufsize, pos, comm))) return info;
  if ((info=MPI_Pack(&B->isreduced, 1, MPI_INT, buf, bufsize, pos, comm))) return info;
  if (B->isreduced == 1)
    {
      const int n = B->nb_variates;
      if ((info=MPI_Pack(&B->center, n, MPI_DOUBLE, buf, bufsize, pos, comm))) return info;
      if ((info=MPI_Pack(&B->scale, n, MPI_DOUBLE, buf, bufsize, pos, comm))) return info;
    }
  return(info);
}

/**
 * Packs a PnlRng
 *
 * @param Obj a PnlObject containing a PnlRng
 * @param buf an already allocated buffer of length bufsize used to store Obj
 * @param bufsize the size of the buffer. It must be at least equal to the value
 * computed by size_matrix
 * @param comm an MPI Communicator
 * @param pos (in/out) the current position in buf
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int pack_rng (const PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  PnlRng *rng = PNL_RNG_OBJECT(Obj);
  
  if ((info=MPI_Pack(&rng->type, 1, MPI_INT, buf, bufsize, pos, comm))) return info;
  if ((info=MPI_Pack(&rng->dimension, 1, MPI_INT, buf, bufsize, pos, comm))) return info;
  if ((info=MPI_Pack(&rng->counter, 1, MPI_INT, buf, bufsize, pos, comm))) return info;
  if ((info=MPI_Pack(&rng->has_gauss, 1, MPI_INT, buf, bufsize, pos, comm))) return info;
  if (rng->has_gauss == 1)
    {
      if ((info=MPI_Pack(&rng->gauss, 1, MPI_DOUBLE, buf, bufsize, pos, comm))) return info;
    }
  if ((info=pnl_rng_state_mpi_pack(rng, buf, bufsize, pos, comm))) return info;
  return (info);
}

/**
 * Packs a PnlList
 *
 * @param Obj a PnlObject containing a PnlList
 * @param buf an already allocated buffer of length bufsize used to store Obj
 * @param bufsize the size of the buffer. It must be at least equal to the value
 * computed by size_matrix
 * @param comm an MPI Communicator
 * @param pos (in/out) the current position in buf
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int pack_list (const PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info, i;
  PnlCell *C;
  PnlList *L = PNL_LIST_OBJECT(Obj);
  
  if ((info=MPI_Pack(&L->len, 1, MPI_INT, buf, bufsize, pos, comm))) return info;
  C = L->first;
  for ( i=0 ; i<L->len ; i++ )
    {
      if ((info=MPI_Pack(&(PNL_GET_TYPE(C->self)), 1, MPI_INT, buf, bufsize, pos, comm))) return info;
      if ((info=pnl_object_mpi_pack(C->self, buf, bufsize, pos, comm))) return info;
      C = C->next;
    }
  return (info);
}

/**
 * Unpacks a PnlVectObject
 *
 * @param Obj a PnlObject containing a PnlVectObject
 * @param buf an already allocated buffer of length bufsize used to store Obj
 * @param bufsize the size of the buffer. It must be at least equal to the value
 * computed by size_vector
 * @param comm an MPI Communicator
 * @param pos (in/out) the current position in buf
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int unpack_vector (PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int n, id, info;
  MPI_Datatype t = MPI_DATATYPE_NULL;
  PnlVectObject *V = PNL_VECT_OBJECT(Obj);

  if ((info=MPI_Unpack(buf,bufsize,pos,&id,1,MPI_INT,comm))) return info;
  if ( id != PNL_GET_TYPE (V) )
    {
      printf ("Wrong type of object in unpack_vector.\n");
      return MPI_ERR_TYPE;
    }
  if ((info=MPI_Unpack(buf,bufsize,pos,&n,1,MPI_INT,comm))) return info;
  pnl_vect_object_resize(V, n);
  switch (PNL_GET_TYPE (V))
    {
    case PNL_TYPE_VECTOR_DOUBLE : t = MPI_DOUBLE;
      break;
    case PNL_TYPE_VECTOR_COMPLEX : n *= 2; t = MPI_DOUBLE;
      break;
    case PNL_TYPE_VECTOR_INT : t = MPI_INT;
      break;
    case PNL_TYPE_VECTOR_COMPACT :
      printf ("Packing is not implemented for %s\n",PNL_GET_TYPENAME(V)); abort ();
      break;
    }
  info=MPI_Unpack(buf,bufsize,pos,V->array,n,t,comm);
  return info;
}

/**
 * Unpacks a PnlMatObject
 *
 * @param Obj a PnlObject containing a PnlMatObject
 * @param buf an already allocated buffer of length bufsize used to store Obj
 * @param bufsize the size of the buffer. It must be at least equal to the value
 * computed by size_vector
 * @param comm an MPI Communicator
 * @param pos (in/out) the current position in buf
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int unpack_matrix (PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int m, n, mn, id, info;
  MPI_Datatype t = MPI_DATATYPE_NULL;
  PnlMatObject *M = PNL_MAT_OBJECT(Obj);

  if ((info=MPI_Unpack(buf,bufsize,pos,&id,1,MPI_INT,comm))) return info;
  if ( id != PNL_GET_TYPE (M) )
    {
      printf ("Wrong type of object in unpack_matrix.\n");
      return MPI_ERR_TYPE;
    }
  if ((info=MPI_Unpack(buf,bufsize,pos,&m,1,MPI_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&n,1,MPI_INT,comm))) return info;
  pnl_mat_object_resize(M, m, n);
  mn = M->mn;
  switch (PNL_GET_TYPE (M))
    {
    case PNL_TYPE_MATRIX_DOUBLE : t = MPI_DOUBLE;
      break;
    case PNL_TYPE_MATRIX_COMPLEX : mn *= 2; t = MPI_DOUBLE;
      break;
    case PNL_TYPE_MATRIX_INT : t = MPI_INT;
      break;
    }
  if((info=MPI_Unpack(buf,bufsize,pos,M->array,mn,t,comm)))return(info);
  return info;
}

/**
 * Unpacks a PnlTridiagMatObject
 *
 * @param Obj a PnlObject containing a PnlTridiagMatObject
 * @param buf an already allocated buffer of length bufsize used to store Obj
 * @param bufsize the size of the buffer. It must be at least equal to the value
 * computed by size_vector
 * @param comm an MPI Communicator
 * @param pos (in/out) the current position in buf
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int unpack_tridiag_matrix (PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int id,n,info;
  PnlTridiagMat *M = (PnlTridiagMat *) Obj;
  if ((info=MPI_Unpack(buf,bufsize,pos,&id,1,MPI_INT,comm))) return info;
  if (PNL_GET_TYPE(M)!=PNL_TYPE_TRIDIAG_MATRIX_DOUBLE)
    {
      PNL_ERROR ("Unknown type", "unpack_tridiag_matrix");
    }
  if ((info=MPI_Unpack(buf,bufsize,pos,&n, 1, MPI_INT, comm))) return info;
  pnl_tridiag_mat_resize (M, n);
  if((info=MPI_Unpack(buf,bufsize,pos,M->D,n,MPI_DOUBLE,comm)))return(info);
  if((info=MPI_Unpack(buf,bufsize,pos,M->DL,n-1,MPI_DOUBLE,comm)))return(info);
  if((info=MPI_Unpack(buf,bufsize,pos,M->DU,n-1,MPI_DOUBLE,comm)))return(info);
  return(info);
}

/**
 * Unpacks a PnlBandMatObject
 *
 * @param Obj a PnlObject containing a PnlBandMatObject
 * @param buf an already allocated buffer of length bufsize used to store Obj
 * @param bufsize the size of the buffer. It must be at least equal to the value
 * computed by size_vector
 * @param comm an MPI Communicator
 * @param pos (in/out) the current position in buf
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int unpack_band_matrix (PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int id,info;
  int m, n, nu, nl;
  PnlBandMat *M = (PnlBandMat *) Obj;
  if ((info=MPI_Unpack(buf,bufsize,pos,&id,1,MPI_INT,comm))) return info;
  if (PNL_GET_TYPE(M)!=PNL_TYPE_BAND_MATRIX_DOUBLE)
    {
      PNL_ERROR ("Unknown type", "pack_band_matrix");
    }

  if ((info=MPI_Unpack(buf,bufsize,pos,&m, 1, MPI_INT, comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&n, 1, MPI_INT, comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&nu, 1, MPI_INT, comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&nl, 1, MPI_INT, comm))) return info;
  pnl_band_mat_resize (M, m, n, nl, nu);
  if((info=MPI_Unpack(buf,bufsize,pos,M->array,M->m_band*M->n_band,MPI_DOUBLE,comm)))return(info);
  return(info);
}

/**
 * Unpacks a PnlHmatObject
 *
 * @param Obj a PnlObject containing a PnlHmatObject
 * @param buf an already allocated buffer of length bufsize used to store Obj
 * @param bufsize the size of the buffer. It must be at least equal to the value
 * computed by size_vector
 * @param comm an MPI Communicator
 * @param pos (in/out) the current position in buf
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int unpack_hmatrix (PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int id, info, ndim, *dims, mn;
  MPI_Datatype t = MPI_DATATYPE_NULL;
  PnlHmatObject *M = PNL_HMAT_OBJECT(Obj);
  if ((info=MPI_Unpack(buf,bufsize,pos,&id,1,MPI_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&ndim, 1, MPI_INT, comm))) return info;

  dims = MALLOC_INT (ndim);
  if ((info=MPI_Unpack(buf,bufsize,pos,dims,ndim,MPI_INT,comm))) return info;
  pnl_hmat_object_resize (M, ndim, dims);
  switch (PNL_GET_TYPE (M))
    {
    case PNL_TYPE_HMATRIX_DOUBLE : mn = M->mn; t = MPI_DOUBLE;
      break;
    case PNL_TYPE_HMATRIX_COMPLEX : mn = 2*M->mn; t = MPI_DOUBLE;
      break;
    case PNL_TYPE_HMATRIX_INT : mn = M->mn; t = MPI_INT;
      break;
    }
  if((info=MPI_Unpack(buf,bufsize,pos,M->array,M->mn,t,comm)))return(info);
  return(info);
}

/**
 * Unpacks a PnlBasis
 *
 * @param Obj a PnlObject containing a PnlBasis
 * @param buf an already allocated buffer of length bufsize used to store Obj
 * @param bufsize the size of the buffer. It must be at least equal to the value
 * computed by size_vector
 * @param comm an MPI Communicator
 * @param pos (in/out) the current position in buf
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int unpack_basis (PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int id, info, basis_id;
  PnlBasis *B = PNL_BASIS_OBJECT(Obj);
  PnlMatInt *T = pnl_mat_int_new();
  if ((info=MPI_Unpack(buf,bufsize,pos,&id,1,MPI_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&basis_id,1,MPI_INT,comm))) return info;
  if ((info=pnl_object_mpi_unpack(PNL_OBJECT(T),buf,bufsize,pos,comm))) return info;
  pnl_basis_set_from_tensor (B, basis_id, T);
  /* do not free T, it will be done by pnl_basis_free later */
  if ((info=MPI_Unpack(buf,bufsize,pos,&B->isreduced,1,MPI_INT,comm))) return info;
  if (B->isreduced == 1)
    {
      const int n = T->n;
      B->center = malloc (n * sizeof(double));
      if ((info=MPI_Unpack(buf,bufsize,pos,B->center,n,MPI_DOUBLE,comm))) return info;
      B->scale = malloc (n * sizeof(double));
      if ((info=MPI_Unpack(buf,bufsize,pos,B->scale,n,MPI_DOUBLE,comm))) return info;
    }
  return (MPI_SUCCESS);
}

/**
 * Unpacks a PnlRng
 *
 * @param Obj a PnlObject containing a PnlRng
 * @param buf an already allocated buffer of length bufsize used to store Obj
 * @param bufsize the size of the buffer. It must be at least equal to the value
 * computed by size_vector
 * @param comm an MPI Communicator
 * @param pos (in/out) the current position in buf
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int unpack_rng (PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int type, info, id;
  PnlRng *rng = PNL_RNG_OBJECT(Obj);

  /* unpacking Obj->object.id */
  if ((info=MPI_Unpack(buf,bufsize,pos,&id,1,MPI_INT,comm))) return info;
  
  if ((info=MPI_Unpack(buf,bufsize,pos,&type,1,MPI_INT,comm))) return info;
  if ( rng->state != NULL )
    {
      free (rng->state); rng->state = NULL;
    }
  pnl_rng_init (rng, type);
  
  if ((info=MPI_Unpack(buf,bufsize,pos,&rng->dimension,1,MPI_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&rng->counter,1,MPI_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&rng->has_gauss,1,MPI_INT,comm))) return info;
  if (rng->has_gauss == 1)
    {
      if ((info=MPI_Unpack(buf,bufsize,pos,&rng->gauss,1,MPI_DOUBLE,comm))) return info;
    }
  if ((info=pnl_rng_state_mpi_unpack(rng,buf,bufsize,pos,comm))) return info;
  return (info);
}

/**
 * Unpacks a PnlList
 *
 * If the length of the list > 0, we try to unpack the buffer into the objects
 * already stored in list. An error occurs if the parent types of the stored
 * objects do not match the ones of the packed objects.
 *
 * @param Obj a PnlObject containing a PnlList
 * @param buf an already allocated buffer of length bufsize used to store Obj
 * @param bufsize the size of the buffer. It must be at least equal to the value
 * computed by size_vector
 * @param comm an MPI Communicator
 * @param pos (in/out) the current position in buf
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
static int unpack_list (PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info, len, i, id;
  PnlList *L = PNL_LIST_OBJECT(Obj);

  /* unpacking Obj->object.id */
  if ((info=MPI_Unpack(buf,bufsize,pos,&id,1,MPI_INT,comm))) return info;
  /* length of packed list */
  if ((info=MPI_Unpack(buf,bufsize,pos,&len,1,MPI_INT,comm))) return info;

  /*
   * if the list is already of the correct length, we assume that we want to
   * unpack buf into the objects already stored in the list. It is only
   * possible if the parent types of the stored objects match the ones the
   * packed objects
   */
  if ( len == L->len )
    {
      PnlCell *C;
      C = L->first;
      for ( i=0 ; i<len ; i++ )
        {
          int ptype;
          if ((info=MPI_Unpack(buf,bufsize,pos,&ptype,1,MPI_INT,comm))) return info;
          if ((info=pnl_object_mpi_unpack (C->self, buf, bufsize, pos, comm))) return info;
          C = C->next;
        }
      return info;
    }
  if ( L->len != 0 && len != L->len )
    {
      PNL_ERROR ( "size_mismatched", "unpack_list");
    }
  for ( i=0 ; i<len ; i++ )
    {
      PnlObject *O;
      int ptype;
      if ((info=MPI_Unpack(buf,bufsize,pos,&ptype,1,MPI_INT,comm))) return info;
      O = pnl_object_create (ptype);
      if ((info=pnl_object_mpi_unpack (O, buf, bufsize, pos, comm))) return info;
      pnl_list_insert_last (L, O);
    }
  return (info);
}


/*
 * Exported wrappers for handling PnlObjects
 */


/**
 * Computes the length of the buffer needed to pack the PnlObject and add this
 * value to size.
 *
 * @param Obj a PnlObject
 * @param comm an MPI Communicator
 * @param size the upper bound on the number of bytes needed to pack Obj
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
int pnl_object_mpi_pack_size (const PnlObject *Obj, MPI_Comm comm, int *size)
{
  int id, count, info;
  *size = 0;
  /* space for type and parent_type */
  if ((info=MPI_Pack_size(1,MPI_INT,comm,&count))) return info;
  *size += count;
  if ((info=MPI_Pack_size(1,MPI_INT,comm,&count))) return info;
  *size += count;

  id = PNL_GET_PARENT_TYPE(Obj);
  switch (id)
    {
    case PNL_TYPE_VECTOR:
      info = size_vector (Obj, comm, &count);
      break;
    case PNL_TYPE_MATRIX:
      info = size_matrix (Obj, comm, &count);
      break;
    case PNL_TYPE_TRIDIAG_MATRIX:
      info = size_tridiag_matrix (Obj, comm, &count);
      break;
    case PNL_TYPE_BAND_MATRIX:
      info = size_band_matrix (Obj, comm, &count);
      break;
    case PNL_TYPE_HMATRIX:
      info = size_hmatrix (Obj, comm, &count);
      break;
    case PNL_TYPE_BASIS:
      info = size_basis (Obj, comm, &count);
      break;
    case PNL_TYPE_RNG:
      info = size_rng (Obj, comm, &count);
      break;
    case PNL_TYPE_LIST:
      info = size_list (Obj, comm, &count);
      break;
    default:
      printf("Computing size of %s is not implemented yet", Obj->label);
      return MPI_ERR_TYPE;
    }
  *size += count; return info;
}


/**
 * Packs a PnlOject
 *
 * @param Obj a PnlObject
 * @param buf an already allocated buffer of length bufsize used to store Obj
 * @param bufsize the size of the buffer. It must be at least equal to the value
 * computed by pnl_object_mpi_Pack_size
 * @param comm an MPI Communicator
 * @param pos (in/out) the current position in buf
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
int pnl_object_mpi_pack (const PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int id, parent_id, info;
  parent_id = PNL_GET_PARENT_TYPE(Obj);
  if ((info=MPI_Pack(&parent_id, 1, MPI_INT, buf, bufsize, pos, comm))) return info;
  id = PNL_GET_TYPE(Obj);
  if ((info=MPI_Pack(&id, 1, MPI_INT, buf, bufsize, pos, comm))) return info;
  switch (parent_id)
    {
    case PNL_TYPE_VECTOR:
      return pack_vector (Obj, buf, bufsize, pos, comm);
      break;
    case PNL_TYPE_MATRIX:
      return pack_matrix (Obj, buf, bufsize, pos, comm);
      break;
    case PNL_TYPE_TRIDIAG_MATRIX:
      return pack_tridiag_matrix (Obj, buf, bufsize, pos, comm);
      break;
    case PNL_TYPE_BAND_MATRIX:
      return pack_band_matrix (Obj, buf, bufsize, pos, comm);
      break;
    case PNL_TYPE_HMATRIX:
      return pack_hmatrix (Obj, buf, bufsize, pos, comm);
      break;
    case PNL_TYPE_BASIS:
      return pack_basis (Obj, buf, bufsize, pos, comm);
      break;
    case PNL_TYPE_RNG:
      return pack_rng (Obj, buf, bufsize, pos, comm);
      break;
    case PNL_TYPE_LIST:
      return pack_list (Obj, buf, bufsize, pos, comm);
      break;
    default:
      printf("Packing for type %s is not implemented yet", Obj->label);
      return MPI_ERR_TYPE;
    }
}

/**
 * Unpacks a PnlOject
 *
 * @param Obj a PnlObject
 * @param buf a buffer created by pnl_object_mpi_Pack of length bufsize used
 * to store Obj
 * @param bufsize the size of the buffer. It must be at least equal to the value
 * computed by pnl_object_mpi_Pack_size
 * @param comm an MPI Communicator
 * @param pos (in/out) the current position in buf
 *
 * @return an error value equal to MPI_SUCCESS when everything is OK
 */
int pnl_object_mpi_unpack (PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int parent_id, info;

  if ((info=MPI_Unpack(buf, bufsize, pos, &parent_id, 1, MPI_INT, comm))) return info;
  if (parent_id != PNL_GET_PARENT_TYPE(Obj))
    {
      printf("Expected object type does not match the type of the received one.\n");
      return MPI_ERR_TYPE;
    }

  switch (parent_id)
    {
    case PNL_TYPE_VECTOR:
      return unpack_vector (Obj, buf, bufsize, pos, comm);
      break;
    case PNL_TYPE_MATRIX:
      return unpack_matrix (Obj, buf, bufsize, pos, comm);
      break;
    case PNL_TYPE_TRIDIAG_MATRIX:
      return unpack_tridiag_matrix (Obj, buf, bufsize, pos, comm);
      break;
    case PNL_TYPE_BAND_MATRIX:
      return unpack_band_matrix (Obj, buf, bufsize, pos, comm);
      break;
    case PNL_TYPE_HMATRIX:
      return unpack_hmatrix (Obj, buf, bufsize, pos, comm);
      break;
    case PNL_TYPE_BASIS:
      return unpack_basis (Obj, buf, bufsize, pos, comm);
      break;
    case PNL_TYPE_RNG:
      return unpack_rng (Obj, buf, bufsize, pos, comm);
      break;
    case PNL_TYPE_LIST:
      return unpack_list (Obj, buf, bufsize, pos, comm);
      break;
    default:
      printf("Unpacking for type %s is not implemented yet", Obj->label);
      return MPI_ERR_TYPE;
    }
}

/**
 * Performs a blocking send of a PnlObject
 *
 * @param Obj a PnlObject to send
 * @param dest the rank of destination (integer)
 * @param tag the message tag (integer)
 * @param comm a Communicator (handle)
 *
 * @return an error value, it should be MPI_SUCCESS when everything went OK
 */
int pnl_object_mpi_send (const PnlObject *Obj, int dest, int tag, MPI_Comm comm)
{
  int info, size=0, pos=0;
  char *buf = NULL;
  info = pnl_object_mpi_pack_size (Obj, comm, &size);
  PNL_MPI_MESSAGE(info, "error in computing size.\n");
  if ((buf = malloc(size)) == NULL) return MPI_ERR_BUFFER;
  info = pnl_object_mpi_pack(Obj, buf, size, &pos, comm);
  PNL_MPI_MESSAGE( info, "error in packing.\n");
  info = MPI_Send (buf, pos, MPI_PACKED, dest, tag, comm);
  free (buf);
  return info;
}

/**
 * Performs a blocking standard synchronous send of a PnlObject
 *
 * @param Obj a PnlObject to send
 * @param dest the rank of destination (integer)
 * @param tag the message tag (integer)
 * @param comm a Communicator (handle)
 *
 * @return an error value, it should be MPI_SUCCESS when everything went OK
 */
int pnl_object_mpi_ssend (const PnlObject *Obj, int dest, int tag, MPI_Comm comm)
{
  int info, size=0, pos=0;
  char *buf = NULL;
  info = pnl_object_mpi_pack_size (Obj, comm, &size);
  PNL_MPI_MESSAGE(info, "error in computing size.\n");
  if ((buf = malloc(size)) == NULL) return MPI_ERR_BUFFER;
  info = pnl_object_mpi_pack(Obj, buf, size, &pos, comm);
  PNL_MPI_MESSAGE( info, "error in packing.\n");
  info = MPI_Ssend (buf, pos, MPI_PACKED, dest, tag, comm);
  free (buf);
  return info;
}


/**
 * Performs a blocking receive of a PnlObject
 *
 * @param Obj a PnlObject used to store the received object. It must have
 * already been allocated by a call to the proper pnl_xxx_new function. Note
 * that the true type of Obj must match the type of the received object.
 * @param src the rank of source (integer)
 * @param tag the message tag (integer)
 * @param comm a Communicator (handle)
 * @param status a status object (status)
 *
 * @return an error value, it should be MPI_SUCCESS when everything went OK
 */
int pnl_object_mpi_recv (PnlObject *Obj, int src, int tag, MPI_Comm comm, MPI_Status *status)
{
  int info, size=0, pos=0;
  char *buf = NULL;
  info = MPI_Probe (src, tag, comm, status);
  PNL_MPI_MESSAGE( info, "error in Probe.\n");
  info = MPI_Get_count (status, MPI_PACKED, &size);
  PNL_MPI_MESSAGE( info, "error in Get_count.\n");
  buf = malloc(size);
  info = MPI_Recv (buf, size, MPI_PACKED, src, tag, comm, status);
  PNL_MPI_MESSAGE( info, "error in MPI_Recv.\n");
  info = pnl_object_mpi_unpack(Obj, buf, size, &pos, comm);
  free (buf); buf=NULL;
  return info;
}


/**
 * Broadcasts  a PnlObject from the process with rank root to
 * all other processes of the group.
 * Note that this function calls MPI_Bcast twice because before receiving a
 * packed PnlObject every process must know the size of the packed object.
 * One call to Bcast for the size and an other one for truly sending the object.
 *
 * @param Obj a PnlObject to broadcast
 * @param root the rank of the root process (integer)
 * @param comm a Communicator (handle)
 *
 * @return an error value, it should be MPI_SUCCESS when everything went OK
 */
int pnl_object_mpi_bcast (PnlObject *Obj, int root, MPI_Comm comm)
{
  int rank, info, size=0, pos=0;
  char *buf = NULL;
  MPI_Comm_rank (comm, &rank);
  if ( rank == 0 )
    {
      info = pnl_object_mpi_pack_size (Obj, comm, &size);
      PNL_MPI_MESSAGE(info, "error in computing size.\n");
      if ((buf = malloc(size)) == NULL) return MPI_ERR_BUFFER;
      info = pnl_object_mpi_pack(Obj, buf, size, &pos, comm);
      PNL_MPI_MESSAGE( info, "error in packing.\n");
    }
  /* Send the size of the packed object to everybody */
  MPI_Bcast (&size, 1, MPI_INT, root, comm);
  if ( rank != 0 ) buf = malloc(size);
  info = MPI_Bcast (buf, size, MPI_PACKED, root, comm);
  PNL_MPI_MESSAGE( info, "error in MPI_Bcast.\n");
  pos = 0;
  pnl_object_mpi_unpack(Obj, buf, size, &pos, comm);
  if ( rank != 0 ) free (buf); buf=NULL;
  return info;
}

/**
 * Performs a non-blocking send of a PnlObject
 *
 * @param Obj a PnlObject to send
 * @param dest the rank of destination (integer)
 * @param tag the message tag (integer)
 * @param comm a Communicator (handle)
 * @param request a communication request (handle)
 *
 * @return an error value, it should be MPI_SUCCESS when everything went OK
 */
int pnl_object_mpi_isend (const PnlObject *Obj, int dest, int tag, MPI_Comm comm, MPI_Request *request)
{
  int info, size=0, pos=0;
  char *buf = NULL;
  info = pnl_object_mpi_pack_size (Obj, comm, &size);
  PNL_MPI_MESSAGE(info, "error in computing size.\n");
  if ((buf = malloc(size)) == NULL) return MPI_ERR_BUFFER;
  info = pnl_object_mpi_pack(Obj, buf, size, &pos, comm);
  PNL_MPI_MESSAGE( info, "error in packing.\n");
  info = MPI_Isend (buf, pos, MPI_PACKED, dest, tag, comm, request);
  free (buf); buf=NULL;
  return info;
}

/**
 * Performs a non-blocking receive of a PnlObject. Unlike all the other
 * functions, this one does not unpack the object because it would require to
 * issue a call to MPI_Wait which would inhibit the effect of using Irecv
 * instead of Recv
 *
 * @param buf (output) the buffer containing the packed object
 * @param size (output) the size of the packed object (needed to later unpack it) 
 * @param src the rank of source (integer)
 * @param tag the message tag (integer)
 * @param comm a Communicator (handle)
 * @param flag (output) set to FALSE if no matching message is ready for being
 * received, set to TRUE otherwise.
 * @param request (output) a communication request (handle)
 *
 * @return an error value, it should be MPI_SUCCESS when everything went OK
 */
int pnl_object_mpi_irecv (void **buf, int *size, int src, int tag, MPI_Comm comm, int *flag, MPI_Request *request)
{
  MPI_Status status;
  int info;
  info = MPI_Iprobe (src, tag, comm, flag, &status);
  PNL_MPI_MESSAGE( info, "error in IProbe.\n");
  if ( *flag == FALSE )
    {
      PNL_MESSAGE (*flag==FALSE, "Cannot find a matching message.\n");
      return MPI_SUCCESS;
    }
  info = MPI_Get_count (&status, MPI_PACKED, size);
  PNL_MPI_MESSAGE( info, "error in Get_count.\n");
  if ((*buf = malloc (*size)) == NULL) return MPI_ERR_BUFFER;
  info = MPI_Irecv (*buf, *size, MPI_PACKED, src, tag, comm, request);
  PNL_MPI_MESSAGE( info, "error in MPI_IRecv.\n");
  return info;
}

