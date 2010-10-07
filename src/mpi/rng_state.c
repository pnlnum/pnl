
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

#include <math.h>
#include <stdlib.h>
#include <mpi.h>

#include <pnl/pnl_object.h>
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_random.h"


static int size_knuth_state (const PnlRng *rng, MPI_Comm comm, int *size);
static int size_mt_state (const PnlRng *rng, MPI_Comm comm, int *size);
static int size_dcmt_state (const PnlRng *rng, MPI_Comm comm, int *size);
static int size_mrgk3_state (const PnlRng *rng, MPI_Comm comm, int *size);
static int size_mrgk5_state (const PnlRng *rng, MPI_Comm comm, int *size);
static int size_shufl_state (const PnlRng *rng, MPI_Comm comm, int *size);

static int pack_knuth_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_mt_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_dcmt_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_mrgk3_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_mrgk5_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_shufl_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);

static int unpack_knuth_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_mt_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_dcmt_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_mrgk3_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_mrgk5_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_shufl_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);


static int size_knuth_state (const PnlRng *rng, MPI_Comm comm, int *size)
{
  int info, count;
  *size = 0;

  if((info=MPI_Pack_size(1,MPI_LONG, comm,&count))) return(info);
  *size += count;
  if((info=MPI_Pack_size(2,MPI_INT, comm,&count))) return(info);
  *size += count;
  if((info=MPI_Pack_size(56,MPI_LONG, comm,&count))) return(info);
  *size += count;
  return info;
}

static int size_mrgk3_state (const PnlRng *rng, MPI_Comm comm, int *size)
{
  int info, count;
  *size = 0;

  if((info=MPI_Pack_size(6,MPI_DOUBLE, comm,&count))) return(info);
  *size += count;
  return info;
}

static int size_mrgk5_state (const PnlRng *rng, MPI_Comm comm, int *size)
{
  int info, count;
  *size = 0;

  if((info=MPI_Pack_size(10,MPI_DOUBLE, comm,&count))) return(info);
  *size += count;
  return info;
}

static int size_shufl_state (const PnlRng *rng, MPI_Comm comm, int *size)
{
  int info, count;
  *size = 0;

  if((info=MPI_Pack_size(34,MPI_LONG, comm,&count))) return(info);
  *size += count;
  return info;
}


static int size_mt_state (const PnlRng *rng, MPI_Comm comm, int *size)
{
  int info, count;
  *size = 0;

  if((info=MPI_Pack_size(MT_N,MPI_UNSIGNED_LONG, comm,&count))) return(info);
  *size += count;
  if((info=MPI_Pack_size(1,MPI_INT, comm,&count))) return(info);
  *size += count;
  return info;
}

static int size_dcmt_state (const PnlRng *rng, MPI_Comm comm, int *size)
{
  int info, count;
  *size = 0;

  if((info=MPI_Pack_size(6,MPI_UNSIGNED_LONG, comm,&count))) return(info);
  *size += count;
  if((info=MPI_Pack_size(9,MPI_INT, comm,&count))) return(info);
  *size += count;
  if((info=MPI_Pack_size(DCMT_N,MPI_UNSIGNED_LONG, comm,&count))) return(info);
  *size += count;
  return info;
}

static int pack_knuth_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  knuth_state *s = (knuth_state *)(rng->state);
  if ((info=MPI_Pack(&(s->SEED),1,MPI_LONG,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->inc1),1,MPI_INT,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->inc2),1,MPI_INT,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(s->t_alea,56,MPI_LONG,buf,bufsize,pos,comm))) return info;
  return info;
}

static int pack_mrgk3_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  mrgk3_state *s = (mrgk3_state *)(rng->state);
  if ((info=MPI_Pack(s,6,MPI_DOUBLE,buf,bufsize,pos,comm))) return info;
  return info;
}

static int pack_mrgk5_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  mrgk5_state *s = (mrgk5_state *)(rng->state);
  if ((info=MPI_Pack(s,10,MPI_DOUBLE,buf,bufsize,pos,comm))) return info;
  return info;
}

static int pack_shufl_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  shufl_state *s = (shufl_state *)(rng->state);
  if ((info=MPI_Pack(s,34,MPI_LONG,buf,bufsize,pos,comm))) return info;
  return info;
}


static int pack_mt_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  mt_state *s = (mt_state *)(rng->state);
  if ((info=MPI_Pack(s->mt,MT_N,MPI_UNSIGNED_LONG,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->mti),1,MPI_INT,buf,bufsize,pos,comm))) return info;
  return info;
}

static int pack_dcmt_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  dcmt_state *s = (dcmt_state *)(rng->state);
  if ((info=MPI_Pack(&(s->aaa),1,MPI_UNSIGNED_LONG,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->mm),1,MPI_INT,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->nn),1,MPI_INT,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->rr),1,MPI_INT,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->ww),1,MPI_INT,buf,bufsize,pos,comm))) return info; 
  if ((info=MPI_Pack(&(s->wmask),1,MPI_UNSIGNED_LONG,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->umask),1,MPI_UNSIGNED_LONG,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->lmask),1,MPI_UNSIGNED_LONG,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->shift0),1,MPI_INT,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->shift1),1,MPI_INT,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->shiftB),1,MPI_INT,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->shiftC),1,MPI_INT,buf,bufsize,pos,comm))) return info; 
  if ((info=MPI_Pack(&(s->maskB),1,MPI_UNSIGNED_LONG,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->maskC),1,MPI_UNSIGNED_LONG,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->i),1,MPI_INT,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(s->state,DCMT_N,MPI_UNSIGNED_LONG,buf,bufsize,pos,comm))) return info;
 return info;
}

static int unpack_knuth_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  knuth_state *s = (knuth_state *)(rng->state);
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->SEED),1,MPI_LONG,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->inc1),1,MPI_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->inc2),1,MPI_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,s->t_alea,56,MPI_LONG,comm))) return info;
  return info;
}

static int unpack_mrgk3_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  mrgk3_state *s = (mrgk3_state *)(rng->state);
  if ((info=MPI_Unpack(buf,bufsize,pos,s,6,MPI_DOUBLE,comm))) return info;
  return info;
}

static int unpack_mrgk5_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  mrgk5_state *s = (mrgk5_state *)(rng->state);
  if ((info=MPI_Unpack(buf,bufsize,pos,s,10,MPI_DOUBLE,comm))) return info;
  return info;
}

static int unpack_shufl_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  shufl_state *s = (shufl_state *)(rng->state);
  if ((info=MPI_Unpack(buf,bufsize,pos,s,34,MPI_LONG,comm))) return info;
  return info;
}

static int unpack_mt_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  mt_state *s = (mt_state *)(rng->state);
  if ((info=MPI_Unpack(buf,bufsize,pos,s->mt,MT_N,MPI_UNSIGNED_LONG,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->mti),1,MPI_INT,comm))) return info;
  return info;
}

static int unpack_dcmt_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  dcmt_state *s = (dcmt_state *)(rng->state);
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->aaa),1,MPI_UNSIGNED_LONG,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->mm),1,MPI_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->nn),1,MPI_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->rr),1,MPI_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->ww),1,MPI_INT,comm))) return info; 
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->wmask),1,MPI_UNSIGNED_LONG,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->umask),1,MPI_UNSIGNED_LONG,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->lmask),1,MPI_UNSIGNED_LONG,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->shift0),1,MPI_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->shift1),1,MPI_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->shiftB),1,MPI_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->shiftC),1,MPI_INT,comm))) return info; 
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->maskB),1,MPI_UNSIGNED_LONG,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->maskC),1,MPI_UNSIGNED_LONG,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->i),1,MPI_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,s->state,DCMT_N,MPI_UNSIGNED_LONG,comm))) return info;
 return info;
}



int pnl_rng_state_mpi_pack_size (const PnlRng *rng, MPI_Comm comm, int *size)
{
  switch (rng->type)
    {
    case PNL_RNG_KNUTH:
      return size_knuth_state (rng, comm, size);
      break;
    case PNL_RNG_MRGK3:
      return size_mrgk3_state (rng, comm, size);
      break;
    case PNL_RNG_MRGK5:
      return size_mrgk5_state (rng, comm, size);
      break;
    case PNL_RNG_SHUFL:
      return size_shufl_state (rng, comm, size);
      break;
    case PNL_RNG_MERSENNE:
      return size_mt_state (rng, comm, size);
      break;
    case PNL_RNG_DCMT:
      return size_dcmt_state (rng, comm, size);
      break;
    }
  return MPI_ERR_TYPE;
}

int pnl_rng_state_mpi_pack (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  switch (rng->type)
    {
    case PNL_RNG_KNUTH:
      return pack_knuth_state (rng, buf, bufsize, pos, comm);
      break;
    case PNL_RNG_MRGK3:
      return pack_mrgk3_state (rng, buf, bufsize, pos, comm);
      break;
    case PNL_RNG_MRGK5:
      return pack_mrgk5_state (rng, buf, bufsize, pos, comm);
      break;
    case PNL_RNG_SHUFL:
      return pack_shufl_state (rng, buf, bufsize, pos, comm);
      break;
    case PNL_RNG_MERSENNE:
      return pack_mt_state (rng, buf, bufsize, pos, comm);
      break;
    case PNL_RNG_DCMT:
      return pack_dcmt_state (rng, buf, bufsize, pos, comm);
      break;
    }
  return MPI_ERR_TYPE;
}

int pnl_rng_state_mpi_unpack (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  switch (rng->type)
    {
    case PNL_RNG_KNUTH:
      return unpack_knuth_state (rng, buf, bufsize, pos, comm);
      break;
    case PNL_RNG_MRGK3:
      return unpack_mrgk3_state (rng, buf, bufsize, pos, comm);
      break;
    case PNL_RNG_MRGK5:
      return unpack_mrgk5_state (rng, buf, bufsize, pos, comm);
      break;
    case PNL_RNG_SHUFL:
      return unpack_shufl_state (rng, buf, bufsize, pos, comm);
      break;
    case PNL_RNG_MERSENNE:
      return unpack_mt_state (rng, buf, bufsize, pos, comm);
      break;
    case PNL_RNG_DCMT:
      return unpack_dcmt_state (rng, buf, bufsize, pos, comm);
      break;
    }
  return MPI_ERR_TYPE;
}
