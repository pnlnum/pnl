
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
static int size_tausworthe_state (const PnlRng *rng, MPI_Comm comm, int *size);
static int size_lecuyer_state (const PnlRng *rng, MPI_Comm comm, int *size);
static int size_sqrt_state (const PnlRng *rng, MPI_Comm comm, int *size);
static int size_halton_state (const PnlRng *rng, MPI_Comm comm, int *size);
static int size_faure_state (const PnlRng *rng, MPI_Comm comm, int *size);
static int size_nied_state (const PnlRng *rng, MPI_Comm comm, int *size);
static int size_sobol_i4_state (const PnlRng *rng, MPI_Comm comm, int *size);
static int size_sobol_i8_state (const PnlRng *rng, MPI_Comm comm, int *size);

static int pack_knuth_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_mt_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_dcmt_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_mrgk3_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_mrgk5_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_shufl_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_lecuyer_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_tausworthe_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_sqrt_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_halton_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_faure_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_nied_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_sobol_i4_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int pack_sobol_i8_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);

static int unpack_knuth_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_mt_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_dcmt_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_mrgk3_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_mrgk5_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_shufl_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_lecuyer_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_tausworthe_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_sqrt_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_halton_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_faure_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_nied_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_sobol_i4_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
static int unpack_sobol_i8_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);

typedef int(pack_size_func)(const PnlRng *rng, MPI_Comm comm, int *size);
typedef int(pack_func)(const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
typedef int(unpack_func)(PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);

typedef struct
{
  PnlRngType t;
  pack_size_func *pack_size;
  pack_func      *pack;
  unpack_func    *unpack;
} PnlRngMPIFunc;

#define MAKE_PROPERTY(TYPE, str) {PNL_RNG_##TYPE, size_##str##_state, pack_##str##_state, unpack_##str##_state}

PnlRngMPIFunc rng_pack_func[] =
  {
    MAKE_PROPERTY(KNUTH,knuth),
    MAKE_PROPERTY(MERSENNE,mt),
    MAKE_PROPERTY(DCMT,dcmt),
    MAKE_PROPERTY(MRGK3,mrgk3),
    MAKE_PROPERTY(MRGK5,mrgk5),
    MAKE_PROPERTY(SHUFL,shufl),
    MAKE_PROPERTY(LECUYER,lecuyer),
    MAKE_PROPERTY(TAUSWORTHE,tausworthe),
    MAKE_PROPERTY(SQRT,sqrt),
    MAKE_PROPERTY(HALTON,halton),
    MAKE_PROPERTY(FAURE,faure),
    MAKE_PROPERTY(NIEDERREITER,nied),
    MAKE_PROPERTY(SOBOL_I4,sobol_i4),
    MAKE_PROPERTY(SOBOL_I8,sobol_i8),
    {PNL_RNG_NULL, NULL, NULL, NULL}
  };


static PnlRngMPIFunc* lookup (PnlRngType t)
{
  int i = 0;
  while ( rng_pack_func[i].t != PNL_RNG_NULL )
    {
      if (  rng_pack_func[i].t == t ) return &(rng_pack_func[i]);
      i++;
    }
  return NULL;
}


int pnl_rng_state_mpi_pack_size (const PnlRng *rng, MPI_Comm comm, int *size)
{
  PnlRngMPIFunc *property;
  if ((property = lookup (rng->type)) == NULL) return MPI_ERR_TYPE;
  return (*(property->pack_size)) (rng, comm, size);
}

int pnl_rng_state_mpi_pack (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  PnlRngMPIFunc *property;
  if ((property = lookup (rng->type)) == NULL) return MPI_ERR_TYPE;
  return (*(property->pack)) (rng, buf, bufsize, pos, comm);
}

int pnl_rng_state_mpi_unpack (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  PnlRngMPIFunc *property;
  if ((property = lookup (rng->type)) == NULL) return MPI_ERR_TYPE;
  return (*(property->unpack)) (rng, buf, bufsize, pos, comm);
}


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

static int size_lecuyer_state (const PnlRng *rng, MPI_Comm comm, int *size)
{
  int info, count;
  *size = 0;

  if((info=MPI_Pack_size(35,MPI_LONG, comm,&count))) return(info);
  *size += count;
  return info;
}

static int size_tausworthe_state (const PnlRng *rng, MPI_Comm comm, int *size)
{
  int info, count;
  *size = 0;

  if((info=MPI_Pack_size(4,MPI_UNSIGNED_LONG, comm,&count))) return(info);
  *size += count;
  return info;
}

static int size_mt_state (const PnlRng *rng, MPI_Comm comm, int *size)
{
  int info, count;
  *size = 0;

  if((info=MPI_Pack_size(PNL_MT_N,MPI_UNSIGNED_LONG, comm,&count))) return(info);
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
  if((info=MPI_Pack_size(PNL_DCMT_N,MPI_UNSIGNED_LONG, comm,&count))) return(info);
  *size += count;
  return info;
}

static int size_sqrt_state (const PnlRng *rng, MPI_Comm comm, int *size)
{
  int info, count;
  *size = 0;

  if((info=MPI_Pack_size(PNL_DIM_MAX_QMC,MPI_INT, comm,&count))) return(info);
  *size += count;
  if((info=MPI_Pack_size(PNL_DIM_MAX_QMC,MPI_DOUBLE, comm,&count))) return(info);
  *size += count;
  return info;
}

static int size_halton_state (const PnlRng *rng, MPI_Comm comm, int *size)
{
  int info, count;
  *size = 0;

  if((info=MPI_Pack_size(PNL_DIM_MAX_QMC,MPI_INT, comm,&count))) return(info);
  *size += count;
  return info;
}

static int size_faure_state (const PnlRng *rng, MPI_Comm comm, int *size)
{
  int info, count;
  *size = 0;

  if((info=MPI_Pack_size(1,MPI_INT, comm,&count))) return(info);
  *size += count;
  return info;
}

static int size_nied_state (const PnlRng *rng, MPI_Comm comm, int *size)
{
  int info, count;
  *size = 0;

  if((info=MPI_Pack_size(1,MPI_LONG, comm,&count))) return(info);
  *size += 2 * count;
  if((info=MPI_Pack_size(1,MPI_DOUBLE, comm,&count))) return(info);
  *size += count;
  if((info=MPI_Pack_size(1,MPI_UNSIGNED_LONG, comm,&count))) return(info);
  *size += count;
  if((info=MPI_Pack_size(PNL_DIM_MAX_NIED+1,MPI_UNSIGNED_LONG, comm,&count))) return(info);
  *size += count;
  return info;
}

static int size_sobol_i4_state (const PnlRng *rng, MPI_Comm comm, int *size)
{
  int info, count;
  *size = 0;

  if((info=MPI_Pack_size(1,MPI_INT, comm,&count))) return(info);
  *size += ( PNL_SOBOL_I4_LOG_MAX * PNL_SOBOL_I4_DIM_MAX2 + PNL_SOBOL_I4_DIM_MAX2 + 1) * count;
  if((info=MPI_Pack_size(1,MPI_FLOAT, comm,&count))) return(info);
  *size += count;
  return info;
}

static int size_sobol_i8_state (const PnlRng *rng, MPI_Comm comm, int *size)
{
  int info, count;
  *size = 0;

  if((info=MPI_Pack_size(1,MPI_LONG_LONG_INT, comm,&count))) return(info);
  *size += ( PNL_SOBOL_I8_LOG_MAX * PNL_SOBOL_I8_DIM_MAX2 + PNL_SOBOL_I8_DIM_MAX2 + 1) * count;
  if((info=MPI_Pack_size(1,MPI_DOUBLE, comm,&count))) return(info);
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

static int pack_lecuyer_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  lecuyer_state *s = (lecuyer_state *)(rng->state);
  if ((info=MPI_Pack(s,35,MPI_LONG,buf,bufsize,pos,comm))) return info;
  return info;
}

static int pack_tausworthe_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  tausworthe_state *s = (tausworthe_state *)(rng->state);
  if ((info=MPI_Pack(s,4,MPI_UNSIGNED_LONG,buf,bufsize,pos,comm))) return info;
  return info;
}

static int pack_mt_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  mt_state *s = (mt_state *)(rng->state);
  if ((info=MPI_Pack(s->mt,PNL_MT_N,MPI_UNSIGNED_LONG,buf,bufsize,pos,comm))) return info;
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
  if ((info=MPI_Pack(s->state,PNL_DCMT_N,MPI_UNSIGNED_LONG,buf,bufsize,pos,comm))) return info;
 return info;
}

static int pack_sqrt_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  sqrt_state *s = (sqrt_state *)(rng->state);
  if ((info=MPI_Pack(s->prime,PNL_DIM_MAX_QMC,MPI_INT,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(s->alpha,PNL_DIM_MAX_QMC,MPI_DOUBLE,buf,bufsize,pos,comm))) return info;
  return info;
}

static int pack_halton_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  halton_state *s = (halton_state *)(rng->state);
  if ((info=MPI_Pack(s->prime,PNL_DIM_MAX_QMC,MPI_INT,buf,bufsize,pos,comm))) return info;
  return info;
}

static int pack_faure_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  faure_state *s = (faure_state *)(rng->state);
  if ((info=MPI_Pack(&(s->r),1,MPI_INT,buf,bufsize,pos,comm))) return info;
  return info;
}

static int pack_nied_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  nied_state *s = (nied_state *)(rng->state);
  if ((info=MPI_Pack(&(s->saut),1,MPI_LONG,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->gray),1,MPI_LONG,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->facteur),1,MPI_DOUBLE,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->initial_d),1,MPI_UNSIGNED_LONG,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->initialX_n),PNL_DIM_MAX_NIED+1,MPI_UNSIGNED_LONG,buf,bufsize,pos,comm))) return info;
  return info;
}

static int pack_sobol_i4_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  sobol_i4_state *s = (sobol_i4_state *)(rng->state);
  if ((info=MPI_Pack(&(s->v),PNL_SOBOL_I4_DIM_MAX2*PNL_SOBOL_I4_LOG_MAX,MPI_INT,
                     buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->lastq),PNL_SOBOL_I4_DIM_MAX2,MPI_INT,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->recipd),1,MPI_FLOAT,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->maxcol),1,MPI_INT,buf,bufsize,pos,comm))) return info;
  return info;
}

static int pack_sobol_i8_state (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  sobol_i8_state *s = (sobol_i8_state *)(rng->state);
  if ((info=MPI_Pack(&(s->v),PNL_SOBOL_I8_DIM_MAX2*PNL_SOBOL_I8_LOG_MAX,MPI_LONG_LONG_INT,
                     buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->lastq),PNL_SOBOL_I8_DIM_MAX2,MPI_LONG_LONG_INT,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->recipd),1,MPI_DOUBLE,buf,bufsize,pos,comm))) return info;
  if ((info=MPI_Pack(&(s->maxcol),1,MPI_LONG_LONG_INT,buf,bufsize,pos,comm))) return info;
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

static int unpack_lecuyer_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  lecuyer_state *s = (lecuyer_state *)(rng->state);
  if ((info=MPI_Unpack(buf,bufsize,pos,s,35,MPI_LONG,comm))) return info;
  return info;
}

static int unpack_tausworthe_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  tausworthe_state *s = (tausworthe_state *)(rng->state);
  if ((info=MPI_Unpack(buf,bufsize,pos,s,4,MPI_UNSIGNED_LONG,comm))) return info;
  return info;
}

static int unpack_mt_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  mt_state *s = (mt_state *)(rng->state);
  if ((info=MPI_Unpack(buf,bufsize,pos,s->mt,PNL_MT_N,MPI_UNSIGNED_LONG,comm))) return info;
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
  if ((info=MPI_Unpack(buf,bufsize,pos,s->state,PNL_DCMT_N,MPI_UNSIGNED_LONG,comm))) return info;
 return info;
}

static int unpack_sqrt_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  sqrt_state *s = (sqrt_state *)(rng->state);
  if ((info=MPI_Unpack(buf,bufsize,pos,s->prime,PNL_DIM_MAX_QMC,MPI_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,s->alpha,PNL_DIM_MAX_QMC,MPI_DOUBLE,comm))) return info;
  return info;
}

static int unpack_halton_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  halton_state *s = (halton_state *)(rng->state);
  if ((info=MPI_Unpack(buf,bufsize,pos,s->prime,PNL_DIM_MAX_QMC,MPI_INT,comm))) return info;
  return info;
}

static int unpack_faure_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  faure_state *s = (faure_state *)(rng->state);
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->r),1,MPI_INT,comm))) return info;
  return info;
}

static int unpack_nied_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  nied_state *s = (nied_state *)(rng->state);
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->saut),1,MPI_LONG,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->gray),1,MPI_LONG,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->facteur),1,MPI_DOUBLE,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->initial_d),1,MPI_UNSIGNED_LONG,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->initialX_n),PNL_DIM_MAX_NIED+1,MPI_UNSIGNED_LONG,comm))) return info;
  return info;
}

static int unpack_sobol_i4_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  sobol_i4_state *s = (sobol_i4_state *)(rng->state);
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->v),PNL_SOBOL_I4_DIM_MAX2*PNL_SOBOL_I4_LOG_MAX,
                       MPI_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->lastq),PNL_SOBOL_I4_DIM_MAX2,MPI_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->recipd),1,MPI_FLOAT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->maxcol),1,MPI_INT,comm))) return info;
  return info;
}

static int unpack_sobol_i8_state (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm)
{
  int info;
  sobol_i8_state *s = (sobol_i8_state *)(rng->state);
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->v),PNL_SOBOL_I8_DIM_MAX2*PNL_SOBOL_I8_LOG_MAX,
                       MPI_LONG_LONG_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->lastq),PNL_SOBOL_I8_DIM_MAX2,MPI_LONG_LONG_INT,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->recipd),1,MPI_DOUBLE,comm))) return info;
  if ((info=MPI_Unpack(buf,bufsize,pos,&(s->maxcol),1,MPI_LONG_LONG_INT,comm))) return info;
  return info;
}
