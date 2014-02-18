#ifndef __PNL_MPI_
#define __PNL_MPI_

#include <mpi.h>
#include "pnl/pnl_object.h"
#include <pnl/pnl_list.h>
#include "pnl/pnl_random.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**
 * \defgroup PnlMpi MPI bindings
 *
 * We provide MPI bindings for most of the MPI Communication functions to
 * natively handle PnlObjects by Packing and Unpacking them transparently.
 */
/* @{ */


extern int pnl_object_mpi_pack_size (const PnlObject *Obj, MPI_Comm comm, int *size);
extern int pnl_object_mpi_pack (const PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm);
extern int pnl_object_mpi_unpack (PnlObject *Obj, void *buf, int bufsize, int *pos, MPI_Comm comm);
extern int pnl_object_mpi_send (const PnlObject *Obj, int dest, int tag, MPI_Comm comm);
extern int pnl_object_mpi_ssend (const PnlObject *Obj, int dest, int tag, MPI_Comm comm);
extern int pnl_object_mpi_recv (PnlObject *Obj, int src, int tag, MPI_Comm comm, MPI_Status *status);
extern int pnl_object_mpi_bcast (PnlObject *Obj, int root, MPI_Comm comm);
extern int pnl_object_mpi_isend (const PnlObject *Obj, int dest, int tag, MPI_Comm comm, MPI_Request *request);
extern int pnl_object_mpi_irecv (void **buf, int *size, int src, int tag, MPI_Comm comm, int *flag, MPI_Request *request);
extern int pnl_object_mpi_reduce (PnlObject *Sendbuf, PnlObject *Recvbuf, MPI_Op op, int root, MPI_Comm comm);

extern int pnl_rng_state_mpi_pack_size (const PnlRng *rng, MPI_Comm comm, int *size);
extern int pnl_rng_state_mpi_pack (const PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);
extern int pnl_rng_state_mpi_unpack (PnlRng *rng, void *buf, int bufsize, int *pos, MPI_Comm comm);

/*
 * Save load of objects
 */
extern PnlRng ** pnl_rng_create_from_file (char *str, int n);
extern int pnl_rng_save_to_file (PnlRng **rngtab, int n, char *str);
extern int pnl_object_save (PnlObject *obj, FILE *stream);
extern PnlObject* pnl_object_load (FILE *stream);
extern PnlList* pnl_object_load_into_list (FILE *stream);

/* @} */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PNL_MPI_ */
