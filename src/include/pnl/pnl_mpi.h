#ifndef __PNL_MPI_
#define __PNL_MPI_

#include <mpi.h>
#include "pnl_object.h"

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
extern int pnl_object_mpi_recv (PnlObject *Obj, int src, int tag, MPI_Comm comm, MPI_Status *status);
extern int pnl_object_mpi_bcast (PnlObject *Obj, int root, MPI_Comm comm);
extern int pnl_object_mpi_isend (const PnlObject *Obj, int dest, int tag, MPI_Comm comm, MPI_Request *request);
extern int pnl_object_mpi_irecv (void **buf, int *size, int src, int tag, MPI_Comm comm, int *flag, MPI_Request *request);

/* @} */

#endif /* __PNL_MPI_ */
