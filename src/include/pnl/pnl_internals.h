#ifndef _PNL_INTERNALS_H
#define _PNL_INTERNALS_H

#include "pnl/pnl_object.h"
#include "pnl/pnl_complex.h"
#include "pnl/pnl_basis.h"

/**
 * \defgroup PnlInternals Internal macros
 *
 * It contains a few macros provided for conveniency when developing some
 * code inside the Pnl
 */
/*@{*/

#define FREE(x) if (x != NULL) { free (x); x = NULL; }
#define MALLOC_DOUBLE(n) malloc ((n) * sizeof (double))
#define MALLOC_COMPLEX(n) malloc ((n) * sizeof (dcomplex))
#define MALLOC_INT(n) malloc ((n) * sizeof (int))
#define PNL_MESSAGE(cond, msg)                  \
  if ( cond && pnl_message_is_on() )            \
    {                                           \
      printf (msg);                             \
    }

/* Check comment character when reading matrices from file */
#define iscomment(c) (((c) == '#') || ((c) == '%'))

extern double pnl_d1mach (int i);
extern double pnl_dlamch (char *);

#define NULLINT -1

#define BASIS_HAS_TENSOR_REP(basis) \
  (basis->id != PNL_BASIS_LOCAL)

/*@}*/

#endif /* _PNL_INTERNALS_H */
