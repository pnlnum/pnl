#ifndef _PNL_PERM_H
#define _PNL_PERM_H


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl/pnl_vector.h"

/**
 * \defgroup PnlPermutation  Permutation 
 */
/*@{*/
typedef PnlVectInt PnlPermutation;

extern PnlPermutation* pnl_permutation_new ();
extern PnlPermutation* pnl_permutation_create (int n);
extern void pnl_permutation_free (PnlPermutation **p);
extern void pnl_vect_permute (PnlVect *px, const PnlVect *x, const PnlPermutation *p);
extern void pnl_vect_permute_inplace (PnlVect *x, const PnlPermutation *p);
extern void pnl_permutation_fprint (FILE *fic, const PnlPermutation *p);
extern void pnl_permutation_print (const PnlPermutation *p);
/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_PERM_H */
