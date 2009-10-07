#ifndef PNL_PERM_H
#define PNL_PERM_H


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl_vector.h"

/**
 * \defgroup PnlPermutation  Permutation structure for Premia
 */
/*@{*/
typedef  struct {
  int size;
  int *array;
} PnlPermutation;

extern PnlPermutation* pnl_permutation_create (int n);
extern PnlPermutation* pnl_permutation_create_from_ptr (int n, const int *);
extern void pnl_permutation_init (PnlPermutation *p);
extern void pnl_permutation_swap (PnlPermutation *p, int i, int j);
extern void pnl_permutation_free (PnlPermutation **p);
extern void pnl_vect_permute (PnlVect *px, const PnlVect *x, const PnlPermutation *p);
extern void pnl_vect_permute_inplace (PnlVect *x, const PnlPermutation *p);
extern void pnl_permutation_fprint (FILE *fic, const PnlPermutation *p);
extern void pnl_permutation_print (const PnlPermutation *p);
/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* PNL_PERM_H */
