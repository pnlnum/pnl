#ifndef _PNL_PERM_H
#define _PNL_PERM_H

#include "pnl/pnl_vector.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**
 * \defgroup PnlPermutation  Permutation 
 */
/*@{*/
typedef PnlVectInt PnlPermutation;

extern PnlPermutation* pnl_permutation_new ();
extern PnlPermutation* pnl_permutation_create (int n);
extern void pnl_permutation_free (PnlPermutation **p);
extern void pnl_permutation_inverse (PnlPermutation *inv, const PnlPermutation *p);
extern void pnl_vect_permute (PnlVect *px, const PnlVect *x, const PnlPermutation *p);
extern void pnl_vect_permute_inplace (PnlVect *x, const PnlPermutation *p);
extern void pnl_vect_permute_inverse (PnlVect *px, const PnlVect *x, const PnlPermutation *p);
extern void pnl_vect_permute_inverse_inplace (PnlVect *x, const PnlPermutation *p);
extern void pnl_mat_col_permute (PnlMat *pX, const PnlMat *X, const PnlPermutation *p);
extern void pnl_mat_row_permute (PnlMat *pX, const PnlMat *X, const PnlPermutation *p);
extern void pnl_permutation_fprint (FILE *fic, const PnlPermutation *p);
extern void pnl_permutation_print (const PnlPermutation *p);
/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_PERM_H */
