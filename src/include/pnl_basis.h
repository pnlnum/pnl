#ifndef _BASIS_ND
#define _BASIS_ND

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


#include "pnl_types.h"
#include "pnl_matrix.h"
/**
 * \defgroup reg_basis Polynomial Bases
 */

/*@{*/

/* basis indices */
enum {CANONICAL, HERMITIAN, TCHEBYCHEV }; 

typedef struct PnlBasis_t PnlBasis;

struct PnlBasis_t {
  int         id;
  const char *label; /*!< string to label the basis */
  int         space_dim;  /*!< dimension of the space in which we are working */
  int         nb_func; /*!< number of elements in the basis */
  PnlMatInt  *T; /*!< Tensor matrix */
  double    (*f)(double    *x, int ind); /*!< the evaluation function */
  double    (*Df)(double   *x, int ind); /*!< the first derivative function */
  double    (*D2f)(double  *x, int ind); /*!< the second derivative function */
};



extern enum_members PnlBases ;
extern PnlBasis* pnl_basis_init ( int index, int nb_func, int space_dim);
extern void pnl_basis_free (PnlBasis **basis);
extern int pnl_fit_least_squares (PnlVect *coef, PnlMat *x, PnlVect *y,
                                  PnlBasis *f);
extern double pnl_basis_eval (PnlVect *coef, double *x, PnlBasis *basis);
extern double pnl_basis_eval_D (PnlVect *coef, double *x, PnlBasis *basis, int i);
extern double pnl_basis_eval_D2 (PnlVect *coef, double *x, PnlBasis *basis, int i, int j);

/*@}*/
#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif
