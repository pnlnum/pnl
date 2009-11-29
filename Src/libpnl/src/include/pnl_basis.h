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
#define CANONICAL 0
#define HERMITIAN 1
#define TCHEBYCHEV 2

typedef struct {
  char * label; /*!< string to label the basis */
  int    space_dim; /*!< dimension of the space in which we are working */
  int    max_dim; /*!< maximum number of elements in the basis */
  double (*f)(double *x, int ind); /*!< the evaluation function */
  double (*Df)(double *x, int ind, int i); /*!< the first derivative function */
  double (*D2f)(double *x, int ind, int i, int j); /*!< the second derivative function */
} PnlBasis;



extern enum_members PnlBases ;
extern PnlBasis* pnl_basis_init ( int index, int nb_func, int space_dim);
extern int pnl_fit_least_squares (PnlVect *coef, PnlMat *x, PnlVect *y,
                                  PnlBasis *f, int dim_base);
extern double pnl_basis_eval (PnlVect *coef, double *x, PnlBasis *basis);
extern double pnl_basis_eval_D (PnlVect *coef, double *x, PnlBasis *basis, int i);
extern double pnl_basis_eval_D2 (PnlVect *coef, double *x, PnlBasis *basis, int i, int j);

/*@}*/
#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif
