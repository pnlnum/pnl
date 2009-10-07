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

typedef double(*PnlBasis)(double *x, int index) ;

typedef struct {
  char * label; /*!< string to label the basis */
  int    space_dim; /*!< dimension of the space in which we are working */
  int    max_dim; /*!< maximum number of elements in the basis */
  double (*Compute)(double *x, int ind); /*!< the evaluation funtion */
} reg_basis;



extern enum_members PnlBases ;
extern PnlBasis pnl_init_basis ( int index, int nb_func, int space_dim);
extern int pnl_fit_least_squares (PnlVect *coef, PnlMat *x, PnlVect *y,
                                  PnlBasis *f, int dim_base);
extern double pnl_basis_eval (PnlVect *coef, double *x, PnlBasis *f);

/*@}*/
#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif
