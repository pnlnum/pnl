#ifndef _PNL_BASIS_H
#define _PNL_BASIS_H

#include "pnl/pnl_matrix.h"
#include "pnl/pnl_object.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/**
 * \defgroup reg_basis Polynomial Bases
 */

/*@{*/

/* basis indices */
enum {PNL_BASIS_CANONICAL, PNL_BASIS_HERMITIAN, PNL_BASIS_TCHEBYCHEV };
/* synonymous for compatibility purposes */
enum {CANONICAL, HERMITIAN, TCHEBYCHEV };

typedef struct _PnlBasis PnlBasis;

struct _PnlBasis
{
  /**
   * Must be the first element in order for the object mechanism to work
   * properly. This allows any PnlBasis pointer to be cast to a PnlObject
   */
  PnlObject object;
  int         id; /*<! basis type */
  const char *label; /*!< string to label the basis */
  int         nb_variates;  /*!< number of variates */
  int         nb_func; /*!< number of elements in the basis */
  PnlMatInt  *T; /*!< Tensor matrix */
  double    (*f)(double x, int i); /*!< Computes the i-th element of the one
                                            dimensional basis */
  double    (*Df)(double x, int i); /*!< Computes the first derivative of i-th element
                                            of the one dimensional basis */
  double    (*D2f)(double x, int i); /*!< Computes the second derivative of the i-th 
                                           element of the one dimensional basis */
  int         isreduced; /* TRUE if the basis is reduced */
  double     *center; /*!< center of the domain */
  double     *scale; /*<! inverse of the scaling factor to map the domain 
                          to [-1, 1]^nb_variates */
};

extern PnlBasis* pnl_basis_new ();
extern PnlBasis* pnl_basis_create ( int index, int nb_func, int space_dim);
extern PnlBasis* pnl_basis_create_from_degree ( int index, int degree, int space_dim);
extern PnlBasis*  pnl_basis_create_from_hyperbolic_degree (int index, double degree, double q, int n);
extern void  pnl_basis_set_from_tensor (PnlBasis *b, int index, const PnlMatInt *T);
extern PnlBasis* pnl_basis_create_from_tensor ( int index, const PnlMatInt *T);
extern void pnl_basis_set_domain (PnlBasis *B, const PnlVect *xmin, const PnlVect *xmax);
extern void pnl_basis_set_reduced (PnlBasis *B, const PnlVect *center, const PnlVect *scale);
extern void pnl_basis_free (PnlBasis **basis);
extern void pnl_basis_print (const PnlBasis *B);
extern int pnl_basis_fit_ls (const PnlBasis *f, PnlVect *coef, const PnlMat *x, const PnlVect *y);
extern double pnl_basis_i ( const PnlBasis *b, const double *x, int i );
extern double pnl_basis_i_D ( const PnlBasis *b, const double *x, int i, int j );
extern double pnl_basis_i_D2 (const PnlBasis *b, const double *x, int i, int j1, int j2);
extern double pnl_basis_eval (const PnlBasis *basis, const PnlVect *coef, const double *x);
extern double pnl_basis_eval_D (const PnlBasis *basis, const PnlVect *coef, const double *x, int i);
extern double pnl_basis_eval_D2 (const PnlBasis *basis, const PnlVect *coef, const double *x, int i, int j);
extern void pnl_basis_eval_derivs (const PnlBasis *basis, const PnlVect *coef, const double *x,
                                   double *val, PnlVect *grad, PnlMat *hes);


/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_BASIS_H */
