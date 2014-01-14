#ifndef _PNL_BASIS_H
#define _PNL_BASIS_H

#include "pnl/pnl_matrix.h"
#include "pnl/pnl_sp_matrix.h"
#include "pnl/pnl_object.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/**
 * \defgroup reg_basis Polynomial Bases
 */

/*@{*/

/* basis indices must start from 0 because they serve an index for the
 * PnlBasisTypeTab array */
enum {PNL_BASIS_NULL=-1, PNL_BASIS_CANONICAL=0, PNL_BASIS_HERMITE=1, PNL_BASIS_TCHEBYCHEV=2 };
/* synonymous for compatibility purposes */
#define CANONICAL PNL_BASIS_CANONICAL
#define HERMITIAN PNL_BASIS_HERMITE
#define TCHEBYCHEV PNL_BASIS_TCHEBYCHEV
#define PNL_BASIS_HERMITIAN PNL_BASIS_HERMITE

typedef struct _PnlBasis PnlBasis;

struct _PnlBasis
{
  /**
   * Must be the first element in order for the object mechanism to work
   * properly. This allows any PnlBasis pointer to be cast to a PnlObject
   */
  PnlObject     object;
  int           id;                     /*<! basis type */
  const char   *label;                  /*!< string to label the basis */
  int           nb_variates;            /*!< number of variates */
  int           nb_func;                /*!< number of elements in the basis */
  PnlMatInt    *T;                      /*!< Tensor matrix */
  PnlSpMatInt  *SpT;                    /*!< Sparse Tensor matrix */
  double      (*f)(double    x, int i); /*!< Computes the i-th element of the one dimensional basis.
                                          As a convention, (*f)(x, 0) MUST be equal to 1 */
  double      (*Df)(double   x, int i); /*!< Computes the first derivative of i-th element of
                                          the one dimensional basis */
  double      (*D2f)(double  x, int i); /*!< Computes the second derivative of the i-th element
                                          of the one dimensional basis */
  int           isreduced;              /* TRUE if the basis is reduced */
  double       *center;                 /*!< center of the domain */
  double       *scale;                  /*<! inverse of the scaling factor to map the
                                          domain to [-1, 1]^nb_variates */
};

extern int pnl_basis_type_register (const char *name, double (*f)(double, int), 
                             double (*Df)(double, int), double (*D2f)(double, int));
extern PnlBasis* pnl_basis_new ();
extern PnlBasis* pnl_basis_create ( int index, int nb_func, int space_dim);
extern PnlBasis* pnl_basis_create_from_degree ( int index, int degree, int space_dim);
extern PnlBasis* pnl_basis_create_from_prod_degree (int index, int degree, int nb_variates);
extern PnlBasis* pnl_basis_create_from_hyperbolic_degree (int index, double degree, double q, int n);
extern void pnl_basis_clone (PnlBasis *dest, const PnlBasis *src);
extern PnlBasis* pnl_basis_copy (const PnlBasis *B);
extern void  pnl_basis_set_from_tensor (PnlBasis *b, int index, const PnlMatInt *T);
extern PnlBasis* pnl_basis_create_from_tensor ( int index, const PnlMatInt *T);
extern void pnl_basis_del_elt (PnlBasis *B, const PnlVectInt *d);
extern void pnl_basis_del_elt_i (PnlBasis *B, int i);
extern void pnl_basis_add_elt (PnlBasis *B, const PnlVectInt *d);
extern void pnl_basis_set_domain (PnlBasis *B, const PnlVect *xmin, const PnlVect *xmax);
extern void pnl_basis_set_reduced (PnlBasis *B, const PnlVect *center, const PnlVect *scale);
extern void pnl_basis_free (PnlBasis **basis);
extern void pnl_basis_print (const PnlBasis *B);
extern int pnl_basis_fit_ls (const PnlBasis *f, PnlVect *coef, const PnlMat *x, const PnlVect *y);
extern double pnl_basis_ik_vect (const PnlBasis *b, const PnlVect *x, int i, int k);
extern double pnl_basis_i_vect (const PnlBasis *b, const PnlVect *x, int i);
extern double pnl_basis_i_D_vect (const PnlBasis *b, const PnlVect *x, int i, int j);
extern double pnl_basis_i_D2_vect (const PnlBasis *b, const PnlVect *x, int i, int j1, int j2);
extern double pnl_basis_eval_vect (const PnlBasis *basis, const PnlVect *coef, const PnlVect *x);
extern double pnl_basis_eval_D_vect (const PnlBasis *basis, const PnlVect *coef, const PnlVect *x, int i);
extern double pnl_basis_eval_D2_vect (const PnlBasis *basis, const PnlVect *coef, const PnlVect *x, int i, int j);
extern void pnl_basis_eval_derivs_vect (const PnlBasis *b, const PnlVect *coef, const PnlVect *x,
                                        double *val, PnlVect *grad, PnlMat *hes);
extern double pnl_basis_ik (const PnlBasis *b, const double *x, int i, int k);
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
