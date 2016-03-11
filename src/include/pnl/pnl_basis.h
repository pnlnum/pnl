#ifndef _PNL_BASIS_H
#define _PNL_BASIS_H

#include "pnl/pnl_matrix.h"
#include "pnl/pnl_sp_matrix.h"
#include "pnl/pnl_object.h"
#include "pnl/pnl_mathtools.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/**
 * \defgroup reg_basis Function Bases
 * It stores multivariate functions.
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

/**
 * @struct _PnlBasis
 *
 * It was originally designed to exclusively hold tensored functions with the exponents of
 * the tensor product stored in PnlBasis#T. This way, we could easily build multivariate polynomials as
 * tensor products of single variate ones. Because we are using polynomials, we advise to map the
 * original domain to [-1,1]^d to avoid numerical blow up; this is the purpose of
 * PnlBasis#isreduced, PnlBasis#center, PnlBasis#scale. These tensored functions can easily be
 * differentiated twice, see pnl_basis_eval_derivs_vect(), pnl_basis_eval_D_vect() and
 * pnl_basis_eval_D2_vect().
 *
 * To make this toolbox more complete, it is now possible to add some extra functions, which are not
 * tensor functions. They are stored using an independent mechanism in PnlBasis#func_list and typed as
 * #PnlRnFuncR.
 * These additional functions are only taken into account by the methods pnl_basis_i(),
 * pnl_basis_i_vect(), pnl_basis_eval() and pnl_basis_eval_vect(). Note in particular that it is not
 * possible to differentiate these functions.
 */
struct _PnlBasis
{
  /**
   * Must be the first element in order for the object mechanism to work
   * properly. This allows any PnlBasis pointer to be cast to a PnlObject
   */
  PnlObject     object;
  /** The basis type */
  int           id;
  /** The string to label the basis */
  const char   *label;
  /** The number of variates */
  int           nb_variates;
  /** The total number of elements in the basis */
  int           nb_func;
  /** The tensor matrix */
  PnlMatInt    *T;
  /** The sparse Tensor matrix */
  PnlSpMatInt  *SpT;
  /** The number of functions in the tensor #T */
  int           len_T;
  /** Compute the i-th element of the one dimensional basis.  As a convention, (*f)(x, 0) MUST be equal to 1 */
  double      (*f)(double    x, int i);
  /** Compute the first derivative of i-th element of the one dimensional basis */
  double      (*Df)(double   x, int i);
  /** Compute the second derivative of the i-th element of the one dimensional basis */
  double      (*D2f)(double  x, int i);
  /** TRUE if the basis is reduced */
  int           isreduced;
  /** The center of the domain */
  double       *center;
  /** The inverse of the scaling factor to map the domain to [-1, 1]^nb_variates */
  double       *scale;
  /** An array of additional functions */
  PnlRnFuncR  *func_list;
  /** The number of functions in #func_list */
  int          len_func_list;
};

extern int pnl_basis_type_register(const char *name, double (*f)(double, int), double (*Df)(double, int), double (*D2f)(double, int));
extern PnlBasis* pnl_basis_new();
extern PnlBasis* pnl_basis_create(int index, int nb_func, int space_dim);
extern PnlBasis* pnl_basis_create_from_degree(int index, int degree, int space_dim);
extern PnlBasis* pnl_basis_create_from_prod_degree(int index, int degree, int nb_variates);
extern PnlBasis* pnl_basis_create_from_hyperbolic_degree(int index, double degree, double q, int n);
extern void pnl_basis_clone(PnlBasis *dest, const PnlBasis *src);
extern PnlBasis* pnl_basis_copy(const PnlBasis *B);
extern void  pnl_basis_set_from_tensor(PnlBasis *b, int index, const PnlMatInt *T);
extern PnlBasis* pnl_basis_create_from_tensor( int index, const PnlMatInt *T);
extern void pnl_basis_del_elt(PnlBasis *B, const PnlVectInt *d);
extern void pnl_basis_del_elt_i(PnlBasis *B, int i);
extern void pnl_basis_add_elt(PnlBasis *B, const PnlVectInt *d);
extern void pnl_basis_set_domain(PnlBasis *B, const PnlVect *xmin, const PnlVect *xmax);
extern void pnl_basis_set_reduced(PnlBasis *B, const PnlVect *center, const PnlVect *scale);
extern void pnl_basis_free(PnlBasis **basis);
extern void pnl_basis_print(const PnlBasis *B);
extern int pnl_basis_fit_ls(const PnlBasis *f, PnlVect *coef, const PnlMat *x, const PnlVect *y);
extern double pnl_basis_ik_vect(const PnlBasis *b, const PnlVect *x, int i, int k);
extern double pnl_basis_i_vect(const PnlBasis *b, const PnlVect *x, int i);
extern double pnl_basis_i_D_vect(const PnlBasis *b, const PnlVect *x, int i, int j);
extern double pnl_basis_i_D2_vect(const PnlBasis *b, const PnlVect *x, int i, int j1, int j2);
extern double pnl_basis_eval_vect(const PnlBasis *basis, const PnlVect *coef, const PnlVect *x);
extern double pnl_basis_eval_D_vect(const PnlBasis *basis, const PnlVect *coef, const PnlVect *x, int i);
extern double pnl_basis_eval_D2_vect(const PnlBasis *basis, const PnlVect *coef, const PnlVect *x, int i, int j);
extern void pnl_basis_eval_derivs_vect(const PnlBasis *b, const PnlVect *coef, const PnlVect *x, double *val, PnlVect *grad, PnlMat *hes);
extern double pnl_basis_ik(const PnlBasis *b, const double *x, int i, int k);
extern double pnl_basis_i_D(const PnlBasis *b, const double *x, int i, int j);
extern double pnl_basis_i_D2(const PnlBasis *b, const double *x, int i, int j1, int j2);
extern double pnl_basis_eval(const PnlBasis *basis, const PnlVect *coef, const double *x);
extern double pnl_basis_eval_D(const PnlBasis *basis, const PnlVect *coef, const double *x, int i);
extern double pnl_basis_eval_D2(const PnlBasis *basis, const PnlVect *coef, const double *x, int i, int j);
extern void pnl_basis_eval_derivs(const PnlBasis *basis, const PnlVect *coef, const double *x, double *val, PnlVect *grad, PnlMat *hes);
extern void pnl_basis_add_function(PnlBasis *b, PnlRnFuncR *f);

#include "pnl/pnl_config.h"
#ifdef PNL_HAVE_INLINE
/**
 * Evaluate the i-th element of the basis b at the point x
 *
 * @deprecated Use the function pnl_basis_i_vect(const PnlBasis *b, const PnlVect *x, int i)
 *
 * @param b a PnlBasis
 * @param x a C array containing the coordinates of the point at which to
 * evaluate the basis
 * @param i an integer describing the index of the element of the basis to
 * considier
 *
 * @return f_i(x) where f is the i-th basis function
 */
PNL_INLINE_FUNC double pnl_basis_i(const PnlBasis *b, const double *x, int i)
{
  int k;
  double aux = 1.;
  if (i > b->len_T - 1)
    {
      PnlVect view = pnl_vect_wrap_array(x, b->nb_variates);
      return PNL_EVAL_RNFUNCR(&(b->func_list[i-b->len_T]), &view);
    }
  if (b->isreduced == 1)
    {
      for (k = b->SpT->I[i] ; k < b->SpT->I[i + 1] ; k++)
        {
          const int j = b->SpT->J[k];
          const int Tij = b->SpT->array[k];
          aux *= (b->f)((x[j] - b->center[j]) * b->scale[j], Tij);
        }
    }
  else
    {
      for (k = b->SpT->I[i] ; k < b->SpT->I[i + 1] ; k++)
        {
          const int j = b->SpT->J[k];
          const int Tij = b->SpT->array[k];
          aux *= (b->f)(x[j], Tij);
        }
    }
  return aux;
}
#endif
PNL_INLINE_DECL double pnl_basis_i(const PnlBasis *b, const double *x, int i);

/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_BASIS_H */
