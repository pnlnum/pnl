/************************************************************************/
/* Copyright Jérôme Lelong <jerome.lelong@gmail.com>                    */
/*                                                                      */
/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as       */
/* published by the Free Software Foundation, either version 3 of the   */
/* License, or (at your option) any later version.                      */
/*                                                                      */
/* This program is distributed in the hope that it will be useful, but  */
/* WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    */
/* Lesser General Public License for more details.                      */
/*                                                                      */
/* You should have received a copy of the GNU Lesser General Public     */
/* License along with this program.  If not, see                        */
/* <http://www.gnu.org/licenses/>.                                      */
/************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pnl/pnl_config.h"
#include "pnl/pnl_basis.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_sp_matrix.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_specfun.h"

#ifndef PNL_RANGE_CHECK_OFF
#define CHECK_NB_FUNC(coef, basis)                                  \
  if ( coef->size != basis->nb_func )                               \
    {                                                               \
      PNL_ERROR("Wrong number of coefficients", "pnl_basis_eval");  \
    }
#define CHECK_NB_VARIATES(x,basis)                           \
  if ( x->size != basis->nb_variate )                        \
    {                                                        \
      PNL_ERROR("Dimension mismatch for the evaluation point", "pnl_basis_eval"); \
    }
#define CHECK_BASIS_TYPE(type, basis)                       \
  if ( type != basis->id )                                  \
    {                                                       \
      PNL_ERROR("Unexpected basis type", "pnl_basis_eval"); \
    }
#else
#define CHECK_NB_FUNC(coef, basis)
#define CHECK_NB_VARIATES(x, basis)
#define CHECK_BASIS_TYPE(type, basis)
#endif

#define CHECK_IS_DIFFERENTIABLE(basis)                    \
  if (basis->Df == NULL)                                  \
    {                                                     \
      fprintf(stderr, "Basis is not differentiable.");    \
      abort();                                            \
    }

#define CHECK_IS_CONSTRUCTIBLE_FROM_TENSOR(index, msg) \
  if (index == PNL_BASIS_LOCAL)                         \
    {                                                   \
      fprintf(stderr, msg);                             \
      return NULL;                                      \
    }

/**
 * Compute the number of elements of a polynomial with product degree
 *
 * The total degree function is the product of the partial degrees
 *
 * @param degree total degree
 * @param nb_variates number of variates of the polynomial
 */
static int nb_degrees_freedom_prod(int degree, int nb_variates)
{
  int partial_degree, nb_elements;
  if (degree == 0)
    {
      return 0;
    }
  if (nb_variates == 1)
    {
      return degree + 1;
    }
  nb_elements = 0;
  for (partial_degree = 0; partial_degree <= degree; partial_degree++)
    {
      nb_elements += nb_degrees_freedom_prod(degree / MAX(1, partial_degree), nb_variates - 1);
    }
  return nb_elements;
}

/**
 * Return the total degree (sum of the partial degrees) of the polynomial
 * represented T
 *
 * @param T an integer array of size @p n holding the partial degrees
 * @param n size of T
 * @return the total degree of the polynomials represented by T
 */
static double count_sum_degree(const int* T, int n, void* _params)
{
  int j, deg;
  deg = 0;

  for (j = 0 ; j < n ; j++)
    {
      deg += T[j];
    }
  return (double)deg;
}

/**
 * Return the total degree (sum of the partial degrees) of the polynomial
 * represented T
 *
 * The total degree is the product of MAX(1, d_i) where d_i is the partial
 * degree. If all d_i are zeros, the total degree is 0.
 *
 * @param T an integer array of size @p n holding the partial degrees
 * @param n size of T
 * @return the total degree of the polynomials represented by T
 */
static double count_prod_degree(const int* T, int n, void* _params)
{
  int j, deg, all_zero = PNL_TRUE;
  deg = 1;

  for (j = 0 ; j < n ; j++)
    {
      const int power = T[j];
      if (all_zero == PNL_TRUE && power > 0) all_zero = PNL_FALSE;
      deg *= MAX(power, 1);
    }
  if (all_zero == PNL_TRUE) return 0;
  return (double)deg;
}

/**
 * Return the total hyperbolic degree at the power q of the polynomial
 * represented by T
 *
 * @param T an integer array of size @p n holding the partial degrees
 * @param n size of T
 * @param params a pointer to a float containing the hyperbolic index
 * @return the total degree of the polynomials represented by T
 */
static double count_hyperbolic_degree(const int* T, int n, void *params)
{
  int j;
  double deg_q;
  double q = ((double *)params)[0];

  for (j = 0, deg_q = 0. ; j < n ; j++)
    {
      deg_q += pow(T[j], q);
    }
  return pow(deg_q, 1. / q);
}


/**
 * Compute the tensor matrix of the nb_variates variate basis with a total degree less or
 * equal than degree. The total degree is defined by the function count_degree
 *
 * @param degree the total degree
 * @param nb_variates the number of variates of the basis.
 * @param nb_elements the maximum number of elements in the tensor matrix
 * @param params extra parameters to pass to count_degree
 * @param count_degree a function to compute the total of a given line in a tensor
 *
 * @return the tensor matrix of the nb_variates variate basis with a total degree less or
 * equal than degree
 */
static PnlMatInt *compute_tensor_from_degree_function(int degree, int nb_variates, int nb_elements,
    void *params, double (*count_degree)(const int*, int, void*)
)
{
  int i, j;
  int *partial_degrees = calloc(nb_variates, sizeof(int));
  PnlMatInt *T = pnl_mat_int_create_from_zero(nb_elements, nb_variates);
  for (i = 1; i < nb_elements; i++)
    {
      for (j = nb_variates - 1; j >= 0; j--)
        {
          partial_degrees[j]++;
          if (count_degree(partial_degrees, nb_variates, params) <= degree) { break; }
          partial_degrees[j] = 0;
        }
      /* We could not find the next element */
      if (j == -1) { break; }
      pnl_mat_int_set_row_from_ptr(T, partial_degrees, i);
    }
  free(partial_degrees);
  pnl_mat_int_resize(T, i, nb_variates);
  return T;
}


static PnlMatInt *compute_tensor_from_sum_degree(int degree, int nb_variates)
{
  int nb_elements = (int) pnl_round(pnl_sf_choose(nb_variates + degree, degree));
  return compute_tensor_from_degree_function(degree, nb_variates, nb_elements, NULL, count_sum_degree);
}

/**
 * Compute the tensor matrix of the n-variate basis with a
 * hyperbolic degree of order q less or equal than degree
 *
 * @param q the hyperbolic index
 * @param nb_variates the number of variates of the basis.
 * @param degree the total hyperbolic degree
 *
 * @return the tensor matrix of the n-variate basis with a total degree less or equal than degree
 */
static PnlMatInt *compute_tensor_from_hyperbolic_degree(double degree, double q, int nb_variates)
{
  int nb_elements = (int) pnl_round(pnl_sf_choose(nb_variates + degree, degree));
  return compute_tensor_from_degree_function(degree, nb_variates, nb_elements, &q, count_hyperbolic_degree);
}

static PnlMatInt *compute_tensor_from_prod_degree(int degree, int nb_variates)
{
  int nb_elements = nb_degrees_freedom_prod(degree, nb_variates);
  return compute_tensor_from_degree_function(degree, nb_variates, nb_elements, NULL, count_prod_degree);
}

/**
 *  Canonical polynomials
 *  @param x the address of a real number
 *  @param l the index of the polynomial to be evaluated
 *  @param dim the index of the component on which the function is applied
 *  @param params extra parameters
 */
static double CanonicalD1(double x, int l, int dim, void *params)
{
  return pnl_pow_i(x, l);
}

/**
 *  First derivative of the Canonical polynomials
 *  @param x the address of a real number
 *  @param l the index of the polynomial whose first derivative is to be evaluated
 *  @param dim the index of the component on which the function is applied
 *  @param params extra parameters
 */
static double DCanonicalD1(double x, int l, int dim, void *params)
{
  if (l == 0) return 0.;
  return l * pnl_pow_i(x, l - 1);
}

/**
 *  Second derivative of the Canonical polynomials
 *  @param x the address of a real number
 *  @param l the index of the polynomial whose second derivative is to be evaluated
 *  @param dim the index of the component on which the function is applied
 *  @param params extra parameters
 */
static double D2CanonicalD1(double x, int l, int dim, void *params)
{
  if (l <= 1) return 0.;
  return l * (l - 1) * pnl_pow_i(x, l - 2);
}

/**
 * The terminal recursive function to compute Hermite polynomials of any
 * order. This function is only used for order > 7
 *  @param x the address of a real number
 *  @param n the order of the polynomial to be evaluated
 *  @param n0 rank of initialization
 *  @param f_n used to store the polynomial of order n0.
 *  @param f_n_1 used to store the polynomial of order n0 - 1.
 */
static double Hermite_rec(double x, int n, int n0, double *f_n, double *f_n_1)
{
  if (n == n0)
    {
      return *f_n;
    }
  else
    {
      double save = *f_n;
      *f_n = (x) * (*f_n) - n0 * (*f_n_1);
      *f_n_1 = save;
      return Hermite_rec(x, n, n0 + 1, f_n, f_n_1);
    }
}

/**
 *  Hermite polynomials
 *  @param x the address of a real number
 *  @param n the index of the polynomial to be evaluated
 *  @param dim the index of the component on which the function is applied
 *  @param params extra parameters
 */
static double HermiteD1(double x, int n, int dim, void *params)
{
  double val = x;
  double val2;
  double f_n, f_n_1;
  switch (n)
    {
    case 0 :
      return 1;
    case 1 :
      return val;
    case 2 :
      return val * val - 1.;
    case 3 :
      return (val * val - 3.) * val;
    case 4 :
      val2 = val * val;
      return (val2 - 6.) * val2 + 3;
    case 5 :
      val2 = val * val;
      return ((val2 - 10) * val2 + 15.) * val;
    case 6 :
      val2 = val * val;
      return ((val2 - 15.) * val2 + 45.) * val2 - 15.;
    case 7:
      val2 = val * val;
      return (((val2 - 21.) * val2 + 105.) * val2 - 105) * val;
    default:
      f_n = HermiteD1(x, 7, dim, params);
      f_n_1 = HermiteD1(x, 6, dim, params);
      return Hermite_rec(x, n, 7, &f_n, &f_n_1);
    }
}

/**
 *  First derivative of the Hermite polynomials
 *  @param x the address of a real number
 *  @param n the index of the polynomial whose derivative is to be evaluated
 *  @param params extra parameters
 */
static double DHermiteD1(double x, int n, int dim, void *params)
{
  if (n == 0) return 0.;
  else return n * HermiteD1(x, n - 1, dim, params);

}

/**
 *  Second derivative of the Hermite polynomials
 *  @param x the address of a real number
 *  @param n the index of the polynomial whose second derivative is to be evaluated
 *  @param dim the index of the component on which the function is applied
 *  @param params extra parameters
 */
static double D2HermiteD1(double x, int n, int dim, void *params)
{
  if (n == 0 || n == 1) return 0.;
  return n * (n - 1) * HermiteD1(x, n - 2, dim, params);
}

/**
 * The terminal recursive function to compute Tchebychev polynomials of any order.
 * @param x the address of a real number
 * @param n the order of the polynomial to be evaluated
 * @param n0 rank of initialization
 * @param f_n0 used to store the polynomial of order n0
 * @param f_n1 used to store the polynomial of order n0 - 1
 */
static double Tchebychev_rec(double x, int n, int n0, double *f_n0, double *f_n1)
{
  if (n == 7)
    {
      return *f_n0;
    }
  else
    {
      double save = *f_n0;
      *f_n0 = 2 * (x) * (*f_n0) - (*f_n1);
      *f_n1 = save;
      return Tchebychev_rec(x, n, n0 + 1, f_n0, f_n1);
    }
}

/**
 *  Tchebytchev polynomials of any order
 *  @param x the address of a real number
 *  @param n the order of the polynomial to be evaluated
 *  @param dim the index of the component on which the function is applied
 *  @param params extra parameters
 */
static double TchebychevD1(double x, int n, int dim, void *params)
{
  double val = x;
  double val2, val3, val4;
  double f_n, f_n_1;
  switch (n)
    {
    case 0 :
      return 1.;
    case 1 :
      return val;
    case 2 :
      return 2. * val * val - 1.;
    case 3 :
      return (4. * val * val - 3.) * val;
    case 4 :
      val2 = val * val;
      return 8. * val2 * val2 - 8. * val2 + 1.;
      break;
    case 5 :
      val2 = val * val;
      val3 = val2 * val;
      return 16. * val3 * val2 - 20. * val3 + 5.* val;
    case 6 :
      val2 = val * val;
      val4 = val2 * val2;
      return 32. * val4 * val2 - 48. * val4 + 18. * val2 - 1;
    case 7 :
      val2 = val * val;
      val3 = val2 * val;
      val4 = val2 * val2;
      return (64. * val4 - 112. * val2 + 56) * val3 - 7. * val;
    default :
      f_n = TchebychevD1(x, 7, dim, params);
      f_n_1 = TchebychevD1(x, 6, dim, params);
      return Tchebychev_rec(x, n, 7, &f_n, &f_n_1);
    }
}

/**
 * The terminal recursive function to compute the first derivative of the
 * Tchebychev polynomials of any order.
 *
 *  @param x the address of a real number
 *  @param n the order of the polynomial to be evaluated
 *  @param n0 rank of initialization
 *  @param f_n used to store the derivative of the polynomial of order n0.
 *  @param f_n_1 used to store the derivative of the polynomial of order
 *  n0-1.
 */
static double DTchebychev_rec(double x, int n, int n0, double *f_n, double *f_n_1)
{
  if (n == n0)
    {
      return *f_n;
    }
  else
    {
      double save = *f_n;
      *f_n = 2 * x * (double)(n0 + 1.0) / (double)n0 * (*f_n) - (double)(n0 + 1.) / (double)(n0 - 1.) * (*f_n_1);
      *f_n_1 = save;
      return DTchebychev_rec(x, n, n0 + 1, f_n, f_n_1);
    }
}

/**
 *  First derivative of the Tchebytchev polynomials
 *  @param x the address of a real number
 *  @param n the index of the polynomial whose first derivative is to be evaluated
 *  @param dim the index of the component on which the function is applied
 *  @param params extra parameters
 */
static double DTchebychevD1(double x, int n, int dim, void *params)
{
  double val = x;
  double val2, val4;
  double f_n, f_n_1;
  switch (n)
    {
    case 0 :
      return 0.;
    case 1 :
      return 1.;
    case 2 :
      return 4. * val;
    case 3 :
      return (12. * val * val - 3.);
    case 4 :
      return (32. * val * val - 16.) * val;
    case 5 :
      val2 = val * val;
      return 80. * val2 * val2 - 60. * val2 + 5.;
    case 6 :
      val2 = val * val;
      val4 = val2 * val2;
      return (192. * val4 - 192. * val2 + 36.) * val;
    case 7 :
      val2 = val * val;
      val4 = val2 * val2;
      return (448. * val4 - 560. * val2 + 168) * val2 - 7.;
    default :
      f_n = DTchebychevD1(x, 7, dim, params);
      f_n_1 = DTchebychevD1(x, 6, dim, params);
      return DTchebychev_rec(x, n, 7, &f_n, &f_n_1);
    }
}

/**
 * The terminal recursive function to compute the second derivative of the
 * Tchebychev polynomials of any order.
 *
 *  @param x the address of a real number
 *  @param n the order of the polynomial to be evaluated
 *  @param n0 rank of initialization
 *  @param f_n used to store the derivative of the polynomial of order n0
 *  @param f_n_1 used to store the derivative of  the polynomial of order n0 - 1
 *  @param dim the index of the component on which the function is applied
 *  @param params extra parameters
 */
static double D2Tchebychev_rec(double x, int n, int n0, double *f_n, double *f_n_1, int dim, void *params)
{
  if (n == n0)
    {
      return *f_n;
    }
  else
    {
      double save = *f_n;
      *f_n = 2 * x * (*f_n) - (*f_n_1) + 4 * DTchebychevD1(x, n0, dim, params);
      *f_n_1 = save;
      return D2Tchebychev_rec(x, n, n0 + 1, f_n, f_n_1, dim, params);
    }
}

/**
 *  Second derivative of the Tchebytchev polynomials
 *  @param x the address of a real number
 *  @param n the index of the polynomial whose second derivative is to be evaluated
 *  @param dim the index of the component on which the function is applied
 *  @param params extra parameters
 */
static double D2TchebychevD1(double x, int n, int dim, void *params)
{
  double val = x;
  double val2, val4;
  double f_n, f_n_1;
  switch (n)
    {
    case 0 :
      return 0.;
    case 1 :
      return 0.;
    case 2 :
      return 4.;
    case 3 :
      return 24. * val;
    case 4 :
      return (96. * val * val - 16.);
    case 5 :
      val2 = val * val;
      return 320. * val2 * val - 120. * val;
    case 6 :
      val2 = val * val;
      val4 = val2 * val2;
      return (960. * val4 - 576. * val2 + 36.);
    case 7 :
      val2 = val * val;
      val4 = val2 * val2;
      return (2688. * val4 - 2240. * val2 + 336) * val;
    default :
      f_n = D2TchebychevD1(x, 7, dim, params);
      f_n_1 = D2TchebychevD1(x, 6, dim, params);
      return D2Tchebychev_rec(x, n, 7, &f_n, &f_n_1, dim, params);
    }
}

static double LocalD1(double x, int n, int dim, void *params)
{
  /* For the tensor mechanism to work, the first interval must have index 1 and not 0. */
  n--;
  int *n_intervals = (int*) params;
  if (-1. + 2. * n / (double) n_intervals[dim] <= x && x < -1. + 2. * (n + 1) / (double) n_intervals[dim])
    {
      return 1.;
    }
  else
    {
      return 0.;
    }
}

static double DLocalD1(double x, int n, int dim, void *params)
{
  printf("Differentiating a local basis is not implemented.\n");
  abort();
}

static double D2LocalD1(double x, int n, int dim, void *params)
{
  printf("Differentiating a local basis is not implemented.\n");
  abort();
}

static double f_reduction_map(const PnlBasis *B, double x, int i, int dim)
{
  double map_x = x;
  if (B->map)
    {
      map_x = B->map(x, i, B->map_params);
    }
  if (B->isreduced)
    {
      return (B->f)((map_x - B->center[dim]) * B->scale[dim], i, dim, B->f_params);
    }
  return (B->f)(map_x, i, dim, B->f_params);
}

static double Df_reduction_map(const PnlBasis *B, double x, int i, int dim)
{
  double map_x = x;
  double y;
  if (B->map)
    {
      map_x = B->map(x, i, B->map_params);
    }
  if (B->isreduced)
    {
      y = (B->Df)((map_x - B->center[dim]) * B->scale[dim], i, dim, B->f_params) * B->scale[dim];
    }
  else
    {
      y = (B->Df)(map_x, i, dim, B->f_params);
    }
  if (B->Dmap)
    {
      y *= B->Dmap(x, dim, B->map_params);
    }
  return y;
}

static double D2f_reduction_map(const PnlBasis *B, double x, int i, int dim)
{
  double map_x = x;
  double y, y1, y2;
  if (B->map)
    {
      map_x = B->map(x, dim, B->map_params);
    }
  if (B->isreduced)
    {
      y2 = (B->D2f)((map_x - B->center[dim]) * B->scale[dim], i, dim, B->f_params) * B->scale[dim] * B->scale[dim];
    }
  else
    {
      y2 = (B->D2f)(map_x, i, dim, B->f_params);
    }
  if (B->Dmap && B->D2map)
    {
      double Dmap = B->Dmap(x, dim, B->map_params);
      if (B->isreduced)
        {
          y1 = (B->Df)((map_x - B->center[dim]) * B->scale[dim], i, dim, B->f_params) * B->scale[dim];
        }
      else
        {
          y1 = (B->Df)(map_x, i, dim, B->f_params);
        }
      y = y2 * Dmap * Dmap + y1 * (B->D2map)(x, dim, B->map_params);
    }
  else
    {
      y = y2;
    }
  return y;
}

/**
 * Struture used to describe the type of a basis
 */
typedef struct PnlBasisType_t PnlBasisType;
struct PnlBasisType_t
{
  int id;
  const char *label;
  double (*f)(double x, int l, int dim, void *params);
  double (*Df)(double x, int l, int dim, void *params);
  double (*D2f)(double x, int l, int dim, void *params);
  int is_orthogonal;
};

#define PNL_BASIS_MAX_TYPE 10
/**
 * The array holding the different basis types registered so far
 */
static PnlBasisType *PnlBasisTypeTab = NULL;
static int pnl_basis_type_next = 0; /*!< next available id for a basis type */
static int pnl_basis_type_tab_length = PNL_BASIS_MAX_TYPE; /*!< length of PnlBasisTypeTab */

/**
 * Register a new type of basis with a given index
 *
 * @param id index with which  the basis should be registered
 * @param label a string identifier
 * @param f the generating function in dimension 1
 * @param Df the first derivative of the generating function in dimension 1
 * @param D2f the second derivative of the generating function in dimension 1
 * @param is_orthogonal a boolean
 *
 * @return PNL_OK or PNL_FAIL
 */
static int pnl_basis_type_register_with_id(int id, const char *label,
    double (*f)(double x, int l, int dim, void *params),
    double (*Df)(double x, int l, int dim, void *params),
    double (*D2f)(double x, int l, int dim, void *params),
    int is_orthogonal
)
{
  /*
   * Enlarge the array if needed
   */
  if (id >= pnl_basis_type_tab_length)
    {
      pnl_basis_type_tab_length *= 2;
      PnlBasisTypeTab = realloc(PnlBasisTypeTab, pnl_basis_type_tab_length * sizeof(PnlBasisType));
    }
  if (pnl_basis_type_next != id) return PNL_FAIL;

  PnlBasisTypeTab[id].id = id;
  PnlBasisTypeTab[id].label = label;
  PnlBasisTypeTab[id].f = f;
  PnlBasisTypeTab[id].Df = Df;
  PnlBasisTypeTab[id].D2f = D2f;
  PnlBasisTypeTab[id].is_orthogonal = is_orthogonal;

  return PNL_OK;
}

/**
 * Initialize the array of basis types with
 *
 * @return
 */
static int pnl_basis_type_init()
{
  if (PnlBasisTypeTab != NULL)  return PNL_OK;
  PnlBasisTypeTab = malloc(PNL_BASIS_MAX_TYPE * sizeof(PnlBasisType));

  if (pnl_basis_type_register_with_id(PNL_BASIS_CANONICAL, "Canonical", CanonicalD1, DCanonicalD1, D2CanonicalD1, PNL_FALSE) != PNL_OK) return PNL_FAIL;
  pnl_basis_type_next++;
  if (pnl_basis_type_register_with_id(PNL_BASIS_HERMITE, "Hermite", HermiteD1, DHermiteD1, D2HermiteD1, PNL_FALSE) != PNL_OK) return PNL_FAIL;
  pnl_basis_type_next++;
  if (pnl_basis_type_register_with_id(PNL_BASIS_TCHEBYCHEV, "Tchebychev", TchebychevD1, DTchebychevD1, D2TchebychevD1, PNL_FALSE) != PNL_OK) return PNL_FAIL;
  pnl_basis_type_next++;
  if (pnl_basis_type_register_with_id(PNL_BASIS_LOCAL, "Local", LocalD1, DLocalD1, D2LocalD1, PNL_TRUE) != PNL_OK) return PNL_FAIL;
  pnl_basis_type_next++;

  return PNL_OK;
}

/**
 * Register a new type of basis
 *
 * @param name a string identifier
 * @param f the generating function in dimension 1
 * @param Df the first derivative of the generating function in dimension 1
 * @param D2f the second derivative of the generating function in dimension 1
 * @param is_orthogonal a boolean
 *
 * @return the next available index or PNL_BASIS_NULL if an error occurred
 */
int pnl_basis_type_register(const char *name,
  double (*f)(double x, int l, int dim, void *params),
  double (*Df)(double x, int l, int dim, void *params),
  double (*D2f)(double x, int l, int dim, void *params),
  int is_orthogonal
)
{
  int id;
  pnl_basis_type_init();
  id = pnl_basis_type_next;
  if (pnl_basis_type_register_with_id(id, name, f, Df, D2f, is_orthogonal) == PNL_FAIL)
    return PNL_BASIS_NULL;
  return id;
}

static char pnl_basis_label[] = "PnlBasis";
/**
 * Create an empty PnlBasis
 *
 * @return a PnlBasis
 */
PnlBasis *pnl_basis_new()
{
  PnlBasis *o;

  pnl_basis_type_init();
  if ((o = malloc(sizeof(PnlBasis))) == NULL) return NULL;
  o->nb_func = 0;
  o->id = 0;
  o->label = "";
  o->nb_variates = 0;
  o->SpT = NULL;
  o->isreduced = 0;
  o->center = NULL;
  o->scale = NULL;
  o->func_list = NULL;
  o->len_func_list = 0;
  o->len_T = 0;
  o->f_params = NULL;
  o->f_params_size = 0;
  o->map = o->Dmap = o->D2map = NULL;
  o->map_params = NULL;
  o->map_params_size = 0;
  o->object.type = PNL_TYPE_BASIS;
  o->object.parent_type = PNL_TYPE_BASIS;
  o->object.label = pnl_basis_label;
  o->object.nref = 0;
  o->object.destroy = (DestroyFunc *) pnl_basis_free;
  o->object.constructor = (NewFunc *) pnl_basis_new;
  o->object.clone = (CloneFunc *) pnl_basis_clone;
  o->object.copy = (CopyFunc *) pnl_basis_copy;
  return o;
}

/**
 * Set an existing basis to a given full tensor representation
 *
 * @param b an already allocated basis (as returned by pnl_basis_new for instance)
 * @param T the full tensor of the multi-dimensional basis.
 */
void pnl_basis_set_from_tensor(PnlBasis *b, const PnlMatInt *T)
{
  PnlSpMatInt *SpT = pnl_sp_mat_int_create_from_mat(T);
  pnl_basis_set_from_sparse_tensor(b, SpT);
}

/**
 * Set an existing basis to a given sparse tensor representation
 *
 * @param b an already allocated basis (as returned by pnl_basis_new for instance)
 * @param SpT the sparse tensor of the multi-dimensional basis. No copy of SpT is done, so
 * do not free it. It will be freed transparently by pnl_basis_free
 */
void pnl_basis_set_from_sparse_tensor(PnlBasis *b, const PnlSpMatInt *SpT)
{
  b->nb_func = SpT->m;
  b->nb_variates = SpT->n;
  if(b->func_list) free(b->func_list);
  b->len_func_list = 0;
  b->len_T = SpT->m;

  /* Not sure this is the right place to put it */
  pnl_sp_mat_int_free(&(b->SpT));
  b->SpT = (PnlSpMatInt *) SpT;
}

/**
 * Set the type of basis
 *
 * @param B a PnlBasis
 * @param index the index of the family to be used
 */
void pnl_basis_set_type(PnlBasis *B, int index)
{
  B->id = index;
  B->label = PnlBasisTypeTab[index].label;
  B->f = PnlBasisTypeTab[index].f;
  B->Df = PnlBasisTypeTab[index].Df;
  B->D2f = PnlBasisTypeTab[index].D2f;
}

/**
 * Return a  PnlBasis
 *
 * @param index the index of the family to be used
 * @param T the tensor of the multi-dimensional basis. No copy of T is done, so
 * do not free T. It will be freed transparently by pnl_basis_free
 * @return a PnlBasis
 */
PnlBasis *pnl_basis_create_from_tensor(int index, const PnlMatInt *T)
{
  PnlBasis *b;
  if ((b = pnl_basis_new()) == NULL) return NULL;
  pnl_basis_set_from_tensor(b, T);
  pnl_basis_set_type(b, index);
  return b;
}

/**
 * Return a  PnlBasis
 *
 * @param index the index of the family to be used
 * @param nb_func the maximum number of functions which may be used
 * @param nb_variates the size of the space in which the basis functions are
 * defined
 * @return a PnlBasis
 */
PnlBasis *pnl_basis_create(int index, int nb_func, int nb_variates)
{
  int degree = 0;
  PnlMatInt *T;
  CHECK_IS_CONSTRUCTIBLE_FROM_TENSOR(index, "Use pnl_basis_local_create to create a local basis");
  for (degree = 0; 1; degree++)
    {
      if (pnl_round(pnl_sf_choose(nb_variates + degree, degree)) > nb_func)
        {
          break;
        }
    }
  degree--;
  T = compute_tensor_from_degree_function(degree, nb_variates, nb_func, NULL, count_sum_degree);
  return pnl_basis_create_from_tensor(index, T);
}

/**
 * Return a  PnlBasis
 *
 * @param index the index of the family to be used
 * @param degree the maximum total degree of the elements in the basis
 * @param nb_variates the size of the space in which the basis functions are
 * defined
 * @return a PnlBasis
 */
PnlBasis *pnl_basis_create_from_degree(int index, int degree, int nb_variates)
{
  PnlMatInt *T;
  CHECK_IS_CONSTRUCTIBLE_FROM_TENSOR(index, "Use pnl_basis_local_create to create a local basis");
  T = compute_tensor_from_sum_degree(degree, nb_variates);
  return pnl_basis_create_from_tensor(index, T);
}

/**
 * Return a  PnlBasis
 *
 * @param index the index of the family to be used
 * @param degree the maximum total degree of the elements in the basis
 * @param nb_variates the size of the space in which the basis functions are
 * defined
 * @return a PnlBasis
 */
PnlBasis *pnl_basis_create_from_prod_degree(int index, int degree, int nb_variates)
{
  PnlMatInt *T;
  CHECK_IS_CONSTRUCTIBLE_FROM_TENSOR(index, "Use pnl_basis_local_create to create a local basis");
  T = compute_tensor_from_prod_degree(degree, nb_variates);
  return pnl_basis_create_from_tensor(index, T);
}

/**
 * Return a  PnlBasis built using an hyperbolic set of indices
 *
 * @param index the index of the family to be used
 * @param degree the hyperbolic maximum degree of the
 * elements in the basis
 * @param q the hyperbolic exponent (0 < q <= 1)
 * @param n the size of the space in which the basis functions are
 * defined
 * @return a PnlBasis such that every element prod_{i=1}^n f_i^(a_i) satisfies
 * (sum_{i=1}^n (a_i^q))^(1/q) <= degree
 */
PnlBasis *pnl_basis_create_from_hyperbolic_degree(int index, double degree, double q, int n)
{
  PnlMatInt *T;
  CHECK_IS_CONSTRUCTIBLE_FROM_TENSOR(index, "Use pnl_basis_local_create to create a local basis");
  T = compute_tensor_from_hyperbolic_degree(degree, q, n);
  return pnl_basis_create_from_tensor(index, T);
}

/**
 * Return a PnlBasis with local and orthogonal function
 *
 * @param n_intervals this is an array of size @p space_dim describing the number of intervals for every dimension.
 * @param space_dim the dimension of the state space
 * @return PnlBasis
 */
PnlBasis* pnl_basis_local_create(const int *n_intervals, int space_dim)
{
  PnlBasis *B;
  int i, prod = 1;
  if ((B = pnl_basis_new()) == NULL) return NULL;
  pnl_basis_set_type(B, PNL_BASIS_LOCAL);
  B->f_params_size = space_dim * sizeof(int);
  B->f_params = malloc(B->f_params_size);
  memcpy(B->f_params, n_intervals, B->f_params_size);
  B->nb_variates = space_dim;
  for (i = 0; i < space_dim; i++)
    {
      prod *= n_intervals[i];
    }
  B->nb_func = prod;
  return B;
}

/**
 * Return a PnlBasis with local and orthogonal function
 *
 * @param n_intervals number of intervals per dimension
 * @param space_dim the dimension of the state space
 * @return PnlBasis
 */
PnlBasis* pnl_basis_local_create_regular(int n_intervals, int space_dim)
{
  PnlBasis *B;
  int i;
  int *intervals_tab = malloc(space_dim * sizeof(int));
  for (i = 0; i < space_dim; i++) intervals_tab[i] = n_intervals;
  B = pnl_basis_local_create(intervals_tab, space_dim);
  free(intervals_tab);
  return B;
}

/**
 * Add an extra function to an existing basis.
 *
 * Note that it is not possible to make a deep copy of f->params, so make sure that its address
 * remains valid as long as b is used. If this does not sound clear, just remember it can be granted
 * by never making a function returning a PnlBasis, nor encapsulating it into a struct.
 *
 * @param b an existing basis
 * @param f a function to be added to b->func_list
 */
void pnl_basis_add_function(PnlBasis *b, PnlRnFuncR *f)
{
  ++(b->len_func_list);
  ++(b->nb_func);
  b->func_list = realloc(b->func_list, b->len_func_list * sizeof(PnlRnFuncR));
  (b->func_list[b->len_func_list - 1]).F = f->F;
  (b->func_list[b->len_func_list - 1]).params = f->params;
}

/**
 * Clone a basis
 *
 * @param dest destination
 * @param src source
 */
void pnl_basis_clone(PnlBasis *dest, const PnlBasis *src)
{
  pnl_basis_set_from_sparse_tensor(dest, src->SpT);
  dest->id = src->id;
  dest->label = src->label;
  dest->f = src->f;
  dest->Df = src->Df;
  dest->D2f = src->D2f;
  dest->f_params = realloc(dest->f_params, src->f_params_size);
  memcpy(dest->f_params, src->f_params, src->f_params_size);
  dest->f_params_size = src->f_params_size;
  pnl_basis_reset_reduced(dest);
  if (src->isreduced == 1)
    {
      int n = dest->nb_variates;
      dest->isreduced = 0;
      dest->center = malloc(n * sizeof(double));
      dest->scale = malloc(n * sizeof(double));
      memcpy(dest->center, src->center, n * sizeof(double));
      memcpy(dest->scale, src->scale, n * sizeof(double));
    }
  /* Copy the extra functions if any*/
  dest->func_list = malloc(src->len_func_list * sizeof(PnlRnFuncR));
  memcpy(dest->func_list, src->func_list, src->len_func_list * sizeof(PnlRnFuncR));
  dest->len_func_list = src->len_func_list;
}

/**
 * Delete the i-th basis function, ie. remove line i from the tensor
 *
 * @param B a PnlBasis
 * @param i an integer between 0 and B->nb_func-1
 */
void pnl_basis_del_elt_i(PnlBasis *B, int i)
{
  PNL_CHECK(i >= B->len_T, "size mismatch",  "basis_del_elt_i");
  pnl_sp_mat_int_del_row(B->SpT, i);
  --(B->nb_func);
  --(B->len_T);
}

/**
 * Delete a basis function described by its tensor decomposition
 *
 * @param B a PnlBasis
 * @param d a PnlVectInt describing the element to remove
 */
void pnl_basis_del_elt(PnlBasis *B, const PnlVectInt *d)
{
  int i;
  PNL_CHECK(B->nb_variates != d->size, "size mismatch",  "basis_del_elt");
  for (i = 0 ; i < B->nb_func ; i++)
    {
      int j;
      int test = 1;
      for (j = 0; j < d->size; j++)
        {
          if (pnl_sp_mat_int_get(B->SpT, i, j) != pnl_vect_int_get(d, j))
            {
              test = 0;
              break;
            }
        }
      if (test) break;
    }
  if (i < B->nb_func) pnl_basis_del_elt_i(B, i);
}

/**
 * Add a function described by its tensor decomposition to a basis
 *
 * @param B a PnlBasis
 * @param d a PnlVectInt describing the function to add
 */
void pnl_basis_add_elt(PnlBasis *B, const PnlVectInt *d)
{
  PNL_CHECK(B->nb_variates != d->size, "size mismatch",  "basis_del_elt");
  pnl_sp_mat_int_add_row(B->SpT, B->SpT->m, d);
  ++(B->nb_func);
  ++(B->len_T);
}

/**
 * Create a copy of a basis
 *
 * @param B a basis
 *
 * @return
 */
PnlBasis *pnl_basis_copy(const PnlBasis *B)
{
  PnlBasis *BB;
  BB = pnl_basis_new();
  pnl_basis_clone(BB, B);
  return BB;
}

/**
 * Set the center and scale field of a PnlBasis using the domain on which
 * the basis will be used
 *
 * @param B
 * @param xmin lower bounds of the domain
 * @param xmax upper bounds of the domain
 */
void pnl_basis_set_domain(PnlBasis *B, const PnlVect *xmin, const PnlVect *xmax)
{
  int i, n;
  PNL_CHECK(xmin->size != xmax->size || xmin->size != B->nb_variates,
            "size mismatch", "pnl_basis_set_domain");

  n = B->nb_variates;
  B->center = realloc(B->center, n * sizeof(double));
  B->scale = realloc(B->scale, n * sizeof(double));
  B->isreduced = 1;
  for (i = 0 ; i < n ; i++)
    {
      const double low = PNL_GET(xmin, i);
      const double high = PNL_GET(xmax, i);
      B->center[i] = (low + high) / 2.;
      B->scale[i] = 2. / (high - low);
    }
}

/**
 * Set the center and scale field of a PnlBasis
 *
 * @param B
 * @param center center of the domain
 * @param scale width of the domain in each direction
 */
void pnl_basis_set_reduced(PnlBasis *B, const PnlVect *center, const PnlVect *scale)
{
  int i, n;
  PNL_CHECK(center->size != scale->size || center->size != B->nb_variates,
            "size mismatch", "pnl_basis_set_reduced");

  n = B->nb_variates;
  B->center = realloc(B->center, n * sizeof(double));
  B->scale = realloc(B->scale, n * sizeof(double));
  B->isreduced = 1;
  memcpy(B->center, center->array, n * sizeof(double));
  for (i = 0 ; i < n ; i++)
    {
      B->scale[i] = 1. / PNL_GET(scale, i);
    }
}

/**
 * @brief Reset center and scale. @p is_reduced is et to 0.
 *
 * @param B a PnlBasis
 */
void pnl_basis_reset_reduced(PnlBasis *B)
{
  if (B->isreduced == 1)
    {
      free(B->center);
      B->center = NULL;
      free(B->scale);
      B->scale = NULL;
      B->isreduced = 0;
    }
}

/**
 * @brief Define the non linear map to apply to every input variable before evaluating the basis functions. First, the centering/rescaling is applied, second this mapping and third the basis functions are evaluated
 * 
 * @param map The mapping
 * @param Dmap The first derivative of the mapping
 * @param D2map The second derivative of the mapping
 * @param params The extra parameters to be passed to the mapping function
 * @param size_params The size in bytes of @p params
 */
void pnl_basis_set_map(
  PnlBasis *B,
  double (*map)(double, int, void*),
  double (*Dmap)(double, int, void*),
  double (*D2map)(double, int, void*),
  void *params, size_t size_params
)
{
  B->map = map;
  B->Dmap = Dmap;
  B->D2map = D2map;
  B->map_params = realloc(B->map_params, size_params);
  memcpy(B->map_params, params, size_params);
  B->map_params_size = size_params;
}

/**
 * Free a PnlBasis
 *
 * @param B
 */
void pnl_basis_free(PnlBasis **B)
{
  if (*B == NULL) return;
  pnl_sp_mat_int_free(&(*B)->SpT);
  if ((*B)->f_params_size > 0)
    {
      free((*B)->f_params);
      (*B)->f_params = NULL;
      (*B)->f_params_size = 0;
    }
  if ((*B)->map_params_size > 0)
    {
      free((*B)->map_params);
      (*B)->map_params = NULL;
      (*B)->map_params_size = 0;
    }
  if ((*B)->isreduced == 1)
    {
      free((*B)->center);
      (*B)->center = NULL;
      free((*B)->scale);
      (*B)->scale = NULL;
    }
  if ((*B)->len_func_list > 0) free((*B)->func_list);
  free(*B);
  *B = NULL;
}

/**
 * Print a PnlBasis
 * @param B a basis
 */
void pnl_basis_print(const PnlBasis *B)
{
  printf("Basis Name: %s\n", B->label);
  printf("\tNumber of variates: %d\n", B->nb_variates);
  printf("\tExtra parameters size: %zu\n", B->f_params_size);
  printf("\tExtra parameters size: %zu\n", B->f_params_size);
  printf("\tNumber of functions in tensor: %d\n", B->len_T);;
  printf("\tNumber of extra functions: %d\n", B->len_func_list);
  printf("\tTotal number of functions: %d\n", B->nb_func);
  printf("\tIs reduced = %d\n", B->isreduced);
  if (B->isreduced)
    {
      int i;
      printf("\tcenter = ");
      for (i = 0 ; i < B->nb_variates ; i++) printf("%f ", B->center[i]);
      printf("\n\tscale = ");
      for (i = 0 ; i < B->nb_variates ; i++) printf("%f ", B->scale[i]);
      printf("\n");
    }
  if (B->SpT)
    {
      printf("\tTensor matrix: \n");
      pnl_sp_mat_int_print(B->SpT);
      printf("\n");
    }
}

/**
 * An element of a basis writes as a product
 *      p_1(x_1) p2(x_2) .... p_n(x_n)
 * for a polynomial with n variates. Each p_k is a polynomial with only
 * one variate.
 *
 * @deprecated Use the function pnl_basis_ik_vect(const PnlBasis *b, const PnlVect *x, int i, int k)
 *
 * This functions evaluates the term p_k of the i-th element of the
 * basis b at the point x
 *
 * @param b a PnlBasis
 * @param x a C array containing the coordinates of the point at which to
 * evaluate the basis
 * @param i an integer describing the index of the element of the basis to
 * consider
 * @param k the index of the term to be evaluated within element i of the
 * basis
 *
 * @return (f_i)_k (x)
 */
double pnl_basis_ik(const PnlBasis *b, const double *x, int i, int k)
{
  const int Tik = pnl_sp_mat_int_get(b->SpT, i, k);
  if (Tik == 0) return 1.;
  return f_reduction_map(b, x[k], Tik, k);
}

/**
 * Evaluate the i-th element of the basis b at the point x
 *
 * @deprecated Use the function pnl_basis_i_vect(const PnlBasis *b, const PnlVect *x, int i)
 *
 * @param b a PnlBasis
 * @param x a C array containing the coordinates of the point at which to
 * evaluate the basis
 * @param i an integer describing the index of the element of the basis to consider
 *
 * @return f_i(x) where f is the i-th basis function
 */
double pnl_basis_i(const PnlBasis *b, const double *x, int i)
{
  int k;
  double aux = 1.;
  if (i > b->len_T - 1)
    {
      PnlVect view = pnl_vect_wrap_array(x, b->nb_variates);
      return PNL_EVAL_RNFUNCR(&(b->func_list[i-b->len_T]), &view);
    }
    for (k = b->SpT->I[i] ; k < b->SpT->I[i + 1] ; k++)
      {
        const int j = b->SpT->J[k];
        const int Tij = b->SpT->array[k];
        aux *= f_reduction_map(b, x[j], Tij, j);
      }
  return aux;
}

/**
 * First order derivative
 *
 * @deprecated Use the function pnl_basis_i_D_vect(const PnlBasis *b, const PnlVect *x, int i, int j)
 *
 * @param b a basis
 * @param x the point at which to evaluate the first derivative
 * @param i the index of the basis element to differentiate
 * @param j the index of the variable w.r.t which we differentiate
 *
 * @return (D(b_i)/Dj)(x)
 */
double pnl_basis_i_D(const PnlBasis *b, const double *x, int i, int j)
{
  int l;
  double aux = 1;

  CHECK_IS_DIFFERENTIABLE(b);
  /* Test if the partial degree is small enough to return 0 */
  if (pnl_sp_mat_int_get(b->SpT, i, j) == 0) return 0.;

  for (l = b->SpT->I[i] ; l < b->SpT->I[i + 1] ; l++)
    {
      const int k = b->SpT->J[l];
      const int Tik = b->SpT->array[l];
      if (k == j)
        aux *= Df_reduction_map(b, x[k], Tik, k);
      else
        aux *= f_reduction_map(b, x[k], Tik, k);
    }
  return aux;
}

/**
 * Second order derivative
 *
 * @deprecated Use the function pnl_basis_i_D2_vect(const PnlBasis *b, const PnlVect *x, int i, int j1, int j2)
 *
 * @param b a basis
 * @param x the point at which to evaluate the first derivative
 * @param i the index of the basis element to differentiate
 * @param j1 the index of the first variable w.r.t which we differentiate
 * @param j2 the index of the second variable w.r.t which we differentiate
 *
 * @return (D(b_i)/(Dj1 Dj2))(x)
 */
double pnl_basis_i_D2(const PnlBasis *b, const double *x, int i, int j1, int j2)
{
  int l;
  double aux = 1;
  CHECK_IS_DIFFERENTIABLE(b);
  /* Test if the partial degree is small enough to return 0 */
  if ((j1 == j2) && pnl_sp_mat_int_get(b->SpT, i, j1) <= 1) return 0.;
  if (pnl_sp_mat_int_get(b->SpT, i, j1) == 0 || pnl_sp_mat_int_get(b->SpT, i, j2) == 0) return 0.;

  if (j1 == j2)
    {
      for (l = b->SpT->I[i] ; l < b->SpT->I[i + 1] ; l++)
        {
          const int k = b->SpT->J[l];
          const int Tik = b->SpT->array[l];
          if (k == j1)
            aux *= D2f_reduction_map(b, x[k], Tik, k);
          else
            aux *= f_reduction_map(b, x[k], Tik, k);
      }
    }
  else
    {
      for (l = b->SpT->I[i] ; l < b->SpT->I[i + 1] ; l++)
        {
          const int k = b->SpT->J[l];
          const int Tik = b->SpT->array[l];
          if (k == j1 || k == j2)
            aux *= Df_reduction_map(b, x[k], Tik, k);
          else
            aux *= f_reduction_map(b, x[k], Tik, k);
        }
    }
  return aux;
}

/**
 * @brief Compute the index of the cell containing x
 *
 * @param basis A local basis
 * @param x
 * @return int an integer between -1 and basis->nb_func - 1. The value -1 means that x lies outside of the domain.
 */
int pnl_basis_local_get_index(const PnlBasis *basis, const double *x)
{
  int i;
  int index_per_dim, global_index, n_intervals_prod;
  int *n_intervals = (int*) basis->f_params;
  global_index = 0;
  n_intervals_prod = 1;
  for (i = 0; i < basis->nb_variates; i++)
    {
      double map_xi = basis->map ? basis->map(x[i], i, basis->map_params) : x[i];
      if (basis->isreduced)
        {
          index_per_dim = (int) (((map_xi - basis->center[i]) * basis->scale[i] + 1) * n_intervals[i] / 2.);
        }
      else
        {
          index_per_dim = (int) ((map_xi + 1) * n_intervals[i] / 2.);
        }
      if (index_per_dim < 0 || index_per_dim >= n_intervals[i])
        {
          return -1;
        }
      global_index += index_per_dim * n_intervals_prod;
      n_intervals_prod *= n_intervals[i];
    }
    return global_index;
}

/**
 * Evaluate a local basis at x
 *
 * @param coef a vector of weights
 * @param x the coordinates of the point at which to evaluate the function
 * @param basis a PnlBasis
 *
 * @return sum (coef .* f(x))
 */
static double pnl_basis_local_eval(const PnlBasis *basis, const PnlVect *coef, const double *x)
{
  CHECK_BASIS_TYPE(PNL_BASIS_LOCAL, basis);
  int global_index = pnl_basis_local_get_index(basis, x);
  if (global_index < 0)
    {
      return 0.;
    }
  return GET(coef, global_index);
}

/**
 * Evaluate a linear combination of basis functions at x
 *
 * @deprecated Use the function pnl_basis_eval_vect(const PnlBasis *basis, const PnlVect *coef, const PnlVect *x)
 *
 * @param coef a vector of weights
 * @param x the coordinates of the point at which to evaluate the function
 * @param basis a PnlBasis
 *
 * @return sum (coef .* f(x))
 */
double pnl_basis_eval(const PnlBasis *basis, const PnlVect *coef, const double *x)
{
  int i;
  double y;

  CHECK_NB_FUNC(coef, basis);

  if (basis->id == PNL_BASIS_LOCAL)
    {
      return pnl_basis_local_eval(basis, coef, x);
    }
  y = 0.;
  for (i = 0 ; i < coef->size ; i++)
    {
      const double a = PNL_GET(coef, i);
      if (a != 0)
        {
          y += a * pnl_basis_i(basis, x, i);
        }
    }
  return y;
}

/**
 * Evaluate the first derivative with respect to x[i] of a linear combination
 * of basis functions at x
 *
 * @deprecated Use the function pnl_basis_eval_D_vect(const PnlBasis *basis, const PnlVect *coef, const PnlVect *x, int i)
 *
 * @param coef a vector typically computed by pnl_basis_fit_ls
 * @param x the coordinates of the point at which to evaluate the function
 * @param basis a PnlBasis
 * @param i the index with respect to which the derivative is computed
 *
 * @return sum (coef .* D_i f(x))
 */
double pnl_basis_eval_D(const PnlBasis *basis, const PnlVect *coef, const double *x, int i)
{
  int k;
  double y;

  CHECK_IS_DIFFERENTIABLE(basis);
  CHECK_NB_FUNC(coef, basis);
  y = 0.;
  for (k = 0 ; k < coef->size ; k++)
    {
      const double a = PNL_GET(coef, k);
      if (a != 0.)
        {
          y += a * pnl_basis_i_D(basis, x, k, i);
        }
    }
  return y;
}

/**
 * Evaluate the second derivative with respect to x[i] and x[j] of a linear
 * combination of basis functions at x
 *
 * @deprecated Use the function pnl_basis_eval_D2_vect(const PnlBasis *basis, const PnlVect *coef, const PnlVect *x, int i, int j)
 *
 * @param coef a vector typically computed by pnl_basis_fit_ls
 * @param x the coordinates of the point at which to evaluate the function
 * @param basis a PnlBasis
 * @param i the index with respect to which the derivative is computed
 * @param j the index with respect to which the derivative is computed
 *
 * @return sum (coef .* D2_{i,j} f(x))
 */
double pnl_basis_eval_D2(const PnlBasis *basis, const PnlVect *coef, const double *x, int i, int j)
{
  int k;
  double y;

  CHECK_IS_DIFFERENTIABLE(basis);
  CHECK_NB_FUNC(coef, basis);
  y = 0.;
  for (k = 0 ; k < coef->size ; k++)
    {
      const double a = PNL_GET(coef, k);
      if (a != 0.)
        {
          y += a * pnl_basis_i_D2(basis, x, k, i, j);
        }
    }
  return y;
}

/**
 * Evaluate the function, its gradient and Hessian matrix at x. The function is
 * defined by the linear combination sum (coef .* f(x))
 *
 * @deprecated Use the function pnl_basis_eval_derivs_vect(const PnlBasis *b, const PnlVect *coef, const PnlVect *x, double *val, PnlVect *grad, PnlMat *hes)
 *
 * @param coef a vector typically computed by pnl_basis_fit_ls
 * @param x the coordinates of the point at which to evaluate the function
 * @param b PnlBasis
 * @param val contains the value of sum (coef .* f(x)) on exit
 * @param grad contains the value of sum (coef .* Df(x)) on exit
 * @param hes contains the value of sum (coef .* D2f(x)) on exit
 *
 */
void pnl_basis_eval_derivs(const PnlBasis *b, const PnlVect *coef, const double *x,
                           double *val, PnlVect *grad, PnlMat *hes)
{
  int i, j, l, n, m;
  double y, *f, *Df, D2f;

  CHECK_IS_DIFFERENTIABLE(b);
  CHECK_NB_FUNC(coef, b);
  y = 0.;
  n = b->nb_variates;
  f = malloc(sizeof(double) * n);
  Df = malloc(sizeof(double) * n);
  pnl_vect_resize(grad, n);
  pnl_mat_resize(hes, n, n);
  pnl_vect_set_zero(grad);
  pnl_mat_set_zero(hes);
  for (i = 0 ; i < coef->size ; i++)
    {
      double auxf;
      const double a = PNL_GET(coef, i);
      if (a == 0.) continue;
      auxf = 1;
      /*
       * computation of val
       */
      for (l = b->SpT->I[i] ; l < b->SpT->I[i + 1] ; l++)
        {
          const int k = b->SpT->J[l];
          const int Tik = b->SpT->array[l];
          f[k] = f_reduction_map(b, x[k], Tik, k);
          auxf *= f[k];
        }
      y += a * auxf;
      /*
       * computation of the gradient and the Hessian matrix
       */
      for (j = 0 ; j < n ; j++)
        {
          auxf = 1;
          for (l = b->SpT->I[i] ; l < b->SpT->I[i + 1] ; l++)
            {
              const int k = b->SpT->J[l];
              if (k != j) auxf *= f[k];
            }
            const int Tij = pnl_sp_mat_int_get(b->SpT, i, j);
            /* gradient */
            if (Tij >= 1)
              {
                Df[j] = Df_reduction_map(b, x[j], Tij, j);
                PNL_LET(grad, j) += a * auxf * Df[j];
              }
            else
              Df[j] = 0.;

            /* diagonal terms of the Hessian matrix */
            if (Tij >= 2)
              {
                D2f = D2f_reduction_map(b, x[j], Tij, j);
                PNL_MLET(hes, j, j) += a * auxf * D2f;
              }
          /* non diagonal terms of the Hessian matrix */
          for (m = 0 ; m < j ; m++)
            {
              if (Df[j] == 0. || Df[m] == 0.) continue;
              auxf = 1;
              for (l = b->SpT->I[i] ; l < b->SpT->I[i + 1] ; l++)
                {
                  const int k = b->SpT->J[l];
                  if ((k != j) && (k != m)) auxf *= f[k];
                }

              PNL_MLET(hes, j, m) += a * auxf * Df[j] * Df[m];
              PNL_MLET(hes, m, j) = PNL_MGET(hes, j, m);
            }
        }
    }
  *val = y;
  free(f);
  free(Df);
}

/**
 * Find the best approximation of the function defined by f(x(i,:)) = y(i)
 * General purpose function
 *
 * @param basis a PnlBasis
 * @param x the matrix of points at which we know the value of the function. One line
 * of the matrix is the vector of the coordinates of one point
 * @param y the values of the function f at the points defined by x
 * @param coef contains on exit the coefficients of the regression
 *
 * @return PNL_OK or PNL_FAIL
 */
static int pnl_basis_fit_ls_general(const PnlBasis *basis, PnlVect *coef, const PnlMat *x, const PnlVect *y)
{
  int N, i, k;
  double b_k;
  PnlMat *A;
  PnlVect *phi_k;

  N = y->size;
  pnl_vect_resize(coef, basis->nb_func);
  pnl_vect_set_all(coef, 0.);
  phi_k = pnl_vect_create_from_scalar(basis->nb_func, 0.);
  A = pnl_mat_create_from_scalar(basis->nb_func, basis->nb_func, 0.);

  /* construct A and b*/
  for (i = 0 ; i < N ; i++)
    {
      for (k = 0 ; k < basis->nb_func ; k++)
        {
          PnlVect xi = pnl_vect_wrap_mat_row(x, i);
          const double tmp = pnl_basis_i_vect(basis, &xi, k);
          b_k =  PNL_GET(coef, k);
          b_k += tmp * PNL_GET(y, i);
          PNL_LET(coef, k) = b_k;
          PNL_LET(phi_k, k) = tmp;
        }
      /* A += phi_k' * phi_k */
      pnl_mat_dger(1., phi_k, phi_k, A);
    }

  /* Because A often comes from simulation, A is not >0. So we use a
   * least-square approach
   */
  pnl_mat_ls(A, coef);

  pnl_vect_free(&phi_k);
  pnl_mat_free(&A);

  return PNL_OK;
}

/**
 * Find the best approximation of the function defined by f(x(i,:)) = y(i)
 * For local basis only
 *
 * @param basis a PnlBasis
 * @param x the matrix of points at which we know the value of the function. One line
 * of the matrix is the vector of the coordinates of one point
 * @param y the values of the function f at the points defined by x
 * @param coef contains on exit the coefficients of the regression
 *
 * @return PNL_OK or PNL_FAIL
 */
static int pnl_basis_fit_ls_local(const PnlBasis *basis, PnlVect *coef, const PnlMat *x, const PnlVect *y)
{
  int i, k;
  int *count;
  CHECK_BASIS_TYPE(PNL_BASIS_LOCAL, basis);
  pnl_vect_resize(coef, basis->nb_func);
  pnl_vect_set_all(coef, 0.);
  count = calloc(coef->size, sizeof(int));


  for (i = 0; i < x->m; i++)
    {
      PnlVect xi = pnl_vect_wrap_mat_row(x, i);
      int x_index = pnl_basis_local_get_index(basis, xi.array);
      if (x_index >= 0)
        {
          LET(coef, x_index) += GET(y, i);
          count[x_index] += 1;
        }
    }

  for (k = 0 ; k < coef->size; k++)
    {
      PNL_LET(coef, k) /= count[k];
    }
  free(count);
  return PNL_OK;
}

/**
 * Find the best approximation of the function defined by f(x(i,:)) = y(i)
 * For orthogonal basis only
 *
 * @param basis a PnlBasis
 * @param x the matrix of points at which we know the value of the function. One line
 * of the matrix is the vector of the coordinates of one point
 * @param y the values of the function f at the points defined by x
 * @param coef contains on exit the coefficients of the regression
 *
 * @return PNL_OK or PNL_FAIL
 */
static int pnl_basis_fit_ls_orthogonal(const PnlBasis *basis, PnlVect *coef, const PnlMat *x, const PnlVect *y)
{
  int N, i, k;
  N = y->size;
  pnl_vect_resize(coef, basis->nb_func);
  pnl_vect_set_all(coef, 0.);

  for (k = 0 ; k < basis->nb_func ; k++)
    {
      double sum, norm;
      sum = norm = 0.;
      for (i = 0 ; i < N ; i++)
        {
          PnlVect xi = pnl_vect_wrap_mat_row(x, i);
          const double val = pnl_basis_i_vect(basis, &xi, k);
          sum += val * PNL_GET(y, i);
          norm += val * val;
        }
      if (norm != 0.)
        {
          PNL_LET(coef, k) = sum / norm;
        }
    }
  return PNL_OK;
}

/**
 * Find the best approximation of the function defined by f(x(i,:)) = y(i)
 *
 * @param basis a PnlBasis
 * @param x the matrix of points at which we know the value of the function. One line
 * of the matrix is the vector of the coordinates of one point
 * @param y the values of the function f at the points defined by x
 * @param coef contains on exit the coefficients of the regression
 *
 * @return PNL_OK or PNL_FAIL
 */
int pnl_basis_fit_ls(const PnlBasis *basis, PnlVect *coef, const PnlMat *x, const PnlVect *y)
{
  if (basis->id == PNL_BASIS_LOCAL)
    {
      return pnl_basis_fit_ls_local(basis, coef, x, y);
    }
  if (PnlBasisTypeTab[basis->id].is_orthogonal == PNL_TRUE)
    {
      return pnl_basis_fit_ls_orthogonal(basis, coef, x, y);
    }
  return pnl_basis_fit_ls_general(basis, coef, x, y);
}

/**
 * An element of a basis writes as a product
 *      p_1(x_1) p2(x_2) .... p_n(x_n)
 * for a polynomial with n variates. Each p_k is a polynomial with only one variate.
 *
 * This functions evaluates the term p_k of the i-th element of the basis b at the point x
 *
 * @param b a PnlBasis
 * @param x a PnlVect containing the coordinates of the point at which to evaluate the basis
 * @param i an integer describing the index of the element of the basis to consider
 * @param k the index of the term to be evaluated with element i of the basis
 *
 * @return (f_i)_k (x)
 */
double pnl_basis_ik_vect(const PnlBasis *b, const PnlVect *x, int i, int k)
{
  return pnl_basis_ik(b, x->array, i, k);
}

/**
 * Evaluate the i-th element of the basis b at the point x
 *
 * @param b a PnlBasis
 * @param x a PnlVect containing the coordinates of the point at which to evaluate the basis
 * @param i an integer describing the index of the element of the basis to consider
 *
 * @return f_i(x) where f is the i-th basis function
 */
double pnl_basis_i_vect(const PnlBasis *b, const PnlVect *x, int i)
{
  if (i > b->len_T - 1)
    {
      return PNL_EVAL_RNFUNCR(&(b->func_list[i-b->len_T]), x);
    }
  return pnl_basis_i(b, x->array, i);
}

/**
 * First order derivative
 *
 * @param b a basis
 * @param x the point at which to evaluate the first derivative
 * @param i the index of the basis element to differentiate
 * @param j the index of the variable w.r.t which we differentiate
 *
 * @return (D(b_i)/Dj)(x)
 */
double pnl_basis_i_D_vect(const PnlBasis *b, const PnlVect *x, int i, int j)
{
  return pnl_basis_i_D(b, x->array, i, j);
}

/**
 * Second order derivative
 *
 * @param b a basis
 * @param x the point at which to evaluate the first derivative
 * @param i the index of the basis element to differentiate
 * @param j1 the index of the first variable w.r.t which we differentiate
 * @param j2 the index of the second variable w.r.t which we differentiate
 *
 * @return (D(b_i)/(Dj1 Dj2))(x)
 */
double pnl_basis_i_D2_vect(const PnlBasis *b, const PnlVect *x, int i, int j1, int j2)
{
  return pnl_basis_i_D2(b, x->array, i, j1, j2);
}

/**
 * Evaluate a linear combination of basis functions at x
 *
 * @param coef a vector typically computed by pnl_basis_fit_ls
 * @param x the coordinates of the point at which to evaluate the function
 * @param basis a PnlBasis
 *
 * @return sum (coef .* f(x))
 */
double pnl_basis_eval_vect(const PnlBasis *basis, const PnlVect *coef, const PnlVect *x)
{
  return pnl_basis_eval(basis, coef, x->array);
}

/**
 * Evaluate the first derivative with respect to x[i] of a linear combination
 * of basis functions at x
 *
 * @param coef a vector typically computed by pnl_basis_fit_ls
 * @param x the coordinates of the point at which to evaluate the function
 * @param basis a PnlBasis
 * @param i the index with respect to which the derivative is computed
 *
 * @return sum (coef .* D_i f(x))
 */
double pnl_basis_eval_D_vect(const PnlBasis *basis, const PnlVect *coef, const PnlVect *x, int i)
{
  return pnl_basis_eval_D(basis, coef, x->array, i);
}

/**
 * Evaluate the second derivative with respect to x[i] and x[j] of a linear
 * combination of basis functions at x
 *
 * @param coef a vector typically computed by pnl_basis_fit_ls
 * @param x the coordinates of the point at which to evaluate the function
 * @param basis a PnlBasis
 * @param i the index with respect to which the derivative is computed
 * @param j the index with respect to which the derivative is computed
 *
 * @return sum (coef .* D2_{i,j} f(x))
 */
double pnl_basis_eval_D2_vect(const PnlBasis *basis, const PnlVect *coef, const PnlVect *x, int i, int j)
{
  return pnl_basis_eval_D2(basis, coef, x->array, i, j);
}

/**
 * Evaluate the function, its gradient and Hessian matrix at x. The function is
 * defined by  linear combination sum (coef .* f(x))
 *
 * @param coef a vector typically computed by pnl_basis_fit_ls
 * @param x the coordinates of the point at which to evaluate the function
 * @param b PnlBasis
 * @param val contains the value of sum (coef .* f(x)) on exit
 * @param grad contains the value of sum (coef .* Df(x)) on exit
 * @param hes contains the value of sum (coef .* D2 f(x)) on exit
 *
 */
void pnl_basis_eval_derivs_vect(const PnlBasis *b, const PnlVect *coef, const PnlVect *x,
                                double *val, PnlVect *grad, PnlMat *hes)
{
  pnl_basis_eval_derivs(b, coef, x->array, val, grad, hes);
}

