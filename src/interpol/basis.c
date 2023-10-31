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
#define CHECK_NB_VARIATES(x, basis)                                               \
  if ( x->size != basis->nb_variates )                                            \
    {                                                                             \
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

#if 0
/**
 * Compute the maximum degree of the last component.
 * if the total degree of the other elements is @p partial.
 *
 * The total degree function is the sum of the partial degrees
 *
 * @param total total degree
 * @param partial partial degree already used
 */
static int freedom_degree_sum(int total, int partial)
{
  if (partial > total) return -1;
  return total - partial;
}
#endif

/**
 * Compute the maximum degree of the last component
 * if the total degree of the other elements is @p partial.
 *
 * The total degree function is the product of the partial degrees
 *
 * @param total total degree
 * @param partial partial degree already used
 */
static int freedom_degree_prod(int total, int partial)
{
  if (partial > total) return -1;
  return total / MAX(partial, 1);
}

/**
 * Return the total degree (sum of the partial degrees) of the polynomial
 * represented by line i of T
 *
 * @param T a matrix of integers representing the decomposition of the mutli-d
 * polynomials
 * @param i the index of the element to be considered in the basis (i.e. the row
 * of T to consider)
 * @return  the total degree of the polynomials represented by
 * line i of T
 */
static int count_sum_degree(const PnlMatInt *T, int i)
{
  int j, deg;
  deg = 0;

  for (j = 0 ; j < T->n ; j++)
    {
      deg += PNL_MGET(T, i, j);
    }
  return deg;
}

/**
 * Return the total degree (product of the partial degrees) of the polynomial
 * represented by line i of T
 *
 * The total degree is the product of MAX(1, d_i) where d_i is the partial
 * degree. If all d_i are zeros, the total degree is 0.
 *
 * @param T a matrix of integers representing the decomposition of the mutli-d
 * polynomials
 * @param i the index of the element to be considered in the basis (i.e. the row
 * of T to consider)
 * @return  the total degree of the polynomials represented by
 * line i of T
 */
static int count_prod_degree(const PnlMatInt *T, int i)
{
  int j, deg, all_zero = PNL_TRUE;
  deg = 1;

  for (j = 0 ; j < T->n ; j++)
    {
      const int power = PNL_MGET(T, i, j);
      if (all_zero == PNL_TRUE && power > 0) all_zero = PNL_FALSE;
      deg *= MAX(power, 1);
    }
  if (all_zero == PNL_TRUE) return 0;
  return deg;
}

/**
 * Return the total hyperbolic degree at the power q of the polynomial
 * represented by line i of T
 *
 * @param T a matrix of integers representing the decomposition of the mutli-d
 * polynomials
 * @param i the index of the element to be considered in the basis (i.e. the row
 * of T to consider)
 * @param q the hyperbolic index
 * @return  the total degree of the polynomials represented by
 * line i of T
 */
static double count_hyperbolic_degree(const PnlMatInt *T, int i, double q)
{
  int j;
  double deg_q;

  for (j = 0, deg_q = 0. ; j < T->n ; j++)
    {
      deg_q += pow(PNL_MGET(T, i, j), q);
    }
  return deg_q;
}

/**
 * Copy T_prev(T_previ, :) into T(Ti,:)
 *
 * @param T an integer matrix with more columns than that T_prev
 * @param T_prev an integer matrix containing the last tensor computed
 * @param Ti the index of the line to consider in T
 * @param T_previ the index of the line to consider in T_prev
 *
 */
static void copy_previous_tensor(PnlMatInt *T, PnlMatInt *T_prev, int Ti, int T_previ)
{
  int j;
  for (j = 0 ; j < T_prev->n ; j++)
    {
      PNL_MLET(T, Ti, j) = PNL_MGET(T_prev, T_previ, j);
    }
}

/**
 * Compute the integer matrix representing the decomposition as a tensor
 * product of the elements of * the multi-b basis onto the 1-d basis
 *
 * @param nb_variates the number of variates
 * @param nb_func the number of elements of the basis
 *
 * @return a tensor
 */
static PnlMatInt *compute_tensor(int nb_func, int nb_variates)
{
  PnlMatInt *T;
  T = pnl_mat_int_create(nb_func, nb_variates);
  if (nb_variates == 1)
    {
      int i;
      for (i = 0 ; i < nb_func ; i++)
        {
          PNL_MLET(T, i, 0) = i;
        }
      return T;
    }
  else
    {
      int current_degree, deg;
      int block, block_start;
      int i, j, k;
      PnlMatInt *T_prev;
      current_degree = 0;
      T_prev = compute_tensor(nb_func, nb_variates - 1);

      i = 0;
      while (1) /* loop on global degree */
        {
          /* determining the last element of global degree <= current_degree */
          j = 0;
          while (j < T_prev->m && count_sum_degree(T_prev, j) < current_degree + 1)
            {
              j++;
            }
          j--;
          for (k = j, deg = current_degree ; k >= 0 ; deg--)
            {
              block_start = k;
              if (deg > 0)
                {
                  /* Find the beginning of the block with global degree deg */
                  while (count_sum_degree(T_prev, block_start) == deg)
                    {
                      block_start--;
                    }
                  block_start++;
                }
              /* Loop in an increasing order over the block of global degree deg */
              for (block = block_start ; block <= k ; block++ , i++)
                {
                  if (i == nb_func)
                    {
                      pnl_mat_int_free(&T_prev);
                      return T;
                    }
                  copy_previous_tensor(T, T_prev, i, block);
                  PNL_MLET(T, i, nb_variates - 1) = current_degree - deg;
                }
              /* Consider the previous block */
              k = block_start - 1;
            }
          current_degree++;
        }
    }
}

/**
 * Compute the number of elements with total degree less or equal than degree
 * in the basis with (T->n + 1) variates
 *
 * @param T the tensor matrix of the basis with n-1 variates
 * @param degree the maximum total degree requested
 * @param count_degree a function to compute the total of a given line in a
 * tensor
 * @param freedom_degree a function to compute the number of degrees of
 * freedom
 *
 * @return the number of elements with total degree less or equal than degree in
 * the basis with n variates
 */
static int compute_nb_elements(const PnlMatInt *T, int degree,
                               int (*count_degree)(const PnlMatInt *, int),
                               int (*freedom_degree)(int, int))
{
  int i;
  int total_elements; /* Number of elements of total degree smaller than degree */
  total_elements = 0;
  for (i = 0 ; i < T->m ; i++)
    {
      int deg = count_degree(T, i);
      total_elements += freedom_degree(degree, deg) + 1;
    }
  return total_elements;
}

/**
 * Compute the tensor matrix of the nb_variates variate basis with a total degree less or
 * equal than degree. The total degree is defined by the function
 * count_degree
 *
 * @param nb_variates the number of variates of the basis.
 * @param degree the total degree
 * @param count_degree a function to compute the total of a given line in a
 * tensor
 * @param freedom_degree a function to compute the number of degrees of
 * freedom
 *
 * @return the tensor matrix of the nb_variates variate basis with a total degree less or
 * equal than degree
 */
static PnlMatInt *compute_tensor_from_degree_function(int degree, int nb_variates,
    int (*count_degree)(const PnlMatInt *, int),
    int (*freedom_degree)(int total, int partial))
{
  PnlMatInt *T;
  if (nb_variates <= 0)
    {
      printf("Nb of variates must be strictly positive in compute_tensor_from_degree\n");
      abort();
    }
  if (nb_variates == 1)
    {
      int i;
      T = pnl_mat_int_create(degree + 1, nb_variates);
      for (i = 0 ; i < degree + 1 ; i++)
        {
          PNL_MLET(T, i, 0) = i;
        }
      return T;
    }
  else
    {
      int nb_elements;
      int current_degree, deg;
      int block_end, block_start;
      int line, i, k;
      PnlMatInt *T_prev;
      current_degree = 0;
      /* Compute the tensor with one variate less */
      T_prev = compute_tensor_from_degree_function(degree, nb_variates - 1, count_degree, freedom_degree);
      /* Compute the number of rows of T */
      nb_elements = compute_nb_elements(T_prev, degree, count_degree, freedom_degree);
      T = pnl_mat_int_create(nb_elements, nb_variates);
      line = 0;
      /* loop on global degree */
      for (current_degree = 0 ; current_degree <= degree ; current_degree++)
        {
          /* loop on  partial degree between 0 and current_degree */
          for (deg = 0 ; deg <= current_degree ; deg++)
            {
              block_start = 0;
              block_end = 0;
              /* Determine the first element with global degree = deg */
              while (block_start < T_prev->m && count_degree(T_prev, block_start) < deg)
                {
                  block_start++;
                }
              block_end = block_start;
              /* Find the last element of the block with global degree = deg */
              while (block_end < T_prev->m && count_degree(T_prev, block_end) == deg)
                {
                  block_end++;
                }
              block_end--;

              /* loop on the degree of the extra dimension */
              for (i = freedom_degree(current_degree - 1, deg)  + 1; i <= freedom_degree(current_degree, deg) ; i++)
                {

                  for (k = block_start ; k <= block_end ; k++, line++)
                    {
                      copy_previous_tensor(T, T_prev, line, k);
                      PNL_MLET(T, line, nb_variates - 1) = i;
                    }
                }
            }
        }
      pnl_mat_int_free(&T_prev);
      return T;
    }
}

/**
 * Compute the tensor representation
 *
 * @param[in,out] T
 * @param degree
 * @param nb_variates
 *
 * @return the number of rows already computed in T
 */
static int compute_tensor_from_sum_degree_rec(PnlMatInt *T, int degree, int nb_variates)
{
  if (nb_variates <= 0)
    {
      printf("Nb of variates must be strictly positive in compute_tensor_from_degree\n");
      abort();
    }
  if (nb_variates == 1)
    {
      int i;
      pnl_mat_int_set_zero(T);
      for (i = 0 ; i < degree + 1 ; i++)
        {
          PNL_MLET(T, i, 0) = i;
        }
      return degree + 1;
    }
  else
    {
      int l, deg, current_row;
      int rows_rec = compute_tensor_from_sum_degree_rec(T, degree, nb_variates - 1);
      current_row = rows_rec;
      for (l = 0; l < rows_rec; l++)
        {
          int partial_degree = count_sum_degree(T, l);
          for (deg = 1; deg <= degree - partial_degree; deg++)
            {
              PnlVectInt row_l = pnl_vect_int_wrap_mat_row(T, l);
              pnl_mat_int_set_row(T, &row_l, current_row);
              PNL_MLET(T, current_row, nb_variates - 1) = deg;
              current_row ++;
            }
        }
      return current_row;
    }
}

static PnlMatInt *compute_tensor_from_sum_degree(int degree, int nb_variates)
{
  int nb_elements = (int) pnl_round(pnl_sf_choose(nb_variates + degree, degree));
  PnlMatInt *T = pnl_mat_int_create(nb_elements, nb_variates);
  compute_tensor_from_sum_degree_rec(T, degree, nb_variates);
  return T;
}

static PnlMatInt *compute_tensor_from_prod_degree(int degree, int nb_variates)
{
  return compute_tensor_from_degree_function(degree, nb_variates, count_prod_degree, freedom_degree_prod);
}

/**
 * Compute the tensor matrix of the n-variate basis with a
 * hyperbolic degree of order q less or equal than degree
 *
 * @param n the number of variates of the basis.
 * @param q the hyperbolic index
 * @param degree the total hyperbolic degree
 *
 * @return the tensor matrix of the n-variate basis with a total degree less or
 * equal than degree
 */
static PnlMatInt *compute_tensor_from_hyperbolic_degree(double degree, double q, int n)
{
  int i, i_sparse;
  double degree_q;
  PnlMatInt *T;
  T = compute_tensor_from_sum_degree((int)ceil(degree), n);
  degree_q = pow(degree, q);
  for (i = 0, i_sparse = 0 ; i < T->m ; i++)
    {
      if (count_hyperbolic_degree(T, i, q) > degree_q) continue;
      if (i_sparse < i)
        {
          memcpy(T->array + i_sparse * n, T->array + i * n, n * sizeof(int));
        }
      i_sparse++;
    }
  pnl_mat_int_resize(T, i_sparse, n);
  return T;
}

/**
 * Compute the tensor for a tensor local basis
 * 
 * @param space_dim the dimension of the state space
 * @param n_intervals this is an array of size \a space_dim describing the number of intervals for every dimension
 * @return PnlMatInt
 */
static PnlMatInt* compute_tensor_permutation(int space_dim, int *n_elements)
{
  int permutation_length, i, col;
  int inner_loop_length, outer_loop_length;
  PnlMatInt *T;
  permutation_length = 1;
  for (i = 0; i < space_dim; i++)
    {
      permutation_length *= n_elements[i];
    }
  T = pnl_mat_int_create(permutation_length, space_dim);
  outer_loop_length = permutation_length;
  inner_loop_length = 1;
  for (col = 0; col < space_dim; col++)
    {
      int outer_index;
      outer_loop_length /= n_elements[col];
      for (outer_index = 0; outer_index < outer_loop_length; outer_index++)
        {
          int middle_loop;
          for (middle_loop = 0; middle_loop < n_elements[col]; middle_loop++)
            {
              int inner_loop;
              for (inner_loop = 0; inner_loop < inner_loop_length; inner_loop++)
                {
                  int row = outer_index * inner_loop_length * n_elements[col] + middle_loop * inner_loop_length + inner_loop;
                  PNL_MLET(T, row, col) = middle_loop + 1;
                }
            }
        }
      inner_loop_length *= n_elements[col];
    }
  return T;
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
 * The terminal recursive function to compute Tchebychev polynomials of any
 * order.
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
 *  @param f_n_1 used to store the derivative of  the polynomial of order
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
  o->T = NULL;
  o->SpT = NULL;
  o->isreduced = 0;
  o->center = NULL;
  o->scale = NULL;
  o->func_list = NULL;
  o->len_func_list = 0;
  o->len_T = 0;
  o->params = NULL;
  o->params_size = 0;
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
 * Create a PnlBasis and stores it into its first argument
 *
 * @param b an already allocated basis (as returned by pnl_basis_new for instance)
 * @param T the tensor of the multi-dimensional basis. No copy of T is done, so
 * do not free T. It will be freed transparently by pnl_basis_free
 */
void  pnl_basis_set_from_tensor(PnlBasis *b, const PnlMatInt *T)
{
  b->nb_func = T->m;
  b->nb_variates = T->n;
  if(b->func_list) free(b->func_list);
  b->len_func_list = 0;
  b->len_T = T->m;

  /* Not sure this is the right place to put it */
  pnl_mat_int_free(&(b->T));
  pnl_sp_mat_int_free(&(b->SpT));
  if (b->isreduced == 1)
    {
      b->isreduced = 0;
      free(b->center);
      b->center = NULL;
      free(b->scale);
      b->scale = NULL;
    }

  b->T = (PnlMatInt *) T;
  b->SpT = pnl_sp_mat_int_create_from_mat(T);
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
  PnlMatInt *T;
  T = compute_tensor(nb_func, nb_variates);
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
  T = compute_tensor_from_hyperbolic_degree(degree, q, n);
  return pnl_basis_create_from_tensor(index, T);
}

/**
 * Return a PnlBasis with local and orthogonal function
 *
 * @param n_intervals this is an array of size \a space_dim describing the number of intervals for every dimension
 * @param space_dim the dimension of the state space
 * @return PnlBasis
 */
PnlBasis* pnl_basis_local_create(int *n_intervals, int space_dim)
{
  PnlMatInt *T;
  PnlBasis *B;
  T = compute_tensor_permutation(space_dim, n_intervals);
  B = pnl_basis_create_from_tensor(PNL_BASIS_LOCAL, T);
  B->params_size = space_dim * sizeof(int);
  B->params = n_intervals;
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
  pnl_basis_set_from_tensor(dest, src->T);
  dest->id = src->id;
  dest->label = src->label;
  dest->f = src->f;
  dest->Df = src->Df;
  dest->D2f = src->D2f;
  dest->params = realloc(dest->params, src->params_size);
  dest->params_size = src->params_size;
  if (src->isreduced == 1)
    {
      int n = dest->nb_variates;
      dest->isreduced = 0;
      /* No need to free, basis_set_from_tensor does it */
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
  pnl_mat_int_del_row(B->T, i);
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
      PnlVectInt Ti = pnl_vect_int_wrap_mat_row(B->T, i);
      if (pnl_vect_int_eq(&Ti, d) == PNL_TRUE) break;
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
  pnl_mat_int_add_row(B->T, B->T->m, d);
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

  if (B->center != NULL) free(B->center);
  if (B->scale != NULL) free(B->scale);
  n = B->nb_variates;
  B->center = malloc(n * sizeof(double));
  B->scale = malloc(n * sizeof(double));
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

  if (B->center != NULL) free(B->center);
  if (B->scale != NULL) free(B->scale);
  n = B->nb_variates;
  B->center = malloc(n * sizeof(double));
  B->scale = malloc(n * sizeof(double));
  B->isreduced = 1;
  memcpy(B->center, center->array, n * sizeof(double));
  for (i = 0 ; i < n ; i++)
    {
      B->scale[i] = 1. / PNL_GET(scale, i);
    }
}

/**
 * Free a PnlBasis
 *
 * @param B
 */
void pnl_basis_free(PnlBasis **B)
{
  if (*B == NULL) return;
  pnl_mat_int_free(&((*B)->T));
  pnl_sp_mat_int_free(&(*B)->SpT);
  if ((*B)->params_size > 0) {
    free((*B)->params);
    (*B)->params = NULL;
    (*B)->params_size = 0;
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
  printf("\tExtra parameters size: %zu\n", B->params_size);
  printf("\tNumber of functions in tensor: %d\n", B->len_T);;
  printf("\tNumber of extra functions: %d\n", B->len_func_list);
  printf("\tTotal number of functions: %d\n", B->nb_func);
  printf("\tisreduced = %d\n", B->isreduced);
  if (B->isreduced)
    {
      int i;
      printf("\tcenter = ");
      for (i = 0 ; i < B->nb_variates ; i++) printf("%f ", B->center[i]);
      printf("\n\tscale = ");
      for (i = 0 ; i < B->nb_variates ; i++) printf("%f ", B->scale[i]);
      printf("\n");
    }
  printf("\tTensor matrix : \n");
  pnl_mat_int_print(B->T);
  printf("\n");
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
  const int Tik = PNL_MGET(b->T, i, k);
  if (Tik == 0) return 1.;
  if (b->isreduced == 1)
    {
      return (b->f)((x[k] - b->center[k]) * b->scale[k], Tik, k, b->params);
    }
  else
    {
      return (b->f)(x[k], Tik, k, b->params);
    }
}

/*
 * pnl_basis_i is moved to the header pnl_basis.h to be decalred as inlined.
 * For low order polynomials, the gain is about 20%
 */

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

  /* Test if the partial degree is small enough to return 0 */
  if (PNL_MGET(b->T, i, j) == 0) return 0.;

  if (b->isreduced == 1)
    {
      for (l = b->SpT->I[i] ; l < b->SpT->I[i + 1] ; l++)
        {
          const int k = b->SpT->J[l];
          const int Tik = b->SpT->array[l];
          if (k == j)
            aux *= b->scale[k] * (b->Df)((x[k] - b->center[k]) * b->scale[k], Tik, k, b->params);
          else
            aux *= (b->f)((x[k] - b->center[k]) * b->scale[k], Tik, k, b->params);
        }
    }
  else
    {
      for (l = b->SpT->I[i] ; l < b->SpT->I[i + 1] ; l++)
        {
          const int k = b->SpT->J[l];
          const int Tik = b->SpT->array[l];
          if (k == j)
            aux *= (b->Df)(x[k], Tik, k, b->params);
          else
            aux *= (b->f)(x[k], Tik, k, b->params);
        }
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
  /* Test if the partial degree is small enough to return 0 */
  if ((j1 == j2) && PNL_MGET(b->T, i, j1) <= 1) return 0.;
  if (PNL_MGET(b->T, i, j1) == 0 || PNL_MGET(b->T, i, j2) == 0) return 0.;

  if (b->isreduced == 1)
    {
      if (j1 == j2)
        {
            for (l = b->SpT->I[i] ; l < b->SpT->I[i + 1] ; l++)
              {
                const int k = b->SpT->J[l];
                const int Tik = b->SpT->array[l];
              if (k == j1)
                aux *= b->scale[k] * b->scale[k] * (b->D2f)((x[k] - b->center[k]) * b->scale[k], Tik, k, b->params);
              else
                aux *= (b->f)((x[k] - b->center[k]) * b->scale[k], Tik, k, b->params);
            }
        }
      else
        {
          for (l = b->SpT->I[i] ; l < b->SpT->I[i + 1] ; l++)
            {
              const int k = b->SpT->J[l];
              const int Tik = b->SpT->array[l];
              if (k == j1 || k == j2)
                aux *= b->scale[k] * (b->Df)((x[k] - b->center[k]) * b->scale[k], Tik, k, b->params);
              else
                aux *= (b->f)((x[k] - b->center[k]) * b->scale[k], Tik, k, b->params);
            }
        }
    }
  else
    {
      if (j1 == j2)
        {
          for (l = b->SpT->I[i] ; l < b->SpT->I[i + 1] ; l++)
            {
              const int k = b->SpT->J[l];
              const int Tik = b->SpT->array[l];
              if (k == j1)
                aux *= (b->D2f)(x[k], Tik, k, b->params);
              else
                aux *= (b->f)(x[k], Tik, k, b->params);
            }
        }
      else
        {
          for (l = b->SpT->I[i] ; l < b->SpT->I[i + 1] ; l++)
            {
              const int k = b->SpT->J[l];
              const int Tik = b->SpT->array[l];
              if (k == j1 || k == j2)
                aux *= (b->Df)(x[k], Tik, k, b->params);
              else
                aux *= (b->f)(x[k], Tik, k, b->params);
            }
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
  int *n_intervals = (int*) basis->params;
  global_index = 0;
  n_intervals_prod = 1;
  for (i = 0; i < basis->nb_variates; i++)
    {
      if (basis->isreduced)
        {
          index_per_dim = (int) (((x[i] - basis->center[i]) * basis->scale[i] + 1) * n_intervals[i] / 2.);
        }
      else
        {
          index_per_dim = (int) ((x[i] + 1) * n_intervals[i] / 2.);
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
      if (b->isreduced == 1)
        {
          for (l = b->SpT->I[i] ; l < b->SpT->I[i + 1] ; l++)
            {
              const int k = b->SpT->J[l];
              const int Tik = b->SpT->array[l];
              f[k] = (b->f)((x[k] - b->center[k]) * b->scale[k], Tik, k, b->params);
              auxf *= f[k];
            }
        }
      else
        {
          for (l = b->SpT->I[i] ; l < b->SpT->I[i + 1] ; l++)
            {
              const int k = b->SpT->J[l];
              const int Tik = b->SpT->array[l];
              f[k] = (b->f)(x[k], Tik, k, b->params);
              auxf *= f[k];
            }
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
          if (b->isreduced == 1)
            {
              const int Tij = PNL_MGET(b->T, i, j);
              /* gradient */
              if (Tij >= 1)
                {
                  Df[j] = b->scale[j] * (b->Df)((x[j] - b->center[j]) * b->scale[j], Tij, j, b->params);
                  PNL_LET(grad, j) += a * auxf * Df[j];
                }
              else
                Df[j] = 0.;

              /* diagonal terms of the Hessian matrix */
              if (Tij >= 2)
                {
                  D2f = b->scale[j] * b->scale[j] * (b->D2f)((x[j] - b->center[j]) * b->scale[j], Tij, j, b->params);
                  PNL_MLET(hes, j, j) += a * auxf * D2f;
                }
            }
          else
            {
              const int Tij = PNL_MGET(b->T, i, j);
              /* gradient */
              if (Tij >= 1)
                {
                  Df[j] = (b->Df)(x[j], Tij, j, b->params);
                  PNL_LET(grad, j) += a * auxf * Df[j];
                }
              else
                Df[j] = 0.;

              /* diagonal terms of the Hessian matrix */
              if (Tij >= 2)
                {
                  D2f = (b->D2f)(x[j], Tij, j, b->params);
                  PNL_MLET(hes, j, j) += a * auxf * D2f;
                }
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

