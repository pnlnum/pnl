
/*************************************************************************/
/* Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>            */
/*                                                                       */
/* This program is free software; you can redistribute it and/or modify  */
/* it under the terms of the GNU General Public License as published by  */
/* the Free Software Foundation; either version 3 of the License, or     */
/* (at your option) any later version.                                   */
/*                                                                       */
/* This program is distributed in the hope that it will be useful,       */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/* GNU General Public License for more details.                          */
/*                                                                       */
/* You should have received a copy of the GNU General Public License     */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pnl_matrix.h"
#include "pnl_vector.h"
#include "pnl_basis.h"
#include "pnl_mathtools.h"

#define function exp


/*
 * example of how to use  pnl_fit_least_squares to regress on a basis.
 * regression of the exponential function on the grid [0:0.05:5]
 */
static void exp_regression2()
{
  int n, basis_name, basis_dim, space_dim;
  int i;
  double a, b, h, err;
  PnlMat *t;
  PnlVect *y;
  PnlVect *alpha;

  PnlBasis *basis;

  printf("\nRegression of the exponential function on [0:0.05:5]\n");

  alpha = pnl_vect_create (0);

  /* creating the grid */
  a=0.0; b=5.0;
  n = 100;
  h = (b-a)/n;
  t = pnl_mat_create_from_double (n+1, 1, h);
  pnl_mat_set (t, 0, 0, 0.0);
  pnl_mat_cumsum (t, 'r');

  /* creating the values of exp on the grid */
  y = pnl_vect_create (n+1);
  for ( i=0 ; i<n+1 ; i++ )
    {
      pnl_vect_set (y, i, function(pnl_mat_get(t, i, 0)));
    }


  basis_name = HERMITIAN; /* TCHEBYCHEV; */
  space_dim = 1; /* real valued basis */
  basis_dim = 5; /* number of elements in the basis */

  basis = pnl_basis_init (basis_name, basis_dim, space_dim);

  pnl_fit_least_squares (alpha, t, y, basis, basis_dim);
  printf("coefficients of the decomposition : ");
  pnl_vect_print (alpha);
  /* computing the infinity norm of the error */
  err = 0.;
  for (i=0; i<t->m; i++)
    {
      double tmp = function(pnl_mat_get(t, i, 0)) -
        pnl_basis_eval (alpha, pnl_mat_lget(t, i, 0), basis);
      if (fabs(tmp) > err) err = fabs(tmp);
    }
  printf ("L^infty error : %f \n", err);

  pnl_mat_free (&t);
  pnl_vect_free (&y);
  pnl_vect_free (&alpha);
}

static double function2d( double *x )
{
  /* return x[0]*x[0] + x[1]*x[1]; */
  return log (1+x[0]*x[0] + x[1]*x[1]);
}

static void regression_multid()
{
  int n, basis_name, basis_dim, space_dim;
  int i, j;
  double a, b, h, err;
  PnlMat *t;
  PnlVect *y;
  PnlVect *alpha;

  PnlBasis *basis;
  printf("\nMulti dimensional regression on of log (1+x[0]*x[0] + x[1]*x[1]\n");
  alpha = pnl_vect_create (0);

  /* creating the grid */
  a=-2.0; b=2.0;
  n = 100;
  h = (b-a)/n;
  t = pnl_mat_create ((n+1)*(n+1), 2);

  /* creating the values of exp on the grid */
  y = pnl_vect_create ((n+1)*(n+1));
  for (i=0; i<n+1; i++)
    {
      for (j=0; j<n+1; j++)
        {
          pnl_mat_set (t, i*(n+1)+j, 0, a + i * h);
          pnl_mat_set (t, i*(n+1)+j, 1, a + j * h);
          pnl_vect_set (y, i*(n+1)+j, function2d(pnl_mat_lget(t, i*(n+1)+j, 0)));
        }
    }


  basis_name = HERMITIAN;
  space_dim = 2; /* real valued basis */
  basis_dim = 10; /* number of elements in the basis */

  basis = pnl_basis_init (basis_name, basis_dim, space_dim);

  pnl_fit_least_squares (alpha, t, y, basis, basis_dim);
  printf("coefficients of the decomposition : ");
  pnl_vect_print (alpha);

  /* computing the infinity norm of the error */
  err = 0.;
  for (i=0; i<t->m; i++)
    {
      double tmp = function2d(pnl_mat_lget(t, i, 0)) -
        pnl_basis_eval (alpha, pnl_mat_lget(t, i, 0), basis);
      if (fabs(tmp) > err) err = fabs(tmp);
    }

  printf ("L^infty error : %f \n", err);
  pnl_mat_free (&t);
  pnl_vect_free (&y);
  pnl_vect_free (&alpha);
}


void basis_test ()
{
  /* exp_regression(); */
  exp_regression2();
  regression_multid();
}
