
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pnl_matrix.h"
#include "pnl_vector.h"
#include "pnl_basis.h"
#include "pnl_mathtools.h"
#include "pnl_random.h"

#define function exp


/*
 * example of how to use  pnl_basis_fit_ls to regress on a basis.
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

  printf("\n\n** Regression of the exponential function on [0:0.05:5] ** \n\n");

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

  basis = pnl_basis_create (basis_name, basis_dim, space_dim);

  pnl_basis_fit_ls (basis, alpha, t, y);
  printf("coefficients of the decomposition : ");
  pnl_vect_print (alpha);
  /* computing the infinity norm of the error */
  err = 0.;
  for (i=0; i<t->m; i++)
    {
      double tmp = function(pnl_mat_get(t, i, 0)) -
        pnl_basis_eval (basis,alpha, pnl_mat_lget(t, i, 0));
      if (fabs(tmp) > err) err = fabs(tmp);
    }
  printf ("L^infty error : %f \n", err);

  pnl_basis_free (&basis);
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
  int n, basis_name, nb_func, nb_variates, degree;
  int i, j;
  double a, b, h, err;
  PnlMat *t;
  PnlVect *y;
  PnlVect *alpha;

  PnlBasis *basis;
  printf("\n\n** Multi dimensional regression on log (1+x[0]*x[0] + x[1]*x[1]) **\n\n");
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
  nb_variates = 2; /* functions with values in R^2 */
  printf ("Creating basis by specifying the number of functions\n");
  nb_func = 15; /* number of elements in the basis */

  basis = pnl_basis_create (basis_name, nb_func, nb_variates);

  pnl_basis_fit_ls (basis, alpha, t, y);
  printf("coefficients of the decomposition : ");
  pnl_vect_print (alpha);

  /* computing the infinity norm of the error */
  err = 0.;
  for (i=0; i<t->m; i++)
    {
      double tmp = function2d(pnl_mat_lget(t, i, 0)) -
        pnl_basis_eval (basis,alpha, pnl_mat_lget(t, i, 0));
      if (fabs(tmp) > err) err = fabs(tmp);
    }

  printf ("L^infty error : %f \n", err);
  pnl_basis_free (&basis);

  printf ("\nCreating basis by specifying the total degree\n");
  degree = 4; /* total degree */
  basis = pnl_basis_create_with_degree (basis_name, degree, nb_variates);

  pnl_basis_fit_ls (basis, alpha, t, y);
  printf("coefficients of the decomposition : ");
  pnl_vect_print (alpha);

  /* computing the infinity norm of the error */
  err = 0.;
  for (i=0; i<t->m; i++)
    {
      double tmp = function2d(pnl_mat_lget(t, i, 0)) -
        pnl_basis_eval (basis,alpha, pnl_mat_lget(t, i, 0));
      if (fabs(tmp) > err) err = fabs(tmp);
    }

  printf ("L^infty error : %f \n", err);
  pnl_basis_free (&basis);

  pnl_mat_free (&t);
  pnl_vect_free (&y);
  pnl_vect_free (&alpha);
}


static double fonction_a_retrouver(double t, double x)
{
  return (1-t)*x*x/5;
}

static double derive_x_fonction_a_retrouver(double t, double x)
{
  return (1-t)*2*x/5;
}

static double derive_xx_fonction_a_retrouver(double t, double x)
{
  return (1-t)*2/5;
}

static double derive_t_fonction_a_retrouver(double t, double x)
{
  return -x*x/5;
}

static double derive_x_approx_fonction(PnlBasis *B, PnlVect *alpha, double t, double x)
{
  double arg[2];
  arg[0] = t; arg[1] = x;
  return pnl_basis_eval_D (B, alpha, arg, 2); 
}

static void derive_approx_fonction(PnlBasis *B, PnlVect *D, PnlVect *alpha, double t, double x)
{
  double sum0, sum1, sum2, sum3;
  double arg[2];
  arg[0] = t; arg[1] = x;
  sum0=0.0;//calcule la valeur de la fonction
  sum1=0.0;//calcule la valeur de sa derivee en x
  sum2=0.0;//calcule la valeur de sa derivee seconde en x
  sum3=0.0;//calcule la valeur de sa derivee en t

  sum0 = pnl_basis_eval (B, alpha, arg);
  sum1 = pnl_basis_eval_D (B, alpha, arg, 2);
  sum2 = pnl_basis_eval_D2 (B, alpha, arg, 2, 2);
  sum3 = pnl_basis_eval_D (B, alpha, arg, 1);

  LET(D,0)=sum0;
  LET(D,1)=sum1;
  LET(D,2)=sum2;
  LET(D,3)=sum3;
}

static void pnl_basis_eval_test ()
{
  PnlMat *X;
  PnlVect *V,*x,*t,*D, *alpha;
  PnlBasis *basis;
  int j,m,n;
  double t0,x0;

  printf ("\n\n** Differentation of the regression **\n\n");
  m=19;//nombre de polynomes
  n=50;
  D=pnl_vect_create(4);
  x=pnl_vect_create(n);
  t=pnl_vect_create(n);
  t0=0.5;
  x0=2.5;
  //on tire aléatoirement les points auxquels on va évaluer la fonction à
  //retrouver
  pnl_vect_rand_uni(x,n,-5,4,7);
  pnl_vect_rand_uni(t,n,0,1,7);
  basis = pnl_basis_create (HERMITIAN, m, 2);
  alpha = pnl_vect_create (m);
  X = pnl_mat_create (n, 2);
  for(j=0;j<n;j++)
    {
      MLET (X, j, 0) = GET(t,j);
      MLET (X, j, 1) = GET(x,j);
    }
  V=pnl_vect_create(n);
  //création du vecteur contenant les valeurs de la fonction à retrouver aux
  //points t(j),x(j)
  for(j=0;j<n;j++)
    {
      LET(V,j)=fonction_a_retrouver(GET(t,j),GET(x,j));
    }
  pnl_basis_fit_ls (basis, alpha, X, V);
  /* calcul des valeurs approchées des dérivées de la fonction à retrouver (1
     en temps et 2 en espace) */
  derive_approx_fonction(basis, D, alpha,t0,x0);
  printf("valeur exacte de la fonction : %f\n",
         fonction_a_retrouver(t0,x0));
  printf("valeur approchee de la fonction :%f \n\n",
         pnl_vect_get(D,0));

  printf("valeur exacte de la derivee de la fonction : %f\n",
         derive_x_fonction_a_retrouver(t0,x0));
  printf("valeur approchee de la derivee de la fonction :%f \n\n",
         derive_x_approx_fonction(basis, alpha, t0, x0));

  printf("valeur exacte de la derivee seconde en espace fonction : %f\n",
         derive_xx_fonction_a_retrouver(t0,x0));
  printf("valeur approchee de la derivee seconde en espace de la fonction  :%f \n\n",
         pnl_vect_get(D,2));

  printf("valeur exacte de la derivee en temps fonction :%f\n",
         derive_t_fonction_a_retrouver(t0,x0));
  printf("valeur approchee de la derivee en temps de la fonction  :%f \n",
         pnl_vect_get(D,3));
  pnl_vect_free(&x);
  pnl_vect_free(&t);
  pnl_vect_free(&V);
  pnl_vect_free(&D);
  pnl_mat_free(&X);
}


void basis_test ()
{
  exp_regression2 ();
  regression_multid ();
  pnl_basis_eval_test ();
}
