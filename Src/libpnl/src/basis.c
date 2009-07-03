
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pnl_basis.h"
#include "pnl_vector.h"
#include "pnl_matrix.h"
#include "pnl_mathtools.h"

/*
 * maximal dimensions of the basis
 */
#define DimBasisDefaultHD1 6
#define DimBasisDefaultHD2 21
#define DimBasisDefaultHD3 20
#define DimBasisDefaultHD4 21
#define DimBasisDefaultHD5 20
#define DimBasisDefaultHD6 19
#define DimBasisDefaultHD7 20
#define DimBasisDefaultHD8 23
#define DimBasisDefaultHD9 26
#define DimBasisDefaultHD10 29
#define DimBasisDefaultTD1 6
#define DimBasisDefaultTD2 21
#define DimBasisDefaultTD3 20
#define DimBasisDefaultTD4 21
#define DimBasisDefaultTD5 20
#define DimBasisDefaultTD6 19
#define DimBasisDefaultTD7 20
#define DimBasisDefaultTD8 23
#define DimBasisDefaultTD9 26
#define DimBasisDefaultTD10 29
#define DimBasisDefaultD1 20
#define DimBasisDefaultD2 21
#define DimBasisDefaultD3 20
#define DimBasisDefaultD4 21
#define DimBasisDefaultD5 20
#define DimBasisDefaultD6 19
#define DimBasisDefaultD7 20
#define DimBasisDefaultD8 23
#define DimBasisDefaultD9 26
#define DimBasisDefaultD10 29



static double HermiteD1(double *x, int ind);
static double HermiteD2(double *x, int ind); 
static double HermiteD3(double *x, int ind);
static double HermiteD4(double *x, int ind);
static double HermiteD5(double *x, int ind);
static double HermiteD6(double *x, int ind);
static double HermiteD7(double *x, int ind); 
static double HermiteD8(double *x, int ind);
static double HermiteD9(double *x, int ind); 
static double HermiteD10(double *x, int ind);
static double CanoniqueD1(double *x, int ind);
static double CanoniqueD2(double *x, int ind); 
static double CanoniqueD3(double *x, int ind);
static double CanoniqueD4(double *x, int ind);
static double CanoniqueD5(double *x, int ind);
static double CanoniqueD6(double *x, int ind);
static double CanoniqueD7(double *x, int ind); 
static double CanoniqueD8(double *x, int ind);
static double CanoniqueD9(double *x, int ind); 
static double CanoniqueD10(double *x, int ind);
static double TchebychevD1(double *x, int ind);
static double TchebychevD2(double *x, int ind); 
static double TchebychevD3(double *x, int ind);
static double TchebychevD4(double *x, int ind);
static double TchebychevD5(double *x, int ind);
static double TchebychevD6(double *x, int ind);
static double TchebychevD7(double *x, int ind); 
static double TchebychevD8(double *x, int ind);
static double TchebychevD9(double *x, int ind); 
static double TchebychevD10(double *x, int ind);


static reg_basis Bases_tab[]=
  {
    {"Hermite", 1,  DimBasisDefaultHD1,  &HermiteD1 },
    {"Hermite", 2,  DimBasisDefaultHD2,  &HermiteD2 },
    {"Hermite", 3,  DimBasisDefaultHD3,  &HermiteD3 },
    {"Hermite", 4,  DimBasisDefaultHD4,  &HermiteD4 },
    {"Hermite", 5,  DimBasisDefaultHD5,  &HermiteD5 },
    {"Hermite", 6,  DimBasisDefaultHD6,  &HermiteD6 },
    {"Hermite", 7,  DimBasisDefaultHD7,  &HermiteD7 },
    {"Hermite", 8,  DimBasisDefaultHD8,  &HermiteD8 },
    {"Hermite", 9,  DimBasisDefaultHD9,  &HermiteD9 },
    {"Hermite", 10,  DimBasisDefaultHD10,  &HermiteD10 },
    {"Canonical", 1,  DimBasisDefaultD1,  &CanoniqueD1 },
    {"Canonical", 2,  DimBasisDefaultD2,  &CanoniqueD2 },
    {"Canonical", 3,  DimBasisDefaultD3,  &CanoniqueD3 },
    {"Canonical", 4,  DimBasisDefaultD4,  &CanoniqueD4 },
    {"Canonical", 5,  DimBasisDefaultD5,  &CanoniqueD5 },
    {"Canonical", 6,  DimBasisDefaultD6,  &CanoniqueD6 },
    {"Canonical", 7,  DimBasisDefaultD7,  &CanoniqueD7 },
    {"Canonical", 8,  DimBasisDefaultD8,  &CanoniqueD8 },
    {"Canonical", 9,  DimBasisDefaultD9,  &CanoniqueD9 },
    {"Canonical", 10,  DimBasisDefaultD10,  &CanoniqueD10 },
    {"Tchebychev", 1,  DimBasisDefaultD1,  &TchebychevD1 },
    {"Tchebychev", 2,  DimBasisDefaultD2,  &TchebychevD2 },
    {"Tchebychev", 3,  DimBasisDefaultD3,  &TchebychevD3 },
    {"Tchebychev", 4,  DimBasisDefaultD4,  &TchebychevD4 },
    {"Tchebychev", 5,  DimBasisDefaultD5,  &TchebychevD5 },
    {"Tchebychev", 6,  DimBasisDefaultD6,  &TchebychevD6 },
    {"Tchebychev", 7,  DimBasisDefaultD7,  &TchebychevD7 },
    {"Tchebychev", 8,  DimBasisDefaultD8,  &TchebychevD8 },
    {"Tchebychev", 9,  DimBasisDefaultD9,  &TchebychevD9 },
    {"Tchebychev", 10,  DimBasisDefaultD10,  &TchebychevD10 },
    {NULL, NULLINT, NULLINT, NULL},
  };

enum_member RegBasis [] =
  {
    { "Canonical", CANONICAL},
    { "Hermite", HERMITIAN},
    { "Tchebychev", TCHEBYCHEV},
    { NULL, NULLINT},
  };

DEFINE_ENUM( PnlBases, RegBasis)

/**
 * identify_basis returns a pointer of function
 * double fake (double *, int)
 *
 * @param index the index of the family to be used
 * @param nb_func the maximum number of functions which may be used
 * @param space_dim the size of the space in which the basis functions are
 * defined
 * @return a PnlBasis
 */
PnlBasis pnl_init_basis ( int index, int nb_func, int space_dim)
{
  enum_member *e = RegBasis;
  reg_basis *b = Bases_tab;

  while (e->label != NULL && e->key != index) { e++; }
  if (e->label == NULL )
    {
      printf ("No basis found : index exceeded\n"); abort();
    }

  while (b->label != NULL)
    {
      if (strcmp (e->label, b->label) == 0 &&
          space_dim == b->space_dim &&
          nb_func <= b->max_dim)
        return b->Compute;
      else
        b++;
    }
  printf ("No basis found : size exceeded \n"); abort();
  return NULL;
}

/* multidimensional bases are obtained as tensor products of
 * one-dimensional ones
 * the numbers inside the braces are the indices of the
 * underlying one-dimensional basis polynomials
 *
 * example : {2,1,4} = p2(x1)*p1(x2)*p4(x3) where p2,p1 and
 * p4 are the 2nd, 3rd, 4th elements of the underlying
 * one-dimensional basis
 *
 * the values come from Premia
 */
static int TensorBasisD2[DimBasisDefaultD2][2]=
  {{0,0},
                                        
   {1,0},{0,1},
                                        
   {1,1},{2,0},{0,2},
                                             
   {2,1},{1,2},{3,0},{0,3},
                                        
   {2,2},{1,3},{3,1},{4,0},{0,4},
                                        
   {1,4},{4,1},{3,2},{2,3},{5,0},{0,5}};
                                        
static int TensorBasisD3[DimBasisDefaultD3][3]=
  {{0,0,0},
   
   {1,0,0},{0,1,0},{0,0,1},
   
   {2,0,0},{0,2,0},{0,0,2},{1,1,0},{1,0,1},{0,1,1},
   
   {1,1,1},{2,1,0},{1,2,0},{0,1,2},{0,2,1},{1,0,2},{2,0,1},{3,0,0},{0,3,0},{0,0,3}};
                                        
static int TensorBasisD4[DimBasisDefaultD4][4]=
  {{0,0,0,0},
   
   {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1},
                                        
   {2,0,0,0},{0,2,0,0},{0,0,2,0},{0,0,0,2},{1,1,0,0},{0,1,1,0},{0,0,1,1},
                                        
   {1,1,1,0},{0,1,1,1},{1,0,1,1},{1,1,0,1},{3,0,0,0},{0,3,0,0},{0,0,3,0},{0,0,0,3},
                                        
   {1,1,1,1}};
                                        
static int TensorBasisD5[DimBasisDefaultD5][5]=
  {{0,0,0,0,0},
                                        
   {1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1},
                                        
   {2,0,0,0,0},{0,2,0,0,0},{0,0,2,0,0},{0,0,0,2,0},{0,0,0,0,2},{1,1,0,0,0},
   {1,1,0,0,0},{1,1,0,0,0},{0,1,1,0,0},{0,0,1,1,0},{0,0,0,1,1},
   
   {1,1,1,0,0},{0,1,1,1,0},{0,0,1,1,1}};
                                        
static int TensorBasisD6[DimBasisDefaultD6][6]=
  {{0,0,0,0,0,0},

   {1,0,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,1,0,0},{0,0,0,0,1,0},{0,0,0,0,0,1},
                                             
   {2,0,0,0,0,0},{0,2,0,0,0,0},{0,0,2,0,0,0},{0,0,0,2,0,0},{0,0,0,0,2,0},{0,0,0,0,0,2},

   {1,1,1,0,0,0},{0,1,1,1,0,0},{0,0,1,1,1,0},{0,0,0,1,1,1},

   {1,1,1,1,1,1}};
                                        
static int TensorBasisD7[DimBasisDefaultD7][7]=
  {{0,0,0,0,0,0,0},
                                        
   {1,0,0,0,0,0,0},{0,1,0,0,0,0,0},{0,0,1,0,0,0,0},{0,0,0,1,0,0,0},{0,0,0,0,1,0,0},
   {0,0,0,0,0,1,0},{0,0,0,0,0,0,1},
                                        
   {2,0,0,0,0,0,0},{0,2,0,0,0,0,0},{0,0,2,0,0,0,0},{0,0,0,2,0,0,0},{0,0,0,0,2,0,0},
   {0,0,0,0,0,2,0},{0,0,0,0,0,0,2},

   {1,1,1,1,0,0,0},{0,1,1,1,1,0,0},{0,0,1,1,1,1,0},{0,0,0,1,1,1,1},

   {1,1,1,1,1,1,1}};

                                        
static int TensorBasisD8[DimBasisDefaultD8][8]=
  {{0,0,0,0,0,0,0,0},

   {1,0,0,0,0,0,0,0},{0,1,0,0,0,0,0,0},{0,0,1,0,0,0,0,0},{0,0,0,1,0,0,0,0},
   {0,0,0,0,1,0,0,0},{0,0,0,0,0,1,0,0},{0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,1},
                                             
   {2,0,0,0,0,0,0,0},{0,2,0,0,0,0,0,0},{0,0,2,0,0,0,0,0},{0,0,0,2,0,0,0,0},
   {0,0,0,0,2,0,0,0},{0,0,0,0,0,2,0,0},{0,0,0,0,0,0,2,0},{0,0,0,0,0,0,0,2},

   {1,1,1,1,0,0,0,0},{0,1,1,1,1,0,0,0},{0,0,1,1,1,1,0,0},{0,0,0,1,1,1,1,0},
   {0,0,0,0,1,1,1,1},

   {1,1,1,1,1,1,1,1}};

                                        
static int TensorBasisD9[DimBasisDefaultD9][9]=
  {{0,0,0,0,0,0,0,0,0},
                                             
   {1,0,0,0,0,0,0,0,0},{0,1,0,0,0,0,0,0,0},{0,0,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0},
   {0,0,0,0,1,0,0,0,0},{0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,1,0,0},{0,0,0,0,0,0,0,1,0},
   {0,0,0,0,0,0,0,0,1},
                                        
   {2,0,0,0,0,0,0,0,0},{0,2,0,0,0,0,0,0,0},{0,0,2,0,0,0,0,0,0},{0,0,0,2,0,0,0,0,0},
   {0,0,0,0,2,0,0,0,0},{0,0,0,0,0,2,0,0,0},{0,0,0,0,0,0,2,0,0},{0,0,0,0,0,0,0,2,0},
   {0,0,0,0,0,0,0,0,2},

   {1,1,1,1,0,0,0,0,0},{0,1,1,1,1,0,0,0,0},{0,0,1,1,1,1,0,0,0},{0,0,0,1,1,1,1,0,0},
   {0,0,0,0,1,1,1,1,0},{0,0,0,0,0,1,1,1,1},
                                             
   {1,1,1,1,1,1,1,1,1}};

                                        
static int TensorBasisD10[DimBasisDefaultD10][10]=
  {{0,0,0,0,0,0,0,0,0,0},

   {1,0,0,0,0,0,0,0,0,0},{0,1,0,0,0,0,0,0,0,0},{0,0,1,0,0,0,0,0,0,0},
   {0,0,0,1,0,0,0,0,0,0},{0,0,0,0,1,0,0,0,0,0},{0,0,0,0,0,1,0,0,0,0},
   {0,0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,1,0,0},{0,0,0,0,0,0,0,0,1,0},
   {0,0,0,0,0,0,0,0,0,1},
                                               
   {2,0,0,0,0,0,0,0,0,0},{0,2,0,0,0,0,0,0,0,0},{0,0,2,0,0,0,0,0,0,0},
   {0,0,0,2,0,0,0,0,0,0},{0,0,0,0,2,0,0,0,0,0},{0,0,0,0,0,2,0,0,0,0},
   {0,0,0,0,0,0,2,0,0,0},{0,0,0,0,0,0,0,2,0,0},{0,0,0,0,0,0,0,0,2,0},
   {0,0,0,0,0,0,0,0,0,2},

   {1,1,1,1,0,0,0,0,0,0},{0,1,1,1,1,0,0,0,0,0},{0,0,1,1,1,1,0,0,0,0},
   {0,0,0,1,1,1,1,0,0,0},{0,0,0,0,1,1,1,1,0,0},{0,0,0,0,0,1,1,1,1,0},
   {0,0,0,0,0,0,1,1,1,1},
                                               
   {1,1,1,1,1,1,1,1,1}};


#define DEFINE_BASIS(name, dim)                     \
  static double name##D##dim(double *x, int ind)    \
  {                                                 \
    int i;                                          \
    double aux = 1;                                 \
    for (i=0;i<dim;i++){                            \
      aux*=name##D1(x+i,TensorBasisD##dim[ind][i]); \
    }                                               \
    return aux;                                     \
  }                                                       


/*
 * Canonical basis, dimension=1..10
 */
static double CanoniqueD1(double *x, int ind)
{
  int i;
  double aux=1;
  for (i=0;i<ind;i++){ aux*=(*x); }
  return aux;
}

DEFINE_BASIS(Canonique, 2)
DEFINE_BASIS(Canonique, 3)
DEFINE_BASIS(Canonique, 4)
DEFINE_BASIS(Canonique, 5)
DEFINE_BASIS(Canonique, 6)
DEFINE_BASIS(Canonique, 7)
DEFINE_BASIS(Canonique, 8)
DEFINE_BASIS(Canonique, 9)
DEFINE_BASIS(Canonique, 10)



/*
 * Hermite basis, dimension=1..10
 */
static double HermiteD1(double *x, int ind)
{
  switch (ind){
  case 0 : return 1;
  case 1 : return 1.414213562*(*x);
  case 2 : return 1.414213562*(*x)*(*x)-0.707106781;
  case 3 : return (1.154700538*(*x)*(*x)-1.732050808)*(*x);
  case 4 : return (0.816496581*(*x)*(*x)-2.449489743)*(*x)*(*x)+0.612372436;
  case 5 : return ((0.516397779*(*x)*(*x)-2.581988897)*(*x)*(*x)+1.936491673)*(*x);

  default : return 1;
  }
}

DEFINE_BASIS( Hermite, 2)
DEFINE_BASIS( Hermite, 3)
DEFINE_BASIS( Hermite, 4)
DEFINE_BASIS( Hermite, 5)
DEFINE_BASIS( Hermite, 6)
DEFINE_BASIS( Hermite, 7)
DEFINE_BASIS( Hermite, 8)
DEFINE_BASIS( Hermite, 9)
DEFINE_BASIS( Hermite, 10)


/*
 * Tchebychev basis, dimension=1..10
 */
static double TchebychevD1(double *x, int ind)
{
  double val = *x;
  double val2, val3, val4;
  switch (ind)
    {
    case 0 : return 1.;
    case 1 : return val;
    case 2 : return 2. * val * val - 1.;
    case 3 : return (4. * val * val - 3.) * val;
    case 4 :
      val2 = val * val;
      return 8. * val2 * val2 - 8. * val2 + 1.;
    case 5 :
      val2 = val * val; val3 = val2 * val;
      return 16. * val3 * val2 - 20. * val3 + 5.* val;
    case 6 :
      val2 = val * val; val4 = val2 * val2;
      return 32. * val4 * val2 - 48. * val4 + 18. * val2 - 1;
    default : return 1.;
    }
}

DEFINE_BASIS( Tchebychev, 2)
DEFINE_BASIS( Tchebychev, 3)
DEFINE_BASIS( Tchebychev, 4)
DEFINE_BASIS( Tchebychev, 5)
DEFINE_BASIS( Tchebychev, 6)
DEFINE_BASIS( Tchebychev, 7)
DEFINE_BASIS( Tchebychev, 8)
DEFINE_BASIS( Tchebychev, 9)
DEFINE_BASIS( Tchebychev, 10)


/**
 * Finds the best approximation of function defined by f(x(i,:)) = y(i)
 * 
 * @param f a PnlBasis
 * @param x the matrix of points at which we know the value of the function. One line
 * of the matrix is the vector of the coordinates of one point
 * @param y the values of the function f at the points defined by x
 * @param nb_func the number of elements to regress upon
 * @param coef contains on exit the coefficients of the regression
 *
 * @return OK or FAIL
 */
int pnl_fit_least_squares (PnlVect *coef, PnlMat *x, PnlVect *y, PnlBasis *f, int nb_func)
{
  int N, i, k;
  double b_k;
  PnlMat *A;
  PnlVect *phi_k;
  
  N = y->size;
  pnl_vect_resize_from_double (coef, nb_func, 0.);
  phi_k = pnl_vect_create_from_double (nb_func, 0.);
  A = pnl_mat_create_from_double (nb_func, nb_func, 0.);

  /* construct A and b*/
  for (i=0; i<N; i++)
    {
      for (k=0; k<nb_func; k++)
        {
          double tmp = (*f)(pnl_mat_lget(x, i, 0), k);
          b_k =  pnl_vect_get(coef, k);
          b_k += tmp * pnl_vect_get (y, i);
          pnl_vect_set (coef, k, b_k);
          pnl_vect_set (phi_k, k, tmp);
        }
      /* A += phi_k' * phi_k */
      pnl_mat_dger(1., phi_k, phi_k, A);
    }

  /* Solve A x = b, with A >0 symmetric */
  pnl_mat_chol (A);
  pnl_mat_chol_syslin_inplace (A, coef);

  pnl_vect_free (&phi_k);
  pnl_mat_free (&A);

  return OK;
}

/**
 * Evaluates a linear combination of basis functions at x
 *
 * @param coef a vector typically computed by pnl_fit_least_squares
 * @param x the coordinates of the point at which to evaluate the function
 * @param f a PnlBasis
 *
 * @return sum (coef .* f(x))
 */
double pnl_basis_eval (PnlVect *coef, double *x, PnlBasis *f)
{
  int i;
  double y;

  y = 0;
  for (i=0; i<coef->size; i++)
    {
      y += pnl_vect_get (coef, i) * (*f)(x, i);
    }
  return y;
}


/*
 * these undef must be at the end of the file
 */
   
#undef DimBasisDefaultHD1
#undef DimBasisDefaultHD2
#undef DimBasisDefaultHD3
#undef DimBasisDefaultHD4
#undef DimBasisDefaultHD5
#undef DimBasisDefaultHD6
#undef DimBasisDefaultHD7
#undef DimBasisDefaultHD8
#undef DimBasisDefaultHD9
#undef DimBasisDefaultHD10
#undef DimBasisDefaultTD1
#undef DimBasisDefaultTD2
#undef DimBasisDefaultTD3
#undef DimBasisDefaultTD4
#undef DimBasisDefaultTD5
#undef DimBasisDefaultTD6
#undef DimBasisDefaultTD7
#undef DimBasisDefaultTD8
#undef DimBasisDefaultTD9
#undef DimBasisDefaultTD10
#undef DimBasisDefaultD1
#undef DimBasisDefaultD2
#undef DimBasisDefaultD3
#undef DimBasisDefaultD4
#undef DimBasisDefaultD5
#undef DimBasisDefaultD6
#undef DimBasisDefaultD7
#undef DimBasisDefaultD8
#undef DimBasisDefaultD9
#undef DimBasisDefaultD10

