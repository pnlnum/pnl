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


#define DEFINE_BASIS_COMPUTES(name, dim)                            \
  static double name##D##dim(double *x, int ind)                    \
{                                                                   \
  int i;                                                            \
  double aux = 1;                                                   \
  for (i = 0 ; i < dim ; i++)                                       \
    {                                                               \
      aux *= name##D1 (x + i, TensorBasisD##dim[ind][i]);           \
    }                                                               \
  return aux;                                                       \
}                                                                   \
                                                                    \
static double D##name##D##dim(double *x, int ind, int k)            \
{                                                                   \
  int i;                                                            \
  double aux = 1;                                                   \
  for ( i = 0 ; i < dim ; i++ )                                     \
    {                                                               \
      if ( i == k-1 )                                               \
      aux *= _D##name##D1 (x + i, TensorBasisD##dim[ind][i]);       \
      else                                                          \
      aux *= name##D1 (x + i, TensorBasisD##dim[ind][i]);           \
    }                                                               \
  return aux;                                                       \
}                                                                   \
static double DD##name##D##dim(double *x, int ind, int k1, int k2)  \
{                                                                   \
  int i;                                                            \
  double aux = 1;                                                   \
  if (k1 == k2)                                                     \
    {                                                               \
      for ( i = 0 ; i < dim ; i++ )                                 \
        {                                                           \
          if ( i == k1-1 )                                          \
          aux *= _DD##name##D1 (x + i, TensorBasisD##dim[ind][i]);  \
          else                                                      \
          aux *= name##D1 (x + i, TensorBasisD##dim[ind][i]);       \
        }                                                           \
    }                                                               \
  else                                                              \
    {                                                               \
      for ( i = 0 ; i < dim ; i++ )                                 \
        {                                                           \
          if ( i == k1-1  || i == k2-1)                             \
          aux *= _D##name##D1 (x + i, TensorBasisD##dim[ind][i]);   \
          else                                                      \
          aux *= name##D1 (x + i, TensorBasisD##dim[ind][i]);       \
        }                                                           \
                                                                    \
    }                                                               \
  return aux;                                                       \
}



/*
 * Canonical basis, dimension=1..10
 */

/**
 *  Canonical polynomials
 *  @ param x the address of a real number
 *  @ param ind the index of the polynom to be evaluated
 */
static double CanonicalD1(double *x, int ind)
{
  return pnl_pow_i (*x, ind);
}

/**
 *  First derivative of the Canonical polynomials
 *  @ param x the address of a real number
 *  @ param ind the index of the polynom whose first derivative is to be evaluated
 */
static double _DCanonicalD1(double *x, int ind)
{
  if (ind == 0) return 0.;
  return ind * pnl_pow_i (*x, ind - 1);
}

static double DCanonicalD1(double *x, int ind, int i)
{
  return _DCanonicalD1(x, ind);
}


/**
 *  Second derivative of the Canonical polynomials
 *  @ param x the address of a real number
 *  @ param ind the index of the polynom whose second derivative is to be evaluated
 */
static double _DDCanonicalD1(double *x, int ind)
{
  if (ind <= 1) return 0.;
  return ind * (ind - 1) * pnl_pow_i (*x, ind - 2);
}

static double DDCanonicalD1(double *x, int ind, int i, int j)
{
  return _DCanonicalD1(x, ind);
}

DEFINE_BASIS_COMPUTES(Canonical, 2);
DEFINE_BASIS_COMPUTES(Canonical, 3);
DEFINE_BASIS_COMPUTES(Canonical, 4);
DEFINE_BASIS_COMPUTES(Canonical, 5);
DEFINE_BASIS_COMPUTES(Canonical, 6);
DEFINE_BASIS_COMPUTES(Canonical, 7);
DEFINE_BASIS_COMPUTES(Canonical, 8);
DEFINE_BASIS_COMPUTES(Canonical, 9);
DEFINE_BASIS_COMPUTES(Canonical, 10);

/*
 * Hermite basis, dimension=1..10
 */
/**
 *  Hermite polynomials
 *  @ param x the address of a real number
 *  @ param ind the index of the polynom to be evaluated
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

/**
 *  First derivative of the Hermite polynomials
 *  @ param x the address of a real number
 *  @ param ind the index of the polynom whose derivative is to be evaluated
 */
static double _DHermiteD1(double *x, int ind)
{
  switch (ind){
      case 0 : return 0;
      case 1 : return 1.414213562;
      case 2 : return 2*1.414213562*(*x);
      case 3 : return 3*1.154700538*(*x)*(*x)-1.732050808;
      case 4 : return 4*0.816496581*(*x)*(*x)*(*x)-2*2.449489743*(*x);
      case 5 : return 5*0.516397779*(*x)*(*x)*(*x)*(*x)-3*2.581988897*(*x)*(*x)+1.936491673;

      default : return 1;
  }
}

static double DHermiteD1(double *x, int ind, int i)
{
  return _DHermiteD1(x, ind);
}
/**
 *  Second derivative of the Hermite polynomials
 *  @ param x the address of a real number
 *  @ param ind the index of the polynom whose second derivative is to be evaluated
 */
static double _DDHermiteD1(double *x, int ind)
{
  switch (ind){
      case 0 : return 0;
      case 1 : return 0;
      case 2 : return 2*1.414213562;
      case 3 : return 6*1.154700538*(*x);
      case 4 : return 12*0.816496581*(*x)*(*x)-2*2.449489743;
      case 5 : return 20*0.516397779*(*x)*(*x)*(*x)-6*2.581988897*(*x);

      default : return 1;
  }
}

static double DDHermiteD1(double *x, int ind, int i, int j)
{
  return _DDHermiteD1(x, ind);
}

DEFINE_BASIS_COMPUTES( Hermite, 2);
DEFINE_BASIS_COMPUTES( Hermite, 3);
DEFINE_BASIS_COMPUTES( Hermite, 4);
DEFINE_BASIS_COMPUTES( Hermite, 5);
DEFINE_BASIS_COMPUTES( Hermite, 6);
DEFINE_BASIS_COMPUTES( Hermite, 7);
DEFINE_BASIS_COMPUTES( Hermite, 8);
DEFINE_BASIS_COMPUTES( Hermite, 9);
DEFINE_BASIS_COMPUTES( Hermite, 10);


/*
 * Tchebychev basis, dimension=1..10
 */
/**
 *  Tchebytchev polynomials
 *  @ param x the address of a real number
 *  @ param ind the index of the polynom to be evaluated
 */
static double TchebychevD1(double *x, int ind)
{
  double val = *x;
  double val2, val3, val4;
  switch (ind)
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
        val2 = val * val; val3 = val2 * val;
        return 16. * val3 * val2 - 20. * val3 + 5.* val;
      case 6 :
        val2 = val * val; val4 = val2 * val2;
        return 32. * val4 * val2 - 48. * val4 + 18. * val2 - 1;
      case 7 :
        val2 = val * val; val3 = val2 * val; val4 = val2 * val2;
        return (64. * val4 - 112. * val2 + 56) * val3 - 7. * val;
      default : PNL_ERROR ("order exceeded", "TchebytchevD1");
    }
}

/**
 *  First derivative of the Tchebytchev polynomials
 *  @ param x the address of a real number
 *  @ param ind the index of the polynom whose first derivative is to be evaluated
 */
static double _DTchebychevD1(double *x, int ind)
{
  double val = *x;
  double val2, val3, val4;
  switch (ind)
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
        val2 = val * val; val4 = val2 * val2;
        return (192. * val4 - 192. * val2 + 36.) * val;
      case 7 :
        val2 = val * val; val3 = val2 * val; val4 = val2 * val2;
        return (448. * val4 - 560. * val2 + 168) * val2 - 7.;
      default : PNL_ERROR ("order exceeded", "DTchebytchevD1");
    }
}

static double DTchebychevD1(double *x, int ind, int i)
{
  return _DTchebychevD1(x, ind);
}
/**
 *  Second derivative of the Tchebytchev polynomials
 *  @ param x the address of a real number
 *  @ param ind the index of the polynom whose second derivative is to be evaluated
 */
static double _DDTchebychevD1(double *x, int ind)
{
  double val = *x;
  double val2, val3, val4;
  switch (ind)
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
        return 240. * val2 * val - 160. * val;
      case 6 :
        val2 = val * val; val4 = val2 * val2;
        return (960. * val4 - 576. * val2 + 36.);
      case 7 :
        val2 = val * val; val3 = val2 * val; val4 = val2 * val2;
        return (2688. * val4 - 2240. * val2 + 336) * val;
      default : PNL_ERROR ("order exceeded", "DDTchebytchevD1");
    }
}

static double DDTchebychevD1(double *x, int ind, int i, int j)
{
  return _DDTchebychevD1(x, ind);
}

DEFINE_BASIS_COMPUTES( Tchebychev, 2);
DEFINE_BASIS_COMPUTES( Tchebychev, 3);
DEFINE_BASIS_COMPUTES( Tchebychev, 4);
DEFINE_BASIS_COMPUTES( Tchebychev, 5);
DEFINE_BASIS_COMPUTES( Tchebychev, 6);
DEFINE_BASIS_COMPUTES( Tchebychev, 7);
DEFINE_BASIS_COMPUTES( Tchebychev, 8);
DEFINE_BASIS_COMPUTES( Tchebychev, 9);
DEFINE_BASIS_COMPUTES( Tchebychev, 10);


#define DEFINE_BASIS(name, d, maxdim) \
{ #name, d, maxdim, name##D##d, D##name##D##d, DD##name##D##d }

static PnlBasis Bases_tab[]=
{
  DEFINE_BASIS (Hermite, 1, DimBasisDefaultHD1),
  DEFINE_BASIS (Hermite, 2, DimBasisDefaultHD2),
  DEFINE_BASIS (Hermite, 3, DimBasisDefaultHD3),
  DEFINE_BASIS (Hermite, 4, DimBasisDefaultHD4),
  DEFINE_BASIS (Hermite, 5, DimBasisDefaultHD5),
  DEFINE_BASIS (Hermite, 6, DimBasisDefaultHD6),
  DEFINE_BASIS (Hermite, 7, DimBasisDefaultHD7),
  DEFINE_BASIS (Hermite, 8, DimBasisDefaultHD8),
  DEFINE_BASIS (Hermite, 9, DimBasisDefaultHD9),
  DEFINE_BASIS (Hermite, 10, DimBasisDefaultHD10),
  DEFINE_BASIS (Canonical, 1, DimBasisDefaultD1),
  DEFINE_BASIS (Canonical, 2, DimBasisDefaultD2),
  DEFINE_BASIS (Canonical, 3, DimBasisDefaultD3),
  DEFINE_BASIS (Canonical, 4, DimBasisDefaultD4),
  DEFINE_BASIS (Canonical, 5, DimBasisDefaultD5),
  DEFINE_BASIS (Canonical, 6, DimBasisDefaultD6),
  DEFINE_BASIS (Canonical, 7, DimBasisDefaultD7),
  DEFINE_BASIS (Canonical, 8, DimBasisDefaultD8),
  DEFINE_BASIS (Canonical, 9, DimBasisDefaultD9),
  DEFINE_BASIS (Canonical, 10, DimBasisDefaultD10),
  DEFINE_BASIS (Tchebychev, 1, DimBasisDefaultD1),
  DEFINE_BASIS (Tchebychev, 2, DimBasisDefaultD2),
  DEFINE_BASIS (Tchebychev, 3, DimBasisDefaultD3),
  DEFINE_BASIS (Tchebychev, 4, DimBasisDefaultD4),
  DEFINE_BASIS (Tchebychev, 5, DimBasisDefaultD5),
  DEFINE_BASIS (Tchebychev, 6, DimBasisDefaultD6),
  DEFINE_BASIS (Tchebychev, 7, DimBasisDefaultD7),
  DEFINE_BASIS (Tchebychev, 8, DimBasisDefaultD8),
  DEFINE_BASIS (Tchebychev, 9, DimBasisDefaultD9),
  DEFINE_BASIS (Tchebychev, 10, DimBasisDefaultD10),
  {NULL, NULLINT, NULLINT, NULL, NULL, NULL},
};

enum_member _reg_basis [] =
{
    { "Canonical", CANONICAL},
    { "Hermite", HERMITIAN},
    { "Tchebychev", TCHEBYCHEV},
    { NULL, NULLINT},
};

DEFINE_ENUM( PnlBases, _reg_basis);

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
PnlBasis*  pnl_basis_init ( int index, int nb_func, int space_dim)
{
  enum_member *e = _reg_basis;
  PnlBasis *b = Bases_tab;

  while (e->label != NULL && e->key != index) { e++; }
  if (e->label == NULL )
    {
      printf ("No basis found : index exceeded\n"); abort();
    }

  while (b->label != NULL)
    {
      if ( strcmp (e->label, b->label) == 0 &&
          space_dim == b->space_dim &&
          nb_func <= b->max_dim )
        return b;
      else
        b++;
    }
  printf ("No basis found : size exceeded \n"); abort();
  return NULL;
}


/**
 * Finds the best approximation of function defined by f(x(i,:)) = y(i)
 *
 * @param basis a PnlBasis
 * @param x the matrix of points at which we know the value of the function. One line
 * of the matrix is the vector of the coordinates of one point
 * @param y the values of the function f at the points defined by x
 * @param nb_func the number of elements to regress upon
 * @param coef contains on exit the coefficients of the regression
 *
 * @return OK or FAIL
 */
int pnl_fit_least_squares (PnlVect *coef, PnlMat *x, PnlVect *y, PnlBasis *basis, int nb_func)
{
  int N, i, k;
  double b_k;
  PnlMat *A;
  PnlVect *phi_k;

  N = y->size;
  pnl_vect_resize (coef, nb_func);
  pnl_vect_set_double (coef, 0.);
  phi_k = pnl_vect_create_from_double (nb_func, 0.);
  A = pnl_mat_create_from_double (nb_func, nb_func, 0.);

  /* construct A and b*/
  for (i=0; i<N; i++)
    {
      for (k=0; k<nb_func; k++)
        {
          double tmp = (basis->f)(pnl_mat_lget(x, i, 0), k);
          b_k =  pnl_vect_get(coef, k);
          b_k += tmp * pnl_vect_get (y, i);
          pnl_vect_set (coef, k, b_k);
          pnl_vect_set (phi_k, k, tmp);
        }
      /* A += phi_k' * phi_k */
      pnl_mat_dger(1., phi_k, phi_k, A);
    }

  /* Solve A x = b, with A >0 symmetric */
  /* pnl_mat_chol (A);
     pnl_mat_chol_syslin_inplace (A, coef); */
  /* Because A often comes from simulation, A is not >0. So we use a QR
     approach */
#ifdef HAVE_LAPACK
  pnl_mat_ls (A, coef);
#else
  pnl_mat_syslin_inplace (A, coef);
#endif
  
  pnl_vect_free (&phi_k);
  pnl_mat_free (&A);

  return OK;
}

/**
 * Evaluates a linear combination of basis functions at x
 *
 * @param coef a vector typically computed by pnl_fit_least_squares
 * @param x the coordinates of the point at which to evaluate the function
 * @param basis a PnlBasis
 *
 * @return sum (coef .* f(x))
 */
double pnl_basis_eval (PnlVect *coef, double *x, PnlBasis *basis)
{
  int i;
  double y;

  y = 0;
  for (i=0; i<coef->size; i++)
    {
      y += pnl_vect_get (coef, i) * (basis->f)(x, i);
    }
  return y;
}

/**
 * Evaluates the first derivative with respect to x[i] of a linear combination
 * of basis functions at x
 *
 * @param coef a vector typically computed by pnl_fit_least_squares
 * @param x the coordinates of the point at which to evaluate the function
 * @param basis a PnlBasis
 * @param i the index with respect to which the derivative is computed
 *
 * @return sum (coef .* D_i f(x))
 */
double pnl_basis_eval_D (PnlVect *coef, double *x, PnlBasis *basis, int i)
{
  int k;
  double y;

  y = 0;
  for (k=0; k<coef->size; k++)
    {
      y += pnl_vect_get (coef, k) * (basis->Df)(x, k, i);
    }
  return y;
}

/**
 * Evaluates the second derivative with respect to x[i] and x[j] of a linear
 * combination of basis functions at x
 *
 * @param coef a vector typically computed by pnl_fit_least_squares
 * @param x the coordinates of the point at which to evaluate the function
 * @param basis a PnlBasis
 * @param i the index with respect to which the derivative is computed
 * @param j the index with respect to which the derivative is computed
 *
 * @return sum (coef .* D2_{i,j} f(x))
 */
double pnl_basis_eval_D2 (PnlVect *coef, double *x, PnlBasis *basis, int i, int j)
{
  int k;
  double y;

  y = 0;
  for (k=0; k<coef->size; k++)
    {
      y += pnl_vect_get (coef, k) * (basis->D2f)(x, k, i, j);
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

