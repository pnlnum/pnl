
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

#include "pnl/pnl_matrix.h"
#include "pnl/pnl_random.h"
#include "tests_utils.h"

#define Ndim 3

/** 
 * Fills in x with normal random numbers
 * 
 * @param x an array
 * @param n size of x
 * @param rng a PnlRng
 */
static void random_array (double *x, int n, PnlRng *rng)
{
  int i;
  for ( i=0 ; i<n ; i++ )
    {
      x[i] = pnl_rng_normal (rng);
    }
}

static void pnl_hmat_copy_test()
{
  PnlHmat *M1;
  PnlHmat *M2;
  PnlRng *rng;
  int dims[Ndim] = { 3, 4, 5};
  rng = pnl_rng_create (PNL_RNG_MERSENNE);
  pnl_rng_sseed (rng, 4172);
  M2 = pnl_hmat_create (Ndim, dims);
  random_array (M2->array, M2->mn, rng);
  M1 = pnl_hmat_copy(M2);
  pnl_test_hmat_eq_abs (M1, M2, 1E-18, "hmat_copy", "");
  pnl_hmat_free(&M1);
  pnl_hmat_free(&M2);
  pnl_rng_free (&rng);
}

static void pnl_hmat_clone_test()
{
  PnlHmat *M1;
  PnlHmat *M2;
  PnlRng *rng;
  int dims[Ndim] = { 3, 4, 5};
  rng = pnl_rng_create (PNL_RNG_MERSENNE);
  pnl_rng_sseed (rng, 4172);
  M2 = pnl_hmat_create (Ndim, dims);
  random_array (M2->array, M2->mn, rng);
  M1 = pnl_hmat_new();
  pnl_hmat_clone (M1, M2);
  pnl_test_hmat_eq_abs (M1, M2, 1E-18, "hmat_copy", "");
  pnl_hmat_free(&M1);
  pnl_hmat_free(&M2);
  pnl_rng_free (&rng);
}

/*
 * Test the linear bijection for storage
 */
static void pnl_hmat_set_test ()
{
  PnlHmat *M;
  int i, j, k;
  int t[Ndim];
  int dims[Ndim] = { 3, 4, 5};

  M = pnl_hmat_create (Ndim, dims);
  for ( i=0 ; i<dims[0] ; i++ )
    for ( j=0 ; j<dims[1] ; j++ )
      for ( k=0 ; k<dims[2] ; k++ )
        {
          t[0] = i; t[1] = j; t[2] = k;
          pnl_hmat_set (M, t, pnl_pow_i (2, i) * pnl_pow_i (3, j));
        }

  for ( i=0 ; i<dims[0] ; i++ )
    for ( j=0 ; j<dims[1] ; j++ )
      {
        double expected = pnl_pow_i (2, i) * pnl_pow_i (3, j);
        for ( k=0 ; k<dims[2] ; k++ )
          {
            double x = M->array[i*M->pdims[0]+j*M->pdims[1]+k];
            if ( x != expected )
              {
                pnl_test_set_fail ( "hmat_set", x, expected);
                goto J1;
              }
          }
      }
  pnl_test_set_ok ("hmat_set");
J1:
  pnl_hmat_free(&M);
}

/*
 * Test the linear bijection for storage
 */
static void pnl_hmat_get_test ()
{
  PnlHmat *M;
  int i, j, k;
  int t[Ndim];
  int dims[Ndim] = { 3, 4, 5};

  M = pnl_hmat_create (Ndim, dims);
  for ( i=0 ; i<dims[0] ; i++ )
    for ( j=0 ; j<dims[1] ; j++ )
      {
        double x = pnl_pow_i (2, i) * pnl_pow_i (3, j);
        for ( k=0 ; k<dims[2] ; k++ )
          {
            M->array[i*M->pdims[0]+j*M->pdims[1]+k] = x;
          }
      }

  for ( i=0 ; i<dims[0] ; i++ )
    for ( j=0 ; j<dims[1] ; j++ )
      {
        double expected = pnl_pow_i (2, i) * pnl_pow_i (3, j);
        for ( k=0 ; k<dims[2] ; k++ )
          {
            double x;
            t[0] = i; t[1] = j; t[2] = k;
            x =  pnl_hmat_get (M, t);
            if ( x != expected )
              {
                pnl_test_set_fail ( "hmat_get", x, expected);
                goto J1;
              }
          }
      }

  pnl_test_set_ok ("hmat_get");
J1:
  pnl_hmat_free(&M);
}

static void pnl_mat_wrap_hmat_test ()
{
  PnlHmat *H;
  PnlMat M;
  PnlRng *rng;
  int i, j;
  int dims[Ndim] = { 3, 4, 5};
  int t = 1;
  rng = pnl_rng_create (PNL_RNG_MERSENNE);
  pnl_rng_sseed (rng, 4172);
  H = pnl_hmat_create (Ndim, dims);
  random_array (H->array, H->mn, rng);

  M = pnl_mat_wrap_hmat (H, &t);
  for ( i=0 ; i<dims[1] ; i++ )
    for ( j=0 ; j<dims[2] ; j++ )
      {
        int tab[Ndim] = { t, i, j};
        double expected = pnl_hmat_get( H, tab);
        double x = MGET(&M, i, j);
        if ( x != expected )
          {
            pnl_test_set_fail ( "mat_wrap_hmat", x, expected);
            goto J1;
          }

      }
  pnl_test_set_ok ("mat_wrap_hmat");
J1:
  pnl_hmat_free (&H);
  pnl_rng_free (&rng);
}

static void pnl_vect_wrap_hmat_test ()
{
  PnlHmat *H;
  PnlVect V;
  PnlRng *rng;
  int i;
  int dims[Ndim] = { 3, 4, 5};
  int t[2] = {1, 2};
  rng = pnl_rng_create (PNL_RNG_MERSENNE);
  pnl_rng_sseed (rng, 4172);
  H = pnl_hmat_create (Ndim, dims);
  random_array (H->array, H->mn, rng);

  V = pnl_vect_wrap_hmat (H, t);
  for ( i=0 ; i<dims[2] ; i++ )
    {
      int tab[Ndim] = { t[0], t[1], i};
      double expected = pnl_hmat_get( H, tab);
      double x = GET(&V, i);
      if ( x != expected )
        {
          pnl_test_set_fail ( "vect_wrap_hmat", x, expected);
          goto J1;
        }

    }
  pnl_test_set_ok ("vect_wrap_hmat");
J1:
  pnl_hmat_free (&H);
  pnl_rng_free (&rng);
}



int main (int argc, char **argv)
{
  pnl_test_init (argc, argv);
  pnl_hmat_copy_test ();
  pnl_hmat_clone_test ();
  pnl_hmat_set_test ();
  pnl_hmat_get_test ();
  pnl_mat_wrap_hmat_test ();
  pnl_vect_wrap_hmat_test ();
  exit (pnl_test_finalize ("Hmatrix"));
}
