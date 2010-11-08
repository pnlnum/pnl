
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
#include <time.h>

#include "pnl/pnl_matrix.h"
#include "pnl/pnl_basis.h"
#include "pnl/pnl_random.h"

#define XMIN -4.
#define XMAX 4.
#define NEVAL 1E6

/*
 * Normal density in dimension d
 */
static double f (const PnlVect *x)
{
  int i, d;
  double y, xi;
  d = x->size;
  y = 0;
  for ( i=0 ; i<d ; i++ ) 
    {
      xi = PNL_GET(x,i);
      y += exp (- xi*xi / 2.);
    }
  return y / pow(M_2PI, d/2.);
}

static void create_data (PnlMat *X, PnlVect *y, int d, int n)
{
  int i;
  PnlRng *rng;
  rng = pnl_rng_create (PNL_RNG_MERSENNE);
  pnl_rng_sseed (rng, time(NULL));
  pnl_mat_rng_uni2 (X, n, d, XMIN, XMAX, rng);
  for ( i=0 ; i<n ; i++ )
    {
      const PnlVect Xi = pnl_vect_wrap_mat_row (X, i);
      PNL_LET(y, i) = f (&Xi);
    }
  pnl_rng_free (&rng);
}

void speed_basis_test ()
{
  PnlMat *X;
  PnlVect *y, *x, *xmin, *xmax;
  PnlVect *alpha, *alphar;
  PnlBasis *b, *br;
  double start, end;
  int i, n, d;
  PnlRng *rng;

  printf ("--> basis evaluation \n");
  rng = pnl_rng_create (PNL_RNG_MERSENNE);
  pnl_rng_sseed (rng, time(NULL));
  n = 1000;
  d = 5;

  alpha = pnl_vect_new ();
  alphar = pnl_vect_new ();
  x = pnl_vect_new ();
  xmin = pnl_vect_create_from_double (d, XMIN);
  xmax = pnl_vect_create_from_double (d, XMAX);
  X = pnl_mat_create (n, d);
  y = pnl_vect_create (n);

  b = pnl_basis_create_from_degree (PNL_BASIS_HERMITIAN, 3, d);
  br = pnl_basis_create_from_degree (PNL_BASIS_HERMITIAN, 3, d);

  pnl_basis_set_domain (br, xmin, xmax); 

  create_data (X, y, d, n);
  pnl_basis_fit_ls(b, alpha, X, y);
  pnl_basis_fit_ls(br, alphar, X, y);

  pnl_vect_rng_uni (x, d, XMIN, XMAX, rng);
  printf ("   Evalution using non reduced bases\n");
  start = clock();
  for ( i=0 ; i<NEVAL ; i++ )
    {
      pnl_basis_eval (b, alpha, x->array);
    }
  end = clock();
  printf ("      CPU time : %f\n", (end-start) / CLOCKS_PER_SEC);
  printf("   Evalution using reduced bases\n");
  start = clock();
  for ( i=0 ; i<NEVAL ; i++ )
    {
      pnl_basis_eval (br, alphar, x->array);
    }
  end = clock();
  printf ("      CPU time : %f\n", (end-start) / CLOCKS_PER_SEC);

  pnl_basis_free (&b);
  pnl_basis_free (&br);
  pnl_vect_free (&y);
  pnl_vect_free (&x);
  pnl_vect_free (&xmin);
  pnl_vect_free (&xmax);
  pnl_vect_free (&alpha);
  pnl_vect_free (&alphar);
  pnl_mat_free (&X);
  pnl_rng_free (&rng);

}
