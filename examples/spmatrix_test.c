
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

#include "pnl/pnl_sp_matrix.h"
#include "pnl/pnl_random.h"
#include "tests_utils.h"

static PnlSpMat* create_random_sp (int m, int n, PnlRng *rng) 
{
  int it;
  /* Sparsity rate between 20% and 60% */
  int nz = floor (((pnl_rng_uni(rng) * 0.4) + 0.2) * m * n);
  PnlSpMat *M = pnl_sp_mat_create (m, n, nz);

  for ( it=0 ; it<nz ; it++ )
    {
      int i = (int) floor (pnl_rng_uni (rng) * m);
      int j = (int) floor (pnl_rng_uni (rng) * n);
      double x = pnl_rng_normal (rng);
      pnl_sp_mat_set (M, i, j, x);
    }
  return M;
}

static void sp_get_set_test ()
{
  PnlSpMat *Sp;
  PnlMat *M;
  int i, j, nb_fail;
  PnlRng * rng = pnl_rng_create(PNL_RNG_MERSENNE);
  int m = 9, n = 11;
  double abserr = 1E-18;
  pnl_rng_sseed (rng, 0);
  Sp = create_random_sp (m, n, rng);
  M = pnl_mat_create_from_sp_mat (Sp);

  nb_fail = 0;
  for ( i=0 ; i<m ; i++ )
    {
      for ( j=0 ; j<n ; j++ )
        {
          nb_fail += fabs(pnl_sp_mat_get (Sp, i, j) -  MGET(M, i, j)) > abserr;
        }
    }
  if ( nb_fail == 0 )
    pnl_test_set_ok ("sp_mat_get_set");
  else
    pnl_test_set_fail ("sp_mat_get_set", nb_fail, 0);
  pnl_sp_mat_free (&Sp);
  pnl_mat_free (&M);
  pnl_rng_free (&rng);
}

static void sp_create_test ()
{
  PnlSpMat *Sp, *Sp2;
  PnlMat *M;
  PnlRng * rng = pnl_rng_create(PNL_RNG_MERSENNE);
  int m = 9, n = 11;
  pnl_rng_sseed (rng, 0);
  Sp = create_random_sp (m, n, rng);
  M = pnl_mat_create_from_sp_mat (Sp);
  Sp2 = pnl_sp_mat_create_from_mat (M);
  if ( pnl_sp_mat_eq (Sp2, Sp) == TRUE )
    pnl_test_set_ok ("sp_mat_create");
  else
    pnl_test_set_fail ("sp_mat_create", 0, 0);
  pnl_sp_mat_free (&Sp);
  pnl_sp_mat_free (&Sp2);
  pnl_mat_free (&M);
  pnl_rng_free (&rng);
}

static void sp_clone_test ()
{
  PnlSpMat *Sp, *Sp2;
  PnlRng * rng = pnl_rng_create(PNL_RNG_MERSENNE);
  int m = 9, n = 11;
  pnl_rng_sseed (rng, 0);
  Sp = create_random_sp (m, n, rng);
  Sp2 = pnl_sp_mat_new ();
  pnl_sp_mat_clone (Sp2, Sp);
   
  if ( pnl_sp_mat_eq (Sp2, Sp) == TRUE )
    pnl_test_set_ok ("sp_mat_clone");
  else
    pnl_test_set_fail ("sp_mat_clone", 0, 0);
  pnl_sp_mat_free (&Sp);
  pnl_sp_mat_free (&Sp2);
  pnl_rng_free (&rng);
}

#define NEQ_ERR(a,b) (fabs(a-b) > err)
static void sp_scalar_ops ()
{
  PnlSpMat *Sp, *Sp2;
  PnlRng * rng = pnl_rng_create(PNL_RNG_MERSENNE);
  int status, i, m = 9, n = 11;
  double x = 3.6;
  double err = 1E-12;
  pnl_rng_sseed (rng, 0);
  Sp = create_random_sp (m, n, rng);
  Sp2 = pnl_sp_mat_new ();

  /* plus op */
  status = OK;
  pnl_sp_mat_clone (Sp2, Sp);
  pnl_sp_mat_plus_scalar (Sp2, x);
  for ( i=0 ; i<Sp->nz ; i++ ) 
    {
      if ( NEQ_ERR(Sp->array[i] + x, Sp2->array[i]) ) { status = FAIL; break; }
    }
  if ( status == FAIL ) 
    pnl_test_set_fail ("sp_mat_plus_scalar", 0, 0);
  else
    pnl_test_set_ok ("sp_mat_plus_scalar");

  /* minus op */
  status = OK;
  pnl_sp_mat_clone (Sp2, Sp);
  pnl_sp_mat_minus_scalar (Sp2, x);
  for ( i=0 ; i<Sp->nz ; i++ ) 
    {
      if ( NEQ_ERR(Sp->array[i] - x, Sp2->array[i]) ) { status = FAIL; break; }
    }
  if ( status == FAIL ) 
    pnl_test_set_fail ("sp_mat_minus_scalar", 0, 0);
  else
    pnl_test_set_ok ("sp_mat_minus_scalar");

  /* mult op */
  status = OK;
  pnl_sp_mat_clone (Sp2, Sp);
  pnl_sp_mat_mult_scalar (Sp2, x);
  for ( i=0 ; i<Sp->nz ; i++ ) 
    {
      if ( NEQ_ERR(Sp->array[i] * x, Sp2->array[i]) ) { status = FAIL; break; }
    }
  if ( status == FAIL ) 
    pnl_test_set_fail ("sp_mat_mult_scalar", 0, 0);
  else
    pnl_test_set_ok ("sp_mat_mult_scalar");

  /* div op */
  status = OK;
  pnl_sp_mat_clone (Sp2, Sp);
  pnl_sp_mat_div_scalar (Sp2, x);
  for ( i=0 ; i<Sp->nz ; i++ ) 
    {
      if ( NEQ_ERR(Sp->array[i] / x, Sp2->array[i]) ) { status = FAIL; break; }
    }
  if ( status == FAIL ) 
    pnl_test_set_fail ("sp_mat_div_scalar", 0, 0);
  else
    pnl_test_set_ok ("sp_mat_div_scalar");


  pnl_sp_mat_free (&Sp);
  pnl_sp_mat_free (&Sp2);
  pnl_rng_free (&rng);
}

static void sp_mat_mult_vect ()
{
  PnlSpMat *Sp;
  PnlMat *M;
  PnlVect *x, *y1, *y2;
  PnlRng * rng = pnl_rng_create(PNL_RNG_MERSENNE);
  int m = 9, n = 11;
  double abserr = 1E-12;
  double lambda = 0.5, b = 2.;
  pnl_rng_sseed (rng, 0);
  Sp = create_random_sp (m, n, rng);
  M = pnl_mat_create_from_sp_mat (Sp);
  x = pnl_vect_new ();
  y1 = pnl_vect_new ();
  y2 = pnl_vect_new ();
  pnl_vect_rng_normal (x, n, rng);

  pnl_sp_mat_mult_vect (y1, Sp, x);
  pnl_mat_mult_vect_inplace (y2, M, x);
  pnl_test_vect_eq_abs (y1, y2, abserr, "sp_mat_mult_vect", "");

  pnl_vect_rng_normal (y1, m, rng);
  pnl_vect_clone (y2, y1);
  pnl_sp_mat_lAxpby (lambda, Sp, x, b, y2);
  pnl_mat_lAxpby (lambda, M, x, b, y1);
  pnl_test_vect_eq_abs (y1, y2, abserr, "sp_mat_lAxpby", "");

  pnl_sp_mat_free (&Sp);
  pnl_mat_free (&M);
  pnl_vect_free (&x);
  pnl_vect_free (&y1);
  pnl_vect_free (&y2);
  pnl_rng_free (&rng);
}

static void sp_mat_del_row ()
{
  PnlSpMat *Sp;
  PnlMat *M, *Mdel;
  PnlRng * rng = pnl_rng_create(PNL_RNG_MERSENNE);
  int m = 9, n = 11;
  double abserr = 1E-12;
  pnl_rng_sseed (rng, 0);

  Sp = create_random_sp (m, n, rng);
  M = pnl_mat_create_from_sp_mat (Sp);
  pnl_mat_del_row (M, 0);
  pnl_sp_mat_del_row (Sp, 0);
  Mdel = pnl_mat_create_from_sp_mat (Sp);
  pnl_test_mat_eq (M, Mdel, abserr, "sp_mat_del_row -- row 0", "");
  pnl_sp_mat_free (&Sp);
  pnl_mat_free (&M);
  pnl_mat_free (&Mdel);

  Sp = create_random_sp (m, n, rng);
  M = pnl_mat_create_from_sp_mat (Sp);
  pnl_mat_del_row (M, m/2);
  pnl_sp_mat_del_row (Sp, m/2);
  Mdel = pnl_mat_create_from_sp_mat (Sp);
  pnl_test_mat_eq (M, Mdel, abserr, "sp_mat_del_row -- middle", "");
  pnl_sp_mat_free (&Sp);
  pnl_mat_free (&M);
  pnl_mat_free (&Mdel);

  Sp = create_random_sp (m, n, rng);
  M = pnl_mat_create_from_sp_mat (Sp);
  pnl_mat_del_row (M, m-1);
  pnl_sp_mat_del_row (Sp, m-1);
  Mdel = pnl_mat_create_from_sp_mat (Sp);
  pnl_test_mat_eq (M, Mdel, abserr, "sp_mat_del_row -- last row", "");
  pnl_sp_mat_free (&Sp);
  pnl_mat_free (&M);
  pnl_mat_free (&Mdel);

  pnl_rng_free (&rng);
}

static double set_zero_neg (double x) { return MAX(0., x); }

static void sp_mat_add_row ()
{
  PnlSpMat *Sp;
  PnlMat *M, *Madd;
  PnlVect *d;
  PnlRng * rng = pnl_rng_create(PNL_RNG_MERSENNE);
  int m = 9, n = 11;
  double abserr = 1E-12;
  pnl_rng_sseed (rng, 0);
  d = pnl_vect_new ();
  pnl_vect_rng_normal (d, n, rng);
  pnl_vect_map_inplace (d, set_zero_neg);


  Sp = create_random_sp (m, n, rng);
  M = pnl_mat_create_from_sp_mat (Sp);
  pnl_mat_add_row (M, 0, d);
  pnl_sp_mat_add_row (Sp, 0, d);
  Madd = pnl_mat_create_from_sp_mat (Sp);
  pnl_test_mat_eq (M, Madd, abserr, "sp_mat_add_row -- row 0", "");
  pnl_sp_mat_free (&Sp);
  pnl_mat_free (&M);
  pnl_mat_free (&Madd);

  Sp = create_random_sp (m, n, rng);
  M = pnl_mat_create_from_sp_mat (Sp);
  pnl_mat_add_row (M, m/2, d);
  pnl_sp_mat_add_row (Sp, m/2, d);
  Madd = pnl_mat_create_from_sp_mat (Sp);
  pnl_test_mat_eq (M, Madd, abserr, "sp_mat_add_row -- middle", "");
  pnl_sp_mat_free (&Sp);
  pnl_mat_free (&M);
  pnl_mat_free (&Madd);

  Sp = create_random_sp (m, n, rng);
  M = pnl_mat_create_from_sp_mat (Sp);
  pnl_mat_add_row (M, m, d);
  pnl_sp_mat_add_row (Sp, m, d);
  Madd = pnl_mat_create_from_sp_mat (Sp);
  pnl_test_mat_eq (M, Madd, abserr, "sp_mat_add_row -- last row", "");
  pnl_sp_mat_free (&Sp);
  pnl_mat_free (&M);
  pnl_mat_free (&Madd);

  pnl_rng_free (&rng);
  pnl_vect_free (&d);
}


int main (int argc, char *argv[])
{
  pnl_test_init (argc, argv);
  sp_get_set_test ();
  sp_create_test ();
  sp_clone_test ();
  sp_scalar_ops ();
  sp_mat_mult_vect ();
  sp_mat_del_row ();
  sp_mat_add_row ();
  /* sp_add_test (); */
  exit (pnl_test_finalize("Sparse matrices"));
}
