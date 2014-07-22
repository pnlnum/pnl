
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

#include "pnl/pnl_vector.h"
#include "pnl/pnl_random.h"
#include "tests_utils.h"

/* static double function_prod(double x, double y) {return x*y;} */

static void pnl_vect_create_from_ptr_test()
{
  PnlVect *v1, *v2;
  int size;
  double x[4]={1.0, 5.0, 3.0, 8.0};
  size=4;
  v1=pnl_vect_create_from_ptr(size,x);
  v2=pnl_vect_create_from_list(size,1.0, 5.0, 3.0, 8.0);
  pnl_test_vect_eq_abs (v1, v2, 1E-18, "vect_create_from_ptr / vect_create_from_list", "");
  pnl_vect_free(&v1);
  pnl_vect_free(&v2);

}

static void pnl_vect_copy_test()
{
  PnlVect *v1;
  PnlVect *v2;
  v2=pnl_vect_create_from_scalar(4,3.0);
  v1=pnl_vect_copy(v2);
  pnl_test_vect_eq_abs (v1, v2, 1E-18, "vect_copy", "");
  pnl_vect_free(&v1);
  pnl_vect_free(&v2);
}

static void pnl_vect_clone_test()
{
  PnlVect *v1;
  PnlVect *v2;
  v2=pnl_vect_create_from_scalar(4,3.0);
  v1=pnl_vect_create(0);
  pnl_vect_clone(v1,v2);
  pnl_test_vect_eq_abs (v1, v2, 1E-18, "vect_clone", "");
  pnl_vect_free(&v1);
  pnl_vect_free(&v2); 
}

static void pnl_vect_extract_test()
{
  double x[10]={0.0,1.0, 2.0, 3.0,4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  PnlVect v1;
  PnlVect *v2, *v3, *v4;
  v2=pnl_vect_create_from_ptr(10,x);
  v1=pnl_vect_wrap_subvect_with_last(v2,3,5);
  v3 = pnl_vect_create_from_ptr (3, x + 3);
  pnl_test_vect_eq_abs (&v1, v3, 1E-18, "wrap_subvect_with_last", "");
  v1=pnl_vect_wrap_subvect(v2,5,4);
  v4 = pnl_vect_create_from_ptr (4, x + 5);
  pnl_test_vect_eq_abs (&v1, v4, 1E-18, "wrap_subvect_with_last", "");
  pnl_vect_free(&v2); 
  pnl_vect_free(&v3); 
  pnl_vect_free(&v4); 
}

static int ispos (double *x) { return *x >= 0; }
static int islarger (double *t) { return t[0] >= t[1]; }

static void pnl_vect_subvect_test ()
{
  int i, gen=PNL_RNG_MERSENNE_RANDOM_SEED;
  PnlVect *v1, *v2, *v3;
  PnlVectInt *ind;
  pnl_rand_init (gen, 10, 1);
  v1 = pnl_vect_create (10);
  v3 = pnl_vect_create (10);
  pnl_vect_rand_normal(v1, 10, gen);
  pnl_vect_rand_normal(v3, 10, gen);
  ind = pnl_vect_int_create (10);

  pnl_vect_find (ind, "v", ispos, v1);
  v2 = pnl_vect_create_subvect_with_ind (v1, ind);
  for ( i=0 ; i<v2->size ; i++ )
    {
      if ( ! ispos (& v2->array[i] ) )
        {
          pnl_test_set_fail ("vect_find (ispos)", GET(v2, i), 0);
          goto J1;
        }
    }
  pnl_test_set_ok ("vect_find (ispos)");
J1:
  pnl_vect_free (&v2);

  pnl_vect_find (ind, "vv", islarger, v1, v3);
  v2 = pnl_vect_create_subvect_with_ind (v1, ind);
  for ( i=0 ; i<v2->size ; i++ )
    {
      if ( GET(v2, i) != GET(v1, PNL_GET(ind, i))  ||
           ! (GET(v1, PNL_GET(ind,i)) >= GET(v3, PNL_GET(ind,i))) )
        {
          pnl_test_set_fail ("vect_find (islarger)", GET(v2, i), 0);
          goto J2;
        }
    }
  pnl_test_set_ok ("vect_find (larger)");
J2:
  pnl_vect_free (&v2);
  pnl_vect_int_free (&ind);
  pnl_vect_free (&v1);
  pnl_vect_free (&v3);
}

static void pnl_vect_plus_vect_test()
{
  int i;
  PnlVect *v1;
  PnlVect *v2;
  double x[4]={1.0, 5.0, 3.0, 8.0};
  v1=pnl_vect_create_from_ptr(4,x);
  v2=pnl_vect_create_from_scalar(4,3.0);
  pnl_vect_plus_vect(v1,v2);
  for ( i=0; i<4; i++ )
    {
      if ( GET(v1,i) != x[i] + 3. )
        {
          pnl_test_set_fail ("vect_plus_vect", GET(v1,i), x[i] + 3.);
          goto J1;
        }
    }
  pnl_test_set_ok ("vect_plus_vect");
J1:
  pnl_vect_free(&v1);
  pnl_vect_free(&v2); 
}

static void pnl_vect_map_inplace_test()
{
  int i;
  PnlVect *v;
  v=pnl_vect_create_from_scalar(4,3.0);
  pnl_vect_map_inplace(v,exp);
  for ( i=0; i<4; i++ )
    {
      if ( GET(v,i) != exp(3.) )
        {
          pnl_test_set_fail ("vect_map_inplace", GET(v,i), exp(3.));
          goto J1;
        }
    }
  pnl_test_set_ok ("vect_map_inplace");
J1:
  pnl_vect_free(&v);
}

static void pnl_vect_plus_scalar_test()
{
  int i;
  PnlVect *v;
  v=pnl_vect_create_from_scalar(4,3.0);
  pnl_vect_plus_scalar(v,0.5);
  for ( i=0; i<4; i++ )
    {
      if ( GET(v,i) != 3.5 )
        {
          pnl_test_set_fail ("vect_plus_scalar", GET(v,i), 3.5);
          goto J1;
        }
    }
  pnl_test_set_ok ("vect_plus_scalar");
J1:
  pnl_vect_free(&v);
}

static void pnl_vect_mult_scalar_test()
{
  int i;
  PnlVect *v;
  v=pnl_vect_create_from_scalar(4,3.0);
  pnl_vect_mult_scalar(v,0.5);
  for ( i=0; i<4; i++ )
    {
      if ( GET(v,i) != 1.5 )
        {
          pnl_test_set_fail ("vect_mult_scalar", GET(v,i), 1.5);
          goto J1;
        }
    }
  pnl_test_set_ok ("vect_mult_scalar");
J1:
  pnl_vect_free(&v);
}

static void pnl_vect_inv_term_test()
{
  int i;
  PnlVect *v;
  v=pnl_vect_create_from_scalar(4,3.0);
  pnl_vect_inv_term(v);
  for ( i=0; i<4; i++ )
    {
      if ( GET(v,i) != 1. / 3. )
        {
          pnl_test_set_fail ("vect_mult_scalar", GET(v,i), 1 / 3.);
          goto J1;
        }
    }
  pnl_test_set_ok ("vect_inv_term");
J1:
  pnl_vect_free(&v);
}

static void pnl_vect_div_vect_term_test()
{
  int i;
  PnlVect *v1;
  PnlVect *v2;
  double x[4]={1.0, 5.0, 3.0, 8.0};
  double abserr = 1E-12;
  v1=pnl_vect_create_from_ptr(4,x);
  v2=pnl_vect_create_from_scalar(4,3.0);
  pnl_vect_div_vect_term(v1,v2);
  for ( i=0; i<4; i++ )
    {
      if ( pnl_cmp_eq_abs(GET(v1,i), x[i] / 3., abserr) )
        {
          pnl_test_set_fail ("vect_div_vect_term", GET(v1,i), x[i] / 3.);
          goto J1;
        }
    }
  pnl_test_set_ok ("vect_div_vect_term");
J1:
  pnl_vect_free(&v1);
  pnl_vect_free(&v2);
}

static void pnl_vect_mult_vect_term_test()
{
  int i;
  PnlVect *v1;
  PnlVect *v2;
  double x[4]={1.0, 5.0, 3.0, 8.0};
  v1=pnl_vect_create_from_ptr(4,x);
  v2=pnl_vect_create_from_scalar(4,3.0);
  pnl_vect_mult_vect_term(v1,v2);
  for ( i=0; i<4; i++ )
    {
      if ( GET(v1,i) != x[i] * 3. )
        {
          pnl_test_set_fail ("vect_mult_vect_term", GET(v1,i), x[i] * 3.);
          goto J1;
        }
    }
  pnl_test_set_ok ("vect_mult_vect_term");
J1:
  pnl_vect_free(&v1);
  pnl_vect_free(&v2);
}

static void  pnl_vect_sum_test()
{
  int i;
  PnlVect *v;
  PnlVect *cumsumv, *cumsumx;
  double sumx, sumv;
  double x[4]={1.0, 5.0, 3.0, 8.0};
  cumsumx = pnl_vect_create (4);
  v=pnl_vect_create_from_ptr(4,x);

  sumv =  pnl_vect_sum(v);
  cumsumv = pnl_vect_copy (v);
  pnl_vect_cumsum(cumsumv);
  sumx = 0;
  for ( i=0 ; i<4 ; i++ )
    {
      sumx += x[i];
      LET(cumsumx, i) = sumx;
    }
  pnl_test_eq_abs (sumv, sumx, 1E-18, "vect_sum", "");
  pnl_test_vect_eq_abs (cumsumv, cumsumx, 1E-18, "vect_cumsum", "");
  pnl_vect_free(&v); 
  pnl_vect_free(&cumsumv); 
  pnl_vect_free(&cumsumx); 
}

static void pnl_vect_map_test()
{
  int i;
  PnlVect *v1;
  PnlVect *v2;
  v1=pnl_vect_create(0);
  v2=pnl_vect_create_from_scalar(4,3.0);
  pnl_vect_map(v1,v2,exp);
  for ( i=0; i<4; i++ )
    {
      if ( GET(v1,i) != exp(GET(v2,i)) )
        {
          pnl_test_set_fail ("vect_map", GET(v1,i), exp(GET(v2,i)));
          goto J1;
        }
    }
  pnl_test_set_ok ("vect_map");
J1:
  pnl_vect_free(&v1);
  pnl_vect_free(&v2);
}

static void pnl_vect_scalar_prod_test()
{
  PnlVect *v1;
  PnlVect *v2;
  int i;
  double scalar, res;
  double x[4]={1.0, 5.0, 3.0, 8.0};
  v1=pnl_vect_create_from_ptr(4,x);
  v2=pnl_vect_create_from_scalar(4,3.0);
  scalar = pnl_vect_scalar_prod(v1,v2);
  res = 0.;
  for ( i=0 ; i<4 ; i++ )
    {
      res += GET(v1, i) * GET(v2, i);
    }
  pnl_test_eq_abs (scalar, res, 1E-18, "vect_scalar_prod", "");
  pnl_vect_free(&v1);
  pnl_vect_free(&v2);
}

static void pnl_vect_prod_test()
{
  PnlVect *v, *cumprodv, *cumprodx, *vcopy;
  double prodx, prodv;
  int i;
  double x[4]={1.0, 5.0, 3.0, 8.0};
  v = pnl_vect_create_from_ptr (4, x);
  cumprodx = pnl_vect_copy (v);
  cumprodv = pnl_vect_copy (v);
  vcopy = pnl_vect_copy (v);
  prodv = pnl_vect_prod(v);

  pnl_vect_cumprod(cumprodv);
  prodx = 1;
  for ( i=0 ; i<4 ; i++ )
    {
      prodx *= x[i];
      LET(cumprodx, i) = prodx;
    }
  pnl_test_eq_abs (prodv, prodx, 1E-18, "vect_prod", "");
  pnl_test_vect_eq_abs (cumprodv, cumprodx, 1E-18, "vect_cumprod", "");
  pnl_vect_free(&v);
  pnl_vect_free(&vcopy);
  pnl_vect_free(&cumprodv);
  pnl_vect_free(&cumprodx);
}

static void pnl_vect_minmax_test()
{
  PnlVect *v;
  double min, max, vmin, vmax, abserr;
  int imin, imax;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  abserr = 1E-18;
  pnl_rand_init (gen, 1, 5);
  v = pnl_vect_create (0);
  pnl_vect_rand_normal (v, 5, gen);
  pnl_vect_minmax (&min, &max, v);
  vmin = pnl_vect_min (v);
  vmax = pnl_vect_max (v);
  pnl_test_eq_abs (min, vmin, abserr, "vect_min", "");
  pnl_test_eq_abs (max, vmax, abserr, "vect_max", "");

  pnl_vect_min_index (&min, &imin, v);
  pnl_test_eq_abs (GET(v,imin), vmin, abserr, "vect_min_index (imin)", "");
  pnl_test_eq_abs (min, vmin, abserr, "vect_min_index (min)", "");

  pnl_vect_max_index (&max, &imax, v);
  pnl_test_eq_abs (GET(v,imax), vmax, abserr, "vect_max_index (imax)", "");
  pnl_test_eq_abs (max, vmax, abserr, "vect_max_index (max)", "");

  pnl_vect_minmax_index (&min, &max, &imin, &imax, v);
  pnl_test_eq_abs (GET(v,imin), vmin, abserr, "vect_minmax_index (imin)", "");
  pnl_test_eq_abs (min, vmin, abserr, "vect_minmax_index (min)", "");
  pnl_test_eq_abs (GET(v,imax), vmax, abserr, "vect_minmax_index (imax)", "");
  pnl_test_eq_abs (max, vmax, abserr, "vect_minmax_index (max)", "");

  pnl_vect_free(&v);
}


static void pnl_vect_qsort_test ()
{
  int i;
  PnlVect *v = pnl_vect_create (0);
  PnlVect *vclone = pnl_vect_create (0);
  PnlVectInt *t = pnl_vect_int_create (0);
  pnl_rand_init (7, 1, 1);
  pnl_vect_rand_uni (v, 20, 0., 1., 7);
  pnl_vect_qsort (v, 'i');

  for ( i=1 ; i<v->size ; i++ )
    {
      if ( GET (v,i-1) > GET(v,i) ) 
        {
          pnl_test_set_fail ("vect_qsort", 0, 0);
          goto J1;
        }
    }
  pnl_test_set_ok ("vect_qsort");

J1:
  pnl_vect_rand_uni (v, 20, 0., 1., 7);
  pnl_vect_clone (vclone, v);
  pnl_vect_qsort_index(v, t, 'i');

  for ( i=1 ; i<v->size ; i++ )
    {
      if ( GET (v,i-1) > GET(v,i) || GET(v,i) != GET(vclone, PNL_GET(t,i)) ) 
        {
          pnl_test_set_fail ("vect_qsort_index (increasing)", 0, 0);
          goto J2;
        }
    }
  pnl_test_set_ok ("vect_qsort_index (increasing)");

J2:
  pnl_vect_clone (v, vclone);
  pnl_vect_qsort_index(v, t, 'd');
  for ( i=1 ; i<v->size ; i++ )
    {
      if ( GET (v,i-1) < GET(v,i) || GET(v,i) != GET(vclone, PNL_GET(t,i)) ) 
        {
          pnl_test_set_fail ( "vect_qsort_index (decreasing)", 0, 0);
          goto J3;
        }
    }
  pnl_test_set_ok ("vect_qsort_index (decreasing)");

J3:
  pnl_vect_free (&v);
  pnl_vect_free (&vclone);
  pnl_vect_int_free (&t);
}


int main (int argc, char *argv[])
{
  pnl_test_init (argc, argv);
  pnl_vect_create_from_ptr_test();
  pnl_vect_copy_test();
  pnl_vect_clone_test();
  pnl_vect_plus_vect_test();
  pnl_vect_map_inplace_test();
  pnl_vect_plus_scalar_test();
  pnl_vect_mult_scalar_test();
  pnl_vect_inv_term_test();
  pnl_vect_div_vect_term_test();
  pnl_vect_mult_vect_term_test();
  pnl_vect_sum_test();
  pnl_vect_map_test();
  pnl_vect_scalar_prod_test();
  pnl_vect_prod_test();
  pnl_vect_minmax_test();
  pnl_vect_qsort_test();
  pnl_vect_extract_test();
  pnl_vect_subvect_test();
  exit (pnl_test_finalize ("Vector"));
}
