
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
#include "pnl/pnl_vector_complex.h"
#include "pnl/pnl_random.h"
#include "tests.h"

/* static double function_prod(double x, double y) {return x*y;} */

static void pnl_vect_set_test()
{
  PnlVect *v;
  printf("test de la fonction 'set' : ");
  v=pnl_vect_create_from_double(4,2.0);
  pnl_vect_set(v,1,3.0);
  LET(v,2) = 4.0;
  pnl_vect_print(v);
  pnl_vect_free(&v);
}

static void pnl_vect_get_test()
{
  PnlVect *v;
  printf("test de la fonction 'get' : ");
  v=pnl_vect_create_from_double(4,2.0);
  printf("v[1]=%f \n",pnl_vect_get(v,1));
  printf("v[2]=%f \n",GET(v,2));    
  pnl_vect_free(&v);
}

static void pnl_vect_lget_test()
{
  PnlVect *v;
  printf("test de la fonction 'lget' : ");
  v=pnl_vect_create_from_double(4,2.0);
  printf("v[1]=%f \n",*(pnl_vect_lget(v,1)));
  pnl_vect_free(&v);
}


static void pnl_vect_create_from_double_test()
{
  PnlVect *v;
  int size;
  double x;
  printf("test de la fonction 'pnl_vect_create_from_double' : ");
  size=4;
  x=2.5;
  v=pnl_vect_create_from_double(size,x);
  pnl_vect_print(v);
  pnl_vect_free(&v);
}

static void pnl_vect_create_from_ptr_test()
{
  PnlVect *v;
  int size;
  double x[4]={1.0, 5.0, 3.0, 8.0};
  printf("test de la fonction 'pnl_vect_create_from_ptr' : ");
  size=4;
  v=pnl_vect_create_from_ptr(size,x);
  pnl_vect_print(v);
  pnl_vect_free(&v);
  printf("test de la fonction 'pnl_vect_create_from_list' : ");
  v=pnl_vect_create_from_list(size,1.0, 5.0, 3.0, 8.0);
  pnl_vect_print(v);
  pnl_vect_free(&v);
  
}

static void pnl_vect_copy_test()
{
  PnlVect *v1;
  PnlVect *v2;
  printf("test de la fonction 'pnl_vect_copy' : ");
  v2=pnl_vect_create_from_double(4,3.0);
  v1=pnl_vect_copy(v2);
  pnl_vect_print(v1);
  pnl_vect_free(&v1);
  pnl_vect_free(&v2);
}

static void pnl_vect_clone_test()
{
  PnlVect *v1;
  PnlVect *v2;
  printf("test de la fonction 'pnl_vect_clone' : ");
  v2=pnl_vect_create_from_double(4,3.0);
  v1=pnl_vect_create(0);
  pnl_vect_clone(v1,v2);
  pnl_vect_print(v1);
  pnl_vect_free(&v1);
  pnl_vect_free(&v2); 
}

static void pnl_vect_resize_test ()
{
  PnlVect *v;
  printf ("Test of function pnl_vect_resize : \n");
  v = pnl_vect_create (0);
  printf ("resize 5 : ");
  pnl_vect_resize (v, 5);
  pnl_vect_set_double (v, 0.2);
  pnl_vect_print (v);
  printf ("resize 3 : ");
  pnl_vect_resize (v, 3);
  pnl_vect_print (v);
  printf ("resize 4 : ");
  pnl_vect_resize (v, 4);
  pnl_vect_print (v);
  pnl_vect_free (&v);
}


static void pnl_vect_extract_test()
{
  double x[10]={0.0,1.0, 2.0, 3.0,4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  PnlVect v1;
  PnlVect *v2;
  printf("test de la fonction 'pnl_vect_extract' : ");
  v2=pnl_vect_create_from_ptr(10,x);
  pnl_vect_print(v2);
  v1=pnl_vect_wrap_subvect_with_last(v2,3,5);
  printf("i=3, j=5 \n");
  pnl_vect_print(&v1);
  v1=pnl_vect_wrap_subvect(v2,5,4);
  printf("i=5, size = 4 \n");
  pnl_vect_print(&v1);
  pnl_vect_free(&v2); 
}

static int ispos (double *x) { return *x >= 0; }
static int islarger (double *t) { return t[0] >= t[1]; }

static void pnl_vect_subvect_test ()
{
  int gen=PNL_RNG_MERSENNE_RANDOM_SEED;
  PnlVect *v1, *v2, *v3;
  PnlVectInt *ind;
  printf ("\nTest of extract_subvect function : \n");
  pnl_rand_init (gen, 10, 1);
  v1 = pnl_vect_create (10);
  v3 = pnl_vect_create (10);
  pnl_vect_rand_normal(v1, 10, gen);
  pnl_vect_rand_normal(v3, 10, gen);
  ind = pnl_vect_int_create (10);
  
  pnl_vect_find (ind, "v", ispos, v1);
  v2 = pnl_vect_create_subvect (v1, ind);
  printf ("v1 = "); pnl_vect_print_nsp (v1);
  printf ("find (v1 >= 0)\n");
  printf ("ind = "); pnl_vect_int_print_nsp (ind);
  printf ("sub = "); pnl_vect_print_nsp (v2);
  pnl_vect_free (&v2);

  pnl_vect_find (ind, "vv", islarger, v1, v3);
  v2 = pnl_vect_create_subvect (v1, ind);
  printf ("v1 = "); pnl_vect_print_nsp (v1);
  printf ("v3 = "); pnl_vect_print_nsp (v3);
  printf ("find (v1 >= v3)\n");
  printf ("ind = "); pnl_vect_int_print_nsp (ind);
  printf ("sub = "); pnl_vect_print_nsp (v2);
  pnl_vect_int_free (&ind);
  pnl_vect_free (&v1);
  pnl_vect_free (&v2);
  pnl_vect_free (&v3);
}

static void pnl_vect_plus_vect_test()
{
  PnlVect *v1;
  PnlVect *v2;
  double x[4]={1.0, 5.0, 3.0, 8.0};
  printf("\nTest of pnl_vect_plus_vect : \n");
  v1=pnl_vect_create_from_ptr(4,x);
  v2=pnl_vect_create_from_double(4,3.0);
  pnl_vect_plus_vect(v1,v2);
  pnl_vect_print(v1);
  pnl_vect_free(&v1);
  pnl_vect_free(&v2); 
}

static void pnl_vect_map_inplace_test()
{
  PnlVect *v;
  printf("test de la fonction 'pnl_vect_map_inplace' : ");
  v=pnl_vect_create_from_double(4,3.0);
  pnl_vect_map_inplace(v,exp);
  pnl_vect_print(v);
  pnl_vect_free(&v);
}

static void pnl_vect_plus_double_test()
{
  PnlVect *v;
  printf("test de la fonction 'pnl_vect_plus_double' : ");
  v=pnl_vect_create_from_double(4,3.0);
  pnl_vect_plus_double(v,0.5);
  pnl_vect_print(v);
  pnl_vect_free(&v);
}

static void pnl_vect_mult_double_test()
{
  PnlVect *v;
  printf("test de la fonction 'pnl_vect_mult_double' : ");
  v=pnl_vect_create_from_double(4,3.0);
  pnl_vect_mult_double(v,0.5);
  pnl_vect_print(v);
  pnl_vect_free(&v);
}

static void pnl_vect_inv_term_test()
{
  PnlVect *v;
  printf("test de la fonction 'pnl_vect_inv_term' : ");
  v=pnl_vect_create_from_double(4,3.0);
  pnl_vect_inv_term(v);
  pnl_vect_print(v);
  pnl_vect_free(&v);
}

static void pnl_vect_div_vect_term_test()
{
  PnlVect *v1;
  PnlVect *v2;
  double x[4]={1.0, 5.0, 3.0, 8.0};
  printf("test de la fonction 'pnl_vect_div_term' : ");
  v1=pnl_vect_create_from_ptr(4,x);
  v2=pnl_vect_create_from_double(4,3.0);
  pnl_vect_div_vect_term(v1,v2);
  pnl_vect_print(v1);
  pnl_vect_free(&v1);
  pnl_vect_free(&v2);
}

static void pnl_vect_mult_vect_term_test()
{
  PnlVect *v1;
  PnlVect *v2;
  double x[4]={1.0, 5.0, 3.0, 8.0};
  printf("test de la fonction 'pnl_vect_mult_vect_term' : ");
  v1=pnl_vect_create_from_ptr(4,x);
  v2=pnl_vect_create_from_double(4,3.0);
  pnl_vect_mult_vect_term(v1,v2);
  pnl_vect_print(v1);
  pnl_vect_free(&v1);
  pnl_vect_free(&v2);
}

static void pnl_vect_set_double_test()
{
  PnlVect *v;
  printf("test de la fonction 'pnl_vect_set_double' : ");
  v=pnl_vect_create_from_double(4,3.0);
  pnl_vect_set_double(v,2.0);
  pnl_vect_print(v);
  pnl_vect_free(&v);
}

static void  pnl_vect_sum_test()
{
  PnlVect *v;
  double x[4]={1.0, 5.0, 3.0, 8.0};
  printf("test de la fonction 'pnl_vect_sum' : ");
  v=pnl_vect_create_from_ptr(4,x);
  printf ("v = "); pnl_vect_print (v);
  printf("sum of v = %f \n", pnl_vect_sum(v));
  printf("cumsum of v = ");
  pnl_vect_cumsum(v); pnl_vect_print (v);
  printf("\n");
  pnl_vect_free(&v); 
}

static void pnl_vect_map_test()
{
  PnlVect *v1;
  PnlVect *v2;
  printf("test de la fonction 'pnl_vect_map' : ");
  v1=pnl_vect_create(0);
  v2=pnl_vect_create_from_double(4,3.0);
  pnl_vect_map(v1,v2,exp);
  pnl_vect_print(v1);
  pnl_vect_free(&v1);
  pnl_vect_free(&v2);
}

/* static void pnl_vect_map_vect_test()
 * {
 *   printf("test de la fonction 'pnl_vect_map_vect' : ");
 *   PnlVect *v1;
 *   PnlVect *v2;
 *   double x[4]={1.0, 5.0, 3.0, 8.0};
 *   v1=pnl_vect_create_from_ptr(4,x);
 *   v2=pnl_vect_create_from_double(4,3.0);
 *   pnl_vect_map_vect(v1,v2,function_prod);
 *   pnl_vect_print(v1);
 *   pnl_vect_free(&v1);
 *   pnl_vect_free(&v2);
 * } */

static void pnl_vect_scalar_prod_test()
{
  PnlVect *v1;
  PnlVect *v2;
  double x[4]={1.0, 5.0, 3.0, 8.0};
  printf("test de la fonction 'pnl_vect_scalar_prod' : ");
  v1=pnl_vect_create_from_ptr(4,x);
  v2=pnl_vect_create_from_double(4,3.0);
  printf("produit scalaire de v1 et v2 = %f \n",pnl_vect_scalar_prod(v1,v2));
  pnl_vect_free(&v1);
  pnl_vect_free(&v2);
}

static void pnl_vect_prod_test()
{
  PnlVect *v;
  double x[4]={1.0, 5.0, 3.0, 8.0};
  printf("test de la fonction 'pnl_vect_prod' : ");
  v=pnl_vect_create_from_ptr(4,x);
  printf ("v = "); pnl_vect_print (v);
  printf("prod of v = %f \n", pnl_vect_prod(v));
  printf("cumprod of v = ");
  pnl_vect_cumprod(v); pnl_vect_print (v);
  printf("\n");
  pnl_vect_free(&v);
}

static void pnl_vect_minmax_test()
{
  PnlVect *v;
  double min, max;
  int imin, imax;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  printf("test de la fonction 'pnl_vect_minmax' : ");
  pnl_rand_init (gen, 1, 5);
  v = pnl_vect_create (0);
  printf ("v = ");  pnl_vect_rand_normal (v, 5, gen);
  pnl_vect_print_nsp (v);
  printf("composante maximale du vecteur v = %f \n", pnl_vect_max(v));
  printf("composante minimale du vecteur v = %f \n", pnl_vect_min(v));
  pnl_vect_minmax (v, &min, &max);
  if (min != pnl_vect_min(v) || max != pnl_vect_max(v))
    {
      printf ("error in pnl_vect_{min, max} or pnl_vect_minmax\n");
    }

  pnl_vect_min_index (v, &min, &imin);
  printf("Minimum : indice , valeur = %d, %f \n", imin, min);
  pnl_vect_max_index (v, &max, &imax);
  printf("Maximum : indice , valeur = %d, %f \n", imax, max);
  pnl_vect_minmax_index (v, &min, &max, &imin, &imax);
  printf("(Minimum | maximum) : indice , valeur = (%d, %f |  %d, %f ) \n", imin, min, imax, max);
  
  pnl_vect_free(&v);
}

/* TEST DES FONCTIONS PnlVectComplex */
 
static void pnl_vect_complex_set_test()
{
  PnlVectComplex *v;
  printf("test de la fonction 'set' : ");
  v=pnl_vect_complex_create_from_dcomplex(4,Complex(2.0,1.0));
  pnl_vect_complex_set(v,1,Complex(3.0,4.0));
  pnl_vect_complex_print(v);
  pnl_vect_complex_free(&v);

  printf("test de la fonction 'set_{real,imag}' : ");
  v=pnl_vect_complex_create_from_dcomplex(4,Complex(2.0,1.0));
  pnl_vect_complex_set_real (v,1,3.0);
  pnl_vect_complex_set_imag (v,1,4.0);
  pnl_vect_complex_print(v);
  pnl_vect_complex_free(&v);
}


static void pnl_vect_complex_get_test()
{
  PnlVectComplex *v;
  dcomplex x=Complex(2.0,1.0);
  printf("test de la fonction 'get' : ");
  v=pnl_vect_complex_create_from_dcomplex(4,x);
  printf("v[1]=%f +i %f\n",Creal(pnl_vect_complex_get(v,1)),Cimag(pnl_vect_complex_get(v,1)));
  printf("test de la fonction 'get_{real,imag}' : ");
  printf("v[1]=%f +i %f\n",pnl_vect_complex_get_real(v,1),pnl_vect_complex_get_imag(v,1));
  pnl_vect_complex_free(&v);
}

static void pnl_vect_complex_lget_test()
{
  PnlVectComplex *v;
  printf("test de la fonction 'lget' : ");
  v=pnl_vect_complex_create_from_dcomplex(4,Complex(2.0,1.0));
  printf("v[1]=%f +i %f \n",Creal(*(pnl_vect_complex_lget(v,1))),Cimag(*(pnl_vect_complex_lget(v,1))));
  pnl_vect_complex_free(&v);
}


static void pnl_vect_complex_create_from_dcomplex_test()
{
  PnlVectComplex *v;
  int size;
  dcomplex x;
  printf("test de la fonction 'pnl_vect_complex_create_from_dcomplex' : ");
  size=4;
  x=Complex(2.0,1.0);
  v=pnl_vect_complex_create_from_dcomplex(size,x);
  pnl_vect_complex_print(v);
  pnl_vect_complex_free(&v);
}

static void pnl_vect_complex_create_from_ptr_test()
{
  PnlVectComplex *v;
  int size;
  dcomplex x[4]={{1.0,4.0}, {5.0,3.0}, {6.0,-1.0}, {8.0,2.0}};
  printf("test de la fonction 'pnl_vect_complex_create_from_ptr' : ");
  size=4;
  v=pnl_vect_complex_create_from_ptr(size,x);
  pnl_vect_complex_print(v);
  pnl_vect_complex_free(&v);
}

static void pnl_vect_complex_copy_test()
{
  PnlVectComplex *v1;
  PnlVectComplex *v2;
  printf("test de la fonction 'pnl_vect_complex_copy' : ");
  v2=pnl_vect_complex_create_from_dcomplex(4,Complex(3.0,1.0));
  v1=pnl_vect_complex_copy(v2);
  pnl_vect_complex_print(v1);
  pnl_vect_complex_free(&v1);
  pnl_vect_complex_free(&v2);
}

static void pnl_vect_complex_clone_test()
{
  PnlVectComplex *v1;
  PnlVectComplex *v2;
  printf("test de la fonction 'pnl_vect_complex_clone' : ");
  v2=pnl_vect_complex_create_from_dcomplex(4,Complex(3.0,1.0));
  v1=pnl_vect_complex_create(0);
  pnl_vect_complex_clone(v1,v2);
  pnl_vect_complex_print(v1);
  pnl_vect_complex_free(&v1);
  pnl_vect_complex_free(&v2); 
}

static void pnl_vect_complex_plus_vect_test()
{
  PnlVectComplex *v1;
  PnlVectComplex *v2;
  dcomplex x[4]={{1.0,4.0}, {5.0,3.0}, {6.0,-1.0}, {8.0,2.0}};
  printf("test de la fonction 'pnl_vect_complex_plus_vect' : ");
  v1=pnl_vect_complex_create_from_ptr(4,x);
  v2=pnl_vect_complex_create_from_dcomplex(4,Complex(3.0,2.0));
  pnl_vect_complex_plus_vect(v1,v2);
  pnl_vect_complex_print(v1);
  pnl_vect_complex_free(&v1);
  pnl_vect_complex_free(&v2); 
}

static void pnl_vect_complex_map_inplace_test()
{
  PnlVectComplex *v;
  printf("test de la fonction 'pnl_vect_complex_map_inplace' : ");
  v=pnl_vect_complex_create_from_dcomplex(4,Complex(3.0,2.0));
  pnl_vect_complex_map_inplace(v,Cexp);
  pnl_vect_complex_print(v);
  pnl_vect_complex_free(&v);
}

static void pnl_vect_complex_plus_complex_test()
{
  PnlVectComplex *v;
  printf("test de la fonction 'pnl_vect_complex_plus_complex' : ");
  v=pnl_vect_complex_create_from_dcomplex(4,Complex(3.0,2.0));
  pnl_vect_complex_plus_dcomplex(v,Complex(0.5,0.1));
  pnl_vect_complex_print(v);
  pnl_vect_complex_free(&v);
}

static void pnl_vect_complex_mult_complex_test()
{
  PnlVectComplex *v;
  printf("test de la fonction 'pnl_vect_complex_mult_complex' : ");
  v=pnl_vect_complex_create_from_dcomplex(4,Complex(3.0,2.0));
  pnl_vect_complex_mult_dcomplex(v,Complex(0.5,0.1));
  pnl_vect_complex_print(v);
  pnl_vect_complex_free(&v);
}

static void pnl_vect_complex_inv_term_test()
{
  PnlVectComplex *v;
  printf("test de la fonction 'pnl_vect_complex_inv_term' : ");
  v=pnl_vect_complex_create_from_dcomplex(4,Complex(3.0,2.0));
  pnl_vect_complex_inv_term(v);
  pnl_vect_complex_print(v);
  pnl_vect_complex_free(&v);
}

static void pnl_vect_complex_div_vect_term_test()
{
  PnlVectComplex *v1;
  PnlVectComplex *v2;
  dcomplex x[4]={{1.0,4.0}, {5.0,3.0}, {6.0,-1.0}, {8.0,2.0}};
  printf("test de la fonction 'pnl_vect_complex_div_term' : ");
  v1=pnl_vect_complex_create_from_ptr(4,x);
  v2=pnl_vect_complex_create_from_dcomplex(4,Complex(3.0,2.0));
  pnl_vect_complex_div_vect_term(v1,v2);
  pnl_vect_complex_print(v1);
  pnl_vect_complex_free(&v1);
  pnl_vect_complex_free(&v2);
}

static void pnl_vect_complex_mult_vect_term_test()
{
  PnlVectComplex *v1;
  PnlVectComplex *v2;
  dcomplex x[4]={{1.0,4.0}, {5.0,3.0}, {6.0,-1.0}, {8.0,2.0}};
  printf("test de la fonction 'pnl_vect_complex_mult_vect_term' : ");
  v1=pnl_vect_complex_create_from_ptr(4,x);
  v2=pnl_vect_complex_create_from_dcomplex(4,Complex(3.0,2.0));
  pnl_vect_complex_mult_vect_term(v1,v2);
  pnl_vect_complex_print(v1);
  pnl_vect_complex_free(&v1);
  pnl_vect_complex_free(&v2);
}

static void pnl_vect_complex_set_complex_test()
{
  PnlVectComplex *v;
  printf("test de la fonction 'pnl_vect_complex_set_complex' : ");
  v=pnl_vect_complex_create_from_dcomplex(4,Complex(3.0,2.0));
  pnl_vect_complex_set_dcomplex(v,Complex(2.0,1.0));
  pnl_vect_complex_print(v);
  pnl_vect_complex_free(&v);
}

static void  pnl_vect_complex_sum_test()
{
  PnlVectComplex *v;
  dcomplex x;
  printf("test de la fonction 'pnl_vect_complex_sum' : ");
  v=pnl_vect_complex_create_from_dcomplex(4,Complex(3.0,2.0));
  x=pnl_vect_complex_sum(v);
  printf("somme de v = %f +i %f \n", Creal(x),Cimag(x));
  pnl_vect_complex_free(&v); 
}

static void pnl_vect_complex_map_test()
{
  PnlVectComplex *v1;
  PnlVectComplex *v2;
  printf("test de la fonction 'pnl_vect_complex_map' : ");
  v1=pnl_vect_complex_create(0);
  v2=pnl_vect_complex_create_from_dcomplex(4,Complex(3.0,2.0));
  pnl_vect_complex_map(v1,v2,Cexp);
  pnl_vect_complex_print(v1);
  pnl_vect_complex_free(&v1);
  pnl_vect_complex_free(&v2);
}

/* static void pnl_vect_complex_map_vect_test()
 * {
 *   printf("test de la fonction 'pnl_vect_complex_map_vect' : ");
 *   PnlVectComplex *v1;
 *   PnlVectComplex *v2;
 *   dcomplex x[4]={Complex(1.0,4.0), Complex(5.0,3.0), Complex(6.0,-1.0), Complex(8.0,2.0)};
 *   v1=pnl_vect_complex_create_from_ptr(4,x);
 *   v2=pnl_vect_complex_create_from_dcomplex(4,Complex(3.0,2.0));
 *   pnl_vect_complex_map_vect(v1,v2,Cmul);
 *   pnl_vect_complex_print(v1);
 *   pnl_vect_complex_free(&v1);
 *   pnl_vect_complex_free(&v2);
 * } */

static void pnl_vect_complex_scalar_prod_test()
{
  PnlVectComplex *v1;
  PnlVectComplex *v2;
  dcomplex y;
  dcomplex x[4]={{1.0,4.0}, {5.0,3.0}, {6.0,-1.0}, {8.0,2.0}};
  printf("test de la fonction 'pnl_vect_complex_scalar_prod' : ");
  v1=pnl_vect_complex_create_from_ptr(4,x);
  v2=pnl_vect_complex_create_from_dcomplex(4,Complex(3.0,2.0));
  y=pnl_vect_complex_scalar_prod(v1,v2);
  printf("produit scalaire de v1 et v2 = %f +i %f \n",Creal(y),Cimag(y));
  pnl_vect_complex_free(&v1);
  pnl_vect_complex_free(&v2);
}

static void pnl_vect_complex_prod_test()
{
  PnlVectComplex *v1;
  dcomplex y;
  dcomplex x[4]={{1.0,4.0}, {5.0,3.0}, {6.0,-1.0}, {8.0,2.0}};
  printf("test de la fonction 'pnl_vect_complex_prod' : ");
  v1=pnl_vect_complex_create_from_ptr(4,x);
  y=pnl_vect_complex_prod(v1);
  printf("produit des composantes de v1 = %f+i %f  \n",Creal(y),Cimag(y) );
  pnl_vect_complex_free(&v1);
}

static void pnl_vect_qsort_test ()
{
  PnlVect *v = pnl_vect_create (0);
  PnlVect *vclone = pnl_vect_create (0);
  PnlVectInt *t = pnl_vect_int_create (0);
  printf("test de la fonction 'pnl_vect_qsort' : \n");
  pnl_rand_init (7, 1, 1);
  pnl_vect_rand_uni (v, 20, 0., 1., 7);
  pnl_vect_print_nsp (v);
  pnl_vect_qsort (v, 'i');
  printf("vecteur trie\n");
  pnl_vect_print_nsp (v);
  printf("\n");
  pnl_vect_rand_uni (v, 20, 0., 1., 7);
  pnl_vect_print_nsp (v);
  pnl_vect_clone (vclone, v);
  pnl_vect_qsort_index(v, t, 'i');
  printf("increasing sort\n");
  pnl_vect_print_nsp (v);
  pnl_vect_int_print_nsp (t);
  pnl_vect_clone (v, vclone);
  pnl_vect_qsort_index(v, t, 'd');
  printf("decreasing sort\n");
  pnl_vect_print_nsp (v);
  pnl_vect_int_print_nsp (t);
    
  pnl_vect_free (&v);
  pnl_vect_free (&vclone);
  pnl_vect_int_free (&t);
}

static void all_test ();
static tst_list vect_tests[] =
  {
    MAKE_ENUM(all_test),
    MAKE_ENUM(pnl_vect_set_test),
    MAKE_ENUM(pnl_vect_get_test),
    MAKE_ENUM(pnl_vect_lget_test),
    MAKE_ENUM(pnl_vect_create_from_double_test),
    MAKE_ENUM(pnl_vect_create_from_ptr_test),
    MAKE_ENUM(pnl_vect_copy_test),
    MAKE_ENUM(pnl_vect_clone_test),
    MAKE_ENUM(pnl_vect_resize_test),
    MAKE_ENUM(pnl_vect_plus_vect_test),
    MAKE_ENUM(pnl_vect_map_inplace_test),
    MAKE_ENUM(pnl_vect_plus_double_test),
    MAKE_ENUM(pnl_vect_mult_double_test),
    MAKE_ENUM(pnl_vect_inv_term_test),
    MAKE_ENUM(pnl_vect_div_vect_term_test),
    MAKE_ENUM(pnl_vect_mult_vect_term_test),
    MAKE_ENUM(pnl_vect_set_double_test),
    MAKE_ENUM(pnl_vect_sum_test),
    MAKE_ENUM(pnl_vect_map_test),
    /*MAKE_ENUM(pnl_vect_map_vect_test),*/
    MAKE_ENUM(pnl_vect_scalar_prod_test),
    MAKE_ENUM(pnl_vect_prod_test),
    MAKE_ENUM(pnl_vect_minmax_test),
    MAKE_ENUM(pnl_vect_complex_set_test),
    MAKE_ENUM(pnl_vect_complex_get_test),
    MAKE_ENUM(pnl_vect_complex_lget_test),
    MAKE_ENUM(pnl_vect_complex_create_from_dcomplex_test),
    MAKE_ENUM(pnl_vect_complex_create_from_ptr_test),/* */
    MAKE_ENUM(pnl_vect_complex_copy_test),
    MAKE_ENUM(pnl_vect_complex_clone_test),
    MAKE_ENUM(pnl_vect_complex_plus_vect_test),/* */
    MAKE_ENUM(pnl_vect_complex_map_inplace_test),
    MAKE_ENUM(pnl_vect_complex_plus_complex_test),
    MAKE_ENUM(pnl_vect_complex_mult_complex_test),
    MAKE_ENUM(pnl_vect_complex_inv_term_test),
    MAKE_ENUM(pnl_vect_complex_div_vect_term_test),
    MAKE_ENUM(pnl_vect_complex_mult_vect_term_test),
    MAKE_ENUM(pnl_vect_complex_set_complex_test),
    MAKE_ENUM(pnl_vect_complex_sum_test),
    MAKE_ENUM(pnl_vect_complex_map_test),
    /*MAKE_ENUM(pnl_vect_complex_map_vect_test),*/
    MAKE_ENUM(pnl_vect_complex_scalar_prod_test),
    MAKE_ENUM(pnl_vect_complex_prod_test),
    MAKE_ENUM(pnl_vect_qsort_test),
    MAKE_ENUM(pnl_vect_extract_test),
    MAKE_ENUM(pnl_vect_subvect_test),
    MAKE_ENUM(NULL)
  };

static void all_test ()
{
  run_all_test (vect_tests);
}

void vector_test()
{
  menu_test (vect_tests);
}
