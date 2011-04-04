
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

#include "pnl/pnl_complex.h"
#include "pnl/pnl_random.h"
#include "tests_utils.h"

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

int main (int argc, char *argv[])
{
  pnl_test_init (argc, argv);
  pnl_vect_complex_set_test();
  pnl_vect_complex_get_test();
  pnl_vect_complex_lget_test();
  pnl_vect_complex_create_from_dcomplex_test();
  pnl_vect_complex_create_from_ptr_test();
  pnl_vect_complex_copy_test();
  pnl_vect_complex_clone_test();
  pnl_vect_complex_plus_vect_test();
  pnl_vect_complex_map_inplace_test();
  pnl_vect_complex_plus_complex_test();
  pnl_vect_complex_mult_complex_test();
  pnl_vect_complex_inv_term_test();
  pnl_vect_complex_div_vect_term_test();
  pnl_vect_complex_mult_vect_term_test();
  pnl_vect_complex_set_complex_test();
  pnl_vect_complex_sum_test();
  pnl_vect_complex_map_test();
  pnl_vect_complex_scalar_prod_test();
  pnl_vect_complex_prod_test();
  exit (pnl_test_finalize ("Complex Vector"));
}

