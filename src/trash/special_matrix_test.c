
/************************************************************************/
/* Copyright David Pommier <david.pommier@gmail.com>                    */
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
#include "pnl/pnl_perm.h"
#include "pnl/pnl_band_matrix.h"
#include "pnl/pnl_spec_matrix.h"
#include "pnl/pnl_random.h"
#include "tests.h"


void Test_Morse_Matrix_algebra(void )
{  
  /*
    Short program to test Morse matrix ...
  */
  int m=6;
  int n=5;
  int Max_size=2;
  int RC=1;
  PnlMorseMat * M;
  PnlMat * FM;
  M=pnl_morse_mat_create(m,n,Max_size,RC);
  pnl_morse_mat_print (M);

  pnl_morse_mat_set(M,2,3,5.0);
  pnl_morse_mat_set(M,1,3,5.0);
  pnl_morse_mat_set(M,2,2,5.0);
  pnl_morse_mat_print (M);
  pnl_morse_mat_set(M,0,4,2.0);
  pnl_morse_mat_set(M,5,4,2.0);
  pnl_morse_mat_freeze(M);
  pnl_morse_mat_print (M);
  
  FM=pnl_morse_mat_full(M);
  pnl_mat_print(FM);

  pnl_morse_mat_free(&M);
  pnl_mat_free(&FM);
}

static void Test_Morse_Matrix_LU(void )
{
  PnlMat *B;
  PnlVect *b,*x,*v;
  int size=15;
  PnlSparseMat * M,*M2;
  PnlMorseMat * MM;
  int gen = 7;
  PnlPermutation *p;
  PnlSparseFactorization *F,*F2;
  pnl_rand_init (gen, size, size);
  b = pnl_vect_create (size);
  v = pnl_vect_create (size);
  B = pnl_mat_create (size, size);
  x = pnl_vect_create(0);
  pnl_mat_set_id (B);
  pnl_vect_rand_normal (v, size, gen);
  pnl_vect_rand_normal (b, size, gen);
  pnl_mat_dger (1., v, v, B);

  M=pnl_sparse_mat_create_fromfull(B);
  printf("test de la fonction 'pnl_sparse_mat_lu_solve' \n");
  p = pnl_permutation_create (size);
  pnl_permutation_init (p);
  F=pnl_sparse_factorization_lu_create (M,0.1);
  
  printf("& test de la fonction 'pnl_sparse_mat_create_frommorse'\n");
  MM=pnl_morse_mat_create_fromfull(B,1);
  M2=pnl_sparse_mat_create_frommorse(MM);
  F2=pnl_sparse_factorization_lu_create (M2,0.1);
  
  pnl_mat_lu (B, p);
  pnl_mat_lu_syslin (x, B, p, b);
  pnl_vect_print(x);
  pnl_vect_clone(x,b);
  pnl_sparse_factorization_lu_syslin (F,x);
  pnl_vect_print(x);
  pnl_vect_clone(x,b);
  pnl_sparse_factorization_lu_syslin (F2,x);
  pnl_vect_print(x);

  
 


  pnl_vect_free(&x);
  pnl_vect_free(&v);
  pnl_vect_free(&b);
  pnl_permutation_free (&p);
  pnl_mat_free (&B);
  pnl_morse_mat_free(&MM);
  pnl_sparse_mat_free(&M);
  pnl_sparse_factorization_free(&F);
  pnl_sparse_mat_free(&M2);
  pnl_sparse_factorization_free(&F2);
 
}

void special_matrix_test()
{
  Test_Morse_Matrix_algebra ();
  Test_Morse_Matrix_LU ();
}
