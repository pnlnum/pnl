
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

#include "config.h"
#include "pnl_matrix.h"
#include "pnl_matrix_complex.h"
#include "pnl_random.h"
#include "tests.h"

/* static double function_prod(double x, double y) {return x*y;} */


/** 
 * Sets the upper (uplo='U') or lower (uplo='L') triangle to 0
 * 
 * @param uplo 'U' or 'L'
 * @param A a real matrix
 */
static void set_triangular_to_zero (PnlMat *A, char uplo)
{
  int i, j;

  if ( uplo == 'L')
    {
      for ( i=0 ; i <A->m ; i++ )
        {
          for ( j=0 ; j<i ; j++ )
            {
              PNL_MLET (A, i, j) = 0;
            }
        }
    }
  else
    {
      for ( i=0 ; i <A->m ; i++ )
        {
          for ( j=i+1 ; j<A->n ; j++ )
            {
              PNL_MLET (A, i, j) = 0;
            }
        }
    }
}


/** 
 * Creates an invertible matrix of size nxn using the exponential function
 * 
 * @param A Output parameter, contains a n x n intervible matrix.
 * @param n size of the matrix
 * @param gen index of the generator
 */
void create_invertible_matrix (PnlMat *A, int n, int gen)
{
  PnlMat *B;

  B = pnl_mat_create (n, n);
  pnl_mat_rand_uni2 (B, n, n, 0., 1., gen);

  pnl_mat_exp (A, B);
  pnl_mat_free (&B);
}

/** 
 * Creates a symmetric positive definite matrix of size nxn
 * 
 * @param S Output parameter, contains a n x n symmetric positive definite matrix.
 * @param n size of the matrix
 * @param gen index of the generator
 */
void create_sym_pos_matrix (PnlMat *S, int n, int gen)
{
  PnlVect *b;
  double g;
  int i;

  b = pnl_vect_create (n);
  pnl_vect_rand_uni (b, n, 0., 1., gen);

  /* S is set to a diagonal matrix with positive eigenvalues */
  pnl_mat_set_double (S, 0.);
  for ( i=0 ; i<S->n ; i++ )
    {
      do
        {
          g = fabs (pnl_rand_normal (gen));
        }
      while (g < 1E-4);
      pnl_mat_set (S, i, i, g);
    }

  pnl_mat_dger (1., b, b, S);
  pnl_vect_free (&b);
}


/* Test of the PnlMat functions */

static void pnl_mat_set_test()
{
  PnlMat *M;
  printf("test de la fonction 'set' : \n");
  M=pnl_mat_create_from_double(4,2,2.0);
  MLET (M,3,1) = 3.0;
  pnl_mat_print(M);
  pnl_mat_free(&M);
}

static void pnl_mat_get_test()
{
  PnlMat *M;
  printf("test de la fonction 'get' : ");
  M=pnl_mat_create_from_double(4,2,2.0);
  printf("M[3,1]=%f \n",MGET(M,3,1));
  pnl_mat_free(&M);
}

static void pnl_mat_lget_test()
{
  PnlMat *M;
  printf("test de la fonction 'lget' : ");
  M=pnl_mat_create_from_double(4,2,2.0);
  printf("M[3,1]=%f \n",*(pnl_mat_lget(M,3,1)));
  pnl_mat_free(&M);
}

static void pnl_mat_create_from_double_test()
{
  PnlMat *M;
  int rows;
  int cols;
  double x;
  printf("test de la fonction 'pnl_mat_create_from_double' : \n");
  rows=4;
  cols=2;
  x=2.5;
  M=pnl_mat_create_from_double(rows,cols,x);
  pnl_mat_print(M);
  pnl_mat_free(&M);
}

static void pnl_mat_create_from_ptr_test()
{
  PnlMat *M;
  int rows;
  int cols;
  double x[8]={1.0, 5.0, 3.0, 8.0, 3.0, 7.0, 6.0, 9.0};
  printf("test de la fonction 'pnl_mat_create_from_ptr' : \n");
  rows=4;
  cols=2;
  M=pnl_mat_create_from_ptr(rows,cols,x);
  pnl_mat_print(M);
  pnl_mat_free(&M);
  M=pnl_mat_create_from_list(rows,cols,1.0, 5.0, 3.0, 8.0, 3.0, 7.0, 6.0, 9.0);
  printf("test de la fonction 'pnl_mat_create_from_list' : \n");
  pnl_mat_print(M);
  pnl_mat_free(&M);
}

static int ispos (double x[]) { return x[0] >= 0; }
static int islarger (double x[]) { return x[0] >= x[1]; }
static int cmp (double x[]) { return (x[0] >= x[1]) && (x[1] < x[2]); }

static void pnl_mat_submat_test ()
{
  double gen=PNL_RNG_MERSENNE_RANDOM_SEED;
  int m = 5, n = 7;
  PnlMat *M1, *M2;
  PnlVect *v1;
  PnlVectInt *indi, *indj;
  printf ("\nTest of extract_submat function : \n");
  pnl_rand_init (gen, m, n);
  M1 = pnl_mat_create (m ,n);
  M2 = pnl_mat_create (m, n);
  pnl_mat_rand_normal(M1, m, n, gen);
  pnl_mat_rand_normal(M2, m, n, gen);
  indi = pnl_vect_int_create (0);
  indj = pnl_vect_int_create (0);
  
  if ( pnl_mat_find (indi, indj, "m", ispos, M1) == FAIL ) 
    {
      printf ("Error in pnl_mat_find\n");
    }
  v1 = pnl_vect_create_submat (M1, indi, indj);
  printf ("M1 = "); pnl_mat_print_nsp (M1);
  printf ("[i, j] = find (M1 >= 0)\n");
  printf ("indi = "); pnl_vect_int_print_nsp (indi);
  printf ("indj = "); pnl_vect_int_print_nsp (indj);
  printf ("sub = "); pnl_vect_print_nsp (v1);
  printf ("\n");
  pnl_vect_free (&v1);

  if ( pnl_mat_find (indi, NULL, "m", ispos, M1) == FAIL )
    {
      printf ("Error in pnl_mat_find\n");
    }
  printf ("M1 = "); pnl_mat_print_nsp (M1);
  printf ("[i] = find (M1 >= 0)\n");
  printf ("indi = "); pnl_vect_int_print_nsp (indi);
  printf ("\n");

  if ( pnl_mat_find (indi, indj, "mm", islarger, M1, M2) == FAIL ) 
    {
      printf ("Error in pnl_mat_find\n");
    }
  v1 = pnl_vect_create_submat (M1, indi, indj);
  printf ("M1 = "); pnl_mat_print_nsp (M1);
  printf ("M2 = "); pnl_mat_print_nsp (M2);
  printf ("[i, j] = find (M1 >= M2)\n");
  printf ("indi = "); pnl_vect_int_print_nsp (indi);
  printf ("indj = "); pnl_vect_int_print_nsp (indj);
  printf ("sub = "); pnl_vect_print_nsp (v1);
  printf ("\n");
  pnl_vect_free (&v1);

  pnl_mat_find (indi, indj, "mmr", cmp, M1, M2, 0.);
  v1 = pnl_vect_create_submat (M1, indi, indj);
  printf ("M1 = "); pnl_mat_print_nsp (M1);
  printf ("M2 = "); pnl_mat_print_nsp (M2);
  printf ("[i, j] = find ((M1 >= M2) && (M2 < 0))\n");
  printf ("indi = "); pnl_vect_int_print_nsp (indi);
  printf ("indj = "); pnl_vect_int_print_nsp (indj);
  printf ("sub = "); pnl_vect_print_nsp (v1);
  printf ("\n");

  pnl_vect_free (&v1);
  pnl_vect_int_free (&indi);
  pnl_vect_int_free (&indj);
  pnl_mat_free (&M1);
  pnl_mat_free (&M2);
}

static void pnl_mat_set_diag_test ()
{
  int d, n;
  double x;
  PnlMat *M;

  n = 8;
  x = 5.;
  M = pnl_mat_create (n, n);

  d = 0;
  pnl_mat_set_double (M, 0.);
  pnl_mat_set_diag (M, x, d);
  pnl_mat_print (M);
  printf ("\n");

  d = 3;
  pnl_mat_set_double (M, 0.);
  pnl_mat_set_diag (M, x, d);
  pnl_mat_print (M);
  printf ("\n");

  d = -2;
  pnl_mat_set_double (M, 0.);
  pnl_mat_set_diag (M, x, d);
  pnl_mat_print (M);
  pnl_mat_free (&M);
  printf ("\n");
}

static void pnl_mat_copy_test()
{
  PnlMat *M1;
  PnlMat *M2;
  printf("test de la fonction 'pnl_mat_copy' : \n");
  M2=pnl_mat_create_from_double(4,2,3.0);
  M1=pnl_mat_copy(M2);
  pnl_mat_print(M1);
  pnl_mat_free(&M1);
  pnl_mat_free(&M2);
}

static void pnl_mat_clone_test()
{
  PnlMat *M1;
  PnlMat *M2;
  printf("test de la fonction 'pnl_mat_clone' : \n");
  M2=pnl_mat_create_from_double(4,2,3.0);
  M1=pnl_mat_create(0,0);
  pnl_mat_clone(M1,M2);
  pnl_mat_print(M1);
  pnl_mat_free(&M1);
  pnl_mat_free(&M2); 
}


static void pnl_mat_resize_test ()
{
  PnlMat *v;
  printf ("Test of function pnl_mat_resize : \n");
  v = pnl_mat_create (0, 0);
  printf ("resize 5x3 : ");
  pnl_mat_resize (v, 5, 3);
  pnl_mat_set_double (v, 0.2);
  pnl_mat_print (v);
  printf ("resize 3x5 : ");
  pnl_mat_resize (v, 3, 5);
  pnl_mat_print (v);
  printf ("resize 2x3 : ");
  pnl_mat_resize (v, 2, 3);
  pnl_mat_print (v);
  printf ("resize 3x4 : ");
  pnl_mat_resize (v, 3, 4);
  pnl_mat_print (v);
  pnl_mat_free (&v);
}

static void pnl_mat_map_inplace_test()
{
  PnlMat *M;
  printf("test de la fonction 'pnl_mat_map_inplace' : \n");
  M=pnl_mat_create_from_double(4,2,3.0);
  pnl_mat_map_inplace(M,exp);
  pnl_mat_print(M);
  pnl_mat_free(&M);
}

static void pnl_mat_plus_mat_test()
{
  PnlMat *M1;
  PnlMat *M2;
  double rows;
  double cols;
  double x[8]={1.0, 5.0, 3.0, 8.0, 3.0, 7.0, 6.0, 9.0};
  printf("test de la fonction 'pnl_mat_plus_mat' : \n");
  rows=4;
  cols=2;
  M1=pnl_mat_create_from_ptr(rows,cols,x);
  M2=pnl_mat_create_from_double(4,2,3.0);
  pnl_mat_plus_mat(M1,M2);
  pnl_mat_print(M1);
  pnl_mat_free(&M1);
  pnl_mat_free(&M2); 
}

static void pnl_mat_plus_double_test()
{
  PnlMat *M;
  printf("test de la fonction 'pnl_mat_plus_double' : \n");
  M=pnl_mat_create_from_double(4,2,3.0);
  pnl_mat_plus_double(M,0.5);
  pnl_mat_print(M);
  pnl_mat_free(&M);
}

static void pnl_mat_mult_mat_test()
{
  PnlMat *M, *M1, *M2;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  printf("test de la fonction 'pnl_mat_mult_mat' : \n");
  pnl_rand_init (gen, 5, 5);
  M1 = pnl_mat_create (5, 4);
  M2 = pnl_mat_create (4, 5);
  pnl_mat_rand_normal (M1, 5, 4, gen);
  pnl_mat_rand_normal (M2, 4, 5, gen);
  printf ("A = "); pnl_mat_print_nsp (M1);
  printf ("B = "); pnl_mat_print_nsp (M2);
  M=pnl_mat_mult_mat(M1,M2);
  pnl_mat_print(M);
  pnl_mat_free(&M);
  pnl_mat_free(&M1);
  pnl_mat_free(&M2);
}

static void pnl_mat_mult_mat_inplace_test()
{
  PnlMat *M, *M1, *M2;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  printf("test de la fonction 'pnl_mat_mult_mat' : \n");
  pnl_rand_init (gen, 5, 5);
  M1 = pnl_mat_create (5, 4);
  M2 = pnl_mat_create (4, 5);
  M = pnl_mat_create (0, 0);
  pnl_mat_rand_normal (M1, 5, 4, gen);
  pnl_mat_rand_normal (M2, 4, 5, gen);
  printf ("A = "); pnl_mat_print_nsp (M1);
  printf ("B = "); pnl_mat_print_nsp (M2);
  pnl_mat_mult_mat_inplace(M,M1,M2);
  pnl_mat_print(M);
  pnl_mat_free(&M);
  pnl_mat_free(&M1);
  pnl_mat_free(&M2);
}

static void pnl_mat_mult_vect_test()
{
  PnlMat *A;
  PnlVect *x, *y;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  printf("test de la fonction pnl_mat_mult_vect : \n");
  pnl_rand_init (gen, 5, 4);
  A = pnl_mat_create (5, 4);
  printf ("A = ");  pnl_mat_rand_normal (A, 5, 4, gen);
  pnl_mat_print_nsp (A);
  printf("\n");
  x = pnl_vect_create (4);
  pnl_vect_rand_normal(x, 4, gen);
  printf ("x = "); pnl_vect_print_nsp (x);
  printf("\n");

  
  y = pnl_vect_create (0);
  pnl_mat_mult_vect_inplace (y, A, x);
  pnl_vect_print_nsp (y);
  printf("\n");
  pnl_vect_free (&y);
  y = pnl_mat_mult_vect (A, x);
  pnl_vect_print_nsp (y);
  printf("\n");

  pnl_vect_free (&x);
  pnl_vect_free (&y);
  pnl_mat_free (&A);
}

static void pnl_mat_mult_double_test()
{
  PnlMat *M;
  printf("test de la fonction 'pnl_mat_mult_double' : \n");
  M=pnl_mat_create_from_double(4,2,3.0);
  pnl_mat_mult_double(M,0.5);
  pnl_mat_print(M);
  pnl_mat_free(&M);
}

static void pnl_mat_set_double_test()
{
  PnlMat *M;
  printf("test de la fonction 'pnl_mat_set_double' : \n");
  M=pnl_mat_create_from_double(4,2,3.0);
  pnl_mat_set_double(M,2.0);
  pnl_mat_print(M);
  pnl_mat_free(&M);
}


static void pnl_mat_mult_mat_term_test()
{
  PnlMat *M1, *M2;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  printf("test de la fonction 'pnl_mat_mult_mat_term' : \n");
  pnl_rand_init (gen, 5, 5);
  M1 = pnl_mat_create (5, 4);
  M2 = pnl_mat_create (5, 4);
  pnl_mat_rand_normal (M1, 5, 4, gen);
  pnl_mat_rand_normal (M2, 5, 4, gen);
  printf ("A = "); pnl_mat_print_nsp (M1);
  printf ("B = "); pnl_mat_print_nsp (M2);
  pnl_mat_mult_mat_term(M1,M2);
  printf ("A .* B = "); pnl_mat_print_nsp (M1);
  pnl_mat_free(&M1);
  pnl_mat_free(&M2);
}

static void pnl_mat_chol_test()
{
  PnlMat *S;
  PnlVect *x;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  printf("test de la fonction 'pnl_mat_chol' : \n");
  pnl_rand_init (gen, 5, 5);
  x = pnl_vect_create (5);
  S = pnl_mat_create (5, 5);
  pnl_mat_set_id (S);
  pnl_vect_rand_normal (x, 5, gen);
  pnl_mat_dger (1., x, x, S);
  printf ("S = "); pnl_mat_print_nsp (S);
  pnl_mat_chol(S);
  printf ("\nchol (S) = "); pnl_mat_print_nsp (S);
  printf ("\n");
  pnl_vect_free(&x);
  pnl_mat_free(&S);
}

static void pnl_mat_sq_transpose_test()
{
  PnlMat *M;
  double x[9]={3.0, 1.0, 4.0, 5.0, 3.0, 6.0, 2.0, 9.0, 3.0};
  printf("test de la fonction 'pnl_mat_sq_transpose' : \n");
  M=pnl_mat_create_from_ptr(3,3,x);
  pnl_mat_sq_transpose(M);
  pnl_mat_print(M);
  pnl_mat_free(&M);
}

static void pnl_mat_transpose_test()
{
  PnlMat *M,*M1;
  double x[8]={1.0, 5.0, 3.0, 8.0, 3.0, 7.0, 6.0, 9.0};
  printf("test de la fonction 'pnl_mat_transpose' : \n");
  M=pnl_mat_create_from_ptr(4,2,x);
  M1=pnl_mat_transpose(M);
  pnl_mat_print(M1);
  pnl_mat_free(&M);
  pnl_mat_free(&M1);
}

static void pnl_vect_wrap_mat_row_test()
{
  PnlMat *M;
  PnlVect V;
  double x[8]={1.0, 5.0, 3.0, 8.0, 3.0, 7.0, 6.0, 9.0};
  printf("test de la fonction 'pnl_vect_wrap_mat_row : \n");
  M=pnl_mat_create_from_ptr(4,2,x);
  V=pnl_vect_wrap_mat_row(M,2);
  pnl_vect_print(&V);
  pnl_mat_free(&M);
}



static void pnl_mat_get_row_test()
{
  PnlMat *M;
  PnlVect *V;
  double x[8]={1.0, 5.0, 3.0, 8.0, 3.0, 7.0, 6.0, 9.0};
  printf("test de la fonction 'pnl_mat_row_to_vect_inplace : \n");
  M=pnl_mat_create_from_ptr(4,2,x);
  V=pnl_vect_create(0);
  pnl_mat_get_row(V,M,2);
  pnl_vect_print(V);
  pnl_vect_free(&V);
  pnl_mat_free(&M);
}

static void pnl_mat_create_diag_test()
{
  PnlMat *M;
  PnlVect *V;
  printf("test de la fonction 'pnl_mat_create_diag : \n");
  V=pnl_vect_create_from_double(3,1.0);
  M=pnl_mat_create_diag(V);
  pnl_mat_print(M);
  pnl_vect_free(&V);
  pnl_mat_free(&M);
}

static void pnl_mat_sum_test()
{
  PnlMat *M, *Mcopy;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  printf("test de la fonction 'pnl_mat_sum : ");
  pnl_rand_init (gen, 5, 5);
  M = pnl_mat_create (5, 4);
  pnl_mat_rand_normal (M, 5, 4, gen);
  printf ("M = "); pnl_mat_print_nsp (M);
  printf("la somme des elements de la matrice est %f \n", pnl_mat_sum(M));
  pnl_mat_free(&M);

  printf("test de la fonction 'pnl_mat_cumsum : \n");
  M = pnl_mat_create (5, 4);
  pnl_mat_rand_normal (M, 5, 4, gen);
  Mcopy = pnl_mat_copy (M);
  printf ("M = "); pnl_mat_print_nsp (M);
  pnl_mat_cumsum (M, 'r');
  printf("\ncumsum (M, 'r') = \n"); pnl_mat_print_nsp (M);
  pnl_mat_cumsum (Mcopy, 'c');
  printf("\ncumsum (M, 'c') = \n"); pnl_mat_print_nsp (Mcopy);
  pnl_mat_free(&M);
  pnl_mat_free(&Mcopy);
}

static void pnl_mat_prod_test()
{
  PnlMat *M, *Mcopy;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  printf("test de la fonction 'pnl_mat_prod : ");
  pnl_rand_init (gen, 5, 5);
  M = pnl_mat_create (5, 4);
  pnl_mat_rand_normal (M, 5, 4, gen);
  printf ("M = "); pnl_mat_print_nsp (M);
  printf("prod = %f \n", pnl_mat_prod(M));
  pnl_mat_free(&M);

  printf("test de la fonction 'pnl_mat_cumprod : \n");
  M = pnl_mat_create (5, 4);
  pnl_mat_rand_normal (M, 5, 4, gen);
  Mcopy = pnl_mat_copy (M);
  printf ("M = "); pnl_mat_print_nsp (M);
  pnl_mat_cumprod (M, 'r');
  printf("\ncumprod (M, 'r') = \n"); pnl_mat_print_nsp (M);
  pnl_mat_cumprod (Mcopy, 'c');
  printf("\ncumprod (M, 'c') = \n"); pnl_mat_print_nsp (Mcopy);
  pnl_mat_free(&M);
  pnl_mat_free(&Mcopy);

}


static void pnl_mat_sum_vect_test()
{
  PnlMat *M;
  PnlVect *V;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  printf("test de la fonction 'pnl_mat_sum_test : \n");
  pnl_rand_init (gen, 5, 5);
  M = pnl_mat_create (5, 4);
  V = pnl_vect_create (0);
  pnl_mat_rand_normal (M, 5, 4, gen);
  printf ("A = "); pnl_mat_print_nsp (M);
  pnl_mat_sum_vect(V, M,'c');
  printf("\nsum (A, 'c') = "); pnl_vect_print(V);
  pnl_mat_sum_vect(V, M,'r');
  printf("\nsum (A, 'r') = "); pnl_vect_print(V);
  pnl_vect_free(&V);
  pnl_mat_free(&M);
}

static void pnl_mat_prod_vect_test()
{
  PnlMat *M;
  PnlVect *V;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  printf("test de la fonction 'pnl_mat_prod_test : \n");
  pnl_rand_init (gen, 5, 5);
  M = pnl_mat_create (5, 4);
  V = pnl_vect_create (0);
  pnl_mat_rand_normal (M, 5, 4, gen);
  printf ("A = "); pnl_mat_print_nsp (M);
  pnl_mat_prod_vect(V, M,'c');
  printf("\nprod (A, 'c') = "); pnl_vect_print(V);
  pnl_mat_prod_vect(V, M,'r');
  printf("\nprod (A, 'r') = "); pnl_vect_print(V);
  pnl_vect_free(&V);
  pnl_mat_free(&M);
}

static void pnl_mat_minmax_test()
{
  PnlMat *A;
  PnlVect *min, *max;
  PnlVectInt *imin, *imax;
  int m, n;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  
  
  printf("test de la fonction 'pnl_vect_minmax' : ");
  m = 5; n = 6;
  pnl_rand_init (gen, 1, m*n);
  A = pnl_mat_create (m, n);
  min = pnl_vect_create (0);
  max = pnl_vect_create (0);
  imin = pnl_vect_int_create (0);
  imax = pnl_vect_int_create (0);
  
  printf ("A = \n");  pnl_mat_rand_normal (A, m, n, gen);
  pnl_mat_print_nsp (A);
  printf ("--> pnl_mat_min \n");
  pnl_mat_min(A, min, 'c');
  printf("min(A, 'c') = "); pnl_vect_print_nsp (min);
  pnl_mat_min(A, min, 'r');
  printf("min(A, 'r') = "); pnl_vect_print_nsp (min);

  printf ("--> pnl_mat_max \n");
  pnl_mat_max(A, max, 'c');
  printf("max(A, 'c') = "); pnl_vect_print_nsp (max);
  pnl_mat_max(A, max, 'r');
  printf("max(A, 'r') = "); pnl_vect_print_nsp (max);

  printf ("--> pnl_mat_minmax \n");
  pnl_mat_minmax (A, min, max, 'r');
  printf("min(A, 'r') = "); pnl_vect_print_nsp (min);
  printf("max(A, 'r') = "); pnl_vect_print_nsp (max);

  printf ("--> pnl_mat_min_index \n");
  pnl_mat_min_index(A, min, imin, 'c');
  printf("min(A, 'c') = "); pnl_vect_print_nsp (min);
  printf("\tindex = "); pnl_vect_int_print_nsp (imin);
  pnl_mat_min_index(A, min, imin, 'r');
  printf("min(A, 'r') = "); pnl_vect_print_nsp (min);
  printf("\tindex = "); pnl_vect_int_print_nsp (imin);

  printf ("--> pnl_mat_max_index \n");
  pnl_mat_max_index(A, max, imax, 'c');
  printf("max(A, 'c') = "); pnl_vect_print_nsp (max);
  printf("\tindex = "); pnl_vect_int_print_nsp (imax);
  pnl_mat_max_index(A, max, imax, 'r');
  printf("max(A, 'r') = "); pnl_vect_print_nsp (max);
  printf("\tindex = "); pnl_vect_int_print_nsp (imax);

  printf ("--> pnl_mat_minmax_index \n");
  pnl_mat_minmax_index (A, min, max, imin, imax, 'r');
  printf("max(A, 'r') = "); pnl_vect_print_nsp (max);
  printf("\tindex = "); pnl_vect_int_print_nsp (imax);
  printf("min(A, 'r') = "); pnl_vect_print_nsp (min);
  printf("\tindex = "); pnl_vect_int_print_nsp (imin);

  pnl_mat_free(&A);
  pnl_vect_free (&min);   pnl_vect_free (&max);
  pnl_vect_int_free (&imin);   pnl_vect_int_free (&imax);
}

static void pnl_mat_qsort_test ()
{
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  PnlMat *M = pnl_mat_create (0, 0);
  PnlMat *Mclone = pnl_mat_create (0, 0);
  PnlMatInt *t = pnl_mat_int_create (0, 0);
  printf("test de la fonction 'pnl_mat_qsort' : \n");
  pnl_rand_init (gen, 1, 1);
  pnl_mat_rand_normal (M, 11, 9, gen);
  pnl_mat_print_nsp (M);
  pnl_mat_clone (Mclone, M);

  pnl_mat_qsort_index (M, t, 'r', 'i');
  printf("sort(M, 'r', 'i') \n");
  pnl_mat_print_nsp (M);
  pnl_mat_int_print_nsp (t);
  printf("\n");

  pnl_mat_clone (M, Mclone);
  pnl_mat_qsort_index (M, t, 'r', 'd');
  printf("sort(M, 'r', 'd') \n");
  pnl_mat_print_nsp (M);
  pnl_mat_int_print_nsp (t);
  printf("\n");

  pnl_mat_clone (M, Mclone);
  pnl_mat_qsort_index (M, t, 'c', 'd');
  printf("sort(M, 'c', 'd') \n");
  pnl_mat_print_nsp (M);
  pnl_mat_int_print_nsp (t);
  printf("\n");

  pnl_mat_clone (M, Mclone);
  pnl_mat_qsort_index (M, t, 'c', 'i');
  printf("sort(M, 'c', 'i') \n");
  pnl_mat_print_nsp (M);
  pnl_mat_int_print_nsp (t);
  printf("\n");

  pnl_mat_free (&M);
  pnl_mat_free (&Mclone);
  pnl_mat_int_free (&t);
}

static void pnl_mat_map_test()
{
  PnlMat *M1,*M2;
  double x[8]={1.0, 5.0, 3.0, 8.0, 3.0, 7.0, 6.0, 9.0};
  printf("test de la fonction 'pnl_mat_map': \n ");
  M1=pnl_mat_create(0,0);
  M2=pnl_mat_create_from_ptr(4,2,x);
  pnl_mat_map(M1,M2,exp);
  pnl_mat_print(M1);
  pnl_mat_free(&M1);
  pnl_mat_free(&M2);
}

static void pnl_mat_set_row_test()
{
  PnlMat *M;
  PnlVect *V;
  double x[8]={1.0, 5.0, 3.0, 8.0, 3.0, 7.0, 6.0, 9.0};
  printf("test de la fonction 'pnl_mat_set_row': \n ");
  M=pnl_mat_create_from_ptr(4,2,x);
  V=pnl_vect_create_from_double(2,0.0);
  pnl_mat_set_row(M,V,2);
  pnl_mat_print(M);
  pnl_vect_free(&V);
  pnl_mat_free(&M);
}

static void pnl_mat_triangular_inverse_test ()
{
  PnlMat *A, *B;
  int gen;
  int n = 5;

  gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init (gen, n, n);
  A = pnl_mat_create (n, n);
  B = pnl_mat_create (n, n);
  create_invertible_matrix (B, n, gen);

  set_triangular_to_zero (B, 'L');
  printf("test de la fonction 'pnl_mat_upper_inverse' : \n");
  printf ("B = "); pnl_mat_print_nsp (B); printf ("\n");
  pnl_mat_upper_inverse(A,B);
  printf ("invB = "); pnl_mat_print_nsp (A); printf ("\n");

  printf("test de la fonction 'pnl_mat_lower_inverse' : \n");
  create_invertible_matrix (B, n, gen);
  set_triangular_to_zero (B, 'U');
  printf ("B = "); pnl_mat_print_nsp (B); printf ("\n");
  pnl_mat_lower_inverse(A,B);
  printf ("invB = "); pnl_mat_print_nsp (A); printf ("\n");
  pnl_mat_free(&A);
  pnl_mat_free(&B);
}

static void pnl_mat_inverse_test ()
{
  PnlMat *A, *invA;
  int gen;
  int n = 5;

  gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init (gen, n, n);
  A = pnl_mat_create (n, n);
  invA = pnl_mat_create (n, n);

  printf("test de la fonction 'pnl_mat_inverse_with_chol' : \n");
  create_sym_pos_matrix (A, n, gen);
  pnl_mat_inverse_with_chol (invA, A);
  printf ("A = "); pnl_mat_print_nsp (A);
  printf ("invA = "); pnl_mat_print_nsp (invA);

  printf("test de la fonction 'pnl_mat_inverse' : \n");
  create_invertible_matrix (A, n, gen);
  pnl_mat_inverse (invA, A);
  printf ("A = "); pnl_mat_print_nsp (A);
  printf ("invA = "); pnl_mat_print_nsp (invA);

  pnl_mat_free (&A);
  pnl_mat_free (&invA);
}

static void pnl_mat_syslin_test ()
{
  PnlMat *S, *Scopy, *B, *Bcopy;
  PnlVect *b, *x;
  PnlVectInt *p;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  printf("test de la fonction 'pnl_mat_chol' : \n");
  pnl_rand_init (gen, 5, 5);
  b = pnl_vect_create (5);
  S = pnl_mat_create (5, 5);
  pnl_vect_rand_normal (b, 5, gen);
  create_sym_pos_matrix (S, 5, gen);

  Scopy = pnl_mat_copy (S);
  printf ("S = "); pnl_mat_print_nsp (S);
  printf ("\nb = "); pnl_vect_print_nsp (b);
  printf ("\n");
  x = pnl_vect_create(0);
  printf("test de la fonction 'pnl_mat_upper_syslin' : \n");
  pnl_mat_upper_syslin(x,S,b);
  pnl_vect_print(x);
    
  printf("test de la fonction 'pnl_mat_lower_syslin' : \n");
  pnl_mat_lower_syslin(x,S,b);
  pnl_vect_print(x);

  printf("test de la fonction 'pnl_mat_chol_syslin' : \n");
  pnl_mat_chol (S);
  pnl_mat_chol_syslin(x, S, b);
  pnl_vect_print(x);

  printf("test de la fonction 'pnl_mat_lu_syslin' (symmetric matrix) : \n");
  pnl_mat_clone (S, Scopy);
  p = pnl_vect_int_create (5);
  pnl_mat_lu (S, p);
  pnl_mat_lu_syslin (x, S, p, b);
  pnl_vect_print(x);

  printf("test de la fonction 'pnl_mat_syslin' (symmetric matrix) : \n");
  pnl_mat_clone (S, Scopy);
  pnl_mat_syslin (x, S, b);
  pnl_vect_print(x);

  B = pnl_mat_create (5,5);
  pnl_mat_rand_normal (B, 5, 5, gen);
  Bcopy = pnl_mat_copy (B);
  printf("test de la fonction 'pnl_mat_chol_syslin_mat' : \n");
  pnl_mat_clone (S, Scopy);
  printf ("B = "); pnl_mat_print_nsp (B); printf ("\n");
  pnl_mat_chol (S);
  pnl_mat_chol_syslin_mat (S, B);
  printf ("X = "); pnl_mat_print_nsp (B); printf ("\n");

  printf("test de la fonction 'pnl_mat_lu_syslin_mat' (symmetric matrix): \n");
  pnl_mat_clone (S, Scopy);
  pnl_mat_clone (B, Bcopy);
  printf ("B = "); pnl_mat_print_nsp (B); printf ("\n");
  pnl_mat_lu (S, p);
  pnl_mat_lu_syslin_mat (S, p, B);
  printf ("X = "); pnl_mat_print_nsp (B); printf ("\n");

  printf("test de la fonction 'pnl_mat_syslin_mat' (symmetric matrix): \n");
  pnl_mat_clone (S, Scopy);
  pnl_mat_clone (B, Bcopy);
  printf ("B = "); pnl_mat_print_nsp (B); printf ("\n");
  pnl_mat_syslin_mat (S, B);
  printf ("X = "); pnl_mat_print_nsp (B); printf ("\n");


  printf("test de la fonction 'pnl_mat_lu_syslin' : \n");
  create_invertible_matrix (S, 5, gen);
  pnl_mat_clone (Scopy, S);
  printf ("A = "); pnl_mat_print_nsp (S);
  pnl_mat_lu (S, p);
  pnl_mat_lu_syslin (x, S, p, b);
  printf("x = "); pnl_vect_print_nsp(x); printf("\n");

  printf("test de la fonction 'pnl_mat_syslin' : \n");
  pnl_mat_clone (S, Scopy);
  pnl_mat_syslin (x, S, b);
  printf("x = "); pnl_vect_print_nsp(x); printf("\n");

  printf("test de la fonction 'pnl_mat_syslin_mat' : \n");
  pnl_mat_clone (S, Scopy);
  pnl_mat_clone (B, Bcopy);
  printf ("B = "); pnl_mat_print_nsp (B); printf ("\n");
  pnl_mat_lu (S, p);
  pnl_mat_lu_syslin_mat (S, p, B);
  printf ("X = "); pnl_mat_print_nsp (B); printf ("\n");

  printf("test de la fonction 'pnl_mat_syslin_mat' : \n");
  pnl_mat_clone (S, Scopy);
  pnl_mat_clone (B, Bcopy);
  printf ("B = "); pnl_mat_print_nsp (B); printf ("\n");
  pnl_mat_syslin_mat (S, B);
  printf ("X = "); pnl_mat_print_nsp (B); printf ("\n");

  pnl_mat_free (&B);
  pnl_mat_free (&Bcopy);
  pnl_mat_free (&S);
  pnl_mat_free (&Scopy);
  pnl_vect_free (&b);
  pnl_vect_free (&x);
  pnl_vect_int_free (&p);
}

static void pnl_mat_create_from_file_test ()
{
  PnlMat *M;
  PnlMatComplex *C;
  printf("test de la fonction 'pnl_mat_create_from_file' : \n");
  M = pnl_mat_create_from_file ("Data/Test_read_mat.dat");
  pnl_mat_print (M);
  pnl_mat_free (&M);

  C = pnl_mat_complex_create_from_file ("Data/Test_read_mat.dat");
  pnl_mat_complex_print (C);
  pnl_mat_complex_free (&C);
  
}

static void pnl_mat_dgemm_test()
{
  PnlMat *A, *B, *C;
  PnlVect *Vd;
  double alpha, beta, *d;
  int type_gen;
  type_gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  pnl_rand_init (type_gen, 20, 1);
  Vd = pnl_vect_create (20);
  pnl_vect_rand_normal (Vd, 20, type_gen);
  d = Vd->array;
  printf("test de la fonction dgemm : \n");

  printf("d = "); pnl_vect_print_nsp (Vd); printf ("\n");
  A=pnl_mat_create_from_ptr (4,5,d);
  B=pnl_mat_create_from_ptr (5,4,d);
  C=pnl_mat_create_from_ptr (4, 4, d);
  printf ("A = matrix (d, 5, 4)'; \n");
  printf ("B = matrix (d, 4, 5)'; \n");
  printf ("C = matrix (d(1:16), 4, 4)'; \n");
  alpha = 7.;
  beta = -2.;
  printf ("alpha = %f; beta = %f;\n", alpha, beta);
  pnl_mat_dgemm ('N', 'N', alpha, A, B, beta, C);
  printf("alpha * A * B + beta * C = \n"); 
  pnl_mat_print(C);
  pnl_mat_free(&A);
  pnl_mat_free(&B);
  pnl_mat_free(&C);

  printf ("A = matrix (d, 4, 5)'; \n");
  A=pnl_mat_create_from_ptr (5,4,d);
  B=pnl_mat_create_from_ptr (5,4,d);
  C=pnl_mat_create_from_ptr (4, 4, d);
  pnl_mat_dgemm ('T', 'N', alpha, A, B, beta, C);
  printf("alpha * A' * B + beta * C = \n"); 
  pnl_mat_print(C);
  pnl_mat_free(&A);
  pnl_mat_free(&B);
  pnl_mat_free(&C);

  printf ("A = matrix (d, 5, 4)'; \n");
  printf ("B = matrix (d, 5, 4)'; \n");
  A=pnl_mat_create_from_ptr (4,5,d);
  B=pnl_mat_create_from_ptr (4,5,d);
  C=pnl_mat_create_from_ptr (4, 4, d);
  pnl_mat_dgemm ('N', 'T', alpha, A, B, beta, C);
  printf("alpha * A * B' + beta * C = \n"); 
  pnl_mat_print(C);
  pnl_mat_free(&A);
  pnl_mat_free(&B);
  pnl_mat_free(&C);

  printf ("A = matrix (d, 4, 5)'; \n");
  printf ("B = matrix (d, 5, 4)'; \n");
  A=pnl_mat_create_from_ptr (5, 4,d);
  B=pnl_mat_create_from_ptr (4, 5,d);
  C=pnl_mat_create_from_ptr (4, 4, d);
  pnl_mat_dgemm ('T', 'T', alpha, A, B, beta, C);
  printf("alpha * A' * B' + beta * C = \n"); 
  pnl_mat_print(C);
  pnl_mat_free(&A);
  pnl_mat_free(&B);
  pnl_mat_free(&C);
  pnl_vect_free (&Vd);
}

static void pnl_mat_dgemv_test ()
{
  double   alpha, beta;
  int      type_gen;
  PnlMat  *A;
  PnlVect *x, *y;
  int      N, m, n;
  m = 4; n = 5; N = MAX (m, n);
  alpha = 2.;
  beta = 3.;
  A = pnl_mat_create (m,n);
  x = pnl_vect_create (n);
  y = pnl_vect_create (m);
  type_gen = 6;
  pnl_rand_init (type_gen, N, N);
  pnl_mat_rand_normal (A, m, n, type_gen);
  pnl_vect_rand_normal (x, n, type_gen);
  pnl_vect_rand_normal (y, m, type_gen);

  printf ("test de la fonction pnl_mat_dgemv : \n");
  printf ("A = "); pnl_mat_print_nsp (A); printf ("\n");
  printf ("x = "); pnl_vect_print_nsp (x); printf ("\n");
  printf ("y = "); pnl_vect_print_nsp (y); printf ("\n");
  pnl_mat_dgemv ('N', alpha, A, x, beta, y);
  printf ("%f * A * x + %f * y = ", alpha, beta);
  pnl_vect_print_nsp (y);printf ("\n");

  pnl_mat_rand_normal (A, n, m, type_gen);
  pnl_vect_rand_normal (x, n, type_gen);
  pnl_vect_rand_normal (y, m, type_gen);
  printf ("A = "); pnl_mat_print_nsp (A); printf ("\n");
  printf ("x = "); pnl_vect_print_nsp (x); printf ("\n");
  printf ("y = "); pnl_vect_print_nsp (y); printf ("\n");
  pnl_mat_dgemv ('T', alpha, A, x, beta, y); 
  printf ("%f * A' * x + %f * y = ", alpha, beta); 
  pnl_vect_print_nsp (y);printf ("\n");

  pnl_mat_free (&A);
  pnl_vect_free (&x);
  pnl_vect_free (&y);
}

static void pnl_mat_mult_vect_transpose_test ()
{
  PnlMat *A;
  PnlVect *x, *y;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  printf("test de la fonction pnl_mat_mult_vect_transpose : \n");
  pnl_rand_init (gen, 5, 4);
  A = pnl_mat_create (5, 4);
  printf ("A = ");  pnl_mat_rand_normal (A, 5, 4, gen);
  pnl_mat_print_nsp (A);
  printf("\n");
  x = pnl_vect_create (5);
  pnl_vect_rand_normal(x, 5, gen);
  printf ("x = "); pnl_vect_print_nsp (x);
  printf("\n");

  y = pnl_vect_create (0);
  pnl_mat_mult_vect_transpose_inplace (y, A, x);
  pnl_vect_print_nsp (y);
  printf("\n");

  pnl_vect_free (&x);
  pnl_vect_free (&y);
  pnl_mat_free (&A);
}

static void pnl_mat_scalar_prod_A_test ()
{
  PnlMat *A;
  PnlVect *x, *y;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  printf("test de la fonction pnl_mat_scalar_prod_A : \n");
  pnl_rand_init (gen, 5, 4);
  A = pnl_mat_create (5, 4);
  pnl_mat_rand_normal (A, 5, 4, gen);
  printf ("A = "); pnl_mat_print_nsp (A);
  printf("\n");
  x = pnl_vect_create (4);
  pnl_vect_rand_normal(x, 4, gen);
  printf ("x = "); pnl_vect_print_nsp (x);
  printf("\n");

  y = pnl_vect_create (5);
  pnl_vect_rand_normal(y, 5, gen);
  printf ("y = "); pnl_vect_print_nsp (y);
  printf("\n");
  printf ("y' A x = %f\n", pnl_mat_scalar_prod_A(A, x, y));

  pnl_vect_free (&x);
  pnl_vect_free (&y);
  pnl_mat_free (&A);
}

static void pnl_mat_exp_test ()
{
  PnlMat *A, *B;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  printf("test de la fonction exp : \n");
  pnl_rand_init (gen, 4, 4);
  A = pnl_mat_create (4, 4);
  pnl_mat_rand_normal (A, 4, 4, gen);
  B = pnl_mat_create (0, 0);
  pnl_mat_exp (B, A);
  printf ("Computing B = exp (A) with A =\n");
  pnl_mat_print_nsp (A);
  printf("and B =\n");
  pnl_mat_print (B);
  pnl_mat_free (&B);
  pnl_mat_free (&A);
}


static void all_test ();
static tst_list mat_tests[] =
  {
    MAKE_ENUM(all_test),
    MAKE_ENUM(pnl_mat_set_test),
    MAKE_ENUM(pnl_mat_get_test),
    MAKE_ENUM(pnl_mat_lget_test),
    MAKE_ENUM(pnl_mat_create_from_double_test),
    MAKE_ENUM(pnl_mat_create_from_ptr_test),
    MAKE_ENUM(pnl_mat_submat_test),
    MAKE_ENUM(pnl_mat_set_diag_test),
    MAKE_ENUM(pnl_mat_copy_test),
    MAKE_ENUM(pnl_mat_clone_test),
    MAKE_ENUM(pnl_mat_resize_test),
    MAKE_ENUM(pnl_mat_map_inplace_test),
    MAKE_ENUM(pnl_mat_plus_mat_test),
    MAKE_ENUM(pnl_mat_plus_double_test),
    MAKE_ENUM(pnl_mat_mult_mat_test),
    MAKE_ENUM(pnl_mat_mult_mat_inplace_test),
    MAKE_ENUM(pnl_mat_mult_vect_test),
    MAKE_ENUM(pnl_mat_mult_double_test),
    MAKE_ENUM(pnl_mat_set_double_test),
    MAKE_ENUM(pnl_mat_mult_mat_term_test),
    MAKE_ENUM(pnl_mat_chol_test),
    MAKE_ENUM(pnl_mat_sq_transpose_test),
    MAKE_ENUM(pnl_mat_transpose_test),
    MAKE_ENUM(pnl_vect_wrap_mat_row_test),
    MAKE_ENUM(pnl_mat_get_row_test),
    MAKE_ENUM(pnl_mat_create_diag_test),
    MAKE_ENUM(pnl_mat_sum_test),
    MAKE_ENUM(pnl_mat_sum_vect_test),
    MAKE_ENUM(pnl_mat_prod_test),
    MAKE_ENUM(pnl_mat_prod_vect_test),
    MAKE_ENUM(pnl_mat_minmax_test),
    MAKE_ENUM(pnl_mat_qsort_test),
    MAKE_ENUM(pnl_mat_map_test),
    MAKE_ENUM(pnl_mat_set_row_test),
    MAKE_ENUM(pnl_mat_triangular_inverse_test),
    MAKE_ENUM(pnl_mat_inverse_test),
    MAKE_ENUM(pnl_mat_syslin_test),
    MAKE_ENUM(pnl_mat_create_from_file_test),
    MAKE_ENUM(pnl_mat_dgemm_test),
    MAKE_ENUM(pnl_mat_exp_test),
    MAKE_ENUM(pnl_mat_mult_vect_transpose_test),
    MAKE_ENUM(pnl_mat_scalar_prod_A_test),
    MAKE_ENUM(pnl_mat_dgemv_test),
    MAKE_ENUM(NULL)
  };

static void all_test ()
{
  run_all_test (mat_tests);
}

void matrix_test()
{
  menu_test (mat_tests);
}
