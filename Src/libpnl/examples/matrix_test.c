
/*************************************************************************/
/* Written and (C) by J�r�me Lelong <jerome.lelong@gmail.com>            */
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "config.h"
#include "pnl_matrix.h"
#include "pnl_matrix_complex.h"
#include "pnl_random.h"
#include "tests.h"

/**
 * Creates a diagonalizable matrix of size n
 * @param n an integer
 * @return a matrix
 */
static PnlMat* pnl_mat_create_diagonalizable (int n)
{
  PnlMat *A;
  PnlMat *P, *invP;
  int i, j;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;

  pnl_rand_init (gen, n, n);
  A = pnl_mat_create (n, n);
  P = pnl_mat_create (n, n);
  invP = pnl_mat_create (n, n);
  pnl_mat_rand_normal (A, n, n, gen);

  /* P = exp (rand) is always invertible */
  pnl_mat_exp(P, A);
  pnl_mat_inverse(invP, P);

  /* P = P * D */
  for ( j=0 ; j<n ; j++ )
    {
      double dj = fabs(pnl_rand_uni (gen));
      for ( i=0 ; i<n ; i++ )
        {
          MLET (P, i, j) *= dj;
        }
    }
  pnl_mat_mult_mat_inplace (A, P, invP);
  pnl_mat_free (&P);
  pnl_mat_free (&invP);
  return A;
}

/* static double function_prod(double x, double y) {return x*y;} */

/*TEST DES FONCTIONS PREMIAMAT */


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

static void pnl_mat_wrap_row_test()
{
  PnlMat *M;
  PnlVect V;
  double x[8]={1.0, 5.0, 3.0, 8.0, 3.0, 7.0, 6.0, 9.0};
  printf("test de la fonction 'pnl_mat_wrap_row : \n");
  M=pnl_mat_create_from_ptr(4,2,x);
  V=pnl_mat_wrap_row(M,2);
  pnl_vect_print(&V);
  pnl_mat_free(&M);
}



static void pnl_mat_row_to_vect_inplace_test()
{
  PnlMat *M;
  PnlVect *V;
  double x[8]={1.0, 5.0, 3.0, 8.0, 3.0, 7.0, 6.0, 9.0};
  printf("test de la fonction 'pnl_mat_row_to_vect_inplace : \n");
  M=pnl_mat_create_from_ptr(4,2,x);
  V=pnl_vect_create(0);
  pnl_mat_row_to_vect_inplace(V,M,2);
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

static void pnl_mat_inverse_test ()
{
  PnlMat *A, *B;
  double x[16] = {1.0, 2.0, 3.0, 4.0, 0.0, 5.0, 6.0, 7.0, 0.0, 0.0, 8.0, 9.0, 0.0, 0.0, 0.0, 10.0};

  B = pnl_mat_create_from_ptr(4,4, x);
  A=pnl_mat_create(0,0);
  printf("test de la fonction 'pnl_mat_upper_inverse' : \n");
  pnl_mat_upper_inverse(A,B);
  pnl_mat_print(A);
  printf("test de la fonction 'pnl_mat_lower_inverse' : \n");
  pnl_mat_sq_transpose (B);
  pnl_mat_lower_inverse(A,B);
  pnl_mat_print(A);
  pnl_mat_free(&A);
  pnl_mat_free(&B);
}

static void pnl_mat_syslin_test ()
{
  PnlMat *S, *Scopy, *B, *Bcopy;
  PnlVect *v, *b, *x;
  PnlPermutation *p;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  printf("test de la fonction 'pnl_mat_chol' : \n");
  pnl_rand_init (gen, 5, 5);
  b = pnl_vect_create (5);
  v = pnl_vect_create (5);
  S = pnl_mat_create (5, 5);
  pnl_mat_set_id (S);
  pnl_vect_rand_normal (v, 5, gen);
  pnl_vect_rand_normal (b, 5, gen);
  pnl_mat_dger (1., v, v, S);
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
  p = pnl_permutation_create (5);
  pnl_permutation_init (p);
  pnl_mat_lu (S, p);
  pnl_mat_lu_syslin (x, S, p, b);
  pnl_vect_print(x);
  pnl_permutation_free (&p);

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
  pnl_mat_chol_syslin_mat (S, B);
  printf ("X = "); pnl_mat_print_nsp (B); printf ("\n");

  printf("test de la fonction 'pnl_mat_syslin_mat' (symmetric matrix): \n");
  pnl_mat_clone (S, Scopy);
  pnl_mat_clone (B, Bcopy);
  printf ("B = "); pnl_mat_print_nsp (B); printf ("\n");
  pnl_mat_syslin_mat (S, B);
  printf ("X = "); pnl_mat_print_nsp (B); printf ("\n");


  printf("test de la fonction 'pnl_mat_lu_syslin' : \n");
  pnl_mat_rand_normal (S, 5, 5, gen);
  pnl_mat_exp (Scopy, S);
  pnl_mat_clone (S, Scopy);
  printf ("A = "); pnl_mat_print_nsp (S);
  p = pnl_permutation_create (5);
  pnl_permutation_init (p);
  pnl_mat_lu (S, p);
  pnl_mat_lu_syslin (x, S, p, b);
  printf("x = "); pnl_vect_print_nsp(x); printf("\n");
  pnl_permutation_free (&p);

  printf("test de la fonction 'pnl_mat_syslin' : \n");
  pnl_mat_clone (S, Scopy);
  pnl_mat_syslin (x, S, b);
  printf("x = "); pnl_vect_print_nsp(x); printf("\n");

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
  pnl_vect_free (&v);
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



static void pnl_mat_eigen_test ()
{
  PnlMat *A, *B, *P;
  PnlVect *v;
  int n = 4;
  int gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  printf("test de la fonction eigen : \n");
  pnl_rand_init (gen, n, n);
  A = pnl_mat_create (n, n);
  B = pnl_mat_create (n, n);
  P = pnl_mat_create (n, n);
  pnl_mat_rand_uni2 (A, n, n, 0., 1., gen);
  pnl_mat_dgemm ('N', 'T', 1., A, A, 0., B);
  v = pnl_vect_create (0);
  printf ("S = "); pnl_mat_print_nsp (B);
  pnl_mat_eigen (v, P, B, TRUE);
  printf ("P = "); pnl_mat_print_nsp (P);
  printf ("v = "); pnl_vect_print_nsp (v);
  pnl_mat_free(&A);
  pnl_mat_free (&B);
  
  A = pnl_mat_create_diagonalizable (n);
  printf ("A = "); pnl_mat_print_nsp (A);
  pnl_mat_eigen (v, P, A, TRUE);
  printf ("P = "); pnl_mat_print_nsp (P);
  printf ("v = "); pnl_vect_print_nsp (v);
  pnl_vect_free (&v);
  pnl_mat_free (&A);
  pnl_mat_free (&P);
}

static void pnl_mat_log_test ()
{
  PnlMat *A, *B;
  int n = 5;
  printf("test de la fonction log : \n");
  A = pnl_mat_create_diagonalizable (n);
  B = pnl_mat_create (n,n);
  printf ("A =");
  pnl_mat_print_nsp (A);
  pnl_mat_log (B, A);
  printf ("A =");
  pnl_mat_print_nsp (A);
  printf("logA =");
  pnl_mat_print_nsp (B);
  printf("\n");
  pnl_mat_free (&B);
  pnl_mat_free (&A);
}

static void all_matrix_test ();
static list mat_tests[] =
  {
    MAKE_ENUM(1, all_matrix_test),
    MAKE_ENUM(2, pnl_mat_set_test),
    MAKE_ENUM(3, pnl_mat_get_test),
    MAKE_ENUM(4, pnl_mat_lget_test),
    MAKE_ENUM(5, pnl_mat_create_from_double_test),
    MAKE_ENUM(6, pnl_mat_create_from_ptr_test),
    MAKE_ENUM(7, pnl_mat_copy_test),
    MAKE_ENUM(8, pnl_mat_clone_test),
    MAKE_ENUM(9, pnl_mat_map_inplace_test),
    MAKE_ENUM(10, pnl_mat_plus_mat_test),
    MAKE_ENUM(11, pnl_mat_plus_double_test),
    MAKE_ENUM(12, pnl_mat_mult_mat_test),
    MAKE_ENUM(13, pnl_mat_mult_mat_inplace_test),
    MAKE_ENUM(14, pnl_mat_mult_vect_test),
    MAKE_ENUM(15, pnl_mat_mult_double_test),
    MAKE_ENUM(16, pnl_mat_set_double_test),
    MAKE_ENUM(17, pnl_mat_mult_mat_term_test),
    MAKE_ENUM(18, pnl_mat_chol_test),
    MAKE_ENUM(19, pnl_mat_sq_transpose_test),
    MAKE_ENUM(20, pnl_mat_transpose_test),
    MAKE_ENUM(21, pnl_mat_wrap_row_test),
    MAKE_ENUM(22, pnl_mat_row_to_vect_inplace_test),
    MAKE_ENUM(23, pnl_mat_create_diag_test),
    MAKE_ENUM(24, pnl_mat_sum_test),
    MAKE_ENUM(25, pnl_mat_sum_vect_test),
    MAKE_ENUM(26, pnl_mat_prod_test),
    MAKE_ENUM(27, pnl_mat_prod_vect_test),
    MAKE_ENUM(28, pnl_mat_minmax_test),
    MAKE_ENUM(29, pnl_mat_qsort_test),
    MAKE_ENUM(30, pnl_mat_map_test),
    MAKE_ENUM(31, pnl_mat_set_row_test),
    MAKE_ENUM(32, pnl_mat_inverse_test),
    MAKE_ENUM(33, pnl_mat_syslin_test),
    MAKE_ENUM(34, pnl_mat_create_from_file_test),
    MAKE_ENUM(35, pnl_mat_dgemm_test),
    MAKE_ENUM(36, pnl_mat_exp_test),
    MAKE_ENUM(37, pnl_mat_mult_vect_transpose_test),
    MAKE_ENUM(38, pnl_mat_scalar_prod_A_test),
    MAKE_ENUM(39, pnl_mat_dgemv_test),
    MAKE_ENUM(40, pnl_mat_eigen_test),
    MAKE_ENUM(41, pnl_mat_log_test),
    MAKE_ENUM(NULL_INT, NULL)
  };

static void all_matrix_test ()
{
  int len=0;
  while (mat_tests[len].id != NULL_INT)
    {
      if (mat_tests[len].func != all_matrix_test) (mat_tests[len].func)();
      len ++;
    }
}

void matrix_test()
{
  int len=0, n=0, choice;

  while (mat_tests[len].id != NULL_INT) len++;
        
    
  while (mat_tests[n].id != NULL_INT)
    {
      printf("%2d. %s\n",  mat_tests[n].id, mat_tests[n].label);
      n ++;
      if (n/10 == (double)n/10.0)
        {
          printf("Which test do you want to run? (type 0 to continue displaying the list)\n");
          while(1)
            {
              scanf("%d", &choice);
              if (choice ==0) break;
              choice --;
              if (choice <0 || choice > len) printf("illegal choice\n");
              else { (mat_tests[choice].func)(); return; }
            }
        }
    }
  printf("Which test do you want to run?\n");
  while(1)
    {
      scanf("%d", &choice); choice --;
      if (choice <0 || choice > len) printf("illegal choice\n");
      else { (mat_tests[choice].func)(); break; }
    }
}
