
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
#include <string.h>
#include <math.h>
#include <time.h>

#include "pnl/pnl_machine.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_random.h"

extern int C2F(dgemm) (char *transa, char *transb, int *m, int *n, int *k, double *alpha,
                       double *a, int *lda, double *b, int *ldb, double *beta,
                       double *c__, int *ldc);

extern int C2F(dgemv) (char *transa, int *m, int *n, double *alpha,
                       double *a, int *lda, double *x, int *incx, double *beta,
                       double *y, int *incy);
extern int C2F(dgesv) (int *n, int *nrhs, double *A, int *lda, int *ipvi, double *B,
                       int *ldb, int *info);
extern int C2F(dposv) (char *uplo, int *n, int *nrhs, double *A, int *lda, double *B,
                       int *ldb, int *info);
extern int C2F(dgetrf) (int *m, int *n, double *A, int *lda, int *ipvi, int *info);
extern int C2F(dgetf2) (int *m, int *n, double *A, int *lda, int *ipvi, int *info);
extern int C2F(dpotrf) (char *uplo, int *n, double *a, int *lda, int *info);


static void speed_access ()
{
  clock_t  start, end;
  double   cpu_time_used;
  int      type_gen, i,j;
  double   sum;
  double  *ptr;
  int      N   = 10000;
  int      dim = 500;
  PnlMat  *M   = pnl_mat_create (0,0);
  type_gen = 6;
  pnl_rand_init (type_gen, N, dim);
  pnl_mat_rand_normal (M, N, dim, type_gen);

  printf ("\n");     
  start = clock();
  sum = 0.0;
  for (i=0; i<N; i++)
    {
      for (j=0; j<dim; j++){
        sum += MGET(M,i,j);
      }
    }
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("%f get version : %f\n", sum, cpu_time_used);
  start = clock();
  sum = 0.0;
  for (i=0, ptr=M->array; i<N*dim; i++, ptr++)
    {
      sum += *ptr;
    }
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("%f ptr version : %f\n", sum, cpu_time_used);
  pnl_mat_free (&M);
}

void speed_copy ()
{
  clock_t  start, end;
  double   cpu_time_used;
  int      type_gen, i;
  PnlMat  *copy;
  double  *d;
  int      N   = 10000;
  int      dim = 500;
  PnlMat  *M   = pnl_mat_create (0,0);
  type_gen = 6;
  pnl_rand_init (type_gen, N, dim);
  pnl_mat_rand_normal (M, N, dim, type_gen);
  printf ("\n");
     
  start = clock();
  copy = pnl_mat_copy (M);
  end = clock();
  pnl_mat_free (&copy);
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("pnl_mat_copy time : %f\n", cpu_time_used);
  start = clock();
  d = malloc (sizeof(double)*N*dim);
  for (i=0; i<N*dim; i++)
    {
      d[i] = M->array[i];
    }
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("brute copy time : %f\n", cpu_time_used);
  pnl_mat_free (&M);
  free (d);
}

void speed_setvalue ()
{
  clock_t  start, end;
  double   cpu_time_used;
  int      n   = 10000;
  int      m = 500;
  PnlMat  *M1   = pnl_mat_create (n, m);
  PnlMat  *M2   = pnl_mat_create (n, m);
     
  start = clock();
  memset (M1->array, 0, n*m*sizeof(double));
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("memset time : %f\n", cpu_time_used);
  start = clock();
  pnl_mat_set_zero (M2);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("pnl_mat_set_zero time : %f\n", cpu_time_used);
  pnl_mat_free (&M1);
  pnl_mat_free (&M2);
}


void dgemm3 (double alpha, const PnlMat *A,
             const PnlMat *B, double beta, PnlMat *C)
{
  int m, n;
  int i, j, k;
  double Cij, Bkj, aibk;
  
  m = A->m; n = B->n;
  if ( beta == 0. )
    {
      pnl_mat_resize (C, m, n);
      pnl_mat_set_all (C, 0.);
    }
  else if ( beta != 1. )
    {
      pnl_mat_mult_scalar (C, beta); 
    }
  if ( alpha == 0.) return;


#ifndef PNL_RANGE_CHECK_OFF
  if (A->n != B->m || (A->m != C->m || B->n != C->n) )
    {
      PNL_ERROR ("size mismatch", "pnl_mat_dgemm");
    }
#endif

  for (k=0; k<A->n; k++) 
    {
      for (i=0; i<C->m; i++)
        {
          double temp = alpha * PNL_MGET (A, i, k);
          for (j=0; j<C->n; j++)
            {
              Bkj = PNL_MGET (B, k, j);
              
              Cij = PNL_MGET(C, i, j);
              aibk = temp * Bkj;
              Cij += aibk;
              PNL_MLET (C, i, j) = Cij;
              
            }
        }
    }
}

void dgemm2 (double alpha, const PnlMat *A,
             const PnlMat *B, double beta, PnlMat *C)
{
  int m, n;
  int i, j, k;
  double Cij, Bkj, aibk;
  
  m = A->m; n = B->n;
  if ( beta == 0. )
    {
      pnl_mat_resize (C, m, n);
      pnl_mat_set_all (C, 0.);
    }
  else if ( beta != 1. )
    {
      pnl_mat_mult_scalar (C, beta); 
    }
  if ( alpha == 0.) return;


#ifndef PNL_RANGE_CHECK_OFF
  if (A->n != B->m || (A->m != C->m || B->n != C->n) )
    {
      PNL_ERROR ("size mismatch", "pnl_mat_dgemm");
    }
#endif

  for (i=0; i<C->m; i++)
    {
      for (k=0; k<A->n; k++) 
        {
          double temp = alpha * PNL_MGET (A, i, k);
          for (j=0; j<C->n; j++)
            {
              Bkj = PNL_MGET (B, k, j);
              
              Cij = PNL_MGET(C, i, j);
              aibk = temp * Bkj;
              Cij += aibk;
              PNL_MLET (C, i, j) = Cij;
              
            }
        }
    }
}

void speed_dgemm ()
{
  clock_t  start, end;
  double   cpu_time_used;
  double alpha, beta;
  int      type_gen;
  PnlMat  *A, *B, *C , *CC;
  int      N   = 1200;
  alpha = 2.; beta = 3.;
  A = pnl_mat_create (N,N);
  B = pnl_mat_create (N,N);
  C = pnl_mat_create (N,N);
  type_gen = 0;
  pnl_rand_init (type_gen, N, N);
  pnl_mat_rand_normal (A, N, N, type_gen);
  pnl_mat_rand_normal (B, N, N, type_gen);
  pnl_mat_rand_normal (C, N, N, type_gen);
  CC = pnl_mat_copy (C);
  
  printf ("\n");

  pnl_mat_clone (C, CC);
  start = clock();
  dgemm_ ("N", "N", &A->m,&A->m,&A->m, &alpha, A->array, &A->m, B->array, &B->m,&beta, C->array,&C->m);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("BLAS dgemm time NN : %f\n", cpu_time_used);
  
  pnl_mat_clone (C, CC);
  start = clock();
  pnl_mat_dgemm ('N', 'N', alpha, A, B, beta, C);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("pnl_mat_dgemm NN time : %f\n", cpu_time_used);


  pnl_mat_clone (C, CC);
  start = clock();
  dgemm2 (alpha, A, B, beta, C);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("dgemm2 time : %f\n", cpu_time_used);

  pnl_mat_clone (C, CC);
  start = clock();
  dgemm3 (alpha, A, B, beta, C);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("dgemm3 time : %f\n", cpu_time_used);

#if 0
  pnl_mat_clone (C, CC);
  start = clock();
  pnl_mat_dgemm ('N', 'T', alpha, A, B, beta, C);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("pnl_mat_dgemm NT time : %f\n", cpu_time_used);

  pnl_mat_clone (C, CC);
  start = clock();
  pnl_mat_dgemm ('T', 'N', alpha, A, B, beta, C);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("pnl_mat_dgemm TN time : %f\n", cpu_time_used);

  pnl_mat_clone (C, CC);
  start = clock();
  pnl_mat_dgemm ('T', 'T', alpha, A, B, beta, C);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("pnl_mat_dgemm TT time : %f\n", cpu_time_used);
#endif

  pnl_mat_free (&A);
  pnl_mat_free (&B);
  pnl_mat_free (&C);
  pnl_mat_free (&CC);
}

void speed_dgemv ()
{
  clock_t  start, end;
  double   cpu_time_used;
  double   alpha, beta;
  int      type_gen, inc;
  PnlMat  *A, *tA;
  PnlVect *x, *y, *yy;
  int      N   = 1000;
  alpha = 1.;
  beta = 2.;
  inc = 1;
  A = pnl_mat_create (N,N);
  x = pnl_vect_create (N);
  y = pnl_vect_create (N);
  type_gen = 6;
  pnl_rand_init (type_gen, N, N);
  pnl_mat_rand_normal (A, N, N, type_gen);
  tA = pnl_mat_copy (A);
  pnl_mat_sq_transpose (tA);
  pnl_vect_rand_normal (x, N, type_gen);
  pnl_vect_rand_normal (y, N, type_gen);
  yy = pnl_vect_copy (y);
  
  printf ("\n");
  start = clock();
  pnl_mat_dgemv ('N', alpha, A, x, beta, y);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("pnl_mat_dgemv 'N' time : %f\n", cpu_time_used);

  pnl_vect_clone (y, yy);
  start = clock();
  dgemv_ ("N", &A->m,&A->m, &alpha, tA->array, &A->m, x->array, &inc, &beta, y->array, &inc);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("BLAS dgemv time 'N': %f\n", cpu_time_used);

  pnl_vect_clone (y, yy);
  start = clock();
  pnl_mat_dgemv ('t', alpha, A, x, beta, y);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("pnl_mat_dgemv 'T' time : %f\n", cpu_time_used);

  pnl_vect_clone (y, yy);
  start = clock();
  dgemv_ ("T", &A->m,&A->m, &alpha, tA->array, &A->m, x->array, &inc, &beta, y->array, &inc);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("BLAS dgemv time 'T': %f\n", cpu_time_used);

  pnl_mat_free (&A);
  pnl_mat_free (&tA);
  pnl_vect_free (&x);
  pnl_vect_free (&y);
  pnl_vect_free (&yy);
}

/**
 * nsp_dsumrows:
 * @x: an array m x n of double with column major order
 * @res: an array of size m with the result
 * @m, @n: sizes of the matrix @x
 * 
 *   res[i] = sum(x[i,j], j=0,n-1)
 * 
 * computes the sum of each row
 *
 * Written by Jean-Philippe Chancelier for Nsp. 
 **/
static void nsp_dsumrows(double *x, double *res, int m, int n)
{
  int i, j, k, k2, k3, p;
  p = n % 3;

  for ( i = 0 ; i < m ; i++ ) 
    res[i] = 0.0;

  k = 0;
  for ( j = 0 ; j < p ; j++ )
    for ( i = 0 ; i < m ; i++, k++)
      res[i] += x[k];

  k2 = k + m;
  k3 = k2 + m;
  for ( j = p ; j < n ; j+=3, k+=2*m, k2+=2*m, k3+=2*m )
    for ( i = 0 ; i < m ; i++, k++, k2++, k3++)
      res[i] += x[k] + x[k2] + x[k3];
}


static void speed_mat_sum_vect ()
{
  PnlMat *A;
  PnlVect *x;
  clock_t  start, end;
  double   cpu_time_used;
  int gen = 7, N=1000;
  pnl_rand_init (gen, N, N);
  A = pnl_mat_create (N, N);
  x = pnl_vect_create (N);
  pnl_mat_rand_normal (A, N, N, gen);
  printf ("\n");
  start = clock();
  pnl_mat_sum_vect (x, A, 'r');
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("pnl_mat_sum_vect 'r' : %f\n", cpu_time_used);

  start = clock();
  pnl_mat_sum_vect (x, A, 'c');
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("pnl_mat_sum_vect 'c' : %f\n", cpu_time_used);

  start = clock();
  nsp_dsumrows (A->array, x->array, A->m, A->n);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("nsp_sum 'c' : %f\n", cpu_time_used);

  pnl_mat_free (&A);
  pnl_vect_free (&x);
}


static void speed_syslin_mat ()
{
  PnlMat *S, *Scopy, *B, *Bcopy;
  PnlVect *v, *vcopy;
  clock_t  start, end;
  double   cpu_time_used;
  int *ipiv, info, nrhs;
  PnlPermutation *p;
  int gen = 7, N=1000;
  pnl_rand_init (gen, N, N);
  v = pnl_vect_create (N);
  S = pnl_mat_create (N, N);
  pnl_mat_set_id (S);
  pnl_vect_rand_normal (v, N, gen);
  pnl_mat_dger (5., v, v, S);
  Scopy = pnl_mat_copy (S);
  vcopy = pnl_vect_copy (v);

  B = pnl_mat_create (N,N);
  pnl_mat_rand_normal (B, N, N, gen);
  Bcopy = pnl_mat_copy (B);
  printf ("\n");
  start = clock();
  pnl_mat_chol_syslin_mat (S, B);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("pnl_mat_chol_syslin_mat time : %f\n", cpu_time_used);

  pnl_mat_clone (S, Scopy);
  pnl_mat_clone (B, Bcopy);
  start = clock();
  pnl_mat_syslin_mat (S, B);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("pnl_mat_syslin_mat time : %f\n", cpu_time_used);

  pnl_mat_clone (S, Scopy);
  pnl_mat_clone (B, Bcopy);
  start = clock();
  ipiv = malloc (sizeof(int) * N);
  dgesv_ (&N, &N, S->array, &N, ipiv, B->array, &N, &info);
  free (ipiv);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("BLAS dgesv mat time : %f \n", cpu_time_used);

  pnl_mat_clone (S, Scopy);
  pnl_mat_clone (B, Bcopy);
  start = clock();
  dposv_ ("L", &N, &N, S->array, &N, B->array, &N, &info);
  if (info != 0) {
    printf ("error in dposv\n"); abort();
  }
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("BLAS dposv mat time : %f \n", cpu_time_used);

  pnl_vect_clone (v, vcopy);
  pnl_mat_clone (S, Scopy);
  start = clock();
  pnl_mat_syslin_inplace (S,v);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("pnl_mat_syslin time : %f\n", cpu_time_used);


  pnl_vect_clone (v, vcopy);
  pnl_mat_clone (S, Scopy);
  start = clock();
  pnl_mat_chol (S);
  pnl_mat_chol_syslin_inplace (S,v);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("pnl_mat_chol_syslin time : %f\n", cpu_time_used);

  nrhs = 1;
  pnl_vect_clone (v, vcopy);
  pnl_mat_clone (S, Scopy);
  start = clock();
  ipiv = malloc (sizeof(int) * N);
  dgesv_ (&N, &nrhs, S->array, &N, ipiv, v->array, &N, &info);
  free (ipiv);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("BLAS dgesv time : %f \n", cpu_time_used);

  nrhs = 1;
  pnl_vect_clone (v, vcopy);
  pnl_mat_clone (S, Scopy);
  start = clock();
  dposv_ ("U", &N, &nrhs, S->array, &N, v->array, &N, &info);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("BLAS dposv time : %f \n", cpu_time_used);
  
  pnl_mat_clone (S, Scopy);
  start = clock();
  p = pnl_permutation_create (N);
  pnl_mat_lu (S,p);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  pnl_permutation_free (&p);
  printf("pnl_mat_lu time : %f\n", cpu_time_used);

  pnl_mat_clone (S, Scopy);
  start = clock();
  pnl_mat_chol (S);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("pnl_mat_chol time : %f\n", cpu_time_used);

  nrhs = 1;
  pnl_mat_clone (S, Scopy);
  start = clock();
  ipiv = malloc (sizeof(int) * N);
  dgetrf_ (&N, &N, S->array, &N, ipiv, &info);
  free (ipiv);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("BLAS LU dgetrf time : %f \n", cpu_time_used);

  nrhs = 1;
  pnl_mat_clone (S, Scopy);
  start = clock();
  dpotrf_ ("U", &N, S->array, &N, &info);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("BLAS Chol dpotrf time : %f \n", cpu_time_used);

  pnl_mat_free (&B);
  pnl_mat_free (&Bcopy);
  pnl_mat_free (&S);
  pnl_mat_free (&Scopy);
  pnl_vect_free (&v);
  pnl_vect_free (&vcopy);

}

void speed_matrix_test ()
{
  speed_access ();
  speed_copy ();
  speed_setvalue ();
  speed_dgemm ();
  speed_dgemv ();
  speed_mat_sum_vect ();
  speed_syslin_mat ();
}
