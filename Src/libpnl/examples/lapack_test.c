
/*************************************************************************/
/* Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>            */
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


static void pnl_mat_ls_test ()
{
  PnlMat *A, *B;
  PnlVect *b;
  int m, n, gen, nrhs;
  m = 5; n = 6; nrhs = 3;
  gen = PNL_RNG_MERSENNE_RANDOM_SEED;
  printf("test de la fonction eigen : \n");
  A = pnl_mat_create (m, n);
  B = pnl_mat_create (m, nrhs);
  b = pnl_vect_create (m);
  pnl_rand_init (gen, m, n);
  pnl_mat_rand_uni2 (A, m, n, 0, 1, gen);
  pnl_mat_rand_uni2 (B, m, nrhs, 0, 1, gen);
  pnl_vect_rand_uni (b, m, 0, 1, gen);

  printf ("A = "); pnl_mat_print_nsp(A);
  printf ("B = "); pnl_mat_print_nsp(B);
  if (pnl_mat_ls_mat (A, B) != OK) printf("error in pnl_mat_ls\n");
  printf ("X = "); pnl_mat_print_nsp(B);


  printf ("b = "); pnl_vect_print_nsp(b);
  if (pnl_mat_ls (A, b) != OK) printf("error in pnl_mat_ls\n");
  printf ("x = "); pnl_vect_print_nsp(b);

  pnl_mat_free (&A);
  pnl_mat_free (&B);
  pnl_vect_free (&b);
}

static void all_lapack_test ();
static list lapack_tests[] =
  {
    MAKE_ENUM(1, all_lapack_test),
    MAKE_ENUM(2, pnl_mat_eigen_test),
    MAKE_ENUM(3, pnl_mat_log_test),
    MAKE_ENUM(4, pnl_mat_ls_test),
    MAKE_ENUM(NULL_INT, NULL)
  };

static void all_lapack_test ()
{
  int len=0;
  while (lapack_tests[len].id != NULL_INT)
    {
      if (lapack_tests[len].func != all_lapack_test) (lapack_tests[len].func)();
      len ++;
    }
}

void lapack_test()
{
  int len=0, n=0, choice;

  while (lapack_tests[len].id != NULL_INT) len++;
        
    
  while (lapack_tests[n].id != NULL_INT)
    {
      printf("%2d. %s\n",  lapack_tests[n].id, lapack_tests[n].label);
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
              else { (lapack_tests[choice].func)(); return; }
            }
        }
    }
  printf("Which test do you want to run?\n");
  while(1)
    {
      scanf("%d", &choice); choice --;
      if (choice <0 || choice > len) printf("illegal choice\n");
      else { (lapack_tests[choice].func)(); break; }
    }
}
