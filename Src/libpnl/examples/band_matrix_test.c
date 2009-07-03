
/*************************************************************************/
/* Written and (C) by David Pommier <david.pommier@gmail.com>            */
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
#include "pnl_matrix.h"
#include"pnl_perm.h"
#include "pnl_band_matrix.h"
#include "pnl_spec_matrix.h"
#include "pnl_random.h"
#include "tests.h"


void Test_Band_Matrix_algebra(void )
{  
  /*
    Short program to test Band matrix ...
  */
  PnlVect *b,*x,*v;
  int i,j,size=100;
  int band=10;//size/3;
  int gen = 7;
  PnlPermutation *p;
  PnlBandMatrix * M;
  PnlMat * FM;
  pnl_rand_init (gen, size, size);
  b = pnl_vect_create (size);
  v = pnl_vect_create (size);
  FM = pnl_mat_create (size, size);
  x = pnl_vect_create(0);
  pnl_mat_set_id (FM);
  pnl_vect_rand_normal (v, size, gen);
  pnl_vect_rand_normal (b, size, gen);
  pnl_mat_dger (1., v, v, FM);

  M=pnl_band_matrix_create_from_full(FM,band);
  p = pnl_permutation_create (size);
  pnl_permutation_init (p);
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
      if(abs(i-j)>band)
        MLET(FM,i,j)=0;
  pnl_mat_lu (FM, p);
  pnl_mat_lu_syslin (x, FM, p, b);
  // x solution of full system  
  pnl_band_matrix_lu(M,0.00000001);
  pnl_band_matrix_solve_syslin_inplace(M,b);
  // b solution of band matrix system
  pnl_vect_axpby(-1,b,1,x);
  printf("error between LU Sys Lin for Full matrix and Band matrix is %7.4f \n",pnl_vect_norm_two(x));  
 
  pnl_mat_free(&FM);
  pnl_band_matrix_free(&M);
  pnl_vect_free(&x);
  pnl_vect_free(&v);
  pnl_vect_free(&b);
  pnl_permutation_free (&p);
}


void band_matrix_test()
{
  Test_Band_Matrix_algebra ();
}
