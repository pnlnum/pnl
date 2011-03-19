
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
#include "pnl/pnl_linalgsolver.h"
#include "tests.h"


/*TEST DES FONCTIONS PREMIAVECT// */
void Test_Solver_sym(void )
{  
  /*
    Short program to test solver CG & shown how to use it...
  */

  int Size;
  PnlVect *x1,*x2,*x3, *b;
  PnlMat *M,*PC;
  PnlCgSolver* Solver;
  PnlBicgSolver* Solver2;
  PnlGmresSolver* Solver3;
  printf(">>>>>> Test iterativ Solver  Symetric matrix \n");
  Size=10;
  b=pnl_vect_create_from_double(Size,1);
  b->array[2]=2.0;
  b->array[1]=0.0;
  x1=pnl_vect_create_from_double(Size,2);
  x2=pnl_vect_create_from_double(Size,2);
  x3=pnl_vect_create_from_double(Size,2);
  M=pnl_mat_create_from_file("Data/Test_mat_sym");
  PC=pnl_mat_create_from_file("Data/Test_PCmat");

  Solver=pnl_cg_solver_create(b->size,20,1e-6);
  Solver2=pnl_bicg_solver_create(b->size,20,1e-6);
  Solver3=pnl_gmres_solver_create(b->size,10,10,1e-6);
  pnl_mat_cg_solver_solve(M,PC,x1,b,Solver);
  pnl_mat_bicg_solver_solve(M,PC,x2,b,Solver2);
  pnl_mat_gmres_solver_solve(M,PC,x3,b,Solver3);
  pnl_vect_axpby(-1.0, x1, 1.0, x2);
  pnl_vect_axpby(-1.0, x1, 1.0, x3);
  printf("==> Error BICG-CG %e \n",pnl_vect_norm_two(x2));
  printf("==> Error GMRES-CG %e \n",pnl_vect_norm_two(x3));
  pnl_cg_solver_free(&Solver);
  pnl_bicg_solver_free(&Solver2);
  pnl_gmres_solver_free(&Solver3);
  pnl_mat_free(&M);
  pnl_mat_free(&PC);
  pnl_vect_free(&x1);
  pnl_vect_free(&x2);
  pnl_vect_free(&x3);
  pnl_vect_free(&b);
  
}
void Test_Solver_no_sym(void )
{  
  /*
    Short program to test solver CG & shown how to use it...
  */
  int Size;
  PnlVect *b,*res, *x1,*x2;
  PnlMat *Q,*PC;
  PnlBicgSolver* Solver2;
  PnlGmresSolver* Solver3;
  printf(">>>>>> Test iterativ Solver Non Symetric matrix \n");
  Size=20;
  b   = pnl_vect_create_from_file ("Data/Test_vect_rhs.dat");
  res = pnl_vect_create_from_file ("Data/Test_res.dat");
  x1=pnl_vect_create_from_double(Size,0.0);
  x2=pnl_vect_create_from_double(Size,0.0);
  Q=pnl_mat_create_from_file("Data/Test_mat_no_sym.dat");
  PC=pnl_mat_create_from_file("Data/Test_PCmat_20.dat");

  Solver2=pnl_bicg_solver_create(b->size,100,1e-6);
  Solver3=pnl_gmres_solver_create(b->size,100,20,1e-6);
  pnl_mat_bicg_solver_solve(Q,PC,x1,b,Solver2);
  pnl_mat_gmres_solver_solve(Q,PC,x2,b,Solver3);
  pnl_vect_axpby(-1.0, res, 1.0, x1);
  pnl_vect_axpby(-1.0, res, 1.0, x2);
  printf("==> Error BICG %g \n",pnl_vect_norm_two(x1));
  printf("==> Error GMRES %g \n",pnl_vect_norm_two(x2));
  pnl_bicg_solver_free(&Solver2);
  pnl_gmres_solver_free(&Solver3);
  pnl_mat_free(&Q);
  pnl_mat_free(&PC);
  pnl_vect_free(&x1);
  pnl_vect_free(&x2);
  pnl_vect_free(&res);
  pnl_vect_free(&b);
}


int main ()
{
  Test_Solver_sym ();
  Test_Solver_no_sym ();
  return OK;
}
