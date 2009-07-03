#ifndef SPARSE_GRID_CONSTRUCTOR_H
#define SPARSE_GRID_CONSTRUCTOR_H
/* Use of C++ code to construct Sparse Grid,  */
/* needed stl container like map ... */
#include "gridsparse_constructor.h"

void complexity_sparse(int d);
int Size_GridSparse(int dim,int lev);
int Size_GridSparse_With_Bnd(int dim,int lev);
void create_grid_sparse_cpp(int dim,
			    int lev,
			    GridSparse  * G);
#endif
