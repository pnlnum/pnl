/************************************************************************/
/* Copyright Jérôme Lelong <jerome.lelong@gmail.com>                    */
/* Copyright David Pommier <pommier.david@gmail.com>                    */
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



#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "pnl/pnl_config.h"
#include "pnl/pnl_tridiag_matrix.h"
#include "pnl/pnl_machine.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_mathtools.h"

static double __op_plus(double a, double b) { return a+b; }
static double __op_minus(double a, double b) { return a-b; }
static double __op_mult(double a, double b) { return a*b; }
static double __op_div(double a, double b) { return a/b; }

extern int C2F(dlagtm)(char *trans, int *n, int *nrhs, double *alpha,
                       double *dl, double *d, double *du, double *x,
                       int *ldx, double *beta, double *b, int *ldb);
extern int C2F(dgtsv)(int *n, int *nrhs, double *dl, double *d, double *du,
                       double *b, int *ldb, int *info);
extern int C2F(dgttrf)(int *n, double *dl, double *d, double *du, double
                       *du2, int *ipiv, int *info);
extern int C2F(dgttrs)(char *trans, int *n, int *nrhs, double *dl, double
                      *d, double *du, double *du2, int *ipiv, double *b,
                      int *ldb, int *info);

/**
 * Solve a tridiagonal system defined by matrix whose three diagonals are
 * defined by three values.
 * @param low the double value of M(i,i-1)
 * @param diag the double value of M(i,i)
 * @param up the double value of M(i,i+1)
 * @param rhs the right and side member of the system
 * @param lhs contains the solution on exit
 */
void pnl_progonka(const double low,
                  const double diag, const double up,
                  const PnlVect * rhs,
                  PnlVect * lhs)
{
  int i;
  double *D,tmp;
  D = malloc((rhs->size)*sizeof(double));
  /* Do Down-Solve L y = rhs  */
  if((D[0]=diag)==0)
    {
      PNL_ERROR("division by zero","pnl_progonka");
    }
  pnl_vect_set(lhs,0,pnl_vect_get(rhs,0));
  for(i=1; i<rhs->size; i++)
    {
      tmp=low/D[i-1];
      if((D[i]=diag-tmp*up)==0)
        {
           PNL_ERROR("division by zero","pnl_progonka");
        }
      pnl_vect_set(lhs,i,pnl_vect_get(rhs,i)-tmp*pnl_vect_get(lhs,i-1));
    }
  /* Do Up-Solve L y = rhs  */
  pnl_vect_set(lhs,rhs->size-1,pnl_vect_get(lhs,rhs->size-1)/D[rhs->size-1]);
  for(i=rhs->size-2;i>=0;i--)
    {
      pnl_vect_set(lhs,i,(pnl_vect_get(lhs,i)-up*pnl_vect_get(lhs,i+1))/D[i]);
    }
  free(D);

}

/*******************************
 *** PnlTridiagmat functions ***
 *******************************/

static char pnl_tridiag_mat_object_label[] = "PnlTridiagMatObject";
static char pnl_tridiag_mat_label[] = "PnlTridiagMat";

PnlTridiagMatObject* pnl_tridiag_mat_object_new ()
{
  PnlTridiagMatObject *o;
  if ( (o = malloc (sizeof(PnlTridiagMatObject))) == NULL) return NULL;
  o->size = 0;
  o->D = o->DU = o->DL = NULL;
  o->object.type = PNL_TYPE_TRIDIAG_MATRIX;
  o->object.parent_type = PNL_TYPE_TRIDIAG_MATRIX;
  o->object.label = pnl_tridiag_mat_object_label;
  o->object.destroy = (DestroyFunc *) pnl_tridiag_mat_object_free;
  return o;
}

/**
 * Free a PnlTridiagMatObject
 *
 * @param v a PnlTridiagMatObject**.
 */
void pnl_tridiag_mat_object_free(PnlTridiagMatObject **v)
{
  if (*v != NULL)
    {
      if ((*v)->D != NULL) free((*v)->D);
      if ((*v)->DU != NULL) free((*v)->DU);
      if ((*v)->DL != NULL) free((*v)->DL);
      free(*v);
      *v=NULL;
    }
}


/**
 * Create an empty PnlTridiagMat
 */
PnlTridiagMat* pnl_tridiag_mat_new()
{
  PnlTridiagMat *o;
  if ( (o = (PnlTridiagMat *) pnl_tridiag_mat_object_new ()) == NULL) return NULL;
  o->object.type = PNL_TYPE_TRIDIAG_MATRIX_DOUBLE;
  o->object.label = pnl_tridiag_mat_label;
  o->object.nref = 0;
  o->object.clone = (CloneFunc *) pnl_tridiag_mat_clone;
  o->object.copy = (CopyFunc *) pnl_tridiag_mat_copy;
  o->object.constructor = (NewFunc *) pnl_tridiag_mat_new;
  return o;
}

/**
 * Create a PnlTridiagMat
 * @param size number of rows
 * @return a PnlTridiagMat pointer
 */
PnlTridiagMat* pnl_tridiag_mat_create(int size)
{
  PnlTridiagMat *o;
  if ( (o = pnl_tridiag_mat_new ()) == NULL) return NULL;
  if (size>0)
    {
      o->size=size;
      if((o->D=malloc(size*sizeof(double)))==NULL) return NULL;
      if((o->DU=malloc((size-1)*sizeof(double)))==NULL) return NULL;
      if((o->DL=malloc((size-1)*sizeof(double)))==NULL) return NULL;
    }
  return o;
}

/**
 * Create a PnlTridiagMat
 * @param size number of rows and colusizes
 * @param x used to fill the diagonal matrix with
 * @param y used to fill the extra diagonal matrix with
 * @return a PnlTridiagMat pointer
 */
PnlTridiagMat* pnl_tridiag_mat_create_from_two_scalar(int size, double x, double y)
{
  PnlTridiagMat *mat;
  int i;

  if ((mat=pnl_tridiag_mat_create(size))==NULL) return NULL;

  for ( i=0 ; i<size-1 ; i++ )
    {
      mat->D[i] = x;
      mat->DL[i] = y;
      mat->DU[i] = y;
    }
  mat->D[size-1] = x;
  return mat;
}

/**
 * Create a PnlTridiagMat
 * @param size number of rows
 * @param x used to fill the three diagonals
 * @return a PnlTridiagMat pointer
 */
PnlTridiagMat* pnl_tridiag_mat_create_from_scalar(int size, double x)
{
  return pnl_tridiag_mat_create_from_two_scalar (size, x, x);
}

/**
 * Create a PnlTridiagMat
 * @param size number of rows
 * @param D principal diagonal
 * @param DU diagonal v[i,i+1]
 * @param DL diagonal v[i,i-1]
 * @return a PnlTridiagMat pointer
 */
PnlTridiagMat* pnl_tridiag_mat_create_from_ptr(int size, const double* DL, const double* D, const double* DU)
{
  PnlTridiagMat *v;

  if ((v=pnl_tridiag_mat_create(size))==NULL) return NULL;

  memcpy(v->D, D, v->size*sizeof(double));
  memcpy(v->DU, DU, (v->size-1)*sizeof(double));
  memcpy(v->DL, DL, (v->size-1)*sizeof(double));
  return v;
}

/**
 * Clone a PnlTridiagMat
 * @param T a PnlTridiagMat
 * @param clone a PnlTridiagMat, on exit contains a copy of T
 */
void pnl_tridiag_mat_clone(PnlTridiagMat *clone, const PnlTridiagMat *T)
{
  pnl_tridiag_mat_resize (clone, T->size);

  memcpy(clone->D, T->D, T->size*sizeof(double));
  memcpy(clone->DU, T->DU, (T->size-1)*sizeof(double));
  memcpy(clone->DL, T->DL, (T->size-1)*sizeof(double));
}

/**
 * Copie a PnlTridiagMat
 * @param T a PnlTridiagMat
 */
PnlTridiagMat* pnl_tridiag_mat_copy(const PnlTridiagMat *T)
{
  PnlTridiagMat *clone;
  clone = pnl_tridiag_mat_create (0);
  pnl_tridiag_mat_clone (clone, T);
  return clone;
}


/**
 * Create a tridiagonal matrix from a standard matrix by only taking into
 * account the three main diagonals.
 *
 * @param mat a PnlMat
 * @return a PnlTridiagMat
 */
PnlTridiagMat* pnl_tridiag_mat_create_from_mat (const PnlMat * mat)
{
  PnlTridiagMat *M;
  int i;
  CheckIsSquare (mat);

  M=pnl_tridiag_mat_create(mat->m);
  for (i=0; i<mat->m-1; i++)
    {
      M->DU[i] = pnl_mat_get (mat, i, i+1);
      M->D[i] = pnl_mat_get (mat, i, i);
      M->DL[i] = pnl_mat_get (mat, i+1, i);
    }
  M->D[mat->m-1] = pnl_mat_get (mat, mat->m-1, mat->m-1);
  return M;
}

/**
 * Create a standard matrix from a tridiagonal matrix.
 * All elements but those of the diagonal, upper diagonal and lower diagonal are set to 0
 * @param T a PnlTridiagMat
 * @return a PnlMat
 */
PnlMat* pnl_tridiag_mat_to_mat (const PnlTridiagMat * T)
{
  PnlMat *M;
  int i;
  M=pnl_mat_create_from_scalar(T->size, T->size,0.);
  for (i=0; i<T->size-1; i++)
    {
      PNL_MLET(M, i, i+1) = T->DU[i];
      PNL_MLET(M, i, i) = T->D[i];
      PNL_MLET(M, i+1, i) = T->DL[i];
    }
  PNL_MLET (M, T->size-1, T->size-1) = T->D[T->size-1];
  return M;
}


/**
 * Free a PnlTridiagMat
 *
 * @param v a PnlTridiagMat**.
 */
void pnl_tridiag_mat_free(PnlTridiagMat **v)
{
  PnlTridiagMatObject *o;
  o = PNL_TRIDIAGMAT_OBJECT (*v);
  pnl_tridiag_mat_object_free (&o);
}

/**
 * Resize a PnlTridiagMat.
 *
 * If the new size is smaller than the current one, no
 * memory is free. If the new size is larger than the
 * current one, more space is allocated. Note that for the
 * sake of efficiency the old data are not copied.
 *
 * @param v : a pointer to an already existing PnlTridiagMat (size
 * must be initialised)
 * @param size : new nb rows
 * @return OK or FAIL. When returns OK, the matrix is changed.
 */
int pnl_tridiag_mat_resize(PnlTridiagMat *v, int size)
{
  if (size < 0) return FAIL;
  if (size==0) /* free array */
    {
      v->size = 0;
      if (v->D != NULL) { free(v->D); v->D = NULL; }
      if (v->DU != NULL) { free(v->DU); v->DU = NULL; }
      if (v->DL != NULL) { free(v->DL); v->DL = NULL; }
      return OK;
    }
  if (v->size >= size) /*nothing to do, just adjust m and n*/
    {
      v->size=size;
      return OK;
    }
  v->size=size;
  if (v->D==NULL)
    {
      if (((v->D = malloc(v->size*sizeof(double)))==NULL)
          || ((v->DU = malloc((v->size-1)*sizeof(double)))==NULL)
          || ((v->DL = malloc((v->size-1)*sizeof(double)))==NULL))
        return FAIL;
    }
  else
    {
      free (v->D);
      free (v->DL);
      free (v->DU);
      if (((v->D = malloc(v->size*sizeof(double)))==NULL)
          ||((v->DU = malloc((v->size-1)*sizeof(double)))==NULL)
          ||((v->DL = malloc((v->size-1)*sizeof(double)))==NULL))
        return FAIL;
    }
  return OK;
}


/**
 * Set the value of self[d,d+up]=x
 *
 * @param self : a PnlTridiagMat
 * @param d : index of element
 * @param up : index of  diagonale
 * @param x : self[d,d+up]=x
 */
void pnl_tridiag_mat_set (PnlTridiagMat *self, int d, int up, double x)
{
  CheckIndexTridiagMat(self,d, up);
  if (up==1) self->DU[d]= x;
  else if (up==-1) self->DL[d-1]= x;
  else if (up==0) self->D[d]= x;
}

/**
 * Get the value of if self[d,d+up]
 *
 * @param self : a PnlTridiagMat
 * @param d : index of element
 * @param up : index of diagonale
 * @return  self[d,d+up]
 */
double pnl_tridiag_mat_get(const PnlTridiagMat *self, int d, int up)
{
  CheckIndexTridiagMat(self,d, up);
  if (up==1) return self->DU[d];
  else if (up==-1) return self->DL[d-1];
  else if (up==0) return self->D[d];
  else return 0.0;
}

/**
 * Return the address of self[d,d+up] to be used as a lvalue
 *
 * @param self : a PnlTridiagMat
 * @param d : index of element
 * @param up : index of diagonale
 * @return  &(self[d,d+up])
 */
double* pnl_tridiag_mat_lget (PnlTridiagMat *self, int d, int up)
{
  CheckIndexTridiagMat(self,d,up);
  if (up==1) return &(self->DU[d]);
  else if (up==-1) return &(self->DL[d-1]);
  else if (up==0) return &(self->D[d]);
  else return NULL;
}

/**
 * in-place map function
 *
 * @param lhs : left hand side vector
 * @param f : the function to be applied term by term
 * @return  lhs = f(lhs)
 */
void pnl_tridiag_mat_map_inplace(PnlTridiagMat *lhs, double(*f)(double))
{
  int i, size;
  size = lhs->size;
  for ( i=0 ; i<size-1 ; i++)
    {
      lhs->D[i] = f(lhs->D[i]);
      lhs->DL[i] = f(lhs->DL[i]);
      lhs->DU[i] = f(lhs->DU[i]);
    }
  lhs->D[size-1] = f(lhs->D[size-1]);
}


/**
 * map matrix componentwise
 *
 * @param lhs : each component lhs(i) contains f(lhs(i),rhs(i))
 * @param rhs : right hand side vector
 * @param f   : real function
 */
void pnl_tridiag_mat_map_tridiag_mat_inplace(PnlTridiagMat *lhs, const PnlTridiagMat *rhs,
                                   double(*f)(double,double))
{
  int i, size;
  CheckTridiagMatMatch(lhs, rhs);
  size = lhs->size;
  for ( i=0 ; i<size-1 ; i++)
    {
      lhs->D[i] = f(lhs->D[i], rhs->D[i]);
      lhs->DL[i] = f(lhs->DL[i], rhs->DL[i]);
      lhs->DU[i] = f(lhs->DU[i], rhs->DU[i]);
    }
  lhs->D[size-1] = f(lhs->D[size-1], rhs->D[size-1]);
}


/**
 * in-place matrix operator application
 *
 * @param lhs : left hand side vector
 * @param x : double arg
 * @param op : a binary operator, given as a function
 * @return  lhs = lhs op x
 */
static void __pnl_tridiag_mat_apply_op(PnlTridiagMat *lhs, double x, double (*op)(double, double))
{
  int i, size;
  size = lhs->size;
  for ( i=0 ; i<size-1 ; i++)
    {
      lhs->D[i] = op(lhs->D[i], x);
      lhs->DL[i] = op(lhs->DL[i], x);
      lhs->DU[i] = op(lhs->DU[i], x);
    }
  lhs->D[size-1] = op(lhs->D[size-1], x);
}

/**
 * in-place matrix matrix addition
 *
 * @param lhs : left hand side matrix
 * @param rhs : rigth hand side matrix
 * @return  lhs = lhs+rhs
 */
void pnl_tridiag_mat_plus_tridiag_mat(PnlTridiagMat *lhs, const PnlTridiagMat *rhs)
{
  pnl_tridiag_mat_map_tridiag_mat_inplace(lhs, rhs, __op_plus);
}

/**
 * in-place matrix matrix substraction
 *
 * @param lhs : left hand side matrix
 * @param rhs : rigth hand side matrix
 * @return  lhs = lhs+rhs
 */
void pnl_tridiag_mat_minus_tridiag_mat(PnlTridiagMat *lhs, const PnlTridiagMat *rhs)
{
  pnl_tridiag_mat_map_tridiag_mat_inplace(lhs, rhs, __op_minus);
}

/**
 * in-place matrix scalar addition
 *
 * @param lhs : left hand side matrix
 * @param  x : scalar
 * @return  lhs = lhs+x
 */
void pnl_tridiag_mat_plus_scalar(PnlTridiagMat *lhs, double x)
{
  __pnl_tridiag_mat_apply_op(lhs, x, __op_plus);
}

/**
 * in-place matrix scalar substraction
 *
 * @param lhs : left hand side matrix
 * @param  x : scalar
 * @return  lhs = lhs-x
 */
void pnl_tridiag_mat_minus_scalar(PnlTridiagMat *lhs, double x)
{
  __pnl_tridiag_mat_apply_op(lhs, x, __op_minus);
}

/**
 * in-place matrix scalar multiplication
 *
 * @param lhs : left hand side matrix
 * @param x : scalar
 * @return  lhs = lhs*x
 */
void pnl_tridiag_mat_mult_scalar(PnlTridiagMat *lhs, double x)
{
  __pnl_tridiag_mat_apply_op(lhs, x, __op_mult);
}


/**
 * in-place matrix scalar division
 *
 * @param lhs : left hand side matrix
 * @param x : scalar
 * @return  lhs = lhs*x
 */
void pnl_tridiag_mat_div_scalar(PnlTridiagMat *lhs, double x)
{
  __pnl_tridiag_mat_apply_op(lhs, x, __op_div);
}


/**
 * in-place matrix matrix term by term product
 *
 * @param lhs : left hand side matrix
 * @param rhs : rigth hand side matrix
 * @return  lhs = lhs.*rhs
 */
void pnl_tridiag_mat_mult_tridiag_mat_term(PnlTridiagMat *lhs, const PnlTridiagMat *rhs)
{
  pnl_tridiag_mat_map_tridiag_mat_inplace(lhs, rhs, __op_mult);
}

/**
 * in-place matrix matrix term by term division
 *
 * @param lhs : left hand side matrix
 * @param rhs : rigth hand side matrix
 * @return  lhs = lhs.*rhs
 */
void pnl_tridiag_mat_div_tridiag_mat_term(PnlTridiagMat *lhs, const PnlTridiagMat *rhs)
{
  pnl_tridiag_mat_map_tridiag_mat_inplace(lhs, rhs, __op_div);
}


/**
 *  in place matrix multiplication
 *
 * @param mat : matrix
 * @param lhs : vector
 * @param rhs : vector
 * @return  lhs=mat*rhs
 */
void pnl_tridiag_mat_mult_vect_inplace(PnlVect *lhs, const PnlTridiagMat *mat, const PnlVect *rhs)
{
  pnl_tridiag_mat_lAxpby (1., mat, rhs, 0., lhs);
}

/**
 *  matrix multiplication
 *
 * @param mat : matrix
 * @param vec : vector
 * @return  mat*vec
 */
PnlVect* pnl_tridiag_mat_mult_vect(const PnlTridiagMat *mat, const PnlVect *vec)
{
  PnlVect *lhs=pnl_vect_create (vec->size);
  pnl_tridiag_mat_lAxpby (1., mat, vec, 0., lhs);
  return lhs;
}


/**
 * Print a tridiagonal matrix to a file
 *
 * @param fic a file descriptor.
 * @param M a PnlTridiagMat pointer.
 */
void pnl_tridiag_mat_fprint (FILE *fic, const PnlTridiagMat *M)
{
  int i;
  for (i=0;i<M->size-1;i++)
    {
      fprintf (fic, "%f" , M->DL[i]);
      if ( i < M->size-2) printf (", ");
    }
  fprintf (fic, "\n");
  for (i=0;i<M->size;i++)
    {
      fprintf (fic, "%f " , M->D[i]);
      if ( i < M->size-1) printf (", ");
    }
  fprintf (fic, "\n");
  for (i=0;i<M->size-1;i++)
    {
      fprintf (fic, "%f " , M->DU[i]);
      if ( i < M->size-2) printf (", ");
    }
  fprintf (fic, "\n");
}

/**
 * Print a tridiagonal matrix to the standard output
 *
 * @param M a PnlTridiagMat pointer.
 */
void pnl_tridiag_mat_print (const PnlTridiagMat *M)
{
  pnl_tridiag_mat_fprint(stdout, M);
}


/**
 * Solve the linear system A x = b
 *
 * @param A a PnlTridiagMat. On exit, A is modified and becomes unusable.
 * @param b right hand side member. On exit b contains the solution x
 * @return FAIL or OK
 */
int pnl_tridiag_mat_syslin_inplace (PnlTridiagMat *A, PnlVect *b)
{
  int n, nrhs, info, ldb;
  n = A->size;
  nrhs = 1;
  ldb = A->size;

  PNL_CHECK (A->size != b->size, "incompatible size", "tridiag_lu");

  C2F(dgtsv)(&n, &nrhs, A->DL, A->D, A->DU, b->array, &ldb, &info);
  if (info == 0) return OK;
  else return FAIL;
}

/**
 * Solve the linear system A x = b
 *
 * @param A a PnlTridiagMat. On exit, A is modified and becomes unusable.
 * @param x contains the solution on exit
 * @param b right hand side member.
 * @return FAIL or OK
 */
int pnl_tridiag_mat_syslin (PnlVect *x, PnlTridiagMat *A, const PnlVect *b)
{
  pnl_vect_clone (x, b);
  return pnl_tridiag_mat_syslin_inplace (A, x);
}


/**
 * Compute the scalar product <x,T * y>
 *
 * @param x a real vector
 * @param T a tridiagonal matrix
 * @param y a real vector
 * @return  x' * A * y
 */
double pnl_tridiag_mat_scalar_prod (const PnlTridiagMat *T, const PnlVect *x, const PnlVect *y)
{
  int n, k;
  double sum;

  PNL_CHECK (T->size != y->size || x->size != T->size, "size mismatched", "tridiag_scalar_prod");
  n = x->size;
  sum = 0.;

  for ( k=1 ; k<n-1 ; k++ )
    {
      sum += PNL_GET(x,k) * (T->DL[k-1] * PNL_GET(y,k-1) + T->D[k] * PNL_GET(y,k) +
                             T->DU[k] * PNL_GET(y,k+1));
    }
  sum += PNL_GET(x,0) * (T->D[0] * PNL_GET(y,0) + T->DU[0] * PNL_GET(y,1) );
  sum += PNL_GET(x, n-1) * (T->DL[n-2] * PNL_GET(y,n-2) + T->D[n-1] * PNL_GET(y, n-1));
  return sum;
}

/**
 * compute y=l * A * x + b * y
 *
 * @param l  : double.
 * @param A : PnlTridiagMat
 * @param x : vector
 * @param b : double
 * @param y : vector
 */
void pnl_tridiag_mat_lAxpby (double l, const PnlTridiagMat *A, const PnlVect *x,
                            double b, PnlVect * y)
{
  int n, nrhs;
  double alpha, beta;
  if ( b==0.0 )
    {
      beta = 0.;
      pnl_vect_resize (y, x->size);
    }
  else
    {
      if (l == 0. || l == 1. || l == 1.)
        {
          alpha = l;
        }
      else
        {
          b /= l;
        }
      if ( b == 1. || b == 1.)
        {
          beta = b;
        }
      else
        {
          beta = 1.;
          pnl_vect_mult_scalar (y, b);
        }
    }

  n = A->size;
  nrhs = 1;
  alpha = 1;
  C2F(dlagtm)("N", &n, &nrhs, &alpha, A->DL, A->D,
              A->DU, x->array, &n, &beta, y->array, &n);
  if (l != 0. && l != 1.&& l != 1.)
    {
      pnl_vect_mult_scalar (y, l);
    }
}

/*
 * PnlTridiagMatLU
 */
static char pnl_tridiag_mat_lu_object_label[] = "PnlTridiagMatLUObject";
static char pnl_tridiag_mat_lu_label[] = "PnlTridiagMatLU";

PnlTridiagMatLUObject* pnl_tridiag_mat_lu_object_new ()
{
  PnlTridiagMatLUObject *o;
  if ( (o = malloc (sizeof(PnlTridiagMatLUObject))) == NULL) return NULL;
  o->size = 0;
  o->D = o->DU = o->DL = o->DU2 = o->ipiv = NULL;
  o->object.type = PNL_TYPE_TRIDIAG_MATRIX_LU;
  o->object.parent_type = PNL_TYPE_TRIDIAG_MATRIX_LU;
  o->object.label = pnl_tridiag_mat_lu_object_label;
  o->object.destroy = (DestroyFunc *) pnl_tridiag_mat_lu_object_free;
  return o;
}

/**
 * Free a PnlTridiagMatLUObject
 *
 * @param v a PnlTridiagMatLUObject**.
 */
void pnl_tridiag_mat_lu_object_free(PnlTridiagMatLUObject **v)
{
  if (*v != NULL)
    {
      if ((*v)->D != NULL) free((*v)->D);
      if ((*v)->DU != NULL) free((*v)->DU);
      if ((*v)->DU2 != NULL) free((*v)->DU2);
      if ((*v)->DL != NULL) free((*v)->DL);
      if ((*v)->ipiv != NULL) free((*v)->ipiv);
      free(*v);
      *v=NULL;
    }
}


/**
 * Create an empty PnlTridiagMatLU
 */
PnlTridiagMatLU* pnl_tridiag_mat_lu_new()
{
  PnlTridiagMatLU *o;
  if ( (o = (PnlTridiagMatLU *) pnl_tridiag_mat_lu_object_new ()) == NULL) return NULL;
  o->object.type = PNL_TYPE_TRIDIAG_MATRIX_LU_DOUBLE;
  o->object.label = pnl_tridiag_mat_lu_label;
  o->object.clone = (CloneFunc *) pnl_tridiag_mat_lu_clone;
  o->object.copy = (CopyFunc *) pnl_tridiag_mat_lu_copy;
  o->object.constructor = (NewFunc *) pnl_tridiag_mat_lu_new;
  return o;
}

/**
 * Create a PnlTridiagMatLU
 * @param size number of rows
 * @return a PnlTridiagMatLU pointer
 */
PnlTridiagMatLU* pnl_tridiag_mat_lu_create(int size)
{
  PnlTridiagMatLU *o;
  if ( (o = pnl_tridiag_mat_lu_new ()) == NULL) return NULL;
  if (size>0)
    {
      o->size=size;
      if((o->D=malloc(size*sizeof(double)))==NULL) return NULL;
      if((o->DU=malloc((size-1)*sizeof(double)))==NULL) return NULL;
      if((o->DU2=malloc((size-2)*sizeof(double)))==NULL) return NULL;
      if((o->DL=malloc((size-1)*sizeof(double)))==NULL) return NULL;
      if((o->ipiv=malloc(size*sizeof(int)))==NULL) return NULL;
    }
  return o;
}

/**
 * Clone a PnlTridiagMatLU
 * @param T a PnlTridiagMatLU
 * @param clone a PnlTridiagMatLU, on exit contains a copy of T
 */
void pnl_tridiag_mat_lu_clone(PnlTridiagMatLU *clone, const PnlTridiagMatLU *T)
{
  pnl_tridiag_mat_lu_resize (clone, T->size);

  memcpy(clone->D, T->D, T->size*sizeof(double));
  memcpy(clone->DU, T->DU, (T->size-1)*sizeof(double));
  memcpy(clone->DU2, T->DU2, (T->size-2)*sizeof(double));
  memcpy(clone->DL, T->DL, (T->size-1)*sizeof(double));
  memcpy(clone->ipiv, T->ipiv, T->size*sizeof(int));
}

/**
 * Copie a PnlTridiagMatLU
 * @param T a PnlTridiagMatLU
 */
PnlTridiagMatLU* pnl_tridiag_mat_lu_copy(const PnlTridiagMatLU *T)
{
  PnlTridiagMatLU *clone;
  clone = pnl_tridiag_mat_lu_create (T->size);
  pnl_tridiag_mat_lu_clone (clone, T);
  return clone;
}

/**
 * Free a PnlTridiagMatLU
 *
 * @param v a PnlTridiagMatLU**.
 */
void pnl_tridiag_mat_lu_free(PnlTridiagMatLU **v)
{
  PnlTridiagMatLUObject *o;
  o = PNL_TRIDIAGMATLU_OBJECT (*v);
  pnl_tridiag_mat_lu_object_free (&o);
}

/**
 * Resize a PnlTridiagMatLU.
 *
 * If the new size is smaller than the current one, no
 * memory is free. If the new size is larger than the
 * current one, more space is allocated. Note that for the
 * sake of efficiency the old data are not copied.
 *
 * @param v : a pointer to an already existing PnlTridiagMatLU (size
 * must be initialised)
 * @param size : new nb rows
 * @return OK or FAIL. When returns OK, the matrix is changed.
 */
int pnl_tridiag_mat_lu_resize(PnlTridiagMatLU *v, int size)
{
  if (size < 0) return FAIL;
  if (size==0) /* free array */
    {
      v->size = 0;
      if (v->D != NULL) { free(v->D); v->D = NULL; }
      if (v->DU != NULL) { free(v->DU); v->DU = NULL; }
      if (v->DU2 != NULL) { free(v->DU2); v->DU2 = NULL; }
      if (v->DL != NULL) { free(v->DL); v->DL = NULL; }
      if (v->ipiv != NULL) { free(v->ipiv); v->ipiv = NULL; }
      return OK;
    }
  if (v->size >= size) /*nothing to do, just adjust m and n*/
    {
      v->size=size;
      return OK;
    }
  v->size=size;
  if (v->size > 0)
    {
      free (v->D);
      free (v->DL);
      free (v->DU);
      free (v->DU2);
      free (v->ipiv);
    }
  if (((v->D = malloc(v->size*sizeof(double)))==NULL)
      || ((v->DU = malloc((v->size-1)*sizeof(double)))==NULL)
      || ((v->DU2 = malloc((v->size-2)*sizeof(double)))==NULL)
      || ((v->DL = malloc((v->size-1)*sizeof(double)))==NULL)
      || ((v->ipiv = malloc((v->size)*sizeof(int)))==NULL))
    return FAIL;
  return OK;
}

/** 
 * Compute the LU factorisation of a tridiagonal matrix
 * 
 * @param LU the LU decomposition. LU must be an already existing
 * PnlTridiagMatLU as returned by pnl_tridiag_mat_lu_new
 * @param A the initial tridiagonal matrix
 * 
 * @return  OK or FAIL
 */
int pnl_tridiag_mat_lu_compute (PnlTridiagMatLU *LU, const PnlTridiagMat *A)
{
  int n, info;
  n = A->size;
  pnl_tridiag_mat_lu_resize (LU, A->size);
  memcpy(LU->D, A->D, A->size*sizeof(double));
  memcpy(LU->DU, A->DU, (A->size-1)*sizeof(double));
  memcpy(LU->DL, A->DL, (A->size-1)*sizeof(double));
  C2F(dgttrf)(&n, LU->DL, LU->D, LU->DU, LU->DU2, LU->ipiv, &info);
  if ( info != 0 )
    {
      PNL_MESSAGE_ERROR ("LU decomposition cannot be computed", "pnl_tridiag_mat_lu_compute");
      return FAIL;
    }
  return OK;
}

/**
 * Solve the linear system A x = b
 *
 * @param LU a PnlTridiagMatLU. A must have been computed by a previous call
 * to pnl_mat_lu_compute
 * @param b right hand side member. On exit b contains the solution x
 * @return FAIL or OK
 */
int pnl_tridiag_mat_lu_syslin_inplace (PnlTridiagMatLU *LU, PnlVect *b)
{
  int n, nrhs, info, ldb;
  n = LU->size;
  nrhs = 1;
  ldb = LU->size;

  PNL_CHECK (LU->size != b->size, "incompatible size", "tridiag_mat_lu_syslin");

  C2F(dgttrs)("N", &n, &nrhs, LU->DL, LU->D, LU->DU, LU->DU2, LU->ipiv, b->array, &ldb, &info);
  if (info == 0) return OK;
  else return FAIL;
}

/**
 * Solve the linear system A x = b
 *
 * @param LU a PnlTridiagMatLU. A must have been computed by a previous call
 * to pnl_mat_lu_compute
 * @param x contains the solution on exit
 * @param b right hand side member.
 * @return FAIL or OK
 */
int pnl_tridiag_mat_lu_syslin (PnlVect *x, PnlTridiagMatLU *LU, const PnlVect *b)
{
  pnl_vect_clone (x, b);
  return pnl_tridiag_mat_lu_syslin_inplace (LU, x);
}

