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

#include "config.h"
#include "pnl_tridiag_matrix.h"
#include "pnl_machine.h"
#include "pnl_matrix.h"
#include "pnl_mathtools.h"

static double __op_plus(double a, double b) { return a+b; }
static double __op_minus(double a, double b) { return a-b; }
static double __op_mult(double a, double b) { return a*b; }
static double __op_div(double a, double b) { return a/b; }

extern int C2F(dlagtm)(char *trans, int *n, int *nrhs, double *alpha, 
                       double *dl, double *d, double *du, double *x,
                       int *ldx, double *beta, double *b, int *ldb);
extern int C2F(dgtsv)(int *n, int *nrhs, double *dl, double *d, double *du, 
                       double *b, int *ldb, int *info);

/**
 * Solves a tridiagonal system defined by matrix whose three diagonals are
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
 *** PnlTriDiagmat functions ***
 *******************************/

/**
 * creates a PnlTriDiagMat
 * @param size number of rows
 * @return a PnlTriDiagMat pointer
 */
PnlTriDiagMat* pnl_tridiagmat_create(int size)
{
  PnlTriDiagMat *v;
  if((v=malloc(sizeof(PnlTriDiagMat)))==NULL) return NULL;
  v->size=size;
  if (size>0)
    {
      if((v->D=malloc(size*sizeof(double)))==NULL) return NULL;
      if((v->DU=malloc((size-1)*sizeof(double)))==NULL) return NULL;
      if((v->DL=malloc((size-1)*sizeof(double)))==NULL) return NULL;
    }
  else
    {
      v->D = (double*)NULL;
      v->DU = (double*)NULL;
      v->DL = (double*)NULL;
    }
  return v;  
}

/**
 * creates a PnlTriDiagMat
 * @param size number of rows and colusizes
 * @param x used to fill the diagonal matrix with
 * @param y used to fill the extra diagonal matrix with
 * @return a PnlTriDiagMat pointer
 */
PnlTriDiagMat* pnl_tridiagmat_create_from_two_double(int size, double x, double y)
{
  PnlTriDiagMat *mat;
  double *ptr,*ptr_up,*ptr_down;
  int i=0;
  
  if ((mat=pnl_tridiagmat_create(size))==NULL) return NULL;

  ptr = mat->D;
  ptr_up = mat->DU;
  ptr_down = mat->DL;
  while(i<mat->size-1)
    {
      *ptr=x;
      *ptr_up=y; 
      *ptr_down=y; 
      ptr++; ptr_up++;   
      ptr_down++; i++;
    }
  *ptr=x;
  return mat;
}

/**
 * creates a PnlTriDiagMat
 * @param size number of rows
 * @param x used to fill the three diagonals 
 * @return a PnlTriDiagMat pointer
 */
PnlTriDiagMat* pnl_tridiagmat_create_from_double(int size, double x)
{
  return pnl_tridiagmat_create_from_two_double (size, x, x);
}

/**
 * creates a PnlTriDiagMat
 * @param size number of rows
 * @param D principal diagonal
 * @param DU diagonal v[i,i+1]
 * @param DL diagonal v[i,i-1]
 * @return a PnlTriDiagMat pointer
 */
PnlTriDiagMat* pnl_tridiagmat_create_from_ptr(int size, const double* DL, const double* D, const double* DU)
{
  PnlTriDiagMat *v;
  
  if ((v=pnl_tridiagmat_create(size))==NULL) return NULL;
  
  memcpy(v->D, D, v->size*sizeof(double));
  memcpy(v->DU, DU, (v->size-1)*sizeof(double));
  memcpy(v->DL, DL, (v->size-1)*sizeof(double));
  return v;
}

/** 
 * Creates a tridiagonal matrix from a standard matrix by only taking into
 * account the three main diagonals.
 *
 * @param mat a PnlMat 
 * @return a PnlTriDiagMat
 */
PnlTriDiagMat* pnl_tridiagmat_create_from_matrix (const PnlMat * mat) 
{
  PnlTriDiagMat *M; 
  int i; 
  CheckIsSquare (mat); 

  M=pnl_tridiagmat_create(mat->m);
  for (i=0; i<mat->m-1; i++) 
    { 
      M->DU[i] = pnl_mat_get (mat, i, i+1); 
      M->D[i] = pnl_mat_get (mat, i, i);
      M->DL[i] = pnl_mat_get (mat, i+1, i); 
    }
  M->D[mat->m] = pnl_mat_get (mat, mat->m-1, mat->m-1); 
  return M; 
}

/**
 * Creates a standard matrix from a tridiagonal matrix. 
 * All elements buth those of the diagonal, upper diagonal and lower diagonal are set to 0
 * @param T a PnlTriDiagMat
 * @return a PnlMat
 */
PnlMat* pnl_tridiagmat_to_matrix (const PnlTriDiagMat * T)
{
  PnlMat *M;
  int i;
  M=pnl_mat_create(T->size, T->size);
  for (i=0; i<T->size-1; i++)
    {
      PNL_MLET(M, i, i+1) = T->DU[i];
      PNL_MLET(M, i, i) = T->D[i];
      PNL_MLET(M, i+1, i) = T->DL[i];
    }
  PNL_MLET (M, i, i) = T->D[T->size-1];
  return M;
}


/**
 * Frees a PnlTriDiagMat
 *
 * @param v a PnlTriDiagMat**. 
 */
void pnl_tridiagmat_free(PnlTriDiagMat **v)
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
 * Resizes a PnlTriDiagMat.
 *
 * If the new size is smaller than the current one, no
 * memory is free. If the new size is larger than the
 * current one, more space is allocated. Note that for the
 * sake of efficiency the old data are not copied.
 *
 * @param v : a pointer to an already existing PnlTriDiagMat (size
 * must be initialised)
 * @param size : new nb rows
 * @return OK or FAIL. When returns OK, the matrix is changed. 
 */
int pnl_tridiagmat_resize(PnlTriDiagMat *v, int size)
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
 * Sets the value of self[d,d+up]=x
 *
 * @param self : a PnlTriDiagMat
 * @param d : index of element
 * @param up : index of  diagonale
 * @param x : self[d,d+up]=x
 */
void pnl_tridiagmat_set (PnlTriDiagMat *self, int d, int up, double x)
{
  CheckIndexTriDiagMat(self,d, up);
  if (up==1) self->DU[d]= x;
  else if (up==-1) self->DL[d-1]= x;
  else if (up==0) self->D[d]= x;
}

/**
 * gets the value of if self[d,d+up]
 *
 * @param self : a PnlTriDiagMat
 * @param d : index of element
 * @param up : index of diagonale
 * @return  self[d,d+up]
 */
double pnl_tridiagmat_get(const PnlTriDiagMat *self, int d, int up)
{
  CheckIndexTriDiagMat(self,d, up);
  if (up==1) return self->DU[d];
  else if (up==-1) return self->DL[d-1];
  else if (up==0) return self->D[d];
  else return 0.0;
}

/**
 * Returns the address of self[d,d+up] to be used as a lvalue
 *
 * @param self : a PnlTriDiagMat
 * @param d : index of element 
 * @param up : index of diagonale
 * @return  &(self[d,d+up])
 */
double* pnl_tridiagmat_lget (PnlTriDiagMat *self, int d, int up)
{
  CheckIndexTriDiagMat(self,d,up);
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
void pnl_tridiagmat_map_inplace(PnlTriDiagMat *lhs, double(*f)(double))
{
  double *lptr = lhs->D;
  double *lptr_up = lhs->DU;
  double *lptr_down = lhs->DL;
  int i=0;
  while (i<lhs->size-1)
    {
      (*lptr) = f(*lptr);  
      (*lptr_up) = f(*lptr_up);  
      (*lptr_down) = f(*lptr_down); 
      lptr++; lptr_up++; lptr_down++; i++;
    }
      (*lptr) = f(*lptr);  
}


/**
 * map matrix componentwise
 *
 * @param lhs : each component lhs(i) contains f(lhs(i),rhs(i))
 * @param rhs : right hand side vector
 * @param f   : real function 
 */
void pnl_tridiagmat_map_tridiagmat(PnlTriDiagMat *lhs, const PnlTriDiagMat *rhs, double(*f)(double,double))
{
  double *lptr = lhs->D;
  double *lptr_up = lhs->DU;
  double *lptr_down = lhs->DL;
  double *rptr = rhs->D;
  double *rptr_up = rhs->DU;
  double *rptr_down = rhs->DL;
  int i=0;
  CheckTriDiagMatMatch(lhs, rhs);
  while (i<lhs->size-1)
    {
      (*lptr)=(*f)(*lptr,*rptr);
      (*lptr_up)=(*f)(*lptr_up,*rptr_up);
      (*lptr_down)=(*f)(*lptr_down,*rptr_down);
      rptr++; rptr_up++; rptr_down++;
      lptr++; lptr_up++; lptr_down++; i++;
    }
  (*lptr)=(*f)(*lptr,*rptr);
}


/**
 * in-place matrix operator application
 *
 * @param lhs : left hand side vector
 * @param x : double arg
 * @param op : a binary operator, given as a function
 * @return  lhs = lhs op x
 */
static void __pnl_tridiagmat_apply_op(PnlTriDiagMat *lhs, double x, double (*op)(double, double))
{
  double *lptr = lhs->D;
  double *lptr_up = lhs->DU;
  double *lptr_down = lhs->DL;
  int i=0;
  while (i<lhs->size-1)
    {
      (*lptr)= op(*lptr,x);  
      (*lptr_up)= op(*lptr_up,x);  
      (*lptr_down)= op(*lptr_down,x); 
      lptr++; lptr_up++;lptr_down++;i++;
    }
  (*lptr)= op(*lptr,x);  
}

/**
 * in-place matrix matrix addition
 *
 * @param lhs : left hand side matrix
 * @param rhs : rigth hand side matrix
 * @return  lhs = lhs+rhs
 */
void pnl_tridiagmat_plus_tridiagmat(PnlTriDiagMat *lhs, const PnlTriDiagMat *rhs)
{
  pnl_tridiagmat_map_tridiagmat(lhs, rhs, __op_plus);
}

/**
 * in-place matrix matrix substraction
 *
 * @param lhs : left hand side matrix
 * @param rhs : rigth hand side matrix
 * @return  lhs = lhs+rhs
 */
void pnl_tridiagmat_minus_tridiagmat(PnlTriDiagMat *lhs, const PnlTriDiagMat *rhs)
{
  pnl_tridiagmat_map_tridiagmat(lhs, rhs, __op_minus);
}
  
/**
 * in-place matrix scalar addition
 *
 * @param lhs : left hand side matrix
 * @param  x : scalar
 * @return  lhs = lhs+x
 */
void pnl_tridiagmat_plus_double(PnlTriDiagMat *lhs, double x)
{
  __pnl_tridiagmat_apply_op(lhs, x, __op_plus);
}

/**
 * in-place matrix scalar substraction
 *
 * @param lhs : left hand side matrix
 * @param  x : scalar
 * @return  lhs = lhs-x
 */
void pnl_tridiagmat_minus_double(PnlTriDiagMat *lhs, double x)
{
  __pnl_tridiagmat_apply_op(lhs, x, __op_minus);
}

/**
 * in-place matrix scalar multiplication
 *
 * @param lhs : left hand side matrix
 * @param x : scalar
 * @return  lhs = lhs*x
 */
void pnl_tridiagmat_mult_double(PnlTriDiagMat *lhs, double x)
{
  __pnl_tridiagmat_apply_op(lhs, x, __op_mult);
}  


/**
 * in-place matrix scalar division
 *
 * @param lhs : left hand side matrix
 * @param x : scalar
 * @return  lhs = lhs*x
 */
void pnl_tridiagmat_div_double(PnlTriDiagMat *lhs, double x)
{
  __pnl_tridiagmat_apply_op(lhs, x, __op_div);
}


/**
 * in-place matrix matrix term by term product
 *
 * @param lhs : left hand side matrix
 * @param rhs : rigth hand side matrix
 * @return  lhs = lhs.*rhs
 */
void pnl_tridiagmat_mult_tridiagmat_term(PnlTriDiagMat *lhs, const PnlTriDiagMat *rhs)
{
  pnl_tridiagmat_map_tridiagmat(lhs, rhs, __op_mult);
}

/**
 * in-place matrix matrix term by term division
 *
 * @param lhs : left hand side matrix
 * @param rhs : rigth hand side matrix
 * @return  lhs = lhs.*rhs
 */
void pnl_tridiagmat_div_tridiagmat_term(PnlTriDiagMat *lhs, const PnlTriDiagMat *rhs)
{
  pnl_tridiagmat_map_tridiagmat(lhs, rhs, __op_div);
}


/**
 *  in place matrix multiplication
 *
 * @param mat : matrix
 * @param lhs : vector
 * @param rhs : vector
 * @return  lhs=mat*rhs
 */
void pnl_tridiagmat_mult_vect_inplace(PnlVect *lhs, const PnlTriDiagMat *mat, const PnlVect *rhs)
{
  pnl_tridiagmat_lAxpby (1., mat, rhs, 0., lhs);
}

/**
 *  matrix multiplication
 *
 * @param mat : matrix
 * @param vec : vector
 * @return  mat*vec
 */
PnlVect* pnl_tridiagmat_mult_vect(const PnlTriDiagMat *mat, const PnlVect *vec)
{
  PnlVect *lhs=pnl_vect_create (vec->size);
  pnl_tridiagmat_lAxpby (1., mat, vec, 0., lhs);
  return lhs;
}


/**
 * Prints a tridiagonal matrix to a file
 *
 * @param fic a file descriptor.
 * @param M a PnlTriDiagMat pointer.
 */
void pnl_tridiagmat_fprint (FILE *fic, const PnlTriDiagMat *M)
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
 * Prints a tridiagonal matrix to the standard output
 *
 * @param M a PnlTriDiagMat pointer.
 */
void pnl_tridiagmat_print (const PnlTriDiagMat *M) 
{
  pnl_tridiagmat_fprint(stdout, M);
}


/**
 * solves the linear system A x = b
 *
 * @param A a PnlTriDiagMat 
 * @param b right hand side member. On exit b contains the solution x
 * @return FAIL or OK
 */
int pnl_tridiagmat_lu_syslin_inplace (const PnlTriDiagMat *A, PnlVect *b)
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
 * solves the linear system A x = b
 *
 * @param A a PnlTriDiagMat 
 * @param x contains the solution on exit
 * @param b right hand side member. 
 * @return FAIL or OK
 */
int pnl_tridiagmat_lu_syslin (PnlVect *x, const PnlTriDiagMat *A, const PnlVect *b)
{
  pnl_vect_clone (x, b);
  return pnl_tridiagmat_lu_syslin_inplace (A, x);
}


/**
 * compute scalar product <lhs,A rhs >
 *
 * @param lhs : vector
 * @param mat : matrix
 * @param rhs : vector
 * @return  =lhs'*mat*rhs
 */
double pnl_tridiagmat_scalar_prod (const PnlVect *lhs, const PnlTriDiagMat *mat,const PnlVect *rhs)
{
  double sum=0.;
  double term;
  double *ptr, *ptr_up,*ptr_down,*lptr,*rptr;
  int i;
  double anc_m1;
  double anc_0;

  lptr = lhs->array;
  rptr = rhs->array;
  ptr  = mat->D;
  ptr_up = mat->DU;
  ptr_down = mat->DL;
  anc_m1=0.0;
  anc_0=(*rptr);
  for (i=0; i<lhs->size-1; i++,ptr++,ptr_up++,ptr_down++)
    {
      term= (*ptr_down)*anc_m1+(*ptr) *anc_0 + (*ptr_up) *(*(rptr+1));
      sum+=(*lptr)*term;
      anc_m1=anc_0;
      lptr++;
      rptr++;
      anc_0=(*rptr);
    }
  term= (*ptr_down)*anc_m1+(*ptr) *anc_0;
  sum+=(*lptr)*term;
  return sum;
}

/**
 * compute y=l * A * x + b * y 
 *
 * @param l  : double.
 * @param A : PnlTriDiagMat
 * @param x : vector
 * @param b : double
 * @param y : vector
 */
void pnl_tridiagmat_lAxpby (double l, const PnlTriDiagMat *A, const PnlVect *x,
                            double b, PnlVect * y)
{
  int n, nrhs, ldb;
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
          pnl_vect_mult_double (y, b);
        }
    }

  n = A->size;
  ldb = A->size;
  nrhs = 1;
  alpha = 1;
  C2F(dlagtm)("N", &n, &nrhs, &alpha, A->DL, A->D, 
              A->DU, x->array, &n, &beta, y->array, &n);
  if (l != 0. && l != 1.&& l != 1.)
    {
      pnl_vect_mult_double (y, l);
    }
}
