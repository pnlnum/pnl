
/*************************************************************************/
/* Written and (C) by David Pommier <pommier.david@gmail.com>            */  
/* 2008                                                                  */
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

/*************************************************************************/
/*These struct and functions are strongly inspired by work of            */
/* - JP Chancelier & al in NSP software                                  */
/*     http://cermics.enpc.fr/~jpc/nsp-tiddly/mine.html                  */ 
/*************************************************************************/

  

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "config.h"
#include "pnl_tridiag_matrix.h"
#include "pnl_matrix.h"
#include "pnl_mathtools.h"

static double __op_plus(double a, double b) { return a+b; }
static double __op_minus(double a, double b) { return a-b; }
static double __op_mult(double a, double b) { return a*b; }
static double __op_div(double a, double b) { return a/b; }
/*static double __op_inv(double a) { return 1.0/a; } */
/*static double __op_square(double a) { return a*a; } */

/***************************
 *** PnlTridiagmat functions ***
 ***************************/

/**
 * creates a PnlTriDiagMat
 * @param size number of rows
 * @return a PnlTriDiagMat pointer
 */
PnlTriDiagMat* pnl_tridiagmat_create(int size)
{
  PnlTriDiagMat *v;
  if((v=malloc(sizeof(PnlTriDiagMat)))==NULL)
    return NULL;
  v->size=size;
  v->owner = 1;
  if (v->size>0)
    {
      if((v->diag=malloc(v->size*sizeof(double)))==NULL)
        return NULL;
      if((v->diag_up=malloc(v->size*sizeof(double)))==NULL)
        return NULL;
      if((v->diag_down=malloc(v->size*sizeof(double)))==NULL)
        return NULL;
    }
  else
    {
      v->diag = (double*)NULL;
      v->diag_up = (double*)NULL;
      v->diag_down = (double*)NULL;
    }
  return v;  
}

/**
 * creates a PnlTriDiagMat
 * @param size number of rows
 * @param x used to fill the matrix with
 * @return a PnlTriDiagMat pointer
 */
PnlTriDiagMat* pnl_tridiagmat_create_from_double(int size, double x)
{
  PnlTriDiagMat *v;
  int i=0;
  double *ptr,*ptr_up,*ptr_down;
  
  if ((v=pnl_tridiagmat_create(size))==NULL)
    return NULL;

  ptr= pnl_tridiagmat_lget (v, 0, 0);
  ptr_up= pnl_tridiagmat_lget (v, 0, 1);
  ptr_down= pnl_tridiagmat_lget (v, 0, -1);
  *ptr=x;
  *ptr_up=x; 
  *ptr_down=0.0; 
  ptr++; ptr_up++;   
  ptr_down++; i++;
  while(i<v->size-1)
    {
      *ptr=x;
      *ptr_up=x; 
      *ptr_down=x; 
      ptr++; ptr_up++;   
      ptr_down++; i++;
    }
  *ptr=x;
  *ptr_up=0.0; 
  *ptr_down=x; 
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
  PnlTriDiagMat *v;
  int i=0;
  double *ptr,*ptr_up,*ptr_down;
  
  if ((v=pnl_tridiagmat_create(size))==NULL)
    return NULL;

  ptr= pnl_tridiagmat_lget (v, 0, 0);
  ptr_up= pnl_tridiagmat_lget (v, 0, 1);
  ptr_down= pnl_tridiagmat_lget (v, 0, -1);
  *ptr=x;
  *ptr_up=y; 
  *ptr_down=0.0; 
  ptr++; ptr_up++;   
  ptr_down++; i++;
  while(i<v->size-1)
    {
      *ptr=x;
      *ptr_up=y; 
      *ptr_down=y; 
      ptr++; ptr_up++;   
      ptr_down++; i++;
    }
  *ptr=x;
  *ptr_up=0.0; 
  *ptr_down=y; 
  return v;
}

/**
 * creates a PnlTriDiagMat
 * @param size number of rows
 * @param diag principal diagonal
 * @param diag_up diagonal v[i,i+1]
 * @param diag_down diagonal v[i,i-1]
 * @return a PnlTriDiagMat pointer
 */
PnlTriDiagMat* pnl_tridiagmat_create_from_ptr(int size, const double* diag, const double* diag_up, const double* diag_down)
{
  PnlTriDiagMat *v;
  
  if ((v=pnl_tridiagmat_create(size))==NULL)
    return NULL;
  
  memcpy(v->diag, diag, v->size*sizeof(double));
  memcpy(v->diag_up, diag_up, v->size*sizeof(double));
  memcpy(v->diag_down, diag_down, v->size*sizeof(double));
  return v;
}

/**
 * Reads a tri-diagonale matrix from a Matrix and creates the corresponding PnlTriDiagMat
 * all element extern to diagonal ,diag up and diag down is put to 0
 * @param mat a PnlMat
 * @return a PnlTriDiagMat
 */
PnlTriDiagMat* pnl_tridiagmat_create_from_matrix (const PnlMat * mat)
{
  PnlTriDiagMat *M;
  int i;
  CheckIsSquare (mat)
    M=pnl_tridiagmat_create(mat->m);
  for (i=0; i<mat->m-1; i++)
    {
      M->diag_up[i] = pnl_mat_get (mat, i, i+1);
      M->diag[i] = pnl_mat_get (mat, i, i);
      M->diag_down[i] = pnl_mat_get (mat, i+1, i);
    }
  M->diag[mat->m] = pnl_mat_get (mat, mat->m, mat->m);
  return M;
}


/**
 * frees a PnlTriDiagMat
 *
 * @param v adress of a PnlTriDiagMat*. v is set to NULL at exit.
 */
void pnl_tridiagmat_free(PnlTriDiagMat **v)
{
  if (*v != NULL)
    {
      if ((*v)->diag != NULL && (*v)->owner == 1) free((*v)->diag);
      if ((*v)->diag_up != NULL && (*v)->owner == 1) free((*v)->diag_up);
      if ((*v)->diag_down != NULL && (*v)->owner == 1) free((*v)->diag_down);
      free(*v);
      *v=NULL;
    }
}

/**
 * resizes a PnlTriDiagMat.
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
  if (v->owner == 0) return FAIL;
  if (size < 0) return FAIL;
  if (size==0) /* free array */
    {
      v->size = 0;
      free(v->diag); v->diag = NULL;
      free(v->diag_up); v->diag_up = NULL;
      free(v->diag_down); v->diag_down = NULL;
      return OK;
    }
  if (v->size >= size) /*nothing to do, just adjust m and n*/
    {
      v->size=size;
      return OK;
    }
  v->size=size;
  if (v->diag==NULL)
    {
      if (((v->diag = malloc(v->size*sizeof(double)))==NULL)
          || ((v->diag_up = malloc(v->size*sizeof(double)))==NULL) || ((v->diag_down = malloc(v->size*sizeof(double)))==NULL))
        return FAIL;
    }else
    {
      if (((v->diag = realloc(v->diag, v->size*sizeof(double)))==NULL)||((v->diag_up = realloc(v->diag_up, v->size*sizeof(double)))==NULL)
          ||((v->diag_down = realloc(v->diag_down, v->size*sizeof(double)))==NULL))
        
        return FAIL;
    }
  return OK;
}


/**
 * sets the value of if up== 1self[d,d+1]  =x, if up== -1 self[d-1,d]  = x if up== 0 self[d,d] 
 *
 * @param self : a PnlTriDiagMat
 * @param d : index of element
 * @param up : index of  diagonale
 * @param x : self[d,d+up]=x
 */
void pnl_tridiagmat_set (PnlTriDiagMat *self, int d, int up, double x)
{
  CheckIndexTriDiagMat(self,d, up);
  if(up==1)
    self->diag_up[d]= x;
  if(up==-1)
    self->diag_down[d]= x;
  if(up==0)
    self->diag[d]= x;
}

/**
 * gets the value of if up== 1self[d,d+1], if up== -1 self[d-1,d] if up== 0 self[d,d] 
 *
 * @param self : a PnlTriDiagMat
 * @param d : index of element
 * @param up : index of diagonale
 * @return  self[d,d+up]
 */
double pnl_tridiagmat_get(const PnlTriDiagMat *self, int d, int up)
{
  CheckIndexTriDiagMat(self,d, up);
  if(up==1)
    return self->diag_up[d];
  if(up==-1)
    return self->diag_down[d];
  if(up==0)
    return self->diag[d];
  return 0.0;
}

/**
 * returns the address of if up== 1self[d,d+1], if up== -1 self[d-1,d] 
 * if up== 0 self[d,d] for use as a lvalue 
 *
 * @param self : a PnlTriDiagMat
 * @param d : index of element 
 * @param up : index of diagonale
 * @return  &(self[d,d+up])
 */
double* pnl_tridiagmat_lget (PnlTriDiagMat *self, int d, int up)
{
  CheckIndexTriDiagMat(self,d,up);
  if(up==1)
    return &(self->diag_up[d]);
  if(up==-1)
    return &(self->diag_down[d]);
  if(up==0)
    return &(self->diag[d]);
  return NULL;
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
  double *lptr = pnl_tridiagmat_lget (lhs, 0, 0);
  double *lptr_up = pnl_tridiagmat_lget (lhs, 0, 1);
  double *lptr_down = pnl_tridiagmat_lget (lhs, 0, -1);
  int i=0;
  while (i<lhs->size)
    {
      (*lptr)= f(*lptr);  
      (*lptr_up)= f(*lptr_up);  
      (*lptr_down)= f(*lptr_down); 
      lptr++; lptr_up++;lptr_down++;i++;
    }
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
  double *lptr = pnl_tridiagmat_lget(lhs, 0, 0);
  double *lptr_up = pnl_tridiagmat_lget (lhs, 0, 1);
  double *lptr_down = pnl_tridiagmat_lget (lhs, 0, -1);
  int i=0;
  CheckTriDiagMatMatch(lhs, rhs);
  while (i<lhs->size)
    {
      (*lptr)=(*f)(*lptr,pnl_tridiagmat_get(rhs,i,0));
      (*lptr_up)=(*f)(*lptr_up,pnl_tridiagmat_get(rhs,i,1));
      (*lptr_down)=(*f)(*lptr_down,pnl_tridiagmat_get(rhs,i,-1));
      lptr++; lptr_up++; lptr_down++; i++;
    }
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
  double *lptr = pnl_tridiagmat_lget(lhs, 0, 0);
  double *lptr_up = pnl_tridiagmat_lget(lhs, 0, 1);
  double *lptr_down = pnl_tridiagmat_lget(lhs, 0, -1);
  int i=0;
  while (i<lhs->size)
    {
      (*lptr)= op(*lptr,x);  
      (*lptr_up)= op(*lptr_up,x);  
      (*lptr_down)= op(*lptr_down,x); 
      lptr++; lptr_up++;lptr_down++;i++;
    }
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
  /* impose Tri mat structure  */
  pnl_tridiagmat_set(lhs, 0, -1,0.0);
  pnl_tridiagmat_set(lhs, lhs->size-1, 1,0.0);
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
  PnlVect *lhs;
  double *rptr, *ptr, *ptr_up,*ptr_down,*lptr;
  int i;
  double anc_m1;
  double anc_0;

  CheckTriDiagMatVectIsCompatible(mat, vec);
  if ((lhs=pnl_vect_create(mat->size))==NULL)
    return NULL;
  rptr = vec->array;
  lptr=lhs->array;
  ptr = mat->diag;
  ptr_up = mat->diag_up;
  ptr_down = mat->diag_down;
  anc_m1=0.0;
  anc_0=(*rptr);
  rptr++;
  for (i=0; i<lhs->size-1; i++)
    {
      *lptr = (*ptr_down)*anc_m1+(*ptr) *anc_0 + (*ptr_up) *(*rptr);
      anc_m1=anc_0;
      anc_0=(*rptr);
      lptr++;rptr++;ptr++;ptr_up++;ptr_down++;
    
    }
  return lhs;
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
  double *rptr, *ptr, *ptr_up,*ptr_down,*lptr;
  int i;
  double anc_m1;
  double anc_0;

  CheckTriDiagMatVectIsCompatible(mat, lhs);
  CheckTriDiagMatVectIsCompatible(mat, rhs);
  rptr = rhs->array;
  lptr = lhs->array;
  ptr  = mat->diag;
  ptr_up = mat->diag_up;
  ptr_down = mat->diag_down;
  anc_m1=0.0;
  anc_0=(*rptr);
  rptr++;
  for (i=0; i<rhs->size-1; i++)
    {
      *lptr = (*ptr_down)*anc_m1+(*ptr) *anc_0 + (*ptr_up) *(*rptr);
      anc_m1=anc_0;
      anc_0=(*rptr);
      lptr++;rptr++;ptr++;ptr_up++;ptr_down++;
    }
  *lptr = (*ptr_down)*anc_m1+(*ptr) *anc_0;
}

/**
 * Prints a tri-diagonale matrix to a file
 *
 * @param fic a file descriptor.
 * @param M a PnlTriDiagMat pointer.
 */
void pnl_tridiagmat_fprint (FILE *fic, const PnlTriDiagMat *M)
{
  int i, up;
  for (i=0;i<M->size;i++)
    {
      for (up=-1;up<=1;up++)
        {
          fprintf (fic, "%7.4f " , pnl_tridiagmat_get (M, i, up));
        }
      fprintf (fic, "\n");
    }
}

/**
 * Prints a tri-diagonale matrix to the standard output
 *
 * @param M a PnlTriDiagMat pointer.
 */
void pnl_tridiagmat_print (const PnlTriDiagMat *M) { pnl_tridiagmat_fprint(stdout, M);}


/**
 * solves the linear system M x = b with M = LU.
 * For a symmetric definite
 * positive system, prefer pnl_tridiagmat_chol_syslin
 *
 * @param lhs existing vector that contains the solution on exit
 * @param M a PnlTriDiagMat 
 * @param rhs right hand side member
 */
void pnl_tridiagmat_lu_syslin (PnlVect *lhs, const PnlTriDiagMat *M,const PnlVect *rhs)
{
  int i;
  double *Diag,tmp;
  Diag = malloc((M->size)*sizeof(double));
  /* Do Down-Solve L y = rhs  */
  if((Diag[0]=pnl_tridiagmat_get(M,0,0))==0)
    {
      PNL_ERROR("division by zero","pnl_tridiagmat_lu_syslin");
    }
  pnl_vect_set(lhs,0,pnl_vect_get(rhs,0));
  for(i=1; i<M->size; i++)
    {
      tmp=pnl_tridiagmat_get(M,i,-1)/Diag[i-1];
      if((Diag[i]=pnl_tridiagmat_get(M,i,0)-tmp*pnl_tridiagmat_get(M,i-1,1))==0)
        {
           PNL_ERROR("division by zero","pnl_tridiagmat_lu_syslin");
        }
      pnl_vect_set(lhs,i,pnl_vect_get(rhs,i)-tmp*pnl_vect_get(lhs,i-1));
    }
  /* Do Up-Solve L y = rhs  */
  pnl_vect_set(lhs,M->size-1,pnl_vect_get(lhs,M->size-1)/Diag[M->size-1]);
  for(i=M->size-2;i>=0;i--)
    {
      pnl_vect_set(lhs,i,(pnl_vect_get(lhs,i)-pnl_tridiagmat_get(M,i,1)*
                          pnl_vect_get(lhs,i+1))/Diag[i]);
    }
  free(Diag);
}

