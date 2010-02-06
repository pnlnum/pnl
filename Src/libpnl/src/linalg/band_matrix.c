/************************************************************************/
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

/*************************************************************************/
/*These struct and functions are strongly inspired by work of            */
/* - F.Hecht & al in RMN class distributed on Freefem Project            */
/*     http://www.freefem.org/                                           */ 
/*************************************************************************/

  

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
  
#include "config.h"
#include "pnl_band_matrix.h"
#include "pnl_matrix.h"
#include "pnl_mathtools.h"
#include "pnl_array.h"


/** 
 * pnl_band_matrix_create 
 * create a void sparse matrix of size n with "band" extra diagonal terms 
 *
 * @param n int size of diagonal square matrix.
 * @param band number of extra-diagonal in matrix.
 * @return adress of PnlBandMatrix.
 */
PnlBandMatrix* pnl_band_matrix_create(const int  n,int band)
{
  int k,i;
  PnlBandMatrix * M = malloc(sizeof(PnlBandMatrix));
  M->n=n;
  M->m=n;/*Square Matrix */
  M->owner=1;
  M->D=malloc(sizeof(double)*n);
  M->pL=malloc(sizeof(int)*(n+1));
  M->pU=M->pL;
  M->typefac=FactorizationNO;
  k=0;
  M->pL[0]=0;
  for (i=1;i<n+1;i++)
    { 
      M->pL[i]=M->pL[i-1]+k;
      k = MIN(k+1,band);
    }
  M->Up=malloc(sizeof(double)*M->pL[n]);
  M->Lo=malloc(sizeof(double)*M->pL[n]);
  for (i =0;i<n;i++) M->D[i] =0;
  for (k =0;k<M->pL[n];k++) M->Lo[k] =0;
  for (k =0;k<M->pU[n];k++) M->Up[k] =0;
  return M;
}

/** 
 * pnl_band_matrix_extract
 * extract a band matrix from PnlMat matrix 
 *
 * @param PM adress of a PnlMat
 * @param band number of extra-diagonal in matrix.
 * @return adress of PnlBandMatrix.
 */
PnlBandMatrix* pnl_band_matrix_create_from_full(const PnlMat *PM,int band)
{
  int k,i,n;
  PnlBandMatrix * M = malloc(sizeof(PnlBandMatrix));
  n=PM->n;
  M->n=n;
  M->m=n;/*Square Matrix */
  M->owner=1;
  M->D=malloc(sizeof(double)*n);
  M->pL=malloc(sizeof(int)*(n+1));
  M->pU=M->pL;
  M->typefac=FactorizationNO;
  k=0;
  M->pL[0]=0;
  for (i=1;i<n+1;i++)
    { 
      M->pL[i]=M->pL[i-1]+k;
      k = MIN(k+1,band);
      /* While Row<band size band = row index  */
    }
  M->Up=malloc(sizeof(double)*M->pL[n]);
  M->Lo=malloc(sizeof(double)*M->pL[n]);
  for (i=0;i<n;i++)
    M->D[i] =MGET(PM,i,i);
  i=0;
  for (k =0;k<M->pL[n];k++) 
    {
      if (k>=M->pL[i+1])
        i++;
      M->Lo[k] =MGET(PM,i,i-M->pL[i+1]+k);
      M->Up[k] =MGET(PM,i-M->pU[i+1]+k,i);
      //printf(" M(%d,%d)= %7.4f \n",i,i-M->pL[i+1]+k,M->Lo[k]);
      //printf(" M(%d,%d)= %7.4f \n",i-M->pL[i+1]+k,i,M->Up[k]);
    }
  return M;
}


/** 
 * frees a PnlBandMatrix
 *
 * @param M adress of a PnlBandMatrix.
 * M is set to NULL at exit.
 */
void pnl_band_matrix_free(PnlBandMatrix **M)
{
  if (*M != NULL)
    {
      if ((*M)->owner==1)
      {
        if ((*M)->D  != NULL )
          free((*M)->D);
        
        if ((*M)->Up != NULL )
          free((*M)->Up);
        if ((*M)->Lo != NULL )
          free((*M)->Lo);
        if ((*M)->pU != (*M)->pL )
          free((*M)->pU);
        if ((*M)->pL != NULL )
          free((*M)->pL);
      }
      free(*M);
      *M=NULL;
    }
}
/**
 * clones a PnlBandMatrix
 *
 * @param clone needs to be already allocated.
 * @param v a constant PnlBandMatrixpointer
 */
void pnl_bnd_matrix_clone(PnlBandMatrix * clone,
                          const PnlBandMatrix * v)
{
  if (clone->owner == 0 && clone->n != v->n && clone->m != v->m)
    {
      PNL_ERROR ("owner == 0 and size mismatch", "pnl_vect_clone");
    }
  memcpy(clone->D, v->D, sizeof(double)*v->n);
  memcpy(clone->pL,v->pL,sizeof(int)*(v->n+1));
  memcpy(clone->Up,v->Up,sizeof(double)*v->pL[v->n]);
  memcpy(clone->Lo,v->Lo,sizeof(double)*v->pL[v->n]);
  clone->pU=clone->pL;
  clone->typefac=v->typefac;
}


/** 
 * pnl_band_matrix_copy_create 
 * create a void sparse matrix with pointer on real sparse Matrix
 * d,u,l is set to NULL in some case.
 *
 * @param n int size of diagonal square matrix .
 * @param d adress of a diagonal.
 * @param u adress of a triangular up matrix.
 * @param pu adress of repartition for triangular up matrix.
 * @param l adress of a triangular down matrix.
 * @param pl adress of repartition for triangular low matrix.
 * @return adress of PnlBandMatrix.
 */
static PnlBandMatrix* pnl_band_matrix_cloned_mat(int      n,
                                                 double  *d,
                                                 double  *u,
                                                 int     *pu,
                                                 double  *l,
                                                 int     *pl)
{
  PnlBandMatrix * M = malloc(sizeof(PnlBandMatrix));
  M->n=n;
  M->m=n;/*Square Matrix */
  M->owner=0;
  M->D=d;
  M->Up=u;
  M->Lo=l;
  M->pL=pl;
  M->pU=pu;
  M->typefac=FactorizationNO;
  return M;
}
 



void pnl_bnd_matrix_store_infull(const PnlBandMatrix *BM,PnlMat *M)
{
  int i,k;
  
  pnl_mat_set_double(M,0.0);
  for (i =0;i<BM->n;i++) MLET(M,i,i)=BM->D[i];
  i=0;
  for (k =0;k<BM->pL[BM->n];k++) 
    {
      if (k>=BM->pL[i+1])
        i++;
      MLET(M,i,i-BM->pL[i+1]+k)=BM->Lo[k] ;
      MLET(M,i-BM->pU[i+1]+k,i)=BM->Up[k] ;
    }
}



/**
 * in-place map function
 *
 * @param lhs left hand side PnlBandMatrix
 * @param f the function to be applied term by term
 * @return  lhs = f(lhs)
 */
void pnl_bnd_matrix_map_inplace(PnlBandMatrix *lhs, 
                                double(*f)(double ))
{
  pnl_array_map_inplace(lhs->D,f,lhs->n);
  pnl_array_map_inplace(lhs->Up,f,lhs->pU[lhs->n]);
  pnl_array_map_inplace(lhs->Lo,f,lhs->pL[lhs->n]);
}

/**
 * in-place PnlBandMatrix scalar addition
 *
 * @param lhs left hand side PnlBandMatrix
 * @param x scalar
 * @return  lhs = lhs+x
 */
void pnl_bnd_matrix_plus_double(PnlBandMatrix *lhs , double x)
{
  pnl_array_plus_double(lhs->D,x,lhs->n);    
  pnl_array_plus_double(lhs->Up,x,lhs->pU[lhs->n]);    
  pnl_array_plus_double(lhs->Lo,x,lhs->pL[lhs->n]);    
 
}

/**
 * in-place PnlBandMatrix scalar substraction
 *
 * @param lhs left hand side PnlBandMatrix
 * @param x scalar
 * @return  lhs = lhs-x
 */
void pnl_bnd_matrix_minus_double(PnlBandMatrix *lhs , double x)
{
  pnl_array_minus_double(lhs->D,x,lhs->n);
  pnl_array_minus_double(lhs->Up,x,lhs->pU[lhs->n]);    
  pnl_array_minus_double(lhs->Lo,x,lhs->pL[lhs->n]);    
}

/**
 * in-place PnlBandMatrix scalar multiplication
 *
 * @param lhs left hand side PnlBandMatrix
 * @param x scalar
 * @return  lhs = lhs*x
 */
void pnl_bnd_matrix_mult_double(PnlBandMatrix *lhs , double x)
{
  pnl_array_mult_double(lhs->D,x,lhs->n);
  pnl_array_mult_double(lhs->Up,x,lhs->pU[lhs->n]);    
  pnl_array_mult_double(lhs->Lo,x,lhs->pL[lhs->n]);    
}

/**
 * in-place PnlBandMatrix scalar division
 *
 * @param lhs left hand side PnlBandMatrix
 * @param x scalar
 * @return  lhs = lhs/x
 */
void pnl_bnd_matrix_div_double(PnlBandMatrix *lhs , double x)
{
  pnl_array_div_double(lhs->D,x,lhs->n);
  pnl_array_div_double(lhs->Up,x,lhs->pU[lhs->n]);    
  pnl_array_div_double(lhs->Lo,x,lhs->pL[lhs->n]);    
}


/**
 * map PnlBandMatrix componentwise
 *
 * @param lhs each component lhs(i) contains f(rhs(i))
 * @param rhs right hand side PnlBandMatrix
 * @param f real function 
 */
void pnl_bnd_matrix_map(PnlBandMatrix *lhs, const PnlBandMatrix *rhs, double(*f)(double))
{
  pnl_bnd_matrix_clone(lhs, rhs);
  pnl_bnd_matrix_map_inplace(lhs, f);
}

/**
 * in-place PnlBandMatrix PnlBandMatrix addition
 *
 * @param lhs left hand side PnlBandMatrix
 * @param rhs rigth hand side PnlBandMatrix
 * @return  lhs = lhs+rhs
 */
void pnl_bnd_matrix_plus_mat(PnlBandMatrix *lhs, const PnlBandMatrix *rhs)
{
  pnl_array_plus_array_term(lhs->D, rhs->D,lhs->n);
  pnl_array_plus_array_term(lhs->Up,rhs->Up,lhs->pU[lhs->n]);    
  pnl_array_plus_array_term(lhs->Lo,rhs->Lo,lhs->pL[lhs->n]);    
}

void pnl_bnd_matrix_minus_mat(PnlBandMatrix *lhs, const PnlBandMatrix *rhs)
{
  pnl_array_minus_array_term(lhs->D,rhs->D,lhs->n);
  pnl_array_minus_array_term(lhs->Up,rhs->Up,lhs->pU[lhs->n]);    
  pnl_array_minus_array_term(lhs->Lo,rhs->Lo,lhs->pL[lhs->n]);
 
}

/**
 * in-place term by term PnlBandMatrix inverse
 *
 * @param lhs left hand side PnlBandMatrix
 * @return  lhs = 1 ./ lhs
 */
void pnl_bnd_matrix_inv_term(PnlBandMatrix *lhs)
{
  pnl_array_inv_term(lhs->D,lhs->n);
  pnl_array_inv_term(lhs->Up,lhs->pU[lhs->n]);    
  pnl_array_inv_term(lhs->Lo,lhs->pL[lhs->n]);    
}

/**
 * in-place term by term PnlBandMatrix inverse
 *
 * @param lhs left hand side PnlBandMatrix
 * @param rhs right hand side PnlBandMatrix
 * @return  lhs = lhs ./ rhs
 */
void pnl_bnd_matrix_div_mat_term(PnlBandMatrix *lhs, const PnlBandMatrix *rhs)
{
  pnl_array_div_array_term(lhs->D,rhs->D,lhs->n);
  pnl_array_div_array_term(lhs->Up,rhs->Up,lhs->pU[lhs->n]);    
  pnl_array_div_array_term(lhs->Lo,rhs->Lo,lhs->pL[lhs->n]);
}

/**
 * in-place PnlBandMatrix term by term multiplication
 *
 * @param lhs left hand side PnlBandMatrix
 * @param rhs right hand side PnlBandMatrix
 * @return  lhs = lhs.*rhs
 */
void pnl_bnd_matrix_mult_mat_term(PnlBandMatrix *lhs, const PnlBandMatrix *rhs)
{
  pnl_array_mult_array_term(lhs->D,rhs->D,lhs->n);
 pnl_array_mult_array_term(lhs->Up,rhs->Up,lhs->pU[lhs->n]);    
 pnl_array_mult_array_term(lhs->Lo,rhs->Lo,lhs->pL[lhs->n]);
}

/** 
 * pnl_band_matrix_transpose
 * return adress of \f$ M^{T} \f$
 *
 * @param M adress of PnlBandMatrix
 * @return adress of PnlBandMatrix.
 */
PnlBandMatrix* pnl_band_matrix_transpose(const PnlBandMatrix* M)
{
  return pnl_band_matrix_cloned_mat(M->n,M->D,M->Lo,M->pL,M->Up,M->pU);
}

/** 
 * pnl_band_matrix_low
 * return adress of \f[A_{i,j} = \left\{   M_{i,j} \textit{ if } i < j;\quad 0
 * \textit{ else } \right.
  \f]
 *
 * @param M adress of PnlBandMatrix
 * @return adress of PnlBandMatrix.
 */
PnlBandMatrix* pnl_band_matrix_low(const PnlBandMatrix* M)
{return pnl_band_matrix_cloned_mat(M->n,NULL,NULL,NULL,M->Lo,M->pL);}

/** 
 * pnl_band_matrix_up
 * return adress of \f[
 A_{i,j} = \left\{
 M_{i,j} \textit{ if } i > j \\
 0 \textit{ else }
 \right.
 \f]
 *
 * @param M adress of PnlBandMatrix
 * @return adress of PnlBandMatrix.
 */
PnlBandMatrix* pnl_band_matrix_up(const PnlBandMatrix* M)
{return pnl_band_matrix_cloned_mat(M->n,NULL,M->Up,M->pU,NULL,NULL);}

/** 
 * pnl_band_matrix_tran_low
 * return adress of
 \f[
 A_{i,j} = \left\{
 M_{j,i} \textit{ if }  j < i; \quad
 0 \textit{ else} \right.
 \f]
 *
 * @param M adress of PnlBandMatrix
 * @return adress of PnlBandMatrix.
 */
PnlBandMatrix* pnl_band_matrix_tran_low(const PnlBandMatrix* M)
{return pnl_band_matrix_cloned_mat(M->n,NULL,M->Lo,M->pL,NULL,NULL);}
/** 
 * pnl_band_matrix_tran_up
 * return adress of
 \f[
 A_{i,j} = \left\{
 M_{j,i} \textit{ if } j > i; \quad
  0 \textit{ else }\right.
 \f]
 *
 * @param M adress of PnlBandMatrix
 * @return adress of PnlBandMatrix.
 */
PnlBandMatrix* pnl_band_matrix_tran_up(const PnlBandMatrix* M)
{return pnl_band_matrix_cloned_mat(M->n,NULL,NULL,NULL,M->Up,M->pU);}

/** 
 * pnl_band_matrix_diag
 * return adress of
 \f[
 A_{i,j} = \left\{ M_{i,j} \textit{ if } i = j; \quad  0 \textit{ else }
 \right.
 \f]
 *
 * @param M adress of PnlBandMatrix
 * @return adress of PnlBandMatrix.
 */
PnlBandMatrix* pnl_band_matrix_diag(const PnlBandMatrix* M)
{return pnl_band_matrix_cloned_mat(M->n,M->D,NULL,NULL,NULL,NULL);}

/** 
 * pnl_band_matrix_low_diag
 * return adress of \f$ A = M_{Diag} + M_{Down} \f$
 *
 * @param M adress of PnlBandMatrix
 * @return adress of PnlBandMatrix.
 */
 PnlBandMatrix* pnl_band_matrix_low_diag(const PnlBandMatrix* M)
{return pnl_band_matrix_cloned_mat(M->n,M->D,NULL,NULL,M->Lo,M->pL);}

/** 
 * pnl_band_matrix_up_diag
 * return adress of \f$ A = M_{Diag} + M_{up} \f$
 *
 * @param M adress of PnlBandMatrix
 * @return adress of PnlBandMatrix.
 */
PnlBandMatrix* pnl_band_matrix_up_diag(const PnlBandMatrix* M)
{return pnl_band_matrix_cloned_mat(M->n,M->D,M->Up,M->pU,NULL,NULL);}

/** 
 * pnl_band_matrix_tran_low_diag
 * return adress of \f$ A = M_{Diag} + M_{Low}^T \f$
 *
 * @param M adress of PnlBandMatrix
 * @return adress of PnlBandMatrix.
 */
PnlBandMatrix* pnl_band_matrix_tran_low_diag(const PnlBandMatrix* M)
{return pnl_band_matrix_cloned_mat(M->n,M->D,M->Lo,M->pL,NULL,NULL);}

/** 
 * pnl_band_matrix_tran_up_diag
 * return adress of \f$ A = M_{Diag} + M_{up}^T \f$
 *
 * @param M adress of PnlBandMatrix
 * @return adress of PnlBandMatrix.
 */
PnlBandMatrix* pnl_band_matrix_tran_up_diag(const PnlBandMatrix* M)
{return pnl_band_matrix_cloned_mat(M->n,M->D,NULL,NULL,M->Up,M->pU);}
 

/** 
 * pnl_band_matrix_get_diag
 * @param i row index
 * @param M adress of PnlBandMatrix
 * @return \f$ M(i,i)\f$.
 */
double pnl_band_matrix_get_diag(PnlBandMatrix * M,int i)
{return M->D[i];}

/** 
 * pnl_band_matrix_get
 * @param M adress of PnlBandMatrix
 * @param i row index 
 * @param j col index
 * @return \f$ M(i,j)\f$.
 */
double pnl_band_matrix_get(PnlBandMatrix * M,int i,int j)
{ 
  if( j<i)
    {
      int ij= M->pL[i+1]-i+j;
#ifndef PNL_RANGE_CHECK_OFF
      assert(M->pL[i]<=ij);
#endif
      return M->Lo[ij];
    }
  else if( j>i)
    {
      int ij= M->pU[j+1]-j+i;
#ifndef PNL_RANGE_CHECK_OFF
      assert(M->pU[j]<=ij);
#endif
      return M->Up[ij];
    }
  else
    return M->D[i];
}

/** 
 * pnl_band_matrix_set
 * @param M adress of PnlBandMatrix
 * @param i row index 
 * @param j col index
 * @param x do \f$ M(i,j) = x \f$.
 */
void pnl_band_matrix_set(PnlBandMatrix * M,int i,int j,double x)
{ 
  if(j<i) {
    int ij= M->pL[i+1]-i+j;
#ifndef PNL_RANGE_CHECK_OFF
    assert(M->pL[i]<=ij);
#endif
    M->Lo[ij]=x;
  }
  else if( j>i)
    { int ij= M->pU[j+1]-j+i;
#ifndef PNL_RANGE_CHECK_OFF
      assert(M->pU[j]<=ij);
#endif
      M->Up[ij]=x;
    }
  else
    M->D[i]=x;
}


/** 
 * pnl_band_matrix_set_double
 * put value x for all entries of M,
 *
 * @param M a bandmatrix.
 * @param x double.
 */
void pnl_band_matrix_set_double(PnlBandMatrix*  M,double x)
{
  int k,i;
  for (i =0;i<M->n;i++) M->D[i] =x;
  for (k =0;k<M->pL[M->n];k++) M->Lo[k] =x;
  for (k =0;k<M->pU[M->n];k++) M->Up[k] =x;
}

void pnl_band_matrix_add(PnlBandMatrix * M,int i,int j,double x)
{ 
  if(j<i) {
    int ij= M->pL[i+1]-i+j;
#ifndef PNL_RANGE_CHECK_OFF
    assert(M->pL[i]<=ij);
#endif
    M->Lo[ij]+=x;
  }
  else if( j>i)
    { int ij= M->pU[j+1]-j+i;
#ifndef PNL_RANGE_CHECK_OFF
      assert(M->pU[j]<=ij);
#endif
      M->Up[ij]+=x;
    }
  else
    M->D[i]+=x;
}


/**
 *  in place matrix multiplication
 *
 * @param lhs : vector
 * @param mat  : Band matrix
 * @param rhs : vector
 * @return  lhs=mat*rhs
 */
void pnl_band_matrix_mult_vect_inplace(PnlVect *lhs, const PnlBandMatrix *mat, const PnlVect *rhs)
{
  int i,k;
  for (i =0;i<mat->n;i++)
    LET(lhs,i)=mat->D[i]*GET(rhs,i);
  i=0;
  for (k =0;k<mat->pL[mat->n];k++) 
    {
      if (k>=mat->pL[i+1])
        i++;
      LET(lhs,i)+=mat->Lo[k]*GET(rhs,i-mat->pL[i+1]+k);
      LET(lhs,i-mat->pU[i+1]+k)+=mat->Up[k]*GET(rhs,i);
    }
}

/**
 *  matrix multiplication
 *
 * @param mat : matrix
 * @param vec : vector
 * @return  mat*vec
 */
PnlVect* pnl_band_matrix_mult_vect(const PnlBandMatrix *mat, const PnlVect *vec)
{
  PnlVect *lhs=pnl_vect_create(vec->size);
  pnl_band_matrix_mult_vect_inplace(lhs,mat,vec);
  return lhs;
}


/**
 * compute scalar product <lhs,A rhs >
 *
 * @param lhs : vector
 * @param mat : matrix
 * @param rhs : vector
 * @return  =lhs'*mat*rhs
 */
double pnl_band_matrixprod_scale(const PnlVect *lhs, const PnlBandMatrix *mat,const PnlVect *rhs)
{
  double res;
  PnlVect *tmp=pnl_band_matrix_mult_vect(mat,rhs);
  res=pnl_vect_scalar_prod(lhs,tmp);
  pnl_vect_free(&tmp);
  return res;
}

/**
 * compute y=l A x + b y 
 *
 * @param l  : double
 * @param A : BandMatrix
 * @param x : vector
 * @param b : double
 * @param y : vector
 */
void pnl_band_matrix_lAxpby(double l, const PnlBandMatrix *A, const PnlVect *x, double b, PnlVect * y)
{
  if(l==0.0)
    {
      pnl_vect_mult_double(y,b);
      return ;
    }
  if(b==0.0)
    {
      pnl_band_matrix_mult_vect_inplace(y,A,x);
      pnl_vect_mult_double(y,l);
      return ;
    }
  {
    int i,k;
    for (i =0;i<A->n;i++)
      LET(y,i)=b*LET(y,i)+l*A->D[i]*GET(x,i);
    i=0;
    for (k =0;k<A->pL[A->n];k++) 
      {
        if (k>=A->pL[i+1])
          i++;
        LET(y,i)+=l*A->Lo[k]*GET(x,i-A->pL[i+1]+k);
        LET(y,i-A->pU[i+1]+k)+=l*A->Up[k]*GET(x,i);
      }
  }
}







/**
 * pnl_band_matrix_solve_syslin_inplace
 * solves the linear system M x = x with M PnlBand Matrix.
 *
 * @param M a PnlBandMatrix 
 * @param x right hand side member, use to store solution
 */
void pnl_band_matrix_solve_syslin_inplace(PnlBandMatrix * M, PnlVect *x)
{
  /* -------------------------------------------------------------------- */
  /*   si La diagonal D n'existe pas alors on suppose 1 dessus (cf crout) */
  /* -------------------------------------------------------------------- */
  const double *ij ,*ii, *ik, *ki;
  double *xk,*xi;
  int i;
  double *v, ss;
  int n;  

  v = &LET(x,0);
  n = M->n;  
  if (x->size != n ) 
    {
      printf("  size of matrix & Vector is not compatible");
      abort();
    }
  switch (M->typefac)
    {
    case FactorizationNO:
      if (M->Up!=NULL && M->Lo!=NULL)
        {
          printf("A PROGRAMMER (Not solve if no factorization)");
          abort();
        }
      if (M->Up!=NULL && M->Lo==NULL)
        { /* Up TriDiag Matrix */
          /*  go up  (M->D ? "DU" : "U")  */
          ki = M->Up + M->pU[n]; 
          i = n;
          while ( i-- )
            {
              ii = M->Up + M->pU[i];
              xi= xk  = v +  i ;
              if (M->D!=NULL)
                *xi /= M->D[i];
              /* pour crout ou LU */
              while ( ki > ii) 
                *--xk  -=  *--ki *  *xi ; 
            }
        }
     else if (M->Up==NULL && M->Lo!=NULL)
      { /* Low TriDiag Matrix */
        /* go down  ( M->D ? "LD" : "L" )  */
        ii = M->Lo;
        for (i=0; i<n; i++)
          { ij = ik = (M->Lo + M->pL[i+1]) ;  /* ii =start,ij=end 1 at row  */
            xk = v + i;
            ss = v[i]; 
            while ( ik > ii) 
              ss -= *--ik * *--xk ; 
            if ( M->D) ss /= M->D[i];/* for crout or LU */
            v[i] = ss ;
            ii = ij;
          }
      }
    else if (M->D!=NULL) 
      { /* Diag Matrix */
        for (i=0;i<n;i++) 
          LET(x,i)/=M->D[i];
      }
      break;
    case FactorizationCholeski:
      {
        PnlBandMatrix *D=pnl_band_matrix_low_diag(M);
        PnlBandMatrix *TL=pnl_band_matrix_tran_low_diag(M);
        pnl_band_matrix_solve_syslin_inplace(D,x);
        pnl_band_matrix_solve_syslin_inplace(TL,x);
        pnl_band_matrix_free(&D);
        pnl_band_matrix_free(&TL);
        break;
      }
    case FactorizationCrout:
      {
        PnlBandMatrix *L=pnl_band_matrix_low(M);
        PnlBandMatrix *D=pnl_band_matrix_diag(M);
        PnlBandMatrix *TL=pnl_band_matrix_tran_low_diag(M);
        pnl_band_matrix_solve_syslin_inplace(L,x);
        pnl_band_matrix_solve_syslin_inplace(D,x);
        pnl_band_matrix_solve_syslin_inplace(TL,x);
        pnl_band_matrix_free(&L);
        pnl_band_matrix_free(&D);
        pnl_band_matrix_free(&TL);
        break;
      }
    case FactorizationLU:
      {
        PnlBandMatrix *L=pnl_band_matrix_low(M);
        PnlBandMatrix *U=pnl_band_matrix_up_diag(M);
        pnl_band_matrix_solve_syslin_inplace(L,x);
        pnl_band_matrix_solve_syslin_inplace(U,x);
        pnl_band_matrix_free(&L);
        pnl_band_matrix_free(&U);
        break;
      }
    default:
      PNL_ERROR(" Error unkwon type of Factorization  ","band matrix.c");
    }
}

/**
 * pnl_band_matrix_solve
 * solves the linear system M x = b with M a PnlBand Matrix.
 * For a symmetric definite Choleski or Crout factorization,
 * else us LU factorization
 *
 * @param M a PnlBandMatrix
 * @param x the solution
 * @param b right hand side member
 */
void pnl_band_matrix_solve(PnlBandMatrix * M, PnlVect *x,const PnlVect *b)
{
  if (M->typefac==0)
    switch(M->typefac)
      {
      case FactorizationCholeski: pnl_band_matrix_cholesky(M,EPSILON) ; break;
      case  FactorizationCrout:   pnl_band_matrix_crout(M,EPSILON); break;
      case  FactorizationLU:      pnl_band_matrix_lu(M,EPSILON); break; 
      default:      ; break; 
        
      }
  pnl_vect_clone(x,b);
  pnl_band_matrix_solve_syslin_inplace(M,x);
} 

/**
 * pnl_band_matrix_cholesky 
 * Cholesky factorization M= U^T U of a Band Matrix (Symetric)
 *
 * @param M a PnlBandMatrix
 * @param eps double use for stability of Pivot coefficient
 */
void pnl_band_matrix_cholesky(PnlBandMatrix * M, double eps)
  {
    double  *ij , *ii  , *ik , *jk , xii;
    int i,j,k;
    if (M->Lo != M->Up)
      {
        printf("Choleski factorization of non symetric matrix");
        abort();
      }  
    M->Up = NULL; /*  */
    M->typefac = FactorizationCholeski;
    M->D[0] = sqrt(M->D[0]); 
    ij = M->Lo ; /* pointer on ij coefficient with j<i  */
    for (i=1;i<M->n;i++) /* loop on row  */
      {
        ii = M->Lo+M->pL[i+1];
        /* pointer on last coefficient in row + one =>  ij < ii; */
        xii = M->D[i] ; 
        for ( ; ij < ii ; ij++) /* for j col index of row i */
          {
            j = i -(ii - ij); 
            k = MAX( j - (M->pL[j+1]-M->pL[j]) ,  i-(M->pL[i+1]-M->pL[i]) ); 
            ik =  ii - (i - k); 
            jk =  M->Lo + M->pL[j+1] -(j - k); 
            k = j - k ; 
            while ( k-- )
              *ij -= *ik++ * *jk++;  
            *ij /=  M->D[j] ;
            xii -= *ij * *ij ;
          }
        if (xii < eps*fabs(M->D[i])) 
          {
            printf(" cholesky : pivot ( %d )=  %7.4f < %7.4f ", i, xii ,eps*fabs(M->D[i]));
            abort();
          }
        M->D[i] = sqrt(xii);
    }
  }

/**
 * pnl_band_matrix_crout
 * Crout factorization M= L U of a Band Matrix (Symetric) 
 * \f[
 L=\left(\begin{array}{cccc}
 l_{11}&0&0&0\\
 l_{21}&l_{21}&0&0\\
 l_{31}&l_{32}&l_{33}&0\\
 l_{41}&l_{42}&l_{43}&l_{44}
 \end{array}\right)
 \quad \quad 
 U=\left(\begin{array}{cccc}
 1&u_{12}&u_{13}&u_{14}\\
 0&1&u_{23}&u_{24}\\
 0&0&1&u_{34}\\
 0&0&&1\\
 \end{array}\right)
 \f]
 * @param M a PnlBandMatrix 
 * @param eps double use for stability of Pivot coefficient
 */
void pnl_band_matrix_crout(PnlBandMatrix * M, double eps)
{
  double  *ij , *ii  , *ik , *jk , xii, *dkk;
  int i,j,k;
  if (M->Lo!= M->Up)
      {
        printf("Crout factorization of non symetric matrix");
        abort();
      }
  M->typefac = FactorizationCrout;
  ij = M->Lo; /* pointeur sur le terme ij de la matrice avec j<i  */
  ij = M->Lo; /* pointer on ij coefficient with j<i  */
  for (i=1;i<M->n;i++) /* loop on row  */
    {
      ii = M->Lo+M->pL[i+1];
      /* pointer on last coefficient in row + one =>  ij < ii; */
      xii = M->D[i] ; 
      for ( ; ij < ii ; ij++) /* for j col index of row i */
        {
          j = i -(ii - ij); 
          k = MAX( j - (M->pL[j+1]-M->pL[j]) ,  i-(M->pL[i+1]-M->pL[i]) ); 
          ik =  ii - (i - k); 
          jk =  M->Lo+ M->pL[j+1] -(j - k); 
          dkk = M->D + k;
          k = j - k ; 
          while ( k-- )
            *ij -= *ik++ * *jk++ * *dkk++;  
          *ij /=  *dkk ; /* k = j ici  */
          xii -= *ij * *ij * *dkk;
        }
    if (fabs(xii) <= MAX(eps*fabs(M->D[i]),1.0e-30))
      {
        printf(" crout : pivot ( %d )=  %7.4f < %7.4f ", i, xii ,eps*fabs(M->D[i]));
        abort();
      }
	M->D[i] = xii;
    }
}

/**
 * pnl_band_matrix_crout
 * Crout factorization M= L U of a Band Matrix (No Symetric) 
 * @param M a PnlBandMatrix 
 * @param eps double use for stability of Pivot coefficient
 */
void pnl_band_matrix_lu(PnlBandMatrix * M, double eps)
{
  double s,uii;
  int i,j,k, k0, j0;
  double *Lik, *Ukj, *Ljk, *Uki;


  if (M->Lo== M->Up &&
      ( M->pL[M->n]  || M->pU[M->n] ) )
    {
      printf(" in LU , symetric matrix ");
      abort();
    }
  M->typefac=FactorizationLU;
  for (i=1;i<M->n;i++)
    /* loop on i rank sub-matrix   */
    { 
      /* for L(i,j)  j=j0,i-1 */
      j0 = i-(M->pL[i+1]-M->pL[i]);
      for ( j = j0; j<i;j++)
        {           
          k0 = MAX(j0,j-(M->pU[j+1]-M->pU[j]));
          Lik = M->Lo+ M->pL[i+1]-i+k0; /* lower */
          Ukj = M->Up + M->pU[j+1]-j+k0; /* upper */
          s =0;
          for (k=k0;k<j;k++) /* k < j < i ; */
            s += *Lik++ * *Ukj++ ;     /* a(i,k)*a(k,j); */
          *Lik -= s;
          *Lik /= M->D[j]; /*  k == j here */
        }
      /* for U(j,i) j=0,i-1         */
      j0=i-M->pU[i+1]+M->pU[i];
      for (j=j0;j<i;j++) 
        {
          s = 0;
          k0 = MAX(j0,j-M->pL[j+1]+M->pL[j]);
          Ljk = M->Lo+ M->pL[j+1]-j+k0;   
          Uki = M->Up + M->pU[i+1]-i+k0;   
          for (k=k0  ;k<j;k++)    /*  */
            s +=  *Ljk++ * *Uki++ ;
          *Uki -= s;  /* k = j here  */
        }
      /* for D (i,i) in last because we need L(i,k) and U(k,i) for k<j */
      k0 = i-MIN(M->pL[i+1]-M->pL[i],M->pU[i+1]-M->pU[i]);
      Lik = M->Lo+ M->pL[i+1]-i+k0; /* lower */
      Uki = M->Up + M->pU[i+1]-i+k0; /* upper */
      s =0;
      for (k=k0;k<i;k++) /* k < i < i ; */
        s += *Lik++ * *Uki++ ;     /* a(i,k)*a(k,i); */
      //printf(" k0 %d, i = %d -> %7.4f \n",k0,i,s);
      uii = M->D[i]-s;
      
      if (fabs(uii) <= MAX(eps*fabs(M->D[i]),1.0e-30))
        {
          printf(" LU : pivot ( %d )=  %7.4f < %7.4f < %7.4f", i,fabs(uii) ,eps*fabs(M->D[i]),eps);
          abort();
        }
      M->D[i] = uii;
      
    }
}

double pnl_band_matrix_conditionning(const PnlBandMatrix *M)
{
  double maxi,mini;
  PnlVect V=pnl_vect_create_wrap_array(M->D,MAX(M->n,M->m));
  pnl_vect_minmax (&V,&mini,&maxi);
  return fabs(mini)/fabs(maxi);
};

