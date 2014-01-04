
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


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <stdarg.h>

#include "pnl/pnl_config.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_internals.h"

static char pnl_matrix_label[] = "PnlMatObject";
static char pnl_hmatrix_label[] = "PnlHmatObject";

/**
 * Create a PnlMatObject which is the parent type of all matrices
 */
PnlMatObject* pnl_mat_object_new ()
{
  PnlMatObject *o;
  if ( (o = malloc (sizeof (PnlMatObject))) == NULL) return NULL;
  o->m = 0;
  o->n = 0;
  o->mn = 0;
  o->mem_size = 0;
  o->owner = 1;
  o->array = NULL;
  o->object.type = PNL_TYPE_MATRIX;
  o->object.parent_type = PNL_TYPE_MATRIX;
  o->object.nref = 0;
  o->object.label = pnl_matrix_label;
  o->object.destroy = (DestroyFunc *) pnl_mat_object_free;
  o->object.constructor = (NewFunc *) pnl_mat_object_new;
  o->object.clone = NULL;
  o->object.copy = NULL;
  return o;
}

/**
 * Free a PnlMatObject
 */
void pnl_mat_object_free (PnlMatObject **o)
{
    if (*o != NULL)
    {
      if ((*o)->array != NULL && (*o)->owner == 1)
        { free((*o)->array); (*o)->array = NULL; }
      (*o)->m = (*o)->n = (*o)->mn = 0;
      free(*o);
    }
}

/**
 * Resize a PnlMatObject.  If the new size is smaller than the current one, no
 * memory is freed. If the new size is larger than the current mem_size, a new
 * pointer is allocated. The old data are kept.
 *
 * @param M : a pointer to an already existing PnlMatObject 
 * @param m : new nb rows
 * @param n : new nb columns
 *
 * @return OK or FAIL. When returns OK, the matrix is changed. 
 */
int pnl_mat_object_resize(PnlMatObject *M, int m, int n)
{
  int mn;
  size_t sizeof_base = 0;
  mn = m*n;
  if (M->owner == 0) return OK;
  if (mn < 0) return FAIL;
  if (mn == 0) /* free array */
    {
      M->m = M->n = M->mn = M->mem_size = 0;
      if (M->array != NULL) { free(M->array); M->array = NULL; }
      return OK;
    }

  if (M->mem_size >= mn) 
    {
      /* If the new size is smaller, we do not reduce the size of the
         allocated block. It may change, but it may allow to grow the matrix
         quicker */
      M->m=m; M->n=n; M->mn = mn;
      return OK;
    }

  /* Now, M->mem_size < mn */
  M->m = m; M->n = n;
  M->mn = M->mem_size = mn;
  switch (PNL_GET_TYPE(M))
    {
    case PNL_TYPE_MATRIX_DOUBLE : sizeof_base = sizeof(double); break;
    case PNL_TYPE_MATRIX_COMPLEX : sizeof_base = sizeof(dcomplex); break;
    case PNL_TYPE_MATRIX_INT : sizeof_base = sizeof(int); break;
    default : PNL_MESSAGE (1, "Unknown type in pnl_mat_object_resize.\n");
      return FAIL;
    }
  if ((M->array=realloc(M->array,mn*sizeof_base)) == NULL) return FAIL;
  return OK;
}

/**
 * Create a PnlHmatObject which is the parent type of all hyper-matrices
 */
PnlHmatObject* pnl_hmat_object_new ()
{
  PnlHmatObject *o;
  if ( (o = malloc (sizeof (PnlHmatObject))) == NULL) return NULL;
  o->ndim = 0;
  o->dims = NULL;
  o->pdims = NULL;
  o->mn = 0;
  o->array = NULL;
  o->object.type = PNL_TYPE_HMATRIX;
  o->object.parent_type = PNL_TYPE_HMATRIX;
  o->object.nref = 0;
  o->object.label = pnl_hmatrix_label;
  o->object.constructor = (NewFunc *) pnl_hmat_object_new;
  o->object.clone = NULL;
  o->object.copy = NULL;
  return o;
}

/**
 * Free a Hmat object
 */
void pnl_hmat_object_free(PnlHmatObject **H)
{
  if (*H != NULL)
    {
      free((*H)->array);
      free((*H)->dims);
      free((*H)->pdims);
      free(*H);
      *H=NULL;
    }
}

/** 
 * Compute the array of leading dimension for Hmat
 * 
 * @param pdims array of size ndim
 * @param ndim number of dimensions of the Hmat
 * @param dims array of the sizes of the Hmat
 */
void pnl_hmat_compute_pdims (int *pdims, int ndim, const int *dims)
{
  int i;
  pdims[ndim-1] = 1.;
  for ( i=ndim-1 ; i>0 ; i-- )
    {
      pdims[i-1] = dims[i] * pdims[i];
    }
}

/** 
 * Compute the value of the linear index (ie. the index in the field
 * array) corresponding to the multi-index defined by tab
 * 
 * @param H a Hmat
 * @param tab multi-index
 * 
 * @return an integer index
 */
int pnl_hmat_compute_linear_index (PnlHmatObject *H, int *tab)
{
  int i, index;
  index = 0;

  for ( i=0 ; i<H->ndim ; i++ ) { index += tab[i] * H->pdims[i]; }
  return index;
}

/**
 * Resize a PnlHmatObject.
 *
 * If the new size is smaller than the current one, no
 * memory is free. If the new size is larger than the
 * current one, more space is allocated. Note that for the
 * sake of efficiency the old data are not copied.
 *
 * @param H : a pointer to an already existing PnlHmatObject
 * @param ndim : new nb dimensions
 * @param dims : new pointer to the dimensions array
 *
 * @return OK or FAIL. When returns OK, the hmatrix is changed. 
 */
int pnl_hmat_object_resize(PnlHmatObject *H, int ndim, const int *dims)
{
  int i, s=1;
  size_t sizeof_base = 0;
  for ( i=0 ; i<ndim ;i++ ) { s *= dims[i]; }

  if (H->mn == s) /* nothing to do, just adjust ndim, dims and pdims */
    {
      H->ndim=ndim;
      if (H->ndim > ndim) 
        {
          if ((H->dims=realloc(H->dims,sizeof(int)*ndim))==NULL) return FAIL;
          if ((H->pdims=realloc(H->pdims,sizeof(int)*ndim))==NULL) return FAIL;
        }
      memcpy(H->dims, dims, ndim*sizeof(int));
      pnl_hmat_compute_pdims (H->pdims, ndim, dims);
      return OK;
    }
  if (s< 0) return FAIL;
  if (s==0) /* free array */
    {
      H->ndim =  H->mn = 0;
      H->dims=NULL;
      H->pdims=NULL;
      free(H->array); H->array = NULL;
      return OK;
    }
  H->ndim=ndim; H->mn=s;
  if ((H->dims=realloc(H->dims,sizeof(int)*ndim))==NULL) return FAIL;
  memcpy(H->dims, dims, ndim*sizeof(int));
  if ((H->pdims=realloc(H->pdims,sizeof(int)*ndim))==NULL) return FAIL;
  pnl_hmat_compute_pdims (H->pdims, ndim, dims);

  switch (PNL_GET_TYPE (H))
    {
    case PNL_TYPE_HMATRIX_DOUBLE : sizeof_base = sizeof(double);
      break;
    case PNL_TYPE_HMATRIX_COMPLEX : sizeof_base = sizeof(dcomplex);
      break;
    case PNL_TYPE_HMATRIX_INT : sizeof_base = sizeof(int);
      break;
    default : PNL_MESSAGE (1, "Unknown type in pnl_hmat_object_resize.\n");
      return FAIL;
    }
  if (H->array==NULL)
    {
      if ((H->array = malloc(H->mn*sizeof_base))==NULL)
        return FAIL;
    }else
    {
      if ((H->array = realloc(H->array,H->mn*sizeof_base))==NULL)
        return FAIL;
    }
  return OK;
}


typedef int(*cmp_func)(const void *, const void *);

#define BASE_DOUBLE
#include "pnl/pnl_templates_on.h"
static char pnl_mat_label[] = "PnlMatrixDouble";
static char pnl_hmat_label[] = "PnlHmatrixDouble";
#include "matrix_source.c"
#include "dgemm_source.c"
#include "pnl/pnl_templates_off.h"
#undef  BASE_DOUBLE


#define BASE_PNL_COMPLEX
#include "pnl/pnl_templates_on.h"
static char pnl_mat_complex_label[] = "PnlMatrixComplex";
static char pnl_hmat_complex_label[] = "PnlHmatrixComplex";
#include "matrix_source.c"
#include "dgemm_source.c"
#include "pnl/pnl_templates_off.h"
#undef  BASE_PNL_COMPLEX

#define BASE_INT
#include "pnl/pnl_templates_on.h"
static char pnl_mat_int_label[] = "PnlMatrixInt";
static char pnl_hmat_int_label[] = "PnlHmatrixInt";
#include "matrix_source.c"
#include "dgemm_source.c"
#include "pnl/pnl_templates_off.h"
#undef  BASE_INT

/** 
 * Compute the cross product of the vectors in A and B.
 * First, a row wise computation is tried; if it fails, the column wise
 * approach is tested
 * 
 * @param[out] lhs Contains (A) x (B) on output
 * @param A must be a matrix with 3 rows or columns
 * @param B must be a matrix with the same size as A
 * 
 * @return FAIL in case of dimension mismatch, OK otherwise
 */
int pnl_mat_cross(PnlMat *lhs, const PnlMat *A, const PnlMat *B)
{
  CheckMatMatch (A, B)
  pnl_mat_resize (lhs, A->m, A->n);  

  if ( A->n == 3 )
    {
      int i;
      for ( i=0 ; i<A->m ; i++ )
        {
          MLET(lhs,i,0) = MGET(A,i,1) * MGET(B,i,2) - MGET(A,i,2) * MGET(B,i,1);
          MLET(lhs,i,1) = MGET(A,i,2) * MGET(B,i,0) - MGET(A,i,0) * MGET(B,i,2);
          MLET(lhs,i,2) = MGET(A,i,0) * MGET(B,i,1) - MGET(A,i,1) * MGET(B,i,0);
        }
      return OK;
    }
  if ( A->m == 3 )
    {
      int j;
      for ( j=0 ; j<A->n ; j++ )
        {
          MLET(lhs,0,j) = MGET(A,1,j) * MGET(B,2,j) - MGET(A,2,j) * MGET(B,1,j);
          MLET(lhs,1,j) = MGET(A,2,j) * MGET(B,0,j) - MGET(A,0,j) * MGET(B,2,j);
          MLET(lhs,2,j) = MGET(A,0,j) * MGET(B,1,j) - MGET(A,1,j) * MGET(B,0,j);
        }
      return OK;
    }
  else
    return FAIL;
}

/**
 * Matrix exponential B = exp( A)
 *
 * @param A a matrix
 * @param B contains exp(A) on return
 *
 * This algorithms has been written in C by Jérôme Lelong
 * following the Fortran implementation given in EXPOKIT
 * Roger B. Sidje Department of Mathematics, University of Queensland 
 * Brisbane, QLD-4072, Australia, (c) 1996-1999 All Rights Reserved
 * The work resulting from EXPOKIT has been published in ACM-Transactions 
 * on Mathematical Software, 24(1):130-156, 1998.
 *
 *  The bibtex record of the citation:
 *
 * \verbatim
 * ARTICLE{EXPOKIT,
 *        AUTHOR  = {Sidje, R. B.},
 *        TITLE   = {{\sc Expokit.} {S}oftware {P}ackage for
 *                 {C}omputing {M}atrix {E}xponentials},
 *        JOURNAL = {ACM Trans. Math. Softw.},
 *        VOLUME  = {24},
 *        NUMBER  = {1},
 *        PAGES   = {130-156}
 *        YEAR    = {1998}
 * }
 * \endverbatim
 */
int pnl_mat_exp (PnlMat *B, const PnlMat *A)
{
  int i,j,k,ih2,ip,iq,iused,ifree,iodd,icoef,ideg,iput,iget,lwork, ns;
  double *work;
  PnlMat Mwork, Mwork1, Mwork2;
  double hnorm,scale,scale2,cp,cq;
  
  CheckIsSquare (A);
  
  /* Mwork is used as a container for pnl_mat_xxx routines */
  pnl_mat_resize (B, A->m, A->n);
  
  icoef = 0;
  ideg = 6;
  ih2 = icoef + (ideg+1);
  ip  = ih2 + A->mn;
  iq  = ip + A->mn;
  ifree = iq + A->mn;
  lwork = 4*A->m*A->m + ideg + 1;
  if ( (work = malloc (lwork * sizeof(double))) == NULL ) abort();

  /*  scaling: seek ns such that ||t*H/2^ns|| < 1/2;  */
  /*  and set scale = t/2^ns ... */
  for (i=0; i<A->m; i++) work[i] = 0.;
  for (j=0; j<A->m; j++)
    {
      for (i=0; i<A->m; i++)
        work[i] += fabs( PNL_MGET(A,i,j) );      
    }

  hnorm = 0.;
  for (i=0; i<A->m; i++)
    {
      hnorm = MAX( hnorm, work[i] );
    }

  if (hnorm == 0.)
    /* matrix is full of zeros */
    {
      pnl_mat_set_id (B);
      free (work); work = NULL;
      return FAIL;
    }
  ns = MAX( 0, (int)(log(hnorm)/log(2.)) + 2 );
  scale = 1. / pnl_pow_i (2., ns);
  scale2 = scale*scale;

  /*  compute Pade coefficients ... */
  i = ideg+1;
  j = 2*ideg+1;
  work[icoef] = 1.0;
  for (k=1; k<ideg+1; k++)
    {
      work[icoef+k] = (work[icoef+k-1]*( i-k )) / ( k*(j-k) );
    }

  Mwork = pnl_mat_wrap_array (&(work[ih2]), A->m, A->n);
  pnl_mat_dgemm ('N', 'N', scale2, A, A, 0., &Mwork);  /* H2 = scale2*H*H */
  
  /* initialize p (numerator) and q (denominator) */
  cp = work[icoef+ideg-1];
  cq = work[icoef+ideg];
  for (j=0; j<A->m; j++)
    {
      for (i=0; i<A->m; i++)
        {
          work[ip + j*A->m + i] = 0.;
          work[iq + j*A->m + i] = 0.;
        }
      work[ip + j*(A->m+1)] = cp; /* sets the diagonal */
      work[iq + j*(A->m+1)] = cq; /* sets the diagonal */
    }

  /* Apply Horner rule */
  iodd = 1;
  k = ideg - 1;
  while (k > 0)
    {
      iused = iodd*iq + (1-iodd)*ip;
      Mwork = pnl_mat_wrap_array (&(work[ifree]), A->m, A->m);
      Mwork1 = pnl_mat_wrap_array (&(work[iused]), A->m, A->m);
      Mwork2 = pnl_mat_wrap_array (&(work[ih2]), A->m, A->m);
      pnl_mat_dgemm('N', 'N', 1., &Mwork1, &Mwork2, 0., &Mwork);
      
      for (j=0; j<A->m; j++)
        {
          /* add work[icoef+k-1]; to the diagonal */
          work[ifree+j*(A->m+1)] += work[icoef+k-1];
        }
      ip = (1-iodd)*ifree + iodd*ip;
      iq = iodd*ifree + (1-iodd)*iq;
      ifree = iused;
      iodd = 1-iodd;
      k--;
    }
 
  /* Obtain (+-)(I + 2 * (p/q))  */
  Mwork = pnl_mat_wrap_array (&(work[ifree]), A->m, A->m);
  Mwork1 = pnl_mat_wrap_array (&(work[ip]), A->m, A->m);
  pnl_mat_dgemm ('N', 'N', scale, &Mwork1, A, 0., &Mwork);
  ip = ifree;

  Mwork1 = pnl_mat_wrap_array (&(work[ip]), A->m, A->m);
  Mwork2 = pnl_mat_wrap_array (&(work[iq]), A->m, A->m);

  pnl_mat_axpy (-1., &Mwork1, &Mwork2 );
  pnl_mat_syslin_mat (&Mwork2, &Mwork1);

  pnl_mat_mult_scalar (&Mwork1, 2.);
  for (j=0; j<A->m; j++) work[ip+j*(A->m+1)] += 1.;
  iput = ip;

  if (ns == 0 && iodd == 1)
    {
      pnl_mat_mult_scalar (&Mwork1, -1.);
    }
  else
    {
      /* squaring : exp(t*H) = (exp(t*H))^(2^ns) */
      iodd = 1;
      for (k=0; k<ns; k++)
        {
          iget = iodd*ip + (1-iodd)*iq;
          iput = (1-iodd)*ip + iodd*iq;
          Mwork = pnl_mat_wrap_array (&(work[iget]), A->m, A->m);
          Mwork1 = pnl_mat_wrap_array (&(work[iput]), A->m, A->m);
          pnl_mat_dgemm ('N', 'N', 1., &Mwork, &Mwork, 0., &Mwork1);
          iodd = 1-iodd;
        }
    }

  /* the solution is located at work[iput] */
  memcpy (B->array, &(work[iput]), A->mn * sizeof(double));
  
  free (work); work = NULL;
  return OK;
}

int pnl_mat_complex_exp (PnlMatComplex *B, const PnlMatComplex *A)
{
  int i,j,k,ih2,ip,iq,iused,ifree,iodd,icoef,ideg,iput,iget,lwork, ns;
  dcomplex *work;
  PnlMatComplex Mwork, Mwork1, Mwork2;
  double hnorm;
  dcomplex scale,scale2,cp,cq;
  
  CheckIsSquare (A);
  
  /* Mwork is used as a container for pnl_mat_xxx routines */
  pnl_mat_complex_resize (B, A->m, A->n);
  
  icoef = 0;
  ideg = 6;
  ih2 = icoef + (ideg+1);
  ip  = ih2 + A->mn;
  iq  = ip + A->mn;
  ifree = iq + A->mn;
  lwork = 4*A->m*A->m + ideg + 1;
  if ( (work = malloc (lwork * sizeof(dcomplex))) == NULL ) abort();

  /*  scaling: seek ns such that ||t*H/2^ns|| < 1/2;  */
  /*  and set scale = t/2^ns ... */
  for (i=0; i<A->m; i++) work[i] = CZERO;
  for (j=0; j<A->m; j++)
    {
      for (i=0; i<A->m; i++)
        work[i].r += Cabs ( PNL_MGET(A,i,j) );      
    }

  hnorm = 0.;
  for (i=0; i<A->m; i++)
    {
      hnorm = MAX( hnorm, work[i].r );
    }

  if (hnorm == 0.)
    /* matrix is full of zeros */
    {
      pnl_mat_complex_set_id (B);
      free (work); work = NULL;
      return FAIL;
    }
  ns = MAX( 0, (int)(log(hnorm)/log(2.)) + 2 );
  scale = Complex(1. / pnl_pow_i (2., ns), 0.);
  scale2 = Cmul(scale,scale);

  /*  compute Pade coefficients ... */
  i = ideg+1;
  j = 2*ideg+1;
  work[icoef] = CONE;
  for (k=1; k<ideg+1; k++)
    {
      work[icoef+k] = CRdiv( CRmul( work[icoef+k-1], ( i-k ) ), k*(j-k) );
    }

  Mwork = pnl_mat_complex_wrap_array (&(work[ih2]), A->m, A->n);
  pnl_mat_complex_dgemm ('N', 'N', scale2, A, A, CZERO, &Mwork);  /* H2 = scale2*H*H */
  
  /* initialize p (numerator) and q (denominator) */
  cp = work[icoef+ideg-1];
  cq = work[icoef+ideg];
  for (j=0; j<A->m; j++)
    {
      for (i=0; i<A->m; i++)
        {
          work[ip + j*A->m + i] = CZERO;
          work[iq + j*A->m + i] = CZERO;
        }
      work[ip + j*(A->m+1)] = cp; /* sets the diagonal */
      work[iq + j*(A->m+1)] = cq; /* sets the diagonal */
    }

  /* Apply Horner rule */
  iodd = 1;
  k = ideg - 1;
  while (k > 0)
    {
      iused = iodd*iq + (1-iodd)*ip;
      Mwork = pnl_mat_complex_wrap_array (&(work[ifree]), A->m, A->m);
      Mwork1 = pnl_mat_complex_wrap_array (&(work[iused]), A->m, A->m);
      Mwork2 = pnl_mat_complex_wrap_array (&(work[ih2]), A->m, A->m);
      pnl_mat_complex_dgemm('N', 'N', CONE, &Mwork1, &Mwork2, CZERO, &Mwork);
      
      for (j=0; j<A->m; j++)
        {
          /* add work[icoef+k-1]; to the diagonal */
          work[ifree+j*(A->m+1)].r += work[icoef+k-1].r;
          work[ifree+j*(A->m+1)].i += work[icoef+k-1].i;
        }
      ip = (1-iodd)*ifree + iodd*ip;
      iq = iodd*ifree + (1-iodd)*iq;
      ifree = iused;
      iodd = 1-iodd;
      k--;
    }
 
  /* Obtain (+-)(I + 2 * (p/q))  */
  Mwork = pnl_mat_complex_wrap_array (&(work[ifree]), A->m, A->m);
  Mwork1 = pnl_mat_complex_wrap_array (&(work[ip]), A->m, A->m);
  pnl_mat_complex_dgemm ('N', 'N', scale, &Mwork1, A, CZERO, &Mwork);
  ip = ifree;

  Mwork1 = pnl_mat_complex_wrap_array (&(work[ip]), A->m, A->m);
  Mwork2 = pnl_mat_complex_wrap_array (&(work[iq]), A->m, A->m);

  pnl_mat_complex_axpy (Complex(-1.,0.), &Mwork1, &Mwork2 );
  pnl_mat_complex_syslin_mat (&Mwork2, &Mwork1);

  pnl_mat_complex_mult_scalar (&Mwork1, Complex(2.,0.));
  for (j=0; j<A->m; j++) work[ip+j*(A->m+1)].r += 1.;
  iput = ip;

  if (ns == 0 && iodd == 1)
    {
      pnl_mat_complex_mult_scalar (&Mwork1, Complex(-1.,0.));
    }
  else
    {
      /* squaring : exp(t*H) = (exp(t*H))^(2^ns) */
      iodd = 1;
      for (k=0; k<ns; k++)
        {
          iget = iodd*ip + (1-iodd)*iq;
          iput = (1-iodd)*ip + iodd*iq;
          Mwork = pnl_mat_complex_wrap_array (&(work[iget]), A->m, A->m);
          Mwork1 = pnl_mat_complex_wrap_array (&(work[iput]), A->m, A->m);
          pnl_mat_complex_dgemm ('N', 'N', CONE, &Mwork, &Mwork, CZERO, &Mwork1);
          iodd = 1-iodd;
        }
    }

  /* the solution is located at work[iput] */
  memcpy (B->array, &(work[iput]), A->mn * sizeof(dcomplex));
  
  free (work); work = NULL;
  return OK;
}

/** 
 * Create a complex matrix from a real one.
 * 
 * @param R a real matrix
 * 
 * @return a c xomplexified version of R
 */
PnlMatComplex* pnl_mat_complex_create_from_mat (const PnlMat *R)
{
  int i;
  PnlMatComplex *C;
  C = pnl_mat_complex_create (R->m, R->n);

  for ( i=0 ; i<C->mn ; i++ )
    {
      C->array[i].r = R->array[i];
      C->array[i].i = 0.;
    }
  return C; 
}


