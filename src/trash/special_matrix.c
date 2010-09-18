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
/* - JP Chancelier & al in NSP software                                  */
/*     http://cermics.enpc.fr/~jpc/nsp-tiddly/mine.html                  */ 
/* - F.Hecht & al in RMN class distributed on Freefem Project            */
/*     http://www.freefem.org/                                           */ 
/*************************************************************************/

  

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
  
#include "pnl/pnl_spec_matrix.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_array.h"

#define GETROWORCOLSIZE(M) ((M)->RC) ? ((M)->n): ((M)->m)
#define GETROWORCOL(M,i,j) ((M)->RC) ? &((M)->array[(j)]): &((M)->array[(i)])
#define GETROWORCOLIND(M,i,j) ((M)->RC) ? (i) : (j)

static void pnl_sprow_init(SpRow *Row,int Max_size)
{
  if(Max_size<0)
    Row = NULL;
  Row->size=0;
  Row->Max_size=Max_size;
  if((Row->Value=malloc(Row->Max_size*sizeof(double)))==NULL)
    Row= NULL;
  if((Row->Index=malloc(Row->Max_size*sizeof(int)))==NULL)
    Row= NULL;
 }

static void pnl_sprow_free(SpRow *row)
{
  if(row!=NULL)
    {
      free(row->Index);
      free(row->Value);
    }
}

static int pnl_sprow_search(SpRow * Row,int j)
{
  int k=0;
  while(k<Row->size && Row->Index[k]!=j){k++;}
  return k;
}

static int pnl_sprow_freeze(SpRow * Row)
{
  if(Row->size==Row->Max_size)
    return OK;
  Row->Max_size=Row->size;
  {
    double * New_Value;
    int * New_Index;
    if((New_Value=malloc(Row->Max_size*sizeof(double)))==NULL)
      return FAIL;
    memcpy (New_Value,Row->Value,Row->Max_size*sizeof(double));
    if((New_Index=malloc(Row->Max_size*sizeof(int)))==NULL)
      return FAIL;
    memcpy (New_Index,Row->Index,Row->Max_size*sizeof(int));
    free(Row->Value);
    free(Row->Index);
    Row->Value=&(New_Value[0]);
    Row->Index=&(New_Index[0]);
  }
  return OK;
}


static int pnl_sprow_set(SpRow * Row,int j,double Val)
{
  int k=pnl_sprow_search(Row,j);
  if(k<Row->size)
    {
      Row->Value[k]=Val;
      return OK;
    }
  if(k<Row->Max_size)
    {
      Row->Index[k]=j;
      Row->Value[k]=Val;
      Row->size++;
      return OK;
    }
  // Reallocation.
  Row->size++;
  Row->Max_size+=10;
  {
    double * New_Value=Row->Value;
    int * New_Index=Row->Index;
    if((Row->Value=malloc(Row->Max_size*sizeof(double)))==NULL)
      return FAIL;
    memcpy (Row->Value,New_Value,Row->size*sizeof(double));
    if((Row->Index=malloc(Row->Max_size*sizeof(int)))==NULL)
      return FAIL;
    memcpy (Row->Index,New_Index,Row->size*sizeof(int));
    free(New_Value);
    free(New_Index);
  }
  Row->Index[k]=j;
  Row->Value[k]=Val;
  return OK;
}

/*
  static double pnl_sprow_get(SpRow * Row,int j)
{
  int k=pnl_sprow_search(Row,j);
  if(k<Row->size)
    {
      return  Row->Value[k];
    }
  return 0.0;
}
*/

static double * pnl_sprow_lget(SpRow * Row,int j)
{
  int k=pnl_sprow_search(Row,j);
  if(k<Row->size)
    {
      return  &(Row->Value[k]);
    }
  if(k<Row->Max_size)
    {
      Row->Index[k]=j;
      Row->size++;
      return &(Row->Value[k]);
    }
  // Reallocation.
  Row->Max_size+=10;
  {
    double * New_Value=Row->Value;
    int * New_Index=Row->Index;
    if((Row->Value=malloc(Row->Max_size*sizeof(double)))==NULL)
      {
        PNL_ERROR(" Error in reallocation ", "pnl_sprow_lget");
        return NULL;
      }
    memcpy (Row->Value,New_Value,Row->size*sizeof(double));
    if((Row->Index=malloc(Row->Max_size*sizeof(int)))==NULL)
      {
        PNL_ERROR(" Error in reallocation ", "pnl_sprow_lget");
        return NULL;
      }
    memcpy (Row->Index,New_Index,Row->size*sizeof(int));
    free(New_Value);
    free(New_Index);
  }
  Row->Index[Row->size]=j;
  Row->size++;
  return &(Row->Value[k]);
  
}


static int pnl_sprow_add(SpRow * Row,int j,double Val)
{
  if(Row->size<Row->Max_size)
    {
      Row->Index[Row->size]=j;
      Row->Value[Row->size]=Val;
      Row->size++;
      return OK;
    }
  // Reallocation.
  Row->Max_size+=10;
  {
    double * New_Value=Row->Value;
    int * New_Index=Row->Index;
    if((Row->Value=malloc(Row->Max_size*sizeof(double)))==NULL)
      return FAIL;
    memcpy (Row->Value,New_Value,Row->size*sizeof(double));
       if((Row->Index=malloc(Row->Max_size*sizeof(int)))==NULL)
      return FAIL;
    memcpy (Row->Index,New_Index,Row->size*sizeof(int));
    free(New_Value);
    free(New_Index);
  }
  Row->Index[Row->size]=j;
  Row->Value[Row->size]=Val;
  Row->size++;
  return OK;
}

/**
 * creates a PnlMorseMat
 * @param m number of rows
 * @param n number of rowumns
 * @param Max_size_row Allocation for row, same for each
 * @param RC int store in row or col
 * @return a PnlMorseMat pointer
 */
PnlMorseMat * pnl_morse_mat_create(int m, int n,int Max_size_row,int RC)
{
  int size;
  PnlMorseMat *Sp;
  if((Sp=malloc(sizeof(PnlMorseMat)))==NULL)
    return NULL;
  Sp->m=m;
  Sp->n=n;
  Sp->RC=RC;
  size=(RC)?Sp->n:Sp->m;
  if (size>0)
    {
      if((Sp->array=malloc(size*sizeof(SpRow)))==NULL) return NULL;
      if(Max_size_row!=0)
        {
          int i=0;
          while(i<size)
            {
              pnl_sprow_init(&(Sp->array[i]),Max_size_row);
              i++;
            }
        }
    }
  else
    Sp->array=(SpRow*) NULL;
  return Sp;
}

/**
 * creates a PnlMorseMat
 * @param m number of rows
 * @param n number of rowumns
 * @param Max_size_row Allocation for row, 
 * @param RC int store in row or col
 * @return a PnlMorseMat pointer
 */
PnlMorseMat * pnl_morse_mat_create_with_size(int m, int n,PnlVectInt *Max_size_row,int RC)
{
  int size;
  PnlMorseMat *Sp;
  if((Sp=malloc(sizeof(PnlMorseMat)))==NULL)
    return NULL;
  Sp->m=m;
  Sp->n=n;
  Sp->RC=RC;
  size=(RC)?Sp->n:Sp->m;
  if (size>0)
    {
      int i=0;
      if((Sp->array=malloc(size*sizeof(SpRow)))==NULL) return NULL;
      if(Max_size_row->size !=size)
        PNL_ERROR(" Error in size ", "pnl_morse_mat_create");
      while(i<size)
        {
          pnl_sprow_init(&(Sp->array[i]),
                         pnl_vect_int_get(Max_size_row,i));
          i++;
        }
    }
  else
    Sp->array=(SpRow*) NULL;
  return Sp;
}


/**
 * create a PnlMat from PnlMorseMat
 * @param M a pointer on PnlMorseMat 
 * @return a PnlMat pointer
 */
PnlMat * pnl_morse_mat_full(PnlMorseMat * M)
{
  int i,k;PnlMat * FM;
  SpRow * Row;
  double *ptr_value;
  int*ptr_index;  
  FM=pnl_mat_create_from_double(M->m,M->n,0.0);
  i=0;
  if(!M->RC)
    while(i<M->m)
    {
      k=0;
      Row=&(M->array[i]);
      ptr_value=&(Row->Value[0]);
      ptr_index=&(Row->Index[0]);
      while(k<Row->size)
        {
          MLET(FM,i,*ptr_index)=*ptr_value;
          ptr_value++;ptr_index++;k++;
        }
      Row++;i++;
    }
  else
    while(i<M->n)
    {
      k=0;
      Row=&(M->array[i]);
      ptr_value=&(Row->Value[0]);
      ptr_index=&(Row->Index[0]);
      while(k<Row->size)
        {
          MLET(FM,*ptr_index,i)=*ptr_value;
          ptr_value++;ptr_index++;k++;
        }
      Row++;i++;
    }
  return FM;
}

/**
 * create a PnlMorseMat from PnlMat
 * @param FM a pointer on PnlMat 
 * @param RC int store in row or col
 * @return a PnlMorseMat pointer
 */
PnlMorseMat * pnl_morse_mat_create_fromfull(PnlMat * FM,int RC)
{
  int i,j;
  PnlMorseMat * M;
  SpRow * Row;
  double val;
  M=pnl_morse_mat_create(FM->m,FM->n,10.0,RC);
  if(!M->RC)
    for(i=0;i<M->m;i++)
    {
      Row=&(M->array[i]);
      for(j=0;j<FM->n;j++)
      {
        if((val=MGET(FM,i,j))!=0)
          pnl_sprow_add(Row,j,val);
      }
      Row++;
    }
  else
    {
      for(j=0;j<M->n;j++)
        {
          Row=&(M->array[j]);
          for(i=0;i<FM->m;i++)
            {
              if((val=MGET(FM,i,j))!=0)
                pnl_sprow_add(Row,i,val);
            }
          Row++;
        }
    }
  pnl_morse_mat_freeze(M);
  return M;
}


/**
 * free a PnlMorseMat pointer and set the data pointer to
 * NULL
 * @param M  address of the pointer to free
 */
void pnl_morse_mat_free(PnlMorseMat ** M)
{
  if((*M)!=NULL)
    {
      int i=0;;
      while(i<(GETROWORCOLSIZE(*M)))
        {
          pnl_sprow_free(&((*M)->array[i]));
          i++;
        }
      free((*M)->array);
      free(*M);
      *M=NULL;
     }
}



/**
 * gets the value of M[i,j]
 *
 * @param M : a PnlMorseMat
 * @param i : index of line
 * @param j : index of col 
 * @return  M[i,j]
 */
double pnl_morse_mat_get(PnlMorseMat* M, int i, int j)
{
  int k;
  SpRow * Row; 
  CheckIndexSparseMat(M,i,j);
  Row = GETROWORCOL(M,i,j);
  k=pnl_sprow_search(Row,GETROWORCOLIND(M,i,j));
  if(k<Row->size)
    return M->array->Value[k];
  PNL_ERROR(" Coefficient not in PnlMorseMat", "pnl_morse_mat_get")
  return 0.0;
}

/**
 * sets the value of M[i,j]=x
 *
 * @param M   : a PnlMorseMat
 * @param i   : index of line
 * @param j   : index of col 
 * @param Val : M[i,j]=x
 * @return  Error Code
 */
int pnl_morse_mat_set(PnlMorseMat* M, int i, int j,double Val)
{
  CheckIndexSparseMat(M,i,j);
  return  pnl_sprow_set(GETROWORCOL(M,i,j),GETROWORCOLIND(M,i,j),Val);
}

/**
 * returns the address of M[i,j] for use as a lvalue.
 *
 * @param M : a PnlMorseMat
 * @param i : index of line
 * @param j : index of col 
 * @return  &(M[i,j])
 */
double* pnl_morse_mat_lget(PnlMorseMat* M, int i, int j)
{
  CheckIndexSparseMat(M,i,j);
  return pnl_sprow_lget(GETROWORCOL(M,i,j),GETROWORCOLIND(M,i,j));
}

/**
 * add entrie M[i,j]=x, must sure i,j doest not exist,
 *
 * @param M   a PnlMorseMat
 * @param i   index of line
 * @param j   index of col 
 * @param Val a double M[i,j]=x
 * @return  Error Code
 */
int pnl_morse_mat_add(PnlMorseMat* M, int i, int j,double Val)
{
  CheckIndexSparseMat(M,i,j);
  return  pnl_sprow_add(GETROWORCOL(M,i,j),GETROWORCOLIND(M,i,j),Val);
}

/**
 * Put Max_size to size and free memory not needed.
 *
 * @param M : a PnlMorseMat
 * @return  Error code
 */
int pnl_morse_mat_freeze(PnlMorseMat* M)
{
  int i;
  if(!M->RC)
    for(i=0;i<M->m;i++)
    if(pnl_sprow_freeze(&(M->array[i]))==FAIL)
      return FAIL;
  else
    for(i=0;i<M->n;i++)
    if(pnl_sprow_freeze(&(M->array[i]))==FAIL)
      return FAIL;
  return OK;
}

/**
 * Prints a matrix to a file
 *
 * @param fic a file descriptor.
 * @param M a PnlMorseMat pointer.
 */
void pnl_morse_mat_fprint (FILE *fic, const PnlMorseMat *M)
{
  int i,k;
  double *ptr_value;
  int *ptr_index;
  SpRow * Row; 
  i=0;
  if(!M->RC)
    while(i<M->m)
      {
        k=0;
        Row=&(M->array[i]);
        ptr_value=&(Row->Value[0]);
        ptr_index=&(Row->Index[0]);
        while(k<Row->size)
          {
            fprintf (fic, "(%d,%d) %7.4f \n" , i,*ptr_index,*ptr_value);
            ptr_value++;ptr_index++;k++;
          }
        if(Row->size>0)
          fprintf (fic, "\n");
        Row++;i++;
      }
  else
    while(i<M->n)
      {
        k=0;
        Row=&(M->array[i]);
        ptr_value=&(Row->Value[0]);
        ptr_index=&(Row->Index[0]);
        while(k<Row->size)
          {
            fprintf (fic, "(%d,%d) %7.4f \n" , *ptr_index,i,*ptr_value);
            ptr_value++;ptr_index++;k++;
          }
        if(Row->size>0)
          fprintf (fic, "\n");
        Row++;i++;
      }

}

/**
 * Prints a matrix to the standard output
 * @param M a PnlMorseMat pointer.
 */
void pnl_morse_mat_print (const PnlMorseMat *M) { pnl_morse_mat_fprint(stdout, M);}


/**
 *  matrix multiplication
 *
 * @param M : matrix
 * @param vec : vector
 * @return  mat*vec
 */
PnlVect* pnl_morse_mat_mult_vect(const PnlMorseMat *M, const PnlVect *vec)
{
  PnlVect *lhs;
  CheckMatVectIsCompatible(M, vec);
  if ((lhs=pnl_vect_create(M->m))==NULL)
    return NULL;
  pnl_morse_mat_mult_vect_inplace(lhs,M,vec);
  return lhs;
}

/**
 *  in place matrix multiplication
 *
 * @param M : matrix
 * @param lhs : vector
 * @param rhs : vector
 * @return  lhs=M*rhs
 */
void pnl_morse_mat_mult_vect_inplace(PnlVect *lhs, const PnlMorseMat *M, const PnlVect *rhs)
{
  double *lptr,*ptr_value;
  int *ptr_index;
  SpRow * Row; 
  int i,k;
  pnl_vect_resize(lhs,M->m);
  CheckMatVectIsCompatible(M, rhs);
  lptr=lhs->array;
  Row =&(M->array[0]);
  if(!M->RC)
    {
      i=0;
      while(i<lhs->size)
        {
          k=0;
          *lptr=0.0;
          ptr_value=&(Row->Value[0]);
          ptr_index=&(Row->Index[0]);
          while (k<Row->size)
            {
              *lptr += *ptr_value * GET(rhs,*ptr_index);
              ptr_value++; ptr_index++; k++;
            }
          lptr++;Row++;i++;
        }
    }
  else
    {
      pnl_vect_set_zero(lhs);
      i=0;
      while(i<rhs->size)
        {
          k=0;
          ptr_value=&(Row->Value[0]);
          ptr_index=&(Row->Index[0]);
          while (k<Row->size)
            {
              LET(lhs,*ptr_index) += *ptr_value * GET(rhs,i);
              ptr_value++; ptr_index++; k++;
            }
          lptr++;Row++;i++;
        }
    }
}


/**
 *  constructor of PnlSparseMat from a PnlMat
 *
 * @param M : matrix
 * @return  pointer on PnlSparseMat
 */
PnlSparseMat *pnl_sparse_mat_create_fromfull(PnlMat * M)
{
    int i, j ;
    cs *T,*A ;
    if (!M) return (NULL) ;				/* check inputs */
    T = cs_spalloc (0, 0, 1, 1, 1) ;			/* allocate result */
    for(i=0;i<M->m;i++)
      for(j=0;j<M->n;j++)
        if (MGET(M,i,j)!=0)
          if (!cs_entry (T, i, j, MGET(M,i,j))) return (cs_spfree (T)) ;
    A=cs_compress (T);
    cs_spfree(T);
    return (A) ;
}

/**
 *  constructor of PnlSparseMat from a PnlMorseMat (RC==1)!
 *
 * @param M : Morse matrix
 * @return  pointer on PnlSparseMat
 */
PnlSparseMat *pnl_sparse_mat_create_frommorse_old(PnlMorseMat * M)
{
    int i, j ;
    cs *A ,*T;
    if (!M) return (NULL) ;				/* check inputs */
    T = cs_spalloc (0, 0, 1, 1, 1) ;			/* allocate result */
    for(i=0;i<M->m;i++)
      for(j=0;j<M->array[i].size;j++)
        if (!cs_entry (T, i, M->array[i].Index[j],  M->array[i].Value[j])) return (cs_spfree (T)) ;
     A=cs_compress (T);
    cs_spfree(T);
    return (A) ;
}

/**
 *  constructor of PnlSparseMat from a PnlMorseMat (RC==1)!
 *
 * @param M : Morse matrix
 * @return  pointer on PnlSparseMat
 */
PnlSparseMat *pnl_sparse_mat_create_frommorse(PnlMorseMat * M)
{
  int i,j,k,ok;
  cs *T;
  if (!M) return (NULL) ;
  if(!M->RC)
    {
      PNL_ERROR(" Morse Mat convertion is not implemented for rox stored morse matrix ", "pnl_sparse_mat_create_frommorse");
      return NULL;
    }
  /* check inputs */
  T = cs_spalloc (0, 0, 1, 1, 1) ;
  T->m = M->m;
  T->n =M->n;
  T->nz=-1;
  /* compressed col*/
  T->nzmax=0;
  for(j=0;j<M->n;j++)
    T->nzmax+=M->array[j].size;
  if (!cs_sprealloc (T,T->nzmax)) return (NULL);
  T->p = cs_realloc (T->p, (T->n)+1, sizeof (int), &ok) ;
  k=0;
  for(j=0;j<M->n;j++)
    {
      T->p[j]=k;
      for(i=0;i<M->array[j].size;i++)
        {
          T->x[k]=M->array[j].Value[i];
          T->i[k]=M->array[j].Index[i];
          k++;
        }
    }
  T->p[M->n]=k;
  return (T) ;
}

/**
 *  destructor for PnlSparseMat
 *
 * @param A adress of pointer on PnlSparseMat
 */
void pnl_sparse_mat_free(PnlSparseMat **A)
{
  *A=cs_spfree(*A);
}

/**
 *  Print function for PnlSparseMat
 *
 * @param A pointer on PnlSparseMat
 */
void pnl_sparse_mat_print(PnlSparseMat *A)
{
  int p, j, m, n, nzmax, nz, *Ap, *Ai ;
  double *Ax ;
    if (!A) { printf ("(null)\n") ; }
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    nzmax = A->nzmax ; nz = A->nz ;
    if (nz < 0)
    {
      printf ("%d-by-%d, nzmax: %d nnz: %d, 1-norm: %g\n", m, n, nzmax,
              Ap [n], cs_norm (A)) ;
      for (j = 0 ; j < n ; j++)
        {
          printf ("    col %d : locations %d to %d\n", j, Ap [j], Ap [j+1]-1);
          for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
              printf ("      %d : %g\n", Ai [p], Ax ? Ax [p] : 1) ;
              if (p > 20) { printf ("  ...\n") ;}
            }
        }
    }
    else
      {
        printf ("triplet: %d-by-%d, nzmax: %d nnz: %d\n", m, n, nzmax, nz) ;
        for (p = 0 ; p < nz ; p++)
          {
            printf ("    %d %d : %g\n", Ai [p], Ap [p], Ax ? Ax [p] : 1) ;
            if (p > 20) { printf ("  ...\n") ;}
          }
      }
}

/**
 *  compute lhs=lhs+ M * rhs
 *
 * @param lhs a PnlVect
 * @param M aPnlSparseMat
 * @param rhs a PnlVect
 * @return Error code
 */
int pnl_sparse_mat_gaxpby(PnlVect *lhs, const PnlSparseMat *M, const PnlVect *rhs)
{
  return cs_gaxpy (M,rhs->array,lhs->array);
}

/**
 *  compute lhs=M * rhs
 *
 * @param lhs a PnlVect
 * @param M aPnlSparseMat
 * @param rhs a PnlVect
 * @return Error code
 */
int pnl_sparse_mat_mult_vect_inplace(PnlVect *lhs, const PnlSparseMat *M, const PnlVect *rhs)
{
  pnl_vect_set_zero(lhs);
  return cs_gaxpy (M,rhs->array,lhs->array);
}



/**
 * in-place map function
 *
 * @param lhs left hand side PnlSparseMat
 * @param f the function to be applied term by term
 * @return  lhs = f(lhs)
 */
void pnl_sparse_matrix_map_inplace(PnlSparseMat *lhs, 
                                double(*f)(double ))
{
  pnl_array_map_inplace(lhs->x,f,lhs->nzmax);
}

/**
 * in-place PnlSparseMat scalar addition
 *
 * @param lhs left hand side PnlSparseMat
 * @param x scalar
 * @return  lhs = lhs+x
 */
void pnl_sparse_matrix_plus_double(PnlSparseMat *lhs , double x)
{
  pnl_array_plus_double(lhs->x,x,lhs->nzmax);     
}

/**
 * in-place PnlSparseMat scalar substraction
 *
 * @param lhs left hand side PnlSparseMat
 * @param x scalar
 * @return  lhs = lhs-x
 */
void pnl_sparse_matrix_minus_double(PnlSparseMat *lhs , double x)
{
  pnl_array_minus_double(lhs->x,x,lhs->nzmax);
}

/**
 * in-place PnlSparseMat scalar multiplication
 *
 * @param lhs left hand side PnlSparseMat
 * @param x scalar
 * @return  lhs = lhs*x
 */
void pnl_sparse_matrix_mult_double(PnlSparseMat *lhs , double x)
{
  pnl_array_mult_double(lhs->x,x,lhs->nzmax);
}

/**
 * in-place PnlSparseMat scalar division
 *
 * @param lhs left hand side PnlSparseMat
 * @param x scalar
 * @return  lhs = lhs/x
 */
void pnl_sparse_matrix_div_double(PnlSparseMat *lhs , double x)
{
  pnl_array_div_double(lhs->x,x,lhs->nzmax);
}

/**
 * in-place PnlSparseMat PnlSparseMat addition
 *
 * @param lhs left hand side PnlSparseMat
 * @param rhs rigth hand side PnlSparseMat
 * @return  lhs = lhs+rhs
 */
void pnl_sparse_matrix_plus_mat(PnlSparseMat *lhs, const PnlSparseMat *rhs)
{
  pnl_array_plus_array_term(lhs->x, rhs->x,lhs->nzmax);
}

void pnl_sparse_matrix_minus_mat(PnlSparseMat *lhs, const PnlSparseMat *rhs)
{
  CheckSparseMatMatch(lhs,rhs);
  pnl_array_minus_array_term(lhs->x,rhs->x,lhs->nzmax);
}

/**
 * in-place term by term PnlSparseMat inverse
 *
 * @param lhs left hand side PnlSparseMat
 * @return  lhs = 1 ./ lhs
 */
void pnl_sparse_matrix_inv_term(PnlSparseMat *lhs)
{
  pnl_array_inv_term(lhs->x,lhs->nzmax);
}

/**
 * in-place term by term PnlSparseMat inverse
 *
 * @param lhs left hand side PnlSparseMat
 * @param rhs right hand side PnlSparseMat
 * @return  lhs = lhs ./ rhs
 */
void pnl_sparse_matrix_div_mat_term(PnlSparseMat *lhs, const PnlSparseMat *rhs)
{
  CheckSparseMatMatch(lhs,rhs);
  pnl_array_div_array_term(lhs->x,rhs->x,lhs->nzmax);
}

/**
 * in-place PnlSparseMat term by term multiplication
 *
 * @param lhs left hand side PnlSparseMat
 * @param rhs right hand side PnlSparseMat
 * @return  lhs = lhs.*rhs
 */
void pnl_sparse_matrix_mult_mat_term(PnlSparseMat *lhs, const PnlSparseMat *rhs)
{
  CheckSparseMatMatch(lhs,rhs);
  pnl_array_mult_array_term(lhs->x,rhs->x,lhs->nzmax);
}

/**
 * compute LU factorization of Sparse Matrix A
 *
 * @param A aPnlSparseMat
 * @param tol double (if need pivoting)
 * @return PnlSparseFactorization pointer with L, U and permutation
 */
/*
  order is the order for fill-reducing.
  => order=0, then no column permutation is used; LU = PA. This is useful if A is
  already known to have a good column ordering - which is the case for a
  M-matrix.
  => order=1, if use if pattern of U and L are identical to pattrens of the
  Choleski factor L and L^T, respectively of a symmetric positive definite
  matrix with the same non zero patterns as A+A^T
  => order=2 or order=3 is to find optimal fill-reducing
  order=2 to is less optimal, less computationnal time thant order=3

  Here we create in interface only for order 0, because Col permuation q don't
  store in PnlSparseFactorization
 */
PnlSparseFactorization * pnl_sparse_factorization_lu_create (const PnlSparseMat*A, double tol)
{
  PnlSparseFactorization *F;
  int order=0;
  if((F=malloc(sizeof(PnlSparseFactorization)))==NULL)
    return NULL;
  F->S = cs_sqr (order, A, 0) ;
  /* ordering and symbolic analysis */
  F->N = cs_lu (A, F->S, tol) ;
  /* numeric LU factorization */
  if (!(F->S && F->N))
    PNL_ERROR(" LU Factorisation in PnlSparse", "pnl_sparse_factorization_lu_create");
  return (F) ;
};


/**
 * free PnlSparseFactorization
 *
 * @param F adress of pointer on PnlSparseFactorization
 */
void pnl_sparse_factorization_free(PnlSparseFactorization ** F)
{
  cs_sfree ((*F)->S) ;
  cs_nfree((*F)->N);
  free(*F);
  (*F)=NULL;
};

/**
 * Solve Lin syslin for LU sparse factorization 
 * Solve LU x = P x
 *
 * @param F a PnlSparseFactorization
 * @param x a PnlVect rhs begin lhs
 */
void pnl_sparse_factorization_lu_syslin(const PnlSparseFactorization * F, PnlVect *x)
{
  PnlVect *b =pnl_vect_copy(x);
  cs_ipvec (F->N->pinv, b->array,x->array, F->N->L->n) ;
  /* x = b(p) */
  cs_lsolve (F->N->L, x->array) ;
  /* x = L\x */
  cs_usolve (F->N->U, x->array) ;
  /* x = U\x */
  cs_ipvec (F->S->q, x->array, b->array, F->N->L->n) ;
  /* b(q) = x */
  pnl_vect_free (&b) ;
  
}


