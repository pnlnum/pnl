#include <stdlib.h>
#include <stdio.h>
#include "pnl_vector_uint.h"
#include "pde_tools.h"


/**
 * creates a PnlPDEBoundary
 * @param X0 down_left point of domain
 * @param X1 Up_right point of domain
 * @return a PnlPDEBoundary pointer
 */
PnlPDEBoundary pnl_pde_boundary_create(double X0,double X1)
{
  PnlPDEBoundary  BP;
  BP.X0=X0;
  BP.H=1/(X1-X0);
  return BP;
}

double pnl_pde_boundary_real_variable(const PnlPDEBoundary BP ,double X)
{return  BP.X0+X/BP.H;}
  
double pnl_pde_boundary_unit_interval(const PnlPDEBoundary BP ,double X)
{
  double res= (X-BP.X0)*BP.H;
  if (abs(2*res-1)<=1.0) {printf("error on boundary Unit_Interval");abort();}
  return res;
}

/**
 * creates a PnlPDEBoundary with
 * Left down corner is \f$ (0,\dots,0)\f$
 * and right up corner is \f$ (1,\dots,1)\f$
 * @return a PnlPDEDimBoundary pointer
 */
PnlPDEDimBoundary* pnl_pde_dim_boundary_create_from_int(int dim)
{
  PnlPDEDimBoundary *TBP;
  PnlPDEBoundary *Tmp;

  int i;
  if((TBP=malloc(sizeof(PnlPDEDimBoundary)))==NULL)
    return NULL;
  if (dim>0)
    {
      if((TBP->array=malloc(dim*sizeof(PnlPDEBoundary)))==NULL)
        return NULL;
    }
  else
    TBP->array = (PnlPDEBoundary*)NULL;
  Tmp=TBP->array;
  i=0;
  while(i<dim)
    {
      (*Tmp)=pnl_pde_boundary_create(0.0,1.0);
      Tmp++;i++;
    }
  return TBP;  
}

/**
 * creates a PnlMat
 * @param X0 left down corner
 * @param X1 right up corner
 * @return a PnlPDEDimBoundary pointer
 */
PnlPDEDimBoundary* pnl_pde_dim_boundary_create(const PnlVect *X0,
						     const PnlVect *X1)
{
  PnlPDEDimBoundary *TBP;
  PnlPDEBoundary *Tmp;
  int i;
  if((TBP=malloc(sizeof(PnlPDEDimBoundary)))==NULL)
    return NULL;
  if (X0->size>0)
    {
      if((TBP->array=malloc(X0->size*sizeof(PnlPDEBoundary)))==NULL)
        return NULL;
    }
  else
    TBP->array = (PnlPDEBoundary*)NULL;
  Tmp=TBP->array;
  i=0;
  while(i<X0->size)
    {
      (*Tmp)=pnl_pde_boundary_create(GET(X0,i),GET(X1,i));
      Tmp++;i++;
    }
  return TBP;  
}

/** 
 * frees a PnlPDEDimBoundary
 *
 * @param v adress of a PnlPDEDimBoundary*. v is set to NULL at exit.
 */
void pnl_pde_dim_boundary_free(PnlPDEDimBoundary **v)
{
  
  if (*v != NULL)
    {
      if ((*v)->array != NULL )
	{
	  free((*v)->array);
	  free(*v);
	  *v=NULL;
	}
    }
}

double pnl_pde_dim_boundary_eval_from_unit(double(*f)(const PnlVect* ),
						     const PnlPDEDimBoundary * BP,
						     const PnlVect * X)  
{	
  
  PnlVect * Res;
  double sol;
  int i;
  Res=pnl_vect_create(X->size);
  for(i=0;i<X->size;i++)
    LET(Res,i)=pnl_pde_boundary_real_variable(BP->array[i],GET(X,i));
  sol=f(Res);
  pnl_vect_free(&Res);
  return  sol;
}

void pnl_pde_dim_boundary_from_unit_to_real_variable(const PnlPDEDimBoundary * BP ,
							PnlVect * X)
{	
  int i;
  for(i=0;i<X->size;i++)
     LET(X,i)=pnl_pde_boundary_real_variable(BP->array[i],GET(X,i));
}

double pnl_pde_dim_boundary_get_step(const PnlPDEDimBoundary * BP,
					       int i) 
{ return BP->array[i].H;}





double standard_time_repartition(int i,int N_T)
{
  return (double) (i)/ (double) N_T;
}
/**
 * creates a PnlPDETimeGrid
 * @param T Terminal time
 * @param N_T number of grids points
 * @param repartition function for repartitions of grids points
 * @return a PnlPDETimeGrid pointer
 */
PnlPDETimeGrid * pnl_pde_time_grid(const double T,
					 const int N_T,
					 double (*repartition)(int i,int NN)
					 )
{
  PnlPDETimeGrid  *TG;
  int i;
  if((TG=malloc(sizeof(PnlPDETimeGrid)))==NULL)
    return NULL;
  if (N_T>0)
    {
      TG->time=pnl_vect_create(N_T+1);
      i=0;
      do{
	LET(TG->time,i)=T*repartition(i,N_T);
	i++;
      }while(i<=N_T);
      TG->is_tuned=1;
      pnl_pde_time_start(TG);
    }
  else
    TG->time= (PnlVect*)NULL;
  return TG;  
}
/**
 * creates a PnlPDETimeGrid
 * @param T Terminal time
 * @param N_T number of grids points
 * @return a PnlPDETimeGrid pointer
 */
PnlPDETimeGrid * pnl_pde_time_homogen_grid(const double T,
						 const int N_T)
{
  PnlPDETimeGrid  *TG=pnl_pde_time_grid(T,N_T,standard_time_repartition);
  TG->is_tuned=0;
  return TG;
}

  

/** 
 * frees a PnlPDETimeGrid
 *
 * @param TG adress of a PnlPDETimeGrid*. TG is set to NULL at exit.
 */
void pnl_pde_time_grid_free(PnlPDETimeGrid **TG)
{
  
  if (*TG != NULL)
    {
      pnl_vect_free(&((*TG)->time));
      free(*TG);
      *TG=NULL;
    }
}

/**
 * initialise PnlPDETimeGrid to the first step.
 *
 * @param TG a PnlPDETimeGrid pointer
 * initialise in first step
 */
void pnl_pde_time_start(PnlPDETimeGrid *TG)
{
  TG->current_index=0;
  TG->current_step=GET(TG->time,1)-GET(TG->time,0);
  /*-GET(TG->time,TG->current_index)+GET(TG->time,++(TG->current_index)); */
}

/**
 * go to the next time step 
 *
 * @param TG a PnlPDETimeGrid pointer
 * increase the current time and compute current step
 * @return 0 if last time step
 */
int pnl_pde_time_grid_increase(PnlPDETimeGrid * TG)
{
      if(TG->is_tuned==1)
	TG->current_step=-GET(TG->time,TG->current_index) + GET(TG->time,TG->current_index+1);
      (TG->current_index)+=1;
      return (TG->current_index<TG->time->size);
}


/**
 * GET function on current step
 *
 * @param TG a PnlPDETimeGrid pointer
 * @return the current step
 */
double pnl_pde_time_grid_step(const PnlPDETimeGrid * TG)
{ return TG->current_step;}

/**
 * GET function on current time
 *
 * @param TG a PnlPDETimeGrid pointer
 * @return the current time
 */
double pnl_pde_time_grid_time(const PnlPDETimeGrid * TG)
{return GET(TG->time,TG->current_index);}


