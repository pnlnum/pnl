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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pnl/pnl_matrix.h"
#include "pnl/pnl_linalgsolver.h"

#define DEBUG_SOLVER 0
#define INFO_SOLVER 0

static char pnl_iteration_label[] = "PnlIterationBase";
static char pnl_cg_solver_label[] = "PnlCgSolver";
static char pnl_bicg_solver_label[] = "PnlBicgSolver";
static char pnl_gmres_solver_label[] = "PnlGmresSolver";

/**
 * Create an empty PnlIterationBase.
 *
 * @return  a PnlIterationBasepointer
 */
static PnlIterationBase* pnl_iteration_base_new() 
{
  PnlIterationBase * o;
  if( (o=malloc(sizeof(PnlIterationBase)))==NULL ) return NULL;
  o->iteration=0;
  o->max_iter=0;
  o->normb=0;
  o->tol_=0;
  o->resid=0.0;
  o->error=0;
  o->object.type = PNL_TYPE_ITERATION_BASE;
  o->object.parent_type = PNL_TYPE_OBJECT;
  o->object.label = pnl_iteration_label;
  return o;
}

/**
 * Create a new PnlIterationBase pointer.
 *
 * @param max_iter_ : maximum number of iterations 
 * @param t : tolerance, parameter of iterative methods, ie relative error
 * @return  a PnlIterationBase pointer
 */
static PnlIterationBase* pnl_iteration_base_create(int max_iter_, double t) 
{
  PnlIterationBase * It;
  if((It = pnl_iteration_base_new ())==NULL) return NULL;
  It->max_iter=max_iter_;
  It->normb=1.0;
  It->tol_=t;
  return It;
}

/**
 * Test if the iterative solution has converged
 *
 * @param it : a PnlIterationBase ptr. 
 * @param rhs : a PnlVect ptr
 * @return   TRUE or FALSE.
 */
static int pnl_iteration_base_converged(PnlIterationBase * it,const PnlVect *rhs) 
{
  it->resid = pnl_vect_norm_two(rhs)/it->normb;
#if DEBUG_SOLVER
  printf(" residu at iteration %d is %7.4f \n",it->iteration,it->resid);
#endif  
#if INFO_SOLVER
  if (it->resid <= it->tol_)
    {
      printf(" converged in %d iterations \n",it->iteration);
      return true;
    }
  return false;
#endif  
  return (it->resid <= it->tol_)?true:false;
}

/**
 * Test if the iterative solution has converged
 *
 * @param it : a PnlIterationBase ptr. 
 * @param Res : a double
 * @return   TRUE or FALSE.
 */
static int pnl_iteration_base_converged_norm(PnlIterationBase * it,double Res) 
{
  it->resid = Res/it->normb;
  return (it->resid <= it->tol_)?true:false;
}

/**
 * Test if the iterative method has finished 
 * ie convergence or iteration index is larger than max_iteration?
 *
 * @param it : a PnlIterationBase ptr. 
 * @param rhs : a PnlVect ptr
 * @return   TRUE or FALSE.
 */
static int pnl_iteration_base_finished(PnlIterationBase * it, const  PnlVect *rhs) 
{
  if (pnl_iteration_base_converged(it,rhs))
    return true;
  /*    printf("converged in %d iterations , with residu = %f max_iter = %d \n" , */
  /*       it->iteration,it->resid,it->max_iter); */
  if ((*it).iteration >(*it).max_iter)
    {
      (*it).error = 1;
      printf("iterative Solver not converged in last iterations with residu = %f iter = %d \n" ,
             it->resid,it->iteration);
      abort();
      return true;
    }
  return false;
       
}

/* /\**
 *  * Test if iterative method is finished 
 *  * ie converged or iteration index is larger than max_iteration 
 *  *
 *  * @param it : a PnlIterationBase ptr. 
 *  * @param Res : a double
 *  * @return   TRUE or FALSE.
 *  *\/
 * static int pnl_iteration_base_finished_norm(PnlIterationBase * it,double Res)
 * {
 *   if (pnl_iteration_base_converged_norm(it,Res))
 *     return true;
 *   /\*       printf("converged in %d iterations , with residu = %f max_iter = %d \n" , *\/
 *   /\*       it->iteration,it->resid,it->max_iter); *\/
 *   if ((*it).iteration > (*it).max_iter)
 *     {
 *       (*it).error = 1;
 *       printf("iterative Solver not converged in last iterations with residu = %f iter = %d \n" ,
 *              it->resid,it->iteration);
 *       abort();
 *       return true;
 *     }
 *   return false;
 * } */

/**
 * Initiatise PnlIterationBase
 * store norm of RHS ohf the linear system
 *
 * @param it : a PnlIterationBase ptr. 
 * @param b  : a PnlVect ptr
 */
static void pnl_iteration_base_initialisation(PnlIterationBase * it,const PnlVect *b ) 
{
  (*it).iteration = 0;
  it->normb=pnl_vect_norm_two(b);
}

/**
 * Initiatise PnlIterationBase
 * store norm of RHS ohf the linear system
 *
 * @param it : a PnlIterationBase ptr. 
 * @param normb : a double
 */
static void pnl_iteration_base_initialisation_norm(PnlIterationBase * it,double normb ) 
{
  (*it).iteration = 0;
  it->normb=fabs(normb);
}

/**
 * increment PnlIterationBase
 * iteration++
 *
 * @param it : a PnlIterationBase ptr. 
 */
static void pnl_iteration_base_increment(PnlIterationBase * it)
{ (*it).iteration++; }

/**
 * Give error code 
 *
 * @param it : a PnlIterationBase ptr. 
 * @return  int
 */
static int pnl_iteration_base_error_code(PnlIterationBase * it) { return it->error; }

/*
  static int pnl_iteration_base_iterate(PnlIterationBase * it) 
  {return (it->iteration < it->max_iter)?true:false;}
  static int pnl_iteration_base_iterations(PnlIterationBase * it) { return it->iteration; }
  static double pnl_iteration_base_residu(PnlIterationBase * it) { return it->resid; }
  static double pnl_iteration_base_norm_rhs(PnlIterationBase * it) { return it->normb; }
  static double pnl_iteration_base_tolerance(PnlIterationBase * it) { return (*it).tol_; }
  static int pnl_iteration_base_iter_max(PnlIterationBase * it){return (*it).max_iter;}    

*/

/**
 * Test first iteration 
 *
 * @param it : a PnlIterationBase ptr. 
 * @return   TRUE or FALSE
 */
static int pnl_iteration_base_first(PnlIterationBase * it) { return it->iteration == 0; }

/**
 * Restart iterationbase, ie iteration=0
 *
 * @param it : a PnlIterationBase ptr. 
 */
static void  pnl_iteration_base_begin(PnlIterationBase * it) { it->iteration = 0; }


/**
 * Store error code if Iteration method fail, ie in BREAKDOWN case
 *
 * @param it : a PnlIterationBase ptr. 
 * @param err_code : a int. 
 */ 
static void pnl_iteration_base_fail(PnlIterationBase * it,int err_code) { it->error = err_code; }

/**
 * Store error code & error message if Iteration method fail, ie in BREAKDOWN case
 * 
 * @param it : a PnlIterationBase ptr. 
 * @param err_code : a int. 
 * @param msg : a char ptr. 
 */ 
static void pnl_iteration_base_failed(PnlIterationBase * it,int err_code, 
                                     const char * msg)
{ (*it).error = err_code;  PNL_ERROR(msg,"iterative solver" );
}


/**
 * Create an empty PnlCgSolver 
 *
 *@return  a PnlCgSolver pointer
 */
PnlCgSolver* pnl_cg_solver_new()
{
  PnlCgSolver * o;
  if((o=malloc(sizeof(PnlCgSolver)))==NULL) return NULL;
  o->r = o->z = o->p = o->q = NULL;
  o->iter=pnl_iteration_base_new();
  o->object.type = PNL_TYPE_CG_SOLVER;
  o->object.parent_type = PNL_TYPE_OBJECT;
  o->object.label = pnl_cg_solver_label;
  o->object.clone = NULL;
  o->object.copy = NULL;
  o->object.constructor = (NewFunc *) pnl_cg_solver_new;
  return o;
}

/**
 * Create a new PnlCgSolver pointer.
 *
 * @param Size : the size of temporary vectors use in the method
 * @param max_iter_  : maximum iteration number 
 * @param tolerance_ : relative error of iterative method
 *@return  a PnlCgSolver pointer
 */
PnlCgSolver* pnl_cg_solver_create(int Size,int max_iter_, double tolerance_)
{
  PnlCgSolver * Solver;
  if((Solver = pnl_cg_solver_new ()) ==NULL) return NULL;
  Solver->r=pnl_vect_create_from_scalar(Size,0.0);
  Solver->z=pnl_vect_create_from_scalar(Size,0.0);
  Solver->p=pnl_vect_create_from_scalar(Size,0.0);
  Solver->q=pnl_vect_create_from_scalar(Size,0.0);
  Solver->iter=pnl_iteration_base_create(max_iter_,tolerance_);
  return Solver;
}

/**
 * initialisation of the solver at the beginning of iterative method.
 *
 * @param Solver : PnlCgSolver  ptr
 * @param b       : PnlVect ptr 
 */
void pnl_cg_solver_initialisation(PnlCgSolver * Solver,const PnlVect * b)
{
  pnl_vect_clone(Solver->r,b);
  Solver->rho=0.0;
  Solver->oldrho=0.0;
  Solver->beta=0.0;
  Solver->alpha=0.0;
  pnl_iteration_base_initialisation(Solver->iter,b);
}

/**
 * destructor of iterative solver 
 *
 * @param Solver : PnlCgSolver  ptr
 */
void pnl_cg_solver_free(PnlCgSolver ** Solver)
{
  if (*Solver != NULL)
    {
      pnl_vect_free(&(*Solver)->r);
      pnl_vect_free(&(*Solver)->z);
      pnl_vect_free(&(*Solver)->p);
      pnl_vect_free(&(*Solver)->q);
      free((*Solver)->iter);
      free(*Solver);
      *Solver=NULL;
    } 
}

/**
 * Solve the linear system
 * matrix_vector_product is the matrix_vector multiplication function
 * matrix_vector_product_PC is the preconditionner function
 * Matrix_Data & PC_Data is data to compute matrix_vector multiplication
 *
 * Solve :
 * Matrix_Data x =  b with preconditionner PC_Data
 *
 * @param matrix_vector_product : ptr function retrun PnlVect ptr and take 
 * void ptr and PnlVect ptr parameters 
 * @param Matrix_Data  : void ptr 
 * @param matrix_vector_product_PC : ptr function retrun PnlVect ptr and take 
 * void ptr and PnlVect ptr parameters 
 * @param PC_Data : void ptr
 * @param x       : PnlVect ptr 
 * @param b       : PnlVect ptr
 * @param Solver  : PnlCgSolver ptr 
 * @return Int : error code.
 */
int pnl_cg_solver_solve(void (* matrix_vector_product )(const void *,const PnlVect*, double, double, PnlVect*), 
                        const void * Matrix_Data,
                        void (*matrix_vector_product_PC)(const void *,const PnlVect*,
                                                         double, double, PnlVect*), 
                        const void * PC_Data,
                        PnlVect * x, 
                        const PnlVect *b,
                        PnlCgSolver * Solver)
{
  /* residual */
  pnl_cg_solver_initialisation(Solver,b);
  /* r = b - A x: initial residual */
  matrix_vector_product(Matrix_Data,x,-1.0,1.0,Solver->r);
  while(! pnl_iteration_base_finished(Solver->iter,Solver->r) )
    {
      matrix_vector_product_PC(PC_Data,Solver->r,1.0,0.0,Solver->z);
      Solver->rho = pnl_vect_scalar_prod(Solver->r,Solver->z);
      if (pnl_iteration_base_first(Solver->iter))
        pnl_vect_axpby(1.0, Solver->z, 0.0, Solver->p);
      else
        {
          Solver->beta = Solver->rho/Solver->oldrho;
          if (fabs(Solver->beta) <1.0e-9)
            {
              pnl_iteration_base_fail(Solver->iter,2);
              printf(">> Error code %d \n",pnl_iteration_base_error_code(Solver->iter));
              return pnl_iteration_base_error_code(Solver->iter);
            }
      
          pnl_vect_axpby(1.0, Solver->z, Solver->beta, Solver->p);
        }
      matrix_vector_product(Matrix_Data,Solver->p,1.0,0.0,Solver->q);
      Solver->alpha = Solver->rho/pnl_vect_scalar_prod(Solver->p,Solver->q);
      /*printf("cg4   %f \n",Solver->alpha); */
      /* update x and r */
      pnl_vect_axpby(Solver->alpha, Solver->p, 1.0, x);     
      pnl_vect_axpby(-Solver->alpha, Solver->q, 1.0, Solver->r);
      Solver->oldrho = Solver->rho;
      pnl_iteration_base_increment(Solver->iter);
    }
  return pnl_iteration_base_error_code(Solver->iter);
}


/**
 * Create an empty PnlBicgSolver 
 *
 *@return  a PnlBicgSolver pointer
 */
PnlBicgSolver* pnl_bicg_solver_new()
{
  PnlBicgSolver * o;
  if((o=malloc(sizeof(PnlBicgSolver)))==NULL) return NULL;
  o->r = NULL;
  o->rtilde = NULL;
  o->p = NULL;
  o->phat = NULL;
  o->s = NULL;
  o->shat = NULL;
  o->t = NULL;
  o->v = NULL;
  o->iter=pnl_iteration_base_new();
  o->object.type = PNL_TYPE_CG_SOLVER;
  o->object.parent_type = PNL_TYPE_OBJECT;
  o->object.label = pnl_bicg_solver_label;
  o->object.clone = NULL;
  o->object.copy = NULL;
  o->object.constructor = (NewFunc *) pnl_bicg_solver_new;
  return o;
}

/**
 * Create a new PnlBicgSolver pointer.
 *
 * @param Size : the size of temporary vectors use in the method
 * @param max_iter_  : maximum iteration number 
 * @param tolerance_ : relative error of iterative method
 *@return  a PnlBicgSolver pointer
 */
PnlBicgSolver* pnl_bicg_solver_create(int Size,int max_iter_, double tolerance_)
{
  PnlBicgSolver * Solver;
  if ((Solver = pnl_bicg_solver_new ()) == NULL) return NULL;
  Solver->r=pnl_vect_create_from_scalar(Size,0.0);
  Solver->rtilde=pnl_vect_create_from_scalar(Size,0.0);
  Solver->p=pnl_vect_create_from_scalar(Size,0.0);
  Solver->phat=pnl_vect_create_from_scalar(Size,0.0);
  Solver->s=pnl_vect_create_from_scalar(Size,0.0);
  Solver->shat=pnl_vect_create_from_scalar(Size,0.0);
  Solver->t=pnl_vect_create_from_scalar(Size,0.0);
  Solver->v=pnl_vect_create_from_scalar(Size,0.0);
  Solver->iter=pnl_iteration_base_create(max_iter_,tolerance_);
  return Solver;
}

  
/**
 * initialisation of the solver at the beginning of iterative method.
 *
 * @param Solver : PnlBicgSolver  ptr
 * @param b       : PnlVect ptr 
 */
void pnl_bicg_solver_initialisation(PnlBicgSolver * Solver,const PnlVect * b)
{
  pnl_vect_clone(Solver->r,b);
  Solver->rho_1=0.0;
  Solver->rho_2=0.0;
  Solver->alpha=0.0;
  Solver->beta=0.0;
  Solver->omega=0.0;
  pnl_iteration_base_initialisation(Solver->iter,b);
}

/**
 * destructor of iterative solver 
 *
 * @param Solver : PnlBicgSolver  ptr
 */
void pnl_bicg_solver_free(PnlBicgSolver ** Solver)
{
  if (*Solver != NULL)
    {
      pnl_vect_free(&(*Solver)->p);
      pnl_vect_free(&(*Solver)->phat);
      pnl_vect_free(&(*Solver)->s);
      pnl_vect_free(&(*Solver)->shat);
      pnl_vect_free(&(*Solver)->t);
      pnl_vect_free(&(*Solver)->v);
      pnl_vect_free(&(*Solver)->r);
      pnl_vect_free(&(*Solver)->rtilde);
      free((*Solver)->iter);
      free(*Solver);
      *Solver=NULL;
    } 
}

/**
 * Solve the linear system
 * matrix_vector_product is the matrix_vector multiplication function
 * matrix_vector_product_PC is the preconditionner function
 * Matrix_Data & PC_Data is data to compute matrix_vector multiplication
 *
 * Solve :
 * Matrix_Data x =  b with preconditionner PC_Data
 *
 * @param matrix_vector_product : ptr function retrun PnlVect ptr and take 
 * void ptr and PnlVect ptr parameters 
 * @param Matrix_Data  : void ptr 
 * @param matrix_vector_product_PC : ptr function retrun PnlVect ptr and take 
 * void ptr and PnlVect ptr parameters 
 * @param PC_Data : void ptr
 * @param x       : PnlVect ptr 
 * @param b       : PnlVect ptr
 * @param Solver  : PnlBicgSolver ptr 
 * @return Int : error code.
 */
int pnl_bicg_solver_solve(void (* matrix_vector_product )(const void *,const PnlVect*,double,double,PnlVect*), 
                          const void * Matrix_Data,
                          void (*matrix_vector_product_PC)(const void *,const PnlVect*,double,double,PnlVect*), 
                          const void * PC_Data,
                          PnlVect * x, 
                          const PnlVect *b,
                          PnlBicgSolver * Solver)
{
  pnl_bicg_solver_initialisation(Solver,b);
  /* r = b - A x: initial residual */
  matrix_vector_product(Matrix_Data,x,-1.0,1.0,Solver->r);
  pnl_vect_clone(Solver->rtilde,Solver->r);            
  /* r~ = r */
  while(! pnl_iteration_base_finished(Solver->iter,Solver->r) )
    {
      Solver->rho_1 = pnl_vect_scalar_prod(Solver->rtilde,Solver->r);
      if (Solver->rho_1==0)
        pnl_iteration_base_failed(Solver->iter,2, "bicg breakdown #1");
      if (pnl_iteration_base_first(Solver->iter))
        pnl_vect_axpby(1.0, Solver->r, 0.0, Solver->p);
      else 
        {
          if (Solver->omega == 0.) 
            {
              pnl_iteration_base_failed(Solver->iter,3, "bicg breakdown #2");
              break;
            }
      
          Solver->beta = (Solver->rho_1 / Solver->rho_2) * (Solver->alpha / Solver->omega);
          pnl_vect_axpby(-1.0*Solver->omega, Solver->v, 1.0, Solver->p);
          pnl_vect_axpby( 1.0, Solver->r, Solver->beta, Solver->p);
          /* p = r + beta * (p - omega * v) */
        }
      matrix_vector_product_PC(PC_Data,Solver->p,1.0,0.0,Solver->phat );
      /*phat=M p */
      /*cout<< "|phat |="<<dot(phat,phat)<<"\n"; */
      matrix_vector_product(Matrix_Data,Solver->phat,1.0,0.0,Solver->v);
      /* v = A * phat */
      Solver->alpha = Solver->rho_1/pnl_vect_scalar_prod(Solver->v, Solver->rtilde);
      pnl_vect_axpby(-1.0*Solver->alpha, Solver->v, 0.0, Solver->s); 
      pnl_vect_axpby(1.0, Solver->r, 1.0, Solver->s);
      /* s = r - alpha * v */
      if(pnl_iteration_base_converged(Solver->iter,Solver->s)) 
        {
          pnl_vect_axpby(Solver->alpha, Solver->phat, 1.0, x); 
          /* x += alpha * phat */
          break;
        }
      matrix_vector_product_PC(PC_Data,Solver->s,1.0,0.0, Solver->shat);
      /* shat = M.Solve(s) */
      matrix_vector_product(Matrix_Data,Solver->shat,1.0,0.0,Solver->t);   
      /* t = A * shat; */
      Solver->omega = pnl_vect_scalar_prod(Solver->t, Solver->s) / pnl_vect_scalar_prod(Solver->t, Solver->t);
      pnl_vect_axpby(Solver->omega, Solver->shat, 1.0, x);
      pnl_vect_axpby(Solver->alpha, Solver->phat, 1.0, x);
      /* x += alpha * phat + omega * shat */
      pnl_vect_axpby(1.0, Solver->s, 0.0, Solver->r);
      pnl_vect_axpby(-1.0*Solver->omega, Solver->t, 1.0, Solver->r);
      /* r = s - omega * t */
      Solver->rho_2 = Solver->rho_1;
      pnl_iteration_base_increment(Solver->iter);
    }
  return pnl_iteration_base_error_code(Solver->iter);
}


/**
 * Create an empty PnlGmresSolver
 *
 *@return  a PnlGmresSolver pointer
 */
PnlGmresSolver* pnl_gmres_solver_new()
{
  int i;
  PnlGmresSolver *o;
  if((o=malloc(sizeof(PnlGmresSolver)))==NULL) return NULL;
  o->restart = 0;
  o->s = NULL;
  o->cs = NULL;
  o->sn = NULL;
  o->w = NULL;
  o->r = NULL;
  o->H = NULL;
  for(i=0;i<MAX_RESTART;i++) { o->v[i] = NULL; }
  o->iter=pnl_iteration_base_new();
  o->iter_inner=pnl_iteration_base_new();
  o->object.type = PNL_TYPE_GMRES_SOLVER;
  o->object.parent_type = PNL_TYPE_OBJECT;
  o->object.label = pnl_gmres_solver_label;
  o->object.clone = NULL;
  o->object.copy = NULL;
  o->object.constructor = (NewFunc *) pnl_gmres_solver_new;
  return o;
}



/**
 * Create a new PnlGmresSolver pointer.
 *
 * @param Size : the size of temporary vectors use in the method
 * @param max_iter_  : maximum iteration number 
 * @param restart_   : restart parameters for gmres method
 * @param tolerance_ : relative error of iterative method
 *@return  a PnlGmresSolver pointer
 */
PnlGmresSolver* pnl_gmres_solver_create(int Size, int max_iter_,
                                        int restart_, double tolerance_)
{
  int i;
  PnlGmresSolver *Solver;
  if( (Solver = pnl_gmres_solver_new ()) == NULL ) return NULL;
  Solver->restart=MIN(MAX_RESTART - 1,restart_);
  Solver->s=pnl_vect_create_from_scalar(Solver->restart+1,0.0);
  Solver->cs=pnl_vect_create_from_scalar(Solver->restart+1,0.0);
  Solver->sn=pnl_vect_create_from_scalar(Solver->restart+1,0.0);
  Solver->w=pnl_vect_create_from_scalar(Size,0.0);
  Solver->r=pnl_vect_create_from_scalar(Size,0.0);
  Solver->H=pnl_mat_create_from_scalar(Solver->restart+1,Solver->restart,0.0);
  for(i=0;i<=Solver->restart+1;i++)
    Solver->v[i]=pnl_vect_create_from_scalar(Size,0.0);
  Solver->iter=pnl_iteration_base_create(max_iter_,tolerance_);
  Solver->iter_inner=pnl_iteration_base_create(restart_,tolerance_);
  return Solver;
}


  
/**
 * initialisation of the solver at the beginning of iterative method.
 *
 * @param Solver  : PnlGmresSolver  ptr
 * @param b       : PnlVect ptr 
 */
void pnl_gmres_solver_initialisation(PnlGmresSolver *Solver,
                                     const PnlVect * b)
{
  double normb=pnl_vect_norm_two(b);
  pnl_vect_clone(Solver->r,b);
  pnl_mat_set_all(Solver->H,0.0);
  Solver->beta=0.0;
  pnl_iteration_base_initialisation_norm(Solver->iter,normb);
  pnl_iteration_base_initialisation_norm(Solver->iter_inner,normb);
}

/**
 * destructor of iterative solver 
 *
 * @param Solver : PnlGmresSolver  ptr
 */
void pnl_gmres_solver_free(PnlGmresSolver ** Solver)
{
  int i;
  if (*Solver != NULL)
    {
      pnl_vect_free(&(*Solver)->s);
      pnl_vect_free(&(*Solver)->cs);
      pnl_vect_free(&(*Solver)->sn);
      pnl_vect_free(&(*Solver)->w);
      pnl_vect_free(&(*Solver)->r);
      pnl_mat_free(&(*Solver)->H);
      for(i=0;i<=(*Solver)->restart+1;i++)
        pnl_vect_free(&(*Solver)->v[i]);
      free((*Solver)->iter);
      free((*Solver)->iter_inner);
      free(*Solver);
      *Solver=NULL;
    } 
}

/**
 * Define rotation for GMRES 
 * @param dx,dy : double  
 * @param cs,sn : double  ptr
 */
static void GeneratePlaneRotation(double dx, double  dy, double *cs, double *sn)
{
  if (dy == 0.0) 
    {
      (*cs) = 1.0;
      (*sn) = 0.0;
    } 
  else if (fabs(dy) > fabs(dx)) 
    {
      double temp = dx / dy;
      (*sn) = 1.0 / sqrt( 1.0 + temp*temp );
      (*cs) = temp * (*sn);
    } 
  else 
    {
      double temp = dy / dx;
      (*cs) = 1.0 / sqrt( 1.0 + temp*temp );
      (*sn) = temp * (*cs);
    }
}
/**
 * Define rotation for GMRES 
 * @param dx,dy : double  ptr 
 * @param cs,sn : double 
 */
static void ApplyPlaneRotation(double *dx, double *dy, double cs, double sn)
{
  double temp  =  cs * (*dx) + sn * (*dy);
  (*dy) = -1.0*sn * (*dx) + cs * (*dy);
  (*dx) = temp;
}

/**
 * Solve the linear system
 * matrix_vector_product is the matrix_vector multiplication function
 * matrix_vector_product_PC is the preconditionner function
 * Matrix_Data & PC_Data is data to compute matrix_vector multiplication
 *
 * Solve :
 * Matrix_Data x =  b with preconditionner PC_Data
 *
 * @param matrix_vector_product : ptr function retrun PnlVect ptr and take 
 * void ptr and PnlVect ptr parameters 
 * @param Matrix_Data  : void ptr 
 * @param matrix_vector_product_PC : ptr function retrun PnlVect ptr and take 
 * void ptr and PnlVect ptr parameters 
 * @param PC_Data : void ptr
 * @param x       : PnlVect ptr 
 * @param b       : PnlVect ptr
 * @param Solver  : PnlGmresSolver ptr 
 * @return Int : error code.
 */
int pnl_gmres_solver_solve(void (* matrix_vector_product )(const void *,const PnlVect*,double,double,PnlVect*), 
                           const void * Matrix_Data,
                           void (*matrix_vector_product_PC)(const void *,const PnlVect*,double,double,PnlVect*), 
                           const void * PC_Data,
                           PnlVect * x, 
                           const PnlVect *b,
                           PnlGmresSolver *Solver)
{
  int i, j, k, p;
  double s_j;
  pnl_gmres_solver_initialisation(Solver,b);
  /* r = b - A x: initial residual */
  matrix_vector_product(Matrix_Data,x,-1.0,1.0,Solver->r);
  matrix_vector_product_PC(PC_Data,Solver->r,1.0,0.0,Solver->w);
  Solver->beta = pnl_vect_norm_two(Solver->w);
  while(!pnl_iteration_base_finished(Solver->iter,Solver->w))
    {
      pnl_vect_axpby(1.0 /Solver->beta, Solver->w, 0.0, Solver->v[0]);    
      pnl_vect_set_all(Solver->s,0.0);
      pnl_vect_set(Solver->s,0,Solver->beta);
      i = 0;
      pnl_iteration_base_begin(Solver->iter_inner);
      /*>>Restart loop */
      do {
        matrix_vector_product(Matrix_Data,Solver->v[i],1.0,0.0,Solver->r);
        matrix_vector_product_PC(PC_Data,Solver->r,1.0,0.0,Solver->w);
        for (k = 0; k <= i; k++) 
          {
            pnl_mat_set(Solver->H,k,i, pnl_vect_scalar_prod(Solver->w, Solver->v[k]));
            pnl_vect_axpby(-1.0*pnl_mat_get(Solver->H,k, i), Solver->v[k], 1.0, Solver->w);
          }
        pnl_mat_set(Solver->H,i+1, i,pnl_vect_norm_two(Solver->w));
        pnl_vect_axpby(1.0/pnl_mat_get(Solver->H,i+1, i), Solver->w, 0., Solver->v[i+1] ); 
        /*Arnoldi process finished */
        /* The   i+1* i matrix H has been computed   */
        for (k = 0; k < i; k++)
          ApplyPlaneRotation(pnl_mat_lget(Solver->H,k,i), 
                             pnl_mat_lget(Solver->H,k+1,i), 
                             pnl_vect_get(Solver->cs,k), 
                             pnl_vect_get(Solver->sn,k));
        GeneratePlaneRotation(pnl_mat_get(Solver->H,i,i), 
                              pnl_mat_get(Solver->H,i+1,i), 
                              pnl_vect_lget(Solver->cs,i), 
                              pnl_vect_lget(Solver->sn,i));
        ApplyPlaneRotation(pnl_mat_lget(Solver->H,i,i), 
                           pnl_mat_lget(Solver->H,i+1,i), 
                           pnl_vect_get(Solver->cs,i), 
                           pnl_vect_get(Solver->sn,i));
        ApplyPlaneRotation(pnl_vect_lget(Solver->s,i), 
                           pnl_vect_lget(Solver->s,i+1), 
                           pnl_vect_get(Solver->cs,i), 
                           pnl_vect_get(Solver->sn,i));
        pnl_iteration_base_increment(Solver->iter_inner);
        pnl_iteration_base_increment(Solver->iter);
        ++i;
      }while((!pnl_iteration_base_converged_norm(Solver->iter_inner,fabs(pnl_vect_get(Solver->s,i))))
             &&(i<Solver->iter_inner->max_iter));
      for(j = i-1;j >= 0;j--) 
        /* loop on the columns of the matrix */
        {
          pnl_vect_set(Solver->s,j,pnl_vect_get(Solver->s,j)/ pnl_mat_get(Solver->H,j,j));
          s_j = pnl_vect_get(Solver->s,j);
          for(p = j-1;p >= 0;p--) 
            /* loop on the lines of the matrix */
            pnl_vect_set(Solver->s,p,pnl_vect_get(Solver->s,p) - pnl_mat_get(Solver->H,p,j)*s_j);
        }
      for(j = 0;j < i;j++)
        pnl_vect_axpby(pnl_vect_get(Solver->s,j), Solver->v[j], 1., x);
      
      matrix_vector_product(Matrix_Data,x,1.0,0.0,Solver->r);
      /*r =A*x;  */
      pnl_vect_axpby(1.0, b, -1.0, Solver->r);
      /*r = b-M*x; */
      matrix_vector_product_PC(PC_Data,Solver->r,1.0,0.0,Solver->w);
      Solver->beta = pnl_vect_norm_two(Solver->w);
    }
  return pnl_iteration_base_error_code(Solver->iter);
}

/**
 * Matrix_Vector Multiplication with void ptr and PnlVect ptr 
 * just do a cast of void ptr to PnlMat ptr
 *
 * @param mat : void ptr 
 * @param vec : PnlVect ptr
 * @param a : PnlVect ptr
 * @param b : PnlVect ptr
 * @param lhs : PnlVect ptr
 * lhs = a*(mat*vec)+b*lhs. 
 */ 
static void pnl_mat_mult_vect_applied(const void *mat, const PnlVect *vec,double a ,double b,PnlVect *lhs)
{pnl_mat_lAxpby(a,(PnlMat*)mat,vec,b,lhs);}

/**
 * Solve the linear system
 * Matrix x =  b with preconditionner PC
 *
 * @param Matrix :PnlMat  ptr 
 * @param PC     :PnlMat  ptr 
 * @param x       : PnlVect ptr 
 * @param b       : PnlVect ptr
 * @param Solver  : PnlCgSolver ptr 
 * @return Int : error code.
 */
int pnl_mat_cg_solver_solve(const PnlMat * Matrix,
                            const PnlMat * PC,
                            PnlVect * x, 
                            const PnlVect *b,
                            PnlCgSolver * Solver)
{ return 
    pnl_cg_solver_solve(pnl_mat_mult_vect_applied,Matrix,pnl_mat_mult_vect_applied,PC,x,b,Solver);}


/**
 * Solve the linear system
 * Matrix x =  b with preconditionner PC
 *
 * @param Matrix :PnlMat  ptr 
 * @param PC     :PnlMat  ptr 
 * @param x       : PnlVect ptr 
 * @param b       : PnlVect ptr
 * @param Solver  : PnlBicgSolver ptr 
 * @return Int : error code.
 */
int pnl_mat_bicg_solver_solve(const PnlMat * Matrix,
                              const PnlMat * PC,
                              PnlVect * x, 
                              const PnlVect *b,
                              PnlBicgSolver * Solver)
{return pnl_bicg_solver_solve(pnl_mat_mult_vect_applied,Matrix,pnl_mat_mult_vect_applied,PC,x,b,Solver);}

/**
 * Solve the linear system
 * Matrix x =  b with preconditionner PC
 *
 * @param Matrix :PnlMat  ptr 
 * @param PC     :PnlMat  ptr 
 * @param x       : PnlVect ptr 
 * @param b       : PnlVect ptr
 * @param Solver  : PnlGmresSolver ptr 
 * @return Int : error code.
 */
int pnl_mat_gmres_solver_solve(const PnlMat * Matrix,
                               const PnlMat * PC,
                               PnlVect * x, 
                               const PnlVect *b,
                               PnlGmresSolver * Solver)
{return pnl_gmres_solver_solve(pnl_mat_mult_vect_applied,Matrix,pnl_mat_mult_vect_applied,PC,x,b,Solver);}


#undef DEBUG_SOLVER
#undef INFO_SOLVER
#undef MAX_RESTART
