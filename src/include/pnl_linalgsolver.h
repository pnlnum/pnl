#ifndef _PNL_LINALGSOLVER_H 
#define _PNL_LINALGSOLVER_H 


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl_mathtools.h"
#include "pnl_matrix.h"

#ifndef MAX_RESTART 
#define MAX_RESTART   50 
#endif

/**
 * \defgroup IterationBase Iterative Solvers
 * \brief Iterative Solvers like Conjugate Gradient, BICGStab and GMRES.
 */
/*@{*/
typedef struct PnlIterationBase{
  int iteration;
  int max_iter;
  double normb;
  double tol_;
  double resid;
  int error;
  /* char *  err_msg; */
} PnlIterationBase;

/* When you repeatedly use iterative solvers, do not malloc each time */
typedef struct PnlCgSolver{
  PnlVect * r;
  PnlVect * z;
  PnlVect * p;
  PnlVect * q;
  double rho;
  double oldrho;
  double beta;
  double alpha;
  PnlIterationBase * iter;
} PnlCgSolver;

typedef struct PnlBicgSolver{ 
  double rho_1, rho_2, alpha, beta, omega;
  PnlVect * p;
  PnlVect * phat;
  PnlVect * s;
  PnlVect * shat;
  PnlVect * t;
  PnlVect * v;
  PnlVect * r;
  PnlVect *  rtilde;
  PnlIterationBase * iter;
} PnlBicgSolver;

typedef struct PnlGmresSolver{ 
  int restart;
  double beta;
  PnlVect * s;
  PnlVect * cs;
  PnlVect * sn;
  PnlVect * w;
  PnlVect * r;
  PnlMat * H;
  PnlVect * v[MAX_RESTART];
  PnlIterationBase *iter;
  PnlIterationBase *iter_inner;
} PnlGmresSolver;


extern PnlCgSolver* pnl_cg_solver_create(int Size,int max_iter_, double tolerance_);
extern void pnl_cg_solver_initialisation(PnlCgSolver * Solver,const PnlVect * b);
extern void pnl_cg_solver_free(PnlCgSolver ** Solver);


extern int pnl_cg_solver_solve(void(* matrix_vector_product )(const void *,const PnlVect*,
                                                              double,double,PnlVect*), 
                               const void * Matrix_Data,
                               void (*matrix_vector_product_PC)(const void *,const PnlVect*,
                                                                double,double,PnlVect*), 
                               const void * PC_Data,
                               PnlVect * x, 
                               const PnlVect * b, 
                               PnlCgSolver * Solver);

extern int pnl_mat_cg_solver_solve(const PnlMat *Matrix,
                                   const PnlMat *PC,
                                   PnlVect *x, 
                                   const PnlVect *b,
                                   PnlCgSolver * Solver);


 
extern PnlBicgSolver* pnl_bicg_solver_create(int Size,int max_iter_, double tolerance_);
extern void pnl_bicg_solver_initialisation(PnlBicgSolver * Solver,const PnlVect * b);
extern void pnl_bicg_solver_free(PnlBicgSolver ** Solver);
extern int
pnl_bicg_solver_solve(void (* matrix_vector_product )(const void *,const PnlVect*,
                                                      double,double,PnlVect*), 
                      const void * Matrix_Data,
                      void (*matrix_vector_product_PC)(const void *,const PnlVect*,
                                                       double,double,PnlVect*), 
                      const void * PC_Data,
                      PnlVect * x, 
                      const PnlVect *b,
                      PnlBicgSolver * Solver);
extern int pnl_mat_bicg_solver_solve(const PnlMat *Matrix,
                                     const PnlMat *PC,
                                     PnlVect *x, 
                                     const PnlVect *b,
                                     PnlBicgSolver * Solver);



extern PnlGmresSolver* pnl_gmres_solver_create(int Size,
                                               int max_iter_, 
                                               int restart_,
                                               double tolerance_);

extern void pnl_gmres_solver_initialisation(PnlGmresSolver *Solver,
                                            const PnlVect * b);
extern void pnl_gmres_solver_free(PnlGmresSolver ** Solver);

extern int
pnl_gmres_solver_solve(void (* matrix_vector_product )(const void *,const PnlVect*,
                                                       double,double,PnlVect*), 
                       const void * Matrix_Data,
                       void (*matrix_vector_product_PC)(const void *,const PnlVect*,
                                                        double,double,PnlVect*), 
                       const void * PC_Data,
                       PnlVect * x, 
                       const PnlVect *b,
                       PnlGmresSolver *Solver);

int pnl_mat_gmres_solver_solve(const PnlMat * Matrix,
                               const PnlMat * PC,
                               PnlVect * x, 
                               const PnlVect *b,
                               PnlGmresSolver * Solver);
/*@}*/
#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_LINALGSOLVER_H  */
