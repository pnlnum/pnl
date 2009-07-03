#ifndef SOLVER_ITERATIV_H
#define SOLVER_ITERATIV_H

#include "pnl_mathtools.h"
#include "pnl_matrix.h"

#ifndef MAX_RESTART 
#define MAX_RESTART   50 
#endif

/*
 * \defgroup IterationBase Iterativ Solvers
 * \author David Pommier 
 * \brief Iterative Solvers like Conjugate Gradient, BICGStab and GMRES.
 * \date October 2008 
 * \version 0.1
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

/* When you repeatedly use Iterative solvers, do not malloc each time */
typedef struct PnlCGSolver{
  PnlVect * r;
  PnlVect * z;
  PnlVect * p;
  PnlVect * q;
  double rho;
  double oldrho;
  double beta;
  double alpha;
  PnlIterationBase * iter;
} PnlCGSolver;

typedef struct PnlBICGSolver{ 
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
} PnlBICGSolver;

typedef struct PnlGMRESSolver{ 
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
} PnlGMRESSolver;


extern PnlCGSolver* pnl_cg_solver_create(int Size,int max_iter_, double tolerance_);
extern void pnl_cg_solver_initialisation(PnlCGSolver * Solver,const PnlVect * b);
extern void pnl_cg_solver_free(PnlCGSolver ** Solver);


extern int pnl_cg_solver_solve(void(* matrix_vector_product )(const void *,const PnlVect*,
                                                              const double,const double,PnlVect*), 
                               const void * Matrix_Data,
                               void (*matrix_vector_product_PC)(const void *,const PnlVect*,
                                                                const double,const double,PnlVect*), 
                               const void * PC_Data,
                               PnlVect * x, 
                               const PnlVect * b, 
                               PnlCGSolver * Solver);

extern int pnl_mat_cg_solver_solve(const PnlMat *Matrix,
                                   const PnlMat *PC,
                                   PnlVect *x, 
                                   const PnlVect *b,
                                   PnlCGSolver * Solver);


 
extern PnlBICGSolver* pnl_bicg_solver_create(int Size,int max_iter_, double tolerance_);
extern void pnl_bicg_solver_initialisation(PnlBICGSolver * Solver,const PnlVect * b);
extern void pnl_bicg_solver_free(PnlBICGSolver ** Solver);
extern int
pnl_bicg_solver_solve(void (* matrix_vector_product )(const void *,const PnlVect*,
                                                      const double,const double,PnlVect*), 
                      const void * Matrix_Data,
                      void (*matrix_vector_product_PC)(const void *,const PnlVect*,
                                                       const double,const double,PnlVect*), 
                      const void * PC_Data,
                      PnlVect * x, 
                      const PnlVect *b,
                      PnlBICGSolver * Solver);
extern int pnl_mat_bicg_solver_solve(const PnlMat *Matrix,
                                     const PnlMat *PC,
                                     PnlVect *x, 
                                     const PnlVect *b,
                                     PnlBICGSolver * Solver);



extern PnlGMRESSolver* pnl_gmres_solver_create(int Size,
                                               int max_iter_, 
                                               int restart_,
                                               double tolerance_);

extern void pnl_gmres_solver_initialisation(PnlGMRESSolver *Solver,
                                            const PnlVect * b);
extern void pnl_gmres_solver_free(PnlGMRESSolver ** Solver);

extern int
pnl_gmres_solver_solve(void (* matrix_vector_product )(const void *,const PnlVect*,
                                                       const double,const double,PnlVect*), 
                       const void * Matrix_Data,
                       void (*matrix_vector_product_PC)(const void *,const PnlVect*,
                                                        const double,const double,PnlVect*), 
                       const void * PC_Data,
                       PnlVect * x, 
                       const PnlVect *b,
                       PnlGMRESSolver *Solver);

int pnl_mat_gmres_solver_solve(const PnlMat * Matrix,
                               const PnlMat * PC,
                               PnlVect * x, 
                               const PnlVect *b,
                               PnlGMRESSolver * Solver);
/*@}*/
#endif /* solver_iterativ_H */
