#ifndef PDE_TOOLS_H
#define PDE_TOOLS_H

#include "pnl_vector.h"
#include "pnl_matrix_uint.h"

/**
 * \defgroup PnlPDEBoundary  to translate domaine from X0,X1 to [0,1]
 */
/*@{*/
typedef struct {
    double X0; /*!< left point */
    double H; /* !< Step */
}PnlPDEBoundary;

extern PnlPDEBoundary pnl_pde_boundary_create(double X0,double X1);
extern double pnl_pde_boundary_real_variable(const PnlPDEBoundary BP ,double X);
extern double pnl_pde_boundary_Unit_interval(const PnlPDEBoundary BP ,double X);
/*@}*/
 

/**
 * \defgroup PnlPDEDimBoundary Vector On boundary
 */
/*@{*/
typedef struct PnlPDEDimBoundary{
  PnlPDEBoundary * array;
/*!< pointer to store the data */
} PnlPDEDimBoundary;

extern PnlPDEDimBoundary* pnl_pde_dim_boundary_create_from_int(int dim);

extern PnlPDEDimBoundary* 
pnl_pde_dim_boundary_create(const PnlVect * X0,
			       const PnlVect * X1);
extern void pnl_pde_dim_boundary_free(PnlPDEDimBoundary **v);

extern double 
pnl_pde_dim_boundary_eval_from_unit(double(*f)(const PnlVect* ),
				       const PnlPDEDimBoundary * BP,
				       const PnlVect * X);  
extern void
pnl_pde_dim_boundary_from_unit_to_real_variable(const PnlPDEDimBoundary * BP,
						   PnlVect * X);

extern double 
pnl_pde_dim_boundary_get_step(const PnlPDEDimBoundary * BP,
				 int i);
/*@}*/

typedef struct {
  double current_step;
  int current_index;
  PnlVect * time; 
  int is_tuned;
}PnlPDETimeGrid;

extern PnlPDETimeGrid * pnl_pde_time_homogen_grid(const double T,
							const int N_T);
extern void pnl_pde_time_grid_free(PnlPDETimeGrid **TG);
extern void pnl_pde_time_start(PnlPDETimeGrid * TG);
extern int pnl_pde_time_grid_increase(PnlPDETimeGrid * TG);
extern double pnl_pde_time_grid_step(const PnlPDETimeGrid * TG);
extern double pnl_pde_time_grid_time(const PnlPDETimeGrid * TG);


#endif

