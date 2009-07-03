#ifndef TRIDIAG_MATRIX_H
#define TRIDIAG_MATRIX_H

#include "pnl_matrix.h"

#define EPSILON 1e-10  
#ifndef PNL_RANGE_CHECK_OFF
#define CheckIndexTriDiagMat(v,i,j){             \
    if (i>=v->size || abs(j)>1  || i<0 )       \
      {perror("index out of range"); abort();\
      }}
#define CheckTriDiagMatMatch(lhs, rhs){          \
    if ((lhs)->size != (rhs)->size )             \
      {perror("non compatible dimensions"); abort();\
}}
#define CheckTriDiagMatVectIsCompatible(mat, vect){            \
    if((mat)->size != (vect)->size)                        \
      {perror("non compatible dimensions"); abort();}}
#else
#define CheckIndexTriDiagMat(v,i,j) {}                          
#define CheckTriDiagMatMatch(lhs, rhs) {}
#define CheckTriDiagMatVectIsCompatible(mat, vect){}
#endif /* PNL_RANGE_CHECK_OFF */


/**
 * \defgroup PnlTriDiagMat Tri diagonal Matrix structure
 * \author David Pommier 
 * \brief A Tridiagonal Matrix structure for Numerical Algorithm to solve PDE,
 * (strongly inspired by work of J P.Chancelier & al in NSP software)
 * \date July 2008 
 * \version 0.1
 * \f[ D[i] = A(i,i),\quad  D_{up}[i] = A(i,i+1),\quad D_{down}[i] = A(i,i-1) \f]
 */

typedef struct PnlTriDiagMat{
  int size; /*!< diagonal dimension product  */
  int owner; /*!< 1 if the structure owns its array pointer */
  double *diag; /*!< pointer to store diagonal*/
  double *diag_up; /*!< pointer to store up diagonal*/
  double *diag_down; /*!< pointer to store down diagonal*/
} PnlTriDiagMat;


extern void pnl_tridiagmat_free(PnlTriDiagMat **m);
extern int pnl_tridiagmat_resize(PnlTriDiagMat *v, int size);
extern PnlTriDiagMat* pnl_tridiagmat_create(int size);
extern PnlTriDiagMat* pnl_tridiagmat_create_from_double(int size, double x);
extern PnlTriDiagMat* pnl_tridiagmat_create_from_two_double(int size, double x, double y);
extern PnlTriDiagMat* pnl_tridiagmat_create_from_ptr(int size, const double* diag, 
						   const double* diag_up, const double* diag_down);
extern PnlTriDiagMat* pnl_tridiagmat_create_from_matrix(const PnlMat * mat);
extern void pnl_tridiagmat_print(const PnlTriDiagMat *M);
extern void pnl_tridiagmat_fprint(FILE *fic, const PnlTriDiagMat *M);
extern void pnl_tridiagmat_plus_tridiagmat(PnlTriDiagMat *lhs, const PnlTriDiagMat *rhs); /*lhs+=rhs*/
extern void pnl_tridiagmat_minus_tridiagmat(PnlTriDiagMat *lhs, const PnlTriDiagMat *rhs); /*lhs-=rhs*/
extern void pnl_tridiagmat_mult_double(PnlTriDiagMat *lhs, double x); /*lhs*=x*/
extern void pnl_tridiagmat_div_double(PnlTriDiagMat *lhs, double x); /*lhs*=x*/
extern void pnl_tridiagmat_mult_tridiagmat_term(PnlTriDiagMat *lhs, 
					   const PnlTriDiagMat *rhs);/*lhs(i,j) *=rhs(i,j)*/
extern void pnl_tridiagmat_div_tridiagmat_term(PnlTriDiagMat *lhs, 
					  const PnlTriDiagMat *rhs);/*lhs(i,j) /=rhs(i,j)*/

extern PnlVect* pnl_tridiagmat_mult_vect(const PnlTriDiagMat *mat,
					   const PnlVect *vec);/*mat*vec*/
extern 
void pnl_tridiagmat_mult_vect_inplace(PnlVect *lhs, const PnlTriDiagMat *mat,
				const PnlVect *rhs);/*lhs=mat*rhs */
extern 
void pnl_tridiagmat_mult_vect_operation(const PnlTriDiagMat *mat, const PnlVect *rhs ,const double a, 
				       const double b, PnlVect * lhs);/*lhs=a*mat*rhs +b*lhs*/
extern 
void pnl_tridiagmat_lu_syslin (PnlVect *lhs, const PnlTriDiagMat *M,const PnlVect *rhs);
/* solve M lhs = rhs */

extern void pnl_tridiagmat_set(PnlTriDiagMat *v, int d, int up, double x);
extern double pnl_tridiagmat_get(const PnlTriDiagMat *self, int d, int up);
extern double* pnl_tridiagmat_lget(PnlTriDiagMat *self, int d, int up);
/*@}*/

#endif /* TRIDIAG_MATRIX_H */


