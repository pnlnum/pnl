#ifndef BAND_MATRIX_H
#define BAND_MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl_matrix.h"

#define EPSILON 1e-10  
#ifndef PNL_RANGE_CHECK_OFF
#define CheckIndexBandMat(v,i,j){             \
    if (i>=v->size || abs(j)>1  || i<0 )       \
      {perror("index out of range"); abort();\
      }}
#define CheckBandMatMatch(lhs, rhs){          \
    if ((lhs)->size != (rhs)->size )             \
      {perror("non compatible dimensions"); abort();\
}}
#define CheckBandMatVectIsCompatible(mat, vect){            \
    if((mat)->size != (vect)->size)                        \
      {perror("non compatible dimensions"); abort();}}
#else
#define CheckIndexBandMat(v,i,j) {}                          
#define CheckBandMatMatch(lhs, rhs) {}
#define CheckBandMatVectIsCompatible(mat, vect){}
#endif /* PNL_RANGE_CHECK_OFF */

/**
 * \defgroup PnlBandMatrix Band Matrix Structure
 * \brief Band Matrix structure for Numerical Algorithm to solve PDE, (strongly inspired by work of F.Hecht on RMN class).
 * \version 0.1
 * \warning Not tested at yet.
 * \f[  D(i) = A(i,i), \f] 
 * \f[ Lo(k) = A(i,j),\quad   j < i \textit{ with: }  \quad  pL(i)<= k < pL(i+1) \textit{ and }
 i-j = pL(i+1)-k. \f]                                                      \
 * \f[ Up(k) = A(i,j),\quad   i < j \textit{ with: }   \quad pU(j)<= k < pU(j+1) \textit{
  and } j-i = pU(i+1)-k. \f]
 * remark:
 * \f$ pL = pU \f$ in most of case 
 * if \f$L = U\f$ then symetric matrix} 
 */
/*@{*/
typedef struct _band_matrix PnlBandMatrix;
typedef enum {
  FactorizationNO=0,
  FactorizationCholeski=1,
  FactorizationCrout=2,
  FactorizationLU=3
}FactorizationType;

struct _band_matrix{
  int n; /*!< size of row */ 
  int m; /*!< size of col */ 
  double *D;  /*!< diagonal vector */  
  double *Up; /*!< Triangular up matrix, only no nul coeffieciens */ 
  double *Lo; /*!< Triangular down matrix, only no nul coeffieciens */
  int    *pU; /*!< Triangular up profil: \f$ Up(k)=A(i,j),\ i<j \f$ with: \f$ pU(j)<=k<pU(j+1)\f$ & \f$i=pU(i+1)-k \f$*/ 
  int    *pL; /*!< Triangular down profil:\f$ Lo(k)=A(i,j), \ j<i \f$
                 with:\f$ pL(i)<=k<pL(i+1)\f$ & \f$ j=pL(i+1)-k\f$ */
  FactorizationType typefac;
  int owner; /*!< 1 if the structure owns its array pointer */
};

/*PnlBandMatrix* pnl_band_matrix_create(const int  n,const double *a); */
extern PnlBandMatrix* pnl_band_matrix_create(const int  n,int band);
extern PnlBandMatrix* pnl_band_matrix_create_from_full(const PnlMat *PM,int band);
extern void pnl_band_matrix_free(PnlBandMatrix ** M);
extern void pnl_band_matrix_add(PnlBandMatrix * M,int i,int j,double x);
extern void pnl_band_matrix_set(PnlBandMatrix * M,int i,int j,double x);
extern void pnl_band_matrix_set_double(PnlBandMatrix*  M,double x);
extern void pnl_band_matrix_map_inplace(PnlBandMatrix *lhs,double(*f)(double ));
extern void pnl_band_matrix_plus_double(PnlBandMatrix *lhs , double x);
extern void pnl_band_matrix_minus_double(PnlBandMatrix *lhs , double x);
extern void pnl_band_matrix_mult_double(PnlBandMatrix *lhs , double x);
extern void pnl_band_matrix_div_double(PnlBandMatrix *lhs , double x);
extern void pnl_band_matrix_plus_mat(PnlBandMatrix *lhs, const PnlBandMatrix *rhs);
extern void pnl_band_matrix_minus_mat(PnlBandMatrix *lhs, const PnlBandMatrix *rhs);
extern void pnl_band_matrix_inv_term(PnlBandMatrix *lhs);
extern void pnl_band_matrix_div_mat_term(PnlBandMatrix *lhs, const PnlBandMatrix *rhs);
extern void pnl_band_matrix_mult_mat_term(PnlBandMatrix *lhs, const PnlBandMatrix *rhs);
extern void pnl_bnd_matrix_clone(PnlBandMatrix * clone,const PnlBandMatrix * v);
extern void pnl_bnd_matrix_store_infull(const PnlBandMatrix *BM,PnlMat *M);

extern PnlVect* pnl_band_matrix_mult_vect(const PnlBandMatrix *mat,const PnlVect *vec);/*mat*vec*/
extern void pnl_band_matrix_mult_vect_inplace(PnlVect *lhs, const PnlBandMatrix *mat, const PnlVect *rhs);/*lhs=mat*rhs */
extern void pnl_band_matrix_lAxpby(double l, const PnlBandMatrix *A, const PnlVect *x, double b, PnlVect * y);/*lhs=l*A*x +b*y*/
extern double pnl_band_matrix_prod_scale(const PnlVect *lhs, const PnlBandMatrix *mat,const PnlVect *rhs);
extern void pnl_band_matrix_lu_syslin (PnlVect *lhs, const PnlBandMatrix *M,const PnlVect *rhs);/* solve M lhs = rhs */


/*,FactorizationType tf = FactorizationNO);  */
extern PnlBandMatrix* pnl_band_matrix_transpose(const PnlBandMatrix* M);
extern PnlBandMatrix* pnl_band_matrix_Low(const PnlBandMatrix* M);
extern PnlBandMatrix* pnl_band_matrix_Up(const PnlBandMatrix* M);
extern PnlBandMatrix* pnl_band_matrix_Tran_Low(const PnlBandMatrix* M);
extern PnlBandMatrix* pnl_band_matrix_Tran_Up(const PnlBandMatrix* M);
extern PnlBandMatrix* pnl_band_matrix_Diag(const PnlBandMatrix* M);
extern PnlBandMatrix* pnl_band_matrix_Low_Diag(const PnlBandMatrix* M);
extern PnlBandMatrix* pnl_band_matrix_Up_Diag(const PnlBandMatrix* M);
extern PnlBandMatrix* pnl_band_matrix_Tran_Low_Diag(const PnlBandMatrix* M);
extern PnlBandMatrix* pnl_band_matrix_Tran_Up_Diag(const PnlBandMatrix* M);

extern double pnl_band_matrix_conditionning(const PnlBandMatrix *M);

extern void pnl_band_matrix_solve_syslin_inplace(PnlBandMatrix * M, PnlVect *x);
extern void pnl_band_matrix_cholesky(PnlBandMatrix * M, double eps);
extern void pnl_band_matrix_crout(PnlBandMatrix * M, double eps);
extern void pnl_band_matrix_lu(PnlBandMatrix * M, double eps);
extern void pnl_band_matrix_solve(PnlBandMatrix * M, PnlVect *x,const PnlVect *b);
/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* BAND_MATRIX_H */
