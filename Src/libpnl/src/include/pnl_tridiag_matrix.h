#ifndef TRID_MATRIX_H
#define TRID_MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl_matrix.h"

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



extern void pnl_progonka(const double low, 
                         const double diag, const double up, 
                         const PnlVect * rhs, PnlVect * lhs);

/**
 * \defgroup PnlTriDiagMat Tridiagonal Matrix
 */

typedef struct PnlTriDiagMat{
  int size; /*!< number of rows, the matrix must be square */
  double *D; /*!< diagonal elements */
  double *DU; /*!< upper diagonal elements */
  double *DL; /*!< lower diagonal elements */
} PnlTriDiagMat;


extern void pnl_tridiagmat_free(PnlTriDiagMat **m);
extern int pnl_tridiagmat_resize(PnlTriDiagMat *v, int size);
extern PnlTriDiagMat* pnl_tridiagmat_create(int size);
extern PnlTriDiagMat* pnl_tridiagmat_create_from_double(int size, double x);
extern PnlTriDiagMat* pnl_tridiagmat_create_from_two_double(int size, double x, double y);
extern PnlTriDiagMat* pnl_tridiagmat_create_from_ptr(int size, const double* DL, 
						   const double* D, const double* DU);
extern PnlTriDiagMat* pnl_tridiagmat_create_from_mat(const PnlMat * mat);
extern PnlMat* pnl_tridiagmat_to_mat(const PnlTriDiagMat * mat);
extern PnlTriDiagMat* pnl_tridiagmat_copy(const PnlTriDiagMat * mat);
extern void pnl_tridiagmat_clone(PnlTriDiagMat *clone, const PnlTriDiagMat * mat);

extern void pnl_tridiagmat_print(const PnlTriDiagMat *M);
extern void pnl_tridiagmat_fprint(FILE *fic, const PnlTriDiagMat *M);

extern void pnl_tridiagmat_plus_tridiagmat(PnlTriDiagMat *lhs, const PnlTriDiagMat *rhs); 
extern void pnl_tridiagmat_minus_tridiagmat(PnlTriDiagMat *lhs, const PnlTriDiagMat *rhs); 
extern void pnl_tridiagmat_plus_double(PnlTriDiagMat *lhs, double x);
extern void pnl_tridiagmat_minus_double(PnlTriDiagMat *lhs, double x);
extern void pnl_tridiagmat_mult_double(PnlTriDiagMat *lhs, double x); 
extern void pnl_tridiagmat_div_double(PnlTriDiagMat *lhs, double x); 
extern void pnl_tridiagmat_mult_tridiagmat_term(PnlTriDiagMat *lhs, const PnlTriDiagMat *rhs);
extern void pnl_tridiagmat_div_tridiagmat_term(PnlTriDiagMat *lhs, const PnlTriDiagMat *rhs);

extern PnlVect* pnl_tridiagmat_mult_vect(const PnlTriDiagMat *mat,const PnlVect *vec);
extern void pnl_tridiagmat_mult_vect_inplace(PnlVect *lhs, const PnlTriDiagMat *mat, const PnlVect *rhs);
extern void pnl_tridiagmat_lAxpby(double l, const PnlTriDiagMat *A, const PnlVect *x, double b, PnlVect * y);
extern double pnl_tridiagmat_scalar_prod(const PnlVect *lhs, const PnlTriDiagMat *mat,const PnlVect *rhs);
extern int pnl_tridiagmat_lu_syslin (PnlVect *lhs, const PnlTriDiagMat *M,const PnlVect *rhs);
extern int pnl_tridiagmat_lu_syslin_inplace (const PnlTriDiagMat *M, PnlVect *rhs);


extern void pnl_tridiagmat_set(PnlTriDiagMat *v, int d, int up, double x);
extern double pnl_tridiagmat_get(const PnlTriDiagMat *self, int d, int up);
extern double* pnl_tridiagmat_lget(PnlTriDiagMat *self, int d, int up);
/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* TRID_MATRIX_H */
