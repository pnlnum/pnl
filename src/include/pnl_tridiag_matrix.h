#ifndef _PNL_TRIDIAG_MATRIX_H 
#define _PNL_TRIDIAG_MATRIX_H 

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl_matrix.h"

#ifndef PNL_RANGE_CHECK_OFF
#define CheckIndexTridiagMat(v,i,j)                                                     \
  if (i<0 || i>=(v)->size || j<-1 || j>1 || (i==0 && j==1) || (i=(v)->size-1 && j==-1)) \
    {                                                                                   \
      perror("index out of range"); abort();                                            \
    }
#define CheckTridiagMatMatch(lhs, rhs)                                                  \
  if ((lhs)->size != (rhs)->size )                                                      \
    {                                                                                   \
      perror("non compatible dimensions"); abort();                                     \
    }
#define CheckTridiagMatVectIsCompatible(mat, vect)                                      \
  if ((mat)->size != (vect)->size)                                                      \
    {                                                                                   \
      perror("non compatible dimensions"); abort();                                     \
    }
#else
#define CheckIndexTridiagMat(v,i,j) {}                          
#define CheckTridiagMatMatch(lhs, rhs) {}
#define CheckTridiagMatVectIsCompatible(mat, vect) {}
#endif /* PNL_RANGE_CHECK_OFF */



extern void pnl_progonka(const double low, 
                         const double diag, const double up, 
                         const PnlVect * rhs, PnlVect * lhs);

/**
 * \defgroup PnlTridiagMat Tridiagonal Matrix
 *
 * Three arrays are used to store the three main diagonals. The arrays used to
 * store the lower and upper diagonals are of size n-1 for a n x n matrix.
 * Tridiagonal matrices must be square.
 */

typedef struct 
{
  int size; /*!< number of rows, the matrix must be square */
  double *D; /*!< diagonal elements */
  double *DU; /*!< upper diagonal elements */
  double *DL; /*!< lower diagonal elements */
} PnlTridiagMat;


extern void pnl_tridiag_mat_free(PnlTridiagMat **m);
extern int pnl_tridiag_mat_resize(PnlTridiagMat *v, int size);
extern PnlTridiagMat* pnl_tridiag_mat_create(int size);
extern PnlTridiagMat* pnl_tridiag_mat_create_from_double(int size, double x);
extern PnlTridiagMat* pnl_tridiag_mat_create_from_two_double(int size, double x, double y);
extern PnlTridiagMat* pnl_tridiag_mat_create_from_ptr(int size, const double* DL, 
						   const double* D, const double* DU);
extern PnlTridiagMat* pnl_tridiag_mat_create_from_mat(const PnlMat * mat);
extern PnlMat* pnl_tridiag_mat_to_mat(const PnlTridiagMat * mat);
extern PnlTridiagMat* pnl_tridiag_mat_copy(const PnlTridiagMat * mat);
extern void pnl_tridiag_mat_clone(PnlTridiagMat *clone, const PnlTridiagMat * mat);

extern void pnl_tridiag_mat_map_inplace(PnlTridiagMat *lhs, double(*f)(double));
extern void pnl_tridiag_mat_map_tridiag_mat(PnlTridiagMat *lhs, const PnlTridiagMat *rhs, double(*f)(double,double));
extern void pnl_tridiag_mat_print(const PnlTridiagMat *M);
extern void pnl_tridiag_mat_fprint(FILE *fic, const PnlTridiagMat *M);

extern void pnl_tridiag_mat_plus_tridiag_mat(PnlTridiagMat *lhs, const PnlTridiagMat *rhs); 
extern void pnl_tridiag_mat_minus_tridiag_mat(PnlTridiagMat *lhs, const PnlTridiagMat *rhs); 
extern void pnl_tridiag_mat_plus_double(PnlTridiagMat *lhs, double x);
extern void pnl_tridiag_mat_minus_double(PnlTridiagMat *lhs, double x);
extern void pnl_tridiag_mat_mult_double(PnlTridiagMat *lhs, double x); 
extern void pnl_tridiag_mat_div_double(PnlTridiagMat *lhs, double x); 
extern void pnl_tridiag_mat_mult_tridiag_mat_term(PnlTridiagMat *lhs, const PnlTridiagMat *rhs);
extern void pnl_tridiag_mat_div_tridiag_mat_term(PnlTridiagMat *lhs, const PnlTridiagMat *rhs);

extern PnlVect* pnl_tridiag_mat_mult_vect(const PnlTridiagMat *mat,const PnlVect *vec);
extern void pnl_tridiag_mat_mult_vect_inplace(PnlVect *lhs, const PnlTridiagMat *mat, const PnlVect *rhs);
extern void pnl_tridiag_mat_lAxpby(double l, const PnlTridiagMat *A, const PnlVect *x, double b, PnlVect * y);
extern double pnl_tridiag_mat_scalar_prod(const PnlVect *lhs, const PnlTridiagMat *mat,const PnlVect *rhs);
extern int pnl_tridiag_mat_syslin (PnlVect *lhs, const PnlTridiagMat *M,const PnlVect *rhs);
extern int pnl_tridiag_mat_syslin_inplace (const PnlTridiagMat *M, PnlVect *rhs);


extern void pnl_tridiag_mat_set(PnlTridiagMat *v, int d, int up, double x);
extern double pnl_tridiag_mat_get(const PnlTridiagMat *self, int d, int up);
extern double* pnl_tridiag_mat_lget(PnlTridiagMat *self, int d, int up);
/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_TRIDIAG_MATRIX_H  */
