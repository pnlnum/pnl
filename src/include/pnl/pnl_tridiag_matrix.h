#ifndef _PNL_TRIDIAG_MATRIX_H 
#define _PNL_TRIDIAG_MATRIX_H 

#include "pnl/pnl_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


#ifndef PNL_RANGE_CHECK_OFF
#define CheckIndexTridiagMat(v,d,up)                                                     \
  if (d<0 || d>=(v)->size || up<-1 || up>1 || (d==0 && up==-1) || (d==(v)->size-1 && up==1)) \
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
/* @{ */
typedef struct _PnlTridiagMatObject PnlTridiagMatObject;
typedef struct _PnlTridiagMat PnlTridiagMat;

struct _PnlTridiagMat
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows any PnlTridiagMat pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int size; /*!< number of rows, the matrix must be square */
  double *D; /*!< diagonal elements */
  double *DU; /*!< upper diagonal elements */
  double *DL; /*!< lower diagonal elements */
};

struct _PnlTridiagMatObject
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows any PnlTridiagMat pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int size; /*!< number of rows, the matrix must be square */
  void *D; /*!< diagonal elements */
  void *DU; /*!< upper diagonal elements */
  void *DL; /*!< lower diagonal elements */
};

extern PnlTridiagMatObject* pnl_tridiag_mat_object_new();
extern void pnl_tridiag_mat_object_free (PnlTridiagMatObject **BM);
extern PnlTridiagMat* pnl_tridiag_mat_new();
extern PnlTridiagMat* pnl_tridiag_mat_create(int size);
extern PnlTridiagMat* pnl_tridiag_mat_create_from_scalar(int size, double x);
extern PnlTridiagMat* pnl_tridiag_mat_create_from_two_scalar(int size, double x, double y);
extern PnlTridiagMat* pnl_tridiag_mat_create_from_ptr(int size, const double* DL, 
						   const double* D, const double* DU);
extern PnlTridiagMat* pnl_tridiag_mat_create_from_mat(const PnlMat * mat);
extern PnlMat* pnl_tridiag_mat_to_mat(const PnlTridiagMat * mat);
extern PnlTridiagMat* pnl_tridiag_mat_copy(const PnlTridiagMat * mat);
extern void pnl_tridiag_mat_clone(PnlTridiagMat *clone, const PnlTridiagMat * mat);
extern void pnl_tridiag_mat_free(PnlTridiagMat **m);
extern int pnl_tridiag_mat_resize(PnlTridiagMat *v, int size);

extern void pnl_tridiag_mat_map_inplace(PnlTridiagMat *lhs, double(*f)(double));
extern void pnl_tridiag_mat_map_tridiag_mat_inplace(PnlTridiagMat *lhs, const PnlTridiagMat *rhs, double(*f)(double,double));
extern void pnl_tridiag_mat_print(const PnlTridiagMat *M);
extern void pnl_tridiag_mat_fprint(FILE *fic, const PnlTridiagMat *M);

extern void pnl_tridiag_mat_plus_tridiag_mat(PnlTridiagMat *lhs, const PnlTridiagMat *rhs); 
extern void pnl_tridiag_mat_minus_tridiag_mat(PnlTridiagMat *lhs, const PnlTridiagMat *rhs); 
extern void pnl_tridiag_mat_plus_scalar(PnlTridiagMat *lhs, double x);
extern void pnl_tridiag_mat_minus_scalar(PnlTridiagMat *lhs, double x);
extern void pnl_tridiag_mat_mult_scalar(PnlTridiagMat *lhs, double x); 
extern void pnl_tridiag_mat_div_scalar(PnlTridiagMat *lhs, double x); 
extern void pnl_tridiag_mat_mult_tridiag_mat_term(PnlTridiagMat *lhs, const PnlTridiagMat *rhs);
extern void pnl_tridiag_mat_div_tridiag_mat_term(PnlTridiagMat *lhs, const PnlTridiagMat *rhs);

extern PnlVect* pnl_tridiag_mat_mult_vect(const PnlTridiagMat *mat,const PnlVect *vec);
extern void pnl_tridiag_mat_mult_vect_inplace(PnlVect *lhs, const PnlTridiagMat *mat, const PnlVect *rhs);
extern void pnl_tridiag_mat_lAxpby(double l, const PnlTridiagMat *A, const PnlVect *x, double b, PnlVect * y);
extern double pnl_tridiag_mat_scalar_prod(const PnlTridiagMat *mat, const PnlVect *lhs, const PnlVect *rhs);
extern int pnl_tridiag_mat_syslin (PnlVect *lhs, PnlTridiagMat *M,const PnlVect *rhs);
extern int pnl_tridiag_mat_syslin_inplace (PnlTridiagMat *M, PnlVect *rhs);


extern void pnl_tridiag_mat_set(PnlTridiagMat *v, int d, int up, double x);
extern double pnl_tridiag_mat_get(const PnlTridiagMat *self, int d, int up);
extern double* pnl_tridiag_mat_lget(PnlTridiagMat *self, int d, int up);

typedef struct _PnlTridiagMatLUObject PnlTridiagMatLUObject;
typedef struct _PnlTridiagMatLU PnlTridiagMatLU;

struct _PnlTridiagMatLU
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows any PnlTridiagMatLU pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int size; /*!< number of rows, the matrix must be square */
  double *D; /*!< diagonal elements */
  double *DU; /*!< upper diagonal elements */
  double *DU2; /*!< second upper diagonal elements */
  double *DL; /*!< lower diagonal elements */
  int *ipiv; /*!< Permutation: row i has been interchanged with row ipiv(i) */
};

struct _PnlTridiagMatLUObject
{
  /** 
   * Must be the first element in order for the object mechanism to work
   * properly. This allows any PnlTridiagMatLU pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int size; /*!< number of rows, the matrix must be square */
  void *D; /*!< diagonal elements */
  void *DU; /*!< upper diagonal elements */
  void *DU2; /*!< second upper diagonal elements */
  void *DL; /*!< lower diagonal elements */
  int *ipiv; /*!< Permutation: row i has been interchanged with row ipiv(i) */
};


extern PnlTridiagMatLUObject* pnl_tridiag_mat_lu_object_new();
extern void pnl_tridiag_mat_lu_object_free (PnlTridiagMatLUObject **BM);
extern PnlTridiagMatLU* pnl_tridiag_mat_lu_new();
extern PnlTridiagMatLU* pnl_tridiag_mat_lu_create(int size);
extern PnlTridiagMatLU* pnl_tridiag_mat_lu_copy(const PnlTridiagMatLU  *mat);
extern void pnl_tridiag_mat_lu_clone(PnlTridiagMatLU *clone, const PnlTridiagMatLU *mat);
extern void pnl_tridiag_mat_lu_free(PnlTridiagMatLU **m);
extern int pnl_tridiag_mat_lu_resize(PnlTridiagMatLU *v, int size);
extern int pnl_tridiag_mat_lu_compute (PnlTridiagMatLU *LU, const PnlTridiagMat *A);
extern int pnl_tridiag_mat_lu_syslin_inplace (PnlTridiagMatLU *LU, PnlVect *b);
extern int pnl_tridiag_mat_lu_syslin (PnlVect *x, PnlTridiagMatLU *LU, const PnlVect *b);

/*@}*/

/*
 * Somes deprecated names
 */
#include "pnl/pnl_deprecated.h"


#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_TRIDIAG_MATRIX_H  */
