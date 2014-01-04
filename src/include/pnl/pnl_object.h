#ifndef _PNL_OBJECT_H
#define _PNL_OBJECT_H


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**
 * \defgroup PnlObject The top level object
 *
 *
 * \brief PnlOject is the top level object used for inheritance
 *
 * All objects of the Pnl must inherit from \a PnlObject and implement a
 * \a pnl_xxx_new function. For instance, a \a PnlFooBar object must implement a
 * function
 * \a pnl_foo_bar_new() which creates an empty \a PnlFooBar object with all its
 * elements properly initialized.
 *
 * Each type type must be added into the function \a pnl_object_create
 */
/*@{*/
typedef struct _PnlObject PnlObject;

/**
 * Casts any object into a PnlObject
 */
#define PNL_OBJECT(o) ((PnlObject *) (o))
/**
 * Casts any object into a PnlVectObject
 */
#define PNL_VECT_OBJECT(o) ((PnlVectObject *) (o))
/**
 * Casts any object into a PnlMatObject
 */
#define PNL_MAT_OBJECT(o) ((PnlMatObject *) (o))
/**
 * Casts any object into a PnlHmatObject
 */
#define PNL_HMAT_OBJECT(o) ((PnlHmatObject *) (o))
/**
 * Casts any object into a PnlBandMatObject
 */
#define PNL_BAND_MAT_OBJECT(o) ((PnlBandMatObject *) (o))
/**
 * Casts any object into a PnlTridiagMatObject
 */
#define PNL_TRIDIAGMAT_OBJECT(o) ((PnlTridiagMatObject *) (o))
/**
 * Casts any object into a PnlTridiagMatLUObject
 */
#define PNL_TRIDIAGMATLU_OBJECT(o) ((PnlTridiagMatLUObject *) (o))
/**
 * Casts any object into a PnlSpMatObject
 */
#define PNL_SP_MAT_OBJECT(o) ((PnlSpMatObject *) (o))
/**
 * Casts any object into a PnlBasis
 */
#define PNL_BASIS_OBJECT(o) ((PnlBasis *) (o))
/**
 * Casts any object into a PnlRng
 */
#define PNL_RNG_OBJECT(o) ((PnlRng *) (o))
/**
 * Casts any object into a PnlList
 */
#define PNL_LIST_OBJECT(o) ((PnlList *) (o))
/**
 * Casts any object into a PnlArray
 */
#define PNL_ARRAY_OBJECT(o) ((PnlArray *) (o))

/**
 * Returns the name of the type of any object inheriting from PnlObject
 */
#define PNL_GET_TYPENAME(o) (((PnlObject *) (o))->label)
/**
 * Returns the type of any object inheriting from PnlObject
 */
#define PNL_GET_TYPE(o) (((PnlObject *) (o))->type)
/**
 * Returns the parent type of any object inheriting from PnlObject
 */
#define PNL_GET_PARENT_TYPE(o) (((PnlObject *) (o))->parent_type)

/**
 * PnlType is used to store the id of all the objects existing in Pnl 
 */
enum {
  PNL_TYPE_NONE,
  PNL_TYPE_OBJECT,
  PNL_TYPE_VECTOR,
  PNL_TYPE_VECTOR_COMPACT,
  PNL_TYPE_VECTOR_DOUBLE,
  PNL_TYPE_VECTOR_INT,
  PNL_TYPE_VECTOR_COMPLEX,
  PNL_TYPE_MATRIX,
  PNL_TYPE_MATRIX_DOUBLE,
  PNL_TYPE_MATRIX_INT,
  PNL_TYPE_MATRIX_COMPLEX,
  PNL_TYPE_TRIDIAG_MATRIX,
  PNL_TYPE_TRIDIAG_MATRIX_DOUBLE,
  PNL_TYPE_TRIDIAG_MATRIX_LU,
  PNL_TYPE_TRIDIAG_MATRIX_LU_DOUBLE,
  PNL_TYPE_BAND_MATRIX,
  PNL_TYPE_BAND_MATRIX_DOUBLE,
  PNL_TYPE_HMATRIX,
  PNL_TYPE_HMATRIX_DOUBLE,
  PNL_TYPE_HMATRIX_INT,
  PNL_TYPE_HMATRIX_COMPLEX,
  PNL_TYPE_SP_MATRIX,
  PNL_TYPE_SP_MATRIX_DOUBLE,
  PNL_TYPE_SP_MATRIX_INT,
  PNL_TYPE_SP_MATRIX_COMPLEX,
  PNL_TYPE_LIST,
  PNL_TYPE_ARRAY,
  PNL_TYPE_BASIS,
  PNL_TYPE_RNG,
  PNL_TYPE_ITERATION_BASE,
  PNL_TYPE_CG_SOLVER,
  PNL_TYPE_BICG_SOLVER, 
  PNL_TYPE_GMRES_SOLVER 
};

/**
 * Because we are planning to write an MPI interface for the Pnl, we do not want
 * to work with enumerations but we would rather consider types which natively
 * exist in MPI.
 */
typedef unsigned int PnlType; 

typedef void (DestroyFunc) (void **);
typedef PnlObject* (CopyFunc) (PnlObject *);
typedef PnlObject* (NewFunc) (PnlObject *);
typedef void (CloneFunc) (PnlObject *dest, const PnlObject *src);
/**
 * The \a PnlObject structure is used to simulate some inheritance between the
 * ojbects of Pnl.  It must be the first element of all the objects existing in
 * Pnl so that casting any object to a PnlObject is legal.
 *
 * To improve the memory management of objects inserted into arrays and
 * lists, we keep track of the number of references on a given object (see
 * nref) to easily free objects which could be inserted in everal lists. 
 * 
 */
struct _PnlObject
{
  PnlType      type; /*!< a unique integer id */
  const char  *label; /*!< a string identifier (for the moment not useful) */
  PnlType      parent_type; /*!< the identifier of the parent object is any,
                                 otherwise parent_type=id */
  int          nref; /*!< number of references to the current object */
  DestroyFunc *destroy; /*!< frees an object */
  NewFunc     *constructor; /*!< New function */
  CopyFunc    *copy; /*!< Copy function */
  CloneFunc   *clone; /*!< Clone function */
};

extern PnlObject* pnl_object_create (PnlType type);

extern void pnl_message_on ();
extern void pnl_message_off ();
extern int pnl_message_is_on ();


/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _PNL_OBJECT_H */
