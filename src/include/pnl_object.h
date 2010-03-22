#ifndef _PNL_OBJECT_H
#define _PNL_OBJECT_H

/**
 * \defgroup PnlObject The top level object used for derivation 
 *
 * All objects of the Pnl must inherit from \a PnlObject and implement a
 * \a pnl_xxx_new function. For instance, a \a PnlFooBar object must implement a
 * function
 * \a pnl_foo_bar_new() which creates an empty \a PnlFooBar object with all its
 * elements properly initialized.
 *
 * Each type type must be added into the function \a pnl_object_new
 */
typedef struct _PnlObject PnlObject;

#define PNL_OBJECT(o) ((PnlObject *) o)
#define PNL_GET_TYPE (o) ( ((PnlObject *) o)->type)
#define PNL_GET_PARENT_TYPE (o) ( ((PnlObject *) o)->parent_type)

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
  PNL_TYPE_BAND_MATRIX,
  PNL_TYPE_BAND_MATRIX_DOUBLE,
  PNL_TYPE_HMATRIX,
  PNL_TYPE_HMATRIX_DOUBLE,
  PNL_TYPE_HMATRIX_INT,
  PNL_TYPE_HMATRIX_COMPLEX,
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

/**
 * The \a PnlObject structure is used to simulate some inheritance between the
 * ojbects of Pnl.  It must be the first element of all the objects existing in
 * Pnl so that casting any object to a PnlObject is legal 
 * 
 */
struct _PnlObject
{
  PnlType type; /*!< a unique integer id */
  const char *label; /*!< a string identifier (for the moment not useful) */
  PnlType parent_type;
};

extern PnlObject* pnl_object_new (PnlType type);


#endif /* _PNL_OBJECT_H */
