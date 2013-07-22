#ifndef _PNL_ARRAY_H_
#define _PNL_ARRAY_H_

#include "pnl/pnl_object.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**
 * \defgroup PnlArray an array object
 * This object is used to store doubly linked arrays.
 *
 * Note that elements are not hard copied into arrays. Arrays only store
 * addresses of #PnlObject. So to easily manage the memory used by the
 * elements inserted into lists, we use the PnlObject#nref field to make
 * sure there is no more reference to an object before deleting it.
 */
/*@{*/

typedef struct _PnlArray PnlArray;
struct _PnlArray
{
  /**
   * Must be the first element in order for the object mechanism to work
   * properly. This allows any PnlArray pointer to be cast to a PnlObject
   */
  PnlObject object; 
  int size;
  PnlObject **array;
  int mem_size;
};

extern PnlArray* pnl_array_new ();
extern PnlArray* pnl_array_create (int n);
extern PnlArray* pnl_array_copy (const PnlArray*);
extern void pnl_array_clone (PnlArray *, const PnlArray*);
extern int pnl_array_resize(PnlArray * v, int size);
extern void pnl_array_free (PnlArray **T);
extern PnlObject* pnl_array_get (const PnlArray *T, int i);
extern void pnl_array_set (PnlArray *T, int i, PnlObject *O);
extern void pnl_array_print (PnlArray *T);

/*@}*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _PNL_ARRAY_H_ */
