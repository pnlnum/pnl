#ifndef _PNL_LIST_H_
#define _PNL_LIST_H_

#include "pnl/pnl_object.h"

/**
 * \defgroup PnlList a List object
 * This object is used to store doubly linked lists.
 */
/*@{*/
typedef struct _PnlCell PnlCell;
struct _PnlCell
{
  struct _PnlCell *prev;  /*!< previous cell or 0 */
  struct _PnlCell *next;  /*!< next cell or 0 */
  PnlObject *self;       /*!< stored object */
};


typedef struct _PnlList PnlList;
struct _PnlList
{
  /**
   * Must be the first element in order for the object mechanism to work
   * properly. This allows any PnlList pointer to be cast to a PnlObject
   */
  PnlObject object; 
  PnlCell *first; /*!< first element of the list */
  PnlCell *last; /*!< last element of the list */
  int len; /*!< length of the list */
};

extern PnlList* pnl_list_new ();
extern PnlCell* pnl_cell_new ();
extern void pnl_list_free (PnlList **L);
extern void pnl_cell_free (PnlCell **c);
extern void pnl_list_insert_first (PnlList *L, PnlObject *o);
extern void pnl_list_insert_last (PnlList *L, PnlObject *o);
extern void pnl_list_remove_last (PnlList *L);
extern void pnl_list_remove_first (PnlList *L);
extern void pnl_list_remove_i (PnlList *L, int i);
extern void pnl_list_concat (PnlList *L1, PnlList *L2);
extern void pnl_list_print (const PnlList *L);

/*@}*/

#endif /* _PNL_LIST_H_ */