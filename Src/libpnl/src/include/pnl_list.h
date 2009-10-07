#ifndef _PNL_LIST_H
#define _PNL_LIST_H


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdlib.h>
#include <stdio.h>
#include "pnl_vector_uint.h"

typedef struct PnlContains{
  int index;
  double value;
}PnlContains;


extern PnlContains * pnl_contains_create(const int ind, double Val);
extern PnlContains * pnl_contains_clone(int ind, double Val);
extern void pnl_contains_fprint(FILE *fic,PnlContains *C);
extern void  pnl_contains_add(PnlContains *C,const PnlContains *C2);
extern int  pnl_contains_less(const PnlContains *C1,const PnlContains *C2);
extern int  pnl_contains_equal(const PnlContains *C1,const PnlContains *C2);
extern PnlContains * pnl_contains_copy(const PnlContains *C2);
extern void pnl_contains_free(PnlContains **C);

typedef struct _PnlNode PnlNode;

struct _PnlNode{
  PnlNode  *previous;
  PnlNode  *next;
  PnlContains *obj;
};

typedef struct PnlSortList {
  int size; /*!< size of the List */
  PnlNode * first;
  PnlNode * last;
  PnlNode * current;
} PnlSortList;

extern PnlSortList * pnl_sort_list_create();
extern void pnl_sort_list_free(PnlSortList ** List);
extern int pnl_sort_list_find(PnlSortList * List,PnlNode **current,int Key, double Val);
extern int pnl_sort_list_find_dicho(PnlSortList * List,PnlNode **current,int Key, double Val);
extern void pnl_sort_list_add(PnlSortList * List,const PnlContains *Val);
extern void pnl_sort_list_add_dicho(PnlSortList * List,const PnlContains *Val);
extern void pnl_sort_list_print(const PnlSortList * List);



typedef struct PnlSparsePoint{
  PnlVectUint *index;
  int value;
}PnlSparsePoint;

extern PnlSparsePoint *pnl_sparse_point_create(const PnlVectUint *ind, int Val);
extern PnlSparsePoint * pnl_sparse_point_clone(PnlVectUint * ind,int val);
extern void pnl_sparse_point_fprint(FILE *fic,PnlSparsePoint *C);
extern void  pnl_sparse_point_add(PnlSparsePoint *C,const PnlSparsePoint *C2);
extern int  pnl_sparse_point_less(const PnlSparsePoint *C1,const PnlSparsePoint *C2);
extern int  pnl_sparse_point_equal(const PnlSparsePoint *C1,const PnlSparsePoint *C2);
extern PnlSparsePoint * pnl_sparse_point_copy(const PnlSparsePoint *C2);
extern void pnl_sparse_point_free(PnlSparsePoint **C);


typedef struct _PnlNodeSparsePoint PnlNodeSparsePoint;

struct _PnlNodeSparsePoint{
  PnlNodeSparsePoint  *previous;
  PnlNodeSparsePoint  *next;
  PnlSparsePoint *obj;
};

extern  void pnl_node_sparse_point_free(PnlNodeSparsePoint **N);
 
typedef struct PnlSortListSparsePoint{
  int size; //!< size of the List 
  PnlNodeSparsePoint * first;
  PnlNodeSparsePoint * last;
  PnlNodeSparsePoint * current;
} PnlSortListSparsePoint;

extern PnlSortListSparsePoint* pnl_sort_list_sparse_point_create();
extern void pnl_sort_list_sparse_point_free(PnlSortListSparsePoint ** List);
extern int pnl_sort_list_sparse_point_find(PnlSortListSparsePoint * List,PnlNodeSparsePoint **current,PnlVectUint *Key, int Val);
extern int pnl_sort_list_sparse_point_find_dicho(PnlSortListSparsePoint * List,PnlNodeSparsePoint **current,PnlVectUint *Key, int Val);
extern void pnl_sort_list_sparse_point_add(PnlSortListSparsePoint * List,const PnlSparsePoint *Val);
extern void pnl_sort_list_sparse_point_add_dicho(PnlSortListSparsePoint * List,const PnlSparsePoint *Val);
extern void pnl_sort_list_sparse_point_print(const PnlSortListSparsePoint * List);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_LIST_H */
