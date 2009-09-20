#ifndef VECTOR_UINT_H
#define VECTOR_UINT_H

#include <stdio.h>
#include "pnl_vector.h"

typedef unsigned int uint;
extern uint log2uint(uint x);
extern double pnl_dyadic_cast(uint i);

/**
 * \defgroup PnlVectUint Uint Vector structure for Premia
 */
/*@{*/
typedef struct PnlVectUint{
  int size;/*!< size of the vector */ 
  int mem_size; /*!< size of the memory block allocated for array */
  uint *array;/*!< pointer to store the data */
  int owner; /*!< 1 if the structure owns its array pointer */
} PnlVectUint;


#ifdef HAVE_INLINE 
extern inline
uint pnl_vect_uint_get (const PnlVectUint *self, int i)
{
  CheckIndexVect(self,i);
  return self->array[i];
}

extern inline
uint* pnl_vect_uint_lget (PnlVectUint *self, int i)
{
  CheckIndexVect(self,i);
  return &(self->array[i]);
}

extern inline
void pnl_vect_uint_set (PnlVectUint *self, int i, uint x)
{
  CheckIndexVect(self,i);
  self->array[i]=x;
}
#endif

extern void pnl_vect_uint_set(PnlVectUint *v, int i, uint x);
extern uint pnl_vect_uint_get(const PnlVectUint *v, int i);
extern uint* pnl_vect_uint_lget(PnlVectUint *v, int i);
extern void pnl_vect_uint_free(PnlVectUint **v);
extern PnlVectUint* pnl_vect_uint_create(int size);
extern PnlVectUint pnl_vect_uint_create_wrap_array(const uint *x, int size);
extern PnlVectUint* pnl_vect_uint_create_from_uint(int size, uint x);
extern PnlVectUint* pnl_vect_uint_create_from_ptr(int size, const uint* x);
extern PnlVectUint* pnl_vect_uint_create_from_list(int size,...);
extern PnlVectUint* pnl_vect_uint_create_from_file (const char * file);
extern int pnl_vect_uint_resize(PnlVectUint *v, int size);
extern int pnl_vect_uint_resize_from_uint(PnlVectUint *v, int size, uint x);
extern int pnl_vect_uint_resize_from_ptr(PnlVectUint *v, int size, const uint *t);
extern PnlVectUint* pnl_vect_uint_copy(const PnlVectUint *v);
extern void pnl_vect_uint_clone(PnlVectUint *clone, const PnlVectUint *v);
extern PnlVectUint pnl_vect_uint_wrap_subvect(const PnlVectUint *V, int i,int s);
extern PnlVectUint pnl_vect_uint_wrap_subvect_with_last(const PnlVectUint *V, int i,int j);

extern void pnl_vect_uint_print(const PnlVectUint *V);
extern void pnl_vect_uint_print_nsp(const PnlVectUint *V);
extern void pnl_vect_uint_fprint(FILE *fic, const PnlVectUint *V);
extern void pnl_vect_uint_plus_vect(PnlVectUint *lhs, const PnlVectUint *rhs); 
/*lhs+=rhs*/
extern void pnl_vect_uint_map_inplace(PnlVectUint *lhs, uint(*f)(uint)); /*lhs=f(lhs)*/
extern void pnl_vect_uint_plus_uint(PnlVectUint *lhs, uint x); /*lhs+=x*/
extern void pnl_vect_uint_minus_uint(PnlVectUint *lhs, uint x); /*lhs-=x*/
extern void pnl_vect_uint_axpby(uint a, const PnlVectUint *x, uint b, PnlVectUint *y); /* y:=a x + b y */
extern void pnl_vect_uint_mult_uint(PnlVectUint *lhs, uint x); /*lhs*=x*/
extern void pnl_vect_uint_div_uint(PnlVectUint *lhs, uint x); /*lhs*=x*/
extern void pnl_vect_uint_inv_term(PnlVectUint *lhs); /* lhs = 1 ./ lhs*/
extern void
pnl_vect_uint_div_vect_term(PnlVectUint *lhs, const PnlVectUint *rhs);/* lhs = lhs ./ rhs*/
extern void
pnl_vect_uint_mult_vect_term(PnlVectUint *lhs, const PnlVectUint *rhs); /* lhs= lhs.*rhs */ 
extern void pnl_vect_uint_set_uint(PnlVectUint *v, uint x);/* v[j]= x */
extern void pnl_vect_uint_set_zero(PnlVectUint *v); /* v[j]= 0 */
extern uint pnl_vect_uint_sum(const PnlVectUint *lhs);/* sum(x) */
extern void pnl_vect_uint_cumsum(PnlVectUint *lhs);
extern void
pnl_vect_uint_map(PnlVectUint *lhs, const PnlVectUint *rhs,
                uint(*f)(uint));/* lhs(i)=f(rhs(i)) */
extern uint
pnl_vect_uint_scalar_prod(const PnlVectUint *rhs1, const PnlVectUint *rhs2); /*rhs1.rhs2*/
extern uint pnl_vect_uint_prod(const PnlVectUint *V); /*res=prod(V(i))*/
extern void pnl_vect_uint_cumprod(PnlVectUint *V); /*res=prod(V(i))*/
extern uint pnl_vect_uint_max(const PnlVectUint *V); /*res=max(V)*/
extern uint pnl_vect_uint_min(const PnlVectUint *V); /*res=min(V)*/
extern void pnl_vect_uint_minmax (const PnlVectUint *, uint *, uint *);
extern void pnl_vect_uint_min_index (const PnlVectUint *, uint *, int *);
extern void pnl_vect_uint_max_index (const PnlVectUint *, uint *, int *);
extern void pnl_vect_uint_minmax_index (const PnlVectUint *, uint *, uint *, int *, int *);
extern void pnl_vect_uint_qsort (PnlVectUint *, char);
extern void pnl_vect_uint_qsort_index (PnlVectUint *, PnlVectInt *, char);

extern double pnl_vect_uint_norm_two(const PnlVectUint *V); /*res=\Vert V \Vert_{l^2} */
extern double pnl_vect_uint_norm_one(const PnlVectUint *V); /*res=\Vert V \Vert_{l^1} */
extern double pnl_vect_uint_norm_infty(const PnlVectUint *V); /*res=\Vert V \Vert_{l^\infty} */
extern double pnl_vect_uint_norm_x(const PnlVectUint *V,double(*f)(uint)); /*res=\Vert V \Vert_{l^X} */

extern uint pnl_vect_uint_level(PnlVectUint *v, int i);
extern uint pnl_vect_uint_level_norm_one(PnlVectUint *v);
extern uint pnl_vect_uint_level_norm_inf(PnlVectUint *v);
extern void pnl_vect_uint_dyadic_cast(const PnlVectUint * v_int,PnlVect * v_out);
extern void pnl_vect_uint_swap_elements(PnlVectUint * v, int i, int j); 
extern void pnl_vect_unit_reverse(PnlVectUint * v);


extern int pnl_vect_uint_less(const PnlVectUint * a,const PnlVectUint * b);
extern int pnl_vect_uint_equal(const PnlVectUint * a,const PnlVectUint * b);

/*extern void pnl_vect_uint_compute_father(const PnlVectUint *v, PnlVectUint * father, int i); */
/*extern PnlVectUint * pnl_vect_uint_create_son(const PnlVectUint *v, int i,boolean LorR); */
/*extern void pnl_vect_uint_compute_son(const PnlVectUint *v,PnlVectUint *son, int i,boolean LorR); */
/*@}*/



#endif /* VECTOR_UINT_H */


