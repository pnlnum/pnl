/* 
 * This file is based on a C library                              
 *  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/DC/dc.html    
 * version 0.6.1 by                                               
 *                                                                
 * Copyright (C) 2001-2009 Makoto Matsumoto and Takuji Nishimura. 
 * Copyright (C) 2009 Mutsuo Saito                                
 * This library is free software; you can redistribute it and/or  
 * modify it under the terms of the GNU Library General Public    
 * License as published by the Free Software Foundation; either   
 * version 2 of the License, or (at your option) any later        
 * version.                                                       
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.           
 * See the GNU Library General Public License for more details.   
 * You should have received a copy of the GNU Library General     
 * Public License along with this library; if not, write to the   
 * Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  
 * 02111-1307  USA                                                
 */

/*
 * The orignal code has been modified for inclusion into PNL
 *
 * Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as 
 * published by  the Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License  along with this program.  If not, see
 * <http://www.gnu.org/licenses/>. 
 *
 * Modifications : 
 *    1. The original implementation of Matsumoto, Nishimura and Saito has
 *    been modified to work with unsigned long instead of uint32_t integers
 *    because uint32_t is not a portable type. A mask has been added to
 *    handle 64-bit machines
 *
 *    2. We needed to ensure that a dcmt_state contains only contiguous
 *    data so we had to fix a few parameters to get rid of any dynamoc
 *    allocation in alloc_dcmt_state. Words are 32-bit long and the period
 *    has been fixed to 2^521-1. The search of new MTs increases very
 *    quickly with the Mersenne exponent.
 *    
 */

#include <string.h>
#include "pnl/pnl_config.h"
#include "pnl/pnl_random.h"
#include "pnl/pnl_internals.h"

#ifdef PNL_HAVE_INLINE
#define __INLINE__ inline
#else
#define __INLINE__
#endif

#define BIT32_MASK 0xffffffffUL
#define MERSENNE_P 521
#define WORDLEN 32
#define LSB 0x1
#define MAX_SEARCH 10000

#define SSS 7
#define TTT 15
#define S00 12
#define S01 18

/* for get_tempering_parameter_hard */
#define LIMIT_V_BEST_OPT 15

#define LIMIT_IRRED_DEG 31
#define NIRREDPOLY 127
#define MAX_IRRED_DEG 9

#define WORD_LEN 32
#define MIN_INFINITE (-2147483647-1)

#define NOT_REJECTED 1
#define REJECTED 0
#define REDU 0
#define IRRED 1
#define NONREDU 1


/*
 * Some name changes for integration in PNL
 * #define free_mt_struct pnl_dcmt_free
 * #define free_mt_struct_array pnl_dcmt_free_array
 * #define alloc_mt_struct alloc_dcmt_state
 */

typedef struct {int *x; int deg;} Polynomial;

typedef struct PRESCR_T {
  int sizeofA; /* parameter size */
  ulong **modlist;
  Polynomial **preModPolys;
} prescr_t;

typedef struct CHECK32_T {
  ulong upper_mask;
  ulong lower_mask;
  ulong word_mask;
} check32_t;

typedef struct EQDEG_T {
  ulong bitmask[32];
  ulong mask_b;
  ulong mask_c;
  ulong upper_v_bits;
  int shift_0;
  int shift_1;
  int shift_s;
  int shift_t;
  int mmm;
  int nnn;
  int rrr;
  int www;
  ulong aaa[2];
  ulong gupper_mask;   /** most significant  (WWW - RRR) bits **/
  ulong glower_mask;	/** least significant RRR bits **/
  ulong greal_mask;	/** upper WWW bitmask **/
  int ggap; /** difference between machine wordsize and dest wordsize **/
  int gcur_maxlengs[32];	/** for optimize_v_hard **/
  ulong gmax_b, gmax_c;
} eqdeg_t;

typedef struct {
  ulong *cf;  /* fraction part */              // status
  int start;     /* beginning of fraction part */ // idx
  int count;	   /* maximum (degree) */
  ulong next; /* (bp) rm (shifted&bitmasked) at the maximum degree */
} Vector;

typedef struct mask_node{
  ulong b,c;
  int v,leng;
  struct mask_node *next;
} MaskNode;

static __INLINE__ ulong trnstmp(eqdeg_t *eq, ulong tmp) {
  tmp &= BIT32_MASK;
  tmp ^= (tmp >> eq->shift_0) & eq->greal_mask;
  return tmp;
}

static __INLINE__ ulong masktmp(eqdeg_t *eq, ulong tmp) {
  tmp &= BIT32_MASK;
  tmp ^= (tmp << eq->shift_s) & eq->mask_b;
  tmp ^= (tmp << eq->shift_t) & eq->mask_c;
  return tmp;
}

static __INLINE__ ulong lsb(eqdeg_t *eq, ulong x) {
  return ((x & BIT32_MASK) >> eq->ggap) & 1;
}

static const unsigned char pivot_calc_tbl[256] = {
  0, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
  4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
  3, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
  4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
  2, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
  4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
  3, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
  4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
  1, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
  4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
  3, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
  4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
  2, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
  4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
  3, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
  4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
};

/* list of irreducible polynomials whose degrees are less than 10 */
static const int irredpolylist[NIRREDPOLY][MAX_IRRED_DEG+1] = {
    {0,1,0,0,0,0,0,0,0,0,},{1,1,0,0,0,0,0,0,0,0,},{1,1,1,0,0,0,0,0,0,0,},
    {1,1,0,1,0,0,0,0,0,0,},{1,0,1,1,0,0,0,0,0,0,},{1,1,0,0,1,0,0,0,0,0,},
    {1,0,0,1,1,0,0,0,0,0,},{1,1,1,1,1,0,0,0,0,0,},{1,0,1,0,0,1,0,0,0,0,},
    {1,0,0,1,0,1,0,0,0,0,},{1,1,1,1,0,1,0,0,0,0,},{1,1,1,0,1,1,0,0,0,0,},
    {1,1,0,1,1,1,0,0,0,0,},{1,0,1,1,1,1,0,0,0,0,},{1,1,0,0,0,0,1,0,0,0,},
    {1,0,0,1,0,0,1,0,0,0,},{1,1,1,0,1,0,1,0,0,0,},{1,1,0,1,1,0,1,0,0,0,},
    {1,0,0,0,0,1,1,0,0,0,},{1,1,1,0,0,1,1,0,0,0,},{1,0,1,1,0,1,1,0,0,0,},
    {1,1,0,0,1,1,1,0,0,0,},{1,0,1,0,1,1,1,0,0,0,},{1,1,0,0,0,0,0,1,0,0,},
    {1,0,0,1,0,0,0,1,0,0,},{1,1,1,1,0,0,0,1,0,0,},{1,0,0,0,1,0,0,1,0,0,},
    {1,0,1,1,1,0,0,1,0,0,},{1,1,1,0,0,1,0,1,0,0,},{1,1,0,1,0,1,0,1,0,0,},
    {1,0,0,1,1,1,0,1,0,0,},{1,1,1,1,1,1,0,1,0,0,},{1,0,0,0,0,0,1,1,0,0,},
    {1,1,0,1,0,0,1,1,0,0,},{1,1,0,0,1,0,1,1,0,0,},{1,0,1,0,1,0,1,1,0,0,},
    {1,0,1,0,0,1,1,1,0,0,},{1,1,1,1,0,1,1,1,0,0,},{1,0,0,0,1,1,1,1,0,0,},
    {1,1,1,0,1,1,1,1,0,0,},{1,0,1,1,1,1,1,1,0,0,},{1,1,0,1,1,0,0,0,1,0,},
    {1,0,1,1,1,0,0,0,1,0,},{1,1,0,1,0,1,0,0,1,0,},{1,0,1,1,0,1,0,0,1,0,},
    {1,0,0,1,1,1,0,0,1,0,},{1,1,1,1,1,1,0,0,1,0,},{1,0,1,1,0,0,1,0,1,0,},
    {1,1,1,1,1,0,1,0,1,0,},{1,1,0,0,0,1,1,0,1,0,},{1,0,1,0,0,1,1,0,1,0,},
    {1,0,0,1,0,1,1,0,1,0,},{1,0,0,0,1,1,1,0,1,0,},{1,1,1,0,1,1,1,0,1,0,},
    {1,1,0,1,1,1,1,0,1,0,},{1,1,1,0,0,0,0,1,1,0,},{1,1,0,1,0,0,0,1,1,0,},
    {1,0,1,1,0,0,0,1,1,0,},{1,1,1,1,1,0,0,1,1,0,},{1,1,0,0,0,1,0,1,1,0,},
    {1,0,0,1,0,1,0,1,1,0,},{1,0,0,0,1,1,0,1,1,0,},{1,0,1,1,1,1,0,1,1,0,},
    {1,1,0,0,0,0,1,1,1,0,},{1,1,1,1,0,0,1,1,1,0,},{1,1,1,0,1,0,1,1,1,0,},
    {1,0,1,1,1,0,1,1,1,0,},{1,1,1,0,0,1,1,1,1,0,},{1,1,0,0,1,1,1,1,1,0,},
    {1,0,1,0,1,1,1,1,1,0,},{1,0,0,1,1,1,1,1,1,0,},{1,1,0,0,0,0,0,0,0,1,},
    {1,0,0,0,1,0,0,0,0,1,},{1,1,1,0,1,0,0,0,0,1,},{1,1,0,1,1,0,0,0,0,1,},
    {1,0,0,0,0,1,0,0,0,1,},{1,0,1,1,0,1,0,0,0,1,},{1,1,0,0,1,1,0,0,0,1,},
    {1,1,0,1,0,0,1,0,0,1,},{1,0,0,1,1,0,1,0,0,1,},{1,1,1,1,1,0,1,0,0,1,},
    {1,0,1,0,0,1,1,0,0,1,},{1,0,0,1,0,1,1,0,0,1,},{1,1,1,1,0,1,1,0,0,1,},
    {1,1,1,0,1,1,1,0,0,1,},{1,0,1,1,1,1,1,0,0,1,},{1,1,1,0,0,0,0,1,0,1,},
    {1,0,1,0,1,0,0,1,0,1,},{1,0,0,1,1,0,0,1,0,1,},{1,1,0,0,0,1,0,1,0,1,},
    {1,0,1,0,0,1,0,1,0,1,},{1,1,1,1,0,1,0,1,0,1,},{1,1,1,0,1,1,0,1,0,1,},
    {1,0,1,1,1,1,0,1,0,1,},{1,1,1,1,0,0,1,1,0,1,},{1,0,0,0,1,0,1,1,0,1,},
    {1,1,0,1,1,0,1,1,0,1,},{1,0,1,0,1,1,1,1,0,1,},{1,0,0,1,1,1,1,1,0,1,},
    {1,0,0,0,0,0,0,0,1,1,},{1,1,0,0,1,0,0,0,1,1,},{1,0,1,0,1,0,0,0,1,1,},
    {1,1,1,1,1,0,0,0,1,1,},{1,1,0,0,0,1,0,0,1,1,},{1,0,0,0,1,1,0,0,1,1,},
    {1,1,0,1,1,1,0,0,1,1,},{1,0,0,1,0,0,1,0,1,1,},{1,1,1,1,0,0,1,0,1,1,},
    {1,1,0,1,1,0,1,0,1,1,},{1,0,0,0,0,1,1,0,1,1,},{1,1,0,1,0,1,1,0,1,1,},
    {1,0,1,1,0,1,1,0,1,1,},{1,1,0,0,1,1,1,0,1,1,},{1,1,1,1,1,1,1,0,1,1,},
    {1,0,1,0,0,0,0,1,1,1,},{1,1,1,1,0,0,0,1,1,1,},{1,0,0,0,0,1,0,1,1,1,},
    {1,0,1,0,1,1,0,1,1,1,},{1,0,0,1,1,1,0,1,1,1,},{1,1,1,0,0,0,1,1,1,1,},
    {1,1,0,1,0,0,1,1,1,1,},{1,0,1,1,0,0,1,1,1,1,},{1,0,1,0,1,0,1,1,1,1,},
    {1,0,0,1,1,0,1,1,1,1,},{1,1,0,0,0,1,1,1,1,1,},{1,0,0,1,0,1,1,1,1,1,},
    {1,1,0,1,1,1,1,1,1,1,},
};

/*******************************************************************/
/* Static functions */
static void MakepreModPolys(prescr_t *pre, int mm, int nn, int rr, int ww);
static Polynomial *make_tntm( int n, int m);
static Polynomial *PolynomialDup(Polynomial *pl);
static void PolynomialMod(Polynomial *wara, const Polynomial *waru);
static Polynomial *PolynomialMult(Polynomial *p0, Polynomial *p1);
static void FreePoly( Polynomial *p);
static Polynomial *NewPoly(int degree);
static int IsReducible(prescr_t *pre, ulong aaa, ulong *polylist);
static ulong word2bit(Polynomial *pl);
static void makemodlist(prescr_t *pre, Polynomial *pl, int nPoly);
static void NextIrredPoly(Polynomial *pl, int nth);


static ulong nextA(mt_state *org, int w);
static ulong nextA_id(mt_state *org, int w, int id, int idw);
static void make_masks(int r, int w, dcmt_state *mts);
static int get_irred_param(check32_t *ck, prescr_t *pre, mt_state *org,
                           dcmt_state *mts,int id, int idw);
static void _get_tempering_parameter_hard_dc(dcmt_state *mts);
static dcmt_state *alloc_dcmt_state(int n);
static dcmt_state *init_mt_search(check32_t *ck, prescr_t *pre, int p);
static void end_mt_search(prescr_t *pre);
static void copy_params_of_dcmt_state(dcmt_state *src, dcmt_state *dst);
static int proper_mersenne_exponent(int p);
static int calc_pivot(ulong v);
static int push_stack(eqdeg_t *eq, ulong b, ulong c,
                      int v, ulong *bbb, ulong *ccc);
static int push_mask(eqdeg_t * eq, int l, int v,
                     ulong b, ulong c, ulong *bbb, ulong *ccc);
static int pivot_reduction(eqdeg_t *eq, int v);
static void init_tempering(eqdeg_t *eq, dcmt_state *mts);
static void free_Vector( Vector *v );
static void free_lattice( Vector **lattice, int v);
static void add(int nnn, Vector *u, Vector *v);
static void optimize_v(eqdeg_t *eq, ulong b, ulong c, int v);
static MaskNode *optimize_v_hard(eqdeg_t *eq, int v, MaskNode *prev);
static Vector *new_Vector(int nnn);
static Vector **make_lattice(eqdeg_t *eq, int v);
static void delete_MaskNodes(MaskNode *head);
static MaskNode *delete_lower_MaskNodes(MaskNode *head, int l);
static MaskNode *cons_MaskNode(MaskNode *head, ulong b, ulong c, int leng);
static void next_state(eqdeg_t *eq, Vector *v, int *count);
/*******************************************************************/

static int _prescreening_dc(prescr_t *pre, ulong aaa)
{

  int i;

  for (i=0; i<NIRREDPOLY; i++) {
    if (IsReducible(pre, aaa,pre->modlist[i])==REDU)
      return REJECTED;
  }
  return NOT_REJECTED;
}

static void _InitPrescreening_dc(prescr_t *pre, int m, int n, int r, int w)
{
  int i;
  Polynomial *pl;

  pre->sizeofA = w;

  pre->preModPolys = (Polynomial **)malloc(
    (pre->sizeofA+1)*(sizeof(Polynomial*)));
  if (NULL == pre->preModPolys) {
    printf ("malloc error in \"InitPrescreening\"\n");
    exit(1);
  }
  MakepreModPolys(pre, m,n,r,w);

  pre->modlist = (ulong**)malloc(NIRREDPOLY * sizeof(ulong*));
  if (NULL == pre->modlist) {
    printf ("malloc error in \"InitPrescreening()\"\n");
    exit(1);
  }
  for (i=0; i<NIRREDPOLY; i++) {
    pre->modlist[i]
      = (ulong*)malloc( (pre->sizeofA + 1) * (sizeof(ulong)) );
    if (NULL == pre->modlist[i]) {
      printf ("malloc error in \"InitPrescreening()\"\n");
      exit(1);
    }
  }


  for (i=0; i<NIRREDPOLY; i++) {
    pl = NewPoly(MAX_IRRED_DEG);
    NextIrredPoly(pl,i);
    makemodlist(pre, pl, i);
    FreePoly(pl);
  }

  for (i=pre->sizeofA; i>=0; i--)
    FreePoly(pre->preModPolys[i]);
  free(pre->preModPolys);

}

static void _EndPrescreening_dc(prescr_t *pre)
{
  int i;

  for (i=0; i<NIRREDPOLY; i++)
    free(pre->modlist[i]);
  free(pre->modlist);
}

static void NextIrredPoly(Polynomial *pl, int nth)
{
  int i, max_deg;

  for (max_deg=0,i=0; i<=MAX_IRRED_DEG; i++) {
    if ( irredpolylist[nth][i] )
      max_deg = i;
    pl->x[i] = irredpolylist[nth][i];
  }

  pl->deg = max_deg;

}

static void makemodlist(prescr_t *pre, Polynomial *pl, int nPoly)
{
  Polynomial *tmpPl;
  int i;

  for (i=0; i<=pre->sizeofA; i++) {
    tmpPl = PolynomialDup(pre->preModPolys[i]);
    PolynomialMod(tmpPl,pl);
    pre->modlist[nPoly][i] = word2bit(tmpPl);
    FreePoly(tmpPl);
  }
}

/* Pack Polynomial into a word */
static ulong word2bit(Polynomial *pl)
{
  int i;
  ulong bx;

  bx = 0;
  for (i=pl->deg; i>0; i--) {
    if (pl->x[i]) bx |= 0x1;
    bx <<= 1;
  }
  if (pl->x[0]) bx |= 0x1;

  return bx;
}

/* REDU -- reducible */
/* aaa = (a_{w-1}a_{w-2}...a_1a_0 */
static int IsReducible(prescr_t *pre, ulong aaa, ulong *polylist)
{
  int i;
  ulong x;

  x = polylist[pre->sizeofA];
  for (i=pre->sizeofA-1; i>=0; i--) {
    if (aaa&0x1)
      x ^= polylist[i];
    aaa >>= 1;
  }

  if ( x == 0 ) return REDU;
  else return NONREDU;
}


/***********************************/
/**   functions for polynomial    **/
/***********************************/
static Polynomial *NewPoly(int degree)
{
  Polynomial *p;

  p = (Polynomial *)calloc( 1, sizeof(Polynomial));
  if( p==NULL ){
    printf("calloc error in \"NewPoly()\"\n");
    exit(1);
  }
  p->deg = degree;

  if (degree < 0) {
    p->x = NULL;
    return p;
  }

  p->x = (int *)calloc( degree + 1, sizeof(int));
  if( p->x == NULL ){
    printf("calloc error\n");
    exit(1);
  }

  return p;
}

static void FreePoly( Polynomial *p)
{
  if (p->x != NULL)
    free( p->x );
  free( p );
}


/** multiplication **/
static Polynomial *PolynomialMult(Polynomial *p0,Polynomial *p1)
{
  int i, j;
  Polynomial *p;

  /* if either p0 or p1 is 0, return 0 */
  if ( (p0->deg < 0) || (p1->deg < 0) ) {
    p = NewPoly(-1);
    return p;
  }

  p = NewPoly(p0->deg + p1->deg);
  for( i=0; i<=p1->deg; i++){
    if( p1->x[i] ){
      for( j=0; j<=p0->deg; j++){
        p->x[i+j] ^= p0->x[j];
      }
    }
  }

  return p;
}

/** wara mod waru **/
/** the result is stored in wara ********/
static void PolynomialMod( Polynomial *wara, const Polynomial *waru)
{
  int i;
  int deg_diff;

  while( wara->deg >= waru->deg  ){
    deg_diff = wara->deg - waru->deg;
    for( i=0; i<=waru->deg; i++){
      wara->x[ i+deg_diff ] ^= waru->x[i];
    }

    for( i=wara->deg; i>=0; i--){
      if( wara->x[i] ) break;
    }
    wara->deg=i;

  }
}

static Polynomial *PolynomialDup(Polynomial *pl)
{
  Polynomial *pt;
  int i;

  pt = NewPoly(pl->deg);
  for (i=pl->deg; i>=0; i--)
    pt->x[i] = pl->x[i];

  return pt;
}

/** make the polynomial  "t**n + t**m"  **/
static Polynomial *make_tntm( int n, int m)
{
  Polynomial *p;

  p = NewPoly(n);
  p->x[n] = p->x[m] = 1;

  return p;
}

static void MakepreModPolys(prescr_t *pre, int mm, int nn, int rr, int ww)
{
  Polynomial *t, *t0, *t1, *s, *s0, *s1;
  int i,j;

  j = 0;
  t = NewPoly(0);
  t->deg = 0;
  t->x[0] = 1;
  pre->preModPolys[j++] = t;

  t = make_tntm (nn, mm);
  t0 = make_tntm (nn, mm);
  s = make_tntm (nn-1, mm-1);

  for( i=1; i<(ww - rr); i++){
    pre->preModPolys[j++] = PolynomialDup(t0);
    t1 = t0;
    t0 = PolynomialMult(t0, t);
    FreePoly(t1);
  }

  pre->preModPolys[j++] = PolynomialDup(t0);

  s0 =PolynomialMult( t0, s);
  FreePoly(t0);	FreePoly(t);
  for( i=(rr-2); i>=0; i--){
    pre->preModPolys[j++] = PolynomialDup(s0);
    s1 = s0;
    s0 = PolynomialMult( s0, s);
    FreePoly(s1);
  }

  pre->preModPolys[j++] = PolynomialDup(s0);

  FreePoly(s0); FreePoly(s);
}

static ulong nextA(mt_state *org, int w)
{
  ulong x, word_mask;

  word_mask = 0xFFFFFFFF;
  word_mask <<= WORDLEN - w;
  word_mask >>= WORDLEN - w;

  x = pnl_mt_genrand(org);
  x &= word_mask;
  x |= (LSB << (w-1));
  x &= BIT32_MASK;

  return x;
}

static ulong nextA_id(mt_state *org, int w, int id, int idw)
{
  ulong x, word_mask;

  word_mask = 0xFFFFFFFF;
  word_mask <<= WORDLEN - w;
  word_mask >>= WORDLEN - w;
  word_mask >>= idw;
  word_mask <<= idw;

  x = pnl_mt_genrand(org);
  x &= word_mask;
  x |= (LSB << (w-1));
  x |= (ulong)id; /* embedding id */
  x &= BIT32_MASK;

  return x;
}

static void _InitCheck32_dc(check32_t *ck, int r, int w)
{
  int i;

  /* word_mask (least significant w bits) */
  ck->word_mask = 0xFFFFFFFFUL;
  ck->word_mask <<= WORDLEN - w;
  ck->word_mask >>= WORDLEN - w;
  /* lower_mask (least significant r bits) */
  for (ck->lower_mask=0,i=0; i<r; ++i) {
    ck->lower_mask <<= 1;
    ck->lower_mask |= LSB;
  }
  /* upper_mask (most significant (w-r) bits */
  ck->upper_mask = (~ck->lower_mask) & ck->word_mask;
}

static int _CheckPeriod_dc(check32_t *ck, mt_state *st,
                    ulong a, int m, int n, int r, int w)
{
  int i, j, p, pp;
  ulong y, *x, *init, mat[2];


  p = n*w-r;
  x = (ulong*) malloc (2*p*sizeof(ulong));
  if (NULL==x) {
    printf("malloc error in \"_CheckPeriod_dc()\"\n");
    exit(1);
  }

  init = (ulong*) malloc (n*sizeof(ulong));
  if (NULL==init) {
    printf("malloc error \"_CheckPeriod_dc()\"\n");
    free(x);
    exit(1);
  }

  /* set initial values */
  for (i=0; i<n; ++i)
    x[i] = init[i] = (ck->word_mask & pnl_mt_genrand(st));
  /* it is better that LSBs of x[2] and x[3] are different */
  if ( (x[2]&LSB) == (x[3]&LSB) ) {
    x[3] ^= 1;
    init[3] ^= 1;
  }

  pp = 2*p-n;
  mat[0] = 0; mat[1] = a;
  for (j=0; j<p; ++j) {

    /* generate */
    for (i=0; i<pp; ++i){
      y = (x[i]&ck->upper_mask) | (x[i+1]&ck->lower_mask);
      x[i+n] = x[i+m] ^ ( (y>>1) ^ mat[y&LSB] );
    }

    /* pick up odd subscritpt elements */
    for (i=2; i<=p; ++i)
      x[i] = x[(i<<1)-1];

    /* reverse generate */
    for (i=p-n; i>=0; --i) {
      y = x[i+n] ^ x[i+m] ^ mat[ x[i+1]&LSB ];
      y <<=1; y |= x[i+1]&LSB;

      x[i+1] = (x[i+1]&ck->upper_mask) | (y&ck->lower_mask);
      x[i] = (y&ck->upper_mask) | (x[i]&ck->lower_mask);
    }

  }

  if ((x[0]&ck->upper_mask)==(init[0]&ck->upper_mask)) {
    for (i=1; i<n; ++i) {
      if (x[i] != init[i])
        break;
    }
    if (i==n) {
      free(x); free(init);
      return IRRED;
    }
  }


  free(x); free(init);
  return REDU;
}

/* When idw==0, id is not embedded into "a" */
#define FOUND 1
#define NOT_FOUND 0
static int get_irred_param(check32_t *ck, prescr_t *pre, mt_state *org,
                           dcmt_state *mts, int id, int idw)
{
  int i;
  ulong a;

  for (i=0; i<MAX_SEARCH; i++) {
    if (idw == 0)
      a = nextA(org, mts->ww);
    else
      a = nextA_id(org, mts->ww, id, idw);
    if (NOT_REJECTED == _prescreening_dc(pre, a) ) {
      if (IRRED
          == _CheckPeriod_dc(ck, org, a,mts->mm,mts->nn,mts->rr,mts->ww)) {
        mts->aaa = a;
        break;
      }
    }
  }

  if (MAX_SEARCH == i) return NOT_FOUND;
  return FOUND;
}



static void make_masks(int r, int w, dcmt_state *mts)
{
  int i;
  ulong ut, wm, um, lm;

  wm = 0xFFFFFFFF;
  wm >>= (WORDLEN - w);

  ut = 0;
  for (i=0; i<r; i++) {
    ut <<= 1;
    ut |= LSB;
  }

  lm = ut;
  um = (~ut) & wm;

  mts->wmask = wm;
  mts->umask = um;
  mts->lmask = lm;
}

static dcmt_state *init_mt_search(check32_t *ck, prescr_t *pre, int p)
{
  int n, m, r, w;
  dcmt_state *mts;

  w = 32; /* we fix the word size to 32 bits */

  if ( !proper_mersenne_exponent(p) ) {
    if (p<521) {
      printf ("\"p\" is too small.\n");
      return NULL;
    }
    else if (p>44497){
      printf ("\"p\" is too large.\n");
      return NULL;
    }
    else {
      printf ("\"p\" is not a Mersenne exponent.\n");
      return NULL;
    }
  }

  n = p/w + 1; /* since p is Mersenne Exponent, w never divids p */
  mts = alloc_dcmt_state(n);
  if (NULL == mts) return NULL;

  m = n/2;
  if (m < 2) m = n-1;
  r = n * w - p;

  make_masks(r, w, mts);
  _InitPrescreening_dc(pre, m, n, r, w);
  _InitCheck32_dc(ck, r, w);

  mts->mm = m;
  mts->nn = n;
  mts->rr = r;
  mts->ww = w;

  return mts;
}

static void end_mt_search(prescr_t *pre)
{
  _EndPrescreening_dc(pre);
}

/*
   w -- word size fixed to 32bits
   p -- Mersenne Exponent
   seed -- seed for original mt19937 to generate parameter.
   */
dcmt_state *get_mt_parameter_st(int p, ulong seed)
{
  dcmt_state *mts;
  prescr_t pre;
  mt_state org;
  check32_t ck;

  pnl_mt_sseed(&org, seed);
  mts = init_mt_search(&ck, &pre, p);
  if (mts == NULL) return NULL;

  if ( NOT_FOUND == get_irred_param(&ck, &pre, &org, mts,0,0) ) 
    {
    pnl_dcmt_free(&mts);
    return NULL;
  }
  _get_tempering_parameter_hard_dc(mts);
  end_mt_search(&pre);

  return mts;
}


/* n : sizeof state vector */
static dcmt_state *alloc_dcmt_state(int n)
{
  dcmt_state *mts;

  mts = (dcmt_state*)malloc(sizeof(dcmt_state));
  if (NULL == mts) return NULL;
  /* mts->state = (ulong*)malloc(n*sizeof(ulong)); */
  /* if (NULL == mts->state) { */
  /*   free(mts); */
  /*   return NULL; */
  /* } */

  return mts;
}


/** 
 * Free a dcmt_state
 * 
 * @param mts
 */
void pnl_dcmt_free(dcmt_state **mts)
{
  if ( *mts != NULL ) free(*mts);
  *mts = NULL;
}

/** 
 * Free an array of dcmt_state
 * 
 * @param mtss an array of dcmt_state
 * @param count number of elements in mtss
 */
void pnl_dcmt_free_array(dcmt_state **mtss, int count)
{
  int i;

  if (mtss == NULL) {
    return;
  }
  for (i=0; i < count; i++) {
    pnl_dcmt_free(&(mtss[i]));
  }
  free(mtss);
}

static void copy_params_of_dcmt_state(dcmt_state *src, dcmt_state *dst)
{
  dst->nn = src->nn;
  dst->mm = src->mm;
  dst->rr = src->rr;
  dst->ww = src->ww;
  dst->wmask = src->wmask;
  dst->umask = src->umask;
  dst->lmask = src->lmask;
}

static int proper_mersenne_exponent(int p)
{
  switch(p) {
  case 521:
  case 607:
  case 1279:
  case 2203:
  case 2281:
  case 3217:
  case 4253:
  case 4423:
  case 9689:
  case 9941:
  case 11213:
  case 19937:
  case 21701:
  case 23209:
  case 44497:
    return 1;
  default:
    return 0;
  }
}

static void _get_tempering_parameter_hard_dc(dcmt_state *mts)
{
  int i;
  MaskNode mn0, *cur, *next;
  eqdeg_t eq;

  init_tempering(&eq, mts);

  for (i=0; i<eq.www; i++)
    eq.gcur_maxlengs[i] = -1;

  mn0.b = mn0.c = mn0.leng = 0;
  mn0.next = NULL;

  cur = &mn0;
  for (i=0; i<LIMIT_V_BEST_OPT; i++) {
    next = optimize_v_hard(&eq, i, cur);
    if (i > 0)
      delete_MaskNodes(cur);
    cur = next;
  }
  delete_MaskNodes(cur);

  optimize_v(&eq, eq.gmax_b, eq.gmax_c,i);
  mts->shift0 = eq.shift_0;
  mts->shift1 = eq.shift_1;
  mts->shiftB = eq.shift_s;
  mts->shiftC = eq.shift_t;
  mts->maskB = eq.mask_b >> eq.ggap;
  mts->maskC = eq.mask_c >> eq.ggap;

  /* show_distrib(mts); */
}

static int calc_pivot(ulong v) 
{
  int p1, p2, p3, p4;

  v &= BIT32_MASK;
  p1 = pivot_calc_tbl[v & 0xff];
  if (p1) {
    return p1 + 24 - 1;
  }
  p2 = pivot_calc_tbl[(v >> 8) & 0xff];
  if (p2) {
    return p2 + 16 - 1;
  }
  p3 = pivot_calc_tbl[(v >> 16) & 0xff];
  if (p3) {
    return p3 + 8 - 1;
  }
  p4 = pivot_calc_tbl[(v >> 24) & 0xff];
  if (p4) {
    return p4 - 1;
  }
  return -1;
}

static int is_zero(int size, Vector *v) 
{
  if (v->cf[0] != 0) {
    return 0;
  } else {
    return (memcmp(v->cf, v->cf + 1, sizeof(ulong) * (size - 1)) == 0);
  }
}

static void init_tempering(eqdeg_t *eq, dcmt_state *mts)
{
  int i;

  eq->mmm = mts->mm;
  eq->nnn = mts->nn;
  eq->rrr = mts->rr;
  eq->www = mts->ww;
  eq->shift_0 = S00;
  eq->shift_1 = S01;
  eq->shift_s = SSS;
  eq->shift_t = TTT;
  eq->ggap = WORD_LEN - eq->www;
  /* bits are filled in mts->aaa from MSB */
  eq->aaa[0] = 0; eq->aaa[1] = (mts->aaa) << eq->ggap;


  for( i=0; i<WORD_LEN; i++)
    eq->bitmask[i] = 0x80000000UL >> i;

  for( i=0, eq->glower_mask=0; i<eq->rrr; i++)
    eq->glower_mask = (eq->glower_mask<<1)| 0x1;

  eq->gupper_mask = ~eq->glower_mask;
  eq->gupper_mask <<= eq->ggap;
  eq->glower_mask <<= eq->ggap;

  eq->greal_mask = (eq->gupper_mask | eq->glower_mask);

#if defined(DEBUG)
  printf ("n=%d m=%d r=%d w=%d\n", eq->nnn, eq->mmm, eq->rrr, eq->www);
  printf ("nw-r=%d\n", eq->nnn * eq->www - eq->rrr);
  printf ("a=%x(%x << %d)\n", eq->aaa[1],mts->aaa,eq->ggap);
  printf ("upper (w-r) bit mask = %x\n", eq->gupper_mask);
  printf ("lower r bit mask     = %x\n", eq->glower_mask);
  printf ("w bit mask           = %x\n", eq->greal_mask);
  fflush(stdout);
#endif
}

/* (v-1) bitmasks of b,c */
static MaskNode *optimize_v_hard(eqdeg_t *eq, int v, MaskNode *prev_masks)
{
  int i, ll, t;
  ulong bbb[8], ccc[8];
  MaskNode *cur_masks;

  cur_masks = NULL;

  while (prev_masks != NULL) {

    ll = push_stack(eq, prev_masks->b,prev_masks->c,v,bbb,ccc);

    for (i=0; i<ll; ++i) {
      eq->mask_b = bbb[i];
      eq->mask_c = ccc[i];
      t = pivot_reduction(eq, v+1);
      if (t >= eq->gcur_maxlengs[v]) {
        eq->gcur_maxlengs[v] = t;
        eq->gmax_b = eq->mask_b;
        eq->gmax_c = eq->mask_c;
        cur_masks = cons_MaskNode(cur_masks, eq->mask_b, eq->mask_c, t);
      }
    }
    prev_masks = prev_masks->next;
  }

  cur_masks = delete_lower_MaskNodes(cur_masks, eq->gcur_maxlengs[v]);

  return cur_masks;
}


/* (v-1) bitmasks of b,c */
static void optimize_v(eqdeg_t *eq, ulong b, ulong c, int v)
{
  int i, max_len, max_i, ll, t;
  ulong bbb[8], ccc[8];

  ll = push_stack(eq, b,c,v,bbb,ccc);

  max_len = max_i = 0;
  if (ll > 1) {
    for (i=0; i<ll; ++i) {
      eq->mask_b = bbb[i];
      eq->mask_c = ccc[i];
      t = pivot_reduction(eq, v+1);
      if (t > max_len) {
        max_len = t;
        max_i = i;
      }
    }
  }

  if ( v >= eq->www-1 ) {
    eq->mask_b = bbb[max_i];
    eq->mask_c = ccc[max_i];
    return;
  }

  optimize_v(eq, bbb[max_i], ccc[max_i], v+1);
}

static int push_stack(eqdeg_t *eq, ulong b, ulong c, int v,
                      ulong *bbb, ulong *ccc)
{
  int i, ll, ncv;
  ulong cv_buf[2];

  ll = 0;

  if( (v+eq->shift_t) < eq->www ){
    ncv = 2; cv_buf[0] = c | eq->bitmask[v]; cv_buf[1] = c;
  }
  else {
    ncv = 1; cv_buf[0] = c;
  }

  for( i=0; i<ncv; ++i)
    ll += push_mask(eq, ll, v, b, cv_buf[i], bbb, ccc);

  return ll;
}

static int push_mask(eqdeg_t *eq, int l, int v, ulong b, ulong c,
                     ulong *bbb, ulong *ccc)
{
  int i, j, k, nbv, nbvt;
  ulong bmask, bv_buf[2], bvt_buf[2];

  k = l;
  if( (eq->shift_s+v) >= eq->www ){
    nbv = 1; bv_buf[0] = 0;
  }
  else if( (v>=eq->shift_t) && (c&eq->bitmask[v-eq->shift_t] ) ){
    nbv = 1; bv_buf[0] = b&eq->bitmask[v];
  }
  else {
    nbv = 2; bv_buf[0] = eq->bitmask[v]; bv_buf[1] = 0;
  }

  if( ((v+eq->shift_t+eq->shift_s) < eq->www) && (c&eq->bitmask[v]) ){
    nbvt = 2; bvt_buf[0] = eq->bitmask[v+eq->shift_t]; bvt_buf[1] = 0;
  }
  else {
    nbvt = 1; bvt_buf[0] = 0;
  }

  bmask = eq->bitmask[v];
  if( (v+eq->shift_t) < eq->www )
    bmask |= eq->bitmask[v+eq->shift_t];
  bmask = ~bmask;
  for( i=0; i<nbvt; ++i){
    for( j=0; j<nbv; ++j){
      bbb[k] = (b&bmask) | bv_buf[j] | bvt_buf[i];
      ccc[k] = c;
      ++k;
    }
  }

  return k-l;
}


/**********************************/
/****  subroutines for lattice ****/
/**********************************/
static int pivot_reduction(eqdeg_t *eq, int v)
{
  Vector **lattice, *ltmp;
  int i;
  int pivot;
  int count;
  int min;

  eq->upper_v_bits = 0;
  for( i=0; i<v; i++) {
    eq->upper_v_bits |= eq->bitmask[i];
  }

  lattice = make_lattice(eq, v );

  for (;;) {
    pivot = calc_pivot(lattice[v]->next);
    if (lattice[pivot]->count < lattice[v]->count) {
      ltmp = lattice[pivot];
      lattice[pivot] = lattice[v];
      lattice[v] = ltmp;
    }
    add(eq->nnn, lattice[v], lattice[pivot]);
    if (lattice[v]->next == 0) {
      count = 0;
      next_state(eq, lattice[v], &count);
      if (lattice[v]->next == 0) {
        if (is_zero(eq->nnn, lattice[v])) {
          break;
        }
        while (lattice[v]->next == 0) {
          count++;
          next_state(eq, lattice[v], &count);
          if (count > eq->nnn * (eq->www-1) - eq->rrr) {
            break;
          }
        }
        if (lattice[v]->next == 0) {
          break;
        }
      }
    }
  }

  min = lattice[0]->count;
  for (i = 1; i < v; i++) {
    if (min > lattice[i]->count) {
      min = lattice[i]->count;
    }
  }
  free_lattice( lattice, v );
  return min;
}


/********************************/
/** allocate momory for Vector **/
/********************************/
static Vector *new_Vector(int nnn)
{
  Vector *v;

  v = (Vector *)malloc( sizeof( Vector ) );
  if( v == NULL ){
    printf("malloc error in \"new_Vector()\"\n");
    exit(1);
  }

  v->cf = (ulong *)calloc( nnn, sizeof( ulong ) );
  if( v->cf == NULL ){
    printf("calloc error in \"new_Vector()\"\n");
    exit(1);
  }

  v->start = 0;

  return v;
}


/************************************************/
/* frees *v which was allocated by new_Vector() */
/************************************************/
static void free_Vector( Vector *v )
{
  if( NULL != v->cf ) free( v->cf );
  if( NULL != v ) free( v );
}

static void free_lattice( Vector **lattice, int v)
{
  int i;

  for( i=0; i<=v; i++)
    free_Vector( lattice[i] );
  free( lattice );
}

/* adds v to u (then u will change) */
static void add(int nnn, Vector *u, Vector *v)
{
  int i;
  int diff = (v->start - u->start + nnn) % nnn;
  for (i = 0; i < nnn - diff; i++) {
    u->cf[i] ^= v->cf[i + diff];
  }
  diff = diff - nnn;
  for (; i < nnn; i++) {
    u->cf[i] ^= v->cf[i + diff];
  }
  u->next ^=  v->next;
}

/* makes a initial lattice */
static Vector **make_lattice(eqdeg_t *eq, int v)
{
  int i;
  int count;
  Vector **lattice, *bottom;

  lattice = (Vector **)malloc( (v+1) * sizeof( Vector *) );
  if( NULL == lattice ){
    printf("malloc error in \"make_lattice\"\n");
    exit(1);
  }

  for( i=0; i<v; i++){ /* from 0th row to v-1-th row */
    lattice[i] = new_Vector(eq->nnn);
    lattice[i]->next = eq->bitmask[i];
    lattice[i]->start = 0;
    lattice[i]->count = 0;
  }

  bottom = new_Vector(eq->nnn); /* last row */
  for(i=0; i< eq->nnn; i++) {
    bottom->cf[i] = 0;
  }
  bottom->cf[eq->nnn -1] = 0xc0000000 & eq->greal_mask;
  bottom->start = 0;
  bottom->count = 0;
  count = 0;
  do {
    next_state(eq, bottom, &count);
  } while (bottom->next == 0);
  //    degree_of_vector(eq, top );
  lattice[v] = bottom;

  return lattice;
}

static void next_state(eqdeg_t *eq, Vector *v, int *count) {
  ulong tmp;

  do {
    tmp = ( v->cf[v->start] & eq->gupper_mask )
      | ( v->cf[(v->start + 1) % eq->nnn] & eq->glower_mask );
    v->cf[v->start] = v->cf[(v->start + eq->mmm) % eq->nnn]
      ^ ( (tmp>>1) ^ eq->aaa[lsb(eq, tmp)] );
    v->cf[v->start] &= eq->greal_mask;
    tmp = v->cf[v->start];
    v->start = (v->start + 1) % eq->nnn;
    v->count++;
    tmp = trnstmp(eq, tmp);
    tmp = masktmp(eq, tmp);
    v->next = tmp & eq->upper_v_bits;
    (*count)++;
    if (*count > eq->nnn * (eq->www-1) - eq->rrr) {
      break;
    }
  } while (v->next == 0);
}

/***********/
static MaskNode *cons_MaskNode(MaskNode *head, ulong b, ulong c, int leng)
{
  MaskNode *t;

  t = (MaskNode*)malloc(sizeof(MaskNode));
  if (t == NULL) {
    printf("malloc error in \"cons_MaskNode\"\n");
    exit(1);
  }

  t->b = b;
  t->c = c;
  t->leng = leng;
  t->next = head;

  return t;
}

static void delete_MaskNodes(MaskNode *head)
{
  MaskNode *t;

  while(head != NULL) {
    t = head->next;
    free(head);
    head = t;
  }
}

static MaskNode *delete_lower_MaskNodes(MaskNode *head, int l)
{
  MaskNode *s, *t, *tail;

  s = head;
  while(1) { /* heading */
    if (s == NULL)
      return NULL;
    if (s->leng >= l)
      break;
    t = s->next;
    free(s);
    s = t;
  }

  head = tail = s;

  while (head != NULL) {
    t = head->next;
    if (head->leng < l) {
      free(head);
    }
    else {
      tail->next = head;
      tail = head;
    }
    head = t;
  }

  tail->next = NULL;
  return s;
}

/*
   w -- word size fixed to 32
   p -- Mersenne Exponent
   */
#define DEFAULT_ID_SIZE 16
/* id <= 0xffff */
static dcmt_state *get_mt_parameter_id_st(int p, int id, ulong seed)
{
  dcmt_state *mts;
  prescr_t pre;
  mt_state org;
  check32_t ck;

  pnl_mt_sseed(&org, seed);
  if (id > 0xffff) {
    printf("\"id\" must be less than 65536\n");
    return NULL;
  }
  if (id < 0) {
    printf("\"id\" must be positive\n");
    return NULL;
  }

  mts = init_mt_search(&ck, &pre, p);
  if (mts == NULL) return NULL;

  if ( NOT_FOUND == get_irred_param(&ck, &pre, &org,
                                    mts, id, DEFAULT_ID_SIZE) ) {
    pnl_dcmt_free(&mts);
    return NULL;
  }
  _get_tempering_parameter_hard_dc(mts);
  end_mt_search(&pre);

  return mts;
}

static dcmt_state **get_mt_parameters_st(int p, int start_id,
                                  int max_id, ulong seed, int *count)
{
  dcmt_state **mtss, *template_mts;
  int i;
  prescr_t pre;
  mt_state org;
  check32_t ck;

  if ((start_id > max_id) || (max_id > 0xffff) || (start_id < 0)) {
    printf("\"id\" error\n");
    return NULL;
  }

  pnl_mt_sseed(&org, seed);
  mtss = (dcmt_state**)malloc(sizeof(dcmt_state*)*(max_id-start_id+1));
  if (NULL == mtss) return NULL;

  template_mts = init_mt_search(&ck, &pre, p);
  if (template_mts == NULL) {
    free(mtss);
    return NULL;
  }
  *count = 0;
  for (i=0; i<=max_id-start_id; i++) {
    mtss[i] = alloc_dcmt_state(template_mts->nn);
    if (NULL == mtss[i]) {
      break;
    }

    copy_params_of_dcmt_state(template_mts, mtss[i]);

    if ( NOT_FOUND == get_irred_param(&ck, &pre, &org, mtss[i],
                                      i+start_id,DEFAULT_ID_SIZE) ) {
      pnl_dcmt_free(&(mtss[i]));
      break;
    }
    _get_tempering_parameter_hard_dc(mtss[i]);
    ++(*count);
  }

  pnl_dcmt_free(&template_mts);
  end_mt_search(&pre);
  if (*count > 0) {
    return mtss;
  } else {
    free(mtss);
    return NULL;
  }
}

/**
 * Create a MT generator. It returns the same generator upon two calls with
 * the same value of seed
 * 
 * @param seed : seed for original mt19937 to generate parameter.
 */
dcmt_state* pnl_dcmt_get_parameter(ulong seed)
{
  /* we use 521 as a Mersenne exponent */
  return get_mt_parameter_st (521, seed);
}

/**
 * Create a MT generator. It returns the same generator upon two calls with
 * the same value of seed
 * 
 * @param id: an integer between 0 and 65536. Choosing different id ensures
 * to get independent generators
 * @param seed : seed for original mt19937 to generate parameter.
 */
dcmt_state* pnl_dcmt_get_parameter_id (int id, ulong seed)
{
  /* we use 521 as a Mersenne exponent */
  return get_mt_parameter_id_st (521, id, seed);
}

/**
 * Create a MT generator. When called several times, the returned generators
 * are independent
 */
int pnl_dcmt_create (dcmt_state *mts)
{
  dcmt_state *template_mts;
  prescr_t pre;
  mt_state org;
  check32_t ck;
  static int id = 0;
  
  pnl_mt_sseed(&org, 4172);
  if ( (template_mts = init_mt_search(&ck, &pre, 521)) == NULL ) return FAIL;

  copy_params_of_dcmt_state(template_mts, mts);
  if ( NOT_FOUND == get_irred_param(&ck,&pre,&org,mts,id,DEFAULT_ID_SIZE) )
    {
      pnl_dcmt_free(&mts);
      return FAIL;
    }
  _get_tempering_parameter_hard_dc(mts);
  end_mt_search(&pre);
  id++;

  return OK;
}

/** 
 * Create a DCMT generator with identifier id
 * 
 * @param id the identifier of generator. Tow generators with different ids
 * are independent
 * @param seed the seed used to initialize the internal MT used to find the
 * parameter set. 
 * 
 * @return a PnlRng with PNL_RNG_DCMT type. Note that this generator must
 * be initialized with pnl_rng_sseed before usage.
 */
PnlRng* pnl_rng_dcmt_create_id (int id, ulong seed)
{
  PnlRng *rng;
  dcmt_state *state;

  state = pnl_dcmt_get_parameter_id (id, seed);
  rng = pnl_rng_create (PNL_RNG_DCMT);
  rng->state = state;
  return rng;
}

/**
 * Create an array of DCMT
 *
 * @param n number of generators to be created
 * @param seed the seed used to initialise the standard MT used to find new
 * DCMT
 * @param count (output) contains the number of generators actually created
 */
dcmt_state** pnl_dcmt_create_array(int n, ulong seed, int *count)
{
  int start_id, max_id, p;

  start_id = 0;
  max_id = n - 1;
  p = 521; /* smallest Mersenne exponent */
  return get_mt_parameters_st (p, start_id, max_id, seed, count);
}

/**
 * Create an array of PnlRng with PNL_RNG_DCMT type. Each PnlRng is
 * associated to a different id.  Note that each generator must
 * be initialized with pnl_rng_sseed before usage.
 *
 *
 * @param start_id smallest id
 * @param max_id largest id
 * @param seed the seed used to initialise the standard MT used to find new
 * DCMT
 * @param count (output) contains the number of generators actually
 * created. It should be max_id - start_id + 1
 */
PnlRng** pnl_rng_dcmt_create_array_id (int start_id, int max_id, ulong seed, int *count)
{
  int i, p, n;
  dcmt_state **states;
  PnlRng **rngs;
  n = max_id - start_id + 1;
  p = 521; /* smallest Mersenne exponent */
  states = get_mt_parameters_st (p, start_id, max_id, seed, count);

  PNL_MESSAGE (*count !=n, "Not all generators could be created");

  rngs = malloc (*count * sizeof (PnlRng *));
  for ( i=0 ; i<*count ; i++ )
    {
      rngs[i] = pnl_rng_create (PNL_RNG_DCMT);
      rngs[i]->state = states[i];
    }

  free (states);
  return rngs;
}


/**
 * Create an array of PnlRng filled with DCMT
 *
 * @param n number of generators to be created
 * @param seed the seed used to initialise the standard MT used to find new
 * DCMT
 * @param count (output) contains the number of generators actually created
 */
PnlRng** pnl_rng_dcmt_create_array (int n, ulong seed, int *count)
{
  return pnl_rng_dcmt_create_array_id (0, n-1, seed, count);
}

/**
 * Set the seed of a DCMT
 *
 * @param s an unsigned integer used as seed
 * @param mts a dcmt
 */
void pnl_dcmt_sseed(dcmt_state *mts, ulong s) 
{
  int i;

  /* when s==0, we fix an arbitrary seed but always the same */
  if ( s == 0 )
    {
      s = 4357;
    }

  mts->state[0] = s & mts->wmask;
  mts->state[0] &= BIT32_MASK;

  for ( i=1 ; i<mts->nn ; i++ ) 
    {
      mts->state[i] = ((1812433253UL * (mts->state[i-1]  ^ (mts->state[i-1] >> 30))) + i) & mts->wmask;
      /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
      /* In the previous versions, MSBs of the seed affect   */
      /* only MSBs of the array mt[].                        */

      /* 32bit mask */
      mts->state[i] &= BIT32_MASK;
    }
  mts->i = mts->nn;
}

/**
 * Generate an unsigned integer using the generator specified by @param mts
 */
ulong pnl_dcmt_genrand(dcmt_state *mts) 
{
  ulong *st, uuu, lll, aa, x;
  int k,n,m,lim;

  if ( mts->i >= mts->nn ) 
    {
      n = mts->nn; m = mts->mm;
      aa = mts->aaa;
      st = mts->state;
      uuu = mts->umask; lll = mts->lmask;

      lim = n - m;
      for (k=0; k<lim; k++) 
        {
          x = (st[k]&uuu)|(st[k+1]&lll);
          st[k] = st[k+m] ^ (x>>1) ^ (x&1U ? aa : 0U);
        }
      lim = n - 1;
      for (; k<lim; k++) 
        {
          x = (st[k]&uuu)|(st[k+1]&lll);
          st[k] = st[k+m-n] ^ (x>>1) ^ (x&1U ? aa : 0U);
        }
      x = (st[n-1]&uuu)|(st[0]&lll);
      st[n-1] = st[m-1] ^ (x>>1) ^ (x&1U ? aa : 0U);
      mts->i=0;
    }

  x = mts->state[mts->i];
  x ^= x >> mts->shift0;
  x ^= (x << mts->shiftB) & mts->maskB;
  x ^= (x << mts->shiftC) & mts->maskC;
  x ^= x >> mts->shift1;

  mts->i++;

  return x;
}

/**
 * Generate a double number between 0 and 1 using the generator specified by @param mts
 */
double pnl_dcmt_genrand_double (dcmt_state *mts)
{
  return pnl_dcmt_genrand(mts) / 4294967296.0; 
}



