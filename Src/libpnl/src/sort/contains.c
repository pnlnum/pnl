/*************************************************************************/
/* Written and (C) by PNL team, 2009                                     */  
/*                                                                       */
/*                                                                       */
/* This program is free software; you can redistribute it and/or modify  */
/* it under the terms of the GNU General Public License as published by  */
/* the Free Software Foundation; either version 3 of the License, or     */
/* (at your option) any later version.                                   */
/*                                                                       */
/* This program is distributed in the hope that it will be useful,       */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/* GNU General Public License for more details.                          */
/*                                                                       */
/* You should have received a copy of the GNU General Public License     */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/*************************************************************************/

#include "pnl_list.h"

/**
 * allocates a contains.
 * @param ind key 
 * @param val value
 * @return a pointeur to PnlContains
 */
PnlContains * pnl_contains_create(const int ind,double val)
{
  PnlContains *C;
  if((C=malloc(sizeof(PnlContains)))==NULL) return NULL;
  C->index=ind;
  C->value=val;
  return C;
}

PnlContains * pnl_contains_clone(int ind,double val)
{
  PnlContains *C;
  if((C=malloc(sizeof(PnlContains)))==NULL) return NULL;
  C->index=ind;
  C->value=val;
  return C;
}

/**
 * allocates a contains - copy constructor.
 * @param C2 contains pointer 
 * @return a pointeur to PnlContains
 */
PnlContains * pnl_contains_copy(const PnlContains *C2)
{
  PnlContains *C;
  if((C=malloc(sizeof(PnlContains)))==NULL) return NULL;
  C->index=C2->index;
  C->value=C2->value;
  return C;
}

/**
 * free a contains
 * @param C address of a contains
 */
void pnl_contains_free(PnlContains **C)
{
  if (*C != NULL)
    {
      free(*C);
      *C=NULL;
    }
}

/**
 * Prints a contains to a file 
 *
 * @param fic a file descriptor.
 * @param C a Contians pointer.
 */
void pnl_contains_fprint(FILE *fic,PnlContains *C)
{fprintf(fic," ( %d,  %7.4f) ;",C->index,C->value);} 


/**
 * Add at value, the value of PnlContains C2  
 *
 * @param C a PnlContains pointer, C.Value+ =C2 .Value.
 * @param C2 a Contians pointer.
 */
void  pnl_contains_add(PnlContains *C,const PnlContains *C2)
{C->value+=C2->value;}

/**
 * Less compute relation C1<C2
 *
 * @param C1 a PnlContains pointer.
 * @param C2 a Contians pointer.
 * @return a int C1<C2
 */
int  pnl_contains_less(const PnlContains *C1,const PnlContains *C2)
{return C1->index<C2->index;}

/**
 * Equal compute relation C1==C2
 *
 * @param C1 a PnlContains pointer.
 * @param C2 a Contians pointer.
 * @return a int C1==C2
 */
int  pnl_contains_equal(const PnlContains *C1,const PnlContains *C2)
{return C1->index==C2->index;}

