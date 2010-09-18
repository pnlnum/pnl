#include "pnl/pnl_internals.h"

/* Function dpmpar                                                           */

/* This function provides double precision machine parameters                */
/* when the appropriate set of data statements is activated (by              */
/* removing the c from column 1) and all other data statements are           */
/* rendered inactive. Most of the parameter values were obtained             */
/* from the corresponding Bell Laboratories Port Library function.           */

/* The function statement is                                                 */

/*   double precision function pnl_minpack_dpmpar(i)                         */

/* where                                                                     */

/*   i is an integer input variable set to 1, 2, or 3 which                  */
/*     selects the desired machine parameter. If the machine has             */
/*     t base b digits and its smallest and largest exponents are            */
/*     emin and emax, respectively, then these parameters are                */

/*     pnl_minpack_dpmpar(1) = b**(1 - t), the machine precision,            */

/*     pnl_minpack_dpmpar(2) = b**(emin - 1), the smallest magnitude,        */

/*     pnl_minpack_dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude. */

/* Argonne National Laboratory. MINPACK Project. November 1996.              */
/* Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More'                    */


double pnl_minpack_dpmpar(int i)
{

  switch (i)
    {
      case 1: return pnl_dlamch ("p"); break;
      case 2: return pnl_dlamch ("u"); break;
      case 3: return pnl_dlamch ("o"); break;
    }
  return 0.0;

}

