
/*************************************************************************/
/* Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>            */
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pnl_bessel.h"

static void complex_bessel_test ()
{
  dcomplex z, c;
  double v = -1.5;
  z = Complex (5., 3.);
  printf("Test des fonctions de bessel complexes \n");
  c = pnl_complex_bessel_i (v, z);
  printf("bessel_i(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_i_scaled (v, z);
  printf("bessel_i_scaled(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_j (v, z);
  printf("bessel_j(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_j_scaled (v, z);
  printf("bessel_j_scaled(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_y (v, z);
  printf("bessel_y(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_y_scaled (v, z);
  printf("bessel_y_scaled(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_k (v, z);
  printf("bessel_k(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_k_scaled (v, z);
  printf("bessel_k_scaled(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_h1 (v, z);
  printf("bessel_h1(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_h1_scaled (v, z);
  printf("bessel_h1_scaled(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_h2 (v, z);
  printf("bessel_h2(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);
  c = pnl_complex_bessel_h2_scaled (v, z);
  printf("bessel_h2_scaled(%f, %f + i %f) = %f + %f i\n", v, z.r, z.i, c.r, c.i);

}

static void real_bessel_test ()
{
  double z, c;
  dcomplex zc;
  double v = -1.5;
  z = 5.;
  printf("Test des fonctions de bessel reelles \n");
  c = pnl_bessel_i (v, z);
  printf("bessel_i(%f, %f) = %f\n", v, z, c);
  c = pnl_bessel_i_scaled (v, z);
  printf("bessel_i_scaled(%f, %f) = %f\n", v, z, c);
  c = pnl_bessel_j (v, z);
  printf("bessel_j(%f, %f) = %f\n", v, z, c);
  c = pnl_bessel_j_scaled (v, z);
  printf("bessel_j_scaled(%f, %f) = %f\n", v, z, c);
  c = pnl_bessel_y (v, z);
  printf("bessel_y(%f, %f) = %f\n", v, z, c);
  c = pnl_bessel_y_scaled (v, z);
  printf("bessel_y_scaled(%f, %f) = %f\n", v, z, c);
  c = pnl_bessel_k (v, z);
  printf("bessel_k(%f, %f) = %f\n", v, z, c);
  c = pnl_bessel_k_scaled (v, z);
  printf("bessel_k_scaled(%f, %f) = %f\n", v, z, c);
  zc = pnl_bessel_h1 (v, z);
  printf("bessel_h1(%f, %f) = %f + %f i\n", v, z, CMPLX(zc));
  zc = pnl_bessel_h1_scaled (v, z);
  printf("bessel_h1_scaled(%f, %f) = %f + %f i\n", v, z, CMPLX(zc));
  zc = pnl_bessel_h2 (v, z);
  printf("bessel_h2(%f, %f) = %f + %f i\n", v, z, CMPLX(zc));
  zc = pnl_bessel_h2_scaled (v, z);
  printf("bessel_h2_scaled(%f, %f) = %f + %f i\n", v, z, CMPLX(zc));
}


void bessel_test()
{
  printf("\n");
  printf("TEST DES FONCTIONS DE BESSEL \n");
  printf("\n");
  complex_bessel_test ();
  real_bessel_test ();
}
