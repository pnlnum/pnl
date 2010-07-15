
/************************************************************************/
/* Copyright Jérôme Lelong <jerome.lelong@gmail.com>                    */
/*                                                                      */
/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as       */
/* published by the Free Software Foundation, either version 3 of the   */
/* License, or (at your option) any later version.                      */
/*                                                                      */
/* This program is distributed in the hope that it will be useful, but  */
/* WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    */
/* Lesser General Public License for more details.                      */
/*                                                                      */
/* You should have received a copy of the GNU Lesser General Public     */
/* License along with this program.  If not, see                        */
/* <http://www.gnu.org/licenses/>.                                      */
/************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl_random.h"

static double square (double x) { return x*x;}

static void test_pnl_vect_rand()
{
    int samples = 10000;
    int type_generator = PNL_RNG_TAUSWORTHE;
    double sum;
    int i;
    PnlVect *G;
    G = pnl_vect_create(0);
    pnl_rand_init(type_generator, 1, samples);

    pnl_vect_rand_uni(G, samples, 0, 1, type_generator);
    pnl_vect_map_inplace(G, square);
    sum = pnl_vect_sum(G)/samples;
    printf("E (U^2) = %f (should be 1/3)\n", sum);

    /* Calling pnl_rand_init again ensures that the next samples will be drawn
       from the beginning of the sequence. Very important for QMC */
    pnl_rand_init(type_generator, 1, samples);    
    sum=0.0;
    for (i=0; i<samples; i++)
    {
        sum += square(pnl_rand_normal(type_generator));
    }
     printf("E (G^2) = %f (should be 1)\n", sum/samples);

     /* equivalently as above */
    pnl_rand_init(type_generator, 1, samples);    
    pnl_vect_rand_normal(G, samples, type_generator);
    pnl_vect_map_inplace(G, square);
    sum = pnl_vect_sum(G)/samples;
    printf("E (G^2) = %f (should be 1)\n", sum);

    pnl_vect_free(&G);
}

static void test_pnl_mat_rand()
{
    PnlMat *M;
    PnlVect *inf, *sup, *Vsum;
    int type_generator = PNL_RNG_KNUTH;
    int samples = 100000;
    int dim = 10;
    int i;
    double sum;
    inf = pnl_vect_create_from_double(dim, 0.0);
    sup = pnl_vect_create_from_double(dim, 1.0);
    M = pnl_mat_create(0,0);
    Vsum = pnl_vect_create (0);
    pnl_rand_init(type_generator, dim, samples);
    pnl_mat_rand_uni(M, samples, dim, inf, sup, type_generator);
    pnl_mat_sum_vect(Vsum, M, 'r');
    pnl_vect_div_double(Vsum, samples);
    pnl_vect_print(Vsum);
    printf("(should be 0.5)\n");
    pnl_vect_free(&Vsum);
    pnl_vect_free(&inf);
    pnl_vect_free(&sup);
    
    
    pnl_rand_init(type_generator, dim, samples);
    pnl_mat_rand_normal(M, samples, dim, type_generator);
    sum = 0.0;
    for (i=0; i<samples; i++)
    {
        sum += MGET(M,i,0)*MGET(M,i,8);
    }
    printf("Cov = %f (should be 0)\n", sum/samples);

    pnl_mat_free(&M);
}

static void test_pnl_rand_gauss()
{
    int dimension=10, samples=100000, type_generator=2;
    int i;
    double g1,g2,sum=0.0;
    pnl_rand_init(type_generator, dimension, samples);
    for (i=0; i<samples; i++)
    {
        g1=pnl_rand_gauss(dimension, CREATE, 0, type_generator);
        g2=pnl_rand_gauss(dimension, RETRIEVE, 8, type_generator);
        sum += g1*g2;
    }
    printf("Cov = %f (should be 0)\n", sum/samples);
}

void random_test(void)
{
    test_pnl_vect_rand();
    test_pnl_mat_rand();
    test_pnl_rand_gauss();
}
