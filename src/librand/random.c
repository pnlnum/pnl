/* 
 * This file contains two parts.  The first one gathers the random generator
 * functions which must never be called directly. These codes are provided by
 * Premia. The second part is an interface to these generators to easily
 * create vectors and matrices of random numbers according to different
 * laws. 
 */



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

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#define _PNL_PRIVATE 1

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_random.h"
#include "pnl/pnl_cdf.h"
#include "randomkit.h"

/* ---------------------------------- */
/* PSEUDO RANDOM NUMBERS GENERATORS.  */
/* ---------------------------------- */

/* ------------------------------------------------------------------- */
/* Random numbers generator of Knuth :                                 */
/* It is  based on MRG and it uses a substractive method               */
/* ------------------------------------------------------------------- */

static void KNUTH(PnlRng *rng,double *sample)
{
  long X_n, y_k;
  int i, ii, l;
  knuth_state *state;

  static const long M    = 1000000000L;
  static const long alea = 1L; 

  state = (knuth_state *) rng->state;
  
  /* First call to the sequence */
  if(rng->counter == 1)
    {
      X_n= state->SEED - alea;
      X_n%= M;
      state->t_alea[55]= X_n;
      y_k= 1;
      /* Initialization of the table */
      for(i= 1; i<= 54; i++)
        {
          ii= (21*i)%55; /* 21 was chosen to alleviate initial
                            nonrandomness problems */
          state->t_alea[ii]= y_k;
          y_k= X_n - y_k;
          if(y_k < 0) 
            y_k+= M;
          X_n= state->t_alea[ii];
        }

      /* Randomization of the elements of the table */   
      for(l=  1; l<= 4; l++)
        {
          for(i= 1; i<= 55; i++)
            {
              state->t_alea[i]-= state->t_alea[1+(i+30)%55];
              if(state->t_alea[i] < 0) 
                state->t_alea[i]+= M;
            }
        }
      state->inc1= 0;
      state->inc2= 31;  /* 31 is a special value of Knuth : 31= 55-24 */
    }

  rng->counter++;

  /* For each call to the sequence, computation of a new point */
  if(++(state->inc1) == 56) 
    state->inc1= 1;
  if(++(state->inc2) == 56) 
    state->inc2= 1;
  /* Substractive method*/
  X_n= state->t_alea[state->inc1] - state->t_alea[state->inc2];

  if(X_n < 0) 
    X_n+= M;
  state->t_alea[state->inc1]= X_n;
  /* Normalized value */
  *sample = (double) X_n / (double) M;
  return;
}

/* ----------------------------------------------------------------------- */
/* Combination of two multiplicative recursive generators of order 3 (k=3) */
/* ----------------------------------------------------------------------- */
static void MRGK3(PnlRng *rng, double *sample)
{
  long k;
  double p1, p2;
  mrgk3_state *s;

  static const double M1   = 4294967087.;
  static const double M2   = 4294944443.;
  static const double A12  = 1403580.;
  static const double A13N = 810728.;
  static const double A21  = 527612.;
  static const double A23N = 1370589.;


  s = (mrgk3_state *) (rng->state);
  rng->counter++; 

  /* First generator */
  p1 = A12*s->x11 - A13N*s->x10;
  k = p1 / M1; 
  p1 -= k*M1;
  if(p1 < 0) p1 += M1;

  s->x10 = s->x11;
  s->x11 = s->x12;
  s->x12 = p1;

  /* Second generator */
  p2 = A21*s->x22 - A23N*s->x20;
  k = p2 / M2;
  p2 -= k*M2;
  if(p2 < 0) p2 += M2;

  s->x20= s->x21;
  s->x21= s->x22;
  s->x22= p2;


  /* Combination of the two generators */
  if (p1< p2) 
    *sample= (p1- p2+ M1) / (double) M1;
  else 
    *sample=(p1- p2) / (double) M1;
  return;
}



/* ----------------------------------------------------------------------- */
/* Combination of two multiplicative recursive generators of order 5 (k=5) */
/* ----------------------------------------------------------------------- */
static void MRGK5(PnlRng *rng,double *sample)
{
  long k;
  double p1, p2;
  mrgk5_state *s;

  static const double M1   = 4294949027.;
  static const double M2   = 4294934327.;
  static const double A12  = 1154721.;
  static const double A14  = 1739991.;
  static const double A15N = 1108499.;
  static const double A21  = 1776413.;
  static const double A23  = 865203.;
  static const double A25N = 1641052.;

  s = (mrgk5_state *)(rng->state);
  
  rng->counter++;

  /* For each call to the sequence, computation of a new point */
  /* First generator with Schrage method */
  p1= A12*s->x13 - A15N*s->x10;
  if(p1> 0) p1 -= A14*M1;
  p1 += A14*s->x11;
  k = p1/M1;
  p1 -= k*M1;

  if (p1< 0) p1+= M1;

  s->x10= s->x11;
  s->x11= s->x12;
  s->x12= s->x13; 
  s->x13= s->x14;
  s->x14= p1;

  /* Second generator with Schrage method */
  p2= A21*s->x24 - A25N*s->x20;
  if (p2> 0) p2-= A23*M2;
  p2 += A23*s->x22;
  k = p2 / M2;
  p2 -= k*M2;

  if (p2 < 0) p2+= M2;

  s->x20= s->x21;
  s->x21= s->x22;
  s->x22= s->x23;
  s->x23= s->x24;
  s->x24= p2;


  /*Combination of the two generators */
  if (p1 <= p2) 
    *sample = (p1 - p2 + M1) / (double) M1;
  else 
    *sample = (p1 - p2) / (double) M1;

  return;
}
/* ------------------------------------------------------------------------ */
/* Random numbers generator of  Park & Miller with Bayes & Durham shuffling
   procedure : the next random number is not obtained from the previous one 
   but we use an intermediate table which contains the 32 precedent random 
   numbers and we choose one of them randomly.  */
/* ----------------------------------------------------------------------- */
static void SHUFL(PnlRng *rng,double *sample)
{
  long N1;
  int j;
  long hi;                 /* high order bit */
  shufl_state *s = (shufl_state *)(rng->state);

  static const long A = 16807;        /* multiplier */
  static const long M = 2147483647;   /* 2**31 - 1  */
  static const long Q = 127773;       /*   M div A  */  
  static const long R = 2836;         /*   M mod A  */


  N1=(M/32);

  /* First call to the sequence */
  if (rng->counter == 1)
    {
      /* After 8 "warm-ups", initialisation of the shuffle table */
      for (j= 39; j>= 0; j--)
        {
          hi= s->x/Q;
          /*Schrage's method to avoid overflows */
          s->x= A*(s->x- hi*Q)- R*hi; 
          if (s->x < 0) s->x+= M;
          if (j< 32) s->t[j]= s->x;
        }
      s->y= s->t[0];                    
    }
  rng->counter++;


  /* For each call to the sequence, computation of a new point */
  hi= s->x/Q;
  s->x= A*(s->x-hi*Q)- R*hi;
  if (s->x < 0) s->x += M;

  /* Shuffling procedure of Bayes & Durham */
  /* Index j dependent on the last point */
  j= s->y/N1;
  /* Next point dependent on j */
  s->y= s->t[j];

  s->t[j]= s->x;
  *sample = (double) s->y / (double) M; 
}



/* ------------------------------------------------------------------------- */
/* Random numbers generator of L'Ecuyer with Bayes & Durham shuffling
   procedure : 
   Combination of two short periods LCG to obtain a longer period generator.
   The period is the least common multiple of the 2 others.  */
/* ------------------------------------------------------------------------- */

static void LECUYER(PnlRng *rng,double *sample)
{ 
  static long A1 = 40014;        /* multiplier of the 1st generator */
  static long A2 = 40692;        /* multiplier of the 2nd generator */
  static long M1 = 2147483647;   /* 2**31 - 1   */ 
  static long M2 = 2147483399;   /* 2**31 - 249 */
  static long Q1 = 53668;        /* m1 div a1   */  
  static long Q2 = 52774;        /* m2 div a2   */ 
  static long R1 = 12221;        /* m1 mod a1   */
  static long R2 = 3791;         /* m2 mod a2   */
  long N1;

  static long x;
  static long y= 978543162;
  int j;
  long hi;               /* high order bit */
  static long z= 0;
  static long t[32];     /* 32 is the size of a computer word */

  N1= (M1/32);

  /* First call to the sequence */
  if (rng->counter == 1)
    {
      x= 437651926;
      y= x;

      /* After 8 "warm-ups", initialisation of the shuffle table */
      for (j= 39; j>= 0; j--)
        {
          /* Park & Miller's generator */
          hi= x/Q1;  
          x= A1*(x-hi*Q1) - R1*hi;
          if (x < 0) 
            x+= M1;
          if (j < 32)
            t[j]= x;
        }
      z= t[0];
    }
  rng->counter++;              

  /* For each call to the sequence, computation of a new point */
  /* First generator */
  hi= x/Q1;
  x= A1*(x-hi*Q1) - R1*hi;
  if (x < 0) 
    x+= M1;

  /* Second generator */
  hi= y/Q2;
  y= A2*(y-hi*Q2) - R2*hi;
  if (y < 0) 
    y+= M2;

  /* Shuffling procedure of Bayes & Durham */
  /* Index j dependent on the last point */
  j= z/N1;
  /* Next point dependent on j */
  z= t[j]- y; 
  t[j]= x;


  /* To avoid 0 value */
  if (z < 1) 
    z+= M1-1;

  *sample = (double) z / (double) M1; 
  return;
}


/* ------------------------ */
/* Tausworthe Algorithm     */
/* ------------------------ */



/* ---------------------------------------------------- */
/* Generation of a random bit
   Algorithm based on a prime polynomial : 
   here we choose x^18 + x^5 + x^2 + x + 1 . */
/* ---------------------------------------------------- */
static int bit_random(void)
{
  static int compt = 1;
  int degre = 18;
  static unsigned long a;
  unsigned long new_bit;

  /* Initialisation for the n first values */
  /* random number over [1, 2^18] ; 2^18= 262144 */
  if(compt == 1)
    {
      a= 176355;
    }
  compt++;

  /* Next bit calculation by the recurrence relation  */
  new_bit= (a & (1<<17)) >> 17
    ^ (a & (1<<4)) >> 4
    ^ (a & (1<<1)) >> 1
    ^ (a & 1);
  a <<= 1;
  /* The new bit is shift on the right */
  a ^= new_bit;
  /* The most left bit is put to 0 */
  a ^= (1 << degre);

  return((int)new_bit);
}



/* --------------------------------------------------------------- */
/* Generation of a word of k random bits. */
/* --------------------------------------------------------------- */
static unsigned long random_word(int k)
{
  int i, bit;
  unsigned long mot;

  mot= 0;
  for(i= 0; i< k; i++)
    {
      bit= bit_random();
      mot= (mot<<1) ^ bit;
    }
  return(mot);
}

/* ---------------------------------------------------------------- */
/*  Tausworthe  Algorithm
    Combination  of 3 Tausworthe generators 
    u(n)[j]= u(n-r)[j] ^ u(n-k)[j], j=1,..,3
    with parameters k, q, r, s and t.
    Generator :
    v= =(u[0] ^ u[1] ^ u[2])/2^32.    

    L= 32 length of a word.
    We use a mask to make the generator 
    (originally designed for 32-bit unsigned int) 
    work on 64 bit machines : this mask enables to drop 
    the 32 highest bits on 64 bit machines*/
/* ---------------------------------------------------------------- */
#define TAUS_MASK 0xffffffffUL

static void TAUS(PnlRng *rng,double *sample)
{
  int i;
  int L= 32;
  static unsigned long u[3];
  static unsigned long c[3];
  unsigned long b;
  static int k[3], q[3];
  static int s[3], r[3], t[3];
  unsigned long v= 0;

  /* First call to the sequence. Initialisation */
  if (rng->counter == 1) 
    {
      /* Choice of the parameters to have ME-CF generators (cf L'Ecuyer) */
      k[0]= 31; q[0]= 13; s[0]= 12;
      r[0]= k[0]- q[0]; t[0]= k[0]- s[0];

      k[1]= 29; q[1]= 2; s[1]= 4;
      r[1]= k[1]- q[1]; t[1]= k[1]- s[1];

      k[2]= 28; q[2]= 3; s[2]= 17;
      r[2]= k[2]- q[2]; t[2]= k[2]-s[2];

      /* constant c : k bits to one and (L-k) bits to zero */
      /* c[j]= 2^32 - 2^(L-k[j]) */
      c[0]= 4294967294ul;
      c[1]= 4294967288ul;
      c[2]= 4294967280ul;

      /* initialisation of each generator */
      u[0]= 0; u[1]= 0; u[2]= 0;
      for(i= 0; i< 3; i++)
        {
          /* The first k bits are chosen randomly */
          u[i]= random_word(k[i]);
          /* The next L-k bits are initially fixed to zero */
          u[i] <<= (L- k[i]);  
          /* Now they are computed with the recurrence on u */
          b= u[i] << q[i];
          b ^= u[i];
          b >>= k[i];
          u[i] ^= b;
        }
    }

  /* For each call to the sequence, computation of a new point */
  for(i= 0; i< 3; i++)
    { 
      /* Calculus of the next points for the 3 generators */
      /* 6 steps explained by L'Ecuyer */
      b=((u[i]<<q[i]) & TAUS_MASK) ^ u[i];         /* Steps 1 and 2 */
      b >>= t[i];                    /* Step 3 */
      u[i]= ((u[i] & c[i]) << s[i]) & TAUS_MASK ;   /* Steps 4 et 5 */
      u[i] ^= b;                     /* Step 6 */
      /* Combination : XOR between the J generators */
      v^= u[i];
    }

  /* Normalization by 2^32 */
  *sample=v / 4294967296.0;
  rng->counter++; 
}

static void MERSENNE(PnlRng *rng,double *sample)
{
  mt_state *state = (mt_state *) (rng->state);

  *sample = pnl_mt_genrand_double (state);
  rng->counter++;
}

static void DYNAMIC_MT (PnlRng *rng,double *sample)
{
  dcmt_state *state = (dcmt_state *) (rng->state);

  *sample = pnl_dcmt_genrand_double (state);
  rng->counter++;
}

/* ------------------------------ */
/*       QUASI MONTE CARLO        */
/* ------------------------------ */

/* ----------------------------------------------------------------*/
/* QUASI RANDOM NUMBERS ; LOW DISCEPANCY SEQUENCES.
   SQRT, Van der Corput - Halton, Faure, Sobol, generalized Faure,
   Niederreiter */
/* ----------------------------------------------------------------*/

#define DIM_MAX_QMC 300

#define MAXI 33
static long Comb[MAXI][MAXI];/*Binomial Coefficients*/

/* -----------------------------------------*/
/* Table for the n first prime numbers */
/* -----------------------------------------*/
static void prime_number(int n, int prime[])
{
  int i, bool_prime;
  int tested_val= 3, nbre= 1;

  prime[0]= 2;

  while(nbre < n)
    {
      bool_prime = 0;
      for(i = 0; i < nbre; i++)
        {
          /* Test if the value is divisible by one of the first prime numbers */
          if((tested_val % prime[i]) == 0)
            {
              bool_prime = 1;
              break;
            }
        }
      /* If no divisor was found, the number is prime */
      if(bool_prime == 0)
        {
          prime[nbre]= tested_val;
          nbre += 1;
        }
      /* Test for odd integers only */
      tested_val+= 2;
    }
}



/* --------------------------------------------------------------------- */
/* Search the smallest element of a sorted table greater than a given
 * threshold.
 * Used for the prime numbers */
/* --------------------------------------------------------------------- */
static int search_value(int seuil, int tab[], int dim)
{
  int min, max, indice;

  min= 0;
  max= dim- 1;
  indice= (int)dim/2;

  while((max- min) >0)
    {
      if(tab[indice]< seuil)
        min= indice+1;
      else
        max= indice;
      indice= (min + max)/2;
    }
  return(tab[indice]);
}



/* -------------------------------------------*/
/* Computation of the binomial coefficients.  */
/* -------------------------------------------*/
static void binomial(int Max)
/*Max should be less than 33 otherwise C[n][p] could exceed LONG_MAX*/
{
  int n, p, i;

  Comb[0][0]= 1;
  for (n= 1; n< Max; n++) 
    {
      Comb[n][0]= 1;
      Comb[n][n]= 1;
      i= n-1;
      for (p= 1; p<= i; p++)
        Comb[n][p]= Comb[i][p]+ Comb[i][p-1];
    }
}

/* ----------------------------------------------------------------*/
/* SQRT Sequence */
/* Computation of the next element for the dim-dimensional sequence. */
/* ----------------------------------------------------------------*/

static void SQRT(PnlRng *rng, double X_m[])
{
  static int dimension=0;
  int i;
  static int prime[DIM_MAX_QMC];
  static double alpha[DIM_MAX_QMC];

  /* Verification of the dimension. It must not change without reinitializing */
  if(dimension != rng->dimension) rng->counter=1;

  /* First call : initialisation */
  if(rng->counter == 1)
    {
      dimension= rng->dimension;
      prime_number(rng->dimension, prime);
      for(i=0; i<rng->dimension; i++)
        {
          alpha[i]= sqrt((double) prime[i]);
        }
    }

  /* For each call to the sequence, computation of a new point */
  for(i= 0; i< rng->dimension; i++)
    X_m[i]= ((rng->counter*alpha[i])-floor(rng->counter*alpha[i]));

  rng->counter++;
}  


/* ------------------------------------------------------------- */
/* Halton sequence. Bases p={p1, p2, ..., pd} of prime numbers ; 
   Computation of the next element of the sequence X_n 
   (index 0 to d-1) */
/* ------------------------------------------------------------- */
static void HALTON(PnlRng *rng, double X_n[])
{
  static int dimension=0;
  int i;
  static int prime[DIM_MAX_QMC];
  int coeff;
  double y;
  int x, puissance;

  /* Verification of the dimension. It must not change without reinitializing */
  if(dimension != rng->dimension)
    rng->counter= 1;

  /* First call : initialization */
  if( rng->counter == 1)
    {
      dimension= rng->dimension;
      prime_number(rng->dimension, prime);
    }

  /* For each call to the sequence, computation of a new point */
  for(i = 0; i< rng->dimension; i++)
    {
      puissance= 1;
      x= 0; y= 0;
      /* radical inverse function */
      while ((rng->counter-x)>0) 
        {
          coeff= ((rng->counter-x)/puissance)%prime[i];
          x += coeff*puissance;
          puissance *= prime[i];
          y += coeff*1.0/puissance;
        }
      X_n[i]= y;
    }
  rng->counter++;
}

/* ----------------------------------------------------------------*/
/* FAURE sequence */
/* r smallest odd prime number greater than d, the dimension of the sequence
   Computation of the next element --> U_n[] (index 0 to d-1) */
/* ----------------------------------------------------------------*/
/* maximal dimension  for the FAURE sequence */
#define DIM_MAX_FAURE 300
/* maximal samples for the FAURE sequence */
#define MAX_SAMPLE_FAURE 10000000

static void FAURE(PnlRng *rng, double U_n[])
{
  int coeff[MAXI];
  int b[MAXI];
  static int dimension=0, r;
  int prime[DIM_MAX_QMC];
  int x= 0, puissance1= 1, puissance2;
  int indice, i, j, k, somm;


  /* Verification of the dimension. It must not change without reinitializing */
  if(dimension != rng->dimension)
    rng->counter= 1;

  /*First call to the sequence */
  if(rng->counter == 1)
    {
      dimension=rng->dimension;
      if((rng->dimension == 2)||(rng->dimension == 1))
        r= 3;
      else
        {
          prime_number(rng->dimension, prime);
          r=search_value(rng->dimension,prime,rng->dimension);
        }
    }

  /* Initialization */
  for (i=0; i< rng->dimension; i++) 
    {
      U_n[i]= 0.;
    }

  indice= 0;
  /* For each call to the sequence, computation of a new point */

  /* r-digit expansion of n --> first term of the sequence */
  while ((rng->counter-x) > 0)

    {
      coeff[indice]= ((rng->counter-x)/puissance1)%r;
      x += coeff[indice]*puissance1;
      puissance1 *= r;
      U_n[0] +=(double)coeff[indice]/(double)puissance1;
      indice +=1;
    }

  /* Other terms of the sequence */
  /* Successive transformations of the r-digit expansion.*/
  for (k=1; k< rng->dimension; k++) 
    {
      puissance2= r;
      for (j=0; j< indice; j++) 
        {
          somm= 0;
          for (i=j; i< indice; i++) 
            somm += Comb[i][j]*coeff[i];
          b[j]= somm%r;
        }
      for (j=0; j< indice; j++) 
        {
          coeff[j]= b[j];
          U_n[k] += (double)coeff[j]/(double)puissance2;
          puissance2 *= r;
        }
    }
  rng->counter++;
}


/* --------------------------------------------------------------- */
/* SOBOL SEQUENCE in dimension d (d<=39)
   X_n[] contains the terms of the n-th element
   The algorithm is based on a relation between X(n) and X(n-1)    */
/* Cf Numerical recipes in C, pages 312-313 */
/* --------------------------------------------------------------- */

/* maximal dimension for the SOBOL sequence */
#define DIM_MAX_SOBOL 39
/* maximal length of bits for the  SOBOL sequence */
#define BIT_MAX_SOBOL 30

static void SOBOL(PnlRng *rng, double X_n[])
{
  int i, j, k, P_j, deg_j;
  unsigned long aux, dim;
  static double facteur;
  static unsigned long initialX_n[DIM_MAX_SOBOL+1];

  /* Degree of the DIM_MAX primitive polynomials */
  static int deg[DIM_MAX_SOBOL+1]={0, 1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8};

  /* Index of each primitive polynomial. It determines the coefficients bi */
  static int P[DIM_MAX_SOBOL+1]={0, 0, 1, 1, 2, 1, 4, 2, 13, 7, 14, 11, 4, 1, 16, 13, 22, 19, 25, 1, 32, 4, 8, 7, 56, 14, 28, 19, 50, 21, 42, 31, 62,37, 41, 55, 59, 14, 56, 21};

  /* Sobol's constants C_i_(j) */
  static unsigned long C[DIM_MAX_SOBOL+1][BIT_MAX_SOBOL+1];
  /* Initial values for C */
  /* The second dimension is the maximum degree of a primitive polynomial +1 */
  static unsigned long C_init[DIM_MAX_SOBOL+1][9]={
    {0,0, 0, 0, 0,  0,  0,  0,   0},
    {0,1, 0, 0, 0,  0,  0,  0,   0},
    {0,1, 1, 0, 0,  0,  0,  0,   0},
    {0,1, 3, 7, 0,  0,  0,  0,   0},
    {0,1, 1, 5, 0,  0,  0,  0,   0},
    {0,1, 3, 1, 1,  0,  0,  0,   0},
    {0,1, 1, 3, 7,  0,  0,  0,   0},
    {0,1, 3, 3, 9,  9,  0,  0,   0},
    {0,1, 3, 7, 13, 3 , 0,  0,   0},
    {0,1, 1, 5, 11, 27, 0,  0,   0},
    {0,1, 3, 5, 1,  15, 0,  0,   0},
    {0,1, 1, 7, 3,  29, 0,  0,   0},
    {0,1, 3, 7, 7,  21, 0,  0,   0},
    {0,1, 1, 1, 9,  23, 37, 0,   0},
    {0,1, 3, 3, 5,  19, 33, 0,   0},
    {0,1, 1, 3, 13, 11, 7,  0,   0},
    {0,1, 1, 7, 13, 25, 5,  0,   0},
    {0,1, 3, 5, 11, 7,  11, 0,   0},
    {0,1, 1, 1, 3,  13, 39, 0,   0},
    {0,1, 3, 1, 15, 17, 63, 13,  0},
    {0,1, 1, 5, 5,  1,  27, 33,  0},
    {0,1, 3, 3, 3,  25, 17, 115, 0},
    {0,1, 1, 3, 15, 29, 15, 41,  0},
    {0,1, 3, 1, 7,  3,  23, 79,  0},
    {0,1, 3, 7, 9,  31, 29, 17,  0},
    {0,1, 1, 5, 13, 11, 3,  29,  0},
    {0,1, 3, 1, 9,  5,  21, 119, 0},
    {0,1, 1, 3, 1,  23, 13, 75,  0},
    {0,1, 3, 3, 11, 27, 31, 73,  0},
    {0,1, 1, 7, 7,  19, 25, 105, 0},
    {0,1, 3, 5, 5,  21, 9,  7,   0},
    {0,1, 1, 1, 15, 5,  49, 59,  0},
    {0,1, 1, 1, 1,  1,  33, 65,  0},
    {0,1, 3, 5, 15, 17, 19, 21,  0},
    {0,1, 1, 7, 11, 13, 29, 3,   0},
    {0,1, 3, 7, 5,  7,  11, 113, 0},
    {0,1, 1, 5, 3,  15, 19, 61,  0},
    {0,1, 3, 1, 1,  9,  27, 89,  7},
    {0,1, 1, 3, 7,  31, 15, 45,  23},
    {0,1, 3, 3, 9,  9,  25, 107, 39}
  };


  /*First call to the sequence */
  if(rng->counter == 1)
    {
      /* Initialization of the full array C[i][j] */
      for (j=1; j<=DIM_MAX_SOBOL; j++) 
        initialX_n[j]= 0;

      facteur=1.0/(1L << BIT_MAX_SOBOL);

      for (j=1; j<=DIM_MAX_SOBOL; j++) 
        {
          deg_j = deg[j];
          /* Values of C_init in C */
          for (i= 0; i<= deg_j; i++) 
            {
              C[j][i]= C_init[j][i];
            }
          /* Recurrence relation to compute the other C[i][j] */
          for (i= deg_j+1; i<= BIT_MAX_SOBOL; i++) 
            {
              P_j= P[j];
              aux= C[j][i-deg_j];
              aux ^= (C[j][i-deg_j] << deg_j);
              for (k= deg_j-1; k>=1; k--) 
                {
                  /* Test for the coefficient b(k) */
                  if (P_j & 1) 
                    aux ^= (C[j][i-k] << k);
                  P_j >>= 1;
                }
              /* Final value for C[i][j] */
              C[j][i]= aux;
            }
        }
    } 

  /* Calculation of a new quasi-random vector on each call */
  dim= rng->counter;
  /* Research of the rightmost 0 bit */
  for (i=1; i<=BIT_MAX_SOBOL; i++) 
    {
      if ((dim & 1) == 0) 
        break;
      dim >>= 1;
    }
  /* Computation of the term n from the term (n-1) */
  for (j=1; j<= rng->dimension; j++) 
    {
      initialX_n[j] = (initialX_n[j])^ (C[j][i]<< (BIT_MAX_SOBOL-i));
      /* normalization */
      X_n[j-1]= initialX_n[j] *facteur;
    }
  rng->counter++;
}


static void SOBOL2(PnlRng *rng, double X_n[])
{
  static rk_sobol_state s;
  static int dimension=0;

  if(rng->counter == 1 || dimension != rng->dimension)
    {
      if (dimension)
        pnl_rk_sobol_free(&s); /* We should free on exit too. */
      dimension=rng->dimension;
      pnl_rk_sobol_init(rng->dimension, &s, NULL, pnl_rk_sobol_Ldirections, NULL);
      /* For randomized QMC, add: rk_sobol_randomshift(&s, NULL); */
    }

  pnl_rk_sobol_double(&s, X_n);
  rng->counter ++;
  return;
}




/* 
 * NIEDERREITER SEQUENCE in dimension d, in base 2 
 * X_n[] contains the terms of the n-th element 
 * The algorithm is based on a relation between X(n) and X(n-1) 
 * ------------------------------------------------ */

/* maximal dimension for the NIEDERREITER sequence */
#define DIM_MAX_NIED 12
/* maximal length of bits for the NIEDERREITER sequence */
#define BIT_MAX_NIED 30

static void NIEDERREITER(PnlRng *rng, double X_n[])
{
  int i, j;
  unsigned long saut, gray, dim;
  static double facteur;
  static unsigned long initial_d, initialX_n[DIM_MAX_NIED+1];

  /* Niederreiter's constants */
  /* Une fonction de calcul de ces coefficients est donnee dans Owen.c */
  static unsigned long C[BIT_MAX_NIED+1][DIM_MAX_NIED+1]={
    {0,1073741824,1073741824,1610612736,1879048192,1879048192,2013265920,
     2013265920,2013265920,2080374784,2080374784,2080374784,2080374784},
    {0,536870912,1610612736,1207959552,1644167168,1644167168,1887436800,
     1887436800,1887436800,2015363072,2015363072,2015363072,2015363072},
    {0,268435456,1342177280,939524096,1174405120,1442840576,1635778560,
     1769996288,1769996288,1885339648,1885339648,1885339648,1952448512},
    {0,134217728,2013265920,2046820352,503316480,771751936,1132462080,
     1400897536,1535115264,1625292800,1692401664,1692401664,1759510528},
    {0,67108864,1140850688,1577058304,742391808,1312817152,260046848,
     796917760,1065353216,1172307968,1306525696,1239416832,1373634560},
    {0,33554432,1711276032,914358272,1522532352,515899392,386400256,
     1602748416,2131230720,266338304,467664896,333447168,601882624},
    {0,16777216,1426063360,1702887424,935329792,1069547520,639107072,
     932708352,1981284352,465633280,937492480,669057024,1205927936},
    {0,8388608,2139095040,1260388352,2106064896,2106064896,1287127040,
     1740111872,1815609344,933429248,1810038784,1338179584,197328896},
    {0,4194304,1077936128,1046478848,1800929280,1767374848,435683328,
     1341652992,1484259328,1869021184,1407647744,461832192,327614464},
    {0,2097152,1616904192,2127036416,1152909312,1387790336,863010816,
     393248768,946896896,1592721408,669974528,858718208,655294464},
    {0,1048576,1347420160,1572339712,456196096,661716992,1851883520,
     644448256,2020179968,973012992,1275002880,1652490240,1243545600},
    {0,524288,2021130240,827981824,639827968,1286275072,1422622720,
     1154711552,1893433344,2011039744,333383680,1090455552,339675136},
    {0,262144,1145307136,1614675968,1548156928,425656320,697794560,
     162496512,1774157824,1807554560,597563392,100603904,679417856},
    {0,131072,1717960704,1216249856,952508416,817766400,1537673216,
     458721280,1543473152,1465595904,1125922816,201209856,1426012160},
    {0,655336,1431633920,945651712,1908695040,1903452160,1069946880,
     1043306496,1073190912,714504192,102332416,402487296,771651584},
    {0,32768,2147450880,2050039808,1669980160,1655300096,2139371520,
     2094512128,2137503744,1359804416,137558016,804976640,1541273600},
    {0,16384,1073758208,1583374336,1192936448,1461445632,2004908032,
     1907357696,1984428032,574156864,342224960,1542844480,865859648},
    {0,8192,1610637312,919095296,511552512,780127232,1736994944,
     1809315968,1812428928,1081139392,686547136,938205376,1664612544},
    {0,4096,1342197760,1706571776,754088960,1320655872,1201168768,
     1471148416,1476817280,79806912,1442300352,1811333568,1114634688},
    {0,2048,2013296640,1264220672,1538054272,522598528,255345536,
     794291072,939811712,161711040,802130880,1473022912,83819456},
    {0,1024,1140868096,1044274688,924431744,1044738432,376967040,
     1580193664,2013284224,325519296,1604198336,731389888,234684352},
    {0,512,1711302144,2130622080,2088397696,2093148032,611882760,
     878683912,1887440776,651102082,993804162,1393639298,536412034},
    {0,256,1426085120,1574823296,1794429712,1774512912,1215931928,
     1614267928,1770004376,1235158790,1987673862,570654470,1005715270},
    {0,128,2139127680,823225120,1168671280,1401549360,410242232,
     1223658552,1535129528,389942798,1825766990,1143404110,1944254158},
    {0,64,1077952576,1610636896,458752112,656139376,812620280,
     299341944,1065379832,777854046,1501951134,208530654,1743054302},
    {0,32,1616928864,1207973576,649592930,1283443810,1758968688,
     590295160,2131249016,1488664830,856416638,484104702,1340654590},
    {0,16,1347440720,939550136,1563492422,414686294,1369962081,
     1180590193,1981288049,827816380,1645658814,966114238,468813820},
    {0,8,2021161080,2046839642,941817886,824588462,726658251,
     205277419,1815616739,1722809210,1143770494,1997176636,935528440},
    {0,4,1145324612,1577074238,1879506988,1879405006,1578619287,
     410522071,1484272071,1298134710,209263358,1846933048,1871120370},
    {0,2,1717986918,914390782,1644568666,1644357534,1009722151,
     947430191,946922383,513795374,483606012,1546380400,1596917668},
    {0,1,1431655765,1702911453,1174691895,1443162927,2019444294,
     2020722270,2020198302,1027523167,1032223673,1014481121,1048512329}};


  /* First call to the sequence */
  if (rng->counter == 1) 
    {

      /* Initialization of initX_n[] */
      for (i=1; i<=DIM_MAX_NIED; i++) 
        initialX_n[i]= 0;

      facteur= 1.0/(double)(1UL << (BIT_MAX_NIED+1));
      saut = 1L << 12;
      initial_d = saut;

      /* Gray code of saut */
      gray = saut^(saut >> 1); 
      for (i=0; i<=BIT_MAX_NIED; i++) 
        {
          if (gray == 0) 
            break;
          for (j= 1; j<= DIM_MAX_NIED; j++) 
            {
              if ((gray & 1) == 0) 
                break;
              /* XOR sum */
              initialX_n[j] ^= C[i][j];

            }
          gray >>= 1;
        }
    }

  /* Calculation of a new quasi-random vector on each call */
  dim= initial_d++;

  /* Research of the rightmost 0 bit */
  for (i=0; i<=BIT_MAX_NIED; i++) 
    {
      if ((dim & 1) == 0) 
        break;
      dim >>= 1;
    }
  /* Computation of the term n from the term n-1 */
  for (j=1; j<= rng->dimension; j++) 
    {
      X_n[j-1]= (double)initialX_n[j]*facteur;
      initialX_n[j] ^= C[i][j];
    }
  rng->counter++;
}


/*
 * Interface for Random Generators
 */

static knuth_state knuth_st;
static mrgk3_state mrgk3_st;
static mrgk5_state mrgk5_st;
static shufl_state shufl_st;
static mt_state mt_st1;
static mt_state mt_st2;

PnlRng PnlRngKnuth = 
  {
    {PNL_TYPE_RNG,pnl_rng_label,PNL_TYPE_RNG, (destroy_func *) pnl_rng_free},
    PNL_RNG_KNUTH,&KNUTH,
    MC,0, 0,0,0,sizeof(knuth_state), &knuth_st
  };
PnlRng PnlRngMrgk3 = 
  {
    {PNL_TYPE_RNG,pnl_rng_label,PNL_TYPE_RNG, (destroy_func *) pnl_rng_free},
    PNL_RNG_MRGK3,&MRGK3,
    MC,0,0,0,0,sizeof(mrgk3_state),&mrgk3_st
  };
PnlRng PnlRngMrgk5 = 
  {
    {PNL_TYPE_RNG,pnl_rng_label,PNL_TYPE_RNG, (destroy_func *) pnl_rng_free},
    PNL_RNG_MRGK5,&MRGK5,
    MC,0,0,0,0,sizeof(mrgk5_state),&mrgk5_st
  };
PnlRng PnlRngShufl = 
  {
    {PNL_TYPE_RNG,pnl_rng_label,PNL_TYPE_RNG, (destroy_func *) pnl_rng_free},
    PNL_RNG_SHUFL,&SHUFL,
    MC,0,0,0,0,sizeof(shufl_state),&shufl_st
  };
PnlRng PnlRngLecuyer =
  {
    {PNL_TYPE_RNG,pnl_rng_label,PNL_TYPE_RNG, (destroy_func *) pnl_rng_free},
    PNL_RNG_LECUYER,&LECUYER,
    MC,0, 0,0,0,0,NULL
  };
PnlRng PnlRngTausworthe = 
  {
    {PNL_TYPE_RNG,pnl_rng_label,PNL_TYPE_RNG, (destroy_func *) pnl_rng_free},
    PNL_RNG_TAUSWORTHE,&TAUS,
    MC,0, 0,0,0,0,NULL
  };
PnlRng PnlRngMersenne = 
  {
    {PNL_TYPE_RNG,pnl_rng_label,PNL_TYPE_RNG, (destroy_func *) pnl_rng_free},
    PNL_RNG_MERSENNE,&MERSENNE,
    MC,0,0, 0,0,sizeof(mt_state),&mt_st1
  };
PnlRng PnlRngMersenneRandomSeed = 
  {
    {PNL_TYPE_RNG,pnl_rng_label,PNL_TYPE_RNG, (destroy_func *) pnl_rng_free},
    PNL_RNG_MERSENNE_RANDOM_SEED,&MERSENNE,
    MC,0,0, 0,0,sizeof(mt_state),&mt_st2
  };
PnlRng PnlRngSqrt = 
  {
    {PNL_TYPE_RNG,pnl_rng_label,PNL_TYPE_RNG, (destroy_func *) pnl_rng_free},
    PNL_RNG_SQRT,&SQRT,
    QMC,0, 0,0,0,0,NULL
  };
PnlRng PnlRngHalton = 
  {
    {PNL_TYPE_RNG,pnl_rng_label,PNL_TYPE_RNG, (destroy_func *) pnl_rng_free},
    PNL_RNG_HALTON,&HALTON,
    QMC,0, 0,0,0,0,NULL
  };
PnlRng PnlRngFaure = 
  {
    {PNL_TYPE_RNG,pnl_rng_label,PNL_TYPE_RNG, (destroy_func *) pnl_rng_free},
    PNL_RNG_FAURE,&FAURE,
    QMC,0, 0,0,0,0,NULL
  };
PnlRng PnlRngSobol = 
  {
    {PNL_TYPE_RNG,pnl_rng_label,PNL_TYPE_RNG, (destroy_func *) pnl_rng_free},
    PNL_RNG_SOBOL,&SOBOL,
    QMC,0, 0,0,0,0,NULL
  };
PnlRng PnlRngSobol2 = 
  {
    {PNL_TYPE_RNG,pnl_rng_label,PNL_TYPE_RNG, (destroy_func *) pnl_rng_free},
    PNL_RNG_SOBOL2,&SOBOL2,
    QMC,0, 0,0,0,0,NULL
  };
PnlRng PnlRngNiederreiter = 
  {
    {PNL_TYPE_RNG,pnl_rng_label,PNL_TYPE_RNG, (destroy_func *) pnl_rng_free},
    PNL_RNG_NIEDERREITER,&NIEDERREITER,
    QMC,0, 0,0,0,0,NULL
  };

/*
 * Random Number Generator Array
 */
PnlRngTypes PnlRngArray[]= 
  {
    {{"KNUTH", PNL_RNG_KNUTH},&PnlRngKnuth},
    {{"MRGK3", PNL_RNG_MRGK3},&PnlRngMrgk3},
    {{"MRGK5", PNL_RNG_MRGK5},&PnlRngMrgk5},
    {{"SHUFL", PNL_RNG_SHUFL},&PnlRngShufl},
    {{"L'ECUYER", PNL_RNG_LECUYER},&PnlRngLecuyer},
    {{"TAUSWORTHE", PNL_RNG_TAUSWORTHE},&PnlRngTausworthe},
    {{"MERSENNE", PNL_RNG_MERSENNE},&PnlRngMersenne},
    {{"MERSENNE (Random Seed)", PNL_RNG_MERSENNE_RANDOM_SEED},&PnlRngMersenneRandomSeed},
    {{"SQRT", PNL_RNG_SQRT},&PnlRngSqrt},
    {{"HALTON", PNL_RNG_HALTON},&PnlRngHalton},
    {{"FAURE", PNL_RNG_FAURE},&PnlRngFaure},
    {{"SOBOL", PNL_RNG_SOBOL},&PnlRngSobol},
    {{"SOBOL2", PNL_RNG_SOBOL2},&PnlRngSobol2},
    {{"NIEDERREITER", PNL_RNG_NIEDERREITER},&PnlRngNiederreiter},
    {{NULL, NULLINT}, NULL}
  };

/*
 * True MC generators do not take into account the parameter dimension in the
 * Compute function.
 */
PnlRngTypes PnlRngMCArray[]= 
  {
    {{"KNUTH", PNL_RNG_KNUTH},&PnlRngKnuth},
    {{"MRGK3", PNL_RNG_MRGK3},&PnlRngMrgk3},
    {{"MRGK5", PNL_RNG_MRGK5},&PnlRngMrgk5},
    {{"SHUFL", PNL_RNG_SHUFL},&PnlRngShufl},
    {{"L'ECUYER", PNL_RNG_LECUYER},&PnlRngLecuyer},
    {{"TAUSWORTHE", PNL_RNG_TAUSWORTHE},&PnlRngTausworthe},
    {{"MERSENNE", PNL_RNG_MERSENNE},&PnlRngMersenne},
    {{"MERSENNE (Random Seed)", PNL_RNG_MERSENNE_RANDOM_SEED},&PnlRngMersenneRandomSeed},
    {{NULL, NULLINT}, NULL}
  };

enum_members RNGs = { sizeof(PnlRngArray[0]), (enum_member*)&PnlRngArray[0], sizeof(PnlRngArray)/sizeof(PnlRngTypes) };
enum_members MC_RNGs = { sizeof(PnlRngMCArray[0]), (enum_member*)&PnlRngMCArray[0], sizeof(PnlRngMCArray)/sizeof(PnlRngTypes) };


PnlRng* pnl_rng_get_from_id (int id)
{
  return PnlRngArray[id].rng;
}

/**
 * Initialises a generator
 * @param type_generator index of the generator to be used
 * @param dimension dimension of the value space to simulate in. Only
 * used for QMC.
 * @param samples maximum number of samples requested. Only used for QMC.
 * @returns OK or WRONG
 */
int pnl_rand_init (int type_generator, int dimension, long samples)
{
  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator); 

  switch (rng->type)
    {
      /*
       * some MC generators which must be initialized
       */
    case PNL_RNG_KNUTH:
      pnl_rand_sseed (type_generator, 161803398L);
      break;
    case PNL_RNG_MRGK3:
      pnl_rand_sseed (type_generator, 0);
      break;
    case PNL_RNG_MRGK5:
      pnl_rand_sseed (type_generator, 0);
      break;
    case PNL_RNG_SHUFL:
      pnl_rand_sseed (type_generator, 1043618065);
      break;
    case PNL_RNG_MERSENNE:
      pnl_rand_sseed (type_generator, 0);
      break;
    case PNL_RNG_MERSENNE_RANDOM_SEED:
      pnl_rand_sseed (type_generator, time(NULL));
      break;
    case PNL_RNG_DCMT:
      pnl_dcmt_sseed ((dcmt_state *)(rng->state), 1234);
      break;
      /*
       * Check if dimension > max_dim for QMC
       */
    case PNL_RNG_FAURE : 
      if (dimension > DIM_MAX_FAURE || samples>MAX_SAMPLE_FAURE) return FAIL;
      binomial(MAXI);
      break;
    case PNL_RNG_SQRT:
    case PNL_RNG_HALTON:
      if ( dimension  > DIM_MAX_QMC ) return FAIL;
      break;
    case PNL_RNG_SOBOL:
      if ( dimension  > DIM_MAX_SOBOL ) return FAIL;
      break;
    case PNL_RNG_NIEDERREITER:
      if ( dimension  > DIM_MAX_NIED ) return FAIL;
      break;
    default:
      break;
    }

  rng->dimension=dimension;
  rng->counter=1;
  rng->has_gauss=0;
  rng->gauss=0.;
  return OK;
}

void pnl_rand_sseed (int type_generator, ulong seed)
{
  PnlRng *rng;

  rng = pnl_rng_get_from_id(type_generator);
  pnl_rng_sseed (rng, seed);
}


/**
 * Initialises a rng of the given type.
 * Note that the fields size_State and state are set to zero, which implies
 * that the created generator is unusable. It is only usefull to receive an
 * already workin generator.
 *
 * @param rng is a rng returned by pnl_rng_new
 * @param type the type of generator to create
 */
void pnl_rng_init (PnlRng *rng, int type)
{

  rng->type = type;
  rng->dimension = 0;
  rng->counter = 0;
  rng->has_gauss = 0;
  rng->gauss = 0;
  rng->size_state = 0;
  rng->state = NULL;
  switch (type)
    {
    case PNL_RNG_KNUTH:
      rng->Compute = KNUTH;
      rng->rand_or_quasi = MC;
      break;
    case PNL_RNG_MRGK3:
      rng->Compute = MRGK3;
      rng->rand_or_quasi = MC;
      break;
    case PNL_RNG_MRGK5:
      rng->Compute = MRGK5;
      rng->rand_or_quasi = MC;
      break;
    case PNL_RNG_SHUFL:
      rng->Compute = SHUFL;
      rng->rand_or_quasi = MC;
      break;
    case PNL_RNG_LECUYER:
      rng->Compute = LECUYER;
      rng->rand_or_quasi = MC;
      break;
    case PNL_RNG_TAUSWORTHE:
      rng->Compute = TAUS;
      rng->rand_or_quasi = MC;
      break;
    case PNL_RNG_MERSENNE:
    case PNL_RNG_MERSENNE_RANDOM_SEED:
      rng->Compute = MERSENNE;
      rng->rand_or_quasi = MC;
      break;
    case PNL_RNG_DCMT:
      rng->Compute = DYNAMIC_MT;
      rng->rand_or_quasi = MC;
      break;
    case PNL_RNG_SQRT:
      rng->Compute = SQRT;
      rng->rand_or_quasi = QMC;
      break;
    case PNL_RNG_HALTON:
      rng->Compute = HALTON;
      rng->rand_or_quasi = QMC;
      break;
    case PNL_RNG_FAURE:
      rng->Compute = FAURE;
      rng->rand_or_quasi = QMC;
      break;
    case PNL_RNG_SOBOL:
      rng->Compute = SOBOL;
      rng->rand_or_quasi = QMC;
      break;
    case PNL_RNG_SOBOL2:
      rng->Compute = SOBOL2;
      rng->rand_or_quasi = QMC;
      break;
    case PNL_RNG_NIEDERREITER:
      rng->Compute = NIEDERREITER;
      rng->rand_or_quasi = QMC;
      break;
    }
}


/**
 * Determines the type of a generator
 * @param type_generator index of the generator
 * @returns MC or QMC
 */
int pnl_rand_or_quasi (int type_generator)
{
  return pnl_rng_get_from_id(type_generator)->rand_or_quasi;
}


/**
 * Frees an rng
 *
 * @param rng the address of an rng
 */
void pnl_rng_free (PnlRng **rng)
{
  if ( *rng == NULL ) return;
  if ( (*rng)->state != NULL )
    {
      free ((*rng)->state); (*rng)->state = NULL;
    }
  free (*rng); *rng = NULL;
}

/**
 * Creates an empty rng
 */
PnlRng* pnl_rng_new ()
{
  PnlRng *rng;
  if ((rng = malloc (sizeof(PnlRng))) == NULL) return NULL;
  
  rng->object.type = PNL_TYPE_RNG;
  rng->object.label = pnl_rng_label;
  rng->object.parent_type = PNL_TYPE_RNG;
  rng->object.destroy = (destroy_func *) pnl_rng_free;

  rng->type = PNL_RNG_NULL;
  rng->Compute = NULL;
  rng->rand_or_quasi = MC;
  rng->dimension = 0;
  rng->counter = 0;
  rng->has_gauss = 0;
  rng->gauss = 0;
  rng->size_state = 0;
  rng->state = NULL;
  return rng;
}


/**
 * Creates a rng of the given type.
 * Note that the fields size_state and state are set to zero, which implies
 * that the created generator is unusable. It is only usefull to receive an
 * already workin generator.
 *
 * @param type the type of generator to create
 * @return a PnlRng or NULL if an error occurred
 */
PnlRng* pnl_rng_create (int type)
{
  PnlRng *rng;
  if ((rng = pnl_rng_new()) == NULL) return NULL;

  rng->type = type;
  rng->dimension = 0;
  rng->counter = 0;
  rng->has_gauss = 0;
  rng->gauss = 0;
  rng->size_state = 0;
  rng->state = NULL;
  switch (type)
    {
    case PNL_RNG_KNUTH:
      rng->Compute = KNUTH;
      rng->rand_or_quasi = MC;
      rng->size_state = sizeof(knuth_state);
      rng->state = malloc(rng->size_state);
      break;
    case PNL_RNG_MRGK3:
      rng->Compute = MRGK3;
      rng->rand_or_quasi = MC;
      rng->size_state = sizeof(mrgk3_state);
      rng->state = malloc(rng->size_state);
      break;
    case PNL_RNG_MRGK5:
      rng->Compute = MRGK5;
      rng->rand_or_quasi = MC;
      rng->size_state = sizeof(mrgk5_state);
      rng->state = malloc(rng->size_state);
      break;
    case PNL_RNG_SHUFL:
      rng->Compute = SHUFL;
      rng->rand_or_quasi = MC;
      rng->size_state = sizeof(shufl_state);
      rng->state = malloc(rng->size_state);
      break;
    case PNL_RNG_MERSENNE:
      rng->Compute = MERSENNE;
      rng->rand_or_quasi = MC;
      rng->size_state = sizeof(mt_state);
      rng->state = malloc(rng->size_state);
      break;
    case PNL_RNG_DCMT:
      rng->Compute = DYNAMIC_MT;
      rng->rand_or_quasi = MC;
      rng->size_state = sizeof(dcmt_state);
      rng->state = pnl_dcmt_create ();
      break;
    default:
      free (rng); rng = NULL;
      printf("Unknown generator type\n");
    }
  return rng;
}


PnlRng** pnl_rng_dcmt_create_array (int n, ulong seed, int *count)
{
  PnlRng **rngtab;
  dcmt_state **mts;
  int i;

  mts = pnl_dcmt_create_array (n, seed, count);
  if (n != *count)
    {
      perror ("Not all generators could be created\n");
      abort ();
    }
  if ( (rngtab = malloc (n * sizeof(PnlRng *))) == NULL)
    {
      pnl_dcmt_free_array (mts, n); return NULL;
    }
  for ( i=0 ; i<n ; i++ )
    {
      rngtab[i] = pnl_rng_new ();
      pnl_rng_init (rngtab[i], PNL_RNG_DCMT);
      rngtab[i]->size_state = sizeof (dcmt_state);
      rngtab[i]->state = mts[i];
    }
  free (mts);
  return rngtab;
}

static void pnl_knuth_sseed (knuth_state *s, ulong seed)
{
  s->SEED = seed;
}

static void pnl_mrgk3_sseed (mrgk3_state *s, ulong seed)
{
  /* Initialization  of the two generators */
  s->x10=231458761.;
  s->x11=34125679.;
  s->x12=45678213.;
  
  s->x20=57964412.;
  s->x21=12365487.;
  s->x22=77221456.;
}

static void pnl_mrgk5_sseed (mrgk5_state *s, ulong seed)
{
  /* Initialization  of the two generators */
  s->x10= 231458761.;
  s->x11= 34125679.;
  s->x12= 45678213.;
  s->x13= 7438902.;
  s->x14= 957345.;

  s->x20= 57964412.;
  s->x21= 12365487.;
  s->x22= 77221456.;
  s->x23= 816403.;
  s->x24= 8488912.;
}

static void pnl_shufl_sseed (shufl_state *s, ulong seed)
{
  s->y = 0;
  s->x= seed;
}

void pnl_rng_sseed (PnlRng *rng, ulong seed)
{
  switch (rng->type)
    {
    case PNL_RNG_KNUTH :
      pnl_knuth_sseed((knuth_state *)(rng->state), seed);
      break;
    case PNL_RNG_MRGK3 :
      pnl_mrgk3_sseed((mrgk3_state *)(rng->state), seed);
      break;
    case PNL_RNG_MRGK5 :
      pnl_mrgk5_sseed((mrgk5_state *)(rng->state), seed);
      break;
    case PNL_RNG_SHUFL :
      pnl_shufl_sseed((shufl_state *)(rng->state), seed);
      break;
    case PNL_RNG_MERSENNE :
    case PNL_RNG_MERSENNE_RANDOM_SEED :
      pnl_mt_sseed((mt_state *)(rng->state), seed);
      break;
    case PNL_RNG_DCMT :
      pnl_dcmt_sseed ((dcmt_state *)(rng->state), seed);
      break;
    }
  rng->counter=1;
  rng->has_gauss=0;
  rng->gauss=0.;
}
