/* This file contains two parts.  The first one gathers the random generator
   functions which must never be called directly. These codes are provided by
   Premia. The second part is an interface to these generators to easily
   create vectors and matrices of random numbers according to different
   laws. */



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


/* ---------------------------------- */
/* PSEUDO RANDOM NUMBERS GENERATORS. */
/* ---------------------------------- */

#include "pnl_mathtools.h"
#include "pnl_random.h"
#include "pnl_cdf.h"
#include "randomkit.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>

#define INV_M (1.0/M)
#define INV_M1 (1.0/M1)

static long counter= 0;
static int draw_new_sample = 0;
#define CheckMaxQMCDim(type_generator, dimension)                       \
  {   if (dimension >= pnl_Random[type_generator].Dimension)         \
      {                                                                 \
        perror("maximum dimension of Monte Carlo exceeded\n"); abort(); \
      }                                                                 \
  }

/********************************/
/******MONTE CARLO STANDARD******/
/********************************/

/* ------------------------------------------------------------------- */
/* Random numbers generator of Knuth : 
   It is  based on MRG and it uses a substractive method  */
/* ------------------------------------------------------------------- */

static void KNUTH(int dimension,double *sample)
{
  static long M= 1000000000;
  static long SEED= 161803398;

  /* Initialize the sequence with a positive seed */
  static long alea= 1; 
  static int inc1, inc2;
  static long t_alea[56];
  long X_n, y_k;
  int i, ii, l;
  
  /* First call to the sequence */
  if(counter == 1)
    {
      X_n= SEED- alea;
      X_n%= M;
      t_alea[55]= X_n;
      y_k= 1;
      /* Initialization of the table */
      for(i= 1; i<= 54; i++)
        {
          ii= (21*i)%55; /* 21 was chosen to alleviate initial
                            nonrandomness problems */
          t_alea[ii]= y_k;
          y_k= X_n - y_k;
          if(y_k < 0) 
            y_k+= M;
          X_n= t_alea[ii];
        }
      
      /* Randomization of the elements of the table */   
      for(l=  1; l<= 4; l++)
        {
          for(i= 1; i<= 55; i++)
            {
              t_alea[i]-= t_alea[1+(i+30)%55];
              if(t_alea[i] < 0) 
                t_alea[i]+= M;
            }
        }
      inc1= 0;
      inc2= 31;  /* 31 is a special value of Knuth : 31= 55-24 */
      alea= 1;
    }
  
  counter++;
  
  /* For each call to the sequence, computation of a new point */
  if(++inc1 == 56) 
    inc1= 1;
  if(++inc2 == 56) 
    inc2= 1;
  /* Substractive method*/
  X_n= t_alea[inc1] - t_alea[inc2];
  
  if(X_n < 0) 
    X_n+= M;
  t_alea[inc1]= X_n;
  /* Normalized value */
  *sample=X_n*INV_M;
  return;
}

/* ----------------------------------------------------------------------- */
/* Combination of two multiplicative recursive generators of order 3 (k=3) */
/* ----------------------------------------------------------------------- */
static void MRGK3(int dimension, double *sample)
{
  static double M1= 4294967087.0;
  static double M2= 4294944443.0;
  static double A12= 1403580.0;
  static double A13N= 810728.0;
  static double A21= 527612.0;
  static double A23N= 1370589.0;
  static double NORM= 2.328306549295728e-10;

  static double x10, x11, x12, x20, x21, x22;
  long k;
  double p1, p2;
  
  /* First call to the sequence */
  if (counter == 1 ) 
    {
      /* Initialization */
      x10=231458761.;
      x11=34125679.;
      x12=45678213.;
      x20=57964412.;
      x21=12365487.;
      x22=77221456.;
      counter+=1; 
    }
  
  /* For each call to the sequence, computation of a new point */
  /* First generator */
  p1= A12*x11 - A13N*x10;
  k= (long)floor(p1/M1); /*TOCHECK*/
  p1-= k*M1;

  if(p1 < 0.0) 
    p1+= M1;

  x10= x11;
  x11= x12;
  x12= p1;
  
  /* Second generator */
  p2= A21*x22 - A23N*x20;
  k= (long)floor(p2/M2);/*TOCHECK*/
  p2-= k*M2;

  if(p2 < 0.0) 
    p2+= M2;

  x20= x21;
  x21= x22;
  x22= p2;
  

  /* Combination of the two generators */
  if (p1< p2) 
    *sample= (p1- p2+ M1)*NORM;
  else 
    *sample=(p1- p2)*NORM;
  return;
}



/* ----------------------------------------------------------------------- */
/* Combination of two multiplicative recursive generators of order 5 (k=5) */
/* ----------------------------------------------------------------------- */
static void MRGK5(int dimension,double *sample)
{
  static double M1= 4294949027.0;
  static double M2=    4294934327.0;
  static double A12=   1154721.0;
  static double A14=   1739991.0;
  static double A15N=  1108499.0;
  static double A21=   1776413.0;
  static double A23=   865203.0;
  static double A25N= 1641052.0;
  static double NORM=  2.3283163396834613e-10;
  
  static double x10, x11, x12, x13, x14, x20, x21, x22, x23, x24;
  long k;
  double p1, p2;
  
  /* First call to the sequence */
  if (counter == 1) 
    {
      /*Initialization*/
      x10= 231458761.;
      x11= 34125679.;
      x12= 45678213.;
      x13= 7438902.;
      x14= 957345.;

      x20= 57964412.;
      x21= 12365487.;
      x22= 77221456.;
      x23= 816403.;
      x24= 8488912.;
      counter++;
    }
  
  /* For each call to the sequence, computation of a new point */
  /* First generator with Schrage method */
  p1= A12*x13 - A15N*x10;
  
  if(p1> 0.0) 
    p1-= A14*M1;

  p1+= A14*x11;
  k= (long)floor(p1/M1);/*TOCHECK*/
  p1-= k*M1;

  if(p1< 0.0) 
    p1+= M1;
  
  x10= x11;
  x11= x12;
  x12= x13; 
  x13= x14;
  x14= p1;
  
  /* Second generator with Schrage method */
  p2= A21*x24 - A25N*x20;

  if(p2> 0.0)
    p2-= A23*M2;

  p2+= A23*x22;
  k= (long)floor(p2/M2);/*TOCHECK*/
  p2-= k*M2;

  if(p2< 0.0)
    p2+= M2;

  x20= x21;
  x21= x22;
  x22= x23;
  x23= x24;
  x24= p2;
  

  /*Combination of the two generators */
  if (p1<= p2) 
    *sample= (p1- p2+ M1)*NORM;
  else 
    *sample=(p1- p2)*NORM;

  return;
}
/* ------------------------------------------------------------------------ */
/* Random numbers generator of  Park & Miller with Bayes & Durham shuffling
   procedure : the next random number is not obtained from the previous one 
   but we use an intermediate table which contains the 32 precedent random 
   numbers and we choose one of them randomly.  */
/* ----------------------------------------------------------------------- */
static void SHUFL(int dimension,double *sample)
{
  static long A= 16807;        /* multiplier */
  static long M= 2147483647L;    /* 2**31 - 1  */
  static long Q= 127773L;        /*   M div A  */  
  static long R= 2836;           /*   M mod A  */
  long N1;
  
  int j;
  long hi;                 /* high order bit */
  static long y= 0;
  static long t[32];       /* 32 refers to the size of a computer word */
  /* Initialisation */
  static long x;
  
  
  N1=(M/32);
  
  /* First call to the sequence */
  if (counter == 1)
    {
      x= 1043618065;
      
      /* After 8 "warm-ups", initialisation of the shuffle table */
      for (j= 39; j>= 0; j--)
        {
          hi= x/Q;
          /*Schrage's method to avoid overflows */
          x= A*(x- hi*Q)- R*hi; 
          if (x < 0) 
            x+= M;
          if (j< 32) 
            t[j]= x;
        }
      y= t[0];                    
    }
  counter++;
  
  
  /* For each call to the sequence, computation of a new point */
  hi= x/Q;
  x= A*(x-hi*Q)- R*hi;
  if (x < 0) 
    x+= M;

  /* Shuffling procedure of Bayes & Durham */
  /* Index j dependent on the last point */
  j= y/N1;
  /* Next point dependent on j */
  y= t[j];
  
  t[j]= x;
  *sample= INV_M*y; 
  
  return;
}



/* ------------------------------------------------------------------------- */
/* Random numbers generator of L'Ecuyer with Bayes & Durham shuffling
   procedure : 
   Combination of two short periods LCG to obtain a longer period generator.
   The period is the least common multiple of the 2 others.  */
/* ------------------------------------------------------------------------- */

static void ECUYER(int dimension,double *sample)
{ 
  static long A1= 40014;        /* multiplier of the 1st generator */
  static long A2= 40692;        /* multiplier of the 2nd generator */
  static long M1= 2147483647;   /* 2**31 - 1   */ 
  static long M2= 2147483399;   /* 2**31 - 249 */
  static long Q1= 53668;        /* m1 div a1   */  
  static long Q2= 52774;        /* m2 div a2   */ 
  static long R1= 12221;        /* m1 mod a1   */
  static long R2= 3791;         /* m2 mod a2   */
  long N1;

  static long x;
  static long y= 978543162;
  int j;
  long hi;               /* high order bit */
  static long z= 0;
  static long t[32];     /* 32 is the size of a computer word */

  N1= (M1/32);
  
  /* First call to the sequence */
  if (counter == 1)
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
  counter++;              
  
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
  
  /* Shuffling procedure of Bqyes & Durham */
  /* Index j dependent on the last point */
  j= z/N1;
  /* Next point dependent on j */
  z= t[j]- y; 
  t[j]= x;
  
  
  /* To avoid 0 value */
  if (z < 1) 
    z+= M1-1;
  
  *sample= INV_M1*z;
  
  return;
}


/* ------------------------- */
/* ------------------------ */
/* Tausworthe Algorithm   */
/* ------------------------ */

/* maximal number of combined generators */
#define TAUS_MAX 10


/* ---------------------------------------------------- */
/* Generation of a random bit
   Algorithm based on a prime polynomial : 
   here we choose x^18 + x^5 + x^2 + x + 1 .
   It is described in 'Numerical Recipes in C' page 296. */
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
    Combination  of J Tausworthe generators 
    u(n)[j]= u(n-r)[j] ^ u(n-k)[j],
    with parameters k, q, r, s and t.
    Generator :
    v= =(u[0] ^ u[1] ... ^ u{J-1])/2^32.    
     
    L= 32 length of a word.    */ 
/* ---------------------------------------------------------------- */
static void TAUS(int dimension,double *sample)
{
  int L= 32;
  int J= 3;
  
  int i;
  static unsigned long u[TAUS_MAX];
  static unsigned long c[TAUS_MAX];
  unsigned long b;
  static int k[TAUS_MAX], q[TAUS_MAX];
  static int s[TAUS_MAX], r[TAUS_MAX], t[TAUS_MAX];
  unsigned long v= 0;

  /* First call to the sequence. Initialisation */
  if (counter == 1) 
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
      for(i= 0; i< J; i++)
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
      
      counter++; 
    }
  
  /* For each call to the sequence, computation of a new point */
  for(i= 0; i< J; i++)
    { 
      /* Calculus of the next point for the J generators */
      /* 6 steps explained by L'Ecuyer */
      b=(u[i]<<q[i]) ^ u[i];         /* Steps 1 and 2 */
      b >>= t[i];                    /* Step 3 */
      u[i]= (u[i] & c[i]) << s[i];   /* Steps 4 et 5 */
      u[i] ^= b;                     /* Step 6 */
      /* Combination : XOR between the J generators */
      v^= u[i];
    }

  /* Normalization by 1/2^32 */
  *sample=v * 2.3283064365e-10;

  return;
}





/* ****************************** */
/* ****************************** */
/* *****QUASI MONTE CARLO   ***** */
/* ****************************** */
/* ----------------------------------------------------------------*/
/* QUASI RANDOM NUMBERS ; LOW DISCEPANCY SEQUENCES.
   SQRT, Van der Corput - Halton, Faure, Sobol, generalized Faure,
   Niederreiter */
/* ----------------------------------------------------------------*/




#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

#define DIM_MAX_QMC 300

#define MAXI 33
long Comb[MAXI][MAXI];/*Binomial Coefficients*/

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
/* Search the smallest element of a table greater than a given threshold.
   Used for the prime numbers */
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
  return;
}

/* ----------------------------------------------------------------*/
/* SQRT Sequence */
/* Computation of the next element for the dim-dimensional sequence. */
/* ----------------------------------------------------------------*/

static void SQRT(int dim, double X_m[])
{
  static int dimension=0;
  int i;
  static int prime[DIM_MAX_QMC];
  static double alpha[DIM_MAX_QMC];

  /* Verification of the dimension. It must not change without reinitializing */
  if(dimension != dim)
    counter= 1;

  /* First call : initialisation */
  if(counter == 1)
    {
      dimension= dim;
      prime_number(dim, prime);
      for(i=0; i<dim; i++)
        {
          alpha[i]= (int) sqrt((double) prime[i]);
        }
    }

  /* For each call to the sequence, computation of a new point */
  for(i= 0; i< dim; i++)
    X_m[i]= ((counter*alpha[i])-floor(counter*alpha[i]));

  counter++;

  return;
}  


/* ------------------------------------------------------------- */
/* Halton sequence. Bases p={p1, p2, ..., pd} of prime numbers ; 
   Computation of the next element of the sequence X_n 
   (index 0 to d-1) */
/* ------------------------------------------------------------- */
static void HALTON(int dim, double X_n[])
{
  static int dimension=0;
  int i;
  static int prime[DIM_MAX_QMC];
  int coeff;
  double y;
  int x, puissance;
  
  /* Verification of the dimension. It must not change without reinitializing */
  if(dimension != dim)
    counter= 1;

  /* First call : initialization */
  if( counter == 1)
    {
      dimension= dim;
      prime_number(dim, prime);
    }

  /* For each call to the sequence, computation of a new point */
  for(i = 0; i< dim; i++)
    {
      puissance= 1;
      x= 0; y= 0;
      /* radical inverse function */
      while ((counter-x)>0) 
        {
          coeff= ((counter-x)/puissance)%prime[i];
          x += coeff*puissance;
          puissance *= prime[i];
          y += coeff*1.0/puissance;
        }
      X_n[i]= y;
    }
  counter++;

  return;
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

static void FAURE(int d, double U_n[])
{
  int coeff[MAXI];
  int b[MAXI];
  static int dimension=0, r;
  int prime[DIM_MAX_QMC];
  int x= 0, puissance1= 1, puissance2;
  int indice, i, j, k, somm;

  
  /* Verification of the dimension. It must not change without reinitializing */
  if(dimension != d)
    counter= 1;
  
  /*First call to the sequence */
  if(counter == 1)
    {
      dimension=d;
      if((d == 2)||(d == 1))
        r= 3;
      else
        {
          prime_number(d, prime);
          r=search_value(d,prime,d);
        }
    }
  
  /* Initialization */
  for (i=0; i< d; i++) 
    {
      U_n[i]= 0.;
    }
  
  indice= 0;
  /* For each call to the sequence, computation of a new point */
  
  /* r-digit expansion of n --> first term of the sequence */
  while ((counter-x) > 0)
    
    {
      coeff[indice]= ((counter-x)/puissance1)%r;
      x += coeff[indice]*puissance1;
      puissance1 *= r;
      U_n[0] +=(double)coeff[indice]/(double)puissance1;
      indice +=1;
    }
  
  /* Other terms of the sequence */
  /* Successive transformations of the r-digit expansion.*/
  for (k=1; k< d; k++) 
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
  counter++;

  return;
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

static void SOBOL(int d, double X_n[])
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
  if(counter == 1)
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
  dim= counter;
  /* Research of the rightmost 0 bit */
  for (i=1; i<=BIT_MAX_SOBOL; i++) 
    {
      if ((dim & 1) == 0) 
        break;
      dim >>= 1;
    }
  /* Computation of the term n from the term (n-1) */
  for (j=1; j<= d; j++) 
    {
      initialX_n[j] = (initialX_n[j])^ (C[j][i]<< (BIT_MAX_SOBOL-i));
      /* normalization */
      X_n[j-1]= initialX_n[j] *facteur;
    }
  counter++;
  
  return;
}


static void SOBOL2(int d, double X_n[])
{
  static rk_sobol_state s;
  static int dimension=0;

  if(counter == 1 || dimension != d)
    {
      if (dimension)
        rk_sobol_free(&s); /* We should free on exit too. */
      dimension=d;
      rk_sobol_init(d, &s, NULL, rk_sobol_Ldirections, NULL);
      /* For randomized QMC, add: rk_sobol_randomshift(&s, NULL); */
    }

  rk_sobol_double(&s, X_n);
  counter ++;
  return;
}


static void MERSENNE(int dimension,double *sample)
{
  static rk_state s;

  if(counter == 1) rk_seed(1, &s); /* Seed is fixed */
  /* For a random seed, use instead: rk_randomseed(&s);  */

  *sample = rk_double(&s);
  counter ++;
  return;
}

static void MERSENNE_RANDOMSEED(int dimension,double *sample)
{
  static rk_state s;

  if(counter == 1)  rk_randomseed(&s);
  
  *sample = rk_double(&s);
  counter ++;
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

static void NIEDERREITER(int d, double X_n[])
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
  if (counter == 1) 
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
  for (j=1; j<= d; j++) 
    {
      X_n[j-1]= (double)initialX_n[j]*facteur;
      initialX_n[j] ^= C[i][j];
    }
  counter++;
}


/* Maximum dimension for random sequences */
#define DIM_MAX 100000
static double ArrayOfRandomNumbers[DIM_MAX];

/*Random Number Generator Array*/
random_generator pnl_Random[]= 
  {
    {{"KNUTH", PNL_RNG_KNUTH},&KNUTH,MC,DIM_MAX},
    {{"MRGK3", PNL_RNG_MRGK3},&MRGK3,MC,DIM_MAX},
    {{"MRGK5", PNL_RNG_MRGK5},&MRGK5,MC,DIM_MAX},
    {{"SHUFL", PNL_RNG_SHUFL},&SHUFL,MC,DIM_MAX},
    {{"L'ECUYER", PNL_RNG_L_ECUYER},&ECUYER,MC,DIM_MAX},
    {{"TAUSWORTHE", PNL_RNG_TAUSWORTHE},&TAUS,MC,DIM_MAX},
    {{"MERSENNE", PNL_RNG_MERSENNE},&MERSENNE,MC,DIM_MAX},
    {{"MERSENNE (Random Seed)", PNL_RNG_MERSENNE_RANDOM_SEED},&MERSENNE_RANDOMSEED,MC,DIM_MAX},
    {{"SQRT", PNL_RNG_SQRT},&SQRT,QMC,DIM_MAX_QMC},
    {{"HALTON", PNL_RNG_HALTON},&HALTON,QMC,DIM_MAX_QMC},
    {{"FAURE", PNL_RNG_FAURE},&FAURE,QMC,DIM_MAX_FAURE},
    {{"SOBOL", PNL_RNG_SOBOL},&SOBOL,QMC,DIM_MAX_SOBOL},
    {{"SOBOL2", PNL_RNG_SOBOL2},&SOBOL2,QMC,DIM_MAX}, 
    {{"NIEDERREITER", PNL_RNG_NIEDERREITER},&NIEDERREITER, QMC,DIM_MAX_NIED},
    {{NULL, NULLINT}, NULL, NULLINT, NULLINT}
  };

/*
 * True MC generators do not take into account the parameter dimension in the
 * Compute function.
 */
random_generator pnl_Random_MC[]= 
  {
    {{"KNUTH", PNL_RNG_KNUTH},&KNUTH,MC,DIM_MAX},
    {{"MRGK3", PNL_RNG_MRGK3},&MRGK3,MC,DIM_MAX},
    {{"MRGK5", PNL_RNG_MRGK5},&MRGK5,MC,DIM_MAX},
    {{"SHUFL", PNL_RNG_SHUFL},&SHUFL,MC,DIM_MAX},
    {{"L'ECUYER", PNL_RNG_L_ECUYER},&ECUYER,MC,DIM_MAX},
    {{"TAUSWORTHE", PNL_RNG_TAUSWORTHE},&TAUS,MC,DIM_MAX},
    {{"MERSENNE", PNL_RNG_MERSENNE},&MERSENNE,MC,DIM_MAX},
    {{"MERSENNE (Random Seed)", PNL_RNG_MERSENNE_RANDOM_SEED},&MERSENNE_RANDOMSEED,MC,DIM_MAX},
    {{NULL, NULLINT}, NULL, NULLINT, NULLINT}
  };


/*The number of generators is kept in the flag GEN_NUMBER in optype.h*/
enum_members RNGs = { sizeof(pnl_Random[0]), (enum_member*)&pnl_Random[0], 13 };
enum_members MC_RNGs = { sizeof(pnl_Random_MC[0]), (enum_member*)&pnl_Random_MC[0], 7 };


/******************************************************
 *  Interface for calling random generators           *
 ******************************************************/


/** 
 * Initialise a generator
 * @param type_generator index of the generator to be used
 * @param simulation_dim dimension of the value space to simulate in. Only
 * used for QMC.
 * @param samples maximum number of samples requested. Only used for QMC.
 * @returns OK or WRONG
 */
int pnl_rand_init (int type_generator, int simulation_dim, long samples)
{
  /*For the Faure sequences we need to compute binomial coefficients
    wich are storable with a Long Type only up to the number 32.
    The corresponding sample value is MAX_SAMPLE_FAURE*/
  if((pnl_Random[type_generator].Compute)==&FAURE)
    {
      if (samples>MAX_SAMPLE_FAURE) return FAIL;
      binomial(MAXI);
    }
  if (pnl_rand_or_quasi (type_generator)==QMC &&  
      simulation_dim >= pnl_Random[type_generator].Dimension)
    return FAIL; 

  counter=1;
  draw_new_sample=1;
  return OK;
}

/** 
 * Determines the type of a generator
 * @param type_generator index of the generator
 * @returns MC or QMC
 */
int pnl_rand_or_quasi (int type_generator)
{
  return pnl_Random[type_generator].RandOrQuasi;
}

/** 
 * Returns the name of a generator
 * @param type_generator index of the generator
 * @returns a pointer to the name
 */
const char * pnl_rand_name (int type_generator)
{
  return pnl_Random[type_generator].base.label;
}


/**
 * Simulates a standard random normal variable using Box Muller's algorithm
 * @param type_generator index of the generator to be used
 * @return a normal random variable
 */    
static double Gauss_BoxMuller(int type_generator)
{
  double xs,ys, g1, g2;
  static double random_number;

  /* do not wast any samples. But be sure to throw away any remaining
     samples when pnl_rand_init is called */
  if (draw_new_sample==1)
    {
      /* draw 2 new samples */
      pnl_Random[type_generator].Compute(1,&xs);
      pnl_Random[type_generator].Compute(1,&ys);
      g1 = sqrt(-2.0*log(xs))*cos(2.0*M_PI*ys);
      g2 = sqrt(-2.0*log(xs))*sin(2.0*M_PI*ys);
      random_number = g2;
      draw_new_sample=0;
      return g1;
    }
  else
    {
      /* use the remaining sample from the last call */
      draw_new_sample=1;
      return random_number;
    }
}


/**
 * Simulation of a Gaussian standard variable.
 *
 * WARNING : A new random variable is drawn at each call.
 *
 * @param dimension size of the vector to simulate
 * @param create_or_retrieve boolean can be CREATE or RETRIEVE. UNUSED
 * @param index UNUSED
 * @param type_generator index of the generator
 */
static double GaussMC(int dimension, int create_or_retrieve, int index, int type_generator)
{
  if (create_or_retrieve == CREATE)
    {
      int i;
      double *ptr = ArrayOfRandomNumbers;
      for (i=0; i<dimension; i++, ptr++)
        *ptr = Gauss_BoxMuller(type_generator);
    }
  return (ArrayOfRandomNumbers[index]);
}


/**
 * Simulation of a Gaussian standard variable for Quasi Monte Carlo Simulation,
 * that is with the pnl_inv_cdfnor function. 
 *  This function can be called for the generation of a n-dimensional vector of 
 *  independent variables: call to a n-dimensional low-discrepancy sequence.
 * @param dimension size of the vector to simulate
 * @param create_or_retrieve boolean can be CREATE or
 * RETRIEVE. if it is CREATE, draw all the dimensions and returns the fisrt
 * one. If it s RETRIEVE, returns the dimension corresponding to index
 * @param index index to be returned
 * @param type_generator index of the generator
 */
static double GaussQMC(int dimension, int create_or_retrieve, int index, int type_generator)
{
  if (dimension >= pnl_Random[type_generator].Dimension)
    {
      perror("maximum dimension of Monte Carlo exceeded\n"); abort();
    }

  if (create_or_retrieve == CREATE)
    pnl_Random[type_generator].Compute(dimension,ArrayOfRandomNumbers);
  return pnl_inv_cdfnor(ArrayOfRandomNumbers[index]);
}


/**
 * Simulation of a Gaussian standard variable for Quasi Monte Carlo Simulation,
 * that is with the pnl_inv_cdfnor function. 
 *  This function can be called for the generation of a n-dimensional vector of 
 *  independent variables: call to a n-dimensional low-discrepancy sequence.
 * @param dimension size of the vector to simulate
 * @param create_or_retrieve boolean can be CREATE or
 * RETRIEVE. if it is CREATE, draw all the dimensions and returns the fisrt
 * one. If it s RETRIEVE, returns the dimension corresponding to index
 * @param index index to be returned
 * @param type_generator index of the generator
 */
double pnl_rand_gauss(int dimension, int create_or_retrieve, int index, int type_generator)
{
  if (pnl_rand_or_quasi (type_generator) == QMC)
    return GaussQMC(dimension, create_or_retrieve, index, type_generator);
  else
    return GaussMC(dimension, create_or_retrieve, index, type_generator);
            
}
/**
 * Simulation of a Bernoulli random variable
 * @param p parameter of the law
 * @param type_generator index of the generator ot be used
 */
int pnl_rand_bernoulli(double p, int type_generator)
{
  double x=0.0;
  pnl_Random[type_generator].Compute(1,&x);
  if (x<p) return 1; else return  0;
}

/**
 * Simulation of a Poisson random variable 
 * @param lambda parameter of the law
 * @param type_generator index of the generator ot be used
 */
long pnl_rand_poisson(double lambda, int type_generator)
{
  double u;
  double a = exp(-lambda);
  long n = 0;
  static double random_number;

  pnl_Random[type_generator].Compute(1,&random_number);
  u = random_number;
  while (u>a){
    pnl_Random[type_generator].Compute(1,&random_number); 
    u *= random_number;
    n++;
  }
  return n;
}

/**
 * Simulation of an exponential random variable
 * @param lambda parameter of the law
 * @param type_generator index of the generator ot be used
 */
double pnl_rand_exp(double lambda,int type_generator)
{
  double x;
  static double random_number;

  do{
    pnl_Random[type_generator].Compute(1,&random_number); 
    x=random_number;
  } while(x==0);
  return (double) (-log(x)/lambda);
}


/**
 * Simulation of a Poisson process 
 * @param lambda parameter of the law
 * @param t time of the simulation
 * @param type_generator index of the generator ot be used
 */
long pnl_rand_poisson1(double lambda, double t, int type_generator)
{
  double S;
  long Nt;
  Nt=0;
  S=0;
  do {
    S=S+pnl_rand_exp (lambda,type_generator);
    Nt=Nt+1;
  } while (S<=t);
  return Nt-1;
}


/**
 * Generate a uniformly distributed number on ]0,1).
 * @param type_generator index ot the generator to be used
 *
 * @see pnl_rand_uni_ab
 */
double pnl_rand_uni (int type_generator)
{
  double u=0.;
  while (u == 0) pnl_Random[type_generator].Compute(1,&u);
  return u;
}

/**
 * Generate a uniformly distributed number on [a,b].
 * @param a lower bound
 * @param b upper bound
 * @param type_generator index ot the generator to be used
 *
 * @see pnl_rand_uni
 */
double pnl_rand_uni_ab (double a, double b, int type_generator)
{
  double u=0.0;
  pnl_Random[type_generator].Compute(1,&u);
  return a+(b-a)*u;
}


/**
 * Generate a normally distributed number.
 * @param type_generator index ot the generator to be used
 */
double pnl_rand_normal (int type_generator)
{
  int mc_or_qmc = pnl_rand_or_quasi (type_generator);
  if (mc_or_qmc == QMC)
    {
      double u;
      pnl_Random[type_generator].Compute(1,&u);
      return pnl_inv_cdfnor(u);
    }
  return Gauss_BoxMuller(type_generator);
}

/**
 * return a vector of uniformly distributed components on [a,b]
 * @param G existing PnlVect containing the random numbers on exit
 * @param samples size of G (number of independent samples requested)
 * @param a lower bound 
 * @param b upper bound 
 * @param type_generator index ot the generator to be used
 *
 * @see pnl_vect_rand_uni_d
 */
void pnl_vect_rand_uni(PnlVect *G, int samples, double a, double b, int type_generator)
{
  int i;
  double random_number;
  double *ptr;
  pnl_vect_resize(G,samples);
  ptr=G->array;
  for(i=0;i<samples;i++)
    {
      pnl_Random[type_generator].Compute(1,&random_number);
      *ptr=a+(b-a)*random_number;
      ptr++;
    }
}


/**
 * return a vector uniformly distributed on [a,b]^dimension
 *
 * if the generator is a true MC generator, no difference between this
 * function and pnl_vect_rand_uni. In case of a QMC generator, this
 * function generator one sample with values in [a,b]^dimension and NOT
 * dimension samples with values in [a, b].
 *
 * @param G existing PnlVect containing the random numbers on exit
 * @param dimension dimension of the state space
 * @param a lower bound 
 * @param b upper bound 
 * @param type_generator index ot the generator to be used
 *
 * @see pnl_vect_rand_uni
 */
void pnl_vect_rand_uni_d (PnlVect *G, int dimension, double a, double b, int type_generator)
{
  int i;
  double *ptr;
  pnl_vect_resize(G,dimension);
  ptr=G->array;
  if (pnl_rand_or_quasi (type_generator) == QMC)
    {
      CheckMaxQMCDim(type_generator, dimension);
      pnl_Random[type_generator].Compute(dimension, ptr);
      for(i=0;i<dimension;i++)
        {
          *ptr=a+(b-a)* (*ptr);  ptr++;
        }
      return;
    }
  for(i=0;i<dimension;i++)
    {
      pnl_Random[type_generator].Compute(1, ptr);
      *ptr=a+(b-a)* (*ptr);  ptr++;
    }
}



/**
 * return a vector of normaly distributed components on R
 *
 * @param samples number of samples
 * @param G : the vector of gaussian numbers, must already be allocated.
 * @param type_generator : the index of the generator to be used 
 *
 * @see pnl_vect_rand_normal_d
 */
void pnl_vect_rand_normal (PnlVect *G, int samples, int type_generator)
{
  int i;
  double *ptr;
  pnl_vect_resize(G,samples);
  ptr = G->array;
  if ( pnl_rand_or_quasi(type_generator) == QMC)
    {
      for (i=0; i<samples; i++)
        {
          pnl_Random[type_generator].Compute(1,ptr);
          *ptr=pnl_inv_cdfnor(*ptr);
          ptr++;
        }
      return;
    }
  for (i=0; i<samples; i++)
    {
      *ptr=Gauss_BoxMuller(type_generator);
      ptr++;
    }
}


/**
 * return a vector normally distributed on R^dimension.
 *
 * if the generator is a true MC generator, no difference between this
 * function and pnl_vect_rand_uni. In case of a QMC generator, this
 * function generator one sample with values in R^dimension and NOT
 * dimension samples with values in R.
 *
 * @param dimension : size of the vector. one sample of a Gaussian vector.
 * @param G : the vector of gaussian numbers, must already be allocated.
 * @param type_generator : the index of the generator to be used 
 *
 * @see pnl_vect_rand_normal
 */
void pnl_vect_rand_normal_d (PnlVect *G, int dimension, int type_generator)
{
  int i;
  double *ptr;
  pnl_vect_resize(G,dimension);
  ptr = G->array;
  if (pnl_rand_or_quasi(type_generator) == QMC)
    {
      CheckMaxQMCDim(type_generator, dimension);
      pnl_Random[type_generator].Compute(dimension,ptr);
      for (i=0; i<dimension; i++)
        {
          *ptr=pnl_inv_cdfnor(*ptr);
          ptr++;
        }
      return;
    }
  for (i=0; i<dimension; i++)
    {
      *ptr=Gauss_BoxMuller(type_generator);
      ptr++;
    }
}

/**
 * return a matrix with its rows uniformly distributed on [a,b]^dimension.
 *
 * the samples have values in [a, b] (space of dimension dimension)
 *
 * @param M : the matrix of gaussian numbers, must already be allocated
 * @param samples : number of Monte Carlo samples (= number of rows of M)
 * @param dimension : dimension of the simulation (= number of columns of M)
 * @param a : lower bound vector of size dimension
 * @param b : upper bound vector of size dimension
 * @param type_generator : index of the generator 
 *
 * WARNING : The rows of M are indenpendent. This is very
 * important if QMC is used 
 */
void pnl_mat_rand_uni(PnlMat *M, int samples, int dimension,
                      const PnlVect *a, const PnlVect *b, int type_generator)
{
  int i, j;
  double *ptr, *aj, *bj;
  pnl_mat_resize(M,samples,dimension);
  ptr=M->array;

  if (pnl_rand_or_quasi (type_generator) == MC)
    {
      for(i=0;i<samples;i++)
        {
          aj=a->array; bj=b->array;
          for (j=0; j<dimension; j++)
            {
              pnl_Random[type_generator].Compute(1, ptr);
              *ptr=(*aj)+(*bj- (*aj))*(*ptr);
              ptr++; aj++; bj++;
            }
        }
      return;
    }
  CheckMaxQMCDim(type_generator, dimension);
  for(i=0;i<samples;i++)
    {
      pnl_Random[type_generator].Compute(dimension, ptr);
      aj=a->array; bj=b->array;
      for (j=0; j<dimension; j++)
        {
          *ptr=(*aj)+(*bj- (*aj))*(*ptr);
          ptr++; aj++; bj++;
        }
    }
}


/**
 * return a matrix with its rows uniformly distributed on [a,b]^dimension.
 *
 * the samples have values in [a, b] (space of dimension dimension)
 *
 * @param M : the matrix of gaussian numbers, must already be allocated
 * @param samples : number of Monte Carlo samples (= number of rows of M)
 * @param dimension : dimension of the simulation (= number of columns of M)
 * @param a : real lower bound 
 * @param b : real upper bound 
 * @param type_generator : index of the generator 
 *
 * WARNING : The rows of M are indenpendent. This is very
 * important if QMC is used 
 */
void pnl_mat_rand_uni2(PnlMat *M, int samples, int dimension,
                      double a, double b, int type_generator)
{
  int i, j;
  double *ptr;
  pnl_mat_resize(M,samples,dimension);
  ptr=M->array;

  if (pnl_rand_or_quasi (type_generator) == MC)
    {
      for(i=0;i<samples;i++)
        {
          for (j=0; j<dimension; j++)
            {
              pnl_Random[type_generator].Compute(1, ptr);
              *ptr = a + (b - a) * (*ptr);
              ptr++; 
            }
        }
      return;
    }
  CheckMaxQMCDim(type_generator, dimension);
  for(i=0;i<samples;i++)
    {
      pnl_Random[type_generator].Compute(dimension, ptr);
      for (j=0; j<dimension; j++)
        {
          *ptr = a + (b - a) * (*ptr);
          ptr++; 
        }
    }
}

/**
 * return a matrix with its rows normally distributed on R^dimension.
 * The samples have values in R^dimension
 *
 * @param M : the matrix of gaussian numbers, must already be allocated
 * @param samples : number of Monte Carlo samples (= number of rows of M)
 * @param dimension : dimension of the simulation (= number of columns of M)
 * @param type_generator : index of the generator 
 *
 * WARNING : The rows of M are indenpendent. This is very important if QMC is
 * used (independent dimensions). Each row represents a sample from the one
 * dimensionnal normal distribution
 */
void pnl_mat_rand_normal(PnlMat *M, int samples, int dimension,
                            int type_generator)
{
  int i, j;
  double *ptr;
  pnl_mat_resize(M,samples,dimension);
  ptr=M->array;
  if (pnl_rand_or_quasi(type_generator) == MC)
    {
      for (i=0; i<M->mn; i++)
        {
          *ptr = Gauss_BoxMuller(type_generator); ptr++;
        }
      return;
    }
  CheckMaxQMCDim(type_generator, dimension);
  for(i=0;i<samples;i++)
    {
      pnl_Random[type_generator].Compute(dimension, ptr);
      for (j=0; j<dimension; j++)
        {
          *ptr=pnl_inv_cdfnor(*ptr); ptr++;
        }
    }
}


/**
 * Simulates Gamma distribution
 *
 * @param a
 * @param b
 * @param gen the generator type
 *
 * New version based on Marsaglia and Tsang, "A Simple Method for
 * generating gamma variables", ACM Transactions on Mathematical
 * Software, Vol 26, No 3 (2000), p363-372.
 */

double pnl_rand_gamma (double a, double b, int gen)
{
  /* assume a > 0 */

  if (a < 1)
    {
      double u = pnl_rand_uni (gen);
      return pnl_rand_gamma ( 1.0 + a, b, gen) * pow (u, 1.0 / a);
    }

  {
    double x, v, u;
    double d = a - 1.0 / 3.0;
    double c = (1.0 / 3.0) / sqrt (d);

    while (1)
      {
        do
          {
            x = pnl_rand_normal (gen);
		     v = 1.0 + c * x;
          }
        while (v <= 0);

        v = v * v * v;
        u = pnl_rand_uni (gen);

        if (u < 1 - 0.0331 * x * x * x * x) 
          break;

        if (log (u) < 0.5 * x * x + d * (1 - v + log (v)))
          break;
      }
    
    return b * d * v;
  }
}

/**
 * Simulates a centered Chi square
 *
 * @param nu a real number, the number of degrees of freedom
 * @param gen the generator type
 *
 * The chisq distribution has the form
 *
 *  p(x) dx = (1/(2*Gamma(nu/2))) (x/2)^(nu/2 - 1) exp(-x/2) dx
 *
 * for x = 0 ... +infty
 */
double pnl_rand_chi2  (double nu, int gen)
{
  return 2. * pnl_rand_gamma ( nu / 2, 1.0, gen);
}



