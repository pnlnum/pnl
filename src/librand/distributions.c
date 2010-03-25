
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

#include "pnl_mathtools.h"
#include "pnl_random.h"
#include "pnl_cdf.h"
#include "randomkit.h"

/* Maximum dimension for random sequences */
#define DIM_MAX 100000
#define CheckMaxQMCDim(type_generator, dimension)                       \
  {   if (dimension >= pnl_Random[type_generator].Dimension)         \
      {                                                                 \
        perror("maximum dimension of QMC exceeded\n"); abort(); \
      }                                                                 \
  }
static double ArrayOfRandomNumbers[DIM_MAX];
static int draw_new_sample = 0;

extern void pnl_rand_reset_all ()
{
  draw_new_sample = 1;
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


