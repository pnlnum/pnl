
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

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_random.h"
#include "pnl/pnl_cdf.h"
#include "randomkit.h"

/* Maximum dimension for random sequences */
#define DIM_MAX 100000

#define CheckMaxQMCDim(rng, dimension)                        \
{   if (dimension >= rng->max_dim)                            \
    {                                                         \
      perror("maximum dimension of QMC exceeded\n"); abort(); \
    }                                                         \
}

static double ArrayOfRandomNumbers[DIM_MAX];


/**
 * Simulates a standard random normal variable using Box Muller's algorithm
 * @param type_generator index of the generator to be used
 * @return a normal random variable
 */    
static double Gauss_BoxMuller(PnlRng *rng)
{
  double xs,ys, g1, g2;

  /* do not wast any samples. But be sure to throw away any remaining
     samples when pnl_rand_init is called */
  if (rng->counter == 1 || rng->has_gauss == 0)
    {
      /* draw 2 new samples */
      rng->Compute(rng,&xs);
      rng->Compute(rng,&ys);
      xs = sqrt( -2. * log(xs) );
      ys = M_2PI * ys;
      g1 = xs * cos(ys);
      g2 = xs * sin(ys);
      rng->gauss = g2;
      rng->has_gauss = 1;
      return g1;
    }
  else
    {
      /* use the remaining sample from the last call */
      rng->has_gauss = 0;
      return rng->gauss;
    }
}


/**
 * Simulation of a Gaussian standard variable.
 *
 * @param dimension size of the vector to simulate
 * @param create_or_retrieve boolean can be CREATE or RETRIEVE. 
 * @param index (unused when calling with CREATE)
 * @param type_generator index of the generator
 */
static double GaussMC(int dimension, int create_or_retrieve, int index, int type_generator)
{
  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);

  if (create_or_retrieve == CREATE)
    {
      int i;
      for (i=0; i<dimension; i++)
        ArrayOfRandomNumbers[i] = Gauss_BoxMuller(rng);
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
  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);
  if (dimension > rng->dimension)
    {
      perror("maximum dimension of Monte Carlo exceeded\n"); abort();
    }

  if (create_or_retrieve == CREATE)
    rng->Compute(rng,ArrayOfRandomNumbers);
  return pnl_inv_cdfnor(ArrayOfRandomNumbers[index]);
}



/*
 * Random number generation interface using the array PnlRngArray[] and the
 * dynamic list PnlRngList
 */

/**
 * Simulation of a Gaussian standard variable in dimension d
 * @param d size od the vector we are simulating
 * @param create_or_retrieve boolean can be CREATE or
 * RETRIEVE. if it is CREATE, draw all the dimensions and returns the fisrt
 * one. If it s RETRIEVE, returns the dimension corresponding to index
 * @param index index to be returned
 * @param type_generator index of the generator
 */
double pnl_rand_gauss(int d, int create_or_retrieve, int index, int type_generator)
{
  if (pnl_rand_or_quasi (type_generator) == QMC)
    return GaussQMC(d, create_or_retrieve, index, type_generator);
  else
    return GaussMC(d, create_or_retrieve, index, type_generator);

}

/**
 * Simulation of a Bernoulli random variable
 * @param p parameter of the law
 * @param type_generator index of the generator ot be used
 */
int pnl_rand_bernoulli(double p, int type_generator)
{
  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);
  return pnl_rng_bernoulli (p, rng);
}

/**
 * Simulation of a Poisson random variable 
 * @param lambda parameter of the law
 * @param type_generator index of the generator ot be used
 */
long pnl_rand_poisson(double lambda, int type_generator)
{
  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);
  return pnl_rng_poisson(lambda, rng);
}

/**
 * Simulation of an exponential random variable
 * @param lambda parameter of the law
 * @param type_generator index of the generator ot be used
 */
double pnl_rand_exp(double lambda,int type_generator)
{
  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);
  return pnl_rng_exp (lambda, rng);

}


/**
 * Simulation of a Poisson process 
 * @param lambda parameter of the law
 * @param t time of the simulation
 * @param type_generator index of the generator ot be used
 */
long pnl_rand_poisson1(double lambda, double t, int type_generator)
{
  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);
  return pnl_rng_poisson1(lambda, t, rng);
}


/**
 * Generate a uniformly distributed number on ]0,1).
 * @param type_generator index ot the generator to be used
 *
 * @see pnl_rand_uni_ab
 */
double pnl_rand_uni (int type_generator)
{
  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);
  return pnl_rng_uni(rng);
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

  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);
  return pnl_rng_uni_ab (a, b, rng);
}


/**
 * Generate a normally distributed number.
 * @param type_generator index ot the generator to be used
 */
double pnl_rand_normal (int type_generator)
{
  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);
  return pnl_rng_normal (rng);
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
  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);
  pnl_vect_rng_uni (G, samples, a, b, rng);
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
  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);
  pnl_vect_rng_uni_d(G, dimension, a, b, rng);
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
  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);
  pnl_vect_rng_normal (G, samples, rng);
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
  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);
  pnl_vect_rng_normal_d (G, dimension, rng);
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
  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);
  pnl_mat_rng_uni (M, samples, dimension, a, b, rng);
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

  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);
  pnl_mat_rng_uni2 (M, samples, dimension, a, b, rng);
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
  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);
  pnl_mat_rng_normal (M, samples, dimension, rng);
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
  PnlRng *rng;
  rng = pnl_rng_get_from_id(gen);
  return pnl_rng_gamma (a, b, rng);
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
  PnlRng *rng;
  rng = pnl_rng_get_from_id(gen);
  return pnl_rng_chi2 (nu, rng);
}

/*
 * Interface using PnlRng
 */

/**
 * Simulation of a Bernoulli random variable
 * @param p parameter of the law
 * @param rng generator to use
 */
int pnl_rng_bernoulli(double p, PnlRng *rng)
{
  double x=0.0;
  rng->Compute(rng,&x);
  if (x<p) return 1; else return  0;
}

/**
 * Simulation of a Poisson random variable 
 * @param lambda parameter of the law
 * @param rng generator to use
 */
long pnl_rng_poisson(double lambda, PnlRng *rng)
{
  double u;
  double a = exp(-lambda);
  long n = 0;
  double random_number;

  rng->Compute(rng,&random_number);
  u = random_number;
  while (u>a)
    {
      rng->Compute(rng,&random_number); 
      u *= random_number;
      n++;
    }
  return n;
}

/**
 * Simulation of an exponential random variable
 * @param lambda parameter of the law
 * @param rng generator to use
 */
double pnl_rng_exp(double lambda,PnlRng *rng)
{
  double x;

  do{
    rng->Compute(rng,&x); 
  } while(x==0);
  return (double) (-log(x)/lambda);
}


/**
 * Simulation of a Poisson process 
 * @param lambda parameter of the law
 * @param t time of the simulation
 * @param rng generator to use
 */
long pnl_rng_poisson1(double lambda, double t, PnlRng *rng)
{
  double S;
  long Nt;
  Nt=0;
  S=0;
  do {
    S=S+pnl_rng_exp (lambda,rng);
    Nt=Nt+1;
  } while (S<=t);
  return Nt-1;
}


/**
 * Generate a uniformly distributed number on ]0,1).
 * @param rng generator to use
 *
 * @see pnl_rng_uni_ab
 */
double pnl_rng_uni (PnlRng *rng)
{
  double u;
  do { rng->Compute(rng,&u); }
  while (u == 0); 
  return u;
}

/**
 * Generate a uniformly distributed number on [a,b].
 * @param a lower bound
 * @param b upper bound
 * @param rng generator to use
 *
 * @see pnl_rng_uni
 */
double pnl_rng_uni_ab (double a, double b, PnlRng *rng)
{
  double u;
  rng->Compute(rng,&u);
  return a+(b-a)*u;
}


/**
 * Generate a normally distributed number.
 * @param rng generator to use
 */
double pnl_rng_normal (PnlRng *rng)
{
  if (rng->rand_or_quasi == QMC)
    {
      double u;
      rng->Compute(rng,&u);
      return pnl_inv_cdfnor(u);
    }
  return Gauss_BoxMuller(rng);
}

/**
 * return a vector of uniformly distributed components on [a,b]
 * @param G existing PnlVect containing the random numbers on exit
 * @param samples size of G (number of independent samples requested)
 * @param a lower bound 
 * @param b upper bound
 * @param rng generator to use 
 *
 * @see pnl_vect_rng_uni_d
 */
void pnl_vect_rng_uni(PnlVect *G, int samples, double a, double b, PnlRng *rng)
{
  int i;
  double u;
  pnl_vect_resize(G,samples);
  for(i=0;i<samples;i++)
    {
      rng->Compute(rng,&u);
      PNL_LET(G, i) = a+(b-a)*u;
    }
}


/**
 * return a vector uniformly distributed on [a,b]^dimension
 *
 * if the generator is a true MC generator, no difference between this
 * function and pnl_vect_rng_uni. In case of a QMC generator, this
 * function generator one sample with values in [a,b]^dimension and NOT
 * dimension samples with values in [a, b].
 *
 * @param G existing PnlVect containing the random numbers on exit
 * @param dimension dimension of the state space
 * @param a lower bound 
 * @param b upper bound 
 * @param rng generator to use
 *
 * @see pnl_vect_rng_uni
 */
void pnl_vect_rng_uni_d (PnlVect *G, int dimension, double a, double b, PnlRng *rng)
{
  int i;
  double u;
  pnl_vect_resize(G,dimension);
  if (rng->rand_or_quasi == QMC)
    {
      CheckMaxQMCDim(rng, dimension);
      for(i=0;i<dimension;i++)
        {
          rng->Compute(rng, &u);
          PNL_LET(G,i) = a+(b-a)*u;
        }
      return;
    }
  for(i=0;i<dimension;i++)
    {
      rng->Compute(rng, &u);
      PNL_LET(G,i) = a+(b-a)*u;
    }
}


/**
 * return a vector of normaly distributed components on R
 *
 * @param samples number of samples
 * @param G : the vector of gaussian numbers, must already be allocated.
 * @param rng generator to use
 *
 * @see pnl_vect_rng_normal_d
 */
void pnl_vect_rng_normal (PnlVect *G, int samples, PnlRng *rng)
{
  int i;
  double u;
  pnl_vect_resize(G,samples);
  if (rng->rand_or_quasi == QMC)
    {
      for (i=0; i<samples; i++)
        {
          rng->Compute(rng,&u);
          PNL_LET(G,i) = pnl_inv_cdfnor(u);
        }
      return;
    }
  for (i=0; i<samples; i++)
    {
      PNL_LET(G,i) = Gauss_BoxMuller(rng);
    }
}


/**
 * return a vector normally distributed on R^dimension.
 *
 * if the generator is a true MC generator, no difference between this
 * function and pnl_vect_rng_uni. In case of a QMC generator, this
 * function generator one sample with values in R^dimension and NOT
 * dimension samples with values in R.
 *
 * @param dimension : size of the vector. one sample of a Gaussian vector.
 * @param G : the vector of gaussian numbers, must already be allocated.
 * @param rng generator to use
 *
 * @see pnl_vect_rng_normal
 */
void pnl_vect_rng_normal_d (PnlVect *G, int dimension, PnlRng *rng)
{
  int i;
  pnl_vect_resize(G,dimension);
  if (rng->rand_or_quasi == QMC)
    {
      CheckMaxQMCDim(rng, dimension);
      rng->Compute(rng,G->array);
      for (i=0; i<dimension; i++)
        {
          PNL_LET(G,i) = pnl_inv_cdfnor(PNL_GET(G,i));
        }
      return;
    }
  for (i=0; i<dimension; i++)
    {
      PNL_LET(G,i) = Gauss_BoxMuller(rng);
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
 * @param rng generator to use
 *
 * WARNING : The rows of M are indenpendent. This is very
 * important if QMC is used 
 */
void pnl_mat_rng_uni(PnlMat *M, int samples, int dimension,
                      const PnlVect *a, const PnlVect *b, PnlRng *rng)
{
  int i, j;
  double u;
  pnl_mat_resize(M,samples,dimension);

  if (rng->rand_or_quasi == MC)
    {
      for(i=0;i<samples;i++)
        {
          for (j=0; j<dimension; j++)
            {
              rng->Compute(rng, &u);
              PNL_MLET(M,i,j)=PNL_GET(a,j)+(PNL_GET(b,j)- PNL_GET(a,j))*u;
            }
        }
      return;
    }
  CheckMaxQMCDim(rng, dimension);
  for(i=0;i<samples;i++)
    {
      rng->Compute(rng, &(PNL_MGET(M, i, 0)));
      for (j=0; j<dimension; j++)
        {
          PNL_MLET(M,i,j)=PNL_GET(a,j)+(PNL_GET(b,j)- PNL_GET(a,j))*PNL_MGET(M,i,j);
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
 * @param rng generator to use
 *
 * WARNING : The rows of M are indenpendent. This is very
 * important if QMC is used 
 */
void pnl_mat_rng_uni2(PnlMat *M, int samples, int dimension,
                       double a, double b, PnlRng *rng)
{
  int i, j;
  double u;
  pnl_mat_resize(M,samples,dimension);

  if (rng->rand_or_quasi == MC)
    {
      for(i=0;i<samples;i++)
        {
          for (j=0; j<dimension; j++)
            {
              rng->Compute(rng, &u);
              PNL_MLET(M,i,j)=a+(b- a)*u;
            }
        }
      return;
    }
  CheckMaxQMCDim(rng, dimension);
  for(i=0;i<samples;i++)
    {
      rng->Compute(rng, &(PNL_MGET(M, i, 0)));
      for (j=0; j<dimension; j++)
        {
          PNL_MLET(M,i,j)=a+(b -a)*PNL_MGET(M,i,j);
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
 * @param rng generator to use
 *
 * WARNING : The rows of M are indenpendent. This is very important if QMC is
 * used (independent dimensions). Each row represents a sample from the one
 * dimensionnal normal distribution
 */
void pnl_mat_rng_normal(PnlMat *M, int samples, int dimension,
                         PnlRng *rng)
{
  int i, j;
  pnl_mat_resize(M,samples,dimension);
  if (rng->rand_or_quasi == MC)
    {
      for (i=0; i<M->mn; i++)
        {
          M->array[i] = Gauss_BoxMuller(rng);
        }
      return;
    }
  CheckMaxQMCDim(rng, dimension);
  for(i=0;i<samples;i++)
    {
      rng->Compute(rng, &(PNL_MGET(M, i, 0)));
      for (j=0; j<dimension; j++)
        {
          PNL_MLET(M,i,j)=pnl_inv_cdfnor(PNL_MGET(M,i,j));
        }
    }
}


/**
 * Simulates Gamma distribution
 *
 * @param a
 * @param b
 * @param rng generator to use
 *
 * New version based on Marsaglia and Tsang, "A Simple Method for
 * generating gamma variables", ACM Transactions on Mathematical
 * Software, Vol 26, No 3 (2000), p363-372.
 */
double pnl_rng_gamma (double a, double b, PnlRng *rng)
{
  /* assume a > 0 */

  if (a < 1)
    {
      double u = pnl_rng_uni (rng);
      return pnl_rng_gamma ( 1.0 + a, b, rng) * pow (u, 1.0 / a);
    }

    {
      double x, v, u;
      double d = a - 1.0 / 3.0;
      double c = (1.0 / 3.0) / sqrt (d);

      while (1)
        {
          do
            {
              x = pnl_rng_normal (rng);
              v = 1.0 + c * x;
            }
          while (v <= 0);

          v = v * v * v;
          u = pnl_rng_uni (rng);

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
 * @param rng generator to use
 *
 * The chisq distribution has the form
 *
 *  p(x) dx = (1/(2*Gamma(nu/2))) (x/2)^(nu/2 - 1) exp(-x/2) dx
 *
 * for x = 0 ... +infty
 */
double pnl_rng_chi2  (double nu, PnlRng *rng)
{
  return 2. * pnl_rng_gamma ( nu / 2, 1.0, rng);
}


