
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
#include "pnl/pnl_specfun.h"
#include "pnl/pnl_cdf.h"

/* Maximum dimension for random sequences */
#define DIM_MAX 100000
static double ArrayOfRandomNumbers[DIM_MAX];

#define CheckQMCDim(rng, dim)                                       \
  if ( rng->dimension != dim )                                      \
    {                                                               \
      perror ("dimension cannot be changed for PNL_QMC\n"); abort ();   \
    }


/**
 * Simulate a standard random normal variable using Box Muller's algorithm
 * @param rng a PnlRng
 * @return a normal random variable
 */
static double Gauss_BoxMuller(PnlRng *rng)
{
  double xs,ys, g1, g2;

  if ( rng->has_gauss == 0 )
    {
      /* draw 2 new samples */
      do
        {
          rng->Compute(rng,&xs);
        }
      while (xs == 0.);
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
 * @param rng a generator
 */
static double GaussMC(int dimension, int create_or_retrieve, int index, PnlRng *rng)
{
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
 * @param rng a generator
 */
static double GaussQMC(int dimension, int create_or_retrieve, int index,PnlRng *rng)
{
  CheckQMCDim(rng, dimension);
  if (create_or_retrieve == CREATE)
    rng->Compute(rng,ArrayOfRandomNumbers);
  return pnl_inv_cdfnor(ArrayOfRandomNumbers[index]);
}


/*
 * Random number generation interface using the array PnlRngArray[]
 */
/**
 * Simulation of a Gaussian standard variable in dimension d
 * @param d size od the vector we are simulating
 * @param create_or_retrieve boolean can be CREATE or
 * RETRIEVE. if it is CREATE, draw all the dimensions and returns the fisrt
 * one. If it s RETRIEVE, returns the dimension corresponding to index
 * @param index index to be returned
 * @param rng a generator
 */
double pnl_rng_gauss(int d, int create_or_retrieve, int index, PnlRng *rng)
{
  if (rng->rand_or_quasi == PNL_QMC)
    return GaussQMC(d, create_or_retrieve, index, rng);
  else
    return GaussMC(d, create_or_retrieve, index, rng);

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
  if (x<p) return 1; else return 0;
}

/**
 * Simulation of a Poisson random variable using the PTRS algorithm
 * @param lambda parameter of the law
 * @param rng generator to use
 */
static long pnl_rng_poisson_ptrs(double lambda, PnlRng *rng)
{
  const double loglambda = log(lambda);
  const double b = 0.931 + 2.53 * sqrt(lambda);
  const double a = -0.059 + 0.02483 * b;
  const double invalpha = 1.1239 + 1.1328 / (b - 3.4);
  const double vr = 0.9277 - 3.6224 / (b - 2.);

  while (1)
    {
      double U = pnl_rng_uni(rng) - 0.5;
      double V = pnl_rng_uni(rng);
      double us = 0.5 - fabs(U);
      long k = (long) floor((2 * a / us + b) * U + lambda + 0.43);
      if ((us >= 0.07) && (V <= vr)) return k;
      if ((k < 0) || ((us < 0.013) && (V > us))) continue;
      if ((log(V) + log(invalpha) - log(a/(us*us)+b)) <= 
          (-lambda + k*loglambda - lgamma(k+1)))
        return k;
    }
}

/**
 * Simulation of a Poisson random variable
 * @param lambda parameter of the law
 * @param rng generator to use
 */
long pnl_rng_poisson(double lambda, PnlRng *rng)
{
  if (rng->rand_or_quasi == PNL_QMC)
    {
      int status, which;
      double x, bound, p, q;
      CheckQMCDim(rng, 1);
      rng->Compute(rng,&p);
      q = 1. - p;
      which = 2;
      pnl_cdf_poi(&which, &p, &q, &x, &lambda, &status, &bound);
      return (long) ceil(x);
    }
  else
    {
      if ( lambda < 10 )
        {
          double random_number;
          double u = 1;
          double elambda = exp(-lambda);
          long n = 0;

          while (1)
            {
              rng->Compute(rng, &random_number);
              u *= random_number;
              if ( u > elambda )
                n++;
              else
                return n;
            }
        }
      else 
        { 
          return pnl_rng_poisson_ptrs(lambda, rng);
        }
    }
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
    rng->Compute(rng, &x);
  } while (x == 0);
  return  (-log(x) / lambda);
}

/** 
 * Generate a non symmetric exponential distribution
 * 
 * @param lambda_p parameter of the positive exponential
 * @param lambda_m parameter of the negative exponential
 * @param p probability of being positive
 * @param rng generator to ise
 * 
 */
double pnl_rng_dblexp (double lambda_p, double lambda_m, double p, PnlRng *rng)
{
  if ( pnl_rng_uni (rng) < p )
    {
      return pnl_rng_exp (lambda_p, rng);
    }
  else
    {
      return - pnl_rng_exp (lambda_m, rng);
    }
}

/**
 * Simulation of a Poisson process
 * @param lambda parameter of the law
 * @param t time of the simulation
 * @param rng generator to use
 */
long pnl_rng_poisson1(double lambda, double t, PnlRng *rng)
{
  return pnl_rng_poisson(lambda * t, rng);
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
  if (rng->rand_or_quasi == PNL_QMC)
    {
      double u;
      CheckQMCDim(rng, 1);
      rng->Compute(rng,&u);
      return pnl_inv_cdfnor(u);
    }
  return Gauss_BoxMuller(rng);
}

/** 
 * Return of log normal random variable
 * 
 * @param m mean of the normal variable
 * @param sigma2 variance of the normal distribution
 * @param rng the generator to use
 */
double pnl_rng_lognormal (double m, double sigma2, PnlRng *rng)
{
  return exp (m + sqrt(sigma2) * pnl_rng_normal (rng));
}

/** 
 * Return an inverse gaussian random sample
 *
 * The algorithm used is based on 
 *
 *     John R. Michael, William R. Schucany and Roy W. Haas
 *     Generating Random Variates Using Transformations with Multiple Roots,
 *     The American Statistician, Vol. 30, No. 2 (May, 1976), pp. 88-90
 * 
 * @param mu mean of the distribution (must be positive)
 * @param lambda shape parameter (must be positive)
 * @param rng a random number generator
 * 
 */
double pnl_rng_invgauss (double mu, double lambda, PnlRng *rng)
{
  double g = pnl_rng_normal (rng);
  double v = g * g;
  double w = mu * v;
  double x1 = mu + mu / (2. * lambda) * (w - sqrt (w * (4 * lambda + w)));
  if ( pnl_rng_uni (rng) < mu / (x1 + mu) )
    {
      return x1;
    }
  else
    {
      return mu * mu / x1;
    }
}

/**
 * Simulate an iid sample of Bernoulli random variables
 * @param p parameter of the law
 * @param a, b the two values taken by the law
 * @param rng generator to use
 */
void pnl_vect_rng_bernoulli(PnlVect *V, int samples, double a, double b, double p, PnlRng *rng)
{
  int i;
  double x = 0.;
  pnl_vect_resize (V, samples);

  for ( i=0 ; i<samples ; i++ )
    {
      rng->Compute(rng,&x);
      PNL_LET (V, i) = (x<p) ? b : a;
    }
}

/**
 * Simulate a sample of a random vector with Bernoulli distribution
 * V(i) is a Bernoulli in {a(i), b(i)} with parameter p(i)
 *
 * @param p vector of parameters of the law
 * @param a, b the vectors of two values taken by the law
 * @param rng generator to use
 */
void pnl_vect_rng_bernoulli_d(PnlVect *V, int dimension, const PnlVect *a, const PnlVect *b, const PnlVect *p, PnlRng *rng)
{
  int i;
  pnl_vect_resize (V, dimension);
  if (rng->rand_or_quasi == PNL_QMC)
    {
      CheckQMCDim(rng, dimension);
      rng->Compute(rng, V->array);
    }
  else
    {
      for ( i=0 ; i<dimension ; i++ ) { rng->Compute(rng, &(PNL_LET(V, i))); }
    }
  for ( i=0 ; i<dimension ; i++ )
    {
      const double ai = PNL_GET(a, i);
      const double bi = PNL_GET(b, i);
      const double pi = PNL_GET(p, i);
      PNL_LET(V,i) = (PNL_GET(V,i) < pi) ? bi : ai; 
    }
}

/** 
 * Sample a vector of i.i.d rv following the Poisson distribution 
 * 
 * @param V contains the rv on output
 * @param samples number of samples
 * @param lambda Poisson parameter
 * @param rng generator to use
 */
void pnl_vect_rng_poisson(PnlVect *V, int samples, double lambda, PnlRng *rng)
{
  int i;
  pnl_vect_resize(V, samples);

  for ( i=0 ; i<samples ; i++ )
    {
      LET(V,i) = (double) pnl_rng_poisson(lambda, rng);
    }
}

/** 
 * Sample a vector from the multidimensional Poisson distribution 
 * 
 * @param V contains the rv on output
 * @param dimension size the vector lambda
 * @param lambda vector of Poisson parameters
 * @param rng generator to use
 */
void pnl_vect_rng_poisson_d(PnlVect *V, int dimension, const PnlVect *lambda, PnlRng *rng)
{
  int i;
  pnl_vect_resize(V, dimension);

  if (rng->rand_or_quasi == PNL_QMC)
    {
      int status, which;
      double bound, p, q, x, lam;
      which = 2;
      CheckQMCDim(rng, dimension);
      rng->Compute(rng, V->array);
      for ( i=0 ; i<dimension ; i++ )
        {
          p = PNL_GET(V, i);
          lam = GET(lambda, i);
          q = 1. - p;
          pnl_cdf_poi(&which, &p, &q, &x, &lam, &status, &bound);
          PNL_LET(V, i) = ceil(x);
        }
    }
  else
    {
      for ( i=0 ; i<dimension ; i++ )
        {
          PNL_LET(V,i) = (double) pnl_rng_poisson(PNL_GET(lambda, i), rng);
        }
    }
}

/**
 * Compute a vector of independent and uniformly distributed r.v. on [a,b]
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
 * Compute a random vector uniformly distributed on [a,b]^dimension
 *
 * if the generator is a true PNL_MC generator, no difference between this
 * function and pnl_vect_rng_uni. In case of a PNL_QMC generator, this
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
  if (rng->rand_or_quasi == PNL_QMC)
    {
      CheckQMCDim(rng, dimension);
      rng->Compute(rng, G->array);
      for(i=0;i<dimension;i++)
        {
          PNL_LET(G,i) = a+(b-a)*PNL_GET(G,i);
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
 * Compute a vector of independent and normaly distributed r.v. on R
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
  if (rng->rand_or_quasi == PNL_QMC)
    {
      CheckQMCDim(rng, 1);
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
 * Compute a random vector normally distributed on R^dimension.
 *
 * if the generator is a true PNL_MC generator, no difference between this
 * function and pnl_vect_rng_uni. In case of a PNL_QMC generator, this
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
  if (rng->rand_or_quasi == PNL_QMC)
    {
      CheckQMCDim(rng, dimension);
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
 * Compute a matrix with independent rows following a Bernoulli distribution
 * onver [a,b] with parameter p (a, b and p are vectors)
 *
 * the samples have values in [a, b] (space of dimension dimension)
 *
 * @param M  the matrix of random numbers, must already be allocated
 * @param samples  number of Monte Carlo samples (= number of rows of M)
 * @param dimension  dimension of the simulation (= number of columns of M)
 * @param a, b values of the Bernoulli distribution
 * @param p parameter of the distribution
 * @param rng generator to use
 *
 * WARNING : The rows of M are indenpendent. This is very
 * important if PNL_QMC is used
 */
void pnl_mat_bernoulli(PnlMat *M, int samples, int dimension, const
                           PnlVect *a, const PnlVect *b, const PnlVect *p,
                           PnlRng *rng)
{
  int i, j;
  pnl_mat_resize(M, samples, dimension);

  if (rng->rand_or_quasi == PNL_QMC)
    {
      CheckQMCDim(rng, dimension);
      for ( i=0 ; i<samples ; i++ ) { rng->Compute(rng, &(PNL_MGET(M, i, 0))); }
    }
  else
    {
      for ( i=0 ; i<M->mn ; i++ ) { rng->Compute(rng, M->array + i); }
    }


  for ( j=0 ; j<dimension ; j++ )
    {
      const double aj = PNL_GET(a, j);
      const double bj = PNL_GET(b, j);
      const double pj = PNL_GET(p, j);
      for ( i=0 ; i<samples ; i++ )
        {
          PNL_MLET(M, i, j) = (PNL_MGET(M, i, j) < pj) ? bj : aj; 
        }
    }
}

/** 
 * Sample a iid vectors from the multidimensional Poisson distribution. The
 * rows are independent.
 * 
 * @param M contains the random matrix on output (samples x dimension)
 * @param dimension size the vector lambda
 * @param lambda vector of Poisson parameters
 * @param rng generator to use
 */
void pnl_mat_rng_poisson(PnlMat *M, int samples, int dimension, const PnlVect *lambda, PnlRng *rng)
{
  int i, j;
  pnl_mat_resize(M, samples, dimension);

  if (rng->rand_or_quasi == PNL_QMC)
    {
      int status, which;
      double bound, p, q, x, lam;
      which = 2;
      CheckQMCDim(rng, dimension);
      for ( i=0 ; i<samples ; i++ )
        {
          rng->Compute(rng, &(PNL_MLET(M, i, 0)));
          for ( j=0 ; j<dimension ; j++ )
            {
              p = PNL_MGET(M, i, j);
              lam = GET(lambda, j);
              q = 1. - p;
              pnl_cdf_poi(&which, &p, &q, &x, &lam, &status, &bound);
              PNL_MLET(M, i, j) = ceil(x);
            }
        }
    }
  else
    {

      for ( i=0 ; i<samples ; i++ )
        {
          for ( j=0 ; j<dimension ; j++ )
            {
              PNL_MLET(M, i, j) = (double) pnl_rng_poisson(PNL_GET(lambda, j), rng);
            }
        }
    }
}

/**
 * Compute a matrix with independent and uniformly distributed rows on [a,b]
 * ( a and b are vectors )
 *
 * the samples have values in [a, b] (space of dimension dimension)
 *
 * @param M : the matrix of random numbers, must already be allocated
 * @param samples : number of Monte Carlo samples (= number of rows of M)
 * @param dimension : dimension of the simulation (= number of columns of M)
 * @param a : lower bound vector of size dimension
 * @param b : upper bound vector of size dimension
 * @param rng generator to use
 *
 * WARNING : The rows of M are indenpendent. This is very
 * important if PNL_QMC is used
 */
void pnl_mat_rng_uni(PnlMat *M, int samples, int dimension,
                     const PnlVect *a, const PnlVect *b, PnlRng *rng)
{
  int i, j;
  double u;
  pnl_mat_resize(M,samples,dimension);

  if (rng->rand_or_quasi == PNL_MC)
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
  CheckQMCDim(rng, dimension);
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
 * Compute a  matrix with independent and uniformly distributed rows on [a,b]^dimension.
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
 * important if PNL_QMC is used
 */
void pnl_mat_rng_uni2(PnlMat *M, int samples, int dimension,
                      double a, double b, PnlRng *rng)
{
  int i, j;
  double u;
  pnl_mat_resize(M,samples,dimension);

  if (rng->rand_or_quasi == PNL_MC)
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
  for(i=0;i<samples;i++)
    {
      CheckQMCDim(rng, dimension);
      rng->Compute(rng, &(PNL_MGET(M, i, 0)));
      for (j=0; j<dimension; j++)
        {
          PNL_MLET(M,i,j)=a+(b -a)*PNL_MGET(M,i,j);
        }
    }
}

/**
 * Compute a matrix with independent and normally distributed rows on R^dimension.
 * The samples have values in R^dimension
 *
 * @param M : the matrix of gaussian numbers, must already be allocated
 * @param samples : number of Monte Carlo samples (= number of rows of M)
 * @param dimension : dimension of the simulation (= number of columns of M)
 * @param rng generator to use
 *
 * WARNING : The rows of M are indenpendent. This is very important if PNL_QMC is
 * used (independent dimensions). Each row represents a sample from the one
 * dimensionnal normal distribution
 */
void pnl_mat_rng_normal(PnlMat *M, int samples, int dimension, PnlRng *rng)
{
  int i, j;
  pnl_mat_resize(M,samples,dimension);
  if (rng->rand_or_quasi == PNL_MC)
    {
      for (i=0; i<M->mn; i++)
        {
          M->array[i] = Gauss_BoxMuller(rng);
        }
      return;
    }
  for(i=0;i<samples;i++)
    {
      CheckQMCDim(rng, dimension);
      rng->Compute(rng, &(PNL_MGET(M, i, 0)));
      for (j=0; j<dimension; j++)
        {
          PNL_MLET(M,i,j)=pnl_inv_cdfnor(PNL_MGET(M,i,j));
        }
    }
}

/**
 * Simulate Gamma distribution
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

  if ( rng->rand_or_quasi == PNL_QMC )
    {
      int status, which;
      double x, bound, p, q;
      CheckQMCDim(rng, 1);
      rng->Compute(rng,&p);
      q = 1. - p;
      which = 2;
      /* The function pnl_cdf_gam uses the rate instead of the scale */
      b = 1. / b; 
      pnl_cdf_gam(&which, &p, &q, &x, &a, &b, &status, &bound);
      return x;
    }
  /* assume a > 0 */
  if (a < 1)
    {
      double u = pnl_rng_uni (rng);
      return pnl_rng_gamma ( 1.0 + a, b, rng) * pow (u, 1.0 / a);
    }

  {
    double x, v, u, d, c;
    d = a - 1.0 / 3.0;
    c = (1.0 / 3.0) / sqrt (d);

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
 * Simulate a centered Chi square
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

/** 
 * Generate a random variable with Bessel distribution with parameters nu
 * and a. We refer to 
 * 
 * article{Devroye2002249,
 * title = "Simulating Bessel random variables",
 * journal = "Statistics & Probability Letters",
 * volume = "57",
 * number = "3",
 * year = "2002",
 * doi = "DOI: 10.1016/S0167-7152(02)00055-X",
 * author = "Luc Devroye"}
 *
 * for a detailed review of the different algorithms to simulate Bessel
 * random variables. We are using the most naive which also turns to be the
 * fastest (probably because we are relying on libamis for evaluating
 * Bessel functions).
 * 
 * @param nu a real number > -1
 * @param a a real number > 0
 * @param rng a PnlRng
 * 
 * @return an integer number
 */
int pnl_rng_bessel (double nu, double a, PnlRng *rng)
{
  double p0;
  double tmp,u;
  int n;

  p0 = pow (a * 0.5, nu) / ( pnl_bessel_i (nu, a) * pnl_sf_gamma (nu + 1.) );
  u = pnl_rng_uni (rng);
  tmp = 0;
  n = 0;
  if ( u <= p0 ) return 0;
  do
    {                                                                        
      tmp += p0;
      p0 = p0 * a * a / (4. * (n + 1) * (n + 1 + nu));
      n++;
    } while ( u> tmp + p0 );
  return n;
}

/*
 * Deprecated funtions using the old rand interface
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
  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);
  return pnl_rng_gauss (d, create_or_retrieve, index, rng);
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
 * Compute a vector of independent and uniformly distributed r.v. on [a,b]
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
 * Compute a random vector uniformly distributed on [a,b]^dimension
 *
 * if the generator is a true PNL_MC generator, no difference between this
 * function and pnl_vect_rand_uni. In case of a PNL_QMC generator, this
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
 * Compute a vector of independent and normaly distributed r.v. on R
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
 * Compute a random vector normally distributed on R^dimension.
 *
 * if the generator is a true PNL_MC generator, no difference between this
 * function and pnl_vect_rand_uni. In case of a PNL_QMC generator, this
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
 * Compute a matrix with independent and uniformly distributed rows on [a,b]
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
 * important if PNL_QMC is used
 */
void pnl_mat_rand_uni(PnlMat *M, int samples, int dimension,
                      const PnlVect *a, const PnlVect *b, int type_generator)
{
  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);
  pnl_mat_rng_uni (M, samples, dimension, a, b, rng);
}

/**
 * Compute a matrix with independent and uniformly distributed rows on [a,b]^dimension.
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
 * important if PNL_QMC is used
 */
void pnl_mat_rand_uni2(PnlMat *M, int samples, int dimension,
                       double a, double b, int type_generator)
{

  PnlRng *rng;
  rng = pnl_rng_get_from_id(type_generator);
  pnl_mat_rng_uni2 (M, samples, dimension, a, b, rng);
}

/**
 * Compute a matrix with independent and normally distributed rows on R^dimension.
 * The samples have values in R^dimension
 *
 * @param M : the matrix of gaussian numbers, must already be allocated
 * @param samples : number of Monte Carlo samples (= number of rows of M)
 * @param dimension : dimension of the simulation (= number of columns of M)
 * @param type_generator : index of the generator
 *
 * WARNING : The rows of M are indenpendent. This is very important if PNL_QMC is
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
 * Simulate Gamma distribution
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
 * Simulate a centered Chi square
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

/** 
 * Generate a random variable with Bessel distribution with parameters nu
 * and a
 * 
 * @param nu a real number > -1
 * @param a a real number > 0
 * @param gen the index of the generator to be used
 * 
 * @return an integer number
 */
int pnl_rand_bessel (double nu, double a, int gen)
{
  PnlRng *rng;
  rng = pnl_rng_get_from_id(gen);
  return pnl_rng_bessel (nu, a, rng);
}

