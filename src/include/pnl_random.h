#ifndef _PNL_RANDOM_H 
#define _PNL_RANDOM_H 


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl_mathtools.h"
#include "pnl_types.h"
#include "pnl_vector.h"
#include "pnl_matrix.h"

/**
 * \defgroup PnlRandom Random generators 
 */
/*@{*/

#define MC 0	
#define QMC 1	
#define CREATE 0 
#define RETRIEVE 1 
#define GEN_NUMBER 20

/* indices of random generators 
 * The first generators must be true MC generators and the last ones are QMC
 * Do NOT mix them otherwise the arrays  pnl_random and pnl_random_MC will be
 * wrong */
#define PNL_RNG_KNUTH 0
#define PNL_RNG_MRGK3 1
#define PNL_RNG_MRGK5 2
#define PNL_RNG_SHUFL 3
#define PNL_RNG_L_ECUYER 4
#define PNL_RNG_TAUSWORTHE 5
#define PNL_RNG_MERSENNE 6
#define PNL_RNG_MERSENNE_RANDOM_SEED 7
/* beginning of QMC sequences */
#define PNL_RNG_SQRT 8
#define PNL_RNG_HALTON 9
#define PNL_RNG_FAURE 10
#define PNL_RNG_SOBOL 11
#define PNL_RNG_SOBOL2 12
#define PNL_RNG_NIEDERREITER 13



/*RandomGenerators*/
typedef struct 
{
  const char * Name;
  void (*Compute)(int, double *);
  int RandOrQuasi;
  int Dimension;
} PnlGenerator;

typedef struct 
{
  enum_member base;
  void (*Compute)(int, double *);
  int RandOrQuasi;
  int Dimension;
} random_generator;

extern random_generator pnl_Random_MC[];
extern random_generator pnl_Random[];
extern enum_members RNGs;
extern enum_members MC_RNGs;

extern int pnl_rand_init(int type_generator, int simulation_dim,long samples);
extern int pnl_rand_or_quasi(int type_generator);
extern const char * pnl_rand_name (int type_generator);
extern double pnl_rand_gauss(int, int, int, int);
extern int pnl_rand_bernoulli(double p, int generator);
extern long pnl_rand_poisson(double lambda, int type_generator);
extern long pnl_rand_poisson1(double lambda,double t, int type_generator);
extern double pnl_rand_exp(double lambda,int type_generator);

extern double pnl_rand_uni (int type_generator);
extern double pnl_rand_uni_ab (double a, double b, int type_generator);
extern double pnl_rand_normal (int type_generator);

extern void
pnl_vect_rand_uni(PnlVect *G, int samples, double a, double b, int type_generator);
extern void
pnl_vect_rand_uni_d(PnlVect *G, int dimension, double a, double b, int type_generator);
extern void pnl_vect_rand_normal(PnlVect *G, int samples, int generator);
extern void pnl_vect_rand_normal_d(PnlVect *G, int dimension, int generator);
extern void
pnl_mat_rand_uni(PnlMat *M, int samples, int dimension, const PnlVect *a,
                 const PnlVect *b, int type_generator);
extern void pnl_mat_rand_uni2(PnlMat *M, int samples, int dimension,
                              double a, double b, int type_generator);
extern void
pnl_mat_rand_normal(PnlMat *M, int samples, int dimension, int type_generator);
extern double pnl_rand_gamma (double a, double b, int gen);
extern double pnl_rand_chi2  (double nu, int gen);

/*@}*/


#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_RANDOM_H */ 
