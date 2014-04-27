#ifndef _PNL_RANDOM_H
#define _PNL_RANDOM_H

#include "pnl/pnl_config.h"
#include "pnl/pnl_object.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**
 * \defgroup PnlRandom Random generators
 */
/*@{*/

#define PNL_MC 0
#define PNL_QMC 1
#define CREATE 0
#define RETRIEVE 1

/* indices of random generators
 * The first generators must be true PNL_MC generators and the last ones are PNL_QMC
 * Do NOT mix them */
typedef enum 
{
  PNL_RNG_NULL                 = -1,
  PNL_RNG_KNUTH                = 0,
  PNL_RNG_MRGK3                = 1,
  PNL_RNG_MRGK5                = 2,
  PNL_RNG_SHUFL                = 3,
  PNL_RNG_LECUYER              = 4,
  PNL_RNG_TAUSWORTHE           = 5,
  PNL_RNG_MERSENNE             = 6,
  PNL_RNG_MERSENNE_RANDOM_SEED = 7,
  /* beginning of PNL_QMC sequences */
  PNL_RNG_SQRT                 = 8,
  PNL_RNG_HALTON               = 9,
  PNL_RNG_FAURE                = 10,
  PNL_RNG_SOBOL_I4             = 11,
  PNL_RNG_SOBOL_I8             = 12,
  PNL_RNG_NIEDERREITER         = 13,
  /* end of static rng */
  PNL_RNG_DCMT                 = 14
} PnlRngType;

/**
 * RandomGenerators
 */
typedef unsigned long ulong;

typedef struct _PnlRng PnlRng;
struct _PnlRng
{
  PnlObject object;
  PnlRngType type; /*!< generator type */
  void (*Compute)(PnlRng *g, double *sample); /*!< the function to compute the
                                                next number in the sequence */
  int rand_or_quasi; /*!< can be PNL_MC or PNL_QMC */
  int dimension; /*!< dimension of the space in which we draw the samples */
  int counter; /*!< counter = number of samples already drawn */
  int has_gauss; /*!< Is a gaussian deviate available? */
  double gauss; /*!< If has_gauss==1, gauss a gaussian sample */
  int size_state; /*!< size in bytes of the state variable */
  void *state; /*!< state of the random generator */
};


extern PnlRng PnlRngKnuth;
extern PnlRng PnlRngMrgk3;
extern PnlRng PnlRngMrgk5;
extern PnlRng PnlRngShufl;
extern PnlRng PnlRngLecuyer;
extern PnlRng PnlRngTausworthe;
extern PnlRng PnlRngMersenne;
extern PnlRng PnlRngMersenneRandomSeed;
extern PnlRng PnlRngSqrt;
extern PnlRng PnlRngHalton;
extern PnlRng PnlRngFaure;
extern PnlRng PnlRngSobolI4;
extern PnlRng PnlRngSobolI8;
extern PnlRng PnlRngNiederreiter;
extern PnlRng *PnlRngArray[];

#ifdef PNL_HAVE_INLINE
PNL_INLINE_FUNC PnlRng* pnl_rng_get_from_id (PnlRngType t)
{
  return PnlRngArray[t];
}
#endif
PNL_INLINE_DECL PnlRng* pnl_rng_get_from_id (PnlRngType t);

extern char* pnl_rng_get_name (PnlRngType t);


/*
 * Rand functions
 */
extern int pnl_rand_init(int type_generator, int simulation_dim,long samples);
extern void pnl_rand_sseed (int type_generator, ulong seed);
extern int pnl_rand_or_quasi(int type_generator);
extern double pnl_rng_gauss(int, int, int, PnlRng*);
extern double pnl_rand_gauss(int, int, int, int);
extern int pnl_rand_bernoulli(double p, int generator);
extern long pnl_rand_poisson(double lambda, int type_generator);
extern long pnl_rand_poisson1(double lambda,double t, int type_generator);
extern double pnl_rand_exp(double lambda,int type_generator);
extern double pnl_rand_uni (int type_generator);
extern double pnl_rand_uni_ab (double a, double b, int type_generator);
extern double pnl_rand_normal (int type_generator);
extern void pnl_vect_rand_uni(PnlVect *G, int samples, double a, double b, int type_generator);
extern void pnl_vect_rand_uni_d(PnlVect *G, int dimension, double a, double b, int type_generator);
extern void pnl_vect_rand_normal(PnlVect *G, int samples, int generator);
extern void pnl_vect_rand_normal_d(PnlVect *G, int dimension, int generator);
extern void pnl_mat_rand_uni(PnlMat *M, int samples, int dimension, const PnlVect *a,
                             const PnlVect *b, int type_generator);
extern void pnl_mat_rand_uni2(PnlMat *M, int samples, int dimension,
                              double a, double b, int type_generator);
extern void pnl_mat_rand_normal(PnlMat *M, int samples, int dimension, int type_generator);
extern double pnl_rand_gamma (double a, double b, int gen);
extern double pnl_rand_chi2  (double nu, int gen);
extern int pnl_rand_bessel (double nu, double a, int gen);

/*
 * Rng interface
 */
extern void pnl_rng_free(PnlRng **);
extern PnlRng* pnl_rng_new ();
extern PnlRng* pnl_rng_create (int type);
extern PnlRng* pnl_rng_copy (const PnlRng *);
extern void pnl_rng_clone (PnlRng *, const PnlRng *);
extern void pnl_rng_init (PnlRng *rng, int type);
extern void pnl_rng_sseed(PnlRng *rng, unsigned long int s);
extern int pnl_rng_sdim (PnlRng *rng, int dim);

extern int pnl_rng_bernoulli(double p, PnlRng *rng);
extern long pnl_rng_poisson(double lambda, PnlRng *rng);
extern long pnl_rng_poisson1(double lambda,double t, PnlRng *rng);
extern double pnl_rng_exp(double lambda,PnlRng *rng);
extern double pnl_rng_dblexp (double lambda_p, double lambda_m, double p, PnlRng *rng);
extern double pnl_rng_uni (PnlRng *rng);
extern double pnl_rng_uni_ab (double a, double b, PnlRng *rng);
extern double pnl_rng_normal (PnlRng *rng);
extern double pnl_rng_lognormal (double m, double sigma2, PnlRng *rng);
extern double pnl_rng_invgauss (double mu, double lambda, PnlRng *rng);
extern double pnl_rng_gamma (double a, double b, PnlRng *rng);
extern double pnl_rng_chi2  (double nu, PnlRng *rng);
extern int pnl_rng_bessel (double nu, double a, PnlRng *rng);
extern void pnl_vect_rng_bernoulli(PnlVect *V, int samples, double a, double b, double p, PnlRng *rng);
extern void pnl_vect_rng_bernoulli_d(PnlVect *V, int dimension, const PnlVect *a, const PnlVect *b, const PnlVect *p, PnlRng *rng);
extern void pnl_vect_rng_poisson(PnlVect *V, int samples, double lambda, PnlRng *rng);
extern void pnl_vect_rng_poisson_d(PnlVect *V, int dimension, const PnlVect *lambda, PnlRng *rng);
extern void pnl_vect_rng_uni(PnlVect *G, int samples, double a, double b, PnlRng *rng);
extern void pnl_vect_rng_uni_d(PnlVect *G, int dimension, double a, double b, PnlRng *rng);
extern void pnl_vect_rng_normal(PnlVect *G, int samples, PnlRng *rng);
extern void pnl_vect_rng_normal_d(PnlVect *G, int dimension, PnlRng *rng);
extern void pnl_mat_rng_uni(PnlMat *M, int samples, int dimension, const PnlVect *a, const PnlVect *b, PnlRng *rng);
extern void pnl_mat_bernoulli(PnlMat *M, int samples, int dimension, const PnlVect *a, const PnlVect *b, const PnlVect *p, PnlRng *rng);
extern void pnl_mat_rng_poisson(PnlMat *M, int samples, int dimension, const PnlVect *lambda, PnlRng *rng);
extern void pnl_mat_rng_uni2(PnlMat *M, int samples, int dimension, double a, double b, PnlRng *rng);
extern void pnl_mat_rng_normal(PnlMat *M, int samples, int dimension, PnlRng *rng);


/*
 * States for different generators
 */
typedef struct
{
  long SEED;
  int inc1, inc2;
  long t_alea[56];
} knuth_state;

typedef struct
{
  double x10, x11, x12;
  double x20, x21, x22;
} mrgk3_state;

typedef struct
{
  double x10, x11, x12, x13, x14;
  double x20, x21, x22, x23, x24;
} mrgk5_state;

typedef struct
{
  long y;
  long t[32];       /* 32 refers to the size of a computer word */
  long x;
} shufl_state;

typedef struct
{
  long x, y, z;
  long t[32];     /* 32 is the size of a computer word */
} lecuyer_state;

typedef struct
{
  ulong u[3];
  ulong a;
} tausworthe_state;

#define PNL_DIM_MAX_QMC 300
typedef struct
{
  int prime[PNL_DIM_MAX_QMC];
  double alpha[PNL_DIM_MAX_QMC];
} sqrt_state;

typedef struct
{
  int prime[PNL_DIM_MAX_QMC];
} halton_state;

typedef struct
{
  int r;
} faure_state;

#define PNL_DIM_MAX_NIED 12
typedef struct
{
  long saut, gray;
  double facteur;
  unsigned long initial_d, initialX_n[PNL_DIM_MAX_NIED+1];
} nied_state;

#define PNL_SOBOL_I4_LOG_MAX 30
#define PNL_SOBOL_I4_DIM_MAX2 1111
typedef struct
{
  int v[PNL_SOBOL_I4_DIM_MAX2][PNL_SOBOL_I4_LOG_MAX];
  int lastq[PNL_SOBOL_I4_DIM_MAX2];
  float recipd;
  int maxcol;
} sobol_i4_state;

#define PNL_SOBOL_I8_LOG_MAX 62
#define PNL_SOBOL_I8_DIM_MAX2 1111
typedef struct
{
  long long int v[PNL_SOBOL_I8_DIM_MAX2][PNL_SOBOL_I8_LOG_MAX];
  long long int lastq[PNL_SOBOL_I8_DIM_MAX2];
  double recipd;
  long long int maxcol;
} sobol_i8_state;

/*
 * MT
 */
#define PNL_MT_N 624
typedef struct
{
  unsigned long mt[PNL_MT_N];
  int mti;
} mt_state;

extern void pnl_mt_sseed(mt_state *state, unsigned long int s);
extern unsigned long pnl_mt_genrand (mt_state *state);
extern double pnl_mt_genrand_double (mt_state *state);



/*
 * DCMT
 */
#define PNL_DCMT_N 17 /* N = p / w + 1; p = 521; w = 32; */

typedef struct
{
  ulong aaa;
  int mm,nn,rr,ww;
  ulong wmask,umask,lmask;
  int shift0, shift1, shiftB, shiftC;
  ulong maskB, maskC;
  int i;
  ulong state[PNL_DCMT_N];
} dcmt_state;

extern dcmt_state* pnl_dcmt_get_parameter(ulong seed);
extern int pnl_dcmt_create (dcmt_state*);
extern dcmt_state** pnl_dcmt_create_array(int n, ulong seed, int *count);
extern void pnl_dcmt_sseed (dcmt_state *mts, ulong s);
extern double pnl_dcmt_genrand_double(dcmt_state *mts);
extern void pnl_dcmt_free(dcmt_state **mts);
extern void pnl_dcmt_free_array(dcmt_state **mts, int count);

extern PnlRng* pnl_rng_dcmt_create_id (int id, ulong seed);
extern PnlRng** pnl_rng_dcmt_create_array (int n, ulong seed, int *count);
extern PnlRng** pnl_rng_dcmt_create_array_id (int start_id, int max_id, ulong seed, int *count);



/*@}*/

#ifdef _PNL_PRIVATE
static char pnl_rng_label[] = "PnlRng";
#endif

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_RANDOM_H */
