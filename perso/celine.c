#include <time.h>
#include "pnl/pnl_random.h"
#include "pnl/pnl_basis.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_cdf.h"
#include "pnl/pnl_optim.h"
#include "pnl/pnl_finance.h"
#include "pnl/pnl_root.h"

# define PRINT_COEFF 1

#define pi 3.14159265358979
/**
 * Characteristics of a product
 */
typedef struct _Product  Product;
struct _Product {
  double r; /*!< interest rate */
  double x0; /*!< spot */
  double T; /*!< maturity */
  double sigma; /*!< volatility */
  double K; /*!< strike */
};

/**
 * Characteristics of a Param_min
 */
typedef struct _Param_min  Param_min;
struct _Param_min {
  PnlVect *lambda_v;
  const Product *P_v;
  PnlVect *a_v;
};

double max(double a, double b)
{
  if(a>b) return a;
  else return b;
}

double gaussian(double x)
{
  return 1.0/sqrt(2*pi)*exp(-pow(x,2)/2.0);
}

double ha(double a, double x, const Product *P)
{
  return exp(-P->r*P->T)*pnl_cdfnor(-(log(x/a)+(P->r-pow(P->sigma,2)/2.0)*P->T)/(P->sigma*sqrt(P->T)));
}

double ha_prime(double a, double x, const Product *P)
{
  return exp(-P->r*P->T)*1.0/(x*P->sigma*sqrt(P->T))*gaussian(-(log(x/a)+(P->r-pow(P->sigma,2)/2.0)*P->T)/(P->sigma*sqrt(P->T)));
}

double g(double x, double K)
{
  return fmax(K-x,0);
}

double g_prime(double x, double K)
{
  if(x<K) return -1;
  else return 0;
}


/**
 * Fonction dont on cherche le zero
 */
double fonction_find_zero(double x, void *params)
{
  Param_min *PM_v = (Param_min *)params;
  PnlVect *a;
  PnlVect *lambda;
  PnlVect *H;
  const Product *P;
  int n,i;
  
  a=PM_v->a_v;
  lambda =PM_v->lambda_v;
  P=PM_v->P_v;
  n=lambda->size;
  H=pnl_vect_create_from_double(n,0);
  for(i=0;i<n;i++)
    {
      pnl_vect_set(H,i,ha_prime(pnl_vect_get(a,i),x,P));
    }
  return pnl_vect_scalar_prod(lambda,H)-g_prime(x,P->K);
}


double solve(double x0, int n, int k, PnlRng *rng, const Product *P)
{
  PnlVect *a;//vecteur d'uniformes 
  PnlVect *C;//vecteur C(i)=h_{a(i)}(x0)
  PnlVect *B;//vecteur de taille j B(i)=-g(x(i));
  PnlVect *H;//vecteur de taille n qui change avec j : H(i)=h_{a(i)}(x(j))
  PnlVect *x;//vecteur solution de taille j
  PnlMat *A;//matrice de taille j*n A(i,j)=-h_{a(i)}(x(j))
  PnlVect *Beq; //vecteur de la fonction pnl_optim_linprog
  PnlVect *xmin;
  PnlVect *lambda;
  PnlFunc func;
  Param_min PM;
  int debug = FALSE;
  double xk;
  double f_lambda, tol;
  int i,j;
  a=pnl_vect_create_from_double(n,0);
  tol=0.001;
  pnl_vect_rng_uni(a,n,0.001,100,rng);
  
  A=pnl_mat_create_from_double(k,n,0);
  x=pnl_vect_create_from_double(k+1,100);
  B=pnl_vect_create_from_double(k,0);
  C=pnl_vect_create_from_double(n,0);
  H=pnl_vect_create_from_double(n,0);
  Beq=pnl_vect_create_from_double(n,0);
  xmin = pnl_vect_create_from_double(n, PNL_NEGINF);
  lambda=pnl_vect_create_from_double(n,0);

  PM.a_v=a;
  PM.P_v=P;
  PM.lambda_v=lambda;
  func.F=fonction_find_zero;
  func.params=(void *) &PM;
  for(i=0;i<n;i++)
    {
      pnl_vect_set(C,i,ha(pnl_vect_get(a,i),x0,P));
    }
  
  for(j=1;j<k;j++)
    {
      printf("j=%d \n",j);
      for(i=0;i<n;i++)
        {
          LET(H,i) = -ha(GET(a,i),GET(x,j-1),P);
        }
      pnl_mat_resize(A,j,n);
      pnl_mat_set_row(A,H,j-1);
      pnl_vect_resize(B,j);
      LET(B,j-1) = -g(GET(x,j-1),P->K);
      pnl_vect_resize(lambda,j);
      pnl_vect_resize(Beq,j);
      pnl_optim_linprog(C,A,B,NULL,NULL,xmin,NULL,debug,lambda,&f_lambda);
      PM.lambda_v=lambda;
      printf("lambda= \n");
      pnl_vect_print_asrow(lambda);
      printf("f_lambda: %f\n", f_lambda);
      xk=pnl_root_brent(&func,0,1000,&tol);
      pnl_vect_resize(x,j+1);
      pnl_vect_set(x,j,xk);
    }
  pnl_vect_free(&a);
  pnl_vect_free(&C);
  pnl_vect_free(&B);
  pnl_vect_free(&H);
  pnl_mat_free(&A);
  pnl_vect_free(&Beq);
  pnl_vect_free(&lambda);

  return f_lambda;
}



int main()
{
  Product P;
  PnlRng *rng= pnl_rng_create (PNL_RNG_MERSENNE);
  int n,k;
  double res;
  P.r=0.06;
  P.x0=100;
  P.T=0.5;
  P.sigma=0.4;
  P.K=110;
  pnl_rng_sseed (rng, 0);
  n=100;
  k=30;
  res=solve(P.x0,n,k,rng,&P);
  pnl_rng_free(&rng);
  return res;
}
