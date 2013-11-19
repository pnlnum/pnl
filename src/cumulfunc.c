/**********************************************************************************
 *                                     DCDFLIB                                    *
 *                                                                                *
 *             Library of Fortran Routines for Cumulative Distribution            *
 *                  Functions, Inverses, and Other Parameters                     *
 *                                                                                *
 *                                 (February, 1994)                               *
 *                                                                                *
 *                                                                                *
 *                                                                                *
 *                                                                                *
 *                                                                                *
 *                                                                                *
 *                     Summary Documentation of Each Routine                      *
 *                                                                                *
 *                                                                                *
 *                                                                                *
 *                                                                                *
 *                                                                                *
 *                                                                                *
 *                                                                                *
 *                                                                                *
 *                             Compiled and Written by:                           *
 *                                                                                *
 *                                  Barry W. Brown                                *
 *                                   James Lovato                                 *
 *                                   Kathy Russell                                *
 *                                                                                *
 *                                                                                *
 *                                                                                *
 *                                                                                *
 *                                                                                *
 *                                                                                *
 *                                                                                *
 *                                                                                *
 *                                                                                *
 *                      Department of Biomathematics, Box 237                     *
 *                      The University of Texas, M.D. Anderson Cancer Center      *
 *                      1515 Holcombe Boulevard                                   *
 *                      Houston, TX      77030                                    *
 *                                                                                *
 *                                                                                *
 *  This work was supported by grant CA-16672 from the National Cancer Institute. *
 *                                                                              *
 *                                                                                *
 *                           SUMMARY OF DCDFLIB                                   *
 *                                                                                *
 * This  library  contains routines  to compute  cumulative  distribution         *
 * functions, inverses, and    parameters  of the  distribution  for  the         *
 * following set of statistical distributions:                                    *
 *                                                                                *
 *     (1) Beta                                                                   *
 *     (2) Binomial                                                               *
 *     (3) Chi-square                                                             *
 *     (4) Noncentral Chi-square                                                  *
 *     (5) F                                                                      *
 *     (6) Noncentral F                                                           *
 *     (7) Gamma                                                                  *
 *     (8) Negative Binomial                                                      *
 *     (9) Normal                                                                 *
 *     (10) Poisson                                                               *
 *     (11) Student's t                                                           *
 *                                                                                *
 * Given values of all but one parameter of a distribution, the other is          *
 * computed.  These calculations are done with  FORTRAN Double Precision          *
 * variables.                                                                     *
 *                                                                                *
 *           -------------------- WARNINGS --------------------                   *
 *                                                                                *
 * The F and  Noncentral F distribution are  not necessarily monotone  in         *
 * either degree  of  freedom argument.  Consequently,  there  may be two         *
* degree of freedom arguments that satisfy the specified condition.  An          *
* arbitrary one of these will be found by the cdf routines.                      *
*                                                                                *
* The  amount of computation  required for  the noncentral chisquare and         *
* noncentral F  distribution    is proportional  to  the  value  of  the         *
* noncentrality   parameter.  Very large values  of   this parameter can         *
* require  immense   numbers of   computation.  Consequently,  when  the         *
* noncentrality parameter is to  be calculated, the upper limit searched         *
* is 10,000.                                                                     *
*                                                                                *
*         -------------------- END WARNINGS --------------------                 *
*                                                                                *
*                             DOCUMENTATION                                      *
*                                                                                *
* This  file  contains an  overview  of the library   and is the primary         *
* documentation.                                                                 *
*                                                                                *
* Other documentation  is  in  directory 'doc'  on  the  distribution as         *
* character  (ASCII) files.  A summary  of all of the available routines         *
* is contained in dcdflib.chs (chs is an abbreviation of 'cheat sheet').         *
* The 'chs'  file  will probably be  the  primary reference.  The  file,         *
* dcdflib.fdoc, contains the  header comments for each  routine intended         *
* for direct use.                                                                *
*                                                                                *
*                              INSTALLATION                                      *
*                                                                                *
* The Fortran source routines are contained in directory src.                    *
*                                                                                *
* A  few  routines use   machine  dependent  constants.  Lists  of  such         *
* constants for different machines are found in ipmpar.f.  Uncomment the         *
* ones  appropriate to your  machine.  The distributed  version uses the         *
* IEEE arithmetic that is used by  the IBM PC,  Macintosh, and most Unix         *
* workstations.                                                                  *
*                                                                                *
* Ignore compilation warnings that lines of code are  not reachable.  We         *
* write in a Fortran structured preprocessor (FLECS)  that is similar in         *
* spirit to  Ratfor.   Sometimes our coding  practices in  FLECS lead to         *
* unreachable lines.   Also, FLECS inserts   lines of code:  STOP  "CODE         *
* FLOWING  INTO   FLECS   PROCEEDURES".  All   such    lines should   be         *
* unreachable.                                                                   *
*                                                                                *
*                                SOURCES                                         *
*                                                                                *
* The following   routines, written  by   others, are  incorporated into         *
* DCDFLIB.                                                                       *
*                                                                                *
*                           Beta Distribution                                    *
*                                                                                *
* DiDinato, A.  R. and Morris, A.  H.   Algorithm 708: Significant Digit         *
* Computation of the Incomplete Beta  Function Ratios.  ACM Trans. Math.         *
* Softw. 18 (1993), 360-373.                                                     *
*                                                                                *
*                  Gamma Distribution and It's Inverse                           *
*                                                                                *
* DiDinato, A. R. and Morris, A.  H. Computation of the Incomplete Gamma         *
* Function  Ratios and  their  Inverse.   ACM  Trans.  Math.   Softw. 12         *
* (1986), 377-393.                                                               *
*                                                                                *
*                          Normal Distribution                                   *
*                                                                                *
* Kennedy and  Gentle, Statistical Computing,  Marcel  Dekker, NY, 1980.         *
* The rational function approximations  from pages 90-95 are used during         *
* the calculation of the inverse normal.                                         *
*                                                                                *
* Cody, W.D.  (1993).  "ALGORITHM  715:  SPECFUN  -  A Portabel  FORTRAN         *
* Package   of  Special  Function   Routines   and  Test  Drivers",  acm         *
* Transactions on Mathematical Software. 19, 22-32.  A slightly modified         *
* version of Cody's function  anorm  is used for the cumultive normal.           *
*                                                                                *
*                              Zero Finder                                       *
*                                                                                *
* J.   C. P.   Bus and  T.  J.  Dekker.   Two Efficient  Algorithms with         *
* Guaranteed Convergence  for Finding a  Zero of a Function.  ACM Trans.         *
* Math. Softw. 4 (1975), 330.                                                    *
*                                                                                *
* We transliterated Algoritm R of this paper from Algol to Fortran.              *
*                                                                                *
*                           General Reference                                    *
*                                                                                *
* Abramowitz,  M. and Stegun,  I. A.  Handbook of Mathematical Functions         *
* With  Formulas, Graphs,  and   Mathematical Tables.   (1964)  National         *
* Bureau of Standards.                                                           *
*                                                                                *
* This book has been reprinted by Dover and others.                              *
*                                                                                *
*                                                                                *
*                               LEGALITIES                                       *
*                                                                                *
* Code that appeared  in an    ACM  publication  is subject  to    their         *
* algorithms policy:                                                             *
*                                                                                *
*      Submittal of  an  algorithm    for publication  in   one of   the  ACM    *
*      Transactions implies that unrestricted use  of the algorithm within  a    *
*      computer is permissible.   General permission  to copy and  distribute    *
*      the algorithm without fee is granted provided that the copies  are not    *
*      made  or   distributed for  direct   commercial  advantage.    The ACM    *
*      copyright notice and the title of the publication and its date appear,    *
*      and  notice is given that copying  is by permission of the Association    *
*      for Computing Machinery.  To copy otherwise, or to republish, requires    *
*      a fee and/or specific permission.                                         *
*                                                                                *
*      Krogh, F.  Algorithms  Policy.  ACM  Tran.   Math.  Softw.   13(1987),    *
*      183-186.                                                                  *
*                                                                                *
* We place the DCDFLIB code that we have written in the public domain.           *
*                                                                                *
*                                  NO WARRANTY                                   *
*                                                                                *
*      WE PROVIDE ABSOLUTELY  NO WARRANTY  OF ANY  KIND  EITHER  EXPRESSED OR    *
*      IMPLIED,  INCLUDING BUT   NOT LIMITED TO,  THE  IMPLIED  WARRANTIES OF    *
*      MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK    *
*      AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS  WITH YOU.  SHOULD    *
*      THIS PROGRAM PROVE  DEFECTIVE, YOU ASSUME  THE COST  OF  ALL NECESSARY    *
*      SERVICING, REPAIR OR CORRECTION.                                          *
*                                                                                *
*      IN NO  EVENT  SHALL THE UNIVERSITY  OF TEXAS OR  ANY  OF ITS COMPONENT    *
*      INSTITUTIONS INCLUDING M. D.   ANDERSON HOSPITAL BE LIABLE  TO YOU FOR    *
*      DAMAGES, INCLUDING ANY  LOST PROFITS, LOST MONIES,   OR OTHER SPECIAL,    *
*      INCIDENTAL   OR  CONSEQUENTIAL DAMAGES   ARISING   OUT  OF  THE USE OR    *
*      INABILITY TO USE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA OR    *
                         *      ITS ANALYSIS BEING  RENDERED INACCURATE OR  LOSSES SUSTAINED  BY THIRD    *
                         *      PARTIES) THE PROGRAM.                                                     *
*                                                                                *
*      (Above NO WARRANTY modified from the GNU NO WARRANTY statement.)          *
*                                                                              *
*                     HOW TO USE THE ROUTINES                                    *
*                                                                                *
* The calling sequence for each routine is of the form:                          *
*                                                                                *
*   SUBROUTINE CDF<name>( WHICH, P, Q, X, <parameters>, STATUS, BOUND ).         *
*                                                                                *
* WHICH and STATUS are INTEGER, all other arguments are DOUBLE PRECISION.        *
*                                                                                *
* <name> is a one to  three character name identifying the distribution.         *
* WHICH  is an input integer value  that identifies what parameter value         *
* is to be calculated from the values of the other parameters.                   *
*                                                                                *
* P is always the cdf evaluated at X, Q is always the compliment of the          *
* cdf evaluated at X, i.e.  1-P, and X is always the value at which the          *
* cdf  is evaluated.   The auxiliary parameters,  <parameters>,  of the          *
* distribution differ by distribution.                                           *
*                                                                                *
* If WHICH is 1, P and  Q are to be calculated, i.e., the cdf; if WHICH          *
* is 2, X is to be calculated, i.e., the inverse cdf.  The value of one          *
* auxiliary parameter in <parameters> can also be the value calculated.          *
*                                                                                *
* STATUS returns 0 if the calculation completes correctly.                       *
*                                                                                *
*            --------------------WARNING--------------------                     *
*                                                                                *
* If STATUS is not 0, no meaningful answer is returned.                          *
*                                                                                *
*         -------------------- END WARNING --------------------                  *
*                                                                                *
* STATUS returns  -I if the I'th  input parameter was  not  in the legal         *
* range (see below).  Parameters are counted  with WHICH being the first         *
* in these return values.                                                        *
*                                                                                *
* A STATUS  value of 1 indicates that  the desired answer was apparently         *
* lower than the lower bound on the search interval.  A return code of 2         *
* indicates that  the answer was  apparently higher than the upper bound         *
* on the search interval.  A return code of 3 indicates that P and Q did         *
* not sum to 1. Other positive codes are routine specific.                       *
*                                                                                *
* BOUND is not  set if STATUS is returned  as 0.  If  STATUS is -I  then         *
* BOUND is   the bound illegally  exceeded by  input  parameter I, where         *
* WHICH  is  counted as 1,  P as 2,  Q as 3,  X as 4, etc.  If STATUS is         *
* returned as 1 or 2 then BOUND  is returned as the lower or upper bound         *
* on the search interval respectively.                                           *
*                                                                                *
*                                 BOUNDS                                         *
*                                                                                *
* Below are  the rules that we used  in determining bounds on quantities         *
* to be  calculated.   Those who don't care   can find a summary  of the         *
* bounds in  dcdflib.chs.   Input bounds  are  checked for  legality  of         *
* input.  The search  range  is  the range   of values searched  for  an         *
* answer.                                                                        *
*                                                                                *
*                              Input Bounds                                      *
*                                                                                *
* Bounds on input parameters are  checked by the  CDF* routines.   These         *
* bounds were set according to the following rules.                              *
*                                                                                *
* P: If the  domain of the cdf (X) extends to  -infinity  then P must be         *
* greater than 0 otherwise P must be greater than or equal to 0.  P must         *
* always be less than or equal to 1.                                             *
*                                                                                *
* Q: If the  domain of the cdf (X) extends to  +infinity  then Q must be         *
* greater than 0 otherwise Q must be greater than or equal to 0.  Q must         *
* always be less than or equal to 1.                                             *
*                                                                                *
* Further, P and Q must sum to 1. The smaller of the two P and Q will be         *
* used in calculations to increase accuracy                                      *
*                                                                                *
* X:  If  the  domain is infinite  in   either the positive  or negative         *
* direction, no check  is performed in  that direction.  If the left end         *
* of the domain is 0, then X is checked to assure non-negativity.                *
*                                                                                *
* DF, SD, etc.:  Some auxiliary parameters must  be positive. The lowest         *
* input values accepted for these parameters is 1E-300.                          *
*                                                                                *
*                                                                                *
*                                 Search Bounds                                  *
*                                                                                *
* These are the  ranges searched for an  answer.   If the domain  of the         *
* parameter in the cdf  is closed at  some  finite value, e.g., 0,  then         *
* this value is the same endpoint of the search range.  If the domain is         *
* open  at  some finite   endpoint (which only  occurs   for  0 --  some         *
                                    * parameters must be strictly positive) then  the endpoint is 1E-300. If         *
* the  domain is infinite in either  direction then +/- 1E300 is used as         *
* the endpoint of the search range.                                              *
*                                                                                *
*                         HOW THE ROUTINES WORK                                  *
*                                                                                *
* The cumulative  distribution   functions are computed  directly.   The         *
* normal, gamma,  and  beta functions use the  code  from the references         *
* cited.  Other  cdfs are calculated  by relating them  to one  of these         *
* distributions.  For example, the  binomial and negative binomial  cdfs         *
* can be converted  to a beta cdf.   This is how fractional observations         *
* are handled.  The  formula from Abramowitz  and Stegun  for converting         *
* the cdfs is cited  in the fdoc file.    (We think the formula  for the         *
                                           * negative binomial in A&S is wrong, but there is a correct one which we         *
                                           * used.)                                                                         *
*                                                                                *
* The inverse normal and gamma are also taken  from the references.  For         *
* all other parameters, a search is made for the value that provides the         *
* desired P.  Initial  values are chosen crudely  for the search  (e.g.,         *
                                                                   * 5).  If the domain  of the cdf for the  parameter being calculated  is         *
* infinite, a step doubling strategy is  used to bound the desired value         *
* then the  zero  finder is  employed  to refine the answer.    The zero         *
* finder attempts to obtain the answer accurately to about eight decimal         *
* places.                                                                        *
**********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_cdf.h"
#include "libamos/amos.h"
#include "pnl/pnl_specfun.h"

#define PRECISION 1.0e-7

#define CHECK_WHICH(tst, b, s)    \
  if (tst)                         \
{                              \
  *bound = b; *status = s;     \
  return;                      \
}


static double algdiv(double *a, double *b);
static double alngam(double *x);
static double alnrel(double *a);
static double apser(double *a,double *b,double *x,double *eps);
static double basym(double *a,double *b,double *lambda,double *eps);
static double bcorr(double *a0,double *b0);
static double betaln(double *a0,double *b0);
static double bfrac(double *a,double *b,double *x,
                    double *y,double *lambda, double *eps);
static void bgrat(double *a,double *b,double *x,
                  double *y,double *w, double *eps,int *ierr);
static double bpser(double *a,double *b,double *x,double *eps);
static void bratio(double *a,double *b,double *x,double *y,double *w,
                   double *w1,int *ierr);
static double brcmp1(int *mu,double *a,double *b,double *x,double *y);
static double brcomp(double *a,double *b,double *x,double *y);
static double bup(double *a,double *b,double *x,double *y,int *n,double *eps);
static double devlpl(double a[],int *n,double *x);
static double dinvnr(double *p,double *q);
static void dinvr(int *status,double *x,double *fx,
                  unsigned long *qleft,unsigned long *qhi);
static void dstinv(double *zsmall,double *zbig,double *zabsst,
                   double *zrelst,double *zstpmu,double *zabsto,
                   double *zrelto);
static double dt1(double *p,double *q,double *df);
static void dzror(int *status,double *x,double *fx,double *xlo,
                  double *xhi,unsigned long *qleft,unsigned long *qhi);
static void dstzr(double *zxlo,double *zxhi,double *zabstl,double *zreltl);
static double erf1(double *x);
static double erfc1(int *ind, double *x);
static double esum(int *mu,double *x);
static double exparg(int *l);
static double fpser(double *a,double *b,double *x,double *eps);
static double gam1(double *a);
static void gaminv(double *a,double *x,double *x0,double *p,double *q,
                   int *ierr);
static double gamln(double *a);
static double gamln1(double *a);
static double Xgamm(double *a);
static void grat1(double *a,double *x,double *r,double *p,double *q,
                  double *eps);
static void gratio(double *a,double *x,double *ans,double *qans,int *ind);
static double gsumln(double *a,double *b);
static double rcomp(double *a,double *x);
static double rexp(double *x);
static double psi(double *xx);
static double rlog(double *x);
static double rlog1(double *x);
static double spmpar(int *i);
static double stvaln(double *p);
static double fifdint(double a);
static double fifdmax1(double a,double b);
static double fifdmin1(double a,double b);
static double fifdsign(double mag,double sign);
static long fifidint(double a);
static void ftnstop(char* msg);
static void cumbet(double *x,double *y,double *a,double *b,double *cum,
                   double *ccum);
static void cumbin(double *s,double *xn,double *pr,double *ompr,
                   double *cum,double *ccum);
static void cumchi(double *x,double *df,double *cum,double *ccum);
static void cumchn(double *x,double *df,double *pnonc,double *cum,
                   double *ccum);
static void cumf(double *f,double *dfn,double *dfd,double *cum,double *ccum);
static void cumfnc(double *f,double *dfn,double *dfd,double *pnonc,
                   double *cum,double *ccum);
static void cumgam(double *x,double *a,double *cum,double *ccum);
static void cumnbn(double *s,double *xn,double *pr,double *ompr,
                   double *result,double *ccum);
static void cumnor(double *arg,double *result,double *ccum);
static void cumpoi(double *s,double *xlam,double *cum,double *ccum);
static void cumt(double *t,double *df,double *cum,double *ccum);


/*
 * COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B .GE. 8
 * 
 * IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY
 * LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*M_PI) + DEL(X).
 */
static double algdiv(double *a, double *b)
{
  double c0 = .833333333333333e-01;
  double c1 = -.277777777760991e-02;
  double c2 = .793650666825390e-03;
  double c3 = -.595202931351870e-03;
  double c4 = .837308034031215e-03;
  double c5 = -.165322962780713e-02;
  static double algdiv,c,d,h,s11,s3,s5,s7,s9,t,u,v,w,x,x2,T1;
  /*
     ..
     .. Executable Statements ..
     */
  if (*a > *b) 
    {
      h = *b/ *a;
      c = 1./(1.+h);
      x = h/(1.+h);
      d = *a+(*b-0.5e0);
    }
  else
    {
      h = *a/ *b;
      c = h/(1.+h);
      x = 1./(1.+h);
      d = *b+(*a-0.5e0);
    }
  /*
     SET SN = (1 - X**N)/(1 - X)
     */
  x2 = x*x;
  s3 = 1.+(x+x2);
  s5 = 1.+(x+x2*s3);
  s7 = 1.+(x+x2*s5);
  s9 = 1.+(x+x2*s7);
  s11 = 1.+(x+x2*s9);
  /*
     SET W = DEL(B) - DEL(A + B)
     */
  t = pow(1./ *b,2.0);
  w = ((((c5*s11*t+c4*s9)*t+c3*s7)*t+c2*s5)*t+c1*s3)*t+c0;
  w *= (c/ *b);
  /*
     COMBINE THE RESULTS
     */
  T1 = *a/ *b;
  u = d*alnrel(&T1);
  v = *a*(log(*b)-1.);
  if (u > v) 
    {
      algdiv = w-v-u;
    }
  else
    {
      algdiv = w-u-v;
    }
  return algdiv;
}

/*
 * double alngam(double *x)
 * double precision LN of the GAMma function
 */
static double alngam(double *x)
{
  return pnl_sf_log_gamma (*x);
}

/*
 * -----------------------------------------------------------------------
 * EVALUATION OF THE FUNCTION LN(1 + A)
 * -----------------------------------------------------------------------
 */
static double alnrel(double *a)
{
  double p1 = -.129418923021993e+01;
  double p2 = .405303492862024e+00;
  double p3 = -.178874546012214e-01;
  double q1 = -.162752256355323e+01;
  double q2 = .747811014037616e+00;
  double q3 = -.845104217945565e-01;
  double t,t2,w,x;
  /*
     ..
     .. Executable Statements ..
     */
  if (fabs(*a) <= 0.375e0) 
    {
      t = *a/(*a+2.);
      t2 = t*t;
      w = (((p3*t2+p2)*t2+p1)*t2+1.)/(((q3*t2+q2)*t2+q1)*t2+1.);
      return 2.*t*w;
    }
  else
    {
      x = 1.e0+*a;
      return log(x);
    }
}

/*
 * -----------------------------------------------------------------------
 * APSER YIELDS THE INCOMPLETE BETA RATIO I(SUB(1-X))(B,A) FOR
 * A .LE. MIN(EPS,EPS*B), B*X .LE. 1, AND X .LE. 0.5. USED WHEN
 * A IS VERY SMALL. USE ONLY IF ABOVE INEQUALITIES ARE SATISFIED.
 * -----------------------------------------------------------------------
 */
static double apser(double *a,double *b,double *x,double *eps)
{
  double g = .577215664901533e0;
  double aj,bx,c,j,s,t,tol;
  /*
     ..
     .. Executable Statements ..
     */
  bx = *b**x;
  t = *x-bx;
  if (*b**eps <= 2.e-2)
    {
      c = log(*x)+psi(b)+g+t;
    }
  else
    {
      c = log(bx)+g+t;
    }
  tol = 5.*(*eps)*fabs(c);
  j = 1.;
  s = 0.;
  do
    {
      j += 1.;
      t *= (*x-bx/j);
      aj = t/j;
      s += aj;
    } while (fabs(aj) > tol);
  return  -(*a*(c+s));
}

/*
 * -----------------------------------------------------------------------
 * ASYMPTOTIC EXPANSION FOR IX(A,B) FOR LARGE A AND B.
 * LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED.
 * IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT
 * A AND B ARE GREATER THAN OR EQUAL TO 15.
 * -----------------------------------------------------------------------
 */
static double basym(double *a,double *b,double *lambda,double *eps)
{
  double e0 = 1.12837916709551e0;
  double e1 = .353553390593274e0;
  int num = 20;
  /*
     ------------------------
   ****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP
   ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN.
   THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1.
   ------------------------
   E0 = 2/SQRT(M_PI)
   E1 = 2**(-3/2)
   ------------------------
   */
  int K3 = 1;
  static double basym_0,bsum,dsum,f,h,h2,hn,j0,j1,r,r0,r1,s,sum,t,t0,t1,u,w,w0,z,z0,
                z2,zn,znm1;
  int i,im1,imj,j,m,mm1,mmj,n,np1;
  double a0[21],b0[21],c[21],d[21],T1,T2;
  /*
     ..
     .. Executable Statements ..
     */
  basym_0 = 0.;
  if (*a < *b) 
    {
      h = *a/ *b;
      r0 = 1./(1.+h);
      r1 = (*b-*a)/ *b;
      w0 = 1./sqrt(*a*(1.+h));
    }
  else
    {
      h = *b/ *a;
      r0 = 1./(1.+h);
      r1 = (*b-*a)/ *a;
      w0 = 1./sqrt(*b*(1.+h));
    }
  T1 = -(*lambda/ *a);
  T2 = *lambda/ *b;
  f = *a*rlog1(&T1)+*b*rlog1(&T2);
  t = exp(-f);
  if (t == 0.) return basym_0;
  z0 = sqrt(f);
  z = 0.5e0*(z0/e1);
  z2 = f+f;
  a0[0] = 2./3.*r1;
  c[0] = -(0.5e0*a0[0]);
  d[0] = -c[0];
  j0 = 0.5e0/e0*erfc1(&K3,&z0);
  j1 = e1;
  sum = j0+d[0]*w0*j1;
  s = 1.;
  h2 = h*h;
  hn = 1.;
  w = w0;
  znm1 = z;
  zn = z2;
  for(n=2; n<=num; n+=2) 
    {
      hn = h2*hn;
      a0[n-1] = 2.*r0*(1.+h*hn)/((double)n+2.);
      np1 = n+1;
      s += hn;
      a0[np1-1] = 2.*r1*s/((double)n+3.);
      for(i=n; i<=np1; i++) 
        {
          r = -(0.5e0*((double)i+1.));
          b0[0] = r*a0[0];
          for(m=2; m<=i; m++) 
            {
              bsum = 0.;
              mm1 = m-1;
              for(j=1; j<=mm1; j++) 
                {
                  mmj = m-j;
                  bsum += (((double)j*r-(double)mmj)*a0[j-1]*b0[mmj-1]);
                }
              b0[m-1] = r*a0[m-1]+bsum/(double)m;
            }
          c[i-1] = b0[i-1]/((double)i+1.);
          dsum = 0.;
          im1 = i-1;
          for(j=1; j<=im1; j++) 
            {
              imj = i-j;
              dsum += (d[imj-1]*c[j-1]);
            }
          d[i-1] = -(dsum+c[i-1]);
        }
      j0 = e1*znm1+((double)n-1.)*j0;
      j1 = e1*zn+(double)n*j1;
      znm1 = z2*znm1;
      zn = z2*zn;
      w = w0*w;
      t0 = d[n-1]*w*j0;
      w = w0*w;
      t1 = d[np1-1]*w*j1;
      sum += (t0+t1);
      if (fabs(t0)+fabs(t1) <= *eps*sum) break;
    }
  u = exp(-bcorr(a,b));
  return e0*t*u*sum;
}

/*
 * -----------------------------------------------------------------------
 * 
 * EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE
 * LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*M_PI) + DEL(A).
 * IT IS ASSUMED THAT A0 .GE. 8 AND B0 .GE. 8.
 * 
 * -----------------------------------------------------------------------
 */
static double bcorr(double *a0,double *b0)
{
  double c0 = .833333333333333e-01;
  double c1 = -.277777777760991e-02;
  double c2 = .793650666825390e-03;
  double c3 = -.595202931351870e-03;
  double c4 = .837308034031215e-03;
  double c5 = -.165322962780713e-02;
  double a,b,c,h,s11,s3,s5,s7,s9,t,w,x,x2;
  /*
     ..
     .. Executable Statements ..
     */
  a = fifdmin1(*a0,*b0);
  b = fifdmax1(*a0,*b0);
  h = a/b;
  c = h/(1.+h);
  x = 1./(1.+h);
  x2 = x*x;
  /*
     SET SN = (1 - X**N)/(1 - X)
     */
  s3 = 1.+(x+x2);
  s5 = 1.+(x+x2*s3);
  s7 = 1.+(x+x2*s5);
  s9 = 1.+(x+x2*s7);
  s11 = 1.+(x+x2*s9);
  /*
     SET W = DEL(B) - DEL(A + B)
     */
  t = pow(1./b,2.0);
  w = ((((c5*s11*t+c4*s9)*t+c3*s7)*t+c2*s5)*t+c1*s3)*t+c0;
  w *= (c/b);
  /*
     COMPUTE  DEL(A) + W
     */
  t = pow(1./a,2.0);
  return (((((c5*t+c4)*t+c3)*t+c2)*t+c1)*t+c0)/a+w;
}

/*
 * -----------------------------------------------------------------------
 * EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION
 * -----------------------------------------------------------------------
 * E = 0.5*LN(2*M_PI)
 * --------------------------
 */
static double betaln(double *a0,double *b0)
{
  double e = .918938533204673e0;
  double a,b,c,h,u,v,w,z;
  int i,n;
  double T1;
  /*
     ..
     .. Executable Statements ..
     */
  a = fifdmin1(*a0,*b0);
  b = fifdmax1(*a0,*b0);
  if (a >= 8.) 
    {
      w = bcorr(&a,&b);
      h = a/b;
      c = h/(1.+h);
      u = -((a-0.5e0)*log(c));
      v = b*alnrel(&h);
      if (u <= v) 
        return -(0.5e0*log(b))+e+w-u-v;
      else
        return -(0.5e0*log(b))+e+w-v-u;
    }
  if (a < 1.) 
    {
      if (b < 8.) 
        {
          T1 = a+b;
          return gamln(&a)+(gamln(&b)-gamln(&T1));
        }
      else
        {
          return gamln(&a)+algdiv(&a,&b);
        }
    }

  if (a > 2.) 
    {
      if (b > 1E3) 
        {
          n = a-1.;
          w = 1.;
          for(i=1; i<=n; i++) 
            {
              a -= 1.;
              w *= (a/(1.+a/b));
            }
          return log(w)-(double)n*log(b)+(gamln(&a)+algdiv(&a,&b));
        }

      n = a-1.;
      w = 1.;
      for(i=1; i<=n; i++) 
        {
          a -= 1.;
          h = a/b;
          w *= (h/(1.+h));
        }
      w = log(w);
      if (b < 8.)
        {
          n = b-1.;
          z = 1.;
          for(i=1; i<=n; i++) 
            {
              b -= 1.;
              z *= (b/(a+b));
            }
          return w+log(z)+(gamln(&a)+(gamln(&b)-gsumln(&a,&b)));
        }
      return w+gamln(&a)+algdiv(&a,&b);
    }
  if (b > 2.) 
    {
      w = 0.;
      if (b < 8.) 
        {
          n = b-1.;
          z = 1.;
          for(i=1; i<=n; i++)
            {
              b -= 1.;
              z *= (b/(a+b));
            }
          return w+log(z)+(gamln(&a)+(gamln(&b)-gsumln(&a,&b)));
        }
      return gamln(&a)+algdiv(&a,&b);
    }
  return gamln(&a)+gamln(&b)-gsumln(&a,&b);
}

/*
 * -----------------------------------------------------------------------
 * CONTINUED FRACTION EXPANSION FOR IX(A,B) WHEN A,B .GT. 1.
 * IT IS ASSUMED THAT  LAMBDA = (A + B)*Y - B.
 * -----------------------------------------------------------------------
 */
static double bfrac(double *a,double *b,double *x,double *y,double *lambda, double *eps)
{
  double bfrac_0,alpha,an,anp1,beta,bn,bnp1,c,c0,c1,e,n,p,r,r0,s,t,w,yp1;
  /*
     ..
     .. Executable Statements ..
     */
  bfrac_0 = brcomp(a,b,x,y);
  if (bfrac_0 == 0.) return bfrac_0;
  c = 1.+*lambda;
  c0 = *b/ *a;
  c1 = 1.+1./ *a;
  yp1 = *y+1.;
  n = 0.;
  p = 1.;
  s = *a+1.;
  an = 0.;
  bn = anp1 = 1.;
  bnp1 = c/c1;
  r = c1/c;
  /*
   *  CONTINUED FRACTION CALCULATION
   */
  while (1)
    {
      n += 1.;
      t = n/ *a;
      w = n*(*b-n)**x;
      e = *a/s;
      alpha = p*(p+c0)*e*e*(w**x);
      e = (1.+t)/(c1+t+t);
      beta = n+w/s+e*(c+n*yp1);
      p = 1.+t;
      s += 2.;
      /*
       * UPDATE AN, BN, ANP1, AND BNP1
       */
      t = alpha*an+beta*anp1;
      an = anp1;
      anp1 = t;
      t = alpha*bn+beta*bnp1;
      bn = bnp1;
      bnp1 = t;
      r0 = r;
      r = anp1/bnp1;
      if (fabs(r-r0) <= *eps*r) 
        {
          bfrac_0 *= r;
          return bfrac_0;
        }

      /*
       * RESCALE AN, BN, ANP1, AND BNP1
       */
      an /= bnp1;
      bn /= bnp1;
      anp1 = r;
      bnp1 = 1.;
    }
}

/*
 * -----------------------------------------------------------------------
 * ASYMPTOTIC EXPANSION FOR IX(A,B) WHEN A IS LARGER THAN B.
 * THE RESULT OF THE EXPANSION IS ADDED TO W. IT IS ASSUMED
 * THAT A .GE. 15 AND B .LE. 1.  EPS IS THE TOLERANCE USED.
 * IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
 * -----------------------------------------------------------------------
 */
static void bgrat(double *a,double *b,double *x,double *y,double *w,
                  double *eps,int *ierr)
{
  double bm1,bp2n,cn,coef,dj,j,l,lnx,n2,nu,p,q,r,s,sum,t,t2,u,v,z;
  int i,n,nm1;
  double c[30],d[30],T1;
  /*
     ..
     .. Executable Statements ..
     */
  bm1 = *b-0.5e0-0.5e0;
  nu = *a+0.5e0*bm1;
  if (*y > 0.375e0) 
    {
      lnx = log(*x);
    }
  else
    {
      T1 = -*y;
      lnx = alnrel(&T1);
    }
  z = -(nu*lnx);
  if (*b*z == 0.)  { *ierr = 1; return; }
  /*
     COMPUTATION OF THE EXPANSION
     SET R = EXP(-Z)*Z**B/GAMMA(B)
     */
  r = *b*(1.+gam1(b))*exp(*b*log(z));
  r *= (exp(*a*lnx)*exp(0.5e0*bm1*lnx));
  u = algdiv(b,a)+*b*log(nu);
  u = r*exp(-u);
  if (u == 0.) { *ierr = 1; return; }
  grat1(b,&z,&r,&p,&q,eps);
  v = 0.25e0*pow(1./nu,2.0);
  t2 = 0.25e0*lnx*lnx;
  l = *w/u;
  j = q/r;
  sum = j;
  t = cn = 1.;
  n2 = 0.;
  for(n=1; n<=30; n++) 
    {
      bp2n = *b+n2;
      j = (bp2n*(bp2n+1.)*j+(z+bp2n+1.)*t)*v;
      n2 += 2.;
      t *= t2;
      cn /= (n2*(n2+1.));
      c[n-1] = cn;
      s = 0.;
      if (n > 1) 
        {
          nm1 = n-1;
          coef = *b-(double)n;
          for(i=1; i<=nm1; i++) {
            s += (coef*c[i-1]*d[n-i-1]);
            coef += *b;
          }
        }
      d[n-1] = bm1*cn+s/(double)n;
      dj = d[n-1]*j;
      sum += dj;
      if (sum <= 0.)  { *ierr = 1; return; }
      if (fabs(dj) <= *eps*(sum+l))
        {
          *ierr = 0;
          *w += (u*sum);
          return;
        }
    }
  *ierr = 1; return; 
}

/*
 * -----------------------------------------------------------------------
 * POWER SERIES EXPANSION FOR EVALUATING IX(A,B) WHEN B .LE. 1
 * OR B*X .LE. 0.7.  EPS IS THE TOLERANCE USED.
 * -----------------------------------------------------------------------
 */
static double bpser(double *a,double *b,double *x,double *eps)
{
  double bpser_0,a0,apb,b0,c,n,sum,t,tol,u,w,z;
  int i,m;
  /*
     ..
     .. Executable Statements ..
     */
  bpser_0 = 0.;
  if (*x == 0.) return bpser_0;
  /*
     -----------------------------------------------------------------------
     COMPUTE THE FACTOR X**A/(A*BETA(A,B))
     -----------------------------------------------------------------------
     */
  a0 = fifdmin1(*a,*b);
  if (a0 < 1.) goto S10;
  z = *a*log(*x)-betaln(a,b);
  bpser_0 = exp(z)/ *a;
  goto S100;
S10:
  b0 = fifdmax1(*a,*b);
  if (b0 >= 8.) goto S90;
  if (b0 > 1.) goto S40;
  /*
     PROCEDURE FOR A0 .LT. 1 AND B0 .LE. 1
     */
  bpser_0 = pow(*x,*a);
  if (bpser_0 == 0.) return bpser_0;
  apb = *a+*b;
  if (apb > 1.)
    {
      u = *a+*b-1.e0;
      z = (1.+gam1(&u))/apb;
    }
  else
    {
      z = 1.+gam1(&apb);
    }
  c = (1.+gam1(a))*(1.+gam1(b))/z;
  bpser_0 *= (c*(*b/apb));
  goto S100;
S40:
  /*
     PROCEDURE FOR A0 .LT. 1 AND 1 .LT. B0 .LT. 8
     */
  u = gamln1(&a0);
  m = b0-1.;
  if (m >= 1) 
    {
      c = 1.;
      for(i=1; i<=m; i++) {
        b0 -= 1.;
        c *= (b0/(a0+b0));
      }
      u = log(c)+u;
    }
  z = *a*log(*x)-u;
  b0 -= 1.;
  apb = a0+b0;
  if (apb > 1.)
    {
      u = a0+b0-1.e0;
      t = (1.+gam1(&u))/apb;
    }
  else
    {
      t = 1.+gam1(&apb);
    }
  bpser_0 = exp(z)*(a0/ *a)*(1.+gam1(&b0))/t;
  goto S100;
S90:
  /*
     PROCEDURE FOR A0 .LT. 1 AND B0 .GE. 8
     */
  u = gamln1(&a0)+algdiv(&a0,&b0);
  z = *a*log(*x)-u;
  bpser_0 = a0/ *a*exp(z);
S100:
  if (bpser_0 == 0. || *a <= 0.1e0**eps) return bpser_0;
  /*
     -----------------------------------------------------------------------
     COMPUTE THE SERIES
     -----------------------------------------------------------------------
     */
  sum = n = 0.;
  c = 1.;
  tol = *eps/ *a;
  do
    {
      n += 1.;
      c *= ((0.5e0+(0.5e0-*b/n))**x);
      w = c/(*a+n);
      sum += w;
    }
  while (fabs(w) > tol) ;
  bpser_0 *= (1.+*a*sum);
  return bpser_0;
}

/*
 * -----------------------------------------------------------------------
 * 
 * EVALUATION OF THE INCOMPLETE BETA FUNCTION IX(A,B)
 * 
 * --------------------
 * 
 * IT IS ASSUMED THAT A AND B ARE NONNEGATIVE, AND THAT X .LE. 1
 * AND Y = 1 - X.  BRATIO ASSIGNS W AND W1 THE VALUES
 * 
 * W  = IX(A,B)
 * W1 = 1 - IX(A,B)
 * 
 * IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
 * IF NO INPUT ERRORS ARE DETECTED THEN IERR IS SET TO 0 AND
 * W AND W1 ARE COMPUTED. OTHERWISE, IF AN ERROR IS DETECTED,
 * THEN W AND W1 ARE ASSIGNED THE VALUE 0 AND IERR IS SET TO
 * ONE OF THE FOLLOWING VALUES ...
 * 
 * IERR = 1  IF A OR B IS NEGATIVE
 * IERR = 2  IF A = B = 0
 * IERR = 3  IF X .LT. 0 OR X .GT. 1
 * IERR = 4  IF Y .LT. 0 OR Y .GT. 1
 * IERR = 5  IF X + Y .NE. 1
 * IERR = 6  IF X = A = 0
 * IERR = 7  IF Y = B = 0
 * 
 * --------------------
 * WRITTEN BY ALFRED H. MORRIS, JR.
 * NAVAL SURFACE WARFARE CENTER
 * DAHLGREN, VIRGINIA
 * REVISED ... NOV 1991
 * -----------------------------------------------------------------------
 */
static void bratio(double *a,double *b,double *x,double *y,double *w,
                   double *w1,int *ierr)
{
  int K1 = 1;
  double a0,b0,eps,lambda,t,x0,y0,z;
  int ierr1,ind,n;
  double T2,T3,T4,T5;
  /*
     ..
     .. Executable Statements ..
     */
  /*
   ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST
   FLOATING POINT NUMBER FOR WHICH 1.0 + EPS .GT. 1.0
   */
  eps = spmpar(&K1);
  *w = *w1 = 0.;
  if (*a < 0. || *b < 0.) goto S270;
  if (*a == 0. && *b == 0.) goto S280;
  if (*x < 0. || *x > 1.) goto S290;
  if (*y < 0. || *y > 1.) goto S300;
  z = *x+*y-0.5e0-0.5e0;
  if (fabs(z) > 3.*eps) goto S310;
  *ierr = 0;
  if (*x == 0.) goto S210;
  if (*y == 0.) goto S230;
  if (*a == 0.) goto S240;
  if (*b == 0.) goto S220;
  eps = fifdmax1(eps,1.e-15);
  if (fifdmax1(*a,*b) < 1.e-3*eps) goto S260;
  ind = 0;
  a0 = *a;
  b0 = *b;
  x0 = *x;
  y0 = *y;
  if (fifdmin1(a0,b0) > 1.) goto S40;
  /*
     PROCEDURE FOR A0 .LE. 1 OR B0 .LE. 1
     */
  if (*x <= 0.5e0) goto S10;
  ind = 1;
  a0 = *b;
  b0 = *a;
  x0 = *y;
  y0 = *x;
S10:
  if (b0 < fifdmin1(eps,eps*a0)) goto S90;
  if (a0 < fifdmin1(eps,eps*b0) && b0*x0 <= 1.) goto S100;
  if (fifdmax1(a0,b0) > 1.) goto S20;
  if (a0 >= fifdmin1(0.2e0,b0)) goto S110;
  if (pow(x0,a0) <= 0.9e0) goto S110;
  if (x0 >= 0.3e0) goto S120;
  n = 20;
  goto S140;
S20:
  if (b0 <= 1.) goto S110;
  if (x0 >= 0.3e0) goto S120;
  if (x0 >= 0.1e0) goto S30;
  if (pow(x0*b0,a0) <= 0.7e0) goto S110;
S30:
  if (b0 > 15.) goto S150;
  n = 20;
  goto S140;
S40:
  /*
     PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1
     */
  if (*a > *b) goto S50;
  lambda = *a-(*a+*b)**x;
  goto S60;
S50:
  lambda = (*a+*b)**y-*b;
S60:
  if (lambda >= 0.) goto S70;
  ind = 1;
  a0 = *b;
  b0 = *a;
  x0 = *y;
  y0 = *x;
  lambda = fabs(lambda);
S70:
  if (b0 < 40. && b0*x0 <= 0.7e0) goto S110;
  if (b0 < 40.) goto S160;
  if (a0 > b0) goto S80;
  if (a0 <= 100.) goto S130;
  if (lambda > 0.03e0*a0) goto S130;
  goto S200;
S80:
  if (b0 <= 100.) goto S130;
  if (lambda > 0.03e0*b0) goto S130;
  goto S200;
S90:
  /*
     EVALUATION OF THE APPROPRIATE ALGORITHM
     */
  *w = fpser(&a0,&b0,&x0,&eps);
  *w1 = 0.5e0+(0.5e0-*w);
  goto S250;
S100:
  *w1 = apser(&a0,&b0,&x0,&eps);
  *w = 0.5e0+(0.5e0-*w1);
  goto S250;
S110:
  *w = bpser(&a0,&b0,&x0,&eps);
  *w1 = 0.5e0+(0.5e0-*w);
  goto S250;
S120:
  *w1 = bpser(&b0,&a0,&y0,&eps);
  *w = 0.5e0+(0.5e0-*w1);
  goto S250;
S130:
  T2 = 15.*eps;
  *w = bfrac(&a0,&b0,&x0,&y0,&lambda,&T2);
  *w1 = 0.5e0+(0.5e0-*w);
  goto S250;
S140:
  *w1 = bup(&b0,&a0,&y0,&x0,&n,&eps);
  b0 += (double)n;
S150:
  T3 = 15.*eps;
  bgrat(&b0,&a0,&y0,&x0,w1,&T3,&ierr1);
  *w = 0.5e0+(0.5e0-*w1);
  goto S250;
S160:
  n = b0;
  b0 -= (double)n;
  if (b0 != 0.) goto S170;
  n -= 1;
  b0 = 1.;
S170:
  *w = bup(&b0,&a0,&y0,&x0,&n,&eps);
  if (x0 > 0.7e0) goto S180;
  *w += bpser(&a0,&b0,&x0,&eps);
  *w1 = 0.5e0+(0.5e0-*w);
  goto S250;
S180:
  if (a0 > 15.) goto S190;
  n = 20;
  *w += bup(&a0,&b0,&x0,&y0,&n,&eps);
  a0 += (double)n;
S190:
  T4 = 15.*eps;
  bgrat(&a0,&b0,&x0,&y0,w,&T4,&ierr1);
  *w1 = 0.5e0+(0.5e0-*w);
  goto S250;
S200:
  T5 = 100.*eps;
  *w = basym(&a0,&b0,&lambda,&T5);
  *w1 = 0.5e0+(0.5e0-*w);
  goto S250;
S210:
  /*
     TERMINATION OF THE PROCEDURE
     */
  if (*a == 0.) goto S320;
S220:
  *w = 0.;
  *w1 = 1.;
  return;
S230:
  if (*b == 0.) goto S330;
S240:
  *w = 1.;
  *w1 = 0.;
  return;
S250:
  if (ind == 0) return;
  t = *w;
  *w = *w1;
  *w1 = t;
  return;
S260:
  /*
     PROCEDURE FOR A AND B .LT. 1.E-3*EPS
     */
  *w = *b/(*a+*b);
  *w1 = *a/(*a+*b);
  return;
S270:
  /*
     ERROR RETURN
     */
  *ierr = 1;
  return;
S280:
  *ierr = 2;
  return;
S290:
  *ierr = 3;
  return;
S300:
  *ierr = 4;
  return;
S310:
  *ierr = 5;
  return;
S320:
  *ierr = 6;
  return;
S330:
  *ierr = 7;
  return;
}

/*
 * -----------------------------------------------------------------------
 * EVALUATION OF  EXP(MU) * (X**A*Y**B/BETA(A,B))
 * -----------------------------------------------------------------------
 */
static double brcmp1(int *mu,double *a,double *b,double *x,double *y)
{
  double Const = .398942280401433e0;
  double brcmp_0,a0,apb,b0,c,e,h,lambda,lnx,lny,t,u,v,x0,y0,z;
  int i,n;
  /*
     -----------------
     = 1/SQRT(2*M_PI)
     -----------------
     */
  double T1,T2,T3,T4;
  /*
     ..
     .. Executable Statements ..
     */
  a0 = fifdmin1(*a,*b);
  if (a0 >= 8.) goto S130;
  if (*x > 0.375e0) goto S10;
  lnx = log(*x);
  T1 = -*x;
  lny = alnrel(&T1);
  goto S30;
S10:
  if (*y > 0.375e0) goto S20;
  T2 = -*y;
  lnx = alnrel(&T2);
  lny = log(*y);
  goto S30;
S20:
  lnx = log(*x);
  lny = log(*y);
S30:
  z = *a*lnx+*b*lny;
  if (a0 < 1.) goto S40;
  z -= betaln(a,b);
  return esum(mu,&z);
S40:
  /*
     -----------------------------------------------------------------------
     PROCEDURE FOR A .LT. 1 OR B .LT. 1
     -----------------------------------------------------------------------
     */
  b0 = fifdmax1(*a,*b);
  if (b0 >= 8.) goto S120;
  if (b0 > 1.) goto S70;
  /*
     ALGORITHM FOR B0 .LE. 1
     */
  brcmp_0 = esum(mu,&z);
  if (brcmp_0 == 0.) return brcmp_0;
  apb = *a+*b;
  if (apb > 1.) goto S50;
  z = 1.+gam1(&apb);
  goto S60;
S50:
  u = *a+*b-1.e0;
  z = (1.+gam1(&u))/apb;
S60:
  c = (1.+gam1(a))*(1.+gam1(b))/z;
  brcmp_0 = brcmp_0*(a0*c)/(1.+a0/b0);
  return brcmp_0;
S70:
  /*
     ALGORITHM FOR 1 .LT. B0 .LT. 8
     */
  u = gamln1(&a0);
  n = b0-1.;
  if (n < 1) goto S90;
  c = 1.;
  for(i=1; i<=n; i++) {
    b0 -= 1.;
    c *= (b0/(a0+b0));
  }
  u = log(c)+u;
S90:
  z -= u;
  b0 -= 1.;
  apb = a0+b0;
  if (apb > 1.) goto S100;
  t = 1.+gam1(&apb);
  goto S110;
S100:
  u = a0+b0-1.e0;
  t = (1.+gam1(&u))/apb;
S110:
  brcmp_0 = a0*esum(mu,&z)*(1.+gam1(&b0))/t;
  return brcmp_0;
S120:
  /*
     ALGORITHM FOR B0 .GE. 8
     */
  u = gamln1(&a0)+algdiv(&a0,&b0);
  T3 = z-u;
  return a0*esum(mu,&T3);
S130:
  /*
     -----------------------------------------------------------------------
     PROCEDURE FOR A .GE. 8 AND B .GE. 8
     -----------------------------------------------------------------------
     */
  if (*a > *b) goto S140;
  h = *a/ *b;
  x0 = h/(1.+h);
  y0 = 1./(1.+h);
  lambda = *a-(*a+*b)**x;
  goto S150;
S140:
  h = *b/ *a;
  x0 = 1./(1.+h);
  y0 = h/(1.+h);
  lambda = (*a+*b)**y-*b;
S150:
  e = -(lambda/ *a);
  if (fabs(e) > 0.6e0) goto S160;
  u = rlog1(&e);
  goto S170;
S160:
  u = e-log(*x/x0);
S170:
  e = lambda/ *b;
  if (fabs(e) > 0.6e0) goto S180;
  v = rlog1(&e);
  goto S190;
S180:
  v = e-log(*y/y0);
S190:
  T4 = -(*a*u+*b*v);
  z = esum(mu,&T4);
  return Const*sqrt(*b*x0)*z*exp(-bcorr(a,b));
}


/*
 * -----------------------------------------------------------------------
 * EVALUATION OF X**A*Y**B/BETA(A,B)
 * -----------------------------------------------------------------------
 */
static double brcomp(double *a,double *b,double *x,double *y)
{
  double Const = .398942280401433e0;
  double brcomp_0,a0,apb,b0,c,e,h,lambda,lnx,lny,t,u,v,x0,y0,z;
  int i,n;
  /*
     -----------------
     = 1/SQRT(2*M_PI)
     -----------------
     */
  double T1,T2;
  /*
     ..
     .. Executable Statements ..
     */
  brcomp_0 = 0.;
  if (*x == 0. || *y == 0.) return brcomp_0;
  a0 = fifdmin1(*a,*b);
  if (a0 >= 8.) goto S130;
  if (*x > 0.375e0) goto S10;
  lnx = log(*x);
  T1 = -*x;
  lny = alnrel(&T1);
  goto S30;
S10:
  if (*y > 0.375e0) goto S20;
  T2 = -*y;
  lnx = alnrel(&T2);
  lny = log(*y);
  goto S30;
S20:
  lnx = log(*x);
  lny = log(*y);
S30:
  z = *a*lnx+*b*lny;
  if (a0 < 1.) goto S40;
  z -= betaln(a,b);
  return exp(z);
S40:
  /*
     -----------------------------------------------------------------------
     PROCEDURE FOR A .LT. 1 OR B .LT. 1
     -----------------------------------------------------------------------
     */
  b0 = fifdmax1(*a,*b);
  if (b0 >= 8.) goto S120;
  if (b0 > 1.) goto S70;
  /*
     ALGORITHM FOR B0 .LE. 1
     */
  brcomp_0 = exp(z);
  if (brcomp_0 == 0.) return brcomp_0;
  apb = *a+*b;
  if (apb > 1.) goto S50;
  z = 1.+gam1(&apb);
  goto S60;
S50:
  u = *a+*b-1.e0;
  z = (1.+gam1(&u))/apb;
S60:
  c = (1.+gam1(a))*(1.+gam1(b))/z;
  brcomp_0 = brcomp_0*(a0*c)/(1.+a0/b0);
  return brcomp_0;
S70:
  /*
     ALGORITHM FOR 1 .LT. B0 .LT. 8
     */
  u = gamln1(&a0);
  n = b0-1.;
  if (n < 1) goto S90;
  c = 1.;
  for(i=1; i<=n; i++) {
    b0 -= 1.;
    c *= (b0/(a0+b0));
  }
  u = log(c)+u;
S90:
  z -= u;
  b0 -= 1.;
  apb = a0+b0;
  if (apb > 1.) goto S100;
  t = 1.+gam1(&apb);
  goto S110;
S100:
  u = a0+b0-1.e0;
  t = (1.+gam1(&u))/apb;
S110:
  return a0*exp(z)*(1.+gam1(&b0))/t;
S120:
  /*
     ALGORITHM FOR B0 .GE. 8
     */
  u = gamln1(&a0)+algdiv(&a0,&b0);
  return a0*exp(z-u);
S130:
  /*
     -----------------------------------------------------------------------
     PROCEDURE FOR A .GE. 8 AND B .GE. 8
     -----------------------------------------------------------------------
     */
  if (*a > *b) goto S140;
  h = *a/ *b;
  x0 = h/(1.+h);
  y0 = 1./(1.+h);
  lambda = *a-(*a+*b)**x;
  goto S150;
S140:
  h = *b/ *a;
  x0 = 1./(1.+h);
  y0 = h/(1.+h);
  lambda = (*a+*b)**y-*b;
S150:
  e = -(lambda/ *a);
  if (fabs(e) > 0.6e0) goto S160;
  u = rlog1(&e);
  goto S170;
S160:
  u = e-log(*x/x0);
S170:
  e = lambda/ *b;
  if (fabs(e) > 0.6e0) goto S180;
  v = rlog1(&e);
  goto S190;
S180:
  v = e-log(*y/y0);
S190:
  z = exp(-(*a*u+*b*v));
  return Const*sqrt(*b*x0)*z*exp(-bcorr(a,b));
}

/*
 * -----------------------------------------------------------------------
 * EVALUATION OF IX(A,B) - IX(A+N,B) WHERE N IS A POSITIVE INTEGER.
 * EPS IS THE TOLERANCE USED.
 * -----------------------------------------------------------------------
 */
static double bup(double *a,double *b,double *x,double *y,int *n,double *eps)
{
  int K1 = 1;
  int K2 = 0;
  double bup_0,ap1,apb,d,l,r,t,w;
  int i,k,kp1,mu,nm1;
  /*
     ..
     .. Executable Statements ..
     */
  /*
     OBTAIN THE SCALING FACTOR EXP(-MU) AND
     EXP(MU)*(X**A*Y**B/BETA(A,B))/A
     */
  apb = *a+*b;
  ap1 = *a+1.;
  mu = 0;
  d = 1.;
  if (*n == 1 || *a < 1.) goto S10;
  if (apb < 1.1e0*ap1) goto S10;
  mu = fabs(exparg(&K1));
  k = exparg(&K2);
  if (k < mu) mu = k;
  t = mu;
  d = exp(-t);
S10:
  bup_0 = brcmp1(&mu,a,b,x,y)/ *a;
  if (*n == 1 || bup_0 == 0.) return bup_0;
  nm1 = *n-1;
  w = d;
  /*
     LET K BE THE INDEX OF THE MAXIMUM TERM
     */
  k = 0;
  if (*b <= 1.) goto S50;
  if (*y > 1.e-4) goto S20;
  k = nm1;
  goto S30;
S20:
  r = (*b-1.)**x/ *y-*a;
  if (r < 1.) goto S50;
  k = t = nm1;
  if (r < t) k = r;
S30:
  /*
     ADD THE INCREASING TERMS OF THE SERIES
     */
  for(i=1; i<=k; i++) {
    l = i-1;
    d = (apb+l)/(ap1+l)**x*d;
    w += d;
  }
  if (k == nm1) goto S70;
S50:
  /*
     ADD THE REMAINING TERMS OF THE SERIES
     */
  kp1 = k+1;
  for(i=kp1; i<=nm1; i++) {
    l = i-1;
    d = (apb+l)/(ap1+l)**x*d;
    w += d;
    if (d <= *eps*w) goto S70;
  }
S70:
  /*
     TERMINATE THE PROCEDURE
     */
  bup_0 *= w;
  return bup_0;
}

/**
 * Cumulative Distribution Function BETa distribution
 *
 * @param which : Integer indicating  which of the  next  parameter
 values is to be calculated using values  of the others
 Legal range: 1..4
 which = 1 : Calculate p and q from x,y,a and b.
 which = 2 : Calculate x and y from p,q,a and b.
 which = 3 : Calculate a from p,q,x,y and b.
 which = 4 : Calculate b from p,q,x,y and a.
 * @param p : The integral from -infinity to x of the beta density.
 Input range: (0,1].
 * @param q : 1-p.
 Input range: (0, 1].
 p+q = 1.0.
 * @param x : Upper limit of integration of the beta
 density. Input range: [0,1]. Search range: [0,1].
 * @param y : 1-x. Input range: [0,1]. Search range: [0,1]. x + y = 1.0.
 * @param a : The first parameter of the beta density. Input range: (0, +infinity).
 * @param b : The second parameter of the beta density. Input range: (0, +infinity).
 * @param status : (0) if calculation completed correctly.
 *              (-I) if input parameter number I is out of range.
 *               (1) if answer appears to be lower than lowest
 *                 search bound.
 *               (2) if answer appears to be higher than greatest
 *                 search bound.
 *               (3) if p + q .ne. 1.
 *               (4) if x+y .ne. 1.
 *
 * @param bound : Undefined if STATUS is 0. Bound exceeded
 by parameter number I if STATUS is negative. Lower search
 bound if STATUS is 1. Upper search bound if STATUS is 2.
 * @return  any one parameter of the beta
 distribution given values for the others
 */


/*
 * 
 * void pnl_cdf_bet(int *which,double *p,double *q,double *x,double *y,
 * double *a,double *b,int *status,double *bound)
 * 
 * Cumulative Distribution Function
 * BETa Distribution
 * 
 * 
 * Function
 * 
 * 
 * Calculates any one parameter of the beta distribution given
 * values for the others.
 * 
 * 
 * Arguments
 * 
 * 
 * WHICH --> Integer indicating which of the next four argument
 * values is to be calculated from the others.
 * Legal range: 1..4
 * iwhich = 1 : Calculate P and Q from X,Y,A and B
 * iwhich = 2 : Calculate X and Y from P,Q,A and B
 * iwhich = 3 : Calculate A from P,Q,X,Y and B
 * iwhich = 4 : Calculate B from P,Q,X,Y and A
 * 
 * P <--> The integral from 0 to X of the chi-square
 * distribution.
 * Input range: [0, 1].
 * 
 * Q <--> 1-P.
 * Input range: [0, 1].
 * P + Q = 1.0.
 * 
 * X <--> Upper limit of integration of beta density.
 * Input range: [0,1].
 * Search range: [0,1]
 * 
 * Y <--> 1-X.
 * Input range: [0,1].
 * Search range: [0,1]
 * X + Y = 1.0.
 * 
 * A <--> The first parameter of the beta density.
 * Input range: (0, +infinity).
 * Search range: [1D-300,1D300]
 * 
 * B <--> The second parameter of the beta density.
 * Input range: (0, +infinity).
 * Search range: [1D-300,1D300]
 * 
 * STATUS <-- 0 if calculation completed correctly
 * -I if input parameter number I is out of range
 * 1 if answer appears to be lower than lowest
 * search bound
 * 2 if answer appears to be higher than greatest
 * search bound
 * 3 if P + Q .ne. 1
 * 4 if X + Y .ne. 1
 * 
 * BOUND <-- Undefined if STATUS is 0
 * 
 * Bound exceeded by parameter number I if STATUS
 * is negative.
 * 
 * Lower search bound if STATUS is 1.
 * 
 * Upper search bound if STATUS is 2.
 * 
 * 
* Method
* 
* 
* Cumulative distribution function  (P)  is calculated directly by
* code associated with the following reference.
* 
* DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant
* Digit Computation of the Incomplete  Beta  Function Ratios.  ACM
* Trans. Math.  Softw. 18 (1993), 360-373.
* 
* Computation of other parameters involve a seach for a value that
* produces  the desired  value  of P.   The search relies  on  the
* monotinicity of P with the other parameter.
* 
* 
* Note
* 
* 
* The beta density is proportional to
* t^(A-1) * (1-t)^(B-1)
*/
void pnl_cdf_bet(int *which,double *p,double *q,double *x,double *y,
                 double *a,double *b,int *status,double *bound)
{
  const double tol = 1.0e-8;
  const double atol = 1.0e-50;
  const double zero = 1.0e-300;
  const double inf = 1.0e300;
  const double one = 1.;
  int K1 = 1;
  double K2 = 0.;
  double K3 = 1.;
  double K8 = 0.5e0;
  double K9 = 5.;
  double fx,xhi,xlo,cum,ccum,xy,pq;
  unsigned long qhi,qleft,qporq;
  double T4,T5,T6,T7,T10,T11,T12,T13,T14,T15;
  /*
     ..
     .. Executable Statements ..
     */
  /*
     Check arguments
     */
  CHECK_WHICH (*which < 1, 1., -1);
  CHECK_WHICH (*which > 4, 4., -1);

  if (*which == 1) goto S70;
  /*
     P
     */
  if (!(*p < 0. || *p > 1.)) goto S60;
  if (!(*p < 0.)) goto S40;
  *bound = 0.;
  goto S50;
S40:
  *bound = 1.;
S50:
  *status = -2;
  return;
S70:
S60:
  if (*which == 1) goto S110;
  /*
     Q
     */
  if (!(*q < 0. || *q > 1.)) goto S100;
  if (!(*q < 0.)) goto S80;
  *bound = 0.;
  goto S90;
S80:
  *bound = 1.;
S90:
  *status = -3;
  return;
S110:
S100:
  if (*which == 2) goto S150;
  /*
     X
     */
  if (!(*x < 0. || *x > 1.)) goto S140;
  if (!(*x < 0.)) goto S120;
  *bound = 0.;
  goto S130;
S120:
  *bound = 1.;
S130:
  *status = -4;
  return;
S150:
S140:
  if (*which == 2) goto S190;
  /*
     Y
     */
  if (!(*y < 0. || *y > 1.)) goto S180;
  if (!(*y < 0.)) goto S160;
  *bound = 0.;
  goto S170;
S160:
  *bound = 1.;
S170:
  *status = -5;
  return;
S190:
S180:
  if (*which == 3) goto S210;
  /*
     A
     */
  if (!(*a <= 0.)) goto S200;
  *bound = 0.;
  *status = -6;
  return;
S210:
S200:
  if (*which == 4) goto S230;
  /*
     B
     */
  if (!(*b <= 0.)) goto S220;
  *bound = 0.;
  *status = -7;
  return;
S230:
S220:
  if (*which == 1) goto S270;
  /*
     P + Q
     */
  pq = *p+*q;
  if (!(fabs(pq-0.5e0-0.5e0) > 3.*spmpar(&K1))) goto S260;
  if (!(pq < 0.)) goto S240;
  *bound = 0.;
  goto S250;
S240:
  *bound = 1.;
S250:
  *status = 3;
  return;
S270:
S260:
  if (*which == 2) goto S310;
  /*
     X + Y
     */
  xy = *x+*y;
  if (!(fabs(xy-0.5e0-0.5e0) > 3.*spmpar(&K1))) goto S300;
  if (!(xy < 0.)) goto S280;
  *bound = 0.;
  goto S290;
S280:
  *bound = 1.;
S290:
  *status = 4;
  return;
S310:
S300:
  if (!(*which == 1)) qporq = *p <= *q;
  /*
     Select the minimum of P or Q
     Calculate ANSWERS
     */
  if (1 == *which) {
    /*
       Calculating P and Q
       */
    cumbet(x,y,a,b,p,q);
    *status = 0;
  }
  else if (2 == *which) {
    /*
       Calculating X and Y
       */
    T4 = atol;
    T5 = tol;
    dstzr(&K2,&K3,&T4,&T5);
    if (!qporq) goto S340;
    *status = 0;
    dzror(status,x,&fx,&xlo,&xhi,&qleft,&qhi);
    *y = one-*x;
S320:
    if (!(*status == 1)) goto S330;
    cumbet(x,y,a,b,&cum,&ccum);
    fx = cum-*p;
    dzror(status,x,&fx,&xlo,&xhi,&qleft,&qhi);
    *y = one-*x;
    goto S320;
S330:
    goto S370;
S340:
    *status = 0;
    dzror(status,y,&fx,&xlo,&xhi,&qleft,&qhi);
    *x = one-*y;
S350:
    if (!(*status == 1)) goto S360;
    cumbet(x,y,a,b,&cum,&ccum);
    fx = ccum-*q;
    dzror(status,y,&fx,&xlo,&xhi,&qleft,&qhi);
    *x = one-*y;
    goto S350;
S370:
S360:
    if (!(*status == -1)) goto S400;
    if (!qleft) goto S380;
    *status = 1;
    *bound = 0.;
    goto S390;
S380:
    *status = 2;
    *bound = 1.;
S400:
S390:
    ;
  }
  else if (3 == *which) {
    /*
       Computing A
       */
    *a = 5.;
    T6 = zero;
    T7 = inf;
    T10 = atol;
    T11 = tol;
    dstinv(&T6,&T7,&K8,&K8,&K9,&T10,&T11);
    *status = 0;
    dinvr(status,a,&fx,&qleft,&qhi);
S410:
    if (!(*status == 1)) goto S440;
    cumbet(x,y,a,b,&cum,&ccum);
    if (!qporq) goto S420;
    fx = cum-*p;
    goto S430;
S420:
    fx = ccum-*q;
S430:
    dinvr(status,a,&fx,&qleft,&qhi);
    goto S410;
S440:
    if (!(*status == -1)) goto S470;
    if (!qleft) goto S450;
    *status = 1;
    *bound = zero;
    goto S460;
S450:
    *status = 2;
    *bound = inf;
S470:
S460:
    ;
  }
  else if (4 == *which) {
    /*
       Computing B
       */
    *b = 5.;
    T12 = zero;
    T13 = inf;
    T14 = atol;
    T15 = tol;
    dstinv(&T12,&T13,&K8,&K8,&K9,&T14,&T15);
    *status = 0;
    dinvr(status,b,&fx,&qleft,&qhi);
S480:
    if (!(*status == 1)) goto S510;
    cumbet(x,y,a,b,&cum,&ccum);
    if (!qporq) goto S490;
    fx = cum-*p;
    goto S500;
S490:
    fx = ccum-*q;
S500:
    dinvr(status,b,&fx,&qleft,&qhi);
    goto S480;
S510:
    if (!(*status == -1)) goto S540;
    if (!qleft) goto S520;
    *status = 1;
    *bound = zero;
    goto S530;
S520:
    *status = 2;
    *bound = inf;
S530:
    ;
  }
S540:
  return;
}


/**
 * Cumulative Distribution Function
 BINa distribution
 *
 * @param which : Integer indicating  which of the  next  parameter
 values is to be calculated using values  of the others
 Legal range: 1..4.
 which = 1 : Calculate p and q from s,xn,pr and ompr.
 which = 2 : Calculate s from p,q,xn,pr and ompr.
 which = 3 : Calculate xn from p,q,s,pr and ompr.
 which = 4 : Calculate pr and ompr from p,q,s and xn.
 * @param p : The cumulation from 0 to s of the binomial distribution.
 (Probability of s or fewer successes in xn trials each
 with probability of success pr.) Input range: [0,1].
 * @param q : 1-p.
 Input range: (0, 1].
 p + q = 1.0.
 * @param s : The number of successes observed. Input range: [0, xn]. Search range: [0, xn].
 density. Input range: [0,1]. Search range: [0,1].
 * @param xn : The number of binomial trials. Input range: (0, +infinity).
 * @param pr : The probability of success in each binomial trial. Input range: [0,1]. Search range: [0,1].
 * @param ompr : 1-pr. Input range: [0,1].
 * @param status : (0) if calculation completed correctly.
 *              (-I) if input parameter number I is out of range.
 *               (1) if answer appears to be lower than lowest
 *                 search bound.
 *               (2) if answer appears to be higher than greatest
 *                 search bound.
 *               (3) if p + q .ne. 1.
 *               (4) if pr+ompr .ne. 1.
 *
 * @param bound : Undefined if STATUS is 0. Bound exceeded
 by parameter number I if STATUS is negative. Lower search
 bound if STATUS is 1. Upper search bound if STATUS is 2.
 * @return  any one parameter of the binomial
 distribution given values for the others
 */

/*
 * 
 * void pnl_cdf_bin(int *which,double *p,double *q,double *s,double *xn,
 * double *pr,double *ompr,int *status,double *bound)
 * 
 * Cumulative Distribution Function
 * BINomial distribution
 * 
 * 
 * Function
 * 
 * 
 * Calculates any one parameter of the binomial
 * distribution given values for the others.
 * 
 * 
 * Arguments
 * 
 * 
 * WHICH --> Integer indicating which of the next four argument
 * values is to be calculated from the others.
 * Legal range: 1..4
 * iwhich = 1 : Calculate P and Q from S,XN,PR and OMPR
 * iwhich = 2 : Calculate S from P,Q,XN,PR and OMPR
 * iwhich = 3 : Calculate XN from P,Q,S,PR and OMPR
 * iwhich = 4 : Calculate PR and OMPR from P,Q,S and XN
 * 
 * P <--> The cumulation from 0 to S of the binomial distribution.
 * (Probablility of S or fewer successes in XN trials each
 * with probability of success PR.)
 * Input range: [0,1].
 * 
 * Q <--> 1-P.
 * Input range: [0, 1].
 * P + Q = 1.0.
 * 
 * S <--> The number of successes observed.
 * Input range: [0, XN]
 * Search range: [0, XN]
 * 
 * XN  <--> The number of binomial trials.
 * Input range: (0, +infinity).
 * Search range: [1E-300, 1E300]
 * 
 * PR  <--> The probability of success in each binomial trial.
 * Input range: [0,1].
 * Search range: [0,1]
 * 
 * OMPR  <--> 1-PR
 * Input range: [0,1].
 * Search range: [0,1]
 * PR + OMPR = 1.0
 * 
 * STATUS <-- 0 if calculation completed correctly
 * -I if input parameter number I is out of range
 * 1 if answer appears to be lower than lowest
 * search bound
 * 2 if answer appears to be higher than greatest
 * search bound
 * 3 if P + Q .ne. 1
 * 4 if PR + OMPR .ne. 1
 * 
 * BOUND <-- Undefined if STATUS is 0
 * 
 * Bound exceeded by parameter number I if STATUS
 * is negative.
 * 
 * Lower search bound if STATUS is 1.
 * 
 * Upper search bound if STATUS is 2.
 * 
* 
* Method
* 
* 
* Formula  26.5.24    of   Abramowitz  and    Stegun,  Handbook   of
* Mathematical   Functions (1966) is   used  to reduce the  binomial
* distribution  to  the  cumulative incomplete    beta distribution.
* 
* Computation of other parameters involve a seach for a value that
* produces  the desired  value  of P.   The search relies  on  the
* monotinicity of P with the other parameter.
* 
* 
*/
void pnl_cdf_bin(int *which,double *p,double *q,double *s,double *xn,
                 double *pr,double *ompr,int *status,double *bound)
{
  const double atol = 1.0e-50;
  const double tol = 1.0e-8;
  const double zero = 1.0e-300;
  const double inf = 1.0e300;
  const double one = 1.;
  int K1 = 1;
  double K2 = 0.;
  double K3 = 0.5e0;
  double K4 = 5.;
  double K11 = 1.;
  double fx,xhi,xlo,cum,ccum,pq,prompr;
  unsigned long qhi,qleft,qporq;
  double T5,T6,T7,T8,T9,T10,T12,T13;
  /*
     ..
     .. Executable Statements ..
     */
  /*
     Check arguments
     */
  CHECK_WHICH (*which < 1, 1., -1);
  CHECK_WHICH (*which > 4, 4., -1);

  if (*which == 1) goto S70;
  /*
     P
     */
  if (!(*p < 0. || *p > 1.)) goto S60;
  if (!(*p < 0.)) goto S40;
  *bound = 0.;
  goto S50;
S40:
  *bound = 1.;
S50:
  *status = -2;
  return;
S70:
S60:
  if (*which == 1) goto S110;
  /*
     Q
     */
  if (!(*q < 0. || *q > 1.)) goto S100;
  if (!(*q < 0.)) goto S80;
  *bound = 0.;
  goto S90;
S80:
  *bound = 1.;
S90:
  *status = -3;
  return;
S110:
S100:
  if (*which == 3) goto S130;
  /*
     XN
     */
  if (!(*xn <= 0.)) goto S120;
  *bound = 0.;
  *status = -5;
  return;
S130:
S120:
  if (*which == 2) goto S170;
  /*
     S
     */
  if (!(*s < 0. || (*which != 3 && *s > *xn))) goto S160;
  if (!(*s < 0.)) goto S140;
  *bound = 0.;
  goto S150;
S140:
  *bound = *xn;
S150:
  *status = -4;
  return;
S170:
S160:
  if (*which == 4) goto S210;
  /*
     PR
     */
  if (!(*pr < 0. || *pr > 1.)) goto S200;
  if (!(*pr < 0.)) goto S180;
  *bound = 0.;
  goto S190;
S180:
  *bound = 1.;
S190:
  *status = -6;
  return;
S210:
S200:
  if (*which == 4) goto S250;
  /*
     OMPR
     */
  if (!(*ompr < 0. || *ompr > 1.)) goto S240;
  if (!(*ompr < 0.)) goto S220;
  *bound = 0.;
  goto S230;
S220:
  *bound = 1.;
S230:
  *status = -7;
  return;
S250:
S240:
  if (*which == 1) goto S290;
  /*
     P + Q
     */
  pq = *p+*q;
  if (!(fabs(pq-0.5e0-0.5e0) > 3.*spmpar(&K1))) goto S280;
  if (!(pq < 0.)) goto S260;
  *bound = 0.;
  goto S270;
S260:
  *bound = 1.;
S270:
  *status = 3;
  return;
S290:
S280:
  if (*which == 4) goto S330;
  /*
     PR + OMPR
     */
  prompr = *pr+*ompr;
  if (!(fabs(prompr-0.5e0-0.5e0) > 3.*spmpar(&K1))) goto S320;
  if (!(prompr < 0.)) goto S300;
  *bound = 0.;
  goto S310;
S300:
  *bound = 1.;
S310:
  *status = 4;
  return;
S330:
S320:
  if (!(*which == 1)) qporq = *p <= *q;
  /*
     Select the minimum of P or Q
     Calculate ANSWERS
     */
  if (1 == *which) {
    /*
       Calculating P
       */
    cumbin(s,xn,pr,ompr,p,q);
    *status = 0;
  }
  else if (2 == *which) {
    /*
       Calculating S
       */
    *s = 5.;
    T5 = atol;
    T6 = tol;
    dstinv(&K2,xn,&K3,&K3,&K4,&T5,&T6);
    *status = 0;
    dinvr(status,s,&fx,&qleft,&qhi);
S340:
    if (!(*status == 1)) goto S370;
    cumbin(s,xn,pr,ompr,&cum,&ccum);
    if (!qporq) goto S350;
    fx = cum-*p;
    goto S360;
S350:
    fx = ccum-*q;
S360:
    dinvr(status,s,&fx,&qleft,&qhi);
    goto S340;
S370:
    if (!(*status == -1)) goto S400;
    if (!qleft) goto S380;
    *status = 1;
    *bound = 0.;
    goto S390;
S380:
    *status = 2;
    *bound = *xn;
S400:
S390:
    ;
  }
  else if (3 == *which) {
    /*
       Calculating XN
       */
    *xn = 5.;
    T7 = zero;
    T8 = inf;
    T9 = atol;
    T10 = tol;
    dstinv(&T7,&T8,&K3,&K3,&K4,&T9,&T10);
    *status = 0;
    dinvr(status,xn,&fx,&qleft,&qhi);
S410:
    if (!(*status == 1)) goto S440;
    cumbin(s,xn,pr,ompr,&cum,&ccum);
    if (!qporq) goto S420;
    fx = cum-*p;
    goto S430;
S420:
    fx = ccum-*q;
S430:
    dinvr(status,xn,&fx,&qleft,&qhi);
    goto S410;
S440:
    if (!(*status == -1)) goto S470;
    if (!qleft) goto S450;
    *status = 1;
    *bound = zero;
    goto S460;
S450:
    *status = 2;
    *bound = inf;
S470:
S460:
    ;
  }
  else if (4 == *which) {
    /*
       Calculating PR and OMPR
       */
    T12 = atol;
    T13 = tol;
    dstzr(&K2,&K11,&T12,&T13);
    if (!qporq) goto S500;
    *status = 0;
    dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
    *ompr = one-*pr;
S480:
    if (!(*status == 1)) goto S490;
    cumbin(s,xn,pr,ompr,&cum,&ccum);
    fx = cum-*p;
    dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
    *ompr = one-*pr;
    goto S480;
S490:
    goto S530;
S500:
    *status = 0;
    dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
    *pr = one-*ompr;
S510:
    if (!(*status == 1)) goto S520;
    cumbin(s,xn,pr,ompr,&cum,&ccum);
    fx = ccum-*q;
    dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
    *pr = one-*ompr;
    goto S510;
S530:
S520:
    if (!(*status == -1)) goto S560;
    if (!qleft) goto S540;
    *status = 1;
    *bound = 0.;
    goto S550;
S540:
    *status = 2;
    *bound = 1.;
S550:
    ;
  }
S560:
  return;
}

/**
 * Cumulative Distribution Function
 CHI-Square distribution
 *
 * @param which : Integer indicating  which of the  next  parameter
 values is to be calculated using values  of the others
 Legal range: 1..3.
 which = 1 : Calculate p and q from x and df.
 which = 2 : Calculate x from p,q and df.
 which = 3 : Calculate df from p,q and x.
 * @param p : The integral from 0 to x of the chi-square
 distribution. Input range: [0, 1].
 * @param q : 1-p.
 Input range: (0, 1].
 p + q = 1.0.
 * @param x : Upper limit of integration of the central
 chi-square distribution.
 Input range: [0, +infinity).
 * @param df : Degrees of freedom of the chi-square distribution. Input range: (0, +infinity).
 * @param status : (0) if calculation completed correctly.
 *              (-I) if input parameter number I is out of range.
 *               (1) if answer appears to be lower than lowest
 *                 search bound.
 *               (2) if answer appears to be higher than greatest
 *                 search bound.
 *               (3) if p + q .ne. 1.
 *               (10) indicates error returned from cumgam.
 *
 * @param bound : Undefined if STATUS is 0. Bound exceeded
 by parameter number I if STATUS is negative. Lower search
 bound if STATUS is 1. Upper search bound if STATUS is 2.
 * @return  any one parameter of the chi-square
 distribution given values for the others
 */


/*
 * 
 * void pnl_cdf_chi(int *which,double *p,double *q,double *x,double *df,
 * int *status,double *bound)
 * 
 * Cumulative Distribution Function
 * CHI-Square distribution
 * 
 * 
 * Function
 * 
 * 
 * Calculates any one parameter of the chi-square
 * distribution given values for the others.
 * 
 * 
 * Arguments
 * 
 * 
 * WHICH --> Integer indicating which of the next three argument
 * values is to be calculated from the others.
 * Legal range: 1..3
 * iwhich = 1 : Calculate P and Q from X and DF
 * iwhich = 2 : Calculate X from P,Q and DF
 * iwhich = 3 : Calculate DF from P,Q and X
 * 
 * P <--> The integral from 0 to X of the chi-square
 * distribution.
 * Input range: [0, 1].
 * 
 * Q <--> 1-P.
 * Input range: (0, 1].
 * P + Q = 1.0.
 * 
 * X <--> Upper limit of integration of the non-central
 * chi-square distribution.
 * Input range: [0, +infinity).
 * Search range: [0,1E300]
 * 
 * DF <--> Degrees of freedom of the
 * chi-square distribution.
 * Input range: (0, +infinity).
 * Search range: [ 1E-300, 1E300]
 * 
 * STATUS <-- 0 if calculation completed correctly
 * -I if input parameter number I is out of range
 * 1 if answer appears to be lower than lowest
 * search bound
 * 2 if answer appears to be higher than greatest
 * search bound
 * 3 if P + Q .ne. 1
 * 10 indicates error returned from cumgam.  See
 * references in pnl_cdf_gam
 * 
 * BOUND <-- Undefined if STATUS is 0
 * 
 * Bound exceeded by parameter number I if STATUS
 * is negative.
 * 
 * Lower search bound if STATUS is 1.
 * 
 * Upper search bound if STATUS is 2.
 * 
 * 
 * Method
 * 
 * 
 * Formula    26.4.19   of Abramowitz  and     Stegun, Handbook  of
 * Mathematical Functions   (1966) is used   to reduce the chisqure
 * distribution to the incomplete distribution.
 * 
* Computation of other parameters involve a seach for a value that
* produces  the desired  value  of P.   The search relies  on  the
* monotinicity of P with the other parameter.
* 
*/
void pnl_cdf_chi(int *which,double *p,double *q,double *x,double *df,
                 int *status,double *bound)
{
  const double tol = 1.0e-8;
  const double atol = 1.0e-50;
  const double zero = 1.0e-300;
  const double inf = 1.0e300;
  int K1 = 1;
  double K2 = 0.;
  double K4 = 0.5e0;
  double K5 = 5.;
  double fx,cum,ccum,pq,porq=0.;
  unsigned long qhi,qleft,qporq;
  double T3,T6,T7,T8,T9,T10,T11;
  /*
     ..
     .. Executable Statements ..
     */
  /*
     Check arguments
     */
  CHECK_WHICH (*which < 1, 1., -1);
  CHECK_WHICH (*which > 4, 4., -1);

  if (*which == 1) goto S70;
  /*
     P
     */
  if (!(*p < 0. || *p > 1.)) goto S60;
  if (!(*p < 0.)) goto S40;
  *bound = 0.;
  goto S50;
S40:
  *bound = 1.;
S50:
  *status = -2;
  return;
S70:
S60:
  if (*which == 1) goto S110;
  /*
     Q
     */
  if (!(*q <= 0. || *q > 1.)) goto S100;
  if (!(*q <= 0.)) goto S80;
  *bound = 0.;
  goto S90;
S80:
  *bound = 1.;
S90:
  *status = -3;
  return;
S110:
S100:
  if (*which == 2) goto S130;
  /*
     X
     */
  if (!(*x < 0.)) goto S120;
  *bound = 0.;
  *status = -4;
  return;
S130:
S120:
  if (*which == 3) goto S150;
  /*
     DF
     */
  if (!(*df <= 0.)) goto S140;
  *bound = 0.;
  *status = -5;
  return;
S150:
S140:
  if (*which == 1) goto S190;
  /*
     P + Q
     */
  pq = *p+*q;
  if (!(fabs(pq-0.5e0-0.5e0) > 3.*spmpar(&K1))) goto S180;
  if (!(pq < 0.)) goto S160;
  *bound = 0.;
  goto S170;
S160:
  *bound = 1.;
S170:
  *status = 3;
  return;
S190:
S180:
  if (*which == 1) goto S220;
  /*
     Select the minimum of P or Q
     */
  qporq = *p <= *q;
  if (!qporq) goto S200;
  porq = *p;
  goto S210;
S200:
  porq = *q;
S220:
S210:
  /*
     Calculate ANSWERS
     */
  if (1 == *which) {
    /*
       Calculating P and Q
       */
    *status = 0;
    cumchi(x,df,p,q);
    if (porq > 1.5e0) {
      *status = 10;
      return;
    }
  }
  else if (2 == *which) {
    /*
       Calculating X
       */
    *x = 5.;
    T3 = inf;
    T6 = atol;
    T7 = tol;
    dstinv(&K2,&T3,&K4,&K4,&K5,&T6,&T7);
    *status = 0;
    dinvr(status,x,&fx,&qleft,&qhi);
S230:
    if (!(*status == 1)) goto S270;
    cumchi(x,df,&cum,&ccum);
    if (!qporq) goto S240;
    fx = cum-*p;
    goto S250;
S240:
    fx = ccum-*q;
S250:
    if (!(fx+porq > 1.5e0)) goto S260;
    *status = 10;
    return;
S260:
    dinvr(status,x,&fx,&qleft,&qhi);
    goto S230;
S270:
    if (!(*status == -1)) goto S300;
    if (!qleft) goto S280;
    *status = 1;
    *bound = 0.;
    goto S290;
S280:
    *status = 2;
    *bound = inf;
S300:
S290:
    ;
  }
  else if (3 == *which) {
    /*
       Calculating DF
       */
    *df = 5.;
    T8 = zero;
    T9 = inf;
    T10 = atol;
    T11 = tol;
    dstinv(&T8,&T9,&K4,&K4,&K5,&T10,&T11);
    *status = 0;
    dinvr(status,df,&fx,&qleft,&qhi);
S310:
    if (!(*status == 1)) goto S350;
    cumchi(x,df,&cum,&ccum);
    if (!qporq) goto S320;
    fx = cum-*p;
    goto S330;
S320:
    fx = ccum-*q;
S330:
    if (!(fx+porq > 1.5e0)) goto S340;
    *status = 10;
    return;
S340:
    dinvr(status,df,&fx,&qleft,&qhi);
    goto S310;
S350:
    if (!(*status == -1)) goto S380;
    if (!qleft) goto S360;
    *status = 1;
    *bound = zero;
    goto S370;
S360:
    *status = 2;
    *bound = inf;
S370:
    ;
  }
S380:
  return;
}

/**
 * Cumulative Distribution Function
 Non-central Chi-Square distribution
 *
 * @param which : Integer indicating  which of the  next  parameter
 values is to be calculated using values  of the others
 Legal range: 1..4.
 which = 1 : Calculate p and q from x and df.
 which = 2 : Calculate x from p,df and pnonc.
 which = 3 : Calculate df from p,x and pnonc.
 which = 4 : Calculate pnonc from p,x and df.
 * @param p : The integral from 0 to x of the non-central chi-square
 distribution.
 * @param q : 1-p. q is not used by this subroutine and is only included
 for similarity with other pnl_cdf_* routines.
 * @param x : Upper limit of integration of the non-central
 chi-square distribution.
 Input range: [0, +infinity).
 * @param df : Degrees of freedom of the
 non-central chi-square distribution. Input range: (0,
 +infinity).
 *@param pnonc : Non-centrality parameter of the non-central chi-square distribution.
 Input range: [0, +infinity).
 * @param status : (0) if calculation completed correctly.
 *              (-I) if input parameter number I is out of range.
 *               (1) if answer appears to be lower than lowest
 *                 search bound.
 *               (2) if answer appears to be higher than greatest
 *                 search bound.
 * @param bound : Undefined if STATUS is 0. Bound exceeded
 by parameter number I if STATUS is negative. Lower search
 bound if STATUS is 1. Upper search bound if STATUS is 2.
 * @return  any one parameter of the non-central chi-square
 distribution given values for the others
 */


/*
 * 
 * void pnl_cdf_chn(int *which,double *p,double *q,double *x,double *df,
 * double *pnonc,int *status,double *bound)
 * 
 * Cumulative Distribution Function
 * Non-central Chi-Square
 * 
 * 
 * Function
 * 
 * 
 * Calculates any one parameter of the non-central chi-square
 * distribution given values for the others.
 * 
 * 
 * Arguments
 * 
 * 
 * WHICH --> Integer indicating which of the next three argument
 * values is to be calculated from the others.
 * Input range: 1..4
 * iwhich = 1 : Calculate P and Q from X and DF
 * iwhich = 2 : Calculate X from P,DF and PNONC
 * iwhich = 3 : Calculate DF from P,X and PNONC
 * iwhich = 3 : Calculate PNONC from P,X and DF
 * 
 * P <--> The integral from 0 to X of the non-central chi-square
 * distribution.
 * Input range: [0, 1-1E-16).
 * 
 * Q <--> 1-P.
 * Q is not used by this subroutine and is only included
 * for similarity with other pnl_cdf_* routines.
 * 
 * X <--> Upper limit of integration of the non-central
 * chi-square distribution.
 * Input range: [0, +infinity).
 * Search range: [0,1E300]
 * 
 * DF <--> Degrees of freedom of the non-central
 * chi-square distribution.
 * Input range: (0, +infinity).
 * Search range: [ 1E-300, 1E300]
 * 
 * PNONC <--> Non-centrality parameter of the non-central
 * chi-square distribution.
 * Input range: [0, +infinity).
 * Search range: [0,1E4]
 * 
 * STATUS <-- 0 if calculation completed correctly
 * -I if input parameter number I is out of range
 * 1 if answer appears to be lower than lowest
 * search bound
 * 2 if answer appears to be higher than greatest
 * search bound
 * 
 * BOUND <-- Undefined if STATUS is 0
 * 
 * Bound exceeded by parameter number I if STATUS
 * is negative.
 * 
 * Lower search bound if STATUS is 1.
 * 
 * Upper search bound if STATUS is 2.
 * 
 * 
 * Method
 * 
 * 
 * Formula  26.4.25   of   Abramowitz   and   Stegun,  Handbook  of
* Mathematical  Functions (1966) is used to compute the cumulative
* distribution function.
* 
* Computation of other parameters involve a seach for a value that
* produces  the desired  value  of P.   The search relies  on  the
* monotinicity of P with the other parameter.
* 
* 
* WARNING
* 
* The computation time  required for this  routine is proportional
* to the noncentrality  parameter  (PNONC).  Very large  values of
* this parameter can consume immense  computer resources.  This is
* why the search range is bounded by 10,000.
* 
*/
void pnl_cdf_chn(int *which,double *p,double *q,double *x,double *df,
                 double *pnonc,int *status,double *bound)
{
  const double tent4 = 1.0e4;
  const double tol = 1.0e-8;
  const double atol = 1.0e-50;
  const double zero = 1.0e-300;
  const double one = 1.-1.0e-16;
  const double inf = 1.0e300;
  double K1 = 0.;
  double K3 = 0.5e0;
  double K4 = 5.;
  double fx,cum,ccum;
  unsigned long qhi,qleft;
  double T2,T5,T6,T7,T8,T9,T10,T11,T12,T13;
  /*
     ..
     .. Executable Statements ..
     */
  /*
     Check arguments
     */
  CHECK_WHICH (*which < 1, 1., -1);
  CHECK_WHICH (*which > 4, 4., -1);

  if (*which == 1) goto S70;
  /*
     P
     */
  if (!(*p < 0. || *p > one)) goto S60;
  if (!(*p < 0.)) goto S40;
  *bound = 0.;
  goto S50;
S40:
  *bound = one;
S50:
  *status = -2;
  return;
S70:
S60:
  if (*which == 2) goto S90;
  /*
     X
     */
  if (!(*x < 0.)) goto S80;
  *bound = 0.;
  *status = -4;
  return;
S90:
S80:
  if (*which == 3) goto S110;
  /*
     DF
     */
  if (!(*df <= 0.)) goto S100;
  *bound = 0.;
  *status = -5;
  return;
S110:
S100:
  if (*which == 4) goto S130;
  /*
     PNONC
     */
  if (!(*pnonc < 0.)) goto S120;
  *bound = 0.;
  *status = -6;
  return;
S130:
S120:
  /*
     Calculate ANSWERS
     */
  if (1 == *which) {
    /*
       Calculating P and Q
       */
    cumchn(x,df,pnonc,p,q);
    *status = 0;
  }
  else if (2 == *which) {
    /*
       Calculating X
       */
    *x = 5.;
    T2 = inf;
    T5 = atol;
    T6 = tol;
    dstinv(&K1,&T2,&K3,&K3,&K4,&T5,&T6);
    *status = 0;
    dinvr(status,x,&fx,&qleft,&qhi);
S140:
    if (!(*status == 1)) goto S150;
    cumchn(x,df,pnonc,&cum,&ccum);
    fx = cum-*p;
    dinvr(status,x,&fx,&qleft,&qhi);
    goto S140;
S150:
    if (!(*status == -1)) goto S180;
    if (!qleft) goto S160;
    *status = 1;
    *bound = 0.;
    goto S170;
S160:
    *status = 2;
    *bound = inf;
S180:
S170:
    ;
  }
  else if (3 == *which) {
    /*
       Calculating DF
       */
    *df = 5.;
    T7 = zero;
    T8 = inf;
    T9 = atol;
    T10 = tol;
    dstinv(&T7,&T8,&K3,&K3,&K4,&T9,&T10);
    *status = 0;
    dinvr(status,df,&fx,&qleft,&qhi);
S190:
    if (!(*status == 1)) goto S200;
    cumchn(x,df,pnonc,&cum,&ccum);
    fx = cum-*p;
    dinvr(status,df,&fx,&qleft,&qhi);
    goto S190;
S200:
    if (!(*status == -1)) goto S230;
    if (!qleft) goto S210;
    *status = 1;
    *bound = zero;
    goto S220;
S210:
    *status = 2;
    *bound = inf;
S230:
S220:
    ;
  }
  else if (4 == *which) {
    /*
       Calculating PNONC
       */
    *pnonc = 5.;
    T11 = tent4;
    T12 = atol;
    T13 = tol;
    dstinv(&K1,&T11,&K3,&K3,&K4,&T12,&T13);
    *status = 0;
    dinvr(status,pnonc,&fx,&qleft,&qhi);
S240:
    if (!(*status == 1)) goto S250;
    cumchn(x,df,pnonc,&cum,&ccum);
    fx = cum-*p;
    dinvr(status,pnonc,&fx,&qleft,&qhi);
    goto S240;
S250:
    if (!(*status == -1)) goto S280;
    if (!qleft) goto S260;
    *status = 1;
    *bound = zero;
    goto S270;
S260:
    *status = 2;
    *bound = tent4;
S270:
    ;
  }
S280:
  return;
}


/**
 * Cumulative Distribution Function
 F distribution
 *
 * @param which : Integer indicating  which of the  next  parameter
 values is to be calculated using values  of the others
 Legal range: 1..4.
 which = 1 : Calculate p and q from f,dfn and dfd.
 which = 2 : Calculate f from p,q,dfn and dfd.
 which = 3 : Calculate dfn from p,q,f and dfd.
 which = 4 : Calculate dfd from p,q,f and dfn.
 * @param p : The integral from 0 to f of the F density. Input range: [0,1].
 * @param q : 1-p. Input range: (0, 1].
 p + q = 1.0.
 * @param f : Upper limit of integration of the F-density.
 Input range: [0, +infinity).
 * @param dfn : Degrees of freedom of the numerator sum of squares.
 Input range: (0, +infinity).
 *@param dfd : Degrees of freedom of the denominator sum of squares.
 Input range: (0, +infinity).
 * @param status : (0) if calculation completed correctly.
 *              (-I) if input parameter number I is out of range.
 *               (1) if answer appears to be lower than lowest
 *                 search bound.
 *               (2) if answer appears to be higher than greatest
 *                 search bound.
 *               (3) if p + q .ne. 1.
 * @param bound : Undefined if STATUS is 0. Bound exceeded
 by parameter number I if STATUS is negative. Lower search
 bound if STATUS is 1. Upper search bound if STATUS is 2.
 * @return  any one parameter of the F distribution given values for the others
 */


/*
 * 
 * void pnl_cdf_f(int *which,double *p,double *q,double *f,double *dfn,
 * double *dfd,int *status,double *bound)
 * 
 * Cumulative Distribution Function
 * F distribution
 * 
 * 
 * Function
 * 
 * 
 * Calculates any one parameter of the F distribution
 * given values for the others.
 * 
 * 
 * Arguments
 * 
 * 
 * WHICH --> Integer indicating which of the next four argument
 * values is to be calculated from the others.
 * Legal range: 1..4
 * iwhich = 1 : Calculate P and Q from F,DFN and DFD
 * iwhich = 2 : Calculate F from P,Q,DFN and DFD
 * iwhich = 3 : Calculate DFN from P,Q,F and DFD
 * iwhich = 4 : Calculate DFD from P,Q,F and DFN
 * 
 * P <--> The integral from 0 to F of the f-density.
 * Input range: [0,1].
 * 
 * Q <--> 1-P.
 * Input range: (0, 1].
 * P + Q = 1.0.
 * 
 * F <--> Upper limit of integration of the f-density.
 * Input range: [0, +infinity).
 * Search range: [0,1E300]
 * 
 * DFN < --> Degrees of freedom of the numerator sum of squares.
 * Input range: (0, +infinity).
 * Search range: [ 1E-300, 1E300]
 * 
 * DFD < --> Degrees of freedom of the denominator sum of squares.
 * Input range: (0, +infinity).
 * Search range: [ 1E-300, 1E300]
 * 
 * STATUS <-- 0 if calculation completed correctly
 * -I if input parameter number I is out of range
 * 1 if answer appears to be lower than lowest
 * search bound
 * 2 if answer appears to be higher than greatest
 * search bound
 * 3 if P + Q .ne. 1
 * 
 * BOUND <-- Undefined if STATUS is 0
 * 
 * Bound exceeded by parameter number I if STATUS
 * is negative.
 * 
 * Lower search bound if STATUS is 1.
 * 
 * Upper search bound if STATUS is 2.
 * 
 * 
 * Method
 * 
 * 
 * Formula   26.6.2   of   Abramowitz   and   Stegun,  Handbook  of
 * Mathematical  Functions (1966) is used to reduce the computation
 * of the  cumulative  distribution function for the  F  variate to
 * that of an incomplete beta.
* 
* Computation of other parameters involve a seach for a value that
* produces  the desired  value  of P.   The search relies  on  the
* monotinicity of P with the other parameter.
* 
* WARNING
* 
* The value of the  cumulative  F distribution is  not necessarily
* monotone in  either degrees of freedom.  There  thus may  be two
* values  that  provide a given CDF  value.   This routine assumes
* monotonicity and will find an arbitrary one of the two values.
* 
*/
void pnl_cdf_f(int *which,double *p,double *q,double *f,double *dfn,
               double *dfd,int *status,double *bound)
{
  const double tol = 1.0e-8;
  const double atol = 1.0e-50;
  const double zero = 1.0e-300;
  const double inf = 1.0e300;
  int K1 = 1;
  double K2 = 0.;
  double K4 = 0.5e0;
  double K5 = 5.;
  double pq,fx,cum,ccum;
  unsigned long qhi,qleft,qporq;
  double T3,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15;
  /*
     ..
     .. Executable Statements ..
     */
  /*
     Check arguments
     */
  CHECK_WHICH (*which < 1, 1., -1);
  CHECK_WHICH (*which > 4, 4., -1);

  if (*which == 1) goto S70;
  /*
     P
     */
  if (!(*p < 0. || *p > 1.)) goto S60;
  if (!(*p < 0.)) goto S40;
  *bound = 0.;
  goto S50;
S40:
  *bound = 1.;
S50:
  *status = -2;
  return;
S70:
S60:
  if (*which == 1) goto S110;
  /*
     Q
     */
  if (!(*q <= 0. || *q > 1.)) goto S100;
  if (!(*q <= 0.)) goto S80;
  *bound = 0.;
  goto S90;
S80:
  *bound = 1.;
S90:
  *status = -3;
  return;
S110:
S100:
  if (*which == 2) goto S130;
  /*
     F
     */
  if (!(*f < 0.)) goto S120;
  *bound = 0.;
  *status = -4;
  return;
S130:
S120:
  if (*which == 3) goto S150;
  /*
     DFN
     */
  if (!(*dfn <= 0.)) goto S140;
  *bound = 0.;
  *status = -5;
  return;
S150:
S140:
  if (*which == 4) goto S170;
  /*
     DFD
     */
  if (!(*dfd <= 0.)) goto S160;
  *bound = 0.;
  *status = -6;
  return;
S170:
S160:
  if (*which == 1) goto S210;
  /*
     P + Q
     */
  pq = *p+*q;
  if (!(fabs(pq-0.5e0-0.5e0) > 3.*spmpar(&K1))) goto S200;
  if (!(pq < 0.)) goto S180;
  *bound = 0.;
  goto S190;
S180:
  *bound = 1.;
S190:
  *status = 3;
  return;
S210:
S200:
  if (!(*which == 1)) qporq = *p <= *q;
  /*
     Select the minimum of P or Q
     Calculate ANSWERS
     */
  if (1 == *which) {
    /*
       Calculating P
       */
    cumf(f,dfn,dfd,p,q);
    *status = 0;
  }
  else if (2 == *which) {
    /*
       Calculating F
       */
    *f = 5.;
    T3 = inf;
    T6 = atol;
    T7 = tol;
    dstinv(&K2,&T3,&K4,&K4,&K5,&T6,&T7);
    *status = 0;
    dinvr(status,f,&fx,&qleft,&qhi);
S220:
    if (!(*status == 1)) goto S250;
    cumf(f,dfn,dfd,&cum,&ccum);
    if (!qporq) goto S230;
    fx = cum-*p;
    goto S240;
S230:
    fx = ccum-*q;
S240:
    dinvr(status,f,&fx,&qleft,&qhi);
    goto S220;
S250:
    if (!(*status == -1)) goto S280;
    if (!qleft) goto S260;
    *status = 1;
    *bound = 0.;
    goto S270;
S260:
    *status = 2;
    *bound = inf;
S280:
S270:
    ;
  }
  else if (3 == *which) {
    /*
       Calculating DFN
       */
    *dfn = 5.;
    T8 = zero;
    T9 = inf;
    T10 = atol;
    T11 = tol;
    dstinv(&T8,&T9,&K4,&K4,&K5,&T10,&T11);
    *status = 0;
    dinvr(status,dfn,&fx,&qleft,&qhi);
S290:
    if (!(*status == 1)) goto S320;
    cumf(f,dfn,dfd,&cum,&ccum);
    if (!qporq) goto S300;
    fx = cum-*p;
    goto S310;
S300:
    fx = ccum-*q;
S310:
    dinvr(status,dfn,&fx,&qleft,&qhi);
    goto S290;
S320:
    if (!(*status == -1)) goto S350;
    if (!qleft) goto S330;
    *status = 1;
    *bound = zero;
    goto S340;
S330:
    *status = 2;
    *bound = inf;
S350:
S340:
    ;
  }
  else if (4 == *which) {
    /*
       Calculating DFD
       */
    *dfd = 5.;
    T12 = zero;
    T13 = inf;
    T14 = atol;
    T15 = tol;
    dstinv(&T12,&T13,&K4,&K4,&K5,&T14,&T15);
    *status = 0;
    dinvr(status,dfd,&fx,&qleft,&qhi);
S360:
    if (!(*status == 1)) goto S390;
    cumf(f,dfn,dfd,&cum,&ccum);
    if (!qporq) goto S370;
    fx = cum-*p;
    goto S380;
S370:
    fx = ccum-*q;
S380:
    dinvr(status,dfd,&fx,&qleft,&qhi);
    goto S360;
S390:
    if (!(*status == -1)) goto S420;
    if (!qleft) goto S400;
    *status = 1;
    *bound = zero;
    goto S410;
S400:
    *status = 2;
    *bound = inf;
S410:
    ;
  }
S420:
  return;
}


/**
 * Cumulative Distribution Function
 Non-central F distribution
 *
 * @param which : Integer indicating  which of the  next  parameter
 values is to be calculated using values  of the others
 Legal range: 1..5.
 which = 1 : Calculate p and p from f,dfn,dfd and pnonc.
 which = 2 : Calculate f from p,q,dfn,dfd and pnonc.
 which = 3 : Calculate dfn from p,q,f,dfd and pnonc.
 which = 4 : Calculate dfd from p,q,f,dfn and pnonc.
 which = 5 : Calculate pnonc from p,q,f,dfn and dfd.
 * @param p : The integral from 0 to f of the non-central F density. Input range: [0,1].
 * @param q : 1-p. Input range: (0, 1].
 p + q = 1.0.
 * @param f : Upper limit of integration of the non-central F-density.
 Input range: [0, +infinity).
 * @param dfn : Degrees of freedom of the numerator sum of squares.
 Input range: (0, +infinity).
 *@param dfd : Degrees of freedom of the denominator sum of squares.
 Input range: (0, +infinity).
 *@param pnonc : The non-centrality parameter
 Input range: [0,infinity)
 * @param status : (0) if calculation completed correctly.
 *              (-I) if input parameter number I is out of range.
 *               (1) if answer appears to be lower than lowest
 *                 search bound.
 *               (2) if answer appears to be higher than greatest
 *                 search bound.
 *               (3) if p + q .ne. 1.
 * @param bound : Undefined if STATUS is 0. Bound exceeded
 by parameter number I if STATUS is negative. Lower search
 bound if STATUS is 1. Upper search bound if STATUS is 2.
 * @return  any one parameter of the non-central F distribution given values for the others
 */

/*
 * 
 * void pnl_cdf_fnc(int *which,double *p,double *q,double *f,double *dfn,
 * double *dfd,double *phonc,int *status,double *bound)
 * 
 * Cumulative Distribution Function
 * Non-central F distribution
 * 
 * 
 * Function
 * 
 * 
 * Calculates any one parameter of the Non-central F
 * distribution given values for the others.
 * 
 * 
 * Arguments
 * 
 * 
 * WHICH --> Integer indicating which of the next five argument
 * values is to be calculated from the others.
 * Legal range: 1..5
 * iwhich = 1 : Calculate P and Q from F,DFN,DFD and PNONC
 * iwhich = 2 : Calculate F from P,Q,DFN,DFD and PNONC
 * iwhich = 3 : Calculate DFN from P,Q,F,DFD and PNONC
 * iwhich = 4 : Calculate DFD from P,Q,F,DFN and PNONC
 * iwhich = 5 : Calculate PNONC from P,Q,F,DFN and DFD
 * 
 * P <--> The integral from 0 to F of the non-central f-density.
 * Input range: [0,1-1E-16).
 * 
 * Q <--> 1-P.
 * Q is not used by this subroutine and is only included
 * for similarity with other pnl_cdf_* routines.
 * 
 * F <--> Upper limit of integration of the non-central f-density.
 * Input range: [0, +infinity).
 * Search range: [0,1E300]
 * 
 * DFN < --> Degrees of freedom of the numerator sum of squares.
 * Input range: (0, +infinity).
 * Search range: [ 1E-300, 1E300]
 * 
 * DFD < --> Degrees of freedom of the denominator sum of squares.
 * Must be in range: (0, +infinity).
 * Input range: (0, +infinity).
 * Search range: [ 1E-300, 1E300]
 * 
 * PNONC <-> The non-centrality parameter
 * Input range: [0,infinity)
 * Search range: [0,1E4]
 * 
 * STATUS <-- 0 if calculation completed correctly
 * -I if input parameter number I is out of range
 * 1 if answer appears to be lower than lowest
 * search bound
 * 2 if answer appears to be higher than greatest
 * search bound
 * 3 if P + Q .ne. 1
 * 
 * BOUND <-- Undefined if STATUS is 0
 * 
 * Bound exceeded by parameter number I if STATUS
 * is negative.
 * 
 * Lower search bound if STATUS is 1.
 * 
 * Upper search bound if STATUS is 2.
 * 
 * 
 * Method
* 
* 
* Formula  26.6.20   of   Abramowitz   and   Stegun,  Handbook  of
* Mathematical  Functions (1966) is used to compute the cumulative
* distribution function.
* 
* Computation of other parameters involve a seach for a value that
* produces  the desired  value  of P.   The search relies  on  the
* monotinicity of P with the other parameter.
* 
* WARNING
* 
* The computation time  required for this  routine is proportional
* to the noncentrality  parameter  (PNONC).  Very large  values of
* this parameter can consume immense  computer resources.  This is
* why the search range is bounded by 10,000.
* 
* WARNING
* 
* The  value  of the  cumulative  noncentral F distribution is not
* necessarily monotone in either degrees  of freedom.  There  thus
* may be two values that provide a given  CDF value.  This routine
* assumes monotonicity  and will find  an arbitrary one of the two
* values.
* 
*/
void pnl_cdf_fnc(int *which,double *p,double *q,double *f,double *dfn,
                 double *dfd,double *pnonc,int *status,double *bound)
{
  const double tent4 = 1.0e4;
  const double tol = 1.0e-13;
  const double atol = 1.0e-50;
  const double  zero = 1.0e-300;
  const double  one = 1.-1.0e-16;
  const double inf = 1.0e300;
  double K1 = 0.;
  double K3 = 0.5;
  double K4 = 5.;
  double fx,cum,ccum;
  unsigned long qhi,qleft;
  double T2,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17;

  /* Check arguments */
  CHECK_WHICH (*which < 1, 1., -1);
  CHECK_WHICH (*which > 5, 5., -1);

  if (*which == 1) goto S70;
  /* P */
  if (!(*p < 0. || *p > one)) goto S60;
  if (!(*p < 0.)) goto S40;
  *bound = 0.;
  goto S50;
S40:
  *bound = one;
S50:
  *status = -2;
  return;
S70:
S60:
  if (*which == 2) goto S90;
  /*
     F
     */
  if (!(*f < 0.)) goto S80;
  *bound = 0.;
  *status = -4;
  return;
S90:
S80:
  if (*which == 3) goto S110;
  /*
     DFN
     */
  if (!(*dfn <= 0.)) goto S100;
  *bound = 0.;
  *status = -5;
  return;
S110:
S100:
  if (*which == 4) goto S130;
  /*
     DFD
     */
  if (!(*dfd <= 0.)) goto S120;
  *bound = 0.;
  *status = -6;
  return;
S130:
S120:
  if (*which == 5) goto S150;
  /*
     PHONC
     */
  if (!(*pnonc < 0.)) goto S140;
  *bound = 0.;
  *status = -7;
  return;
S150:
S140:
  /*
     Calculate ANSWERS
     */
  if (1 == *which) {
    /*
       Calculating P
       */
    cumfnc(f,dfn,dfd,pnonc,p,q);
    *status = 0;
  }
  else if (2 == *which) {
    /*
       Calculating F
       */
    *f = 5.;
    T2 = inf;
    T5 = atol;
    T6 = tol;
    dstinv(&K1,&T2,&K3,&K3,&K4,&T5,&T6);
    *status = 0;
    dinvr(status,f,&fx,&qleft,&qhi);
S160:
    if (!(*status == 1)) goto S170;
    cumfnc(f,dfn,dfd,pnonc,&cum,&ccum);
    fx = cum-*p;
    dinvr(status,f,&fx,&qleft,&qhi);
    goto S160;
S170:
    if (!(*status == -1)) goto S200;
    if (!qleft) goto S180;
    *status = 1;
    *bound = 0.;
    goto S190;
S180:
    *status = 2;
    *bound = inf;
S200:
S190:
    ;
  }
  else if (3 == *which) {
    /*
       Calculating DFN
       */
    *dfn = 5.;
    T7 = zero;
    T8 = inf;
    T9 = atol;
    T10 = tol;
    dstinv(&T7,&T8,&K3,&K3,&K4,&T9,&T10);
    *status = 0;
    dinvr(status,dfn,&fx,&qleft,&qhi);
S210:
    if (!(*status == 1)) goto S220;
    cumfnc(f,dfn,dfd,pnonc,&cum,&ccum);
    fx = cum-*p;
    dinvr(status,dfn,&fx,&qleft,&qhi);
    goto S210;
S220:
    if (!(*status == -1)) goto S250;
    if (!qleft) goto S230;
    *status = 1;
    *bound = zero;
    goto S240;
S230:
    *status = 2;
    *bound = inf;
S250:
S240:
    ;
  }
  else if (4 == *which) {
    /*
       Calculating DFD
       */
    *dfd = 5.;
    T11 = zero;
    T12 = inf;
    T13 = atol;
    T14 = tol;
    dstinv(&T11,&T12,&K3,&K3,&K4,&T13,&T14);
    *status = 0;
    dinvr(status,dfd,&fx,&qleft,&qhi);
S260:
    if (!(*status == 1)) goto S270;
    cumfnc(f,dfn,dfd,pnonc,&cum,&ccum);
    fx = cum-*p;
    dinvr(status,dfd,&fx,&qleft,&qhi);
    goto S260;
S270:
    if (!(*status == -1)) goto S300;
    if (!qleft) goto S280;
    *status = 1;
    *bound = zero;
    goto S290;
S280:
    *status = 2;
    *bound = inf;
S300:
S290:
    ;
  }
  else if (5 == *which) {
    /*
       Calculating PHONC
       */
    *pnonc = 5.;
    T15 = tent4;
    T16 = atol;
    T17 = tol;
    dstinv(&K1,&T15,&K3,&K3,&K4,&T16,&T17);
    *status = 0;
    dinvr(status,pnonc,&fx,&qleft,&qhi);
S310:
    if (!(*status == 1)) goto S320;
    cumfnc(f,dfn,dfd,pnonc,&cum,&ccum);
    fx = cum-*p;
    dinvr(status,pnonc,&fx,&qleft,&qhi);
    goto S310;
S320:
    if (!(*status == -1)) goto S350;
    if (!qleft) goto S330;
    *status = 1;
    *bound = 0.;
    goto S340;
S330:
    *status = 2;
    *bound = tent4;
S340:
    ;
  }
S350:
  return;
}

/**
 * Cumulative Distribution Function
 GAMma distribution
 *
 * @param which : Integer indicating  which of the  next  parameter
 values is to be calculated using values  of the others
 Legal range: 1..4.
 which = 1 : Calculate p and p from x,shape and scale.
 which = 2 : Calculate x from p,q,shape and scale.
 which = 3 : Calculate shape from p,q,x and scale.
 which = 4 : Calculate scale from p,q,x and shape.
 * @param p : The integral from 0 to x of the gamma density. Input range: [0,1].
 * @param q : 1-p. Input range: (0, 1].
 p + q = 1.0.
 * @param x : Upper limit of integration of the gamma density.
 Input range: [0, +infinity).
 * @param shape : The shape parameter of the gamma density.
 Input range: (0, +infinity).
 *@param scale : The scale parameter of the gamma density.
 Input range: (0, +infinity).
 * @param status : (0) if calculation completed correctly.
 *              (-I) if input parameter number I is out of range.
 *               (1) if answer appears to be lower than lowest
 *                 search bound.
 *               (2) if answer appears to be higher than greatest
 *                 search bound.
 *               (3) if p + q .ne. 1.
 *               (10) if the gamma or inverse gamma routine cannot
 compute the answer.  Usually happens only for
 x and shape very large.
 * @param bound : Undefined if STATUS is 0. Bound exceeded
 by parameter number I if STATUS is negative. Lower search
 bound if STATUS is 1. Upper search bound if STATUS is 2.
 * @return  any one parameter of the gamma distribution given values for the others
 */


/*
 * 
 * void pnl_cdf_gam(int *which,double *p,double *q,double *x,double *shape,
 * double *scale,int *status,double *bound)
 * 
 * Cumulative Distribution Function
 * GAMma Distribution
 * 
 * 
 * Function
 * 
 * 
 * Calculates any one parameter of the gamma
 * distribution given values for the others.
 * 
 * 
 * Arguments
 * 
 * 
 * WHICH --> Integer indicating which of the next four argument
 * values is to be calculated from the others.
 * Legal range: 1..4
 * iwhich = 1 : Calculate P and Q from X,SHAPE and SCALE
 * iwhich = 2 : Calculate X from P,Q,SHAPE and SCALE
 * iwhich = 3 : Calculate SHAPE from P,Q,X and SCALE
 * iwhich = 4 : Calculate SCALE from P,Q,X and SHAPE
 * 
 * P <--> The integral from 0 to X of the gamma density.
 * Input range: [0,1].
 * 
 * Q <--> 1-P.
 * Input range: (0, 1].
 * P + Q = 1.0.
 * 
 * X <--> The upper limit of integration of the gamma density.
 * Input range: [0, +infinity).
 * Search range: [0,1E300]
 * 
 * SHAPE <--> The shape parameter of the gamma density.
 * Input range: (0, +infinity).
 * Search range: [1E-300,1E300]
 * 
 * SCALE <--> The scale parameter of the gamma density.
 * Input range: (0, +infinity).
 * Search range: (1E-300,1E300]
 * 
 * STATUS <-- 0 if calculation completed correctly
 * -I if input parameter number I is out of range
 * 1 if answer appears to be lower than lowest
 * search bound
 * 2 if answer appears to be higher than greatest
 * search bound
 * 3 if P + Q .ne. 1
 * 10 if the gamma or inverse gamma routine cannot
 * compute the answer.  Usually happens only for
 * X and SHAPE very large (gt 1E10 or more)
 * 
 * BOUND <-- Undefined if STATUS is 0
 * 
 * Bound exceeded by parameter number I if STATUS
 * is negative.
 * 
 * Lower search bound if STATUS is 1.
 * 
 * Upper search bound if STATUS is 2.
 * 
 * 
 * Method
 * 
 * 
 * Cumulative distribution function (P) is calculated directly by
* the code associated with:
* 
* DiDinato, A. R. and Morris, A. H. Computation of the  incomplete
* gamma function  ratios  and their  inverse.   ACM  Trans.  Math.
* Softw. 12 (1986), 377-393.
* 
* Computation of other parameters involve a seach for a value that
* produces  the desired  value  of P.   The search relies  on  the
* monotinicity of P with the other parameter.
* 
* 
* Note
* 
* 
* 
* The gamma density is proportional to
* T**(SHAPE - 1) * EXP(- SCALE * T)
        * 
*/
void pnl_cdf_gam(int *which,double *p,double *q,double *x,double *shape,
                 double *scale,int *status,double *bound)
{
  const double tol = 1.0e-8;
  const double atol = 1.0e-50;
  const double zero = 1.0e-300;
  const double inf = 1.0e300;
  int K1 = 1;
  double K5 = 0.5e0;
  double K6 = 5.;
  double xx,fx,xscale,cum,ccum,pq,porq=0.;
  int ierr;
  unsigned long qhi,qleft,qporq;
  double T2,T3,T4,T7,T8,T9;
  /*
     ..
     .. Executable Statements ..
     */
  /*
     Check arguments
     */
  CHECK_WHICH (*which < 1, 1., -1);
  CHECK_WHICH (*which > 4, 4., -1);

  if (*which == 1) goto S70;
  /*
     P
     */
  if (!(*p < 0. || *p > 1.)) goto S60;
  if (!(*p < 0.)) goto S40;
  *bound = 0.;
  goto S50;
S40:
  *bound = 1.;
S50:
  *status = -2;
  return;
S70:
S60:
  if (*which == 1) goto S110;
  /*
     Q
     */
  if (!(*q <= 0. || *q > 1.)) goto S100;
  if (!(*q <= 0.)) goto S80;
  *bound = 0.;
  goto S90;
S80:
  *bound = 1.;
S90:
  *status = -3;
  return;
S110:
S100:
  if (*which == 2) goto S130;
  /*
     X
     */
  if (!(*x < 0.)) goto S120;
  *bound = 0.;
  *status = -4;
  return;
S130:
S120:
  if (*which == 3) goto S150;
  /*
     SHAPE
     */
  if (!(*shape <= 0.)) goto S140;
  *bound = 0.;
  *status = -5;
  return;
S150:
S140:
  if (*which == 4) goto S170;
  /*
     SCALE
     */
  if (!(*scale <= 0.)) goto S160;
  *bound = 0.;
  *status = -6;
  return;
S170:
S160:
  if (*which == 1) goto S210;
  /*
     P + Q
     */
  pq = *p+*q;
  if (!(fabs(pq-0.5e0-0.5e0) > 3.*spmpar(&K1))) goto S200;
  if (!(pq < 0.)) goto S180;
  *bound = 0.;
  goto S190;
S180:
  *bound = 1.;
S190:
  *status = 3;
  return;
S210:
S200:
  if (*which == 1) goto S240;
  /*
     Select the minimum of P or Q
     */
  qporq = *p <= *q;
  if (!qporq) goto S220;
  porq = *p;
  goto S230;
S220:
  porq = *q;
S240:
S230:
  /*
     Calculate ANSWERS
     */
  if (1 == *which) {
    /*
       Calculating P
       */
    *status = 0;
    xscale = *x**scale;
    cumgam(&xscale,shape,p,q);
    if (porq > 1.5e0) *status = 10;
  }
  else if (2 == *which) {
    /*
       Computing X
       */
    T2 = -1.;
    gaminv(shape,&xx,&T2,p,q,&ierr);
    if (ierr < 0.) {
      *status = 10;
      return;
    }
    else  {
      *x = xx/ *scale;
      *status = 0;
    }
  }
  else if (3 == *which) {
    /*
       Computing SHAPE
       */
    *shape = 5.;
    xscale = *x**scale;
    T3 = zero;
    T4 = inf;
    T7 = atol;
    T8 = tol;
    dstinv(&T3,&T4,&K5,&K5,&K6,&T7,&T8);
    *status = 0;
    dinvr(status,shape,&fx,&qleft,&qhi);
S250:
    if (!(*status == 1)) goto S290;
    cumgam(&xscale,shape,&cum,&ccum);
    if (!qporq) goto S260;
    fx = cum-*p;
    goto S270;
S260:
    fx = ccum-*q;
S270:
    if ( !( (qporq && cum > 1.5e0) || (!qporq && ccum > 1.5e0))) goto S280;
    *status = 10;
    return;
S280:
    dinvr(status,shape,&fx,&qleft,&qhi);
    goto S250;
S290:
    if (!(*status == -1)) goto S320;
    if (!qleft) goto S300;
    *status = 1;
    *bound = zero;
    goto S310;
S300:
    *status = 2;
    *bound = inf;
S320:
S310:
    ;
  }
  else if (4 == *which) {
    /*
       Computing SCALE
       */
    T9 = -1.;
    gaminv(shape,&xx,&T9,p,q,&ierr);
    if (ierr < 0.) {
      *status = 10;
      return;
    }
    else  {
      *scale = xx/ *x;
      *status = 0;
    }
  }
  return;
}

/**
 * Cumulative Distribution Function
 Negative BiNomial distribution
 *
 * @param which : Integer indicating  which of the  next  parameter
 values is to be calculated using values  of the others
 Legal range: 1..4.
 which = 1 : Calculate p and q from s,xn,pr and ompr.
 which = 2 : Calculate s from p,q,xn,pr and ompr.
 which = 3 : Calculate xn from p,q,s,pr and ompr.
 which = 4 : Calculate pr and ompr from p,q,s and xn.
 * @param p : The cumulation from 0 to s of the  negative
 binomial distribution.
 Input range: [0,1].
 * @param q : 1-p.
 Input range: (0, 1].
 p + q = 1.0.
 * @param s : The upper limit of cumulation of the binomial distribution.
 There are s or fewer failures before the xn-th success.
 Input range: [0, +infinity).           
 * @param xn : The number of successes.
 Input range: [0, +infinity).
 * @param pr : The probability of success in each binomial trial.
 Input range: [0,1].
 * @param ompr : 1-pr
 * @param status :  (0) if calculation completed correctly.
 *              (-I) if input parameter number I is out of range.
 *               (1) if answer appears to be lower than lowest
 *                 search bound.
 *               (2) if answer appears to be higher than greatest
 *                 search bound.
 *               (3) if p + q .ne. 1.
 *               (4) if pr + ompr .ne. 1  
 *
 * @param bound : Undefined if STATUS is 0. Bound exceeded
 by parameter number I if STATUS is negative. Lower search
 bound if STATUS is 1. Upper search bound if STATUS is 2.
 * @return  any one parameter of the normal
 distribution given values for the others
 */


/*
 * 
 * void pnl_cdf_nbn(int *which,double *p,double *q,double *s,double *xn,
 * double *pr,double *ompr,int *status,double *bound)
 * 
 * Cumulative Distribution Function
 * Negative BiNomial distribution
 * 
 * 
 * Function
 * 
 * 
 * Calculates any one parameter of the negative binomial
 * distribution given values for the others.
 * 
 * The  cumulative  negative   binomial  distribution  returns  the
 * probability that there  will be  F or fewer failures before  the
 * XNth success in binomial trials each of which has probability of
 * success PR.
 * 
 * The individual term of the negative binomial is the probability of
 * S failures before XN successes and is
 * Choose( S, XN+S-1 ) * PR^(XN) * (1-PR)^S
 * 
 * 
 * Arguments
 * 
 * 
 * WHICH --> Integer indicating which of the next four argument
 * values is to be calculated from the others.
 * Legal range: 1..4
 * iwhich = 1 : Calculate P and Q from S,XN,PR and OMPR
 * iwhich = 2 : Calculate S from P,Q,XN,PR and OMPR
 * iwhich = 3 : Calculate XN from P,Q,S,PR and OMPR
 * iwhich = 4 : Calculate PR and OMPR from P,Q,S and XN
 * 
 * P <--> The cumulation from 0 to S of the  negative
 * binomial distribution.
 * Input range: [0,1].
 * 
 * Q <--> 1-P.
 * Input range: (0, 1].
 * P + Q = 1.0.
 * 
 * S <--> The upper limit of cumulation of the binomial distribution.
 * There are F or fewer failures before the XNth success.
 * Input range: [0, +infinity).
 * Search range: [0, 1E300]
 * 
 * XN  <--> The number of successes.
 * Input range: [0, +infinity).
 * Search range: [0, 1E300]
 * 
 * PR  <--> The probability of success in each binomial trial.
 * Input range: [0,1].
 * Search range: [0,1].
 * 
 * OMPR  <--> 1-PR
 * Input range: [0,1].
 * Search range: [0,1]
 * PR + OMPR = 1.0
 * 
 * STATUS <-- 0 if calculation completed correctly
 * -I if input parameter number I is out of range
 * 1 if answer appears to be lower than lowest
 * search bound
 * 2 if answer appears to be higher than greatest
 * search bound
 * 3 if P + Q .ne. 1
 * 4 if PR + OMPR .ne. 1
 * 
* BOUND <-- Undefined if STATUS is 0
* 
* Bound exceeded by parameter number I if STATUS
* is negative.
* 
* Lower search bound if STATUS is 1.
* 
* Upper search bound if STATUS is 2.
* 
* 
* Method
* 
* 
* Formula   26.5.26   of   Abramowitz  and  Stegun,  Handbook   of
* Mathematical Functions (1966) is used  to  reduce calculation of
* the cumulative distribution  function to that of  an  incomplete
* beta.
* 
* Computation of other parameters involve a seach for a value that
* produces  the desired  value  of P.   The search relies  on  the
* monotinicity of P with the other parameter.
* 
*/
void pnl_cdf_nbn(int *which,double *p,double *q,double *s,double *xn,
                 double *pr,double *ompr,int *status,double *bound)
{
  const double tol = 1.0e-8;
  const double atol = 1.0e-50;
  const double inf = 1.0e300;
  const double one = 1.;
  int K1 = 1;
  double K2 = 0.;
  double K4 = 0.5e0;
  double K5 = 5.;
  double K11 = 1.;
  double fx,xhi,xlo,pq,prompr,cum,ccum;
  unsigned long qhi,qleft,qporq;
  double T3,T6,T7,T8,T9,T10,T12,T13;
  /*
     ..
     .. Executable Statements ..
     */
  /*
     Check arguments
     */
  CHECK_WHICH (*which < 1, 1., -1);
  CHECK_WHICH (*which > 4, 4., -1);

  if (*which == 1) goto S70;
  /*
     P
     */
  if (!(*p < 0. || *p > 1.)) goto S60;
  if (!(*p < 0.)) goto S40;
  *bound = 0.;
  goto S50;
S40:
  *bound = 1.;
S50:
  *status = -2;
  return;
S70:
S60:
  if (*which == 1) goto S110;
  /*
     Q
     */
  if (!(*q <= 0. || *q > 1.)) goto S100;
  if (!(*q <= 0.)) goto S80;
  *bound = 0.;
  goto S90;
S80:
  *bound = 1.;
S90:
  *status = -3;
  return;
S110:
S100:
  if (*which == 2) goto S130;
  /*
     S
     */
  if (!(*s < 0.)) goto S120;
  *bound = 0.;
  *status = -4;
  return;
S130:
S120:
  if (*which == 3) goto S150;
  /*
     XN
     */
  if (!(*xn < 0.)) goto S140;
  *bound = 0.;
  *status = -5;
  return;
S150:
S140:
  if (*which == 4) goto S190;
  /*
     PR
     */
  if (!(*pr < 0. || *pr > 1.)) goto S180;
  if (!(*pr < 0.)) goto S160;
  *bound = 0.;
  goto S170;
S160:
  *bound = 1.;
S170:
  *status = -6;
  return;
S190:
S180:
  if (*which == 4) goto S230;
  /*
     OMPR
     */
  if (!(*ompr < 0. || *ompr > 1.)) goto S220;
  if (!(*ompr < 0.)) goto S200;
  *bound = 0.;
  goto S210;
S200:
  *bound = 1.;
S210:
  *status = -7;
  return;
S230:
S220:
  if (*which == 1) goto S270;
  /*
     P + Q
     */
  pq = *p+*q;
  if (!(fabs(pq-0.5e0-0.5e0) > 3.*spmpar(&K1))) goto S260;
  if (!(pq < 0.)) goto S240;
  *bound = 0.;
  goto S250;
S240:
  *bound = 1.;
S250:
  *status = 3;
  return;
S270:
S260:
  if (*which == 4) goto S310;
  /*
     PR + OMPR
     */
  prompr = *pr+*ompr;
  if (!(fabs(prompr-0.5e0-0.5e0) > 3.*spmpar(&K1))) goto S300;
  if (!(prompr < 0.)) goto S280;
  *bound = 0.;
  goto S290;
S280:
  *bound = 1.;
S290:
  *status = 4;
  return;
S310:
S300:
  if (!(*which == 1)) qporq = *p <= *q;
  /*
     Select the minimum of P or Q
     Calculate ANSWERS
     */
  if (1 == *which) {
    /*
       Calculating P
       */
    cumnbn(s,xn,pr,ompr,p,q);
    *status = 0;
  }
  else if (2 == *which) {
    /*
       Calculating S
       */
    *s = 5.;
    T3 = inf;
    T6 = atol;
    T7 = tol;
    dstinv(&K2,&T3,&K4,&K4,&K5,&T6,&T7);
    *status = 0;
    dinvr(status,s,&fx,&qleft,&qhi);
S320:
    if (!(*status == 1)) goto S350;
    cumnbn(s,xn,pr,ompr,&cum,&ccum);
    if (!qporq) goto S330;
    fx = cum-*p;
    goto S340;
S330:
    fx = ccum-*q;
S340:
    dinvr(status,s,&fx,&qleft,&qhi);
    goto S320;
S350:
    if (!(*status == -1)) goto S380;
    if (!qleft) goto S360;
    *status = 1;
    *bound = 0.;
    goto S370;
S360:
    *status = 2;
    *bound = inf;
S380:
S370:
    ;
  }
  else if (3 == *which) {
    /*
       Calculating XN
       */
    *xn = 5.;
    T8 = inf;
    T9 = atol;
    T10 = tol;
    dstinv(&K2,&T8,&K4,&K4,&K5,&T9,&T10);
    *status = 0;
    dinvr(status,xn,&fx,&qleft,&qhi);
S390:
    if (!(*status == 1)) goto S420;
    cumnbn(s,xn,pr,ompr,&cum,&ccum);
    if (!qporq) goto S400;
    fx = cum-*p;
    goto S410;
S400:
    fx = ccum-*q;
S410:
    dinvr(status,xn,&fx,&qleft,&qhi);
    goto S390;
S420:
    if (!(*status == -1)) goto S450;
    if (!qleft) goto S430;
    *status = 1;
    *bound = 0.;
    goto S440;
S430:
    *status = 2;
    *bound = inf;
S450:
S440:
    ;
  }
  else if (4 == *which) {
    /*
       Calculating PR and OMPR
       */
    T12 = atol;
    T13 = tol;
    dstzr(&K2,&K11,&T12,&T13);
    if (!qporq) goto S480;
    *status = 0;
    dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
    *ompr = one-*pr;
S460:
    if (!(*status == 1)) goto S470;
    cumnbn(s,xn,pr,ompr,&cum,&ccum);
    fx = cum-*p;
    dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
    *ompr = one-*pr;
    goto S460;
S470:
    goto S510;
S480:
    *status = 0;
    dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
    *pr = one-*ompr;
S490:
    if (!(*status == 1)) goto S500;
    cumnbn(s,xn,pr,ompr,&cum,&ccum);
    fx = ccum-*q;
    dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
    *pr = one-*ompr;
    goto S490;
S510:
S500:
    if (!(*status == -1)) goto S540;
    if (!qleft) goto S520;
    *status = 1;
    *bound = 0.;
    goto S530;
S520:
    *status = 2;
    *bound = 1.;
S530:
    ;
  }
S540:
  return;
}

/**
 * Cumulative Distribution Function
 NORmal distribution
 *
 * @param which : Integer indicating  which of the  next  parameter
 values is to be calculated using values  of the others
 Legal range: 1..4.
 which = 1 : Calculate p and q from x,mean and sd.
 which = 2 : Calculate x from p,q,mean and sd.
 which = 3 : Calculate mean from p,q,x and sd.
 which = 4 : Calculate sd from p,q,x and mean.
 * @param p : The integral from -infinity to x of the normal density.
 Input range: (0,1].
 * @param q : 1-p.
 Input range: (0, 1].
 p + q = 1.0.
 * @param x : Upper limit of integration of the normal-density.
 Input range: ( -infinity, +infinity)
 * @param mean : The mean of the normal density.
 * @param sd : Standard Deviation of the normal density.
 * @param status :  (0) if calculation completed correctly.
 *              (-I) if input parameter number I is out of range.
 *               (1) if answer appears to be lower than lowest
 *                 search bound.
 *               (2) if answer appears to be higher than greatest
 *                 search bound.
 *               (3) if p + q .ne. 1.
 *
 * @param bound : Undefined if STATUS is 0. Bound exceeded
 by parameter number I if STATUS is negative. Lower search
 bound if STATUS is 1. Upper search bound if STATUS is 2.
 * @return  any one parameter of the normal
 distribution given values for the others
 */

/*
 * 
 * A slightly modified version of ANORM from
 * 
 * Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
 * Package of Special Function Routines and Test Drivers"
 * acm Transactions on Mathematical Software. 19, 22-32.
 * 
 * is used to calulate the  cumulative standard normal distribution.
 * 
 * The rational functions from pages  90-95  of Kennedy and Gentle,
 * Statistical  Computing,  Marcel  Dekker, NY,  1980 are  used  as
 * starting values to Newton's Iterations which compute the inverse
 * standard normal.  Therefore no  searches  are necessary for  any
 * parameter.
 * 
 * For X < -15, the asymptotic expansion for the normal is used  as
 * the starting value in finding the inverse standard normal.
 * This is formula 26.2.12 of Abramowitz and Stegun.
 * 
 * 
 * Note
 * 
 * 
 * The normal density is proportional to
 * exp( - 0.5 * (( X - MEAN)/SD)**2)
 * 
 */
void pnl_cdf_nor(int *which,double *p,double *q,double *x,double *mean,
                 double *sd,int *status,double *bound)
{
  int K1 = 1;
  double z,pq;
  /*
     Check arguments
     */
  CHECK_WHICH(*which < 1, 1., -1);
  CHECK_WHICH(*which > 4, 4., -1);


  *status = 0;
  if (*which != 1) 
    {
      /* Check p and q */
      if (*p <= 0.) 
        {
          *bound = 0.; *status = -2;
          return;
        }
      if (*p > 1.)
        {
          *bound = 1.; *status = -2;
          return;
        }
      if (*q <= 0.) 
        {
          *bound = 1.; *status = -3;
          return;
        }
      if (*q > 1.)
        {
          *bound = 0.; *status = -3;
          return;
        }
      pq = *p+*q;
      if ( (fabs(*p + *q - 1) > 3.0*spmpar(&K1)) )
        {
          if ((pq < 0.0)) 
            {
              *bound = 0.0; *status = 3;
              return;
            }
          else
            {
              *bound = 1.0; *status = 3;
              return;
            }
        }
    }
  if ( (*which != 4) && (*sd <= 0.) )
    {
      *bound = 0.;
      *status = -6;
      return;
    }
  /*
     Calculate ANSWERS
     */
  switch (*which)
    {
    case 1:
      /* Computing P */
      z = (*x-*mean)/ *sd;
      cumnor(&z,p,q);
      break;
    case 2:
      /* Computing X */
      z = dinvnr(p,q);
      *x = *sd*z+*mean;
      break;
    case 3:
      /* Computing the MEAN */
      z = dinvnr(p,q);
      *mean = *x-*sd*z;
      break;
    case 4:
      /* Computing SD */
      z = dinvnr(p,q);
      *sd = (*x-*mean)/z;
      break;
    default:
      return;
    }
}

/**
 * Cumulative Distribution Function
 POIsson distribution
 *
 * @param which : Integer indicating  which of the  next  parameter
 * values is to be calculated using values  of the others
 *Legal range: 1..3.
 which = 1 : Calculate p and q from s and xlam.
 which = 2 : Calculate q from p,q and xlam.
 which = 3 : Calculate xlam from p,q and s.
 * @param p : The cumulation from 0 to s of the poisson density.
 Input range: [0,1].
 * @param q : 1-p.
 Input range: (0, 1].
 p + q = 1.0.
 * @param s : Upper limit of cumulation of the Poisson.
 Input range: [0, +infinity).
 * @param xlam : Mean of the Poisson distribution.
 Input range: [0, +infinity).
 * @param status :  (0) if calculation completed correctly.
 *              (-I) if input parameter number I is out of range.
 *               (1) if answer appears to be lower than lowest
 *                 search bound.
 *               (2) if answer appears to be higher than greatest
 *                 search bound.
 *               (3) if p + q .ne. 1.
 *
 * @param bound : Undefined if STATUS is 0. Bound exceeded
 by parameter number I if STATUS is negative. Lower search
 bound if STATUS is 1. Upper search bound if STATUS is 2.
 * @return  any one parameter of the Poisson
 distribution given values for the others
 */

/*
 * 
 * void pnl_cdf_poi(int *which,double *p,double *q,double *s,double *xlam,
 * int *status,double *bound)
 * 
 * Cumulative Distribution Function
 * POIsson distribution
 * 
 * 
 * Function
 * 
 * 
 * Calculates any one parameter of the Poisson
 * distribution given values for the others.
 * 
 * 
 * Arguments
 * 
 * 
 * WHICH --> Integer indicating which  argument
 * value is to be calculated from the others.
 * Legal range: 1..3
 * iwhich = 1 : Calculate P and Q from S and XLAM
 * iwhich = 2 : Calculate A from P,Q and XLAM
 * iwhich = 3 : Calculate XLAM from P,Q and S
 * 
 * P <--> The cumulation from 0 to S of the poisson density.
 * Input range: [0,1].
 * 
 * Q <--> 1-P.
 * Input range: (0, 1].
 * P + Q = 1.0.
 * 
 * S <--> Upper limit of cumulation of the Poisson.
 * Input range: [0, +infinity).
 * Search range: [0,1E300]
 * 
 * XLAM <--> Mean of the Poisson distribution.
 * Input range: [0, +infinity).
 * Search range: [0,1E300]
 * 
 * STATUS <-- 0 if calculation completed correctly
 * -I if input parameter number I is out of range
 * 1 if answer appears to be lower than lowest
 * search bound
 * 2 if answer appears to be higher than greatest
 * search bound
 * 3 if P + Q .ne. 1
 * 
 * BOUND <-- Undefined if STATUS is 0
 * 
 * Bound exceeded by parameter number I if STATUS
 * is negative.
 * 
 * Lower search bound if STATUS is 1.
 * 
 * Upper search bound if STATUS is 2.
 * 
 * 
 * Method
 * 
 * 
 * Formula   26.4.21  of   Abramowitz  and   Stegun,   Handbook  of
 * Mathematical Functions (1966) is used  to reduce the computation
 * of  the cumulative distribution function to that  of computing a
 * chi-square, hence an incomplete gamma function.
 * 
 * Cumulative  distribution function  (P) is  calculated  directly.
 * Computation of other parameters involve a seach for a value that
 * produces  the desired value of  P.   The  search relies  on  the
 * monotinicity of P with the other parameter.
* 
*/
void pnl_cdf_poi(int *which,double *p,double *q,double *s,double *xlam,
                 int *status,double *bound)
{
  const double tol = 1.0e-8;
  const double atol = 1.0e-50;
  const double inf = 1.0e300;
  int K1 = 1;
  double K2 = 0.;
  double K4 = 0.5e0;
  double K5 = 5.;
  double fx,cum,ccum,pq;
  unsigned long qhi,qleft,qporq;
  double T3,T6,T7,T8,T9,T10;
  /*
     ..
     .. Executable Statements ..
     */
  /*
     Check arguments
     */
  CHECK_WHICH (*which < 1, 1., -1);
  CHECK_WHICH (*which > 3, 3., -1);

  if (*which == 1) goto S70;
  /*
     P
     */
  if (!(*p < 0. || *p > 1.)) goto S60;
  if (!(*p < 0.)) goto S40;
  *bound = 0.;
  goto S50;
S40:
  *bound = 1.;
S50:
  *status = -2;
  return;
S70:
S60:
  if (*which == 1) goto S110;
  /*
     Q
     */
  if (!(*q <= 0. || *q > 1.)) goto S100;
  if (!(*q <= 0.)) goto S80;
  *bound = 0.;
  goto S90;
S80:
  *bound = 1.;
S90:
  *status = -3;
  return;
S110:
S100:
  if (*which == 2) goto S130;
  /*
     S
     */
  if (!(*s < 0.)) goto S120;
  *bound = 0.;
  *status = -4;
  return;
S130:
S120:
  if (*which == 3) goto S150;
  /*
     XLAM
     */
  if (!(*xlam < 0.)) goto S140;
  *bound = 0.;
  *status = -5;
  return;
S150:
S140:
  if (*which == 1) goto S190;
  /*
     P + Q
     */
  pq = *p+*q;
  if (!(fabs(pq-0.5e0-0.5e0) > 3.*spmpar(&K1))) goto S180;
  if (!(pq < 0.)) goto S160;
  *bound = 0.;
  goto S170;
S160:
  *bound = 1.;
S170:
  *status = 3;
  return;
S190:
S180:
  if (!(*which == 1)) qporq = *p <= *q;
  /*
     Select the minimum of P or Q
     Calculate ANSWERS
     */
  if (1 == *which) {
    /*
       Calculating P
       */
    cumpoi(s,xlam,p,q);
    *status = 0;
  }
  else if (2 == *which) {
    /*
       Calculating S
       */
    *s = 5.;
    T3 = inf;
    T6 = atol;
    T7 = tol;
    dstinv(&K2,&T3,&K4,&K4,&K5,&T6,&T7);
    *status = 0;
    dinvr(status,s,&fx,&qleft,&qhi);
S200:
    if (!(*status == 1)) goto S230;
    cumpoi(s,xlam,&cum,&ccum);
    if (!qporq) goto S210;
    fx = cum-*p;
    goto S220;
S210:
    fx = ccum-*q;
S220:
    dinvr(status,s,&fx,&qleft,&qhi);
    goto S200;
S230:
    if (!(*status == -1)) goto S260;
    if (!qleft) goto S240;
    *status = 1;
    *bound = 0.;
    goto S250;
S240:
    *status = 2;
    *bound = inf;
S260:
S250:
    ;
  }
  else if (3 == *which) {
    /*
       Calculating XLAM
       */
    *xlam = 5.;
    T8 = inf;
    T9 = atol;
    T10 = tol;
    dstinv(&K2,&T8,&K4,&K4,&K5,&T9,&T10);
    *status = 0;
    dinvr(status,xlam,&fx,&qleft,&qhi);
S270:
    if (!(*status == 1)) goto S300;
    cumpoi(s,xlam,&cum,&ccum);
    if (!qporq) goto S280;
    fx = cum-*p;
    goto S290;
S280:
    fx = ccum-*q;
S290:
    dinvr(status,xlam,&fx,&qleft,&qhi);
    goto S270;
S300:
    if (!(*status == -1)) goto S330;
    if (!qleft) goto S310;
    *status = 1;
    *bound = 0.;
    goto S320;
S310:
    *status = 2;
    *bound = inf;
S320:
    ;
  }
S330:
  return;
}

/**
 * Cumulative Distribution Function
 T distribution
 *
 * @param which : Integer indicating  which of the  next  parameter
 * values is to be calculated using values  of the others
 * Legal range: 1..3.
 iwhich = 1 : Calculate p and q from t and df.
 iwhich = 2 : Calculate t from p,q and df.
 iwhich = 3 : Calculate df from p,q and t.
 * @param p : The integral from -infinity to t of the t-density.
 Input range: (0,1].
 * @param q : 1-p.
 Input range: (0, 1].
 p + q = 1.0.
 * @param t : Upper limit of integration of the t-density.
 Input range: ( -infinity, +infinity).
 * @param df :  Degrees of freedom of the t-distribution.
 Input range: (0 , +infinity).
 * @param status :  (0) if calculation completed correctly.
 *              (-I) if input parameter number I is out of range.
 *               (1) if answer appears to be lower than lowest
 *                 search bound.
 *               (2) if answer appears to be higher than greatest
 *                 search bound.
 *               (3) if P + Q .ne. 1.
 *
 * @param bound : Undefined if STATUS is 0. Bound exceeded
 by parameter number I if STATUS is negative. Lower search
 bound if STATUS is 1. Upper search bound if STATUS is 2.
 * @return  any one parameter of the T
 distribution given values for the others
 */

/*
 * 
 * void pnl_cdf_t(int *which,double *p,double *q,double *t,double *df,
 * int *status,double *bound)
 * 
 * Cumulative Distribution Function
 * T distribution
 * 
 * 
 * Function
 * 
 * 
 * Calculates any one parameter of the t distribution given
 * values for the others.
 * 
 * 
 * Arguments
 * 
 * 
 * WHICH --> Integer indicating which  argument
 * values is to be calculated from the others.
 * Legal range: 1..3
 * iwhich = 1 : Calculate P and Q from T and DF
 * iwhich = 2 : Calculate T from P,Q and DF
 * iwhich = 3 : Calculate DF from P,Q and T
 * 
 * P <--> The integral from -infinity to t of the t-density.
 * Input range: (0,1].
 * 
 * Q <--> 1-P.
 * Input range: (0, 1].
 * P + Q = 1.0.
 * 
 * T <--> Upper limit of integration of the t-density.
 * Input range: ( -infinity, +infinity).
 * Search range: [ -1E300, 1E300 ]
 * 
 * DF <--> Degrees of freedom of the t-distribution.
 * Input range: (0 , +infinity).
 * Search range: [1e-300, 1E10]
 * 
 * STATUS <-- 0 if calculation completed correctly
 * -I if input parameter number I is out of range
 * 1 if answer appears to be lower than lowest
 * search bound
 * 2 if answer appears to be higher than greatest
 * search bound
 * 3 if P + Q .ne. 1
 * 
 * BOUND <-- Undefined if STATUS is 0
 * 
 * Bound exceeded by parameter number I if STATUS
 * is negative.
 * 
 * Lower search bound if STATUS is 1.
 * 
 * Upper search bound if STATUS is 2.
 * 
 * 
 * Method
 * 
 * 
 * Formula  26.5.27  of   Abramowitz   and  Stegun,   Handbook   of
 * Mathematical Functions  (1966) is used to reduce the computation
 * of the cumulative distribution function to that of an incomplete
 * beta.
 * 
 * Computation of other parameters involve a seach for a value that
 * produces  the desired  value  of P.   The search relies  on  the
 * monotinicity of P with the other parameter.
 * 
*/
void pnl_cdf_t(int *which,double *p,double *q,double *t,double *df,
               int *status,double *bound)
{
  const double tol = 1.0e-8;
  const double atol = 1.0e-50;
  const double zero = 1.0e-300;
  const double inf = 1.0e300;
  const double maxdf = 1.0e10;
  int K1 = 1;
  double K4 = 0.5e0;
  double K5 = 5.;
  double fx,cum,ccum,pq;
  unsigned long qhi,qleft,qporq;
  double T2,T3,T6,T7,T8,T9,T10,T11;
  /*
     ..
     .. Executable Statements ..
     */
  /*
     Check arguments
     */
  CHECK_WHICH (*which < 1, 1., -1);
  CHECK_WHICH (*which > 3, 3., -1);

  if (*which == 1) goto S70;
  /* P */
  if (!(*p <= 0. || *p > 1.)) goto S60;
  if (!(*p <= 0.)) goto S40;
  *bound = 0.;
  goto S50;
S40:
  *bound = 1.;
S50:
  *status = -2;
  return;
S70:
S60:
  if (*which == 1) goto S110;
  /* Q */
  if (!(*q <= 0. || *q > 1.)) goto S100;
  if (!(*q <= 0.)) goto S80;
  *bound = 0.;
  goto S90;
S80:
  *bound = 1.;
S90:
  *status = -3;
  return;
S110:
S100:
  if (*which == 3) goto S130;
  /*
     DF
     */
  if (!(*df <= 0.)) goto S120;
  *bound = 0.;
  *status = -5;
  return;
S130:
S120:
  if (*which == 1) goto S170;
  /*
     P + Q
     */
  pq = *p+*q;
  if (!(fabs(pq-0.5e0-0.5e0) > 3.*spmpar(&K1))) goto S160;
  if (!(pq < 0.)) goto S140;
  *bound = 0.;
  goto S150;
S140:
  *bound = 1.;
S150:
  *status = 3;
  return;
S170:
S160:
  if (!(*which == 1)) qporq = *p <= *q;
  /*
     Select the minimum of P or Q
     Calculate ANSWERS
     */
  if (1 == *which) {
    /*
       Computing P and Q
       */
    cumt(t,df,p,q);
    *status = 0;
  }
  else if (2 == *which) {
    /*
       Computing T
       .. Get initial approximation for T
       */
    *t = dt1(p,q,df);
    T2 = -inf;
    T3 = inf;
    T6 = atol;
    T7 = tol;
    dstinv(&T2,&T3,&K4,&K4,&K5,&T6,&T7);
    *status = 0;
    dinvr(status,t,&fx,&qleft,&qhi);
S180:
    if (!(*status == 1)) goto S210;
    cumt(t,df,&cum,&ccum);
    if (!qporq) goto S190;
    fx = cum-*p;
    goto S200;
S190:
    fx = ccum-*q;
S200:
    dinvr(status,t,&fx,&qleft,&qhi);
    goto S180;
S210:
    if (!(*status == -1)) goto S240;
    if (!qleft) goto S220;
    *status = 1;
    *bound = -inf;
    goto S230;
S220:
    *status = 2;
    *bound = inf;
S240:
S230:
    ;
  }
  else if (3 == *which) {
    /*
       Computing DF
       */
    *df = 5.;
    T8 = zero;
    T9 = maxdf;
    T10 = atol;
    T11 = tol;
    dstinv(&T8,&T9,&K4,&K4,&K5,&T10,&T11);
    *status = 0;
    dinvr(status,df,&fx,&qleft,&qhi);
S250:
    if (!(*status == 1)) goto S280;
    cumt(t,df,&cum,&ccum);
    if (!qporq) goto S260;
    fx = cum-*p;
    goto S270;
S260:
    fx = ccum-*q;
S270:
    dinvr(status,df,&fx,&qleft,&qhi);
    goto S250;
S280:
    if (!(*status == -1)) goto S310;
    if (!qleft) goto S290;
    *status = 1;
    *bound = zero;
    goto S300;
S290:
    *status = 2;
    *bound = maxdf;
S300:
    ;
  }
S310:
  return;
}

/*
 * 
 * 
 * void cumbet(double *x,double *y,double *a,double *b,double *cum,
 * double *ccum)
 * 
 * Double precision cUMulative incomplete BETa distribution
 * 
 * 
 * Function
 * 
 * 
 * Calculates the pnl_cdf_ to X of the incomplete beta distribution
 * with parameters a and b.  This is the integral from 0 to x
 * of (1/B(a,b))*f(t)) where f(t) = t**(a-1) * (1-t)**(b-1)
 * 
 * 
 * Arguments
 * 
 * 
 * X --> Upper limit of integration.
 * X is DOUBLE PRECISION
 * 
 * Y --> 1 - X.
 * Y is DOUBLE PRECISION
 * 
 * A --> First parameter of the beta distribution.
 * A is DOUBLE PRECISION
 * 
 * B --> Second parameter of the beta distribution.
 * B is DOUBLE PRECISION
 * 
 * CUM <-- Cumulative incomplete beta distribution.
 * CUM is DOUBLE PRECISION
 * 
 * CCUM <-- Compliment of Cumulative incomplete beta distribution.
 * CCUM is DOUBLE PRECISION
 * 
 * 
 * Method
 * 
 * 
 * Calls the routine BRATIO.
 * 
 * References
 * 
 * Didonato, Armido R. and Morris, Alfred H. Jr. (1992) Algorithim
 * 708 Significant Digit Computation of the Incomplete Beta Function
 * Ratios. ACM ToMS, Vol.18, No. 3, Sept. 1992, 360-373.
 * 
 */
static void cumbet(double *x,double *y,double *a,double *b,double *cum,
                   double *ccum)
{
  int ierr;
  if ((*x <= 0.)) 
    {
      *cum = 0.;
      *ccum = 1.;
    }
  else if (*y <= 0.) 
    {
      *cum = 1.;
      *ccum = 0.;
    }
  else
    {
      bratio(a,b,x,y,cum,ccum,&ierr);
    }
}

/*
 * 
 * void cumbin(double *s,double *xn,double *pr,double *ompr,
 * double *cum,double *ccum)
 * 
 * CUmulative BINomial distribution
 * 
 * 
 * Function
 * 
 * 
 * Returns the probability   of 0  to  S  successes in  XN   binomial
 * trials, each of which has a probability of success, PBIN.
 * 
 * 
 * Arguments
 * 
 * 
 * S --> The upper limit of cumulation of the binomial distribution.
 * S is DOUBLE PRECISION
 * 
 * XN --> The number of binomial trials.
 * XN is DOUBLE PRECISIO
 * 
 * PBIN --> The probability of success in each binomial trial.
 * PBIN is DOUBLE PRECIS
 * 
 * OMPR --> 1 - PBIN
 * OMPR is DOUBLE PRECIS
 * 
 * CUM <-- Cumulative binomial distribution.
 * CUM is DOUBLE PRECISI
 * 
 * CCUM <-- Compliment of Cumulative binomial distribution.
 * CCUM is DOUBLE PRECIS
 * 
 * 
 * Method
 * 
 * 
 * Formula  26.5.24    of   Abramowitz  and    Stegun,  Handbook   of
 * Mathematical   Functions (1966) is   used  to reduce the  binomial
 * distribution  to  the  cumulative    beta distribution.
 * 
 */
static void cumbin(double *s,double *xn,double *pr,double *ompr,
                   double *cum,double *ccum)
{
  if (*s < *xn) 
    {
      double T1 = *s+1.;
      double T2 = *xn-*s;
      cumbet(pr,ompr,&T1,&T2,ccum,cum);
    }
  else
    {
      *cum = 1.;
      *ccum = 0.;
    }
}

/*
 * 
 * void cumchi(double *x,double *df,double *cum,double *ccum)
 * CUMulative of the CHi-square distribution
 * 
 * 
 * Function
 * 
 * 
 * Calculates the cumulative chi-square distribution.
 * 
 * 
 * Arguments
 * 
 * 
 * X       --> Upper limit of integration of the
 * chi-square distribution.
 * X is DOUBLE PRECISION
 * 
 * DF      --> Degrees of freedom of the
 * chi-square distribution.
 * DF is DOUBLE PRECISION
 * 
 * CUM <-- Cumulative chi-square distribution.
 * CUM is DOUBLE PRECISIO
 * 
 * CCUM <-- Compliment of Cumulative chi-square distribution.
 * CCUM is DOUBLE PRECISI
 * 
 * 
 * Method
 * 
 * 
 * Calls incomplete gamma function (CUMGAM)
 * 
 */
static void cumchi(double *x,double *df,double *cum,double *ccum)
{
  double a,xx;
  a = *df*0.5e0;
  xx = *x*0.5e0;
  cumgam(&xx,&a,cum,ccum);
  return;
}

/*
 * 
 * void cumchn(double *x,double *df,double *pnonc,double *cum,
 * double *ccum)
 * 
 * CUMulative of the Non-central CHi-square distribution
 * 
 * 
 * Function
 * 
 * 
 * Calculates     the       cumulative      non-central    chi-square
 * distribution, i.e.,  the probability   that  a   random   variable
 * which    follows  the  non-central chi-square  distribution,  with
 * non-centrality  parameter    PNONC  and   continuous  degrees   of
 * freedom DF, is less than or equal to X.
 * 
 * 
 * Arguments
 * 
 * 
 * X       --> Upper limit of integration of the non-central
 * chi-square distribution.
 * X is DOUBLE PRECISION
 * 
 * DF      --> Degrees of freedom of the non-central
 * chi-square distribution.
 * DF is DOUBLE PRECISION
 * 
 * PNONC   --> Non-centrality parameter of the non-central
 * chi-square distribution.
 * PNONC is DOUBLE PRECIS
 * 
 * CUM <-- Cumulative non-central chi-square distribution.
 * CUM is DOUBLE PRECISIO
 * 
 * CCUM <-- Compliment of Cumulative non-central chi-square distribut
 * CCUM is DOUBLE PRECISI
 * 
 * 
 * Method
 * 
 * 
 * Uses  formula  26.4.25   of  Abramowitz  and  Stegun, Handbook  of
 * Mathematical    Functions,  US   NBS   (1966)    to calculate  the
 * non-central chi-square.
 * 
 * 
 * Variables
 * 
 * 
 * EPS     --- Convergence criterion.  The sum stops when a
 * term is less than EPS*SUM.
 * EPS is DOUBLE PRECISIO
 * 
 * NTIRED  --- Maximum number of terms to be evaluated
 * in each sum.
 * NTIRED is INTEGER
 * 
 * QCONV   --- .TRUE. if convergence achieved -
 * i.e., program did not stop on NTIRED criterion.
 * QCONV is LOGICAL
 * 
 * CCUM <-- Compliment of Cumulative non-central
 * chi-square distribution.
 * CCUM is DOUBLE PRECISI
 * 
 */
static void cumchn(double *x,double *df,double *pnonc,double *cum,
                   double *ccum)
{
  double eps = 1.0e-5;
  int ntired = 1000;
  static double adj,centaj,centwt,chid2,dfd2,lcntaj,lcntwt,lfact,pcent,pterm,sum,
                sumadj,term,wt,xnonc;
  int i,icent,iterb,iterf;
  double T1,T2,T3;
  /*
     ..
     .. Executable Statements ..
     */
  if (*x <= 0.) 
    {
      *cum = 0.;
      *ccum = 1.;
      return;
    }

  if (*pnonc <= 1.0e-10)
    {
      /*
       * When non-centrality parameter is (essentially) zero,
       * use cumulative chi-square distribution
       */
      cumchi(x,df,cum,ccum);
      return;
    }

  xnonc = *pnonc/2.;
  /*
   **********************************************************************
   The following code calcualtes the weight, chi-square, and
   adjustment term for the central term in the infinite series.
   The central term is the one in which the poisson weight is
   greatest.  The adjustment term is the amount that must
   be subtracted from the chi-square to move up two degrees
   of freedom.
   **********************************************************************
   */
  icent = fifidint(xnonc);
  if (icent == 0) icent = 1;
  chid2 = *x/2.;
  /*
   * Calculate central weight term
   */
  T1 = (double)(icent+1);
  lfact = alngam(&T1);
  lcntwt = -xnonc+(double)icent*log(xnonc)-lfact;
  centwt = exp(lcntwt);
  /*
   * Calculate central chi-square
   */
  T2 = *df + 2. * icent;
  cumchi(x,&T2,&pcent,ccum);
  /*
   * Calculate central adjustment term
   */
  dfd2 = (*df + 2. * icent) / 2.;
  T3 = 1.+dfd2;
  lfact = alngam(&T3);
  lcntaj = dfd2*log(chid2)-chid2-lfact;
  centaj = exp(lcntaj);
  sum = centwt*pcent;
  /*
   **********************************************************************
   Sum backwards from the central term towards zero.
   Quit whenever either
   (1) the zero term is reached, or
   (2) the term gets small relative to the sum, or
   (3) More than NTIRED terms are totaled.
   **********************************************************************
   */
  iterb = 0;
  sumadj = 0.;
  adj = centaj;
  wt = centwt;
  i = icent;
  while ( 1 )
    {
      dfd2 = (*df + 2. * i)/2.;
      /*
         Adjust chi-square for two fewer degrees of freedom.
         The adjusted value ends up in PTERM.
         */
      adj = adj*dfd2/chid2;
      sumadj += adj;
      pterm = pcent+sumadj;
      /*
         Adjust poisson weight for J decreased by one
         */
      wt *= ((double)i/xnonc);
      term = wt*pterm;
      sum += term;
      i -= 1;
      iterb += 1;
      if ((iterb > ntired) || (sum < 1.0e-20 || term < eps*sum) || i == 0) break;
    }

  iterf = 0;
  /*
   **********************************************************************
   Now sum forward from the central term towards infinity.
   Quit when either
   (1) the term gets small relative to the sum, or
   (2) More than NTIRED terms are totaled.
   **********************************************************************
   */
  sumadj = adj = centaj;
  wt = centwt;
  i = icent;

  while ( 1 )
    {
      /*
       * Update weights for next higher J
       */
      wt *= (xnonc/(double)(i+1));
      /*
       * Calculate PTERM and add term to sum
       */
      pterm = pcent-sumadj;
      term = wt*pterm;
      sum += term;
      /*
       * Update adjustment term for DF for next iteration
       */
      i += 1;
      dfd2 = (*df + 2. * i)/2.;
      adj = adj*chid2/dfd2;
      sumadj += adj;
      iterf += 1;

      if ((iterf > ntired) || (sum < 1.0e-20 || term < eps*sum))  break;
    }
  *cum = sum;
  *ccum = 0.5e0+(0.5e0-*cum);
}



/* CUMulative F distribution
 * 
 * Computes  the  integral from  0  to  F of  the f-density  with DFN
 * and DFD degrees of freedom.
 * 
 * 
 * Arguments
 * 
 * 
 * F --> Upper limit of integration of the f-density.
 * DFN --> Degrees of freedom of the numerator sum of squares.
 * DFD --> Degrees of freedom of the denominator sum of squares.
 * CUM <-- Cumulative f distribution.
 * CCUM <-- Compliment of Cumulative f distribution.
 * 
 * Formula  26.5.28 of  Abramowitz and   Stegun   is  used to  reduce
 * the cumulative F to a cumulative beta distribution.
 * 
 * Note if F is less than or equal to 0, 0 is returned.
 */
static void cumf(double *f,double *dfn,double *dfd,double *cum,double *ccum)
{
  double dsum, prod, xx, yy, T1, T2;
  int ierr;

  if (*f <= 0.)
    {
      *cum = 0.;
      *ccum = 1.;
      return;
    }

  prod = *dfn * *f;
  /*
     XX is such that the incomplete beta with parameters
     DFD/2 and DFN/2 evaluated at XX is 1 - CUM or CCUM
     YY is 1 - XX
     Calculate the smaller of XX and YY accurately
     */
  dsum = *dfd + prod;
  xx = *dfd / dsum;
  if (xx > 0.5)
    {
      yy = prod / dsum;
      xx = 1. - yy;
    }
  else
    yy = 1. -xx;
  T1 = *dfd * 0.5;
  T2 = *dfn * 0.5;
  bratio (&T1, &T2, &xx, &yy, ccum, cum, &ierr);
  return;
}

/*
 * 
 * F -NON- -C-ENTRAL F DISTRIBUTION
 * 
 * 
 * 
 * Function
 * 
 * 
 * COMPUTES NONCENTRAL F DISTRIBUTION WITH DFN AND DFD
 * DEGREES OF FREEDOM AND NONCENTRALITY PARAMETER PNONC
 * 
 * 
 * Arguments
 * 
 * 
 * X --> UPPER LIMIT OF INTEGRATION OF NONCENTRAL F IN EQUATION
 * 
 * DFN --> DEGREES OF FREEDOM OF NUMERATOR
 * 
 * DFD -->  DEGREES OF FREEDOM OF DENOMINATOR
 * 
 * PNONC --> NONCENTRALITY PARAMETER.
 * 
 * CUM <-- CUMULATIVE NONCENTRAL F DISTRIBUTION
 * 
 * CCUM <-- COMPLIMENT OF CUMMULATIVE
 * 
 * 
 * Method
 * 
 * 
 * USES FORMULA 26.6.20 OF REFERENCE FOR INFINITE SERIES.
 * SERIES IS CALCULATED BACKWARD AND FORWARD FROM J = LAMBDA/2
 * (THIS IS THE TERM WITH THE LARGEST POISSON WEIGHT) UNTIL
 * THE CONVERGENCE CRITERION IS MET.
 * 
 * FOR SPEED, THE INCOMPLETE BETA FUNCTIONS ARE EVALUATED
 * BY FORMULA 26.5.16.
 * 
 * 
 * REFERENCE
 * 
 * 
 * HANDBOOD OF MATHEMATICAL FUNCTIONS
 * EDITED BY MILTON ABRAMOWITZ AND IRENE A. STEGUN
 * NATIONAL BUREAU OF STANDARDS APPLIED MATEMATICS SERIES - 55
 * MARCH 1965
 * P 947, EQUATIONS 26.6.17, 26.6.18
 * 
 * 
 * Note
 * 
 * 
 * THE SUM CONTINUES UNTIL A SUCCEEDING TERM IS LESS THAN EPS
 * TIMES THE SUM (OR THE SUM IS LESS THAN 1.0E-20).  EPS IS
 * SET TO 1.0E-4 IN A DATA STATEMENT WHICH CAN BE CHANGED.
 * 
 */
static void cumfnc(double *f,double *dfn,double *dfd,double *pnonc,
                   double *cum,double *ccum)
{
  double eps = 1.0e-4;
  static double dsum,dummy,prod,xx,yy,adn,aup,b,betdn,betup,centwt,dnterm,sum,
                upterm,xmult,xnonc;
  int i,icent,ierr;
  double T1,T2,T3,T4,T5,T6;

  if ((*f <= 0.)) 
    {
      *cum = 0.;
      *ccum = 1.;
      return;
    }
  if ((*pnonc < 1.0e-10)) 
    {
      /*
       * Handle case in which the non-centrality parameter is
       * (essentially) zero.
       */
      cumf(f,dfn,dfd,cum,ccum);
      return;
    }
  xnonc = *pnonc/2.;
  /*
   * Calculate the central term of the poisson weighting factor.
   */
  icent = xnonc;
  if (icent == 0) icent = 1;
  /*
   * Compute central weight term
   */
  T1 = (double)(icent+1);
  centwt = exp(-xnonc+(double)icent*log(xnonc)-alngam(&T1));
  /*
   * Compute central incomplete beta term
   * Assure that minimum of arg to beta and 1 - arg is computed
   * accurately.
   */
  prod = *dfn**f;
  dsum = *dfd+prod;
  yy = *dfd/dsum;
  if (yy > 0.5) 
    {
      xx = prod/dsum;
      yy = 1.-xx;
    }
  else  xx = 1.-yy;
  T2 = *dfn*0.5+(double)icent;
  T3 = *dfd*0.5;
  bratio(&T2,&T3,&xx,&yy,&betdn,&dummy,&ierr);
  adn = *dfn/2.+(double)icent;
  aup = adn;
  b = *dfd/2.;
  betup = betdn;
  sum = centwt*betdn;
  /*
     Now sum terms backward from icent until convergence or all done
     */
  xmult = centwt;
  i = icent;
  T4 = adn+b;
  T5 = adn+1.;
  dnterm = exp(alngam(&T4)-alngam(&T5)-alngam(&b)+adn*log(xx)+b*log(yy));

  while ( !((sum < 1.0e-20) || ((xmult*betdn) < eps*sum) || i <= 0) )
    {
      xmult *= ((double)i/xnonc);
      i -= 1;
      adn -= 1.0;
      dnterm = (adn+1.0)/((adn+b)*xx)*dnterm;
      betdn += dnterm;
      sum += (xmult*betdn);
    }

  i = icent+1;
  /*
     Now sum forwards until convergence
     */
  xmult = centwt;
  if (aup-1.0+b == 0)
    {
      upterm = exp(-alngam(&aup)-alngam(&b)+(aup-1.0)*log(xx)+
                   b*log(yy));
    }
  else 
    {
      T6 = aup-1.0+b;
      upterm = exp(alngam(&T6)-alngam(&aup)-alngam(&b)+(aup-1.0)*log(xx)+b*
                   log(yy));
    }

  while ( 1 )
    {
      xmult *= (xnonc/(double)i);
      i += 1;
      aup += 1.0;
      upterm = (aup+b-2.)*xx/(aup-1.0)*upterm;
      betup -= upterm;
      sum += (xmult*betup);
      if ( (sum < 1.0e-20) || ((xmult*betup) < eps*sum) ) break;
    }
  *cum = sum;
  *ccum = 0.5e0+(0.5e0-*cum);
}


/* Double precision cUMulative incomplete GAMma distribution
 *  
 * Computes   the  cumulative        of    the     incomplete   gamma
 * distribution, i.e., the integral from 0 to X of
 * (1/GAM(A))*EXP(-T)*T**(A-1) DT
 * where GAM(A) is the complete gamma function of A, i.e.,
 * GAM(A) = integral from 0 to infinity of
 * EXP(-T)*T**(A-1) DT
 *  
 * X --> The upper limit of integration of the incomplete gamma.
 * A --> The shape parameter of the incomplete gamma.
 * CUM <-- Cumulative incomplete gamma distribution.
 * CCUM <-- Compliment of Cumulative incomplete gamma distribution.
 */

/**
 * Gamma(a, 1) cdf
 *
 * @param x the upper bound of the integral
 * @param a the parameter of the distribution
 * @param cum on output = int_0^x (1/gam(a))*exp(-t)*t**(a-1) dt
 * @param ccum on output = 1-cum
 */
static void cumgam(double *x,double *a,double *cum,double *ccum)
{
  int K1 = 0;

  if (*x <= 0)
    {
      *cum = 0.;
      *ccum = 1.;
    }
  else
    gratio(a,x,cum,ccum,&K1);
  return;
}

/*
 * 
 * void cumnbn(double *s,double *xn,double *pr,double *ompr,
 * double *result,double *ccum)
 * 
 * CUmulative Negative BINomial distribution
 * 
 * 
 * Function
 * 
 * 
 * Returns the probability that it there will be S or fewer failures
 * before there are XN successes, with each binomial trial having
 * a probability of success PR.
 * 
 * Prob(# failures = S | XN successes, PR)  =
 * ( XN + S - 1 )
 * (            ) * PR^XN * (1-PR)^S
 * (      S     )
 * 
 * 
 * Arguments
 * 
 * 
 * S --> The number of failures
 * S is DOUBLE PRECISION
 * 
 * XN --> The number of successes
 * XN is DOUBLE PRECISIO
 * 
 * PR --> The probability of success in each binomial trial.
 * PR is DOUBLE PRECISIO
 * 
 * OMPR --> 1 - PR
 * OMPR is DOUBLE PRECIS
 * 
 * RESULT <-- Cumulative negative binomial distribution.
 * RESULT is DOUBLE PRECISI
 * 
 * CCUM <-- Compliment of Cumulative negative binomial distribution.
 * CCUM is DOUBLE PRECIS
 * 
 * 
 * Method
 * 
 * 
 * Formula  26.5.26    of   Abramowitz  and    Stegun,  Handbook   of
 * Mathematical   Functions (1966) is   used  to reduce the  negative
 * binomial distribution to the cumulative beta distribution.
 * 
 */
static void cumnbn(double *s,double *xn,double *pr,double *ompr,
                   double *result,double *ccum)
{
  double T1;
  T1 = *s + 1.;
  cumbet(pr,ompr,xn,&T1,result,ccum);
  return;
}

/*
 * 
 * void cumnor(double *arg,double *result,double *ccum)
 * 
 * 
 * Function
 * 
 * 
 * Computes the cumulative  of    the  normal   distribution,   i.e.,
 * the integral from -infinity to x of
 * (1/sqrt(2*pi)) exp(-u*u/2) du
 * 
 * X --> Upper limit of integration.
 * X is DOUBLE PRECISION
 * 
 * RESULT <-- Cumulative normal distribution.
 * RESULT is DOUBLE PRECISION
 * 
 * CCUM <-- Compliment of Cumulative normal distribution.
 * CCUM is DOUBLE PRECISION
 * 
 * Renaming of function ANORM from:
 * 
 * Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
 * Package of Special Function Routines and Test Drivers"
 * acm Transactions on Mathematical Software. 19, 22-32.
 * 
 * with slight modifications to return ccum and to deal with
 * machine constants.
 * 
 * **********************************************************************
 * Original Comments:
 * ------------------------------------------------------------------
 * 
 * This function evaluates the normal distribution function:
 * 
 * / x
 * 1       |       -t*t/2
 * P(x) = ----------- |      e       dt
 * sqrt(2 pi)  |
 * /-oo
 * 
 * The main computation evaluates near-minimax approximations
 * derived from those in "Rational Chebyshev approximations for
 * the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
 * This transportable program uses rational functions that
 * theoretically approximate the normal distribution function to
 * at least 18 significant decimal digits.  The accuracy achieved
 * depends on the arithmetic system, the compiler, the intrinsic
 * functions, and proper selection of the machine-dependent
 * constants.
 * 
 * *******************************************************************
 * *******************************************************************
 * 
 * Explanation of machine-dependent constants.
 * 
 * MIN   = smallest machine representable number.
 * 
 * EPS   = argument below which anorm(x) may be represented by
 * 0.5  and above which  x*x  will not underflow.
 * A conservative value is the largest machine number X
 * such that   1.0 + X = 1.0   to machine precision.
 * *******************************************************************
 * *******************************************************************
 * 
 * Error returns
 * 
 * The program returns  ANORM = 0     for  ARG .LE. XLOW.
 * 
 * 
* Intrinsic functions required are:
* 
* ABS, AINT, EXP
* 
* 
* Author: W. J. Cody
* Mathematics and Computer Science Division
* Argonne National Laboratory
* Argonne, IL 60439
* 
* Latest modification: March 15, 1992
* 
*/
static void cumnor(double *arg,double *result,double *ccum)
{
  double a[5] = {
    2.2352520354606839287e00,1.6102823106855587881e02,1.0676894854603709582e03,
    1.8154981253343561249e04,6.5682337918207449113e-2
  };
  double b[4] = {
    4.7202581904688241870e01,9.7609855173777669322e02,1.0260932208618978205e04,
    4.5507789335026729956e04
  };
  double c[9] = {
    3.9894151208813466764e-1,8.8831497943883759412e00,9.3506656132177855979e01,
    5.9727027639480026226e02,2.4945375852903726711e03,6.8481904505362823326e03,
    1.1602651437647350124e04,9.8427148383839780218e03,1.0765576773720192317e-8
  };
  double d[8] = {
    2.2266688044328115691e01,2.3538790178262499861e02,1.5193775994075548050e03,
    6.4855582982667607550e03,1.8615571640885098091e04,3.4900952721145977266e04,
    3.8912003286093271411e04,1.9685429676859990727e04
  };
  double half = 0.5;
  double p[6] = {
    2.1589853405795699e-1,1.274011611602473639e-1,2.2235277870649807e-2,
    1.421619193227893466e-3,2.9112874951168792e-5,2.307344176494017303e-2
  };
  double one = 1.;
  double q[5] = {
    1.28426009614491121e00,4.68238212480865118e-1,6.59881378689285515e-2,
    3.78239633202758244e-3,7.29751555083966205e-5
  };
  double sixten = 1.60e0;
  double sqrpi = 3.9894228040143267794e-1;
  double thrsh = 0.66291e0;
  double root32 = 5.656854248e0;
  double zero = 0.;
  int K1 = 1;
  int K2 = 2;
  int i;
  double del,eps,temp,x,xden,xnum,y,xsq,min;
  /*
     Machine dependent constants
     */
  eps = spmpar(&K1)*0.5e0;
  min = spmpar(&K2);
  x = *arg;
  y = fabs(x);
  if (y <= thrsh) 
    {
      /*
         Evaluate  anorm  for  |X| <= 0.66291
         */
      xsq = zero;
      if (y > eps) xsq = x*x;
      xnum = a[4]*xsq;
      xden = xsq;
      for(i=0; i<3; i++) 
        {
          xnum = (xnum+a[i])*xsq;
          xden = (xden+b[i])*xsq;
        }
      *result = x*(xnum+a[3])/(xden+b[3]);
      temp = *result;
      *result = half+temp;
      *ccum = half-temp;
    }
  /*
     Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
     */
  else if (y <= root32) 
    {
      xnum = c[8]*y;
      xden = y;
      for(i=0; i<7; i++) 
        {
          xnum = (xnum+c[i])*y;
          xden = (xden+d[i])*y;
        }
      *result = (xnum+c[7])/(xden+d[7]);
      xsq = fifdint(y*sixten)/sixten;
      del = (y-xsq)*(y+xsq);
      *result = exp(-(xsq*xsq*half))*exp(-(del*half))**result;
      *ccum = one-*result;
      if (x > zero) 
        {
          temp = *result;
          *result = *ccum;
          *ccum = temp;
        }
    }
  /*
     Evaluate  anorm  for |X| > sqrt(32)
     */
  else  
    {
      *result = zero;
      xsq = one/(x*x);
      xnum = p[5]*xsq;
      xden = xsq;
      for(i=0; i<4; i++) 
        {
          xnum = (xnum+p[i])*xsq;
          xden = (xden+q[i])*xsq;
        }
      *result = xsq*(xnum+p[4])/(xden+q[4]);
      *result = (sqrpi-*result)/y;
      xsq = fifdint(x*sixten)/sixten;
      del = (x-xsq)*(x+xsq);
      *result = exp(-(xsq*xsq*half))*exp(-(del*half))**result;
      *ccum = one-*result;
      if (x > zero) 
        {
          temp = *result;
          *result = *ccum;
          *ccum = temp;
        }
    }
  if (*result < min) *result = 0.;
  /*
     Fix up for negative argument, erf, etc.
     ----------Last card of ANORM ----------
     */
  if (*ccum < min) *ccum = 0.;
}

/* CUMulative POIsson distribution
 *  
 * Returns the  probability  of  S   or  fewer events in  a   Poisson
 * distribution with mean XLAM.
 *  
 * S --> Upper limit of cumulation of the Poisson.
 * XLAM --> Mean of the Poisson distribution.
 * CUM <-- Cumulative poisson distribution.
 * CCUM <-- Compliment of Cumulative poisson distribution.
 *  
 * Uses formula  26.4.21   of   Abramowitz and  Stegun,  Handbook  of
 * Mathematical   Functions  to reduce   the   cumulative Poisson  to
 * the cumulative chi-square distribution.
 */
static void cumpoi(double *s,double *xlam,double *cum,double *ccum)
{
  double chi,df;
  df = 2. * (*s + 1.);
  chi = 2. * *xlam;
  cumchi(&chi,&df,ccum,cum);
  return;
}

/* CUMulative T-distribution
 *  
 * Computes the integral from -infinity to T of the t-density.
 *  
 * T --> Upper limit of integration of the t-density.
 * DF --> Degrees of freedom of the t-distribution.
 * CUM <-- Cumulative t-distribution.
 * CCUM <-- Compliment of Cumulative t-distribution.
 *  
 * Formula 26.5.27   of     Abramowitz  and   Stegun,    Handbook  of
 * Mathematical Functions  is   used   to  reduce the  t-distribution
 * to an incomplete beta.
 */
static void cumt(double *t,double *df,double *cum,double *ccum)
{
  double K2 = 0.5e0;
  double xx,a,oma,tt,yy,dfptt,T1;


  tt = *t * *t;
  dfptt = *df + tt;
  xx = *df / dfptt;
  yy = tt / dfptt;
  T1 = 0.5 * *df;
  cumbet(&xx, &yy, &T1, &K2, &a, &oma);

  if (*t <= 0.)
    {
      *cum = 0.5 * a;
      *ccum = oma + *cum;
    }
  else
    {
      *ccum = 0.5 * a;
      *cum = oma + *ccum;
    }
}

/* Double precision EVALuate a PoLynomial at X
 *  
 * returns
 * A(1) + A(2)*X + ... + A(N)*X**(N-1)
 *  
 * A --> Array of coefficients of the polynomial.
 * N --> Length of A, also degree of polynomial - 1.
 * X --> Point at which the polynomial is to be evaluated.
 */
static double devlpl(double a[],int *n,double *x)
{
  double term;
  int i;

  term = a[*n-1];
  for(i= *n-1-1; i>=0; i--) term = a[i]+term**x;
  return term;
}

/*
 * 
 * double dinvnr(double *p,double *q)
 * Double precision NoRmal distribution INVerse
 * 
 * 
 * Function
 * 
 * 
 * Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
 * infinity to X of (1/SQRT(2*M_PI)) EXP(-U*U/2) dU is P
 * 
 * 
 * Arguments
 * 
 * 
 * P --> The probability whose normal deviate is sought.
 * P is DOUBLE PRECISION
 * 
 * Q --> 1-P
 * P is DOUBLE PRECISION
 * 
 * 
 * Method
 * 
 * 
 * The  rational   function   on  page 95    of Kennedy  and  Gentle,
 * Statistical Computing, Marcel Dekker, NY , 1980 is used as a start
 * value for the Newton method of finding roots.
 * 
 * 
 * Note
 * 
 * 
 * If P or Q .lt. machine EPS returns +/- DINVNR(EPS)
 * 
 */
static double dinvnr(double *p,double *q)
{
  const int maxit = 100;
  const double eps = 1.0e-13;
  const double r2pi = 0.3989422804014326e0;
  double dinvnr_0,strtx,xcur,cum,ccum,pp,dx;
  int i;
  unsigned long qporq;
  /* FIND MINIMUM OF P AND Q */
  qporq = *p <= *q;
  if (!qporq) pp = *q; else pp = *p;

  /* INITIALIZATION STEP */
  strtx = stvaln (&pp);
  xcur = strtx;
  /* NEWTON INTERATIONS */
  for (i=1; i<=maxit; i++)
    {
      cumnor (&xcur, &cum, &ccum);
      dx = (cum-pp) / (r2pi*exp(-0.5 * xcur * xcur));
      xcur -= dx;
      if (fabs(dx/xcur) < eps) 
        {
          /* IF WE GET HERE, NEWTON HAS SUCCEDED */
          dinvnr_0 = xcur;
          if (!qporq) dinvnr_0 = -dinvnr_0;
          return dinvnr_0;
        }
    }
  dinvnr_0 = strtx;
  /* IF WE GET HERE, NEWTON HAS FAILED */
  if (!qporq) dinvnr_0 = -dinvnr_0;
  return dinvnr_0;
}

static void E0000(int IENTRY,int *status,double *x,double *fx,
                  unsigned long *qleft,unsigned long *qhi,double *zabsst,
                  double *zabsto,double *zbig,double *zrelst,
                  double *zrelto,double *zsmall,double *zstpmu)
{
  static double absstp,abstol,big,fbig,fsmall,relstp,reltol,small,step,stpmul,xhi,
                xlb,xlo,xsave,xub,yy;
  static int i99999;
  unsigned long qbdd,qdum1,qdum2;
  static unsigned long qcond, qincr;
  switch(IENTRY)
    {
    case 0: goto DINVR; 
    case 1: goto DSTINV;
    }
DINVR:
  if ( *status == 0 )
    {
      if (!( (small <= *x) && (*x <= big) )) ftnstop(" SMALL, X, BIG not monotone in INVR");
      /* See that SMALL and BIG bound the zero and set QINCR */
      xsave = *x; *x = small;
      i99999 = 1; /* GET-FUNCTION-VALUE */
      *status = 1; /* TO GET-FUNCTION-VALUE */
      return;
    }

  switch((int)i99999)
    {
    case 1: 
      fsmall = *fx;
      *x = big;
      i99999 = 2; /* GET-FUNCTION-VALUE */
      *status = 1; /* TO GET-FUNCTION-VALUE */
      return;

    case 2: 
      fbig = *fx;
      qincr = fbig > fsmall;
      if (qincr) 
        {
          if (fsmall > 0.) 
            {
              *status = -1; *qleft = *qhi = 1;
              return;
            }
          if (fbig < 0.) 
            {
              *status = -1; *qleft = *qhi = 0;
              return;
            }
        }
      else
        {
          if (fsmall < 0.) 
            {
              *status = -1; *qleft = 1; *qhi = 0;
              return;
            }
          if (fbig > 0.) 
            {
              *status = -1; *qleft = 0; *qhi = 1;
              return;
            }
        }
      *x = xsave;
      step = MAX(absstp,relstp*fabs(*x));
      i99999 = 3; /* * YY = F(X) - Y ; GET-FUNCTION-VALUE */
      *status = 1; /* TO GET-FUNCTION-VALUE */
      return;

    case 3: 
      yy = *fx;
      if (yy == 0.) 
        {
          *status = 0;
          return;
        }
      /* HANDLE CASE IN WHICH WE MUST STEP HIGHER */
      if (! ((qincr && yy < 0.) || (!qincr && yy > 0.)) ) 
        {
          /* HANDLE CASE IN WHICH WE MUST STEP LOWER */
          xub = xsave;
          xlb = MAX(xub-step,small);
          *x = xlb;
          i99999 = 5; /* GET-FUNCTION-VALUE */
        }
      else
        {
          xlb = xsave;
          xub = MIN(xlb+step,big);
          *x = xub; /* YY = F(XUB) - Y */
          i99999 = 4; /* GET-FUNCTION-VALUE */
        }
      *status = 1; /* TO GET-FUNCTION-VALUE */
      return;

    case 4: 
      break;

    case 5: 
      yy = *fx;
      qbdd = (qincr && yy <= 0.) || (!qincr && yy >= 0.);
      qcond = qbdd || (xlb <= small);
      if (!qcond) 
        {
          step = stpmul*step;
          xub = xlb;
          xlb = MAX(xub-step,small);
          /* YY = F(XLB) - Y */
          *x = xlb;
          /* GET-FUNCTION-VALUE */
          *status = 1; /* TO GET-FUNCTION-VALUE */
          return;
        }

      if (((xlb <= small) && !qbdd))
        {
          *status = -1;
          *qleft = 1;
          *qhi = qincr;
          *x = small;
          return;
        }

    case 6: 
      if (!(*status == 1)) 
        {
          *x = xlo;
          *status = 0;
          return;
        }

    default: 
      if ( !qcond )
        {
          *x = xub; 
          i99999 = 4; /* GET-FUNCTION-VALUE */
          *status = 1; /* TO GET-FUNCTION-VALUE */
          return;
        }
      break;
    }


  if ( (i99999 != 5) && (i99999 != 6) )
    {
      /* Extension of cases 4 and default from the above switch */
      yy = *fx;
      qbdd =( qincr && yy >= 0.) || ( !qincr && yy <= 0.);
      qcond = qbdd || (xub >= big);
      if (!qcond) 
        {
          step = stpmul*step;
          xlb = xub;
          xub = MIN(xlb+step,big);
          *x = xub;
          i99999 = 4; /* GET-FUNCTION-VALUE */
          *status = 1; /* TO GET-FUNCTION-VALUE */
          return;
        }

      if ((xub >= big) && !qbdd) 
        {
          *status = -1;
          *qleft = 0;
          *qhi = !qincr;
          *x = big;
          return;
        }
    }

  if ( i99999 != 6 )
    {
      dstzr(&xlb,&xub,&abstol,&reltol);
      /* IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.  */
      *status = 0;
    }

  dzror(status,x,fx,&xlo,&xhi,&qdum1,&qdum2);
  if (*status == 1) 
    {
      /* GET-FUNCTION-VALUE */
      i99999 = 6;
    }
  else
    {
      *x = xlo;
      *status = 0;
    }
  return;
DSTINV:
  small = *zsmall;
  big = *zbig;
  absstp = *zabsst;
  relstp = *zrelst;
  stpmul = *zstpmu;
  abstol = *zabsto;
  reltol = *zrelto;
  return;
}

/*
 * 
 * void dinvr(int *status,double *x,double *fx,
 * unsigned long *qleft,unsigned long *qhi)
 * 
 * Double precision
 * bounds the zero of the function and invokes zror
 * Reverse Communication
 * 
 * 
 * Function
 * 
 * 
 * Bounds the    function  and  invokes  ZROR   to perform the   zero
 * finding.  STINVR  must  have   been  called  before this   routine
 * in order to set its parameters.
 * 
 * 
 * Arguments
 * 
 * 
 * STATUS <--> At the beginning of a zero finding problem, STATUS
 * should be set to 0 and INVR invoked.  (The value
 * of parameters other than X will be ignored on this cal
 * 
 * When INVR needs the function evaluated, it will set
 * STATUS to 1 and return.  The value of the function
 * should be set in FX and INVR again called without
 * changing any of its other parameters.
 * 
 * When INVR has finished without error, it will return
 * with STATUS 0.  In that case X is approximately a root
 * of F(X).
 * 
 * If INVR cannot bound the function, it returns status
 * -1 and sets QLEFT and QHI.
 * INTEGER STATUS
 * 
 * X <-- The value of X at which F(X) is to be evaluated.
 * DOUBLE PRECISION X
 * 
 * FX --> The value of F(X) calculated when INVR returns with
 * STATUS = 1.
 * DOUBLE PRECISION FX
 * 
 * QLEFT <-- Defined only if QMFINV returns .FALSE.  In that
 * case it is .TRUE. If the stepping search terminated
 * unsucessfully at SMALL.  If it is .FALSE. the search
 * terminated unsucessfully at BIG.
 * QLEFT is LOGICAL
 * 
 * QHI <-- Defined only if QMFINV returns .FALSE.  In that
 * case it is .TRUE. if F(X) .GT. Y at the termination
 * of the search and .FALSE. if F(X) .LT. Y at the
 * termination of the search.
 * QHI is LOGICAL
 * 
 */
static void dinvr(int *status,double *x,double *fx,
                  unsigned long *qleft,unsigned long *qhi)
{
  E0000(0,status,x,fx,qleft,qhi,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
}

/*
 * void dstinv(double *zsmall,double *zbig,double *zabsst,
 * double *zrelst,double *zstpmu,double *zabsto,
 * double *zrelto)
 * 
 * Double Precision - SeT INverse finder - Reverse Communication
 * Function
 * Concise Description - Given a monotone function F finds X
 * such that F(X) = Y.  Uses Reverse communication -- see invr.
 * This routine sets quantities needed by INVR.
 * More Precise Description of INVR -
 * F must be a monotone function, the results of QMFINV are
 * otherwise undefined.  QINCR must be .TRUE. if F is non-
 * decreasing and .FALSE. if F is non-increasing.
 * QMFINV will return .TRUE. if and only if F(SMALL) and
 * F(BIG) bracket Y, i. e.,
 * QINCR is .TRUE. and F(SMALL).LE.Y.LE.F(BIG) or
 * QINCR is .FALSE. and F(BIG).LE.Y.LE.F(SMALL)
 * if QMFINV returns .TRUE., then the X returned satisfies
 * the following condition.  let
 * TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
 * then if QINCR is .TRUE.,
 * F(X-TOL(X)) .LE. Y .LE. F(X+TOL(X))
 * and if QINCR is .FALSE.
 * F(X-TOL(X)) .GE. Y .GE. F(X+TOL(X))
 * Arguments
 * SMALL --> The left endpoint of the interval to be
 * searched for a solution.
 * SMALL is DOUBLE PRECISION
 * BIG --> The right endpoint of the interval to be
 * searched for a solution.
 * BIG is DOUBLE PRECISION
 * ABSSTP, RELSTP --> The initial step size in the search
 * is MAX(ABSSTP,RELSTP*ABS(X)). See algorithm.
 * ABSSTP is DOUBLE PRECISION
 * RELSTP is DOUBLE PRECISION
 * STPMUL --> When a step doesn't bound the zero, the step
 * size is multiplied by STPMUL and another step
 * taken.  A popular value is 2.0
 * DOUBLE PRECISION STPMUL
 * ABSTOL, RELTOL --> Two numbers that determine the accuracy
 * of the solution.  See function for a precise definition.
 * ABSTOL is DOUBLE PRECISION
 * RELTOL is DOUBLE PRECISION
 * Method
 * Compares F(X) with Y for the input value of X then uses QINCR
 * to determine whether to step left or right to bound the
 * desired x.  the initial step size is
 * MAX(ABSSTP,RELSTP*ABS(S)) for the input value of X.
 * Iteratively steps right or left until it bounds X.
 * At each step which doesn't bound X, the step size is doubled.
 * The routine is careful never to step beyond SMALL or BIG.  If
 * it hasn't bounded X at SMALL or BIG, QMFINV returns .FALSE.
 * after setting QLEFT and QHI.
 * If X is successfully bounded then Algorithm R of the paper
 * 'Two Efficient Algorithms with Guaranteed Convergence for
 * Finding a Zero of a Function' by J. C. P. Bus and
 * T. J. Dekker in ACM Transactions on Mathematical
 * Software, Volume 1, No. 4 page 330 (DEC. '75) is employed
 * to find the zero of the function F(X)-Y. This is routine
 * QRZERO.
 */
static void dstinv(double *zsmall,double *zbig,double *zabsst,
                   double *zrelst,double *zstpmu,double *zabsto,
                   double *zrelto)
{
  E0000(1,NULL,NULL,NULL,NULL,NULL,zabsst,zabsto,zbig,zrelst,zrelto,zsmall,
        zstpmu);
}

/* Double precision Initalize Approximation to
 * INVerse of the cumulative T distribution
 *  
 * Returns  the  inverse   of  the T   distribution   function, i.e.,
 * the integral from 0 to INVT of the T density is P. This is an
 * initial approximation
 *  
 * P --> The p-value whose inverse from the T distribution is
 * desired.
 * Q --> 1-P.
 * DF --> Degrees of freedom of the T distribution.
 */
static double dt1(double *p,double *q,double *df)
{
  double denpow,sum,term,x,xx;
  int i;
  double coef[4][5] = {
      {1.,1.,0.,0.,0.},
      {3.,16.,5.,0.,0.},
      {-15.,17., 19.,3.,0.},
      {-945.,-1920.,1482.,776.,79.}
  };
  double denom[4] = { 4.,96.,384.,92160. };
  int ideg[4] = { 2,3,4,5 };

  x = fabs (dinvnr (p, q) );
  xx = x * x;
  sum = x;
  denpow = 1.;
  for(i=0; i<4; i++)
    {
      term = devlpl (&coef[i][0], &ideg[i], &xx) * x;
      denpow *= *df;
      sum += (term / (denpow * denom[i]) );
    }
  if (*p >= 0.5) return sum;
  return -sum;
}

/* DEFINE DZROR */
static void E0001(int IENTRY,int *status,double *x,double *fx,
                  double *xlo,double *xhi,unsigned long *qleft,
                  unsigned long *qhi,double *zabstl,double *zreltl,
                  double *zxhi,double *zxlo)
{
#define ftol(zx) (0.5e0*fifdmax1(abstol,reltol*fabs((zx))))
  static double a,abstol,b,c,d,fa,fb,fc,fd,fda,fdb,m,mb,p,q,reltol,tol,w,xxhi,xxlo;
  static int ext,i99999;
  static unsigned long first,qrzero;
  switch(IENTRY){case 0: goto DZROR; case 1: goto DSTZR;}
DZROR:
  if (*status > 0) goto S280;
  *xlo = xxlo;
  *xhi = xxhi;
  b = *x = *xlo;
  /*
     GET-FUNCTION-VALUE
     */
  i99999 = 1;
  goto S270;
S10:
  fb = *fx;
  *xlo = *xhi;
  a = *x = *xlo;
  /*
     GET-FUNCTION-VALUE
     */
  i99999 = 2;
  goto S270;
S20:
  /*
     Check that F(ZXLO) < 0 < F(ZXHI)  or
     F(ZXLO) > 0 > F(ZXHI)
     */
  if (!(fb < 0.)) goto S40;
  if (!(*fx < 0.)) goto S30;
  *status = -1;
  *qleft = *fx < fb;
  *qhi = 0;
  return;
S40:
S30:
  if (!(fb > 0.)) goto S60;
  if (!(*fx > 0.)) goto S50;
  *status = -1;
  *qleft = *fx > fb;
  *qhi = 1;
  return;
S60:
S50:
  fa = *fx;
  first = 1;
S70:
  c = a;
  fc = fa;
  ext = 0;
S80:
  if (!(fabs(fc) < fabs(fb))) goto S100;
  if (!(c != a)) goto S90;
  d = a;
  fd = fa;
S90:
  a = b;
  fa = fb;
  *xlo = c;
  b = *xlo;
  fb = fc;
  c = a;
  fc = fa;
S100:
  tol = ftol(*xlo);
  m = (c+b)*.5e0;
  mb = m-b;
  if (!(fabs(mb) > tol)) goto S240;
  if (!(ext > 3)) goto S110;
  w = mb;
  goto S190;
S110:
  tol = fifdsign(tol,mb);
  p = (b-a)*fb;
  if (!first) goto S120;
  q = fa-fb;
  first = 0;
  goto S130;
S120:
  fdb = (fd-fb)/(d-b);
  fda = (fd-fa)/(d-a);
  p = fda*p;
  q = fdb*fa-fda*fb;
S130:
  if (!(p < 0.)) goto S140;
  p = -p;
  q = -q;
S140:
  if (ext == 3) p *= 2.;
  if (!(p*1. == 0. || p <= q*tol)) goto S150;
  w = tol;
  goto S180;
S150:
  if (!(p < mb*q)) goto S160;
  w = p/q;
  goto S170;
S160:
  w = mb;
S190:
S180:
S170:
  d = a;
  fd = fa;
  a = b;
  fa = fb;
  b += w;
  *xlo = b;
  *x = *xlo;
  /*
     GET-FUNCTION-VALUE
     */
  i99999 = 3;
  goto S270;
S200:
  fb = *fx;
  if (!(fc*fb >= 0.)) goto S210;
  goto S70;
S210:
  if (!(w == mb)) goto S220;
  ext = 0;
  goto S230;
S220:
  ext += 1;
S230:
  goto S80;
S240:
  *xhi = c;
  qrzero = ( fc >= 0. && fb <= 0.) || ( fc < 0. && fb >= 0.);
  if (!qrzero) goto S250;
  *status = 0;
  goto S260;
S250:
  *status = -1;
S260:
  return;
DSTZR:
  xxlo = *zxlo;
  xxhi = *zxhi;
  abstol = *zabstl;
  reltol = *zreltl;
  return;
S270:
  /*
     TO GET-FUNCTION-VALUE
     */
  *status = 1;
  return;
S280:
  switch((int)i99999){case 1: goto S10;case 2: goto S20;case 3: goto S200;
  default: break;}
#undef ftol
}

/*
 * 
 * 
 * void dzror(int *status,double *x,double *fx,double *xlo,
 * double *xhi,unsigned long *qleft,unsigned long *qhi)
 * 
 * Double precision ZeRo of a function -- Reverse Communication
 * 
 * 
 * Function
 * 
 * 
 * Performs the zero finding.  STZROR must have been called before
 * this routine in order to set its parameters.
 * 
 * 
 * Arguments
 * 
 * 
 * STATUS <--> At the beginning of a zero finding problem, STATUS
 * should be set to 0 and ZROR invoked.  (The value
 * of other parameters will be ignored on this call.)
 * 
 * When ZROR needs the function evaluated, it will set
 * STATUS to 1 and return.  The value of the function
 * should be set in FX and ZROR again called without
 * changing any of its other parameters.
 * 
 * When ZROR has finished without error, it will return
 * with STATUS 0.  In that case (XLO,XHI) bound the answe
 * 
 * If ZROR finds an error (which implies that F(XLO)-Y an
 * F(XHI)-Y have the same sign, it returns STATUS -1.  In
 * this case, XLO and XHI are undefined.
 * INTEGER STATUS
 * 
 * X <-- The value of X at which F(X) is to be evaluated.
 * DOUBLE PRECISION X
 * 
 * FX --> The value of F(X) calculated when ZROR returns with
 * STATUS = 1.
 * DOUBLE PRECISION FX
 * 
 * XLO <-- When ZROR returns with STATUS = 0, XLO bounds the
 * inverval in X containing the solution below.
 * DOUBLE PRECISION XLO
 * 
 * XHI <-- When ZROR returns with STATUS = 0, XHI bounds the
 * inverval in X containing the solution above.
 * DOUBLE PRECISION XHI
 * 
 * QLEFT <-- .TRUE. if the stepping search terminated unsucessfully
 * at XLO.  If it is .FALSE. the search terminated
 * unsucessfully at XHI.
 * QLEFT is LOGICAL
 * 
 * QHI <-- .TRUE. if F(X) .GT. Y at the termination of the
 * search and .FALSE. if F(X) .LT. Y at the
 * termination of the search.
 * QHI is LOGICAL
 * 
 */
static void dzror(int *status,double *x,double *fx,double *xlo,
                  double *xhi,unsigned long *qleft,unsigned long *qhi)
{
  E0001(0,status,x,fx,xlo,xhi,qleft,qhi,NULL,NULL,NULL,NULL);
}

/*
 * 
 * void dstzr(double *zxlo,double *zxhi,double *zabstl,double *zreltl)
 * Double precision SeT ZeRo finder - Reverse communication version
 * Function
 * Sets quantities needed by ZROR.  The function of ZROR
 * and the quantities set is given here.
 * Concise Description - Given a function F
 * find XLO such that F(XLO) = 0.
 * More Precise Description -
 * Input condition. F is a double precision function of a single
 * double precision argument and XLO and XHI are such that
 * F(XLO)*F(XHI)  .LE.  0.0
 * If the input condition is met, QRZERO returns .TRUE.
 * and output values of XLO and XHI satisfy the following
 * F(XLO)*F(XHI)  .LE. 0.
 * ABS(F(XLO)  .LE. ABS(F(XHI)
 * ABS(XLO-XHI)  .LE. TOL(X)
 * where
 * TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
 * If this algorithm does not find XLO and XHI satisfying
 * these conditions then QRZERO returns .FALSE.  This
 * implies that the input condition was not met.
 * Arguments
 * XLO --> The left endpoint of the interval to be
 * searched for a solution.
 * XLO is DOUBLE PRECISION
 * XHI --> The right endpoint of the interval to be
 * for a solution.
 * XHI is DOUBLE PRECISION
 * ABSTOL, RELTOL --> Two numbers that determine the accuracy
 * of the solution.  See function for a
 * precise definition.
 * ABSTOL is DOUBLE PRECISION
 * RELTOL is DOUBLE PRECISION
 * Method
 * Algorithm R of the paper 'Two Efficient Algorithms with
 * Guaranteed Convergence for Finding a Zero of a Function'
 * by J. C. P. Bus and T. J. Dekker in ACM Transactions on
 * Mathematical Software, Volume 1, no. 4 page 330
 * (Dec. '75) is employed to find the zero of F(X)-Y.
 * 
 */
static void dstzr(double *zxlo,double *zxhi,double *zabstl,double *zreltl)
{
  E0001(1,NULL,NULL,NULL,NULL,NULL,NULL,NULL,zabstl,zreltl,zxhi,zxlo);
}

/*
   evaluation of the real error function
   */
static double erf1(double *x)
{
  return pnl_sf_erf (*x);
}

/*
 * EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION
 * 
 * ERFC1(IND,X) = ERFC(X)            IF IND = 0
 * ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE
 */
static double erfc1(int *ind, double *x)
{
  double c = .564189583547756e0;
  double a[5] = {
    .771058495001320e-04,-.133733772997339e-02,.323076579225834e-01,
    .479137145607681e-01,.128379167095513e+00
  };
  double b[3] = {
    .301048631703895e-02,.538971687740286e-01,.375795757275549e+00
  };
  double p[8] = {
    -1.36864857382717e-07,5.64195517478974e-01,7.21175825088309e+00,
    4.31622272220567e+01,1.52989285046940e+02,3.39320816734344e+02,
    4.51918953711873e+02,3.00459261020162e+02
  };
  double q[8] = {
    1.00000000000000e+00,1.27827273196294e+01,7.70001529352295e+01,
    2.77585444743988e+02,6.38980264465631e+02,9.31354094850610e+02,
    7.90950925327898e+02,3.00459260956983e+02
  };
  double r[5] = {
    2.10144126479064e+00,2.62370141675169e+01,2.13688200555087e+01,
    4.65807828718470e+00,2.82094791773523e-01
  };
  double s[4] = {
    9.41537750555460e+01,1.87114811799590e+02,9.90191814623914e+01,
    1.80124575948747e+01
  };
  int K1 = 1;
  static double erfc1_0,ax,bot,e,t,top,w;
  /*
     ABS(X) .LE. 0.5
     */
  ax = fabs(*x);
  if (ax > 0.5e0) goto S10;
  t = *x**x;
  top = (((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4]+1.;
  bot = ((b[0]*t+b[1])*t+b[2])*t+1.;
  erfc1_0 = 0.5e0+(0.5e0-*x*(top/bot));
  if (*ind != 0) erfc1_0 = exp(t)*erfc1_0;
  return erfc1_0;
S10:
  /*
     0.5 .LT. ABS(X) .LE. 4
     */
  if (ax > 4.) goto S20;
  top = ((((((p[0]*ax+p[1])*ax+p[2])*ax+p[3])*ax+p[4])*ax+p[5])*ax+p[6])*ax+p[
    7];
  bot = ((((((q[0]*ax+q[1])*ax+q[2])*ax+q[3])*ax+q[4])*ax+q[5])*ax+q[6])*ax+q[
    7];
  erfc1_0 = top/bot;
  goto S40;
S20:
  /*
     ABS(X) .GT. 4
     */
  if (*x <= -5.6e0) goto S60;
  if (*ind != 0) goto S30;
  if (*x > 100.) goto S70;
  if (*x**x > -exparg(&K1)) goto S70;
S30:
  t = pow(1./ *x,2.0);
  top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t+r[4];
  bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t+1.;
  erfc1_0 = (c-t*top/bot)/ax;
S40:
  /*
     FINAL ASSEMBLY
     */
  if (*ind == 0) goto S50;
  if (*x < 0.) erfc1_0 = 2.*exp(*x**x)-erfc1_0;
  return erfc1_0;
S50:
  w = *x**x;
  t = w;
  e = w-t;
  erfc1_0 = (0.5e0+(0.5e0-e))*exp(-t)*erfc1_0;
  if (*x < 0.) erfc1_0 = 2.-erfc1_0;
  return erfc1_0;
S60:
  /*
     LIMIT VALUE FOR LARGE NEGATIVE X
     */
  erfc1_0 = 2.;
  if (*ind != 0) erfc1_0 = 2.*exp(*x**x);
  return erfc1_0;
S70:
  /*
     LIMIT VALUE FOR LARGE POSITIVE X
     WHEN IND = 0
     */
  return 0.;

}

/*
 * Computes exp(mu + x)
 */
static double esum(int *mu,double *x)
{
  double w;

  if (*x > 0.)
    {
      if (*mu > 0) return exp ((double) *mu) * exp (*x);
      w = (double)*mu+*x;
      if (w < 0.) return exp ((double) *mu) * exp (*x);
      return exp(w);
    }

  if (*mu < 0) return exp ((double) *mu) * exp (*x);
  w = (double)*mu + *x;
  if (w > 0.) return exp ((double) *mu) * exp (*x); 
  return exp (w);
}

/*
 * if l = 0 then  exparg(l) = the largest positive w for which
 * exp(w) can be computed.
 *
 * if l  nonzero then  exparg(l) = the largest negative w for
 * which the computed value of exp(w) is nonzero.
 *
 * note... only an approximate value for exparg(l) is needed.
 */
static double exparg(int *l)
{
  double lnb;
  int b,m;

  b = pnl_ipmpar(4);

  if (b == 2) lnb = .69314718055995;
  else if (b == 8) lnb = 2.0794415416798;
  else if (b == 16) lnb = 2.7725887222398;
  else lnb = log ((double) b);

  if (*l == 0)
    {
      m = pnl_ipmpar(10);
      return 0.99999 * ((double)m * lnb);
    }
  m = pnl_ipmpar(9) - 1;
  return 0.99999 * ((double)m * lnb);
}

/*
 * EVALUATION OF I_X (A,B)
 * 
 * FOR B .LT. MIN(EPS,EPS*A) AND X .LE. 0.5.
 * 
 * SET  FPSER = X**A
 */
static double fpser(double *a,double *b,double *x,double *eps)
{
  int K1 = 1;
  double fpser_0,an,c,s,t,tol;

  fpser_0 = 1.;
  if (*a > 1.e-3**eps) 
    {
      fpser_0 = 0.;
      t = *a*log(*x);
      if (t < exparg(&K1)) return fpser_0;
      fpser_0 = exp(t);
    }
  /*
     NOTE THAT 1/B(A,B) = B
     */
  fpser_0 = *b/ *a*fpser_0;
  tol = *eps/ *a;
  an = *a+1.;
  t = *x;
  s = t/an;

  do
    {
      an += 1.;
      t = *x*t;
      c = t/an;
      s += c;
    }
  while  (fabs(c) > tol) ;
  fpser_0 *= (1.+*a*s);
  return fpser_0;
}

/*
 * COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 .LE. A .LE. 1.5
 */
static double gam1(double *a)
{
  double s1 = .273076135303957e+00;
  double s2 = .559398236957378e-01;
  double p[7] = {
    .577215664901533e+00,-.409078193005776e+00,-.230975380857675e+00,
    .597275330452234e-01,.766968181649490e-02,-.514889771323592e-02,
    .589597428611429e-03
  };
  double q[5] = {
    .100000000000000e+01,.427569613095214e+00,.158451672430138e+00,
    .261132021441447e-01,.423244297896961e-02
  };
  double r[9] = {
    -.422784335098468e+00,-.771330383816272e+00,-.244757765222226e+00,
    .118378989872749e+00,.930357293360349e-03,-.118290993445146e-01,
    .223047661158249e-02,.266505979058923e-03,-.132674909766242e-03
  };
  double bot,d,t,top,w,T1;
  /*
     ..
     .. Executable Statements ..
     */
  t = *a;
  d = *a-0.5e0;
  if (d > 0.) t = d-0.5e0;
  T1 = t;
  if (T1 < 0) 
    {
      top = (((((((r[8]*t+r[7])*t+r[6])*t+r[5])*t+r[4])*t+r[3])*t+r[2])*t+r[1])*t+
        r[0];
      bot = (s2*t+s1)*t+1.;
      w = top/bot;
      if (d > 0.) 
        return t*w/ *a;
      else
        return *a*(w+0.5e0+0.5e0);
    }
  else if (T1 == 0) return 0.;
  else 
    {
      top = (((((p[6]*t+p[5])*t+p[4])*t+p[3])*t+p[2])*t+p[1])*t+p[0];
      bot = (((q[4]*t+q[3])*t+q[2])*t+q[1])*t+1.;
      w = top/bot;
      if (d > 0.) 
        return t/ *a*(w-0.5e0-0.5e0);
      else
        return *a*w;
    }
}

/*
 * INVERSE INCOMPLETE GAMMA RATIO FUNCTION
 * 
 * GIVEN POSITIVE A, AND NONEGATIVE P AND Q WHERE P + Q = 1.
 * THEN X IS COMPUTED WHERE P(A,X) = P AND Q(A,X) = Q. SCHRODER
 * ITERATION IS EMPLOYED. THE ROUTINE ATTEMPTS TO COMPUTE X
 * TO 10 SIGNIFICANT DIGITS IF THIS IS POSSIBLE FOR THE
 * PARTICULAR COMPUTER ARITHMETIC BEING USED.
 * 
 * ------------
 * 
 * X IS A VARIABLE. IF P = 0 THEN X IS ASSIGNED THE VALUE 0,
 * AND IF Q = 0 THEN X IS SET TO THE LARGEST FLOATING POINT
 * NUMBER AVAILABLE. OTHERWISE, GAMINV ATTEMPTS TO OBTAIN
 * A SOLUTION FOR P(A,X) = P AND Q(A,X) = Q. IF THE ROUTINE
 * IS SUCCESSFUL THEN THE SOLUTION IS STORED IN X.
 * 
 * X0 IS AN OPTIONAL INITIAL APPROXIMATION FOR X. IF THE USER
 * DOES NOT WISH TO SUPPLY AN INITIAL APPROXIMATION, THEN SET
 * X0 .LE. 0.
 * 
 * IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
 * WHEN THE ROUTINE TERMINATES, IERR HAS ONE OF THE FOLLOWING
 * VALUES ...
 * 
 * IERR =  0    THE SOLUTION WAS OBTAINED. ITERATION WAS
 * NOT USED.
 * IERR.GT.0    THE SOLUTION WAS OBTAINED. IERR ITERATIONS
 * WERE PERFORMED.
 * IERR = -2    (INPUT ERROR) A .LE. 0
 * IERR = -3    NO SOLUTION WAS OBTAINED. THE RATIO Q/A
 * IS TOO LARGE.
 * IERR = -4    (INPUT ERROR) P + Q .NE. 1
 * IERR = -6    20 ITERATIONS WERE PERFORMED. THE MOST
 * RECENT VALUE OBTAINED FOR X IS GIVEN.
 * THIS CANNOT OCCUR IF X0 .LE. 0.
 * IERR = -7    ITERATION FAILED. NO VALUE IS GIVEN FOR X.
 * THIS MAY OCCUR WHEN X IS APPROXIMATELY 0.
 * IERR = -8    A VALUE FOR X HAS BEEN OBTAINED, BUT THE
 * ROUTINE IS NOT CERTAIN OF ITS ACCURACY.
 * ITERATION CANNOT BE PERFORMED IN THIS
 * CASE. IF X0 .LE. 0, THIS CAN OCCUR ONLY
 * WHEN P OR Q IS APPROXIMATELY 0. IF X0 IS
 * POSITIVE THEN THIS CAN OCCUR WHEN A IS
 * EXCEEDINGLY CLOSE TO X AND A IS EXTREMELY
 * LARGE (SAY A .GE. 1.E20).
 * ----------------------------------------------------------------------
 * WRITTEN BY ALFRED H. MORRIS, JR.
 * NAVAL SURFACE WEAPONS CENTER
 * DAHLGREN, VIRGINIA
 */
static void gaminv(double *a,double *x,double *x0,double *p,double *q,
                   int *ierr)
{
  double a0 = 3.31125922108741e0;
  double a1 = 11.6616720288968e0;
  double a2 = 4.28342155967104e0;
  double a3 = .213623493715853e0;
  double b1 = 6.61053765625462e0;
  double b2 = 6.40691597760039e0;
  double b3 = 1.27364489782223e0;
  double b4 = .036117081018842e0;
  double c = .577215664901533e0;
  double ln10 = 2.302585e0;
  double tol = 1.e-5;
  double amin[2] = {
    500.,100.
  };
  double bmin[2] = {
    1.e-28,1.e-13
  };
  double dmin[2] = {
    1.e-06,1.e-04
  };
  double emin[2] = {
    2.e-03,6.e-03
  };
  double eps0[2] = {
    1.e-10,1.e-08
  };
  int K1 = 1;
  int K2 = 2;
  int K3 = 3;
  int K8 = 0;
  static double am1,amax,ap1,ap2,ap3,apn,b,c1,c2,c3,c4,c5,d,e,e2,eps,g,h,pn,qg,qn,
                r,rta,s,s2,sum,t,u,w,xmax,xmin,xn,y,z;
  int iop;
  double T4,T5,T6,T7,T9;
  /*
     ..
     .. Executable Statements ..
     */
  /*
   ****** E, XMIN, AND XMAX ARE MACHINE DEPENDENT CONSTANTS.
   E IS THE SMALLEST NUMBER FOR WHICH 1.0 + E .GT. 1.0.
   XMIN IS THE SMALLEST POSITIVE NUMBER AND XMAX IS THE
   LARGEST POSITIVE NUMBER.
   */
  e = spmpar(&K1);
  xmin = spmpar(&K2);
  xmax = spmpar(&K3);
  *x = 0.;
  if (*a <= 0.) goto S300;
  t = *p+*q-1.e0;
  if (fabs(t) > e) goto S320;
  *ierr = 0;
  if (*p == 0.) return;
  if (*q == 0.) goto S270;
  if (*a == 1.) goto S280;
  e2 = 2.*e;
  amax = 0.4e-10/(e*e);
  iop = 1;
  if (e > 1.e-10) iop = 2;
  eps = eps0[iop-1];
  xn = *x0;
  if (*x0 > 0.) goto S160;
  /*
     SELECTION OF THE INITIAL APPROXIMATION XN OF X
     WHEN A .LT. 1
     */
  if (*a > 1.) goto S80;
  T4 = *a+1.;
  g = Xgamm(&T4);
  qg = *q*g;
  if (qg == 0.) goto S360;
  b = qg/ *a;
  if (qg > 0.6e0**a) goto S40;
  if (*a >= 0.30e0 || b < 0.35e0) goto S10;
  t = exp(-(b+c));
  u = t*exp(t);
  xn = t*exp(u);
  goto S160;
S10:
  if (b >= 0.45e0) goto S40;
  if (b == 0.) goto S360;
  y = -log(b);
  s = 0.5e0+(0.5e0-*a);
  z = log(y);
  t = y-s*z;
  if (b < 0.15e0) goto S20;
  xn = y-s*log(t)-log(1.+s/(t+1.));
  goto S220;
S20:
  if (b <= 0.01e0) goto S30;
  u = ((t+2.*(3.-*a))*t+(2.-*a)*(3.-*a))/((t+(5.-*a))*t+2.);
  xn = y-s*log(t)-log(u);
  goto S220;
S30:
  c1 = -(s*z);
  c2 = -(s*(1.+c1));
  c3 = s*((0.5e0*c1+(2.-*a))*c1+(2.5e0-1.5e0**a));
  c4 = -(s*(((c1/3.+(2.5e0-1.5e0**a))*c1+((*a-6.)**a+7.))*c1+(
                                                              (11.**a-46.0)**a+47.)/6.));
  c5 = -(s*((((-(c1/4.)+(11.**a-17.)/6.)*c1+((-(3.**a)+13.)*
                                             *a-13.))*c1+0.5e0*(((2.**a-25.)**a+72.)**a-61.))*c1+((
                                                                                                   (25.**a-195.)**a+477.)**a-379.)/12.));
  xn = (((c5/y+c4)/y+c3)/y+c2)/y+c1+y;
  if (*a > 1.) goto S220;
  if (b > bmin[iop-1]) goto S220;
  *x = xn;
  return;
S40:
  if (b**q > 1.e-8) goto S50;
  xn = exp(-(*q/ *a+c));
  goto S70;
S50:
  if (*p <= 0.9e0) goto S60;
  T5 = -*q;
  xn = exp((alnrel(&T5)+gamln1(a))/ *a);
  goto S70;
S60:
  xn = exp(log(*p*g)/ *a);
S70:
  if (xn == 0.) goto S310;
  t = 0.5e0+(0.5e0-xn/(*a+1.));
  xn /= t;
  goto S160;
S80:
  /*
     SELECTION OF THE INITIAL APPROXIMATION XN OF X
     WHEN A .GT. 1
     */
  if (*q <= 0.5e0) goto S90;
  w = log(*p);
  goto S100;
S90:
  w = log(*q);
S100:
  t = sqrt(-(2.*w));
  s = t-(((a3*t+a2)*t+a1)*t+a0)/((((b4*t+b3)*t+b2)*t+b1)*t+1.);
  if (*q > 0.5e0) s = -s;
  rta = sqrt(*a);
  s2 = s*s;
  xn = *a+s*rta+(s2-1.)/3.+s*(s2-7.)/(36.*rta)-((3.*s2+7.)*
                                                s2-16.)/(810.**a)+s*((9.*s2+256.)*s2-433.)/(38880.**a*
                                                                                            rta);
  xn = fifdmax1(xn,0.);
  if (*a < amin[iop-1]) goto S110;
  *x = xn;
  d = 0.5e0+(0.5e0-*x/ *a);
  if (fabs(d) <= dmin[iop-1]) return;
S110:
  if (*p <= 0.5e0) goto S130;
  if (xn < 3.**a) goto S220;
  y = -(w+gamln(a));
  d = fifdmax1(2.,*a*(*a-1.));
  if (y < ln10*d) goto S120;
  s = 1.-*a;
  z = log(y);
  goto S30;
S120:
  t = *a-1.;
  T6 = -(t/(xn+1.));
  xn = y+t*log(xn)-alnrel(&T6);
  T7 = -(t/(xn+1.));
  xn = y+t*log(xn)-alnrel(&T7);
  goto S220;
S130:
  ap1 = *a+1.;
  if (xn > 0.70e0*ap1) goto S170;
  w += gamln(&ap1);
  if (xn > 0.15e0*ap1) goto S140;
  ap2 = *a+2.;
  ap3 = *a+3.;
  *x = exp((w+*x)/ *a);
  *x = exp((w+*x-log(1.+*x/ap1*(1.+*x/ap2)))/ *a);
  *x = exp((w+*x-log(1.+*x/ap1*(1.+*x/ap2)))/ *a);
  *x = exp((w+*x-log(1.+*x/ap1*(1.+*x/ap2*(1.+*x/ap3))))/ *a);
  xn = *x;
  if (xn > 1.e-2*ap1) goto S140;
  if (xn <= emin[iop-1]*ap1) return;
  goto S170;
S140:
  apn = ap1;
  t = xn/apn;
  sum = 1.+t;
S150:
  apn += 1.;
  t *= (xn/apn);
  sum += t;
  if (t > 1.e-4) goto S150;
  t = w-log(sum);
  xn = exp((xn+t)/ *a);
  xn *= (1.-(*a*log(xn)-xn-t)/(*a-xn));
  goto S170;
S160:
  /*
     SCHRODER ITERATION USING P
     */
  if (*p > 0.5e0) goto S220;
S170:
  if (*p <= 1.e10*xmin) goto S350;
  am1 = *a-0.5e0-0.5e0;
S180:
  if (*a <= amax) goto S190;
  d = 0.5e0+(0.5e0-xn/ *a);
  if (fabs(d) <= e2) goto S350;
S190:
  if (*ierr >= 20) goto S330;
  *ierr += 1;
  gratio(a,&xn,&pn,&qn,&K8);
  if (pn == 0. || qn == 0.) goto S350;
  r = rcomp(a,&xn);
  if (r == 0.) goto S350;
  t = (pn-*p)/r;
  w = 0.5e0*(am1-xn);
  if (fabs(t) <= 0.1e0 && fabs(w*t) <= 0.1e0) goto S200;
  *x = xn*(1.-t);
  if (*x <= 0.) goto S340;
  d = fabs(t);
  goto S210;
S200:
  h = t*(1.+w*t);
  *x = xn*(1.-h);
  if (*x <= 0.) goto S340;
  if (fabs(w) >= 1. && fabs(w)*t*t <= eps) return;
  d = fabs(h);
S210:
  xn = *x;
  if (d > tol) goto S180;
  if (d <= eps) return;
  if (fabs(*p-pn) <= tol**p) return;
  goto S180;
S220:
  /*
     SCHRODER ITERATION USING Q
     */
  if (*q <= 1.e10*xmin) goto S350;
  am1 = *a-0.5e0-0.5e0;
S230:
  if (*a <= amax) goto S240;
  d = 0.5e0+(0.5e0-xn/ *a);
  if (fabs(d) <= e2) goto S350;
S240:
  if (*ierr >= 20) goto S330;
  *ierr += 1;
  gratio(a,&xn,&pn,&qn,&K8);
  if (pn == 0. || qn == 0.) goto S350;
  r = rcomp(a,&xn);
  if (r == 0.) goto S350;
  t = (*q-qn)/r;
  w = 0.5e0*(am1-xn);
  if (fabs(t) <= 0.1e0 && fabs(w*t) <= 0.1e0) goto S250;
  *x = xn*(1.-t);
  if (*x <= 0.) goto S340;
  d = fabs(t);
  goto S260;
S250:
  h = t*(1.+w*t);
  *x = xn*(1.-h);
  if (*x <= 0.) goto S340;
  if (fabs(w) >= 1. && fabs(w)*t*t <= eps) return;
  d = fabs(h);
S260:
  xn = *x;
  if (d > tol) goto S230;
  if (d <= eps) return;
  if (fabs(*q-qn) <= tol**q) return;
  goto S230;
S270:
  /*
     SPECIAL CASES
     */
  *x = xmax;
  return;
S280:
  if (*q < 0.9e0) goto S290;
  T9 = -*p;
  *x = -alnrel(&T9);
  return;
S290:
  *x = -log(*q);
  return;
S300:
  /*
     ERROR RETURN
     */
  *ierr = -2;
  return;
S310:
  *ierr = -3;
  return;
S320:
  *ierr = -4;
  return;
S330:
  *ierr = -6;
  return;
S340:
  *ierr = -7;
  return;
S350:
  *x = xn;
  *ierr = -8;
  return;
S360:
  *x = xmax;
  *ierr = -8;
  return;
}

/*
 * Compute the log of the gamma function
 */
static double gamln(double *a)
{
  return pnl_sf_log_gamma (*a);
}

/*
 * -----------------------------------------------------------------------
 * EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 .LE. A .LE. 1.25
 * -----------------------------------------------------------------------
 */
static double gamln1(double *a)
{
  double p0 = .577215664901533e+00;
  double p1 = .844203922187225e+00;
  double p2 = -.168860593646662e+00;
  double p3 = -.780427615533591e+00;
  double p4 = -.402055799310489e+00;
  double p5 = -.673562214325671e-01;
  double p6 = -.271935708322958e-02;
  double q1 = .288743195473681e+01;
  double q2 = .312755088914843e+01;
  double q3 = .156875193295039e+01;
  double q4 = .361951990101499e+00;
  double q5 = .325038868253937e-01;
  double q6 = .667465618796164e-03;
  double r0 = .422784335098467e+00;
  double r1 = .848044614534529e+00;
  double r2 = .565221050691933e+00;
  double r3 = .156513060486551e+00;
  double r4 = .170502484022650e-01;
  double r5 = .497958207639485e-03;
  double s1 = .124313399877507e+01;
  double s2 = .548042109832463e+00;
  double s3 = .101552187439830e+00;
  double s4 = .713309612391000e-02;
  double s5 = .116165475989616e-03;
  double w,x;

  if (*a >= 0.6e0) 
    {
      x = *a-0.5e0-0.5e0;
      w = (((((r5*x+r4)*x+r3)*x+r2)*x+r1)*x+r0)/
        (((((s5*x+s4)*x+s3)*x+s2)*x+s1)*x +1.);
      return x*w;
    }
  else
    {
      w = ((((((p6**a+p5)**a+p4)**a+p3)**a+p2)**a+p1)**a+p0)/
        ((((((q6**a+q5)**a+ q4)**a+q3)**a+q2)**a+q1)**a+1.);
      return -(*a*w);
    }
}

/*
 * Compute the gamma function
 */
static double Xgamm(double *a)
{
  return pnl_sf_gamma (*a);
}


static void grat1(double *a,double *x,double *r,double *p,double *q,
                  double *eps)
{
  int K2 = 0;
  double a2n,a2nm1,am0,an,an0,b2n,b2nm1,c,cma,g,h,j,l,sum,t,tol,w,z,T1,T3;
  /*
     ..
     .. Executable Statements ..
     */
  /*
     -----------------------------------------------------------------------
     EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS
     P(A,X) AND Q(A,X)
     IT IS ASSUMED THAT A .LE. 1.  EPS IS THE TOLERANCE TO BE USED.
     THE INPUT ARGUMENT R HAS THE VALUE E**(-X)*X**A/GAMMA(A).
     -----------------------------------------------------------------------
     */
  if (*a**x == 0.) goto S120;
  if (*a == 0.5e0) goto S100;
  if (*x < 1.1e0) goto S10;
  goto S60;
S10:
  /*
     TAYLOR SERIES FOR P(A,X)/X**A
     */
  an = 3.;
  c = *x;
  sum = *x/(*a+3.);
  tol = 0.1e0**eps/(*a+1.);
S20:
  an += 1.;
  c = -(c*(*x/an));
  t = c/(*a+an);
  sum += t;
  if (fabs(t) > tol) goto S20;
  j = *a**x*((sum/6.-0.5e0/(*a+2.))**x+1./(*a+1.));
  z = *a*log(*x);
  h = gam1(a);
  g = 1.+h;
  if (*x < 0.25e0) goto S30;
  if (*a < *x/2.59e0) goto S50;
  goto S40;
S30:
  if (z > -.13394e0) goto S50;
S40:
  w = exp(z);
  *p = w*g*(0.5e0+(0.5e0-j));
  *q = 0.5e0+(0.5e0-*p);
  return;
S50:
  l = rexp(&z);
  w = 0.5e0+(0.5e0+l);
  *q = (w*j-l)*g-h;
  if (*q < 0.) goto S90;
  *p = 0.5e0+(0.5e0-*q);
  return;
S60:
  /*
     CONTINUED FRACTION EXPANSION
     */
  a2nm1 = a2n = 1.;
  b2nm1 = *x;
  b2n = *x+(1.-*a);
  c = 1.;
S70:
  a2nm1 = *x*a2n+c*a2nm1;
  b2nm1 = *x*b2n+c*b2nm1;
  am0 = a2nm1/b2nm1;
  c += 1.;
  cma = c-*a;
  a2n = a2nm1+cma*a2n;
  b2n = b2nm1+cma*b2n;
  an0 = a2n/b2n;
  if (fabs(an0-am0) >= *eps*an0) goto S70;
  *q = *r*an0;
  *p = 0.5e0+(0.5e0-*q);
  return;
S80:
  /*
     SPECIAL CASES
     */
  *p = 0.;
  *q = 1.;
  return;
S90:
  *p = 1.;
  *q = 0.;
  return;
S100:
  if (*x >= 0.25e0) goto S110;
  T1 = sqrt(*x);
  *p = erf1(&T1);
  *q = 0.5e0+(0.5e0-*p);
  return;
S110:
  T3 = sqrt(*x);
  *q = erfc1(&K2,&T3);
  *p = 0.5e0+(0.5e0-*q);
  return;
S120:
  if (*x <= *a) goto S80;
  goto S90;
}

/*
 * ----------------------------------------------------------------------
 * EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS
 * P(A,X) AND Q(A,X)
 * 
 * ----------
 * 
 * IT IS ASSUMED THAT A AND X ARE NONNEGATIVE, WHERE A AND X
 * ARE NOT BOTH 0.
 * 
 * ANS AND QANS ARE VARIABLES. GRATIO ASSIGNS ANS THE VALUE
 * P(A,X) AND QANS THE VALUE Q(A,X). IND MAY BE ANY INTEGER.
 * IF IND = 0 THEN THE USER IS REQUESTING AS MUCH ACCURACY AS
 * POSSIBLE (UP TO 14 SIGNIFICANT DIGITS). OTHERWISE, IF
 * IND = 1 THEN ACCURACY IS REQUESTED TO WITHIN 1 UNIT OF THE
 * 6-TH SIGNIFICANT DIGIT, AND IF IND .NE. 0,1 THEN ACCURACY
 * IS REQUESTED TO WITHIN 1 UNIT OF THE 3RD SIGNIFICANT DIGIT.
 * 
 * ERROR RETURN ...
 * ANS IS ASSIGNED THE VALUE 2 WHEN A OR X IS NEGATIVE,
 * WHEN A*X = 0, OR WHEN P(A,X) AND Q(A,X) ARE INDETERMINANT.
 * P(A,X) AND Q(A,X) ARE COMPUTATIONALLY INDETERMINANT WHEN
 * X IS EXCEEDINGLY CLOSE TO A AND A IS EXTREMELY LARGE.
 * ----------------------------------------------------------------------
 * WRITTEN BY ALFRED H. MORRIS, JR.
 * NAVAL SURFACE WEAPONS CENTER
 * DAHLGREN, VIRGINIA
 * --------------------
 */
static void gratio(double *a,double *x,double *ans,double *qans,int *ind)
{
  double alog10 = 2.30258509299405e0;
  double d10 = -.185185185185185e-02;
  double d20 = .413359788359788e-02;
  double d30 = .649434156378601e-03;
  double d40 = -.861888290916712e-03;
  double d50 = -.336798553366358e-03;
  double d60 = .531307936463992e-03;
  double d70 = .344367606892378e-03;
  double rt2pin = .398942280401433e0;
  double rtpi = 1.77245385090552e0;
  double third = .333333333333333e0;
  double acc0[3] = {
    5.e-15,5.e-7,5.e-4
  };
  double big[3] = {
    20.,14.,10.
  };
  double d0[13] = {
    .833333333333333e-01,-.148148148148148e-01,.115740740740741e-02,
    .352733686067019e-03,-.178755144032922e-03,.391926317852244e-04,
    -.218544851067999e-05,-.185406221071516e-05,.829671134095309e-06,
    -.176659527368261e-06,.670785354340150e-08,.102618097842403e-07,
    -.438203601845335e-08
  };
  double d1[12] = {
    -.347222222222222e-02,.264550264550265e-02,-.990226337448560e-03,
    .205761316872428e-03,-.401877572016461e-06,-.180985503344900e-04,
    .764916091608111e-05,-.161209008945634e-05,.464712780280743e-08,
    .137863344691572e-06,-.575254560351770e-07,.119516285997781e-07
  };
  double d2[10] = {
    -.268132716049383e-02,.771604938271605e-03,.200938786008230e-05,
    -.107366532263652e-03,.529234488291201e-04,-.127606351886187e-04,
    .342357873409614e-07,.137219573090629e-05,-.629899213838006e-06,
    .142806142060642e-06
  };
  double d3[8] = {
    .229472093621399e-03,-.469189494395256e-03,.267720632062839e-03,
    -.756180167188398e-04,-.239650511386730e-06,.110826541153473e-04,
    -.567495282699160e-05,.142309007324359e-05
  };
  double d4[6] = {
    .784039221720067e-03,-.299072480303190e-03,-.146384525788434e-05,
    .664149821546512e-04,-.396836504717943e-04,.113757269706784e-04
  };
  double d5[4] = {
    -.697281375836586e-04,.277275324495939e-03,-.199325705161888e-03,
    .679778047793721e-04
  };
  double d6[2] = {
    -.592166437353694e-03,.270878209671804e-03
  };
  double e00[3] = {
    .25e-3,.25e-1,.14e0
  };
  double x00[3] = {
    31.,17.,9.7e0
  };
  int K1 = 1;
  int K2 = 0;
  static double a2n,a2nm1,acc,am0,amn,an,an0,apn,b2n,b2nm1,c,c0,c1,c2,c3,c4,c5,c6,
                cma,e,e0,g,h,j,l,r,rta,rtx,s,sum,t,t1,tol,twoa,u,w,x0,y,z;
  int i,iop,m,max,n;
  double wk[20],T3;
  int T4,T5;
  double T6,T7;
  /*
     ..
     .. Executable Statements ..
     */
  /*
     --------------------
   ****** E IS A MACHINE DEPENDENT CONSTANT. E IS THE SMALLEST
   FLOATING POINT NUMBER FOR WHICH 1.0 + E .GT. 1.0 .
   */
  e = spmpar(&K1);
  if (*a < 0. || *x < 0.) goto S430;
  if (*a == 0. && *x == 0.) goto S430;
  if (*a**x == 0.) goto S420;
  iop = *ind+1;
  if (iop != 1 && iop != 2) iop = 3;
  acc = fifdmax1(acc0[iop-1],e);
  e0 = e00[iop-1];
  x0 = x00[iop-1];
  /*
     SELECT THE APPROPRIATE ALGORITHM
     */
  if (*a >= 1.) goto S10;
  if (*a == 0.5e0) goto S390;
  if (*x < 1.1e0) goto S160;
  t1 = *a*log(*x)-*x;
  u = *a*exp(t1);
  if (u == 0.) goto S380;
  r = u*(1.+gam1(a));
  goto S250;
S10:
  if (*a >= big[iop-1]) goto S30;
  if (*a > *x || *x >= x0) goto S20;
  twoa = *a+*a;
  m = fifidint(twoa);
  if (twoa != (double)m) goto S20;
  i = m/2;
  if (*a == (double)i) goto S210;
  goto S220;
S20:
  t1 = *a*log(*x)-*x;
  r = exp(t1)/Xgamm(a);
  goto S40;
S30:
  l = *x/ *a;
  if (l == 0.) goto S370;
  s = 0.5e0+(0.5e0-l);
  z = rlog(&l);
  if (z >= 700./ *a) goto S410;
  y = *a*z;
  rta = sqrt(*a);
  if (fabs(s) <= e0/rta) goto S330;
  if (fabs(s) <= 0.4e0) goto S270;
  t = pow(1./ *a,2.0);
  t1 = (((0.75e0*t-1.)*t+3.5e0)*t-105.)/(*a*1260.);
  t1 -= y;
  r = rt2pin*rta*exp(t1);
S40:
  if (r == 0.) goto S420;
  if (*x <= fifdmax1(*a,alog10)) goto S50;
  if (*x < x0) goto S250;
  goto S100;
S50:
  /*
     TAYLOR SERIES FOR P/R
     */
  apn = *a+1.;
  t = *x/apn;
  wk[0] = t;
  for(n=2; n<=20; n++) {
    apn += 1.;
    t *= (*x/apn);
    if (t <= 1.e-3) goto S70;
    wk[n-1] = t;
  }
  n = 20;
S70:
  sum = t;
  tol = 0.5e0*acc;
S80:
  apn += 1.;
  t *= (*x/apn);
  sum += t;
  if (t > tol) goto S80;
  max = n-1;
  for(m=1; m<=max; m++) {
    n -= 1;
    sum += wk[n-1];
  }
  *ans = r/ *a*(1.+sum);
  *qans = 0.5e0+(0.5e0-*ans);
  return;
S100:
  /*
     ASYMPTOTIC EXPANSION
     */
  amn = *a-1.;
  t = amn/ *x;
  wk[0] = t;
  for(n=2; n<=20; n++) {
    amn -= 1.;
    t *= (amn/ *x);
    if (fabs(t) <= 1.e-3) goto S120;
    wk[n-1] = t;
  }
  n = 20;
S120:
  sum = t;
S130:
  if (fabs(t) <= acc) goto S140;
  amn -= 1.;
  t *= (amn/ *x);
  sum += t;
  goto S130;
S140:
  max = n-1;
  for(m=1; m<=max; m++) {
    n -= 1;
    sum += wk[n-1];
  }
  *qans = r/ *x*(1.+sum);
  *ans = 0.5e0+(0.5e0-*qans);
  return;
S160:
  /*
     TAYLOR SERIES FOR P(A,X)/X**A
     */
  an = 3.;
  c = *x;
  sum = *x/(*a+3.);
  tol = 3.*acc/(*a+1.);
S170:
  an += 1.;
  c = -(c*(*x/an));
  t = c/(*a+an);
  sum += t;
  if (fabs(t) > tol) goto S170;
  j = *a**x*((sum/6.-0.5e0/(*a+2.))**x+1./(*a+1.));
  z = *a*log(*x);
  h = gam1(a);
  g = 1.+h;
  if (*x < 0.25e0) goto S180;
  if (*a < *x/2.59e0) goto S200;
  goto S190;
S180:
  if (z > -.13394e0) goto S200;
S190:
  w = exp(z);
  *ans = w*g*(0.5e0+(0.5e0-j));
  *qans = 0.5e0+(0.5e0-*ans);
  return;
S200:
  l = rexp(&z);
  w = 0.5e0+(0.5e0+l);
  *qans = (w*j-l)*g-h;
  if (*qans < 0.) goto S380;
  *ans = 0.5e0+(0.5e0-*qans);
  return;
S210:
  /*
     FINITE SUMS FOR Q WHEN A .GE. 1
     AND 2*A IS AN INTEGER
     */
  sum = exp(-*x);
  t = sum;
  n = 1;
  c = 0.;
  goto S230;
S220:
  rtx = sqrt(*x);
  sum = erfc1(&K2,&rtx);
  t = exp(-*x)/(rtpi*rtx);
  n = 0;
  c = -0.5e0;
S230:
  if (n == i) goto S240;
  n += 1;
  c += 1.;
  t = *x*t/c;
  sum += t;
  goto S230;
S240:
  *qans = sum;
  *ans = 0.5e0+(0.5e0-*qans);
  return;
S250:
  /*
     CONTINUED FRACTION EXPANSION
     */
  tol = fifdmax1(5.*e,acc);
  a2nm1 = a2n = 1.;
  b2nm1 = *x;
  b2n = *x+(1.-*a);
  c = 1.;
S260:
  a2nm1 = *x*a2n+c*a2nm1;
  b2nm1 = *x*b2n+c*b2nm1;
  am0 = a2nm1/b2nm1;
  c += 1.;
  cma = c-*a;
  a2n = a2nm1+cma*a2n;
  b2n = b2nm1+cma*b2n;
  an0 = a2n/b2n;
  if (fabs(an0-am0) >= tol*an0) goto S260;
  *qans = r*an0;
  *ans = 0.5e0+(0.5e0-*qans);
  return;
S270:
  /*
     GENERAL TEMME EXPANSION
     */
  if (fabs(s) <= 2.*e && *a*e*e > 3.28e-3) goto S430;
  c = exp(-y);
  T3 = sqrt(y);
  w = 0.5e0*erfc1(&K1,&T3);
  u = 1./ *a;
  z = sqrt(z+z);
  if (l < 1.) z = -z;
  T4 = iop-2;
  if (T4 < 0) goto S280;
  else if (T4 == 0) goto S290;
  else  goto S300;
S280:
  if (fabs(s) <= 1.e-3) goto S340;
  c0 = ((((((((((((d0[12]*z+d0[11])*z+d0[10])*z+d0[9])*z+d0[8])*z+d0[7])*z+d0[
              6])*z+d0[5])*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-third;
  c1 = (((((((((((d1[11]*z+d1[10])*z+d1[9])*z+d1[8])*z+d1[7])*z+d1[6])*z+d1[5]
            )*z+d1[4])*z+d1[3])*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
  c2 = (((((((((d2[9]*z+d2[8])*z+d2[7])*z+d2[6])*z+d2[5])*z+d2[4])*z+d2[3])*z+
          d2[2])*z+d2[1])*z+d2[0])*z+d20;
  c3 = (((((((d3[7]*z+d3[6])*z+d3[5])*z+d3[4])*z+d3[3])*z+d3[2])*z+d3[1])*z+
        d3[0])*z+d30;
  c4 = (((((d4[5]*z+d4[4])*z+d4[3])*z+d4[2])*z+d4[1])*z+d4[0])*z+d40;
  c5 = (((d5[3]*z+d5[2])*z+d5[1])*z+d5[0])*z+d50;
  c6 = (d6[1]*z+d6[0])*z+d60;
  t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u+c0;
  goto S310;
S290:
  c0 = (((((d0[5]*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-third;
  c1 = (((d1[3]*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
  c2 = d2[0]*z+d20;
  t = (c2*u+c1)*u+c0;
  goto S310;
S300:
  t = ((d0[2]*z+d0[1])*z+d0[0])*z-third;
S310:
  if (l < 1.) goto S320;
  *qans = c*(w+rt2pin*t/rta);
  *ans = 0.5e0+(0.5e0-*qans);
  return;
S320:
  *ans = c*(w-rt2pin*t/rta);
  *qans = 0.5e0+(0.5e0-*ans);
  return;
S330:
  /*
     TEMME EXPANSION FOR L = 1
     */
  if (*a*e*e > 3.28e-3) goto S430;
  c = 0.5e0+(0.5e0-y);
  w = (0.5e0-sqrt(y)*(0.5e0+(0.5e0-y/3.))/rtpi)/c;
  u = 1./ *a;
  z = sqrt(z+z);
  if (l < 1.) z = -z;
  T5 = iop-2;
  if (T5 < 0) goto S340;
  else if (T5 == 0) goto S350;
  else  goto S360;
S340:
  c0 = ((((((d0[6]*z+d0[5])*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-
    third;
  c1 = (((((d1[5]*z+d1[4])*z+d1[3])*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
  c2 = ((((d2[4]*z+d2[3])*z+d2[2])*z+d2[1])*z+d2[0])*z+d20;
  c3 = (((d3[3]*z+d3[2])*z+d3[1])*z+d3[0])*z+d30;
  c4 = (d4[1]*z+d4[0])*z+d40;
  c5 = (d5[1]*z+d5[0])*z+d50;
  c6 = d6[0]*z+d60;
  t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u+c0;
  goto S310;
S350:
  c0 = (d0[1]*z+d0[0])*z-third;
  c1 = d1[0]*z+d10;
  t = (d20*u+c1)*u+c0;
  goto S310;
S360:
  t = d0[0]*z-third;
  goto S310;
S370:
  /*
     SPECIAL CASES
     */
  *ans = 0.;
  *qans = 1.;
  return;
S380:
  *ans = 1.;
  *qans = 0.;
  return;
S390:
  if (*x >= 0.25e0) goto S400;
  T6 = sqrt(*x);
  *ans = erf1(&T6);
  *qans = 0.5e0+(0.5e0-*ans);
  return;
S400:
  T7 = sqrt(*x);
  *qans = erfc1(&K2,&T7);
  *ans = 0.5e0+(0.5e0-*qans);
  return;
S410:
  if (fabs(s) <= 2.*e) goto S430;
S420:
  if (*x <= *a) goto S370;
  goto S380;
S430:
  /*
     ERROR RETURN
     */
  *ans = 2.;
  return;
}

/*
 * EVALUATION OF THE FUNCTION LN(GAMMA(A + B))
 * FOR 1 .LE. A .LE. 2  AND  1 .LE. B .LE. 2
 */
static double gsumln(double *a,double *b)
{
  double x,T1,T2;
  x = *a+*b-2.e0;
  if (x <= 0.25e0) 
    {
      T1 = 1.+x;
      return gamln1(&T1);
    }
  else if (x > 1.25e0)
    {
      T2 = x-1.;
      return gamln1(&T2)+log(x*(1.+x));
    }
  else
    return gamln1(&x)+alnrel(&x);
}

static double psi(double *xx)
{
  return pnl_sf_psi (*xx);
}

/*
 * EVALUATION OF EXP(-X)*X**A/GAMMA(A)
 * 
 * RT2PIN = 1/SQRT(2*M_PI)
 */
static double rcomp(double *a,double *x)
{
  double rt2pin = .398942280401433e0;
  double rcomp_0,t,t1,u;
  /*
     ..
     .. Executable Statements ..
     */
  rcomp_0 = 0.;
  if (*a >= 20.) goto S20;
  t = *a*log(*x)-*x;
  if (*a >= 1.) goto S10;
  rcomp_0 = *a*exp(t)*(1.+gam1(a));
  return rcomp_0;
S10:
  return exp(t)/Xgamm(a);

S20:
  u = *x/ *a;
  if (u == 0.) return rcomp_0;
  t = pow(1./ *a,2.0);
  t1 = (((0.75e0*t-1.)*t+3.5e0)*t-105.)/(*a*1260.);
  t1 -= (*a*rlog(&u));
  return rt2pin*sqrt(*a)*exp(t1);

}

/*
 * EVALUATION OF THE FUNCTION EXP(X) - 1
 */
static double rexp(double *x)
{
  return pnl_expm1 (*x);
}

/*
 * COMPUTATION OF  X - 1 - LN(X)
 */
static double rlog(double *x)
{
  double a = .566749439387324e-01;
  double b = .456512608815524e-01;
  double p0 = .333333333333333e+00;
  double p1 = -.224696413112536e+00;
  double p2 = .620886815375787e-02;
  double q1 = -.127408923933623e+01;
  double q2 = .354508718369557e+00;
  double r,t,u,w,w1;
  /*
     ..
     .. Executable Statements ..
     */
  if (*x < 0.61e0 || *x > 1.57e0) goto S40;
  if (*x < 0.82e0) goto S10;
  if (*x > 1.18e0) goto S20;
  /*
     ARGUMENT REDUCTION
     */
  u = *x-0.5e0-0.5e0;
  w1 = 0.;
  goto S30;
S10:
  u = *x-0.7e0;
  u /= 0.7e0;
  w1 = a-u*0.3e0;
  goto S30;
S20:
  u = 0.75e0**x-1.e0;
  w1 = b+u/3.;
S30:
  /*
     SERIES EXPANSION
     */
  r = u/(u+2.);
  t = r*r;
  w = ((p2*t+p1)*t+p0)/((q2*t+q1)*t+1.);
  return 2.*t*(1./(1.-r)-r*w)+w1;
S40:
  r = *x-0.5e0-0.5e0;
  return r-log(*x);
}

/*
 * EVALUATION OF THE FUNCTION X - LN(1 + X)
 */
static double rlog1(double *x)
{
  double a = .566749439387324e-01;
  double b = .456512608815524e-01;
  double p0 = .333333333333333e+00;
  double p1 = -.224696413112536e+00;
  double p2 = .620886815375787e-02;
  double q1 = -.127408923933623e+01;
  double q2 = .354508718369557e+00;
  double h,r,t,w,w1;
  /*
     ..
     .. Executable Statements ..
     */
  if (*x < -0.39e0 || *x > 0.57e0) goto S40;
  if (*x < -0.18e0) goto S10;
  if (*x > 0.18e0) goto S20;
  /*
     ARGUMENT REDUCTION
     */
  h = *x;
  w1 = 0.;
  goto S30;
S10:
  h = *x+0.3e0;
  h /= 0.7e0;
  w1 = a-h*0.3e0;
  goto S30;
S20:
  h = 0.75e0**x-0.25e0;
  w1 = b+h/3.;
S30:
  /*
     SERIES EXPANSION
     */
  r = h/(h+2.);
  t = r*r;
  w = ((p2*t+p1)*t+p0)/((q2*t+q1)*t+1.);
  return 2.*t*(1./(1.-r)-r*w)+w1;
S40:
  w = *x+0.5e0+0.5e0;
  return *x-log(w);
}

/*
 *     spmpar provides the single precision machine constants for 
 *     the computer being used. it is assumed that the argument 
 *     i is an int having one of the values 1, 2, or 3. if the 
 *     single precision arithmetic being used has m base b digits and 
 *     its smallest and largest exponents are emin and emax, then 
 *        spmpar(1) = b**(1 - m), the machine precision, 
 *        spmpar(2) = b**(emin - 1), the smallest magnitude, 
 *        spmpar(3) = b**emax*(1 - b**(-m)), the largest magnitude. 
 *     rewriten  by jpc to use lapack dlamch 
 */
static double spmpar (int *i)
{
  switch (*i)
    {
    case 1 : return  pnl_dlamch ("p");
    case 2 : return  pnl_dlamch ("u");
    case 3 : return  pnl_dlamch ("o");
    }
  return 0.0;
}

/*
 * 
 * double stvaln(double *p)
 * STarting VALue for Neton-Raphon
 * calculation of Normal distribution Inverse
 * 
 * Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
 * infinity to X of (1/SQRT(2*M_PI)) EXP(-U*U/2) dU is P
 * 
 * P --> The probability whose normal deviate is sought.
 * P is DOUBLE PRECISION
 * 
 * 
 * Method
 * 
 * The  rational   function   on  page 95    of Kennedy  and  Gentle,
 * Statistical Computing, Marcel Dekker, NY , 1980.
 * 
 */
static double stvaln(double *p)
{
  double xden[5] = {
    0.993484626060e-1,0.588581570495e0,0.531103462366e0,0.103537752850e0,
    0.38560700634e-2
  };
  double xnum[5] = {
    -0.322232431088e0,-1.000000000000e0,-0.342242088547e0,-0.204231210245e-1,
    -0.453642210148e-4
  };
  int K1 = 5;
  double stvaln_0,sign,y,z;

  if ((*p <= 0.5e0))
    {
      sign = -1.;
      z = *p;
    }
  else
    {
      sign = 1.;
      z = 1.-*p;
    }

  y = sqrt(-(2.*log(z)));
  stvaln_0 = y+devlpl(xnum,&K1,&y)/devlpl(xden,&K1,&y);
  return sign*stvaln_0;
}

/*
 * Truncates a double precision number to an integer and returns the
 * value in a double.
 */
static double fifdint(double a)
{
  return (double) ((int) a);
}

/*
 * returns the maximum of two numbers a and b
 */
static double fifdmax1(double a,double b)
{
  if (a < b) return b;
  else return a;
}

/*
 * returns the minimum of two numbers a and b
 */
static double fifdmin1(double a,double b)
{
  if (a < b) return a;
  else return b;
}

/*
 * transfers the sign of the variable "sign" to the variable "mag"
 */
static double fifdsign(double mag,double sign)
{
  if (mag < 0) mag = -mag;
  if (sign < 0) mag = -mag;
  return mag;
}

/*
 * Truncates a double precision number to a long integer
 */
static long fifidint(double a)
{
  if (a < 1.0) return (long) 0;
  else return (long) a;
}

/* msg - error message */
static void ftnstop(char* msg)
{
  if (msg != NULL) fprintf(stderr,"%s\n",msg);
  exit(0);
}

/*
 * Simple useable cumulative functions
 */

/**
 * Computes the CDF of Chi2(df, nc)
 *
 * @param x a real value
 * @param df the number of degrees of freedom
 * @param ncparam the non centrality parameter
 * @return CDF(x) for Chi2(df, nc)
 */
double pnl_cdfchi2n(double x, double df, double ncparam)
{
  double P,Q,r;
  int status;
  int i=1;

  pnl_cdf_chn(&i,&P,&Q,&x,&df,&ncparam,&status,&r);
  return P;
}

/**
 * Computes the cumulative distribution of beta*X where X is non centrally
 * chi squared distributed with nu degrees of freedom and non centrality
 * parameter lambda. Result is saved in variable P
 * 
 * @param x point at which the cdf is computed
 * @param nu the number of degrees of freedom
 * @param lambda the non centrality parameter
 * @param beta a real number
 * @param P contains CDF (x) for beta * X on exit
 */
void pnl_cdfbchi2n(double x, double nu, double lambda, double beta, double *P) 
{
  double bound, p, q, xx;
  int which=1, status;

  xx=x/beta;
  if (xx<=0)
    *P=(beta>0?0:1);
  else 
    {
      pnl_cdf_chn(&which,&p,&q,&xx,&nu,&lambda,&status,&bound);
      *P=(beta>0?p:q);
    }
}

/** One-Dimensional Normal Law. Cumulative distribution function.
 * Abramowitz, Milton et Stegun, Handbook of MathematicalFunctions, 1968, Dover
 * Publication, New York, page 932 (26.2.18).Precision 10-7
 *
 * @param x upper bound of the integral
 * @return the cumulative distribution function at x
 */
double pnl_cdfnor(double x)
{
  double p= 0.2316419;
  double b1= 0.319381530;
  double b2= -0.356563782;
  double b3= 1.781477937;
  double b4= -1.821255978;
  double b5= 1.330274429;
  double one_over_twopi= 0.39894228;

  double t;

  if (x >= 0.0)
    {
      t = 1.0 / ( 1.0 + p * x );
      return (1.0 - one_over_twopi * exp( -x * x / 2.0 ) * t *
              ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    } 
  else /* x < 0 */
    {
      t = 1.0 / ( 1.0 - p * x );
      return ( one_over_twopi * exp( -x * x / 2.0 ) * t *
               ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }
}

double cdf_nor(double x)
{
  return pnl_cdfnor(x);
}

/*Two-Dimensional Normal Law. Cumulative distribution function.
  Drezner-Mathematics of Computation 32(1-1978) pp.277-279*/
static double ff(double x, double y, double aa, double bb, double r)
{
  return exp(aa*(2.*x-aa) + bb*(2.*y - bb) + 2.*r*(x-aa)*(y-bb));     
}

static double NN1( double a, double b, double r)
{   
  double u[4] = {0.3253030,0.4211071,0.1334425,0.006374323};
  double v[4] = {0.1337764,0.6243247,1.3425378,2.2626645};
  double aa,bb,m;
  int i,j;

  aa = a/(sqrt(2.*(1.-r*r)));
  bb = b/(sqrt(2.*(1.-r*r)));
  m = 0.0;

  for ( i = 0; i <= 3; i++)     
    for ( j = 0; j <= 3; j++)           
      m = m + u[i]*u[j]*ff(v[i], v[j], aa, bb, r);

  m = m*(sqrt(1.0- r*r))/M_PI;
  return (m);
}

static double NN2( double a, double b, double r)
{
  double res=0.0;

  if (( a <= 0.0 ) && ( b <= 0.0 ) && ( r <= 0.0 ))
    {
      res =  NN1(a, b, r);
    }
  else if (( a <= 0.0 ) && ( b >= 0.0 ) && ( r >= 0.0 ))
    {
      res =  cdf_nor(a) - NN1(a, -b, -r);
    }
  else if (( a >= 0.0 ) && ( b <= 0.0 ) && ( r >= 0.0 ))
    {
      res = cdf_nor(b) - NN1(-a, b, -r);
    }
  else if (( a >= 0.0 ) && ( b >= 0.0 ) && ( r <= 0.0 ))
    {
      res = cdf_nor(a) + cdf_nor(b) - 1.0 + NN1(-a, -b, r);
    }
  return (res);
}

static double sgn( double x )
{
  if ( x > 0.0 ) return 1.;
  else if (x < 0.0 ) return -1.;
  return 0.;
}

/**
 * Cumulative bivariate normal distribution function
 *
 * @param a a real number
 * @param b a real number
 * @param r the correlation factor between -1 and 1
 *
 * @return  the cumlative bivariate normal distribution function with
 * correlation r at the point (a,b)
 */
double pnl_cdf2nor( double a, double b, double r)
{
  double r1,r2,d,res;

  if (r+PRECISION>=1)
    {
      res=cdf_nor(MIN(a,b));
    }
  else
    if (r-PRECISION<=-1)
      {
        if (a>-b)
          {res=cdf_nor(a)-cdf_nor(-b);}
        else
          {res=0;}
      }
    else
      if ( a*b*r <= 0.0 )
        {
          res = NN2(a,b,r);
        }
      else
        {
          r1 = (r*a - b)*sgn(a)/sqrt(a*a - 2.0*r*a*b + b*b);
          r2 = (r*b - a)*sgn(b)/sqrt(a*a - 2.0*r*a*b + b*b);
          d = (1.0- sgn(a)*sgn(b))/4.0;
          res = NN2( a, 0.0, r1 ) + NN2( b, 0.0, r2) - d;
        }

  return (res);
}

/** One-dimensional normal density function
 * @param x point at which the density is computed
 * @return the value of the density at x
 */
double pnl_normal_density(double x)
{
  return(exp(-SQR(x)/2.) * M_1_SQRT2PI);
}

double InvNormal2P1[] = {
  0.160304955844066229311E2,
  -0.90784959262960326650E2,
  0.18644914861620987391E3,
  -0.16900142734642382420E3,
  0.6545466284794487048E2,
  -0.864213011587247794E1,
  0.1760587821390590
};

double InvNormal2Q1[] = {
  0.147806470715138316110E2,
  -0.91374167024260313396E2,
  0.21015790486205317714E3,
  -0.22210254121855132366E3,
  0.10760453916055123830E3,
  -0.206010730328265443E2,
  0.1E1
};

double InvNormal2P2[] = {
  -0.152389263440726128E-1,
  0.3444556924136125216,
  -0.29344398672542478687E1,
  0.11763505705217827302E2,
  -0.22655292823101104193E2,
  0.19121334396580330163E2,
  -0.5478927619598318769E1,
  0.237516689024448000
};

double InvNormal2Q2[] = {
  -0.108465169602059954E-1,
  0.2610628885843078511,
  -0.24068318104393757995E1,
  0.10695129973387014469E2,
  -0.23716715521596581025E2,
  0.24640158943917284883E2,
  -0.10014376349783070835E2,
  0.1E1
};

double InvNormal2P3[] = {
  0.56451977709864482298E-4,
  0.53504147487893013765E-2,
  0.12969550099727352403,
  0.10426158549298266122E1,
  0.28302677901754489974E1,
  0.26255672879448072726E1,
  0.20789742630174917228E1,
  0.72718806231556811306,
  0.66816807711804989575E-1,
  -0.17791004575111759979E-1,
  0.22419563223346345828E-2
};

double InvNormal2Q3[] = {
  0.56451699862760651514E-4,
  0.53505587067930653953E-2,
  0.12986615416911646934,
  0.10542932232626491195E1,
  0.30379331173522206237E1,
  0.37631168536405028901E1,
  0.38782858277042011263E1,
  0.20372431817412177929E1,
  0.1E1
};

/**
 * Computes the inverse of the normal cumulative distribution function.
 * Comes from Pierre L'Ecuyer (SSJ)
 */
double pnl_inv_cdfnor (double u)
{
  /*
   * Returns the inverse of the cdf of the normal distribution.
   * Rational approximations giving 16 decimals of precision.
   * The authors also give an approximation with 23 decimals of
   * precision.
   * J.M. Blair, C.A. Edwards, J.H. Johnson, "Rational Chebyshev
   * approximations for the Inverse of the Error Function", in
   * Mathematics of Computation, Vol. 30, 136, pp 827, (1976)
   */

  int i;
  int negatif;
  double z, v, numer, denom;

  if (u==0.) return log(0.);
  if (u==1.) return -log(0.);
  if (u <= 0.0 || u >= 1.0)
    {
      printf ("u is out of range\n"); abort();
    }

  /* Transform u as argument of InvErf */
  z = u;
  u = 2.0*u - 1.0;
  if (u >= 1.0)
    return 1000.0;

  if (u < 0.0) {
    u = -u;
    negatif = 1;
  } else
    negatif = 0;

  if (u <= 0.75) {
    v = u*u - 0.5625;
    numer = denom = 0.0;

    /* Evaluation of the 2 polynomials by Horner  */
    for (i = 6; i >= 0; i--) {
      numer = numer*v + InvNormal2P1[i];
      denom = denom*v + InvNormal2Q1[i];
    }

    z = u*numer/denom;
  }
  else if (u <= 0.9375) {
    v = u*u - 0.87890625;
    numer = denom = 0.0;

    for ( i = 7; i >= 0; i-- ) {
      numer = numer*v + InvNormal2P2[i];
      denom = denom*v + InvNormal2Q2[i];
    }

    z = u*numer/denom;
  }
  else {
    if (z > 0.1)
      v = 1.0 / sqrt(-log (1.0 - u));
    else
      v = 1.0 / sqrt (-log (2.0*z));         
    numer = denom = 0.0;

    for (i = 10; i >= 0; i--)
      numer = numer*v + InvNormal2P3[i];

    for (i = 8; i >= 0; i--)
      denom = denom*v + InvNormal2Q3[i];

    z = (1.0/v)*numer/denom;
  }

  if (negatif)
    return -z*M_SQRT2;
  else
    return z*M_SQRT2;
}

