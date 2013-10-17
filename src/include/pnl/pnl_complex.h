#ifndef _PNL_COMPLEX_H
#define _PNL_COMPLEX_H


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**
 * \defgroup PnlComplex  Complex numbers
 * \brief Operations on complex numbers.
 */
/*@{*/

typedef struct
{
    double r; /*!< real part */
    double i; /*!< imaginary part */
} dcomplex;

/** zero complex  0 + i 0  */
#define CZERO (Complex(0.,0.))
/** unitary real complex  1 + i 0 */
#define CUNO (Complex(1.,0.))
#define CONE (Complex(1.,0.))
/** unitary pure imaginary complex  0 + i */
#define CI (Complex(0.,1.)) 
#define CMPLX(z) z.r, z.i

extern double Creal( dcomplex g );
extern double Cimag( dcomplex g );
extern dcomplex Cadd(dcomplex a, dcomplex b);
extern dcomplex CRadd(dcomplex z, double x);
extern dcomplex RCadd(double b, dcomplex z);
extern dcomplex Csub(dcomplex a, dcomplex b);
extern dcomplex CRsub(dcomplex a, double b);
extern dcomplex RCsub(double a, dcomplex b);
extern dcomplex Cminus (dcomplex z);
extern dcomplex Cmul(dcomplex a, dcomplex b);
extern dcomplex RCmul(double a, dcomplex b);
extern dcomplex CRmul(dcomplex b, double a);
extern dcomplex Complex(double re, double im);
extern dcomplex Complex_polar(double r, double theta);
extern void Cprintf( dcomplex z);
extern dcomplex Conj(dcomplex z);
extern dcomplex Cinv(dcomplex a);
extern dcomplex Cdiv(dcomplex a, dcomplex b);
extern dcomplex RCdiv(double a, dcomplex b);
extern dcomplex CRdiv(dcomplex a, double b);
extern double Csqr_norm(dcomplex z);
extern double Cabs(dcomplex z);
extern dcomplex Csqrt(dcomplex z);
extern dcomplex Clog(dcomplex z);
extern dcomplex Cexp(dcomplex z);
extern dcomplex CIexp(double t);
extern double Carg(dcomplex a);
extern dcomplex Ctgamma(dcomplex a);
extern dcomplex Clgamma(dcomplex xx);
extern dcomplex Ccos(dcomplex g);
extern dcomplex Csin(dcomplex g);
extern dcomplex Ctan(dcomplex z);
extern dcomplex Ccotan(dcomplex z);
extern dcomplex Ccosh(dcomplex g);
extern dcomplex Csinh(dcomplex g);
extern dcomplex Ctanh(dcomplex z);
extern dcomplex Ccotanh(dcomplex z);
extern dcomplex Cpow(dcomplex z, dcomplex exp);
extern dcomplex Cpow_real (dcomplex z, double y);
/* Algebirc operation on C : */
/* use i for multiply by i */
/* use c for congugate */
/* use d for double */
/* use a,b for complex */
/* use p or m for plus or minus */
extern dcomplex C_op_apib(dcomplex a, dcomplex b);
extern dcomplex C_op_amib(dcomplex a, dcomplex b);
extern dcomplex C_op_apcb(dcomplex a, dcomplex b);
extern dcomplex C_op_amcb(dcomplex a, dcomplex b);
extern dcomplex C_op_dapb(double d,dcomplex a, dcomplex b);
extern dcomplex C_op_damb(double d,dcomplex a, dcomplex b);
extern dcomplex C_op_dapib(double d,dcomplex a, dcomplex b);
extern dcomplex C_op_damib(double d,dcomplex a, dcomplex b);
extern dcomplex C_op_dapcb(double d,dcomplex a, dcomplex b);
extern dcomplex C_op_damcb(double d,dcomplex a, dcomplex b);
extern dcomplex C_op_idapb(double d,dcomplex a, dcomplex b);
extern dcomplex C_op_idamb(double d,dcomplex a, dcomplex b);
extern dcomplex C_op_idapcb(double d,dcomplex a, dcomplex b);
extern dcomplex C_op_idamcb(double d,dcomplex a, dcomplex b);

/*@}*/

typedef struct
{
  dcomplex (*F) (dcomplex x, void *params);
  void *params;
} PnlCmplxFunc;

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_COMPLEX_H */
