#ifndef __COMPLEX_H__
#define __COMPLEX_H__

/**
 * \defgroup PnlComplex  Complex structure
 * \brief Operations on complex numbers.
 */
/*@{*/

typedef struct {
    double r; /*!< real part */
    double i; /*!< imaginary part */
} fcomplex;

/** zero complex  0 + i 0  */
#define CZERO (Complex(0,0))
/** unitary real complex  1 + i 0 */
#define CUNO (Complex(1,0))
#define CONE (Complex(1,0))
/** unitary pure imaginary complex  0 + i */
#define CI (Complex(0,1)) 
#define CMPLX(z) z.r, z.i

extern double Creal( fcomplex g );
extern double Cimag( fcomplex g );
extern fcomplex Cadd(fcomplex a, fcomplex b);
extern fcomplex CRadd(fcomplex z, double x);
extern fcomplex RCadd(double b, fcomplex z);
extern fcomplex Csub(fcomplex a, fcomplex b);
extern fcomplex CRsub(fcomplex a, double b);
extern fcomplex RCsub(double a, fcomplex b);
extern fcomplex Cminus (fcomplex z);
extern fcomplex Ciadd(fcomplex a, fcomplex b);
extern fcomplex Cisub(fcomplex a, fcomplex b);
extern fcomplex Cmul(fcomplex a, fcomplex b);
extern fcomplex RCmul(double a, fcomplex b);
extern fcomplex CRmul(fcomplex b, double a);
extern fcomplex Complex(double re, double im);
extern fcomplex Complex_polar(double r, double theta);
extern void Cprintf( fcomplex z);
extern fcomplex Conj(fcomplex z);
extern fcomplex Cinv(fcomplex a);
extern fcomplex Cdiv(fcomplex a, fcomplex b);
extern fcomplex RCdiv(double a, fcomplex b);
extern fcomplex CRdiv(fcomplex a, double b);
extern double Csqr_norm(fcomplex z);
extern double Cabs(fcomplex z);
extern fcomplex Csqrt(fcomplex z);
extern fcomplex Clog(fcomplex z);
extern fcomplex Cexp(fcomplex z);
extern fcomplex CIexp(double t);
extern double Carg(fcomplex a);
extern fcomplex Ctgamma(fcomplex a);
extern fcomplex Clgamma(fcomplex xx);
extern fcomplex Ccos(fcomplex g);
extern fcomplex Csin(fcomplex g);
extern fcomplex Ctan(fcomplex z);
extern fcomplex Ccotan(fcomplex z);
extern fcomplex Ccosh(fcomplex g);
extern fcomplex Csinh(fcomplex g);
extern fcomplex Ctanh(fcomplex z);
extern fcomplex Ccotanh(fcomplex z);
extern fcomplex Cpow(fcomplex z, fcomplex exp);
extern fcomplex Cpow_real (fcomplex z, double y);
/* Algebirc operation on C : */
/* use i for multiply by i */
/* use c for congugate */
/* use d for double */
/* use a,b for complex */
/* use p or m for plus or minus */
extern fcomplex C_op_apib(fcomplex a, fcomplex b);
extern fcomplex C_op_amib(fcomplex a, fcomplex b);
extern fcomplex C_op_apcb(fcomplex a, fcomplex b);
extern fcomplex C_op_amcb(fcomplex a, fcomplex b);
extern fcomplex C_op_dapb(double d,fcomplex a, fcomplex b);
extern fcomplex C_op_damb(double d,fcomplex a, fcomplex b);
extern fcomplex C_op_dapib(double d,fcomplex a, fcomplex b);
extern fcomplex C_op_damib(double d,fcomplex a, fcomplex b);
extern fcomplex C_op_dapcb(double d,fcomplex a, fcomplex b);
extern fcomplex C_op_damcb(double d,fcomplex a, fcomplex b);
extern fcomplex C_op_idapb(double d,fcomplex a, fcomplex b);
extern fcomplex C_op_idamb(double d,fcomplex a, fcomplex b);
extern fcomplex C_op_idapcb(double d,fcomplex a, fcomplex b);
extern fcomplex C_op_idamcb(double d,fcomplex a, fcomplex b);

/*@}*/

typedef struct {
  fcomplex (*function) (fcomplex x, void *params);
  void *params;
} PnlCmplxFunc;

#endif /* __COMPLEX_H_ */



