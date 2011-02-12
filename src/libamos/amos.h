#ifndef AMOS_H
#define AMOS_H

#include "pnl/pnl_complex.h"
#include "pnl/pnl_mathtools.h"

#define D_SIGN(a,b) ( b >= 0 ? (a >= 0 ? a : - a) : -(a >= 0 ? a : -a))
#define D_INT(x) ( (x>0) ? floor(x) : -floor(- x) )

extern double pnl_d1mach (int i);
extern double pnl_dlamch (char *cmach);
extern int pnl_ipmpar (int i);
extern double amos_dgamln (double *z__, int *ierr);
extern int amos_dsclmr (void);
extern int amos_fdump (void);
extern int amos_i1mach (int i);
extern int amos_xerror (char *mess, int *nmess, int *l1, int *l2,
			long int mess_len);
extern double amos_azabs (const double *zr,const double *zi);
extern int amos_zacai (double *zr, double *zi, double *fnu,const int *kode,
		       int *mr,const int *n, double *yr, double *yi, int *nz,
		       double *rl, double *tol, double *elim, double *alim);
extern int amos_zacon (double *zr, double *zi, double *fnu, int *kode,
		       int *mr, const int *n, double *yr, double *yi, int *nz,
		       double *rl, double *fnul, double *tol, double *elim,
		       double *alim);
extern int pnl_zairy (double *zr, double *zi,const int *id,const int *kode,
		       double *air, double *aii, int *nz, int *ierr);
extern int amos_zasyi (double *zr, double *zi, double *fnu,const int *kode,const  int *n,
		       double *yr, double *yi, int *nz, double *rl,
		       double *tol, double *elim, double *alim);
extern int pnl_zbesh (double *zr, double *zi, double *fnu, int *kode,const int *m,
		       const int *n, double *cyr, double *cyi, int *nz, int *ierr);
extern int pnl_zbesi (double *zr, double *zi, double *fnu, int *kode,const  int *n,
		       double *cyr, double *cyi, int *nz, int *ierr);
extern int pnl_zbesj (double *zr, double *zi, double *fnu, int *kode,const  int *n,
		       double *cyr, double *cyi, int *nz, int *ierr);
extern int pnl_zbesk (double *zr, double *zi, double *fnu, int *kode,const  int *n,
		       double *cyr, double *cyi, int *nz, int *ierr);
extern int pnl_zbesy (double *zr, double *zi, double *fnu, int *kode,const  int *n,
		       double *cyr, double *cyi, int *nz, double *cwrkr,
		       double *cwrki, int *ierr);
extern int amos_zbinu (double *zr, double *zi, double *fnu, int *kode,const  int *n,
		       double *cyr, double *cyi, int *nz, double *rl,
		       double *fnul, double *tol, double *elim, double *alim);
extern int pnl_zbiry  (double *zr, double *zi, int *id, int *kode,
		       double *bir, double *bii, int *ierr);
extern int amos_zbknu (double *zr, double *zi, double *fnu,const int *kode,const int *n,
		       double *yr, double *yi, int *nz, double *tol,
		       double *elim, double *alim);
extern int amos_zbuni (double *zr, double *zi, double *fnu, int *kode,const  int *n,
		       double *yr, double *yi, int *nz, int *nui, int *nlast,
		       double *fnul, double *tol, double *elim, double *alim);
extern int amos_zbunk (double *zr, double *zi, double *fnu, int *kode,
		       int *mr,const  int *n, double *yr, double *yi, int *nz,
		       double *tol, double *elim, double *alim);

extern int amos_zdiv (const double *ar,const  double *ai, const double *br, 
		      const double *bi, double *cr, double *ci);
extern int amos_azexp (const double *ar,const double *ai, double *br, double *bi);
extern int amos_zkscl (double *zrr, double *zri, double *fnu,const int *n,
		       double *yr, double *yi, int *nz, double *rzr,
		       double *rzi, double *ascle, double *tol, double *elim);
extern int amos_azlog (const double *ar,const double *ai, double *br, double *bi, int *ierr);

extern int amos_zmlri (double *zr, double *zi, double *fnu, const int *kode,const  int *n,
		       double *yr, double *yi, int *nz, double *tol);
extern int amos_zmlt (double *ar, double *ai, double *br, double *bi,
		      double *cr, double *ci);
extern int pnl_zrati (double *zr, double *zi, double *fnu,const  int *n,
		       double *cyr, double *cyi, double *tol);
extern int amos_zs1s2 (double *zrr, double *zri, double *s1r, double *s1i,
		       double *s2r, double *s2i, int *nz, double *ascle,
		       double *alim, int *iuf);
extern int amos_zseri (double *zr, double *zi, double *fnu,const int *kode,const  int *n,
		       double *yr, double *yi, int *nz, double *tol,
		       double *elim, double *alim);
extern int amos_zshch (double *zr, double *zi, double *cshr, double *cshi,
		       double *cchr, double *cchi);
extern int amos_azsqrt (const double *ar,const double *ai, double *br, double *bi);

extern int amos_zuchk (double *yr, double *yi, int *nz, double *ascle,
		       double *tol);
extern int amos_zunhj (double *zr, double *zi, double *fnu,const int *ipmtr,
		       double *tol, double *phir, double *phii, double *argr,
		       double *argi, double *zeta1r, double *zeta1i,
		       double *zeta2r, double *zeta2i, double *asumr,
		       double *asumi, double *bsumr, double *bsumi);
extern int amos_zuni1 (double *zr, double *zi, double *fnu, int *kode,const int *n,
		       double *yr, double *yi, int *nz, int *nlast,
		       double *fnul, double *tol, double *elim, double *alim);
extern int pnl_zuni2 (double *zr, double *zi, double *fnu, int *kode,const int *n,
		       double *yr, double *yi, int *nz, int *nlast,
		       double *fnul, double *tol, double *elim, double *alim);
extern int amos_zunik (double *zrr, double *zri, double *fnu,const int *ikflg,const int *ipmtr,
		       double *tol, int *init, double *phir, double *phii,
		       double *zeta1r, double *zeta1i, double *zeta2r, double *zeta2i,
		       double *sumr, double *sumi, double *cwrkr, double *cwrki);

extern int amos_zunk1 (double *zr, double *zi, double *fnu, int *kode,
		       int *mr,const  int *n, double *yr, double *yi, int *nz,
		       double *tol, double *elim, double *alim);
extern int amos_zunk2 (double *zr, double *zi, double *fnu, int *kode,
		       int *mr,const  int *n, double *yr, double *yi, int *nz,
		       double *tol, double *elim, double *alim);
extern int amos_zuoik (double *zr, double *zi, double *fnu,const int *kode,
		       const int *ikflg,const int *n, double *yr, double *yi, int *nuf,
		       double *tol, double *elim, double *alim);
extern int amos_zwrsk (double *zrr, double *zri, double *fnu, int *kode,
		       const int *n, double *yr, double *yi, int *nz, double *cwr,
		       double *cwi, double *tol, double *elim, double *alim);


#endif


