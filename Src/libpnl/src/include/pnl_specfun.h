#ifndef  _SPECIAL_FUNC_H
#define _SPECIAL_FUNS_H

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

void pnl_gamma_inc(double a,double x,double * Ga,double *P,double *Q);
double pnl_gamma_inc_func(double a,double x);
double pnl_expint_En(int n,double x);
double pnl_expint_Ei(double x);

#endif
