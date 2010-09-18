#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <stdarg.h>
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_complex.h"

void pnl_array_map_inplace(double *lhs, double(*f)(double ),int n);
void pnl_array_plus_double(double *lhs , double x, int n);
void pnl_array_minus_double(double *lhs , double x, int n);
void pnl_array_mult_double(double *lhs , double x, int n);
void pnl_array_div_double(double *lhs , double x, int n);
void pnl_array_map(double *lhs, const double *rhs, double(*f)(double),int n);
void pnl_array_map_array(double *lhs, const double *rhs, double(*f)(double,double),int n);
void pnl_array_inv_term(double *lhs,int n);
void pnl_array_plus_array_term(double *lhs, const double *rhs,int n);
void pnl_array_minus_array_term(double *lhs, const double *rhs,int n);
void pnl_array_mult_array_term(double *lhs, const double *rhs,int n);
void pnl_array_div_array_term(double *lhs, const double *rhs,int n);



void pnl_array_complex_map_inplace(dcomplex *lhs,dcomplex(*f)(dcomplex ),int n);
void pnl_array_complex_plus_dcomplex(dcomplex *lhs , dcomplex x, int n);
void pnl_array_complex_minus_dcomplex(dcomplex *lhs , dcomplex x, int n);
void pnl_array_complex_mult_dcomplex(dcomplex *lhs , dcomplex x, int n);
void pnl_array_complex_div_dcomplex(dcomplex *lhs , dcomplex x, int n);
void pnl_array_complex_map(dcomplex *lhs, const dcomplex *rhs, dcomplex(*f)(dcomplex),int n);
void pnl_array_complex_map_array(dcomplex *lhs, const dcomplex *rhs, dcomplex(*f)(dcomplex,dcomplex),int n);
void pnl_array_complex_plus_array_term(dcomplex *lhs, const dcomplex *rhs,int n);
void pnl_array_complex_minus_array_term(dcomplex *lhs, const dcomplex *rhs,int n);
void pnl_array_complex_mult_array_term(dcomplex *lhs, const dcomplex *rhs,int n);
void pnl_array_complex_div_array_term(dcomplex *lhs, const dcomplex *rhs,int n);
void pnl_array_complex_inv_term(dcomplex *lhs,int n);



void pnl_array_int_map_inplace(int *lhs,int(*f)(int ),int n);
void pnl_array_int_plus_int(int *lhs , int x, int n);
void pnl_array_int_minus_int(int *lhs , int x, int n);
void pnl_array_int_mult_int(int *lhs , int x, int n);
void pnl_array_int_div_int(int *lhs , int x, int n);
void pnl_array_int_map(int *lhs, const int *rhs, int(*f)(int),int n);
void pnl_array_int_map_array(int *lhs, const int *rhs, int(*f)(int,int),int n);
void pnl_array_int_plus_array_term(int *lhs, const int *rhs,int n);
void pnl_array_int_minus_array_term(int *lhs, const int *rhs,int n);
void pnl_array_int_mult_array_term(int *lhs, const int *rhs,int n);
void pnl_array_int_inv_term(int *lhs,int n);
void pnl_array_int_div_array_term(int *lhs, const int *rhs,int n);



void pnl_array_uint_map_inplace(uint *lhs,uint(*f)(uint ),int n);
void pnl_array_uint_plus_uint(uint *lhs , uint x, int n);
void pnl_array_uint_minus_uint(uint *lhs , uint x, int n);
void pnl_array_uint_mult_uint(uint *lhs , uint x, int n);
void pnl_array_uint_div_uint(uint *lhs , uint x, int n);
void pnl_array_uint_map(uint *lhs, const uint *rhs, uint(*f)(uint),int n);
void pnl_array_uint_map_array(uint *lhs, const uint *rhs, uint(*f)(uint,uint),int n);
void pnl_array_uint_plus_array_term(uint *lhs, const uint *rhs,int n);
void pnl_array_uint_mult_array_term(uint *lhs, const uint *rhs,int n);
void pnl_array_uint_minus_array_term(uint *lhs, const uint *rhs,int n);
void pnl_array_uint_div_array_term(uint *lhs, const uint *rhs,int n);




