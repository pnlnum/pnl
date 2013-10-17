#include "pnl/pnl_config.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_cdf.h"
#include "pnl/pnl_finance.h"
#include "pnl/pnl_root.h"


/************************************************************************/
/* Copyright David Pommier <pommier.david@gmail.com>                    */
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

/*
 * This code has been significantly modified by 
 * Jérôme Lelong <jerome.lelong@gmail.com>
 * November 2009
 */

typedef struct Pnl_Data_Vol_Impli_BS{
  int is_call;
  double price, r, divid, spot, strike, T;
}Pnl_Data_Vol_Impli_BS;


/**
 * give the price and delta of a call option in BS model
 *
 * @param s a double value of spot
 * @param k a double, the Strike
 * @param t a double, the Maturity
 * @param r a double, the interest rate
 * @param divid a double, dividend rate
 * @param sigma a doule the volatility
 * @param ptprice a pointer to store price
 * @param ptdelta a pointer to store first derivative against spot
 */
int pnl_cf_call_bs(double s, double k, double t, double r, double divid, 
                   double sigma, double *ptprice, double *ptdelta)
{
  double V_Sqrt_T;
  double D1;
  double D2;

  PNL_CHECK(sigma < 0., "Volatility required to be >= 0", "pnl_cf_call_bs");
  PNL_CHECK(t < 0., "Maturity required to be >= 0", "pnl_cf_call_bs");
  PNL_CHECK(r < 0., "Maturity required to be >= 0", "pnl_cf_call_bs");
  PNL_CHECK(k < 0., "Strike required to be >= 0", "pnl_cf_call_bs");

  V_Sqrt_T = sigma * sqrt(t);
  D1 = (log (s / k) + ((r - divid) *t + 0.5 * sigma * sigma * t)) / V_Sqrt_T; 
  D2 = D1 - V_Sqrt_T;
  /* price */
  *ptprice = s * exp (-divid * t) * cdf_nor(D1) - k * exp (-r * t) * cdf_nor(D2);
  /*Delta*/
  *ptdelta = exp ( -divid * t) * cdf_nor(D1);
  return OK;
}

/**
 * give the price and delta of a put option in BS model
 *
 * @param s a double value of spot
 * @param k a double, the Strike
 * @param t a double, the Maturity
 * @param r a double, the interest rate
 * @param divid a double, dividend rate
 * @param sigma a doule the volatility
 * @param ptprice a pointer to store price
 * @param ptdelta a pointer to store first derivative against spot
 */
int pnl_cf_put_bs(double s, double k, double t, double r, double divid,
                  double sigma, double *ptprice, double *ptdelta)
{
  double V_Sqrt_T;
  double D1;
  double D2;

  PNL_CHECK(sigma < 0., "Volatility required to be >= 0", "pnl_cf_put_bs");
  PNL_CHECK(t < 0., "Maturity required to be >= 0", "pnl_cf_put_bs");
  PNL_CHECK(r < 0., "Maturity required to be >= 0", "pnl_cf_put_bs");
  PNL_CHECK(k < 0., "Strike required to be >= 0", "pnl_cf_put_bs");

  V_Sqrt_T = sigma * sqrt(t);
  D1 = (log (s / k) + ((r - divid) *t + 0.5 * sigma * sigma * t)) / V_Sqrt_T; 
  D2 = D1 - V_Sqrt_T;
  /* price */
  *ptprice = k * exp (-r * t) * cdf_nor(-D2) - s * exp (-divid * t) * cdf_nor(-D1);
  /*Delta*/
  *ptdelta = -exp ( -divid * t) * cdf_nor(-D1);
  return OK;
}

/**
 * give the price of a call option in BS model
 *
 * @param s a double value of spot
 * @param k a double, the Strike
 * @param t a double, the Maturity
 * @param r a double, the interest rate
 * @param divid a double, dividend rate
 * @param sigma a doule the volatility
 */
double pnl_bs_call(double s, double k, double t, double r, double divid, double sigma)
{
  double ptprice, ptdelta;

  pnl_cf_call_bs(s, k, t, r, divid, sigma, &ptprice, &ptdelta);
  return ptprice;
}

/**
 * give the price of a put option in BS model
 *
 * @param s a double value of spot
 * @param k a double, the Strike
 * @param t a double, the Maturity
 * @param r a double, the interest rate
 * @param divid a double, dividend rate
 * @param sigma a doule the volatility
 */
double pnl_bs_put(double s, double k, double t, double r, double divid, double sigma)
{
  double ptprice, ptdelta;

  pnl_cf_put_bs(s, k, t, r, divid, sigma, &ptprice, &ptdelta);
  return ptprice;
}

/**
 * give the price of a call/put option in BS model
 *
 * @param iscall an integer (1 for a call, 0 for a put)
 * @param s a double value of spot
 * @param k a double, the Strike
 * @param t a double, the Maturity
 * @param r a double, the interest rate
 * @param divid a double, dividend rate
 * @param sigma a doule the volatility
 */
double pnl_bs_call_put(int iscall, double s, double k, double t, double r, double divid, double sigma)
{
  if (iscall) 
    return pnl_bs_call(s, k, t, r, divid, sigma);
  else
    return pnl_bs_put(s, k, t, r, divid, sigma);
}


/**
 * give the second derivative w.r.t the spot of a call/put option price in the BS model 
 *
 * @param s a double value of spot
 * @param k a double, the Strike
 * @param t a double, the Maturity
 * @param r a double, the interest rate
 * @param divid a double, dividend rate
 * @param sigma a doule the volatility
 * @return gamma of a call/put option
 */
double pnl_bs_gamma(double s, double k, double t, double r, double divid, double sigma)
{
  double V_Sqrt_T;
  double D1;

  PNL_CHECK(sigma < 0., "Volatility required to be >= 0", "pnl_cf_gamma_bs");
  PNL_CHECK(t < 0., "Maturity required to be >= 0", "pnl_cf_gamma_bs");
  PNL_CHECK(r < 0., "Maturity required to be >= 0", "pnl_cf_gamma_bs");
  PNL_CHECK(k < 0., "Strike required to be >= 0", "pnl_cf_gamma_bs");

  V_Sqrt_T = sigma * sqrt(t);
  D1 = (log (s / k) + ((r - divid) *t + 0.5 * sigma * sigma * t)) / V_Sqrt_T; 
  return exp (-divid * t) * pnl_normal_density (D1) / (s * V_Sqrt_T);
}

/**
 * give the first derivative of the price w.r.t the volatility of a call/put
 * option in a BS model
 *
 * @param s a double value of spot
 * @param k a double, the Strike
 * @param t a double, the Maturity
 * @param r a double, the interest rate
 * @param divid a double, dividend rate
 * @param sigma a doule the volatility
 * @return vega of a call/put option
 */
double pnl_bs_vega (double s, double k, double t, double r, double divid, double sigma)
{
  double V_Sqrt_T;
  double D1;

  PNL_CHECK(sigma < 0., "Volatility required to be >= 0", "pnl_cf_vega_bs");
  PNL_CHECK(t < 0., "Maturity required to be >= 0", "pnl_cf_vega_bs");
  PNL_CHECK(r < 0., "Maturity required to be >= 0", "pnl_cf_vega_bs");
  PNL_CHECK(k < 0., "Strike required to be >= 0", "pnl_cf_vega_bs");

  V_Sqrt_T = sigma * sqrt(t);
  D1 = (log (s / k) + ((r - divid) *t + 0.5 * sigma * sigma * t)) / V_Sqrt_T; 
  return s * exp (-divid * t) * pnl_normal_density (D1) * sqrt (t);

}

static void pnl_bs_increment_call_put_Type(double x, double * fx, double * dfx, const Pnl_Data_Vol_Impli_BS * Data)
{
  *fx = Data->price - pnl_bs_call_put (Data->is_call, Data->spot, Data->strike, Data->T, Data->r, Data->divid, x);
  *dfx = - pnl_bs_vega (Data->spot, Data->strike, Data->T, Data->r, Data->divid, x);
}

static void pnl_bs_increment_call_put(double x, double * fx, double * dfx, void* Data)
{ pnl_bs_increment_call_put_Type(x,fx,dfx,Data);}

/**
 * Computes the implied volatility of an option price
 *
 * @param is_call a int, the option type, 1 for call, 0 for put
 * @param Price a double, option price today
 * @param r the instantaneous interest rate
 * @param divid the instantaneous dividend rate
 * @param spot the initial value of the asset
 * @param Strike a double, for value contract
 * @param T a double, echeance time (T)
 * @param error an integer containing the error code on output (OK or FAIL) 
 * @return implied of a call/put option
 */
double pnl_bs_implicit_vol (int is_call, double Price, double spot, double Strike,
                            double T, double r, double divid, int *error)
{
  double impli_vol;
  Pnl_Data_Vol_Impli_BS data;
  PnlFuncDFunc func;

  data.is_call = is_call;
  data.price = Price;
  data.spot = spot;
  data.r = r;
  data.divid = divid;
  data.strike = Strike;
  data.T = T;
  *error = FAIL;
  if(is_call)
    {
      if (Price <= exp(-r * T) * MAX (spot * exp ( (r-divid) * T) - Strike, 0.0))
        return 0.0;
      else if (Price >= spot * exp (-divid * T))
        return PNL_POSINF;
    }
  else
    {
      if (Price <= exp (-r * T) * MAX (Strike - spot * exp ( (r-divid) * T), 0.0))
        return 0.0;
      else if (Price >= exp (-r * T) * Strike)
        return PNL_POSINF;
    }

  func.F = pnl_bs_increment_call_put;
  func.params = &data;
  *error = pnl_root_newton_bisection(&func,0.001,10.0,0.0001,20,&impli_vol);
  return impli_vol;
}


/**
 * compute implied volatility matrix of a list of options prices
 *
 * @param is_call a matrix of int (1 for call, 0 for put)
 * @param Price a matrix of Prices
 * @param spot the spot price
 * @param r the instantaneous interest rate
 * @param divid the instantaneous dividend rate
 * @param Strike a Vector, for value contract
 * @param Maturity a Vector, echeance time (T)
 * @param Vol a Matrix, to store matrix volatility
 * @return error code  indicating the number of errors (0 if all was OK)
 */
int pnl_bs_matrix_implicit_vol (const PnlMatInt * is_call, const PnlMat * Price,double spot,double r, double divid,
                                const PnlVect * Strike,const PnlVect * Maturity,PnlMat * Vol)
{
  double T, price, k;
  int i, j ,iscall;
  int error, g_error = 0;

  for(j=0;j<Maturity->size;j++)
    {
      T = GET(Maturity,j);
      for(i=0;i<Strike->size;i++)
        {
          price = MGET(Price,i,j);;
          iscall = pnl_mat_int_get(is_call,i,j);
          k = GET(Strike,i);

          MLET(Vol, i, j) = pnl_bs_implicit_vol (iscall, price, spot, k, T, r, divid, &error);
          g_error += (error == FAIL ? 1 : 0);
        }
    }
  return g_error;
}

