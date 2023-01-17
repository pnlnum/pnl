#include "pnl/pnl_finance.h"

int main()
{
  int error;
  double vol;
  vol = pnl_bs_implicit_vol(0,8.040982,100,100,0.5,0.09531,0,&error);
  printf("err: %d, vol: %f\n", error, vol);
  exit(0);
}
