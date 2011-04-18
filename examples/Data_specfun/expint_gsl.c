#include <gsl/gsl_specfunc.h>
#include <stdlib.h>
#include <time.h>




int main ()
{
  int      n;
  FILE *dat;
  srand48(time(NULL));

  dat = fopen ("expint_test.dat", "w");

  for(n=1;n<10;n++)
    {
      double x = drand48() * 20;
      double res = gsl_sf_expint_En (n,x);
      fprintf(dat,"{ \"%s\", %s, %d, %.18f, %.18f },\n", 
              "pnl_sf_expint_En", "pnl_sf_expint_En", n, x, res);
    }
  fclose (dat);

  dat = fopen ("expintEi_test.dat", "w");

  for(n=1;n<10;n++)
    {
      double x = drand48() * 20;
      double res = gsl_sf_expint_Ei (x);
      fprintf(dat,"{ \"%s\", %s, %.18f, %.18f },\n", 
              "pnl_sf_expint_Ei", "pnl_sf_expint_Ei", x, res);
    }
  fclose (dat);


  exit (0);
}

