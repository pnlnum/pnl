#include <gsl/gsl_specfunc.h>
#include <gsl/gsl_sf_hyperg.h>
#include <stdlib.h>
#include <time.h>




int main ()
{
  int      n;
  FILE *dat;
  srand48(time(NULL));

  dat = fopen ("hyperg0F1_test.dat", "w");
  for(n=1;n<10;n++)
    {
      double x = drand48() * 20 - 10.;
      double c = drand48() * 20 - 10.;
      double res = gsl_sf_hyperg_0F1(c,x);
      fprintf(dat,"{ \"%s\", %s, %.18f, %.18f, %.18f },\n", 
              "pnl_sf_hyperg_0F1", "pnl_sf_hyperg_0F1", x, c, res);
    }
  fclose (dat);

  dat = fopen ("hyperg1F1_test.dat", "w");
  for(n=1;n<10;n++)
    {
      double x = drand48() * 20 - 10.;
      double b = drand48() * 20 - 10.;
      double a = drand48() * 20 - 10.;
      double res = gsl_sf_hyperg_1F1(a, b, x);
      fprintf(dat,"{ \"%s\", %s, %.18f, %.18f, %.18f, %.18f },\n", 
              "pnl_sf_hyperg_1F1", "pnl_sf_hyperg_1F1", a, b, x, res);
    }
  fclose (dat);

  dat = fopen ("hypergU_test.dat", "w");
  for(n=1;n<10;n++)
    {
      double x = drand48() * 10; /* x > 0 */
      double b = drand48() * 20 - 10.;
      double a = drand48() * 20 - 10.;
      double res = gsl_sf_hyperg_U(a, b, x);
      fprintf(dat,"{ \"%s\", %s, %.18f, %.18f, %.18f, %.18f },\n", 
              "pnl_sf_hyperg_U", "pnl_sf_hyperg_U", a, b, x, res);
    }
  fclose (dat);

  dat = fopen ("hyperg2F0_test.dat", "w");
  for(n=1;n<10;n++)
    {
      double x = -drand48() * 10; /* x < 0 */
      double b = drand48() * 20 - 10.;
      double a = drand48() * 20 - 10.;
      double res = gsl_sf_hyperg_2F0(a, b, x);
      fprintf(dat,"{ \"%s\", %s, %.18f, %.18f, %.18f, %.18f },\n", 
              "pnl_sf_hyperg_2F0", "pnl_sf_hyperg_2F0", a, b, x, res);
    }
  fclose (dat);

  dat = fopen ("hyperg2F1_test.dat", "w");
  for(n=1;n<10;n++)
    {
      double x = 2 * drand48() - 1; /* |x| < 1 */
      double c = drand48() * 20 - 10.;
      double b = drand48() * 20 - 10.;
      double a = drand48() * 20 - 10.;
      double res = gsl_sf_hyperg_2F1(a, b, c, x);
      fprintf(dat,"{ \"%s\", %s, %.18f, %.18f, %.18f, %.18f, %.18f },\n", 
              "pnl_sf_hyperg_2F1", "pnl_sf_hyperg_2F1", a, b, c, x, res);
    }
  fclose (dat);

  exit (0);
}

