
/************************************************************************/
/* Copyright Jérôme Lelong <jerome.lelong@gmail.com>                    */
/* Copyright David Pommier <david.pommier@gmail.com>                    */
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_specfun.h"
#include "pnl/pnl_random.h"
#include "tests_utils.h"

static double wrap_pnl_sf_choose (double n, double p)
{
  return pnl_sf_choose ((int) n, (int) p);
}

/* 
 * struct for double(*f)(double)
 */
struct d2d_test
{
  char *label;
  double (*f)(double);
  double arg;
  double res;
};

/* 
 * struct for double(*f)(double, double)
 */
struct dd2d_test
{
  char *label;
  double (*f)(double, double);
  double arg;
  double nu;
  double res;
};

/* 
 * struct for double(*f)(int, double)
 */
struct id2d_test
{
  char *label;
  double (*f)(int, double);
  int n;
  double arg;
  double res;
};

/* 
 * struct for dcomplex(*f)(double, double)
 */
struct dd2c_test
{
  char *label;
  dcomplex (*f)(double, double);
  double arg;
  double nu;
  double res_r, res_i;
};

/* 
 * struct for double(*f)(double, double, double)
 */
struct ddd2d_test
{
  char *label;
  double (*f)(double, double, double);
  double a, b, x;
  double res;
};

/* 
 * struct for double(*f)(double, double, double, double)
 */
struct dddd2d_test
{
  char *label;
  double (*f)(double, double, double, double);
  double a, b, c, x;
  double res;
};

/* 
 * struct for dcomplex(*f)(doubla, dcomplex)
 */
struct dc2c_test
{
  char *label;
  dcomplex (*f)(double, dcomplex);
  double arg_r, arg_i;
  double nu;
  double res_r, res_i;
};

struct d2d_test list_gamma_tst [] =
{
#include "Data_specfun/gamma_test.dat"
    { NULL, NULL, 0, 0}
};

struct id2d_test list_expint_tst [] =
{
#include "Data_specfun/expint_test.dat"
    { NULL, NULL, 0, 0, 0}
};

struct dd2d_test list_gammainc_tst [] =
{
#include "Data_specfun/gammainc_test.dat"
    { NULL, NULL, 0, 0, 0}
};

struct dd2d_test list_real_bessel_tst [] =
{
#include "Data_specfun/real_bessel_test.dat"
    { NULL, NULL, 0, 0, 0}
};

struct dd2c_test list_real_besselh_tst [] =
{
#include "Data_specfun/real_besselh_test.dat"
    { NULL, NULL, 0, 0, 0, 0}
};

struct dc2c_test list_complex_bessel_tst [] =
{
#include "Data_specfun/complex_bessel_test.dat"
    { NULL, NULL, 0, 0, 0, 0, 0}
};

struct dd2d_test list_hyper0F1_tst [] =
{
#include "Data_specfun/hyperg0F1_test.dat"
    { NULL, NULL, 0, 0, 0}
};

struct ddd2d_test list_hyper1F1_tst [] =
{
#include "Data_specfun/hyperg1F1_test.dat"
    { NULL, NULL, 0, 0, 0, 0}
};

struct ddd2d_test list_hyperU_tst [] =
{
#include "Data_specfun/hypergU_test.dat"
    { NULL, NULL, 0, 0, 0, 0}
};

struct ddd2d_test list_hyper2F0_tst [] =
{
#include "Data_specfun/hyperg2F0_test.dat"
    { NULL, NULL, 0, 0, 0, 0}
};

struct dddd2d_test list_hyper2F1_tst [] =
{
#include "Data_specfun/hyperg2F1_test.dat"
    { NULL, NULL, 0, 0, 0, 0, 0}
};

void d2d_funcs_test (struct d2d_test *tst)
{
  int i;
  double tol  = 1E-5;
  for ( i=0 ; tst[i].f != NULL ; i++ )
    {
      struct d2d_test t = tst[i];
      double res = (t.f)(t.arg);
      pnl_test_eq (res, t.res, tol, t.label, 
                  "computed at %g", t.arg);
    }
}

void id2d_funcs_test (struct id2d_test *tst)
{
  int i;
  double tol  = 1E-5;
  for ( i=0 ; tst[i].f != NULL ; i++ )
    {
      struct id2d_test t = tst[i];
      double res = (t.f)(t.n, t.arg);
      pnl_test_eq (res, t.res, tol, t.label, 
                  "computed at (%d, %g)", t.n, t.arg);
    }
}

void dd2d_funcs_test (struct dd2d_test *tst)
{
  int i;
  double tol  = 1E-5;
  for ( i=0 ; tst[i].f != NULL ; i++ )
    {
      struct dd2d_test t = tst[i];
      double res = (t.f)(t.nu, t.arg);
      pnl_test_eq (res, t.res, tol, t.label, 
                  "computed at (%.18g, %.18g)", t.nu, t.arg);
    }
}

void ddd2d_funcs_test (struct ddd2d_test *tst)
{
  int i;
  double tol  = 1E-5;
  for ( i=0 ; tst[i].f != NULL ; i++ )
    {
      struct ddd2d_test t = tst[i];
      double res = (t.f)(t.a, t.b, t.x);
      pnl_test_eq (res, t.res, tol, t.label, 
                   "computed at (%g, %g, %g)", t.a, t.b, t.x);
    }
}

void dddd2d_funcs_test (struct dddd2d_test *tst)
{
  int i;
  double tol  = 1E-5;
  for ( i=0 ; tst[i].f != NULL ; i++ )
    {
      struct dddd2d_test t = tst[i];
      double res = (t.f)(t.a, t.b, t.c, t.x);
      pnl_test_eq (res, t.res, tol, t.label, 
                   "computed at (%g, %g, %g, %g)", t.a, t.b, t.c, t.x);
    }
}

void dd2c_funcs_test (struct dd2c_test *tst)
{
  int i;
  double tol  = 1E-5;
  for ( i=0 ; tst[i].f != NULL ; i++ )
    {
      struct dd2c_test t = tst[i];
      dcomplex res = (t.f)(t.nu, t.arg);
      pnl_test_eq (res.r, t.res_r, tol, t.label, 
                  "computed at (%g,%g+i %g)", t.nu, t.arg);
      pnl_test_eq (res.i, t.res_i, tol, t.label, 
                  "computed at (%g,%g+i %g)", t.nu, t.arg);
    }
}

void dc2c_funcs_test (struct dc2c_test *tst)
{
  int i;
  double tol  = 1E-4;
  for ( i=0 ; tst[i].f != NULL ; i++ )
    {
      struct dc2c_test t = tst[i];
      dcomplex arg = Complex(t.arg_r, t.arg_i);
      dcomplex res = (t.f)(t.nu, arg);
      pnl_test_eq (res.r, t.res_r, tol, t.label, 
                  "computed at (%g,%g+i %g)", t.nu, t.arg_r, t.arg_i);
      pnl_test_eq (res.i, t.res_i, tol, t.label, 
                  "computed at (%g,%g+i %g)", t.nu, t.arg_r, t.arg_i);
    }
}

void test_faddeeva()
{
    {
      int i;
#define NTST 57 // define instead of const for C compatibility
      dcomplex z[NTST] = {
        Complex(624.2,-0.26123),
        Complex(-0.4,3.),
        Complex(0.6,2.),
        Complex(-1.,1.),
        Complex(-1.,-9.),
        Complex(-1.,9.),
        Complex(-0.0000000234545,1.1234),
        Complex(-3.,5.1),
        Complex(-53,30.1),
        Complex(0.0,0.12345),
        Complex(11,1),
        Complex(-22,-2),
        Complex(9,-28),
        Complex(21,-33),
        Complex(1e5,1e5),
        Complex(1e14,1e14),
        Complex(-3001,-1000),
        Complex(1e160,-1e159),
        Complex(-6.01,0.01),
        Complex(-0.7,-0.7),
        Complex(2.611780000000000e+01, 4.540909610972489e+03),
        Complex(0.8e7,0.3e7),
        Complex(-20,-19.8081),
        Complex(1e-16,-1.1e-16),
        Complex(2.3e-8,1.3e-8),
        Complex(6.3,-1e-13),
        Complex(6.3,1e-20),
        Complex(1e-20,6.3),
        Complex(1e-20,16.3),
        Complex(9,1e-300),
        Complex(6.01,0.11),
        Complex(8.01,1.01e-10),
        Complex(28.01,1e-300),
        Complex(10.01,1e-200),
        Complex(10.01,-1e-200),
        Complex(10.01,0.99e-10),
        Complex(10.01,-0.99e-10),
        Complex(1e-20,7.01),
        Complex(-1,7.01),
        Complex(5.99,7.01),
        Complex(1,0),
        Complex(55,0),
        Complex(-0.1,0),
        Complex(1e-20,0),
        Complex(0,5e-14),
        Complex(0,51),
        Complex(PNL_INF,0),
        Complex(-PNL_INF,0),
        Complex(0,PNL_INF),
        Complex(0,-PNL_INF),
        Complex(PNL_INF,PNL_INF),
        Complex(PNL_INF,-PNL_INF),
        Complex(PNL_NAN,PNL_NAN),
        Complex(PNL_NAN,0),
        Complex(0,PNL_NAN),
        Complex(PNL_NAN,PNL_INF),
        Complex(PNL_INF,PNL_NAN)
      };
      dcomplex w[NTST] = { /* w(z), computed with WolframAlpha
                              ... note that WolframAlpha is problematic
                              some of the above inputs, so I had to
                              use the continued-fraction expansion
                              in WolframAlpha in some cases, or switch
                              to Maple */
        Complex(-3.78270245518980507452677445620103199303131110e-7,
                0.000903861276433172057331093754199933411710053155),
        Complex(0.1764906227004816847297495349730234591778719532788,
                -0.02146550539468457616788719893991501311573031095617),
        Complex(0.2410250715772692146133539023007113781272362309451,
                0.06087579663428089745895459735240964093522265589350),
        Complex(0.30474420525691259245713884106959496013413834051768,
                -0.20821893820283162728743734725471561394145872072738),
        Complex(7.317131068972378096865595229600561710140617977e34,
                8.321873499714402777186848353320412813066170427e34),
        Complex(0.0615698507236323685519612934241429530190806818395,
                -0.00676005783716575013073036218018565206070072304635),
        Complex(0.3960793007699874918961319170187598400134746631,
                -5.593152259116644920546186222529802777409274656e-9),
        Complex(0.08217199226739447943295069917990417630675021771804,
                -0.04701291087643609891018366143118110965272615832184),
        Complex(0.00457246000350281640952328010227885008541748668738,
                -0.00804900791411691821818731763401840373998654987934),
        Complex(0.8746342859608052666092782112565360755791467973338452,
                0.),
        Complex(0.00468190164965444174367477874864366058339647648741,
                0.0510735563901306197993676329845149741675029197050),
        Complex(-0.0023193175200187620902125853834909543869428763219,
                -0.025460054739731556004902057663500272721780776336),
        Complex(9.11463368405637174660562096516414499772662584e304,
                3.97101807145263333769664875189354358563218932e305),
        Complex(-4.4927207857715598976165541011143706155432296e281,
                -2.8019591213423077494444700357168707775769028e281),
        Complex(2.820947917809305132678577516325951485807107151e-6,
                2.820947917668257736791638444590253942253354058e-6),
        Complex(2.82094791773878143474039725787438662716372268e-15,
                2.82094791773878143474039725773333923127678361e-15),
        Complex(-0.0000563851289696244350147899376081488003110150498,
                -0.000169211755126812174631861529808288295454992688),
        Complex(-5.586035480670854326218608431294778077663867e-162,
                5.586035480670854326218608431294778077663867e-161),
        Complex(0.00016318325137140451888255634399123461580248456,
                -0.095232456573009287370728788146686162555021209999),
        Complex(0.69504753678406939989115375989939096800793577783885,
                -1.8916411171103639136680830887017670616339912024317),
        Complex(0.0001242418269653279656612334210746733213167234822,
                7.145975826320186888508563111992099992116786763e-7),
        Complex(2.318587329648353318615800865959225429377529825e-8,
                6.182899545728857485721417893323317843200933380e-8),
        Complex(-0.0133426877243506022053521927604277115767311800303,
                -0.0148087097143220769493341484176979826888871576145),
        Complex(1.00000000000000012412170838050638522857747934,
                1.12837916709551279389615890312156495593616433e-16),
        Complex(0.9999999853310704677583504063775310832036830015,
                2.595272024519678881897196435157270184030360773e-8),
        Complex(-1.4731421795638279504242963027196663601154624e-15,
                0.090727659684127365236479098488823462473074709),
        Complex(5.79246077884410284575834156425396800754409308e-18,
                0.0907276596841273652364790985059772809093822374),
        Complex(0.0884658993528521953466533278764830881245144368,
                1.37088352495749125283269718778582613192166760e-22),
        Complex(0.0345480845419190424370085249304184266813447878,
                2.11161102895179044968099038990446187626075258e-23),
        Complex(6.63967719958073440070225527042829242391918213e-36,
                0.0630820900592582863713653132559743161572639353),
        Complex(0.00179435233208702644891092397579091030658500743634,
                0.0951983814805270647939647438459699953990788064762),
        Complex(9.09760377102097999924241322094863528771095448e-13,
                0.0709979210725138550986782242355007611074966717),
        Complex(7.2049510279742166460047102593255688682910274423e-304,
                0.0201552956479526953866611812593266285000876784321),
        Complex(3.04543604652250734193622967873276113872279682e-44,
                0.0566481651760675042930042117726713294607499165),
        Complex(3.04543604652250734193622967873276113872279682e-44,
                0.0566481651760675042930042117726713294607499165),
        Complex(0.5659928732065273429286988428080855057102069081e-12,
                0.056648165176067504292998527162143030538756683302),
        Complex(-0.56599287320652734292869884280802459698927645e-12,
                0.0566481651760675042929985271621430305387566833029),
        Complex(0.0796884251721652215687859778119964009569455462,
                1.11474461817561675017794941973556302717225126e-22),
        Complex(0.07817195821247357458545539935996687005781943386550,
                -0.01093913670103576690766705513142246633056714279654),
        Complex(0.04670032980990449912809326141164730850466208439937,
                0.03944038961933534137558064191650437353429669886545),
        Complex(0.36787944117144232159552377016146086744581113103176,
                0.60715770584139372911503823580074492116122092866515),
        Complex(0,
                0.010259688805536830986089913987516716056946786526145),
        Complex(0.99004983374916805357390597718003655777207908125383,
                -0.11208866436449538036721343053869621153527769495574),
        Complex(0.99999999999999999999999999999999999999990000,
                1.12837916709551257389615890312154517168802603e-20),
        Complex(0.999999999999943581041645226871305192054749891144158,
                0),
        Complex(0.0110604154853277201542582159216317923453996211744250,
                0),
        Complex(0,0),
        Complex(0,0),
        Complex(0,0),
        Complex(PNL_INF,0),
        Complex(0,0),
        Complex(PNL_NAN,PNL_NAN),
        Complex(PNL_NAN,PNL_NAN),
        Complex(PNL_NAN,PNL_NAN),
        Complex(PNL_NAN,0),
        Complex(PNL_NAN,PNL_NAN),
        Complex(PNL_NAN,PNL_NAN)
      };
      double errmax = 1E-13;
      for (i = 0; i < NTST; ++i)
        {
          dcomplex fw = pnl_sf_w(z[i]);
          pnl_test_eq_rel(fw.r, w[i].r, errmax, "pnl_sf_w", "Computed at (%g + i %g)", z[i].r, z[i].i);
        }
    }
  {
#undef NTST
#define NTST 41 // define instead of const for C compatibility
    dcomplex z[NTST] = {
      Complex(1,2),
      Complex(-1,2),
      Complex(1,-2),
      Complex(-1,-2),
      Complex(9,-28),
      Complex(21,-33),
      Complex(1e3,1e3),
      Complex(-3001,-1000),
      Complex(1e160,-1e159),
      Complex(5.1e-3, 1e-8),
      Complex(-4.9e-3, 4.95e-3),
      Complex(4.9e-3, 0.5),
      Complex(4.9e-4, -0.5e1),
      Complex(-4.9e-5, -0.5e2),
      Complex(5.1e-3, 0.5),
      Complex(5.1e-4, -0.5e1),
      Complex(-5.1e-5, -0.5e2),
      Complex(1e-6,2e-6),
      Complex(0,2e-6),
      Complex(0,2),
      Complex(0,20),
      Complex(0,200),
      Complex(PNL_INF,0),
      Complex(-PNL_INF,0),
      Complex(0,PNL_INF),
      Complex(0,-PNL_INF),
      Complex(PNL_INF,PNL_INF),
      Complex(PNL_INF,-PNL_INF),
      Complex(PNL_NAN,PNL_NAN),
      Complex(PNL_NAN,0),
      Complex(0,PNL_NAN),
      Complex(PNL_NAN,PNL_INF),
      Complex(PNL_INF,PNL_NAN),
      Complex(1e-3,PNL_NAN),
      Complex(7e-2,7e-2),
      Complex(7e-2,-7e-4),
      Complex(-9e-2,7e-4),
      Complex(-9e-2,9e-2),
      Complex(-7e-4,9e-2),
      Complex(7e-2,0.9e-2),
      Complex(7e-2,1.1e-2)
    };
    dcomplex w[NTST] = { // erf(z[i]), evaluated with Maple
      Complex(-0.5366435657785650339917955593141927494421,
        -5.049143703447034669543036958614140565553),
      Complex(0.5366435657785650339917955593141927494421,
        -5.049143703447034669543036958614140565553),
      Complex(-0.5366435657785650339917955593141927494421,
        5.049143703447034669543036958614140565553),
      Complex(0.5366435657785650339917955593141927494421,
        5.049143703447034669543036958614140565553),
      Complex(0.3359473673830576996788000505817956637777e304,
        -0.1999896139679880888755589794455069208455e304),
      Complex(0.3584459971462946066523939204836760283645e278,
        0.3818954885257184373734213077678011282505e280),
      Complex(0.9996020422657148639102150147542224526887,
        0.00002801044116908227889681753993542916894856),
      Complex(-1, 0),
      Complex(1, 0),
      Complex(0.005754683859034800134412990541076554934877,
        0.1128349818335058741511924929801267822634e-7),
      Complex(-0.005529149142341821193633460286828381876955,
        0.005585388387864706679609092447916333443570),
      Complex(0.007099365669981359632319829148438283865814,
        0.6149347012854211635026981277569074001219),
      Complex(0.3981176338702323417718189922039863062440e8,
        -0.8298176341665249121085423917575122140650e10),
      Complex(-PNL_INF,
        -PNL_INF),
      Complex(0.007389128308257135427153919483147229573895,
        0.6149332524601658796226417164791221815139),
      Complex(0.4143671923267934479245651547534414976991e8,
        -0.8298168216818314211557046346850921446950e10),
      Complex(-PNL_INF,
        -PNL_INF),
      Complex(0.1128379167099649964175513742247082845155e-5,
        0.2256758334191777400570377193451519478895e-5),
      Complex(0,
        0.2256758334194034158904576117253481476197e-5),
      Complex(0,
        18.56480241457555259870429191324101719886),
      Complex(0,
        0.1474797539628786202447733153131835124599e173),
      Complex(0,
        PNL_INF),
      Complex(1,0),
      Complex(-1,0),
      Complex(0,PNL_INF),
      Complex(0,-PNL_INF),
      Complex(PNL_NAN,PNL_NAN),
      Complex(PNL_NAN,PNL_NAN),
      Complex(PNL_NAN,PNL_NAN),
      Complex(PNL_NAN,0),
      Complex(0,PNL_NAN),
      Complex(PNL_NAN,PNL_NAN),
      Complex(PNL_NAN,PNL_NAN),
      Complex(PNL_NAN,PNL_NAN),
      Complex(0.07924380404615782687930591956705225541145,
        0.07872776218046681145537914954027729115247),
      Complex(0.07885775828512276968931773651224684454495,
        -0.0007860046704118224342390725280161272277506),
      Complex(-0.1012806432747198859687963080684978759881,
        0.0007834934747022035607566216654982820299469),
      Complex(-0.1020998418798097910247132140051062512527,
        0.1010030778892310851309082083238896270340),
      Complex(-0.0007962891763147907785684591823889484764272,
        0.1018289385936278171741809237435404896152),
      Complex(0.07886408666470478681566329888615410479530,
        0.01010604288780868961492224347707949372245),
      Complex(0.07886723099940260286824654364807981336591,
        0.01235199327873258197931147306290916629654)
    };
#define TST(f,isc)                                                                                             \
    double errmax = 1E-13;                                                                                     \
    int i;                                                                                                     \
    for (i = 0; i < NTST; ++i) {                                                                               \
      dcomplex fw = pnl_sf_complex_##f(z[i]);                                                                  \
        pnl_test_eq_rel(fw.r, w[i].r, errmax, "pnl_sf_complex_" #f, "Computed at (%g + i %g)", z[i].r, z[i].i); \
    }

    TST(erf, 1e-20);
  }
  {
    // since erfi just calls through to erf, just one test should
    // be sufficient to make sure I didn't screw up the signs or something
#undef NTST
#define NTST 1 // define instead of const for C compatibility
    dcomplex z[NTST] = { Complex(1.234,0.5678) };
    dcomplex w[NTST] = { // erfi(z[i]), computed with Maple
      Complex(1.081032284405373149432716643834106923212,
        1.926775520840916645838949402886591180834)
    };
    TST(erfi, 0);
  }
  {
    // since erfcx just calls through to w, just one test should
    // be sufficient to make sure I didn't screw up the signs or something
#undef NTST
#define NTST 1 // define instead of const for C compatibility
    dcomplex z[NTST] = { Complex(1.234,0.5678) };
    dcomplex w[NTST] = { // erfcx(z[i]), computed with Maple
      Complex(0.3382187479799972294747793561190487832579,
        -0.1116077470811648467464927471872945833154)
    };
    TST(erfcx, 0);
  }
  {
#undef NTST
#define NTST 30 // define instead of const for C compatibility
    dcomplex z[NTST] = {
      Complex(1,2),
      Complex(-1,2),
      Complex(1,-2),
      Complex(-1,-2),
      Complex(9,-28),
      Complex(21,-33),
      Complex(1e3,1e3),
      Complex(-3001,-1000),
      Complex(1e160,-1e159),
      Complex(5.1e-3, 1e-8),
      Complex(0,2e-6),
      Complex(0,2),
      Complex(0,20),
      Complex(0,200),
      Complex(2e-6,0),
      Complex(2,0),
      Complex(20,0),
      Complex(200,0),
      Complex(PNL_INF,0),
      Complex(-PNL_INF,0),
      Complex(0,PNL_INF),
      Complex(0,-PNL_INF),
      Complex(PNL_INF,PNL_INF),
      Complex(PNL_INF,-PNL_INF),
      Complex(PNL_NAN,PNL_NAN),
      Complex(PNL_NAN,0),
      Complex(0,PNL_NAN),
      Complex(PNL_NAN,PNL_INF),
      Complex(PNL_INF,PNL_NAN),
      Complex(88,0)
    };
    dcomplex w[NTST] = { // erfc(z[i]), evaluated with Maple
      Complex(1.536643565778565033991795559314192749442,
        5.049143703447034669543036958614140565553),
      Complex(0.4633564342214349660082044406858072505579,
        5.049143703447034669543036958614140565553),
      Complex(1.536643565778565033991795559314192749442,
        -5.049143703447034669543036958614140565553),
      Complex(0.4633564342214349660082044406858072505579,
        -5.049143703447034669543036958614140565553),
      Complex(-0.3359473673830576996788000505817956637777e304,
        0.1999896139679880888755589794455069208455e304),
      Complex(-0.3584459971462946066523939204836760283645e278,
        -0.3818954885257184373734213077678011282505e280),
      Complex(0.0003979577342851360897849852457775473112748,
        -0.00002801044116908227889681753993542916894856),
      Complex(2, 0),
      Complex(0, 0),
      Complex(0.9942453161409651998655870094589234450651,
        -0.1128349818335058741511924929801267822634e-7),
      Complex(1,
        -0.2256758334194034158904576117253481476197e-5),
      Complex(1,
        -18.56480241457555259870429191324101719886),
      Complex(1,
        -0.1474797539628786202447733153131835124599e173),
      Complex(1, -PNL_INF),
      Complex(0.9999977432416658119838633199332831406314,
        0),
      Complex(0.004677734981047265837930743632747071389108,
        0),
      Complex(0.5395865611607900928934999167905345604088e-175,
        0),
      Complex(0, 0),
      Complex(0, 0),
      Complex(2, 0),
      Complex(1, -PNL_INF),
      Complex(1, PNL_INF),
      Complex(PNL_NAN, PNL_NAN),
      Complex(PNL_NAN, PNL_NAN),
      Complex(PNL_NAN, PNL_NAN),
      Complex(PNL_NAN, 0),
      Complex(1, PNL_NAN),
      Complex(PNL_NAN, PNL_NAN),
      Complex(PNL_NAN, PNL_NAN),
      Complex(0,0)
    };
    TST(erfc, 1e-20);
  }
  {
#undef NTST
#define NTST 48 // define instead of const for C compatibility
    dcomplex z[NTST] = {
      Complex(2,1),
      Complex(-2,1),
      Complex(2,-1),
      Complex(-2,-1),
      Complex(-28,9),
      Complex(33,-21),
      Complex(1e3,1e3),
      Complex(-1000,-3001),
      Complex(1e-8, 5.1e-3),
      Complex(4.95e-3, -4.9e-3),
      Complex(5.1e-3, 5.1e-3),
      Complex(0.5, 4.9e-3),
      Complex(-0.5e1, 4.9e-4),
      Complex(-0.5e2, -4.9e-5),
      Complex(0.5e3, 4.9e-6),
      Complex(0.5, 5.1e-3),
      Complex(-0.5e1, 5.1e-4),
      Complex(-0.5e2, -5.1e-5),
      Complex(1e-6,2e-6),
      Complex(2e-6,0),
      Complex(2,0),
      Complex(20,0),
      Complex(200,0),
      Complex(0,4.9e-3),
      Complex(0,-5.1e-3),
      Complex(0,2e-6),
      Complex(0,-2),
      Complex(0,20),
      Complex(0,-200),
      Complex(PNL_INF,0),
      Complex(-PNL_INF,0),
      Complex(0,PNL_INF),
      Complex(0,-PNL_INF),
      Complex(PNL_INF,PNL_INF),
      Complex(PNL_INF,-PNL_INF),
      Complex(PNL_NAN,PNL_NAN),
      Complex(PNL_NAN,0),
      Complex(0,PNL_NAN),
      Complex(PNL_NAN,PNL_INF),
      Complex(PNL_INF,PNL_NAN),
      Complex(39, 6.4e-5),
      Complex(41, 6.09e-5),
      Complex(4.9e7, 5e-11),
      Complex(5.1e7, 4.8e-11),
      Complex(1e9, 2.4e-12),
      Complex(1e11, 2.4e-14),
      Complex(1e13, 2.4e-16),
      Complex(1e300, 2.4e-303)
    };
    dcomplex w[NTST] = { // dawson(z[i]), evaluated with Maple
      Complex(0.1635394094345355614904345232875688576839,
        -0.1531245755371229803585918112683241066853),
      Complex(-0.1635394094345355614904345232875688576839,
        -0.1531245755371229803585918112683241066853),
      Complex(0.1635394094345355614904345232875688576839,
        0.1531245755371229803585918112683241066853),
      Complex(-0.1635394094345355614904345232875688576839,
        0.1531245755371229803585918112683241066853),
      Complex(-0.01619082256681596362895875232699626384420,
        -0.005210224203359059109181555401330902819419),
      Complex(0.01078377080978103125464543240346760257008,
        0.006866888783433775382193630944275682670599),
      Complex(-0.5808616819196736225612296471081337245459,
        0.6688593905505562263387760667171706325749),
      Complex(PNL_INF,
        -PNL_INF),
      Complex(0.1000052020902036118082966385855563526705e-7,
        0.005100088434920073153418834680320146441685),
      Complex(0.004950156837581592745389973960217444687524,
        -0.004899838305155226382584756154100963570500),
      Complex(0.005100176864319675957314822982399286703798,
        0.005099823128319785355949825238269336481254),
      Complex(0.4244534840871830045021143490355372016428,
        0.002820278933186814021399602648373095266538),
      Complex(-0.1021340733271046543881236523269967674156,
        -0.00001045696456072005761498961861088944159916),
      Complex(-0.01000200120119206748855061636187197886859,
        0.9805885888237419500266621041508714123763e-8),
      Complex(0.001000002000012000023960527532953151819595,
        -0.9800058800588007290937355024646722133204e-11),
      Complex(0.4244549085628511778373438768121222815752,
        0.002935393851311701428647152230552122898291),
      Complex(-0.1021340732357117208743299813648493928105,
        -0.00001088377943049851799938998805451564893540),
      Complex(-0.01000200120119126652710792390331206563616,
        0.1020612612857282306892368985525393707486e-7),
      Complex(0.1000000000007333333333344266666666664457e-5,
        0.2000000000001333333333323199999999978819e-5),
      Complex(0.1999999999994666666666675199999999990248e-5,
        0),
      Complex(0.3013403889237919660346644392864226952119,
        0),
      Complex(0.02503136792640367194699495234782353186858,
        0),
      Complex(0.002500031251171948248596912483183760683918,
        0),
      Complex(0,0.004900078433419939164774792850907128053308),
      Complex(0,-0.005100088434920074173454208832365950009419),
      Complex(0,0.2000000000005333333333341866666666676419e-5),
      Complex(0,-48.16001211429122974789822893525016528191),
      Complex(0,0.4627407029504443513654142715903005954668e174),
      Complex(0,-PNL_INF),
      Complex(0,0),
      Complex(-0,0),
      Complex(0, PNL_INF),
      Complex(0, -PNL_INF),
      Complex(PNL_NAN, PNL_NAN),
      Complex(PNL_NAN, PNL_NAN),
      Complex(PNL_NAN, PNL_NAN),
      Complex(PNL_NAN, 0),
      Complex(0, PNL_NAN),
      Complex(PNL_NAN, PNL_NAN),
      Complex(PNL_NAN, PNL_NAN),
      Complex(0.01282473148489433743567240624939698290584,
        -0.2105957276516618621447832572909153498104e-7),
      Complex(0.01219875253423634378984109995893708152885,
        -0.1813040560401824664088425926165834355953e-7),
      Complex(0.1020408163265306334945473399689037886997e-7,
        -0.1041232819658476285651490827866174985330e-25),
      Complex(0.9803921568627452865036825956835185367356e-8,
        -0.9227220299884665067601095648451913375754e-26),
      Complex(0.5000000000000000002500000000000000003750e-9,
        -0.1200000000000000001800000188712838420241e-29),
      Complex(5.00000000000000000000025000000000000000000003e-12,
        -1.20000000000000000000018000000000000000000004e-36),
      Complex(5.00000000000000000000000002500000000000000000e-14,
        -1.20000000000000000000000001800000000000000000e-42),
      Complex(5e-301, 0)
    };
    TST(dawson, 1e-20);
  }
}

int main (int argc, char **argv)
{
  pnl_test_init (argc, argv);
  d2d_funcs_test (list_gamma_tst);
  id2d_funcs_test (list_expint_tst);
  dd2d_funcs_test (list_gammainc_tst);

  dd2d_funcs_test (list_real_bessel_tst);
  dd2c_funcs_test (list_real_besselh_tst);
  dc2c_funcs_test (list_complex_bessel_tst);

  dd2d_funcs_test (list_hyper0F1_tst);
  ddd2d_funcs_test (list_hyper1F1_tst);
  ddd2d_funcs_test (list_hyperU_tst);
  ddd2d_funcs_test (list_hyper2F0_tst);
  dddd2d_funcs_test (list_hyper2F1_tst);
  test_faddeeva();

  exit (pnl_test_finalize("Special functions"));
}
