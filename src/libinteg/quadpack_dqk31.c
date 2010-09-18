/*
 * This file contains the routine dqk31 from the quadpack distribution
 * available on www.netlib.org
 *
 * This routine has been translated from Fortran using f2c and then slightly
 * modified by Jérôme Lelong for integration into Pnl. All Quadpack names have
 * been prefixed by pnl to avoid possible clash in names. Calls to d1mach have
 * been replaced by calls to dlamch coming from Lapack
 */

#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"

int pnl_dqk31(PnlFunc * f, double *a, double *b, double *
	result, double *abserr, double *resabs, double *resasc)
{
    /* Initialized data */

    double wg[8] = { .030753241996117268354628393577204,
	    .070366047488108124709267416450667,
	    .107159220467171935011869546685869,
	    .139570677926154314447804794511028,
	    .166269205816993933553200860481209,
	    .186161000015562211026800561866423,
	    .198431485327111576456118326443839,
	    .202578241925561272880620199967519 };
    double xgk[16] = { .998002298693397060285172840152271,
	    .987992518020485428489565718586613,
	    .967739075679139134257347978784337,
	    .937273392400705904307758947710209,
	    .897264532344081900882509656454496,
	    .848206583410427216200648320774217,
	    .790418501442465932967649294817947,
	    .724417731360170047416186054613938,
	    .650996741297416970533735895313275,
	    .570972172608538847537226737253911,
	    .485081863640239680693655740232351,
	    .394151347077563369897207370981045,
	    .299180007153168812166780024266389,
	    .201194093997434522300628303394596,
	    .101142066918717499027074231447392,0. };
    double wgk[16] = { .005377479872923348987792051430128,
	    .015007947329316122538374763075807,
	    .025460847326715320186874001019653,
	    .03534636079137584622203794847836,
	    .04458975132476487660822729937328,
	    .05348152469092808726534314723943,
	    .062009567800670640285139230960803,
	    .069854121318728258709520077099147,
	    .076849680757720378894432777482659,
	    .083080502823133021038289247286104,
	    .088564443056211770647275443693774,
	    .093126598170825321225486872747346,
	    .096642726983623678505179907627589,
	    .099173598721791959332393173484603,
	    .10076984552387559504494666261757,
	    .101330007014791549017374792767493 };

    double c_b7 = 1.5;

    /* System generated locals */
    double d__1, d__2, d__3;

    /* Local variables */
    int j;
    double fc, fv1[15], fv2[15];
    int jtw;
    double absc, resg, resk, fsum, fval1, fval2;
    int jtwm1;
    double hlgth, centr, reskh, uflow;
    
    double epmach, dhlgth;

/* ***begin prologue  dqk31 */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a1a2 */
/* ***keywords  31-point gauss-kronrod rules */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  to compute i = integral of f over (a,b) with error */
/*                           estimate */
/*                       j = integral of ABS(f) over (a,b) */
/* ***description */

/*           integration rules */
/*           standard fortran subroutine */
/*           double precision version */

/*           parameters */
/*            on entry */
/*              f      - double precision */
/*                       function subprogram defining the integrand */
/*                       function f(x). the actual name for f needs to be */
/*                       declared e x t e r n a l in the calling program. */

/*              a      - double precision */
/*                       lower limit of integration */

/*              b      - double precision */
/*                       upper limit of integration */

/*            on return */
/*              result - double precision */
/*                       approximation to the integral i */
/*                       result is computed by applying the 31-point */
/*                       gauss-kronrod rule (resk), obtained by optimal */
/*                       addition of abscissae to the 15-point gauss */
/*                       rule (resg). */

/*              abserr - double precison */
/*                       estimate of the modulus of the modulus, */
/*                       which should not exceed ABS(i-result) */

/*              resabs - double precision */
/*                       approximation to the integral j */

/*              resasc - double precision */
/*                       approximation to the integral of ABS(f-i/(b-a)) */
/*                       over (a,b) */

/* ***references  (none) */
/* ***end prologue  dqk31 */


/*           the abscissae and weights are given for the interval (-1,1). */
/*           because of symmetry only the positive abscissae and their */
/*           corresponding weights are given. */

/*           xgk    - abscissae of the 31-point kronrod rule */
/*                    xgk(2), xgk(4), ...  abscissae of the 15-point */
/*                    gauss rule */
/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
/*                    added to the 15-point gauss rule */

/*           wgk    - weights of the 31-point kronrod rule */

/*           wg     - weights of the 15-point gauss rule */


/* gauss quadrature weights and kronron quadrature abscissae and weights */
/* as evaluated with 80 decimal digit arithmetic by l. w. fullerton, */
/* bell labs, nov. 1981. */





/*           list of major variables */
/*           ----------------------- */
/*           centr  - mid point of the interval */
/*           hlgth  - half-length of the interval */
/*           absc   - abscissa */
/*           fval*  - function value */
/*           resg   - result of the 15-point gauss formula */
/*           resk   - result of the 31-point kronrod formula */
/*           reskh  - approximation to the mean value of f over (a,b), */
/*                    i.e. to i/(b-a) */

/*           machine dependent constants */
/*           --------------------------- */
/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */
/* ***first executable statement  dqk31 */
    epmach = pnl_d1mach(4);
    uflow = pnl_d1mach(1);

    centr = (*a + *b) * .5;
    hlgth = (*b - *a) * .5;
    dhlgth = ABS(hlgth);

/*           compute the 31-point kronrod approximation to */
/*           the integral, and estimate the absolute error. */

    fc = PNL_EVAL_FUNC(f,centr);
    resg = wg[7] * fc;
    resk = wgk[15] * fc;
    *resabs = ABS(resk);
    for (j = 1; j <= 7; ++j) {
	jtw = j << 1;
	absc = hlgth * xgk[jtw - 1];
	d__1 = centr - absc;
	fval1 = PNL_EVAL_FUNC(f,d__1);
	d__1 = centr + absc;
	fval2 = PNL_EVAL_FUNC(f,d__1);
	fv1[jtw - 1] = fval1;
	fv2[jtw - 1] = fval2;
	fsum = fval1 + fval2;
	resg += wg[j - 1] * fsum;
	resk += wgk[jtw - 1] * fsum;
	*resabs += wgk[jtw - 1] * (ABS(fval1) + ABS(fval2));
/* L10: */
    }
    for (j = 1; j <= 8; ++j) {
	jtwm1 = (j << 1) - 1;
	absc = hlgth * xgk[jtwm1 - 1];
	d__1 = centr - absc;
	fval1 = PNL_EVAL_FUNC(f,d__1);
	d__1 = centr + absc;
	fval2 = PNL_EVAL_FUNC(f,d__1);
	fv1[jtwm1 - 1] = fval1;
	fv2[jtwm1 - 1] = fval2;
	fsum = fval1 + fval2;
	resk += wgk[jtwm1 - 1] * fsum;
	*resabs += wgk[jtwm1 - 1] * (ABS(fval1) + ABS(fval2));
/* L15: */
    }
    reskh = resk * .5;
    *resasc = wgk[15] * (d__1 = fc - reskh, ABS(d__1));
    for (j = 1; j <= 15; ++j) {
	*resasc += wgk[j - 1] * ((d__1 = fv1[j - 1] - reskh, ABS(d__1)) + (
		d__2 = fv2[j - 1] - reskh, ABS(d__2)));
/* L20: */
    }
    *result = resk * hlgth;
    *resabs *= dhlgth;
    *resasc *= dhlgth;
    *abserr = (d__1 = (resk - resg) * hlgth, ABS(d__1));
    if (*resasc != 0. && *abserr != 0.) {
/* Computing MIN */
	d__3 = *abserr * 200. / *resasc;
	d__1 = 1., d__2 = pow(d__3, c_b7);
	*abserr = *resasc * MIN(d__1,d__2);
    }
    if (*resabs > uflow / (epmach * 50.)) {
/* Computing MAX */
	d__1 = epmach * 50. * *resabs;
	*abserr = MAX(d__1,*abserr);
    }
    return 0;
} /* pnl_dqk31 */

