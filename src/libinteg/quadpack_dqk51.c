/*
 * This file contains the routine dqk51 from the quadpack distribution
 * available on www.netlib.org
 *
 * This routine has been translated from Fortran using f2c and then slightly
 * modified by Jérôme Lelong for integration into Pnl. All Quadpack names have
 * been prefixed by pnl to avoid possible clash in names. Calls to d1mach have
 * been replaced by calls to dlamch coming from Lapack
 */

#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"

int pnl_dqk51(PnlFunc * f, double *a, double *b, double *
	result, double *abserr, double *resabs, double *resasc)
{
    /* Initialized data */

    double wg[13] = { .011393798501026287947902964113235,
	    .026354986615032137261901815295299,
	    .040939156701306312655623487711646,
	    .054904695975835191925936891540473,
	    .068038333812356917207187185656708,
	    .080140700335001018013234959669111,
	    .091028261982963649811497220702892,
	    .100535949067050644202206890392686,
	    .108519624474263653116093957050117,
	    .114858259145711648339325545869556,
	    .119455763535784772228178126512901,
	    .122242442990310041688959518945852,
	    .12317605372671545120390287307905 };
    double xgk[26] = { .999262104992609834193457486540341,
	    .995556969790498097908784946893902,
	    .988035794534077247637331014577406,
	    .976663921459517511498315386479594,
	    .961614986425842512418130033660167,
	    .942974571228974339414011169658471,
	    .920747115281701561746346084546331,
	    .894991997878275368851042006782805,
	    .86584706529327559544899696958834,
	    .83344262876083400142102110869357,
	    .797873797998500059410410904994307,
	    .759259263037357630577282865204361,
	    .717766406813084388186654079773298,
	    .673566368473468364485120633247622,
	    .626810099010317412788122681624518,
	    .577662930241222967723689841612654,
	    .52632528433471918259962377815801,
	    .473002731445714960522182115009192,
	    .417885382193037748851814394594572,
	    .361172305809387837735821730127641,
	    .303089538931107830167478909980339,
	    .243866883720988432045190362797452,
	    .183718939421048892015969888759528,
	    .122864692610710396387359818808037,
	    .061544483005685078886546392366797,0. };
    double wgk[26] = { .001987383892330315926507851882843,
	    .005561932135356713758040236901066,
	    .009473973386174151607207710523655,
	    .013236229195571674813656405846976,
	    .016847817709128298231516667536336,
	    .020435371145882835456568292235939,
	    .024009945606953216220092489164881,
	    .027475317587851737802948455517811,
	    .030792300167387488891109020215229,
	    .034002130274329337836748795229551,
	    .03711627148341554356033062536762,
	    .040083825504032382074839284467076,
	    .042872845020170049476895792439495,
	    .04550291304992178890987058475266,
	    .047982537138836713906392255756915,
	    .05027767908071567196332525943344,
	    .052362885806407475864366712137873,
	    .054251129888545490144543370459876,
	    .055950811220412317308240686382747,
	    .057437116361567832853582693939506,
	    .058689680022394207961974175856788,
	    .059720340324174059979099291932562,
	    .060539455376045862945360267517565,
	    .061128509717053048305859030416293,
	    .061471189871425316661544131965264,
	    .061580818067832935078759824240066 };

    double c_b7 = 1.5;

    /* System generated locals */
    double d__1, d__2, d__3;

    /* Local variables */
    int j;
    double fc, fv1[25], fv2[25];
    int jtw;
    double absc, resg, resk, fsum, fval1, fval2;
    int jtwm1;
    double hlgth, centr, reskh, uflow;
    
    double epmach, dhlgth;

/* ***begin prologue  dqk51 */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a1a2 */
/* ***keywords  51-point gauss-kronrod rules */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math & progr. div. - k.u.leuven */
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
/*                       function subroutine defining the integrand */
/*                       function f(x). the actual name for f needs to be */
/*                       declared e x t e r n a l in the calling program. */

/*              a      - double precision */
/*                       lower limit of integration */

/*              b      - double precision */
/*                       upper limit of integration */

/*            on return */
/*              result - double precision */
/*                       approximation to the integral i */
/*                       result is computed by applying the 51-point */
/*                       kronrod rule (resk) obtained by optimal addition */
/*                       of abscissae to the 25-point gauss rule (resg). */

/*              abserr - double precision */
/*                       estimate of the modulus of the absolute error, */
/*                       which should not exceed ABS(i-result) */

/*              resabs - double precision */
/*                       approximation to the integral j */

/*              resasc - double precision */
/*                       approximation to the integral of ABS(f-i/(b-a)) */
/*                       over (a,b) */

/* ***references  (none) */
/* ***end prologue  dqk51 */



/*           the abscissae and weights are given for the interval (-1,1). */
/*           because of symmetry only the positive abscissae and their */
/*           corresponding weights are given. */

/*           xgk    - abscissae of the 51-point kronrod rule */
/*                    xgk(2), xgk(4), ...  abscissae of the 25-point */
/*                    gauss rule */
/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
/*                    added to the 25-point gauss rule */

/*           wgk    - weights of the 51-point kronrod rule */

/*           wg     - weights of the 25-point gauss rule */


/* gauss quadrature weights and kronron quadrature abscissae and weights */
/* as evaluated with 80 decimal digit arithmetic by l. w. fullerton, */
/* bell labs, nov. 1981. */



/*       note: wgk (26) was calculated from the values of wgk(1..25) */


/*           list of major variables */
/*           ----------------------- */

/*           centr  - mid point of the interval */
/*           hlgth  - half-length of the interval */
/*           absc   - abscissa */
/*           fval*  - function value */
/*           resg   - result of the 25-point gauss formula */
/*           resk   - result of the 51-point kronrod formula */
/*           reskh  - approximation to the mean value of f over (a,b), */
/*                    i.e. to i/(b-a) */

/*           machine dependent constants */
/*           --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */

/* ***first executable statement  dqk51 */
    epmach = pnl_d1mach(4);
    uflow = pnl_d1mach(1);

    centr = (*a + *b) * .5;
    hlgth = (*b - *a) * .5;
    dhlgth = ABS(hlgth);

/*           compute the 51-point kronrod approximation to */
/*           the integral, and estimate the absolute error. */

    fc = PNL_EVAL_FUNC(f,centr);
    resg = wg[12] * fc;
    resk = wgk[25] * fc;
    *resabs = ABS(resk);
    for (j = 1; j <= 12; ++j) {
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
    for (j = 1; j <= 13; ++j) {
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
    *resasc = wgk[25] * (d__1 = fc - reskh, ABS(d__1));
    for (j = 1; j <= 25; ++j) {
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
} /* pnl_dqk51 */

