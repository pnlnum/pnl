/*
 * This file contains the routine dqk41 from the quadpack distribution
 * available on www.netlib.org
 *
 * This routine has been translated from Fortran using f2c and then slightly
 * modified by Jérôme Lelong for integration into Pnl. All Quadpack names have
 * been prefixed by pnl to avoid possible clash in names. Calls to d1mach have
 * been replaced by calls to dlamch coming from Lapack
 */

#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"

int pnl_dqk41(PnlFunc * f, double *a, double *b, double *
	result, double *abserr, double *resabs, double *resasc)
{
    /* Initialized data */

    double wg[10] = { .017614007139152118311861962351853,
	    .040601429800386941331039952274932,
	    .062672048334109063569506535187042,
	    .083276741576704748724758143222046,
	    .10193011981724043503675013548035,
	    .118194531961518417312377377711382,
	    .131688638449176626898494499748163,
	    .142096109318382051329298325067165,
	    .149172986472603746787828737001969,
	    .152753387130725850698084331955098 };
    double xgk[21] = { .998859031588277663838315576545863,
	    .99312859918509492478612238847132,
	    .981507877450250259193342994720217,
	    .963971927277913791267666131197277,
	    .940822633831754753519982722212443,
	    .912234428251325905867752441203298,
	    .878276811252281976077442995113078,
	    .839116971822218823394529061701521,
	    .795041428837551198350638833272788,
	    .746331906460150792614305070355642,
	    .693237656334751384805490711845932,
	    .636053680726515025452836696226286,
	    .575140446819710315342946036586425,
	    .510867001950827098004364050955251,
	    .44359317523872510319999221349264,
	    .373706088715419560672548177024927,
	    .301627868114913004320555356858592,
	    .227785851141645078080496195368575,
	    .152605465240922675505220241022678,
	    .076526521133497333754640409398838,0. };
    double wgk[21] = { .003073583718520531501218293246031,
	    .008600269855642942198661787950102,
	    .014626169256971252983787960308868,
	    .020388373461266523598010231432755,
	    .025882133604951158834505067096153,
	    .031287306777032798958543119323801,
	    .036600169758200798030557240707211,
	    .041668873327973686263788305936895,
	    .046434821867497674720231880926108,
	    .050944573923728691932707670050345,
	    .055195105348285994744832372419777,
	    .059111400880639572374967220648594,
	    .062653237554781168025870122174255,
	    .065834597133618422111563556969398,
	    .068648672928521619345623411885368,
	    .07105442355344406830579036172321,
	    .073030690332786667495189417658913,
	    .074582875400499188986581418362488,
	    .075704497684556674659542775376617,
	    .076377867672080736705502835038061,
	    .076600711917999656445049901530102 };

    double c_b7 = 1.5;

    /* System generated locals */
    double d__1, d__2, d__3;

    /* Local variables */
    int j;
    double fc, fv1[20], fv2[20];
    int jtw;
    double absc, resg, resk, fsum, fval1, fval2;
    int jtwm1;
    double hlgth, centr, reskh, uflow;
    
    double epmach, dhlgth;

/* ***begin prologue  dqk41 */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a1a2 */
/* ***keywords  41-point gauss-kronrod rules */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  to compute i = integral of f over (a,b), with error */
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
/*                       result is computed by applying the 41-point */
/*                       gauss-kronrod rule (resk) obtained by optimal */
/*                       addition of abscissae to the 20-point gauss */
/*                       rule (resg). */

/*              abserr - double precision */
/*                       estimate of the modulus of the absolute error, */
/*                       which should not exceed ABS(i-result) */

/*              resabs - double precision */
/*                       approximation to the integral j */

/*              resasc - double precision */
/*                       approximation to the integal of ABS(f-i/(b-a)) */
/*                       over (a,b) */

/* ***references  (none) */
/* ***end prologue  dqk41 */



/*           the abscissae and weights are given for the interval (-1,1). */
/*           because of symmetry only the positive abscissae and their */
/*           corresponding weights are given. */

/*           xgk    - abscissae of the 41-point gauss-kronrod rule */
/*                    xgk(2), xgk(4), ...  abscissae of the 20-point */
/*                    gauss rule */
/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
/*                    added to the 20-point gauss rule */

/*           wgk    - weights of the 41-point gauss-kronrod rule */

/*           wg     - weights of the 20-point gauss rule */


/* gauss quadrature weights and kronron quadrature abscissae and weights */
/* as evaluated with 80 decimal digit arithmetic by l. w. fullerton, */
/* bell labs, nov. 1981. */





/*           list of major variables */
/*           ----------------------- */

/*           centr  - mid point of the interval */
/*           hlgth  - half-length of the interval */
/*           absc   - abscissa */
/*           fval*  - function value */
/*           resg   - result of the 20-point gauss formula */
/*           resk   - result of the 41-point kronrod formula */
/*           reskh  - approximation to mean value of f over (a,b), i.e. */
/*                    to i/(b-a) */

/*           machine dependent constants */
/*           --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */

/* ***first executable statement  dqk41 */
    epmach = pnl_d1mach(4);
    uflow = pnl_d1mach(1);

    centr = (*a + *b) * .5;
    hlgth = (*b - *a) * .5;
    dhlgth = ABS(hlgth);

/*           compute the 41-point gauss-kronrod approximation to */
/*           the integral, and estimate the absolute error. */

    resg = 0.;
    fc = PNL_EVAL_FUNC(f,centr);
    resk = wgk[20] * fc;
    *resabs = ABS(resk);
    for (j = 1; j <= 10; ++j) {
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
    for (j = 1; j <= 10; ++j) {
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
    *resasc = wgk[20] * (d__1 = fc - reskh, ABS(d__1));
    for (j = 1; j <= 20; ++j) {
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
} /* pnl_dqk41 */

