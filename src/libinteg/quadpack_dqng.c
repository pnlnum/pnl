/*
 * This file contains the routine dqng from the quadpack distribution
 * available on www.netlib.org
 *
 * This routine has been translated from Fortran using f2c and then slightly
 * modified by Jérôme Lelong for integration into Pnl. All Quadpack names have
 * been prefixed by pnl to avoid possible clash in names. Calls to d1mach have
 * been replaced by calls to dlamch coming from Lapack
 */

#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"



int pnl_dqng(PnlFunc *f, double *a, double *b, double *
	epsabs, double *epsrel, double *result, double *abserr, 
	int *neval, int *ier)
{
    /* Initialized data */

    double x1[5] = { .973906528517171720077964012084452,
	    .865063366688984510732096688423493,
	    .679409568299024406234327365114874,
	    .433395394129247190799265943165784,
	    .14887433898163121088482600112972 };
    double w87a[21] = { .00814837738414917290000287844819,
	    .018761438201562822243935059003794,
	    .027347451050052286161582829741283,
	    .033677707311637930046581056957588,
	    .036935099820427907614589586742499,
	    .002884872430211530501334156248695,
	    .013685946022712701888950035273128,
	    .023280413502888311123409291030404,
	    .030872497611713358675466394126442,
	    .035693633639418770719351355457044,
	    9.15283345202241360843392549948e-4,
	    .005399280219300471367738743391053,
	    .010947679601118931134327826856808,
	    .01629873169678733526266570322328,
	    .02108156888920383511243306018819,
	    .02537096976925382724346799983171,
	    .02918969775647575250144615408492,
	    .032373202467202789685788194889595,
	    .034783098950365142750781997949596,
	    .036412220731351787562801163687577,
	    .037253875503047708539592001191226 };
    double w87b[23] = { 2.74145563762072350016527092881e-4,
	    .001807124155057942948341311753254,
	    .00409686928275916486445807068348,
	    .006758290051847378699816577897424,
	    .009549957672201646536053581325377,
	    .01232944765224485369462663996378,
	    .015010447346388952376697286041943,
	    .0175489679862431910996653529259,
	    .019938037786440888202278192730714,
	    .022194935961012286796332102959499,
	    .024339147126000805470360647041454,
	    .026374505414839207241503786552615,
	    .02828691078877120065996800298796,
	    .030052581128092695322521110347341,
	    .031646751371439929404586051078883,
	    .033050413419978503290785944862689,
	    .034255099704226061787082821046821,
	    .035262412660156681033782717998428,
	    .036076989622888701185500318003895,
	    .036698604498456094498018047441094,
	    .037120549269832576114119958413599,
	    .037334228751935040321235449094698,
	    .037361073762679023410321241766599 };
    double w10[5] = { .066671344308688137593568809893332,
	    .149451349150580593145776339657697,
	    .219086362515982043995534934228163,
	    .269266719309996355091226921569469,
	    .295524224714752870173892994651338 };
    double x2[5] = { .995657163025808080735527280689003,
	    .930157491355708226001207180059508,
	    .780817726586416897063717578345042,
	    .562757134668604683339000099272694,
	    .294392862701460198131126603103866 };
    double w21a[5] = { .03255816230796472747881897245939,
	    .07503967481091995276704314091619,
	    .109387158802297641899210590325805,
	    .134709217311473325928054001771707,
	    .147739104901338491374841515972068 };
    double w21b[6] = { .011694638867371874278064396062192,
	    .05475589657435199603138130024458,
	    .093125454583697605535065465083366,
	    .123491976262065851077958109831074,
	    .142775938577060080797094273138717,
	    .149445554002916905664936468389821 };
    double x3[11] = { .999333360901932081394099323919911,
	    .987433402908088869795961478381209,
	    .954807934814266299257919200290473,
	    .900148695748328293625099494069092,
	    .82519831498311415084706673258852,
	    .732148388989304982612354848755461,
	    .622847970537725238641159120344323,
	    .499479574071056499952214885499755,
	    .364901661346580768043989548502644,
	    .222254919776601296498260928066212,
	    .074650617461383322043914435796506 };
    double w43a[10] = { .016296734289666564924281974617663,
	    .037522876120869501461613795898115,
	    .054694902058255442147212685465005,
	    .067355414609478086075553166302174,
	    .073870199632393953432140695251367,
	    .005768556059769796184184327908655,
	    .027371890593248842081276069289151,
	    .046560826910428830743339154433824,
	    .061744995201442564496240336030883,
	    .071387267268693397768559114425516 };
    double w43b[12] = { .001844477640212414100389106552965,
	    .010798689585891651740465406741293,
	    .021895363867795428102523123075149,
	    .032597463975345689443882222526137,
	    .042163137935191811847627924327955,
	    .050741939600184577780189020092084,
	    .058379395542619248375475369330206,
	    .064746404951445885544689259517511,
	    .069566197912356484528633315038405,
	    .072824441471833208150939535192842,
	    .074507751014175118273571813842889,
	    .074722147517403005594425168280423 };
    double x4[22] = { .999902977262729234490529830591582,
	    .99798989598667874542749632236596,
	    .992175497860687222808523352251425,
	    .981358163572712773571916941623894,
	    .965057623858384619128284110607926,
	    .943167613133670596816416634507426,
	    .91580641468550720959182643072005,
	    .883221657771316501372117548744163,
	    .845710748462415666605902011504855,
	    .803557658035230982788739474980964,
	    .75700573068549555832894279343202,
	    .70627320978732181982409427474084,
	    .651589466501177922534422205016736,
	    .593223374057961088875273770349144,
	    .531493605970831932285268948562671,
	    .46676362304202284487196678165927,
	    .399424847859218804732101665817923,
	    .329874877106188288265053371824597,
	    .258503559202161551802280975429025,
	    .185695396568346652015917141167606,
	    .111842213179907468172398359241362,
	    .037352123394619870814998165437704 };

    /* System generated locals */
    double d__1, d__2, d__3, d__4;

    /* Local variables */
    int k, l;
    double fv1[5], fv2[5], fv3[5], fv4[5];
    int ipx=0;
    double absc, fval, res10, res21=0, res43=0, res87, fval1, fval2, hlgth, 
	    centr, reskh, uflow;
    
    double epmach, dhlgth, resabs=0, resasc=0, fcentr, savfun[21];

/* non-adaptive integration */
/* double precision version */

/*           f      - double precision */
/*                    function subprogram defining the integrand function */
/*                    f(x). the actual name for f needs to be declared */
/*                    e x t e r n a l in the driver program. */

/*           a      - double precision */
/*                    lower limit of integration */

/*           b      - double precision */
/*                    upper limit of integration */

/*           epsabs - double precision */
/*                    absolute accuracy requested */
/*           epsrel - double precision */
/*                    relative accuracy requested */
/*                    if  epsabs.le.0 */
/*                    and epsrel.lt.MAX(50*rel.mach.acc.,0.5d-28), */
/*                    the routine will end with ier = 6. */

/*         on return */
/*           result - double precision */
/*                    approximation to the integral i */
/*                    result is obtained by applying the 21-point */
/*                    gauss-kronrod rule (res21) obtained by optimal */
/*                    addition of abscissae to the 10-point gauss rule */
/*                    (res10), or by applying the 43-point rule (res43) */
/*                    obtained by optimal addition of abscissae to the */
/*                    21-point gauss-kronrod rule, or by applying the */
/*                    87-point rule (res87) obtained by optimal addition */
/*                    of abscissae to the 43-point rule. */

/*           abserr - double precision */
/*                    estimate of the modulus of the absolute error, */
/*                    which should equal or exceed ABS(i-result) */

/*           neval  - int */
/*                    number of integrand evaluations */

/*           ier    - ier = 0 normal and reliable terMINation of the */
/*                            routine. it is assumed that the requested */
/*                            accuracy has been achieved. */
/*                    ier.gt.0 abnormal terMINation of the routine. it is */
/*                            assumed that the requested accuracy has */
/*                            not been achieved. */
/*           error messages */
/*                    ier = 1 the MAXimum number of steps has been */
/*                            executed. the integral is probably too */
/*                            difficult to be calculated by dqng. */
/*                        = 6 the input is invalid, because */
/*                            epsabs.le.0 and */
/*                            epsrel.lt.MAX(50*rel.mach.acc.,0.5d-28). */
/*                            result, abserr and neval are set to zero. */

/* ***end prologue  dqng */



/*           the following data statements contain the */
/*           abscissae and weights of the integration rules used. */

/*           x1      abscissae common to the 10-, 21-, 43- and 87- */
/*                   point rule */
/*           x2      abscissae common to the 21-, 43- and 87-point rule */
/*           x3      abscissae common to the 43- and 87-point rule */
/*           x4      abscissae of the 87-point rule */
/*           w10     weights of the 10-point formula */
/*           w21a    weights of the 21-point formula for abscissae x1 */
/*           w21b    weights of the 21-point formula for abscissae x2 */
/*           w43a    weights of the 43-point formula for abscissae x1, x3 */
/*           w43b    weights of the 43-point formula for abscissae x3 */
/*           w87a    weights of the 87-point formula for abscissae x1, */
/*                   x2, x3 */
/*           w87b    weights of the 87-point formula for abscissae x4 */


/* gauss-kronrod-patterson quadrature coefficients for use in */
/* quadpack routine qng.  these coefficients were calculated with */
/* 101 decimal digit arithmetic by l. w. fullerton, bell labs, nov 1981. */


/*           list of major variables */
/*           ----------------------- */

/*           centr  - mid point of the integration interval */
/*           hlgth  - half-length of the integration interval */
/*           fcentr - function value at mid point */
/*           absc   - abscissa */
/*           fval   - function value */
/*           savfun - array of function values which have already been */
/*                    computed */
/*           res10  - 10-point gauss result */
/*           res21  - 21-point kronrod result */
/*           res43  - 43-point result */
/*           res87  - 87-point result */
/*           resabs - approximation to the integral of ABS(f) */
/*           resasc - approximation to the integral of ABS(f-i/(b-a)) */

/*           machine dependent constants */
/*           --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */

/* ***first executable statement  dqng */
    epmach = pnl_d1mach(4);
    uflow = pnl_d1mach(1);

/*           test on validity of parameters */
/*           ------------------------------ */

    *result = 0.;
    *abserr = 0.;
    *neval = 0;
    *ier = 6;
/* Computing MAX */
    d__1 = epmach * 50.;
    if (*epsabs <= 0 && *epsrel < MAX(d__1,5e-29))
      {
        * result = 0;
        * abserr = 0;
        * neval = 0;
        PNL_ERROR ("Wrong values of epsabs and epsrel", "pnl_dqng" );
      }
    hlgth = (*b - *a) * .5;
    dhlgth = ABS(hlgth);
    centr = (*b + *a) * .5;
    fcentr = PNL_EVAL_FUNC(f,centr);
    *neval = 21;
    *ier = 1;

/*          compute the integral using the 10- and 21-point formula. */

    for (l = 1; l <= 3; ++l) {
	switch (l) {
	    case 1:  goto L5;
	    case 2:  goto L25;
	    case 3:  goto L45;
	}
L5:
	res10 = 0.;
	res21 = w21b[5] * fcentr;
	resabs = w21b[5] * ABS(fcentr);
	for (k = 1; k <= 5; ++k) {
	    absc = hlgth * x1[k - 1];
	    d__1 = centr + absc;
	    fval1 = PNL_EVAL_FUNC(f,d__1);
	    d__1 = centr - absc;
	    fval2 = PNL_EVAL_FUNC(f,d__1);
	    fval = fval1 + fval2;
	    res10 += w10[k - 1] * fval;
	    res21 += w21a[k - 1] * fval;
	    resabs += w21a[k - 1] * (ABS(fval1) + ABS(fval2));
	    savfun[k - 1] = fval;
	    fv1[k - 1] = fval1;
	    fv2[k - 1] = fval2;
/* L10: */
	}
	ipx = 5;
	for (k = 1; k <= 5; ++k) {
	    ++ipx;
	    absc = hlgth * x2[k - 1];
	    d__1 = centr + absc;
	    fval1 = PNL_EVAL_FUNC(f,d__1);
	    d__1 = centr - absc;
	    fval2 = PNL_EVAL_FUNC(f,d__1);
	    fval = fval1 + fval2;
	    res21 += w21b[k - 1] * fval;
	    resabs += w21b[k - 1] * (ABS(fval1) + ABS(fval2));
	    savfun[ipx - 1] = fval;
	    fv3[k - 1] = fval1;
	    fv4[k - 1] = fval2;
/* L15: */
	}

/*          test for convergence. */

	*result = res21 * hlgth;
	resabs *= dhlgth;
	reskh = res21 * .5;
	resasc = w21b[5] * (d__1 = fcentr - reskh, ABS(d__1));
	for (k = 1; k <= 5; ++k) {
	    resasc = resasc + w21a[k - 1] * ((d__1 = fv1[k - 1] - reskh, ABS(
		    d__1)) + (d__2 = fv2[k - 1] - reskh, ABS(d__2))) + w21b[k 
		    - 1] * ((d__3 = fv3[k - 1] - reskh, ABS(d__3)) + (d__4 = 
		    fv4[k - 1] - reskh, ABS(d__4)));
/* L20: */
	}
	*abserr = (d__1 = (res21 - res10) * hlgth, ABS(d__1));
	resasc *= dhlgth;
	goto L65;

/*          compute the integral using the 43-point formula. */

L25:
	res43 = w43b[11] * fcentr;
	*neval = 43;
	for (k = 1; k <= 10; ++k) {
	    res43 += savfun[k - 1] * w43a[k - 1];
/* L30: */
	}
	for (k = 1; k <= 11; ++k) {
	    ++ipx;
	    absc = hlgth * x3[k - 1];
	    d__1 = absc + centr;
	    d__2 = centr - absc;
	    fval = PNL_EVAL_FUNC(f,d__1) + PNL_EVAL_FUNC(f,d__2);
	    res43 += fval * w43b[k - 1];
	    savfun[ipx - 1] = fval;
/* L40: */
	}

/*          test for convergence. */

	*result = res43 * hlgth;
	*abserr = (d__1 = (res43 - res21) * hlgth, ABS(d__1));
	goto L65;

/*          compute the integral using the 87-point formula. */

L45:
	res87 = w87b[22] * fcentr;
	*neval = 87;
	for (k = 1; k <= 21; ++k) {
	    res87 += savfun[k - 1] * w87a[k - 1];
/* L50: */
	}
	for (k = 1; k <= 22; ++k) {
	    absc = hlgth * x4[k - 1];
	    d__1 = absc + centr;
	    d__2 = centr - absc;
	    res87 += w87b[k - 1] * (PNL_EVAL_FUNC(f,d__1) + PNL_EVAL_FUNC(f,d__2));
/* L60: */
	}
	*result = res87 * hlgth;
	*abserr = (d__1 = (res87 - res43) * hlgth, ABS(d__1));
L65:
	if (resasc != 0. && *abserr != 0.) {
/* Computing MIN */
	    d__3 = *abserr * 200. / resasc;
	    d__1 = 1., d__2 = pow(d__3, 1.5);
	    *abserr = resasc * MIN(d__1,d__2);
	}
	if (resabs > uflow / (epmach * 50.)) {
/* Computing MAX */
	    d__1 = epmach * 50. * resabs;
	    *abserr = MAX(d__1,*abserr);
	}
/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * ABS(*result);
	if (*abserr <= MAX(d__1,d__2)) {
	    *ier = 0;
	}
/* ***jump out of do-loop */
	if (*ier == 0) {
	    goto L999;
	}
/* L70: */
    }
L999:
    return 0;
} /* pnl_dqng */

