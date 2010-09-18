/*
 * This file contains the routine dqk61 from the quadpack distribution
 * available on www.netlib.org
 *
 * This routine has been translated from Fortran using f2c and then slightly
 * modified by Jérôme Lelong for integration into Pnl. All Quadpack names have
 * been prefixed by pnl to avoid possible clash in names. Calls to d1mach have
 * been replaced by calls to dlamch coming from Lapack
 */

#include "pnl/pnl_internals.h"
#include "pnl/pnl_mathtools.h"

int pnl_dqk61(PnlFunc * f, double *a, double *b, double *
	result, double *abserr, double *resabs, double *resasc)
{
    /* Initialized data */

    double wg[15] = { .007968192496166605615465883474674,
	    .018466468311090959142302131912047,
	    .028784707883323369349719179611292,
	    .038799192569627049596801936446348,
	    .048402672830594052902938140422808,
	    .057493156217619066481721689402056,
	    .065974229882180495128128515115962,
	    .073755974737705206268243850022191,
	    .08075589522942021535469493846053,
	    .086899787201082979802387530715126,
	    .092122522237786128717632707087619,
	    .09636873717464425963946862635181,
	    .099593420586795267062780282103569,
	    .101762389748405504596428952168554,
	    .102852652893558840341285636705415 };
    double xgk[31] = { .999484410050490637571325895705811,
	    .996893484074649540271630050918695,
	    .991630996870404594858628366109486,
	    .983668123279747209970032581605663,
	    .973116322501126268374693868423707,
	    .960021864968307512216871025581798,
	    .944374444748559979415831324037439,
	    .926200047429274325879324277080474,
	    .905573307699907798546522558925958,
	    .882560535792052681543116462530226,
	    .857205233546061098958658510658944,
	    .829565762382768397442898119732502,
	    .799727835821839083013668942322683,
	    .767777432104826194917977340974503,
	    .733790062453226804726171131369528,
	    .69785049479331579693229238802664,
	    .660061064126626961370053668149271,
	    .620526182989242861140477556431189,
	    .57934523582636169175602493217254,
	    .536624148142019899264169793311073,
	    .492480467861778574993693061207709,
	    .447033769538089176780609900322854,
	    .400401254830394392535476211542661,
	    .352704725530878113471037207089374,
	    .304073202273625077372677107199257,
	    .254636926167889846439805129817805,
	    .204525116682309891438957671002025,
	    .153869913608583546963794672743256,
	    .102806937966737030147096751318001,
	    .051471842555317695833025213166723,0. };
    double wgk[31] = { .00138901369867700762455159122676,
	    .003890461127099884051267201844516,
	    .00663070391593129217331982636975,
	    .009273279659517763428441146892024,
	    .011823015253496341742232898853251,
	    .01436972950704580481245143244358,
	    .016920889189053272627572289420322,
	    .019414141193942381173408951050128,
	    .021828035821609192297167485738339,
	    .024191162078080601365686370725232,
	    .026509954882333101610601709335075,
	    .028754048765041292843978785354334,
	    .030907257562387762472884252943092,
	    .032981447057483726031814191016854,
	    .034979338028060024137499670731468,
	    .036882364651821229223911065617136,
	    .038678945624727592950348651532281,
	    .040374538951535959111995279752468,
	    .04196981021516424614714754128597,
	    .043452539701356069316831728117073,
	    .044814800133162663192355551616723,
	    .046059238271006988116271735559374,
	    .047185546569299153945261478181099,
	    .048185861757087129140779492298305,
	    .049055434555029778887528165367238,
	    .049795683427074206357811569379942,
	    .050405921402782346840893085653585,
	    .050881795898749606492297473049805,
	    .051221547849258772170656282604944,
	    .051426128537459025933862879215781,
	    .051494729429451567558340433647099 };

    double c_b7 = 1.5;

    /* System generated locals */
    double d__1, d__2, d__3;

    /* Local variables */
    int j;
    double fc, fv1[30], fv2[30];
    int jtw;
    double resg, resk, fsum, fval1, fval2;
    int jtwm1;
    double dabsc, hlgth, centr, reskh, uflow;
    
    double epmach, dhlgth;

/* ***begin prologue  dqk61 */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a1a2 */
/* ***keywords  61-point gauss-kronrod rules */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  to compute i = integral of f over (a,b) with error */
/*                           estimate */
/*                       j = integral of dABS(f) over (a,b) */
/* ***description */

/*        integration rule */
/*        standard fortran subroutine */
/*        double precision version */


/*        parameters */
/*         on entry */
/*           f      - double precision */
/*                    function subprogram defining the integrand */
/*                    function f(x). the actual name for f needs to be */
/*                    declared e x t e r n a l in the calling program. */

/*           a      - double precision */
/*                    lower limit of integration */

/*           b      - double precision */
/*                    upper limit of integration */

/*         on return */
/*           result - double precision */
/*                    approximation to the integral i */
/*                    result is computed by applying the 61-point */
/*                    kronrod rule (resk) obtained by optimal addition of */
/*                    abscissae to the 30-point gauss rule (resg). */

/*           abserr - double precision */
/*                    estimate of the modulus of the absolute error, */
/*                    which should equal or exceed dABS(i-result) */

/*           resabs - double precision */
/*                    approximation to the integral j */

/*           resasc - double precision */
/*                    approximation to the integral of dABS(f-i/(b-a)) */


/* ***references  (none) */
/* ***end prologue  dqk61 */



/*           the abscissae and weights are given for the */
/*           interval (-1,1). because of symmetry only the positive */
/*           abscissae and their corresponding weights are given. */

/*           xgk   - abscissae of the 61-point kronrod rule */
/*                   xgk(2), xgk(4)  ... abscissae of the 30-point */
/*                   gauss rule */
/*                   xgk(1), xgk(3)  ... optimally added abscissae */
/*                   to the 30-point gauss rule */

/*           wgk   - weights of the 61-point kronrod rule */

/*           wg    - weigths of the 30-point gauss rule */


/* gauss quadrature weights and kronron quadrature abscissae and weights */
/* as evaluated with 80 decimal digit arithmetic by l. w. fullerton, */
/* bell labs, nov. 1981. */




/*           list of major variables */
/*           ----------------------- */

/*           centr  - mid point of the interval */
/*           hlgth  - half-length of the interval */
/*           dabsc  - abscissa */
/*           fval*  - function value */
/*           resg   - result of the 30-point gauss rule */
/*           resk   - result of the 61-point kronrod rule */
/*           reskh  - approximation to the mean value of f */
/*                    over (a,b), i.e. to i/(b-a) */

/*           machine dependent constants */
/*           --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */

    epmach = pnl_d1mach(4);
    uflow = pnl_d1mach(1);

    centr = (*b + *a) * .5;
    hlgth = (*b - *a) * .5;
    dhlgth = ABS(hlgth);

/*           compute the 61-point kronrod approximation to the */
/*           integral, and estimate the absolute error. */

/* ***first executable statement  dqk61 */
    resg = 0.;
    fc = PNL_EVAL_FUNC(f,centr);
    resk = wgk[30] * fc;
    *resabs = ABS(resk);
    for (j = 1; j <= 15; ++j) {
	jtw = j << 1;
	dabsc = hlgth * xgk[jtw - 1];
	d__1 = centr - dabsc;
	fval1 = PNL_EVAL_FUNC(f,d__1);
	d__1 = centr + dabsc;
	fval2 = PNL_EVAL_FUNC(f,d__1);
	fv1[jtw - 1] = fval1;
	fv2[jtw - 1] = fval2;
	fsum = fval1 + fval2;
	resg += wg[j - 1] * fsum;
	resk += wgk[jtw - 1] * fsum;
	*resabs += wgk[jtw - 1] * (ABS(fval1) + ABS(fval2));
/* L10: */
    }
    for (j = 1; j <= 15; ++j) {
	jtwm1 = (j << 1) - 1;
	dabsc = hlgth * xgk[jtwm1 - 1];
	d__1 = centr - dabsc;
	fval1 = PNL_EVAL_FUNC(f,d__1);
	d__1 = centr + dabsc;
	fval2 = PNL_EVAL_FUNC(f,d__1);
	fv1[jtwm1 - 1] = fval1;
	fv2[jtwm1 - 1] = fval2;
	fsum = fval1 + fval2;
	resk += wgk[jtwm1 - 1] * fsum;
	*resabs += wgk[jtwm1 - 1] * (ABS(fval1) + ABS(fval2));
/* L15: */
    }
    reskh = resk * .5;
    *resasc = wgk[30] * (d__1 = fc - reskh, ABS(d__1));
    for (j = 1; j <= 30; ++j) {
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
} /* pnl_dqk61 */

