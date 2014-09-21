/**
Computation of Faddeeva's complex scaled error function,
  w(z) = exp(-z^2) * erfc(-i*z),
nameless function (7.1.3) of Abramowitz&Stegun (1964),
also known as the plasma dispersion function.

This implementation uses a combination of different algorithms.
 
Copyright:
    © 2012 Massachusetts Institute of Technology,    
    © 2013 Forschungszentrum Jülich GmbH,
    © 2014 Ilya Yaroshenko
 
Authors:
    Steven G. Johnson, core author;
    Joachim Wuttke, C package maintainer;
    $(LINK2 http://plus.google.com/+IlyaYaroshenko, Ilya Yaroshenko), D package maintainer

License: 
    Subject to the terms of the MIT license, as written in the included LICENSE.txt file.
*/
module libcerf.w_of_z;

import std.complex;
import core.stdc.math : exp, fabs, sin, cos;
import std.math : isNaN, isInfinity, floor;

import libcerf.erfcx_;
import libcerf.im_w_of_x;

version(unittest)
{
    import libcerf.testutils;
}

///
//import core.stdc.complex : cexp;
extern (C) @trusted nothrow @nogc Complex!double cexp(Complex!double z);

/**
 w_of_z, Faddeeva's scaled complex error function                         

   Computes various error functions (erf, erfc, erfi, erfcx), 
   including the Dawson integral, in the complex plane, based
   on algorithms for the computation of the Faddeeva function 
              w(z) = exp(-z^2) * erfc(-i*z).
   Given w(z), the error functions are mostly straightforward
   to compute, except for certain regions where we have to
   switch to Taylor expansions to avoid cancellation errors
   [e.g. near the origin for erf(z)].

*/
//Steven G. Johnson, October 2012.
Complex!double w_of_z(Complex!double z) @safe @nogc nothrow
{
    if (!z.re) {
        // Purely imaginary input, purely real output.
        // However, use z.re to give correct sign of 0 in cimag(w).
        return Complex!double(erfcx(z.im), z.re);
    }
    if (!z.im) {
        // Purely real input, complex output.
        return Complex!double(exp(-sqr(z.re)),  im_w_of_x(z.re));
    }
    enum double a = 0.518321480430085929872; // pi / sqrt(-log(eps*0.5))
    enum double c = 0.329973702884629072537; // (2/pi) * a;
    enum double a2 = 0.268657157075235951582; // a^2

    immutable x = fabs(z.re);
    immutable y = z.im;
    immutable ya = fabs(y);

    Complex!double ret = 0.; // return value

    double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0;

    if (ya > 7 || (x > 6  // continued fraction is faster
                   /* As pointed out by M. Zaghloul, the continued
                      fraction seems to give a large relative error in
                      Re w(z) for |x| ~ 6 and small |y|, so use
                      algorithm 816 in this region: */
                   && (ya > 0.1 || (x > 8 && ya > 1e-10) || x > 28))) {

        /* Poppe & Wijers suggest using a number of terms
           nu = 3 + 1442 / (26*rho + 77)
           where rho = sqrt((x/x0)^2 + (y/y0)^2) where x0=6.3, y0=4.4.
           (They only use this expansion for rho >= 1, but rho a little less
           than 1 seems okay too.)
           Instead, I did my own fit to a slightly different function
           that avoids the hypotenuse calculation, using NLopt to minimize
           the sum of the squares of the errors in nu with the constraint
           that the estimated nu be >= minimum nu to attain machine precision.
           I also separate the regions where nu == 2 and nu == 1. */
        enum double ispi = M_2_SQRTPI / 2; // 1 / sqrt(pi)
        immutable xs = y < 0 ? -z.re : z.re; // compute for -z if y < 0
        if (x + ya > 4000) { // nu <= 2
            if (x + ya > 1e7) { // nu == 1, w(z) = i/sqrt(pi) / z
                // scale to avoid overflow
                if (x > ya) {
                    immutable double yax = ya / xs; 
                    immutable double denom = ispi / (xs + yax*ya);
                    ret = Complex!double(denom*yax, denom);
                }
                else if (isInfinity(ya)) {
                    return ((isNaN(x) || y < 0) 
                            ? Complex!double.init : Complex!double(0,0));
                }
                else {
                    immutable double xya = xs / ya;
                    immutable double denom = ispi / (xya*xs + ya);
                    ret = Complex!double(denom, denom*xya);
                }
            }
            else { // nu == 2, w(z) = i/sqrt(pi) * z / (z*z - 0.5)
                immutable double dr = xs*xs - ya*ya - 0.5, di = 2*xs*ya;
                immutable double denom = ispi / (dr*dr + di*di);
                ret = Complex!double(denom * (xs*di-ya*dr), denom * (xs*dr+ya*di));
            }
        }
        else { // compute nu(z) estimate and do general continued fraction
            immutable double c0=3.9, c1=11.398, c2=0.08254, c3=0.1421, c4=0.2023; // fit
            double nu = floor(c0 + c1 / (c2*x + c3*ya + c4));
            double wr = xs, wi = ya;
            for (nu = 0.5 * (nu - 1); nu > 0.4; nu -= 0.5) {
                // w <- z - nu/w:
                double denom = nu / (wr*wr + wi*wi);
                wr = xs - wr * denom;
                wi = ya + wi * denom;
            }
            { // w(z) = i/sqrt(pi) / w:
                immutable double denom = ispi / (wr*wr + wi*wi);
                ret = Complex!double(denom*wi, denom*wr);
            }
        }
        if (y < 0) {
            // use w(z) = 2.0*exp(-z*z) - w(-z), 
            // but be careful of overflow in exp(-z*z) 
            //                                = exp(-(xs*xs-ya*ya) -2*i*xs*ya) 
            return 2.0 * cexp(Complex!double((ya-xs)*(xs+ya), 2*xs*y)) - ret;
        }
        else
            return ret;
    }

    /* Note: The test that seems to be suggested in the paper is x <
       sqrt(-log(DBL_MIN)), about 26.6, since otherwise exp(-x^2)
       underflows to zero and sum1,sum2,sum4 are zero.  However, long
       before this occurs, the sum1,sum2,sum4 contributions are
       negligible in double precision; I find that this happens for x >
       about 6, for all y.  On the other hand, I find that the case
       where we compute all of the sums is faster (at least with the
       precomputed expa2n2 table) until about x=10.  Furthermore, if we
       try to compute all of the sums for x > 20, I find that we
       sometimes run into numerical problems because underflow/overflow
       problems start to appear in the various coefficients of the sums,
       below.  Therefore, we use x < 10 here. */
    else if (x < 10) {

        double prod2ax = 1, prodm2ax = 1;
        double expx2;

        if (isNaN(y)) {
            return Complex!double(y,y);
        }

        if (x < 5e-4) { // compute sum4 and sum5 together as sum5-sum4
                        // This special case is needed for accuracy.
            immutable double x2 = x*x;
            expx2 = 1 - x2 * (1 - 0.5*x2); // exp(-x*x) via Taylor
            // compute exp(2*a*x) and exp(-2*a*x) via Taylor, to double precision
            immutable double ax2 = 1.036642960860171859744*x; // 2*a*x
            immutable double exp2ax =
                1 + ax2 * (1 + ax2 * (0.5 + 0.166666666666666666667*ax2));
            immutable double expm2ax =
                1 - ax2 * (1 - ax2 * (0.5 - 0.166666666666666666667*ax2));
            for (uint n = 1; ; n++) {
                immutable double coef = expa2n2[n-1] * expx2 / (a2*(n*n) + y*y);
                prod2ax *= exp2ax;
                prodm2ax *= expm2ax;
                sum1 += coef;
                sum2 += coef * prodm2ax;
                sum3 += coef * prod2ax;
          
                // really = sum5 - sum4
                sum5 += coef * (2*a) * n * sinh_taylor((2*a)*n*x);
          
                // test convergence via sum3
                if (coef * prod2ax < double.epsilon * sum3) 
                    break;
            }
        }
        else { // x > 5e-4, compute sum4 and sum5 separately
            expx2 = exp(-x*x);
            immutable double exp2ax = exp((2*a)*x), expm2ax = 1 / exp2ax;
            for (int n = 1; 1; ++n) {
                immutable double coef = expa2n2[n-1] * expx2 / (a2*(n*n) + y*y);
                prod2ax *= exp2ax;
                prodm2ax *= expm2ax;
                sum1 += coef;
                sum2 += coef * prodm2ax;
                sum4 += (coef * prodm2ax) * (a*n);
                sum3 += coef * prod2ax;
                sum5 += (coef * prod2ax) * (a*n);
                // test convergence via sum5, since this sum has the slowest decay
                if ((coef * prod2ax) * (a*n) < double.epsilon * sum5) 
                    break;
            }
        }
        immutable double expx2erfcxy = // avoid spurious overflow for large negative y
            y > -6 // for y < -6, erfcx(y) = 2*exp(y*y) to double precision
            ? expx2*erfcx(y) : 2*exp(y*y-x*x);
        if (y > 5) { // imaginary terms cancel
            immutable double sinxy = sin(x*y);
            ret = (expx2erfcxy - c*y*sum1) * cos(2*x*y)
                + (c*x*expx2) * sinxy * sinc(x*y, sinxy);
        }
        else {
            double xs = z.re;
            immutable double sinxy = sin(xs*y);
            immutable double sin2xy = sin(2*xs*y), cos2xy = cos(2*xs*y);
            immutable double coef1 = expx2erfcxy - c*y*sum1;
            immutable double coef2 = c*xs*expx2;
            ret = Complex!double(coef1 * cos2xy + coef2 * sinxy * sinc(xs*y, sinxy),
                    coef2 * sinc(2*xs*y, sin2xy) - coef1 * sin2xy);
        }
    }
    else { // x large: only sum3 & sum5 contribute (see above note)    


        if (isNaN(x))
            return Complex!double(x,x);
        if (isNaN(y))
            return Complex!double(y,y);
        ret = exp(-x*x); // |y| < 1e-10, so we only need exp(-x*x) term
        // (round instead of ceil as in original paper; note that x/a > 1 here)
        immutable double n0 = floor(x/a + 0.5); // sum in both directions, starting at n0
        immutable double dx = a*n0 - x;
        sum3 = exp(-dx*dx) / (a2*(n0*n0) + y*y);
        sum5 = a*n0 * sum3;
        immutable double exp1 = exp(4*a*dx);
        double exp1dn = 1;
        int dn;
        for (dn = 1; n0 - dn > 0; ++dn) { // loop over n0-dn and n0+dn terms
            immutable double np = n0 + dn, nm = n0 - dn;
            double tp = exp(-sqr(a*dn+dx));
            double tm = tp * (exp1dn *= exp1); // trick to get tm from tp
            tp /= (a2*(np*np) + y*y);
            tm /= (a2*(nm*nm) + y*y);
            sum3 += tp + tm;
            sum5 += a * (np * tp + nm * tm);
            if (a * (np * tp + nm * tm) < double.epsilon * sum5) 
                goto finish;
        }
        while (1) { // loop over n0+dn terms only (since n0-dn <= 0)
            immutable double np = n0 + dn++;
            immutable double tp = exp(-sqr(a*dn+dx)) / (a2*(np*np) + y*y);
            sum3 += tp;
            sum5 += a * np * tp;
            if (a * np * tp < double.epsilon * sum5) 
                goto finish;
        }
    }
finish:
    return ret + Complex!double((0.5*c)*y*(sum2+sum3), 
                   (0.5*c)*copysign(sum5-sum4, z.re));
}

private:

/**
auxiliary function                                            
*/
double sinc(double x, double sinx) @safe @nogc nothrow
{ 
    // return sinc(x) = sin(x)/x, given both x and sin(x) 
    // [since we only use this in cases where sin(x) has already been computed]
    return fabs(x) < 1e-4 ? 1 - 0.1666666666666666666667*x*x : sinx / x; 
}

///ditto
double sinh_taylor(double x) pure @safe @nogc nothrow
{
    // sinh(x) via Taylor series, accurate to machine precision for |x| < 1e-2
    immutable x2 = x*x;
    return x * (1 + (x2) * (0.1666666666666666666667
                             + 0.00833333333333333333333 * (x2)));
}

///ditto
double sqr(double x) pure @safe @nogc nothrow { return x*x; } 

/**
precomputed table of expa2n2[n-1] = exp(-a2*n*n)                           
for double-precision a2 = 0.26865... in w_of_z, below.                     
*/
immutable double[] expa2n2 = [
    7.64405281671221563e-01,
    3.41424527166548425e-01,
    8.91072646929412548e-02,
    1.35887299055460086e-02,
    1.21085455253437481e-03,
    6.30452613933449404e-05,
    1.91805156577114683e-06,
    3.40969447714832381e-08,
    3.54175089099469393e-10,
    2.14965079583260682e-12,
    7.62368911833724354e-15,
    1.57982797110681093e-17,
    1.91294189103582677e-20,
    1.35344656764205340e-23,
    5.59535712428588720e-27,
    1.35164257972401769e-30,
    1.90784582843501167e-34,
    1.57351920291442930e-38,
    7.58312432328032845e-43,
    2.13536275438697082e-47,
    3.51352063787195769e-52,
    3.37800830266396920e-57,
    1.89769439468301000e-62,
    6.22929926072668851e-68,
    1.19481172006938722e-73,
    1.33908181133005953e-79,
    8.76924303483223939e-86,
    3.35555576166254986e-92,
    7.50264110688173024e-99,
    9.80192200745410268e-106,
    7.48265412822268959e-113,
    3.33770122566809425e-120,
    8.69934598159861140e-128,
    1.32486951484088852e-135,
    1.17898144201315253e-143,
    6.13039120236180012e-152,
    1.86258785950822098e-160,
    3.30668408201432783e-169,
    3.43017280887946235e-178,
    2.07915397775808219e-187,
    7.36384545323984966e-197,
    1.52394760394085741e-206,
    1.84281935046532100e-216,
    1.30209553802992923e-226,
    5.37588903521080531e-237,
    1.29689584599763145e-247,
    1.82813078022866562e-258,
    1.50576355348684241e-269,
    7.24692320799294194e-281,
    2.03797051314726829e-292,
    3.34880215927873807e-304,
    0.0 // underflow (also prevents reads past array end, below)
]; 

unittest 
{
    static immutable Z = [
        Complex!double(624.2,-0.26123),
        Complex!double(-0.4,3.),
        Complex!double(0.6,2.),
        Complex!double(-1.,1.),
        Complex!double(-1.,-9.),
        Complex!double(-1.,9.),
        Complex!double(-0.0000000234545,1.1234),
        Complex!double(-3.,5.1),
        Complex!double(-53,30.1),
        Complex!double(0.0,0.12345),
        Complex!double(11,1),
        Complex!double(-22,-2),
        Complex!double(9,-28),
        Complex!double(21,-33),
        Complex!double(1e5,1e5),
        Complex!double(1e14,1e14),
        Complex!double(-3001,-1000),
        Complex!double(1e160,-1e159),
        Complex!double(-6.01,0.01),
        Complex!double(-0.7,-0.7),
        Complex!double(2.611780000000000e+01, 4.540909610972489e+03),
        Complex!double(0.8e7,0.3e7),
        Complex!double(-20,-19.8081),
        Complex!double(1e-16,-1.1e-16),
        Complex!double(2.3e-8,1.3e-8),
        Complex!double(6.3,-1e-13),
        Complex!double(6.3,1e-20),
        Complex!double(1e-20,6.3),
        Complex!double(1e-20,16.3),
        Complex!double(9,1e-300),
        Complex!double(6.01,0.11),
        Complex!double(8.01,1.01e-10),
        Complex!double(28.01,1e-300),
        Complex!double(10.01,1e-200),
        Complex!double(10.01,-1e-200),
        Complex!double(10.01,0.99e-10),
        Complex!double(10.01,-0.99e-10),
        Complex!double(1e-20,7.01),
        Complex!double(-1,7.01),
        Complex!double(5.99,7.01),
        Complex!double(1,0),
        Complex!double(55,0),
        Complex!double(-0.1,0),
        Complex!double(1e-20,0),
        Complex!double(0,5e-14),
        Complex!double(0,51),
        Complex!double(double.infinity,0),
        Complex!double(-double.infinity,0),
        Complex!double(0,double.infinity),
        Complex!double(0,-double.infinity),
        Complex!double(double.infinity,double.infinity),
        Complex!double(double.infinity,-double.infinity),
        Complex!double(double.nan,double.nan),
        Complex!double(double.nan,0),
        Complex!double(0,double.nan),
        Complex!double(double.nan,double.infinity),
        Complex!double(double.infinity,double.nan)
    ];
    /* w(z), computed with WolframAlpha
    ... note that WolframAlpha is problematic
    some of the above inputs, so I had to
    use the continued-fraction expansion
    in WolframAlpha in some cases, or switch
    to Maple */
    static immutable W = [ 
        Complex!double(-3.78270245518980507452677445620103199303131110e-7,
          0.000903861276433172057331093754199933411710053155),
        Complex!double(0.1764906227004816847297495349730234591778719532788,
          -0.02146550539468457616788719893991501311573031095617),
        Complex!double(0.2410250715772692146133539023007113781272362309451,
          0.06087579663428089745895459735240964093522265589350),
        Complex!double(0.30474420525691259245713884106959496013413834051768,
          -0.20821893820283162728743734725471561394145872072738),
        Complex!double(7.317131068972378096865595229600561710140617977e34,
          8.321873499714402777186848353320412813066170427e34),
        Complex!double(0.0615698507236323685519612934241429530190806818395,
          -0.00676005783716575013073036218018565206070072304635),
        Complex!double(0.3960793007699874918961319170187598400134746631,
          -5.593152259116644920546186222529802777409274656e-9),
        Complex!double(0.08217199226739447943295069917990417630675021771804,
          -0.04701291087643609891018366143118110965272615832184),
        Complex!double(0.00457246000350281640952328010227885008541748668738,
          -0.00804900791411691821818731763401840373998654987934),
        Complex!double(0.8746342859608052666092782112565360755791467973338452,
          0.),
        Complex!double(0.00468190164965444174367477874864366058339647648741,
          0.0510735563901306197993676329845149741675029197050),
        Complex!double(-0.0023193175200187620902125853834909543869428763219,
          -0.025460054739731556004902057663500272721780776336),
        Complex!double(9.11463368405637174660562096516414499772662584e304,
          3.97101807145263333769664875189354358563218932e305),
        Complex!double(-4.4927207857715598976165541011143706155432296e281,
          -2.8019591213423077494444700357168707775769028e281),
        Complex!double(2.820947917809305132678577516325951485807107151e-6,
          2.820947917668257736791638444590253942253354058e-6),
        Complex!double(2.82094791773878143474039725787438662716372268e-15,
          2.82094791773878143474039725773333923127678361e-15),
        Complex!double(-0.0000563851289696244350147899376081488003110150498,
          -0.000169211755126812174631861529808288295454992688),
        Complex!double(-5.586035480670854326218608431294778077663867e-162,
          5.586035480670854326218608431294778077663867e-161),
        Complex!double(0.00016318325137140451888255634399123461580248456,
          -0.095232456573009287370728788146686162555021209999),
        Complex!double(0.69504753678406939989115375989939096800793577783885,
          -1.8916411171103639136680830887017670616339912024317),
        Complex!double(0.0001242418269653279656612334210746733213167234822,
          7.145975826320186888508563111992099992116786763e-7),
        Complex!double(2.318587329648353318615800865959225429377529825e-8,
          6.182899545728857485721417893323317843200933380e-8),
        Complex!double(-0.0133426877243506022053521927604277115767311800303,
          -0.0148087097143220769493341484176979826888871576145),
        Complex!double(1.00000000000000012412170838050638522857747934,
          1.12837916709551279389615890312156495593616433e-16),
        Complex!double(0.9999999853310704677583504063775310832036830015,
          2.595272024519678881897196435157270184030360773e-8),
        Complex!double(-1.4731421795638279504242963027196663601154624e-15,
          0.090727659684127365236479098488823462473074709),
        Complex!double(5.79246077884410284575834156425396800754409308e-18,
          0.0907276596841273652364790985059772809093822374),
        Complex!double(0.0884658993528521953466533278764830881245144368,
          1.37088352495749125283269718778582613192166760e-22),
        Complex!double(0.0345480845419190424370085249304184266813447878,
          2.11161102895179044968099038990446187626075258e-23),
        Complex!double(6.63967719958073440070225527042829242391918213e-36,
          0.0630820900592582863713653132559743161572639353),
        Complex!double(0.00179435233208702644891092397579091030658500743634,
          0.0951983814805270647939647438459699953990788064762),
        Complex!double(9.09760377102097999924241322094863528771095448e-13,
          0.0709979210725138550986782242355007611074966717),
        Complex!double(7.2049510279742166460047102593255688682910274423e-304,
          0.0201552956479526953866611812593266285000876784321),
        Complex!double(3.04543604652250734193622967873276113872279682e-44,
          0.0566481651760675042930042117726713294607499165),
        Complex!double(3.04543604652250734193622967873276113872279682e-44,
          0.0566481651760675042930042117726713294607499165),
        Complex!double(0.5659928732065273429286988428080855057102069081e-12,
          0.056648165176067504292998527162143030538756683302),
        Complex!double(-0.56599287320652734292869884280802459698927645e-12,
          0.0566481651760675042929985271621430305387566833029),
        Complex!double(0.0796884251721652215687859778119964009569455462,
          1.11474461817561675017794941973556302717225126e-22),
        Complex!double(0.07817195821247357458545539935996687005781943386550,
          -0.01093913670103576690766705513142246633056714279654),
        Complex!double(0.04670032980990449912809326141164730850466208439937,
          0.03944038961933534137558064191650437353429669886545),
        Complex!double(0.36787944117144232159552377016146086744581113103176,
          0.60715770584139372911503823580074492116122092866515),
        Complex!double(0,
          0.010259688805536830986089913987516716056946786526145),
        Complex!double(0.99004983374916805357390597718003655777207908125383,
          -0.11208866436449538036721343053869621153527769495574),
        Complex!double(0.99999999999999999999999999999999999999990000,
          1.12837916709551257389615890312154517168802603e-20),
        Complex!double(0.999999999999943581041645226871305192054749891144158,
          0),
        Complex!double(0.0110604154853277201542582159216317923453996211744250,
          0),
        Complex!double(0,0),
        Complex!double(0,0),
        Complex!double(0,0),
        Complex!double(double.infinity,0),
        Complex!double(0,0),
        Complex!double(double.nan,double.nan),
        Complex!double(double.nan,double.nan),
        Complex!double(double.nan,double.nan),
        Complex!double(double.nan,0),
        Complex!double(double.nan,double.nan),
        Complex!double(double.nan,double.nan)
    ];
    commonTest!w_of_z(Z, W);
}
