/**
Computate Dawson, Voigt, and several error functions,
based on erfcx, im_w_of_x, w_of_z as implemented in separate files.

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
module libcerf.err_fcts;

import std.complex;
import std.math : isNaN, copysign, SQRT2, M_1_PI;
import core.stdc.math : sin, cos, erf;

import libcerf.w_of_z;
import libcerf.im_w_of_x;
import libcerf.erfcx_;

version(unittest)
{
    import libcerf.testutils;
}

private enum double spi2 = 0.8862269254527580136490837416705725913990; // sqrt(pi)/2
private enum double s2pi = 2.5066282746310005024157652848110; // sqrt(2*pi)

/**
Compute erfcx(z) = exp(z^2) erfc(z),
the complex underflow-compensated complementary error function,
trivially related to Faddeeva's w_of_z.
*/
Complex!double cerfcx(Complex!double z) @safe @nogc nothrow
{
    return w_of_z(Complex!double(-z.im, z.re));
}

unittest
{
    static immutable Z = [ Complex!double(1.234,0.5678) ];
    static immutable W = [ 
        Complex!double(0.3382187479799972294747793561190487832579,
          -0.1116077470811648467464927471872945833154)
    ];
    commonTest!cerfcx(Z, W);
}

/**
Compute erfi(z) = -i erf(iz), the rotated complex error function.
*/
Complex!double cerfi(Complex!double z) @safe @nogc nothrow
{
    Complex!double e = cerf(Complex!double(-z.im, z.re));
    return Complex!double(e.im, -e.re);
}

unittest
{
    static immutable Z = [ Complex!double(1.234,0.5678) ];
    static immutable W = [ 
        Complex!double(1.081032284405373149432716643834106923212,
          1.926775520840916645838949402886591180834)
    ];
    commonTest!cerfi(Z, W);
}

/**
Compute erfi(x) = -i erf(ix), the imaginary error function.
*/
double erfi(double x) @safe @nogc nothrow
{
    immutable x2 = x*x;
    return x2 > 720 ? copysign(double.infinity, x) : exp(x2) * im_w_of_x(x);
}

unittest {
    double[2][] tests = [
        [0.0, 0],
        [1.0, 1.65042575879754287602533772956136244389],
        [-1.0, -1.65042575879754287602533772956136244389],
        [5.0, 8.29827388067680351614622319074691995187273573754139844564e+9],
        [-5.0, -8.29827388067680351614622319074691995187273573754139844564e+9],
    ];
    foreach(test; tests)
    {
        auto x = test[0];
        auto t = erfi(x);
        auto y = test[1];
        assert(relativeErrorCheck(t, y));
    }
}

/**
Compute dawson(x) = sqrt(pi)/2 * exp(-x^2) * erfi(x),
Dawson's integral for a real argument.
*/
double dawson(double x) @safe @nogc nothrow
{


    return spi2 * im_w_of_x(x);
}

/**
Compute Voigt's convolution of a Gaussian
G(x,sigma) = 1/sqrt(2*pi)/|sigma| * exp(-x^2/2/sigma^2)
and a Lorentzian
L(x,gamma) = |gamma| / pi / ( x^2 + gamma^2 ),
namely
voigt(x,sigma,gamma) =
      \int_{-infty}^{infty} dx' G(x',sigma) L(x-x',gamma)
using the relation
voigt(x,sigma,gamma) = Re{ w(z) } / sqrt(2*pi) / |sigma|
with
z = (x+i*|gamma|) / sqrt(2) / |sigma|.

Reference: Abramowitz&Stegun (1964), formula (7.4.13).
*/
//Joachim Wuttke, January 2013.
double voigt( double x, double sigma, double gamma ) @safe @nogc nothrow
{

    immutable gam = fabs(gamma);
    immutable sig = fabs(sigma);

    if (!gam) {
        if (!sig) {
            // It's kind of a delta function
            return x ? 0 : double.infinity;
        } else {
            // It's a pure Gaussian
            return exp( -x*x/2/(sig*sig) ) / s2pi / sig;
        }
    } else {
        if (!sig) {
            // It's a pure Lorentzian
            return gam * M_1_PI / (x*x + gam*gam);
        } else {
            // Regular case, both parameters are nonzero
            Complex!double z = Complex!double(x,gam) / double(SQRT2) / sig;
            return  w_of_z(z).re / s2pi / sig;
        }
    }
}

/**
Compute erf(z), the complex error function,
using w_of_z except for certain regions.
*/
//Steven G. Johnson, October 2012.
Complex!double cerf(Complex!double z) @safe @nogc nothrow
{
    immutable x = z.re, y = z.im;

    if (!y)
        return Complex!double(erf(x), y); // preserve sign of 0
    if (!x) // handle separately for speed & handling of y = Inf or NaN
        return Complex!double(x, // preserve sign of 0
                 /* handle y -> Inf limit manually, since
                    exp(y^2) -> Inf but Im[w(y)] -> 0, so
                    IEEE will give us a NaN when it should be Inf */
                 y*y > 720 ? copysign(double.infinity, y)
                 : exp(y*y) * im_w_of_x(y));
  
    immutable mRe_z2 = (y - x) * (x + y); // Re(-z^2), being careful of overflow
    immutable mIm_z2 = -2*x*y; // Im(-z^2)
    if (mRe_z2 < -750) // underflow
        return Complex!double(copysign(1.0, x), 0.0);

    /* Handle positive and negative x via different formulas,
       using the mirror symmetries of w, to avoid overflow/underflow
       problems from multiplying exponentially large and small quantities. */
    if (x >= 0) {
        if (x < 8e-2) {
            if (fabs(y) < 1e-2)
                goto taylor;
            else if (fabs(mIm_z2) < 5e-3 && x < 5e-3)
                goto taylor_erfi;
        }
        /* don't use complex exp function, since that will produce spurious NaN
           values when multiplying w in an overflow situation. */
        return 1.0 - exp(mRe_z2) *
            (Complex!double(cos(mIm_z2), sin(mIm_z2))
             * w_of_z(Complex!double(-y,x)));
    }
    else { // x < 0
        if (x > -8e-2) { // duplicate from above to avoid fabs(x) call
            if (fabs(y) < 1e-2)
                goto taylor;
            else if (fabs(mIm_z2) < 5e-3 && x > -5e-3)
                goto taylor_erfi;
        }
        else if (isNaN(x))
            return Complex!double(double.nan, y == 0 ? 0 : double.nan);
        /* don't use complex exp function, since that will produce spurious NaN
           values when multiplying w in an overflow situation. */
        return exp(mRe_z2) *
            (Complex!double(cos(mIm_z2), sin(mIm_z2))
             * w_of_z(Complex!double(y,-x))) - 1.0;
    }

    // Use Taylor series for small |z|, to avoid cancellation inaccuracy
    //   erf(z) = 2/sqrt(pi) * z * (1 - z^2/3 + z^4/10 - z^6/42 + z^8/216 + ...)
taylor:
    Complex!double mz2 = Complex!double(mRe_z2, mIm_z2); // -z^2
    return z * (1.1283791670955125739
                + mz2 * (0.37612638903183752464
                         + mz2 * (0.11283791670955125739
                                  + mz2 * (0.026866170645131251760
                                           + mz2 * 0.0052239776254421878422))));

    /* for small |x| and small |xy|, 
       use Taylor series to avoid cancellation inaccuracy:
       erf(x+iy) = erf(iy)
       + 2*exp(y^2)/sqrt(pi) *
       [ x * (1 - x^2 * (1+2y^2)/3 + x^4 * (3+12y^2+4y^4)/30 + ... 
       - i * x^2 * y * (1 - x^2 * (3+2y^2)/6 + ...) ]
       where:
       erf(iy) = exp(y^2) * Im[w(y)]
    */
taylor_erfi:
    double x2 = x*x, y2 = y*y;
    double expy2 = exp(y2);
    return Complex!double
        (expy2 * x * (1.1283791670955125739
                      - x2 * (0.37612638903183752464
                              + 0.75225277806367504925*y2)
                      + x2*x2 * (0.11283791670955125739
                                 + y2 * (0.45135166683820502956
                                         + 0.15045055561273500986*y2))),
         expy2 * (im_w_of_x(y)
                  - x2*y * (1.1283791670955125739 
                            - x2 * (0.56418958354775628695 
                                    + 0.37612638903183752464*y2))));
}

unittest
{
    static immutable Z = [
        Complex!double(1,2),
        Complex!double(-1,2),
        Complex!double(1,-2),
        Complex!double(-1,-2),
        Complex!double(9,-28),
        Complex!double(21,-33),
        Complex!double(1e3,1e3),
        Complex!double(-3001,-1000),
        Complex!double(1e160,-1e159),
        Complex!double(5.1e-3, 1e-8),
        Complex!double(-4.9e-3, 4.95e-3),
        Complex!double(4.9e-3, 0.5),
        Complex!double(4.9e-4, -0.5e1),
        Complex!double(-4.9e-5, -0.5e2),
        Complex!double(5.1e-3, 0.5),
        Complex!double(5.1e-4, -0.5e1),
        Complex!double(-5.1e-5, -0.5e2),
        Complex!double(1e-6,2e-6),
        Complex!double(0,2e-6),
        Complex!double(0,2),
        Complex!double(0,20),
        Complex!double(0,200),
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
        Complex!double(double.infinity,double.nan),
        Complex!double(1e-3,double.nan),
        Complex!double(7e-2,7e-2),
        Complex!double(7e-2,-7e-4),
        Complex!double(-9e-2,7e-4),
        Complex!double(-9e-2,9e-2),
        Complex!double(-7e-4,9e-2),
        Complex!double(7e-2,0.9e-2),
        Complex!double(7e-2,1.1e-2)
    ];
    static immutable W = [ 
        Complex!double(-0.5366435657785650339917955593141927494421,
          -5.049143703447034669543036958614140565553),
        Complex!double(0.5366435657785650339917955593141927494421,
          -5.049143703447034669543036958614140565553),
        Complex!double(-0.5366435657785650339917955593141927494421,
          5.049143703447034669543036958614140565553),
        Complex!double(0.5366435657785650339917955593141927494421,
          5.049143703447034669543036958614140565553),
        Complex!double(0.3359473673830576996788000505817956637777e304,
          -0.1999896139679880888755589794455069208455e304),
        Complex!double(0.3584459971462946066523939204836760283645e278,
          0.3818954885257184373734213077678011282505e280),
        Complex!double(0.9996020422657148639102150147542224526887,
          0.00002801044116908227889681753993542916894856),
        Complex!double(-1, 0),
        Complex!double(1, 0),
        Complex!double(0.005754683859034800134412990541076554934877,
          0.1128349818335058741511924929801267822634e-7),
        Complex!double(-0.005529149142341821193633460286828381876955,
          0.005585388387864706679609092447916333443570),
        Complex!double(0.007099365669981359632319829148438283865814,
          0.6149347012854211635026981277569074001219),
        Complex!double(0.3981176338702323417718189922039863062440e8,
          -0.8298176341665249121085423917575122140650e10),
        Complex!double(-double.infinity,
          -double.infinity),
        Complex!double(0.007389128308257135427153919483147229573895,
          0.6149332524601658796226417164791221815139),
        Complex!double(0.4143671923267934479245651547534414976991e8,
          -0.8298168216818314211557046346850921446950e10),
        Complex!double(-double.infinity,
          -double.infinity),
        Complex!double(0.1128379167099649964175513742247082845155e-5,
          0.2256758334191777400570377193451519478895e-5),
        Complex!double(0,
          0.2256758334194034158904576117253481476197e-5),
        Complex!double(0,
          18.56480241457555259870429191324101719886),
        Complex!double(0,
          0.1474797539628786202447733153131835124599e173),
        Complex!double(0,
          double.infinity),
        Complex!double(1,0),
        Complex!double(-1,0),
        Complex!double(0,double.infinity),
        Complex!double(0,-double.infinity),
        Complex!double(double.nan,double.nan),
        Complex!double(double.nan,double.nan),
        Complex!double(double.nan,double.nan),
        Complex!double(double.nan,0),
        Complex!double(0,double.nan),
        Complex!double(double.nan,double.nan),
        Complex!double(double.nan,double.nan),
        Complex!double(double.nan,double.nan),
        Complex!double(0.07924380404615782687930591956705225541145,
          0.07872776218046681145537914954027729115247),
        Complex!double(0.07885775828512276968931773651224684454495,
          -0.0007860046704118224342390725280161272277506),
        Complex!double(-0.1012806432747198859687963080684978759881,
          0.0007834934747022035607566216654982820299469),
        Complex!double(-0.1020998418798097910247132140051062512527,
          0.1010030778892310851309082083238896270340),
        Complex!double(-0.0007962891763147907785684591823889484764272,
          0.1018289385936278171741809237435404896152),
        Complex!double(0.07886408666470478681566329888615410479530,
          0.01010604288780868961492224347707949372245),
        Complex!double(0.07886723099940260286824654364807981336591,
          0.01235199327873258197931147306290916629654)
    ];
    commonTest!cerf(Z, W);
}


/**
Compute erfc(z) = 1 - erf(z), the complex complementary error function,
using w_of_z except for certain regions.
*/
//Steven G. Johnson, October 2012.
Complex!double cerfc(Complex!double z) @safe @nogc nothrow
{

    immutable x = z.re, y = z.im;

    if (!x)
        return Complex!double(1,
                 /* handle y -> Inf limit manually, since
                    exp(y^2) -> Inf but Im[w(y)] -> 0, so
                    IEEE will give us a double.nan when it should be Inf */
                 y*y > 720 ? -copysign(double.infinity, y)
                 : -exp(y*y) * im_w_of_x(y));
    if (!y) {
        if (x*x > 750) // underflow
            return Complex!double(x >= 0 ? 0.0 : 2.0,
                     -y); // preserve sign of 0
        return Complex!double(x >= 0 ? exp(-x*x) * erfcx(x) 
                 : 2. - exp(-x*x) * erfcx(-x),
                 -y); // preserve sign of zero
    }

    immutable mRe_z2 = (y - x) * (x + y); // Re(-z^2), being careful of overflow
    immutable mIm_z2 = -2*x*y; // Im(-z^2)
    if (mRe_z2 < -750) // underflow
        return Complex!double(x >= 0 ? 0.0 : 2.0, 0);

    if (x >= 0)
        return cexp(Complex!double(mRe_z2, mIm_z2))
            * w_of_z(Complex!double(-y,x));
    else
        return 2.0 - cexp(Complex!double(mRe_z2, mIm_z2))
            * w_of_z(Complex!double(y,-x));
}

unittest
{
    static immutable Z = [
        Complex!double(1,2),
        Complex!double(-1,2),
        Complex!double(1,-2),
        Complex!double(-1,-2),
        Complex!double(9,-28),
        Complex!double(21,-33),
        Complex!double(1e3,1e3),
        Complex!double(-3001,-1000),
        Complex!double(1e160,-1e159),
        Complex!double(5.1e-3, 1e-8),
        Complex!double(0,2e-6),
        Complex!double(0,2),
        Complex!double(0,20),
        Complex!double(0,200),
        Complex!double(2e-6,0),
        Complex!double(2,0),
        Complex!double(20,0),
        Complex!double(200,0),
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
        Complex!double(double.infinity,double.nan),
        Complex!double(88,0)
    ];
    static immutable W = [
        Complex!double(1.536643565778565033991795559314192749442,
          5.049143703447034669543036958614140565553),
        Complex!double(0.4633564342214349660082044406858072505579,
          5.049143703447034669543036958614140565553),
        Complex!double(1.536643565778565033991795559314192749442,
          -5.049143703447034669543036958614140565553),
        Complex!double(0.4633564342214349660082044406858072505579,
          -5.049143703447034669543036958614140565553),
        Complex!double(-0.3359473673830576996788000505817956637777e304,
          0.1999896139679880888755589794455069208455e304),
        Complex!double(-0.3584459971462946066523939204836760283645e278,
          -0.3818954885257184373734213077678011282505e280),
        Complex!double(0.0003979577342851360897849852457775473112748,
          -0.00002801044116908227889681753993542916894856),
        Complex!double(2, 0),
        Complex!double(0, 0),
        Complex!double(0.9942453161409651998655870094589234450651,
          -0.1128349818335058741511924929801267822634e-7),
        Complex!double(1,
          -0.2256758334194034158904576117253481476197e-5),
        Complex!double(1,
          -18.56480241457555259870429191324101719886),
        Complex!double(1,
          -0.1474797539628786202447733153131835124599e173),
        Complex!double(1, -double.infinity),
        Complex!double(0.9999977432416658119838633199332831406314,
          0),
        Complex!double(0.004677734981047265837930743632747071389108,
          0),
        Complex!double(0.5395865611607900928934999167905345604088e-175,
          0),
        Complex!double(0, 0),
        Complex!double(0, 0),
        Complex!double(2, 0),
        Complex!double(1, -double.infinity),
        Complex!double(1, double.infinity),
        Complex!double(double.nan, double.nan),
        Complex!double(double.nan, double.nan),
        Complex!double(double.nan, double.nan),
        Complex!double(double.nan, 0),
        Complex!double(1, double.nan),
        Complex!double(double.nan, double.nan),
        Complex!double(double.nan, double.nan),
        Complex!double(0,0)
    ];
    commonTest!cerfc(Z, W);
}

/**
Compute Dawson(z) = sqrt(pi)/2  *  exp(-z^2) * erfi(z),
Dawson's integral for a complex argument,
using w_of_z except for certain regions.
*/
//Steven G. Johnson, October 2012.
Complex!double cdawson(Complex!double z) @safe @nogc nothrow
{

    immutable x = z.re, y = z.im;

    // handle axes separately for speed & proper handling of x or y = Inf or double.nan
    if (!y)
        return Complex!double(spi2 * im_w_of_x(x), -y); // preserve sign of 0
    if (!x) {
        immutable y2 = y*y;
        if (y2 < 2.5e-5) { // Taylor expansion
            return Complex!double(x, // preserve sign of 0
                     y * (1.
                          + y2 * (0.6666666666666666666666666666666666666667
                                  + y2 * 0.26666666666666666666666666666666666667)));
        }
        return Complex!double(x, // preserve sign of 0
                 spi2 * (y >= 0 
                         ? exp(y2) - erfcx(y)
                         : erfcx(-y) - exp(y2)));
    }

    immutable mRe_z2 = (y - x) * (x + y); // Re(-z^2), being careful of overflow
    immutable mIm_z2 = -2*x*y; // Im(-z^2)
    Complex!double mz2 = Complex!double(mRe_z2, mIm_z2); // -z^2

    /* Handle positive and negative x via different formulas,
       using the mirror symmetries of w, to avoid overflow/underflow
       problems from multiplying exponentially large and small quantities. */
    if (y >= 0) {
        if (fabs(y) < 5e-3) {
            if (fabs(x) < 5e-3)
                goto taylor;
            else if (fabs(mIm_z2) < 5e-3)
                goto taylor_realaxis;
        }
        Complex!double res = cexp(mz2) - w_of_z(z);
        return spi2 * Complex!double(-res.im, res.re);
    }
    else { // y < 0
        if (y > -5e-3) { // duplicate from above to avoid fabs(x) call
            if (fabs(x) < 5e-3)
                goto taylor;
            else if (fabs(mIm_z2) < 5e-3)
                goto taylor_realaxis;
        }
        else if (isNaN(y))
            return Complex!double(!x ? 0 : double.nan, double.nan);
        Complex!double res = w_of_z(-z) - cexp(mz2);
        return spi2 * Complex!double(-res.im, res.re);
    }

    // Use Taylor series for small |z|, to avoid cancellation inaccuracy
    //     dawson(z) = z - 2/3 z^3 + 4/15 z^5 + ...
taylor:
    return z * (1.
                + mz2 * (0.6666666666666666666666666666666666666667
                         + mz2 * 0.2666666666666666666666666666666666666667));

    /* for small |y| and small |xy|, 
       use Taylor series to avoid cancellation inaccuracy:
       dawson(x + iy)
       = D + y^2 (D + x - 2Dx^2)
       + y^4 (D/2 + 5x/6 - 2Dx^2 - x^3/3 + 2Dx^4/3)
       + iy [ (1-2Dx) + 2/3 y^2 (1 - 3Dx - x^2 + 2Dx^3)
       + y^4/15 (4 - 15Dx - 9x^2 + 20Dx^3 + 2x^4 - 4Dx^5) ] + ...
       where D = dawson(x) 

       However, for large |x|, 2Dx -> 1 which gives cancellation problems in
       this series (many of the leading terms cancel).  So, for large |x|,
       we need to substitute a continued-fraction expansion for D.

       dawson(x) = 0.5 / (x-0.5/(x-1/(x-1.5/(x-2/(x-2.5/(x...))))))

       The 6 terms shown here seems to be the minimum needed to be
       accurate as soon as the simpler Taylor expansion above starts
       breaking down.  Using this 6-term expansion, factoring out the
       denominator, and simplifying with Maple, we obtain:

       Re dawson(x + iy) * (-15 + 90x^2 - 60x^4 + 8x^6) / x
       = 33 - 28x^2 + 4x^4 + y^2 (18 - 4x^2) + 4 y^4
       Im dawson(x + iy) * (-15 + 90x^2 - 60x^4 + 8x^6) / y
       = -15 + 24x^2 - 4x^4 + 2/3 y^2 (6x^2 - 15) - 4 y^4

       Finally, for |x| > 5e7, we can use a simpler 1-term continued-fraction
       expansion for the real part, and a 2-term expansion for the imaginary
       part.  (This avoids overflow problems for huge |x|.)  This yields:
     
       Re dawson(x + iy) = [1 + y^2 (1 + y^2/2 - (xy)^2/3)] / (2x)
       Im dawson(x + iy) = y [ -1 - 2/3 y^2 + y^4/15 (2x^2 - 4) ] / (2x^2 - 1)

    */
taylor_realaxis:
    immutable x2 = x*x;
    if (x2 > 1600) // |x| > 40
    { 
        immutable y2 = y*y;
        if (x2 > 25e14) {// |x| > 5e7
            immutable xy2 = (x*y)*(x*y);
            return Complex!double((0.5 + y2 * (0.5 + 0.25*y2
                                  - 0.16666666666666666667*xy2)) / x,
                     y * (-1 + y2 * (-0.66666666666666666667
                                     + 0.13333333333333333333*xy2
                                     - 0.26666666666666666667*y2))
                     / (2*x2 - 1));
        }
        return (1. / (-15 + x2*(90 + x2*(-60 + 8*x2)))) *
            Complex!double(x * (33 + x2 * (-28 + 4*x2)
                   + y2 * (18 - 4*x2 + 4*y2)),
              y * (-15 + x2 * (24 - 4*x2)
                   + y2 * (4*x2 - 10 - 4*y2)));
    }
    else 
    {
        immutable D = spi2 * im_w_of_x(x);
        immutable y2 = y*y;
        return Complex!double
            (D + y2 * (D + x - 2*D*x2)
             + y2*y2 * (D * (0.5 - x2 * (2 - 0.66666666666666666667*x2))
                        + x * (0.83333333333333333333
                               - 0.33333333333333333333 * x2)),
             y * (1 - 2*D*x
                  + y2 * 0.66666666666666666667 * (1 - x2 - D*x * (3 - 2*x2))
                  + y2*y2 * (0.26666666666666666667 -
                             x2 * (0.6 - 0.13333333333333333333 * x2)
                             - D*x * (1 - x2 * (1.3333333333333333333
                                                - 0.26666666666666666667 * x2)))));
    }
}

unittest
{
    static immutable Z = [
        Complex!double(2,1),
        Complex!double(-2,1),
        Complex!double(2,-1),
        Complex!double(-2,-1),
        Complex!double(-28,9),
        Complex!double(33,-21),
        Complex!double(1e3,1e3),
        Complex!double(-1000,-3001),
        Complex!double(1e-8, 5.1e-3),
        Complex!double(4.95e-3, -4.9e-3),
        Complex!double(5.1e-3, 5.1e-3),
        Complex!double(0.5, 4.9e-3),
        Complex!double(-0.5e1, 4.9e-4),
        Complex!double(-0.5e2, -4.9e-5),
        Complex!double(0.5e3, 4.9e-6),
        Complex!double(0.5, 5.1e-3),
        Complex!double(-0.5e1, 5.1e-4),
        Complex!double(-0.5e2, -5.1e-5),
        Complex!double(1e-6,2e-6),
        Complex!double(2e-6,0),
        Complex!double(2,0),
        Complex!double(20,0),
        Complex!double(200,0),
        Complex!double(0,4.9e-3),
        Complex!double(0,-5.1e-3),
        Complex!double(0,2e-6),
        Complex!double(0,-2),
        Complex!double(0,20),
        Complex!double(0,-200),
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
        Complex!double(double.infinity,double.nan),
        Complex!double(39, 6.4e-5),
        Complex!double(41, 6.09e-5),
        Complex!double(4.9e7, 5e-11),
        Complex!double(5.1e7, 4.8e-11),
        Complex!double(1e9, 2.4e-12),
        Complex!double(1e11, 2.4e-14),
        Complex!double(1e13, 2.4e-16),
        Complex!double(1e300, 2.4e-303)
    ];
    static immutable W = [
        Complex!double(0.1635394094345355614904345232875688576839,
          -0.1531245755371229803585918112683241066853),
        Complex!double(-0.1635394094345355614904345232875688576839,
          -0.1531245755371229803585918112683241066853),
        Complex!double(0.1635394094345355614904345232875688576839,
          0.1531245755371229803585918112683241066853),
        Complex!double(-0.1635394094345355614904345232875688576839,
          0.1531245755371229803585918112683241066853),
        Complex!double(-0.01619082256681596362895875232699626384420,
          -0.005210224203359059109181555401330902819419),
        Complex!double(0.01078377080978103125464543240346760257008,
          0.006866888783433775382193630944275682670599),
        Complex!double(-0.5808616819196736225612296471081337245459,
          0.6688593905505562263387760667171706325749),
        Complex!double(double.infinity,
          -double.infinity),
        Complex!double(0.1000052020902036118082966385855563526705e-7,
          0.005100088434920073153418834680320146441685),
        Complex!double(0.004950156837581592745389973960217444687524,
          -0.004899838305155226382584756154100963570500),
        Complex!double(0.005100176864319675957314822982399286703798,
          0.005099823128319785355949825238269336481254),
        Complex!double(0.4244534840871830045021143490355372016428,
          0.002820278933186814021399602648373095266538),
        Complex!double(-0.1021340733271046543881236523269967674156,
          -0.00001045696456072005761498961861088944159916),
        Complex!double(-0.01000200120119206748855061636187197886859,
          0.9805885888237419500266621041508714123763e-8),
        Complex!double(0.001000002000012000023960527532953151819595,
          -0.9800058800588007290937355024646722133204e-11),
        Complex!double(0.4244549085628511778373438768121222815752,
          0.002935393851311701428647152230552122898291),
        Complex!double(-0.1021340732357117208743299813648493928105,
          -0.00001088377943049851799938998805451564893540),
        Complex!double(-0.01000200120119126652710792390331206563616,
          0.1020612612857282306892368985525393707486e-7),
        Complex!double(0.1000000000007333333333344266666666664457e-5,
          0.2000000000001333333333323199999999978819e-5),
        Complex!double(0.1999999999994666666666675199999999990248e-5,
          0),
        Complex!double(0.3013403889237919660346644392864226952119,
          0),
        Complex!double(0.02503136792640367194699495234782353186858,
          0),
        Complex!double(0.002500031251171948248596912483183760683918,
          0),
        Complex!double(0,0.004900078433419939164774792850907128053308),
        Complex!double(0,-0.005100088434920074173454208832365950009419),
        Complex!double(0,0.2000000000005333333333341866666666676419e-5),
        Complex!double(0,-48.16001211429122974789822893525016528191),
        Complex!double(0,0.4627407029504443513654142715903005954668e174),
        Complex!double(0,-double.infinity),
        Complex!double(0,0),
        Complex!double(-0,0),
        Complex!double(0, double.infinity),
        Complex!double(0, -double.infinity),
        Complex!double(double.nan, double.nan),
        Complex!double(double.nan, double.nan),
        Complex!double(double.nan, double.nan),
        Complex!double(double.nan, 0),
        Complex!double(0, double.nan),
        Complex!double(double.nan, double.nan),
        Complex!double(double.nan, double.nan),
        Complex!double(0.01282473148489433743567240624939698290584,
          -0.2105957276516618621447832572909153498104e-7),
        Complex!double(0.01219875253423634378984109995893708152885,
          -0.1813040560401824664088425926165834355953e-7),
        Complex!double(0.1020408163265306334945473399689037886997e-7,
          -0.1041232819658476285651490827866174985330e-25),
        Complex!double(0.9803921568627452865036825956835185367356e-8,
          -0.9227220299884665067601095648451913375754e-26),
        Complex!double(0.5000000000000000002500000000000000003750e-9,
          -0.1200000000000000001800000188712838420241e-29),
        Complex!double(5.00000000000000000000025000000000000000000003e-12,
          -1.20000000000000000000018000000000000000000004e-36),
        Complex!double(5.00000000000000000000000002500000000000000000e-14,
          -1.20000000000000000000000001800000000000000000e-42),
        Complex!double(5e-301, 0)
    ];
    commonTest!cdawson(Z, W);
}

unittest
{
    specialTest!(cerf, erf, 1e-20);
    specialTest!(cerfi, erfi, 0);
    specialTest!(cerfcx, erfcx, 0);
    specialTest!(cerfc, erfc, 1e-20);
    specialTest!(cdawson, dawson, 1e-20);
}
