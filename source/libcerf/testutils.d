/**
Internal test utilities.
Use flags -unittest -debug=libcerf to receive debug info.

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
module libcerf.testutils;

import std.complex;
import std.math : isNaN, isInfinity, isFinite;
import core.stdc.math : fabs;

debug(libcerf)
{
    import std.traits : fullyQualifiedName;
    import std.stdio : writefln;
}

import core.stdc.math : pow, erfc;

package:

///
bool relativeErrorCheck(double a, double b)
@safe nothrow @nogc
{
    if(a.isNaN)
    {
        return b.isNaN;
    }
    if(a.isInfinity)
    {
        return b == a;
    }
    if(a == 0)
    {
        return fabs(b) < 1e-13;
    }
    return fabs((b-a) / a) < 1e-13;
}

///
bool relativeErrorCheck(Complex!double a, Complex!double b)
@safe nothrow @nogc
{
    return relativeErrorCheck(a.re, b.re) && relativeErrorCheck(a.im, b.im);
}

///
void commonTest(alias cfun)(in Complex!double[] Z, in Complex!double[] W)
{
    debug(libcerf){
        enum name = fullyQualifiedName!cfun;
        writefln("%s common tests...", name);
        scope (success)
        writefln("%s common tests success", name);
    }
    assert(W.length == Z.length);   
    foreach(i; 0..Z.length) {
        immutable z = Z[i];
        immutable w = W[i];
        immutable f = cfun(z);
        auto e = f-w;
        e.re /= w.re;
        e.im /= w.im;
        debug(libcerf) 
            writefln("%s(%s) = %s, vs. %s, rel. err. = %s)", name, z, f, w, e);
        assert(relativeErrorCheck(w, f));
    }
}

///
void specialTest(alias cfun, alias fun, double C)()
{
    debug(libcerf){
        enum name = fullyQualifiedName!cfun;
        writefln("%s special tests...", name);
        scope (success)
        writefln("%s special tests success", name);
    }
    foreach (i; 0..10000) 
    {
        immutable x = pow(10., -300. + i * 600. / (10000 - 1));
        assert(relativeErrorCheck(fun(+x), cfun(Complex!double(+x,x*C)).re));
        assert(relativeErrorCheck(fun(-x), cfun(Complex!double(-x,x*C)).re));
    }
    assert(relativeErrorCheck(fun(double.infinity), cfun(Complex!double(double.infinity,0.)).re));
    assert(relativeErrorCheck(fun(-double.infinity), cfun(Complex!double(-double.infinity,0.)).re));
    assert(relativeErrorCheck(fun(double.nan), cfun(Complex!double(double.nan,0.)).re));
}
