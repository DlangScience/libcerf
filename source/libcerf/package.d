/**
Package module with aliases.

Example:
--------------------
import std.complex, std.math;
import libcerf;

auto c = erfi(complex(1.0, 0.0));
auto f = erfi(1.0);

assert(fabs((c.re - f)/f) < 1e-13);
--------------------

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

See_Also:
    Original $(LINK2 http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package, Faddeeva Package)

    C $(LINK2 http://apps.jcns.fz-juelich.de/libcerf, source code)

    D source code: $(LINK2 http://github.com/9il/libcerf, github)

Version: 1.4

Date: September 21, 2014
*/
module libcerf;

static import libcerf.err_fcts;
static import libcerf.erfcx_;
static import libcerf.w_of_z;
static import libcerf.im_w_of_x;

///Alias for libcerf.w_of_z.w_of_z.
alias fadeeva = libcerf.w_of_z.w_of_z;

///Alias for libcerf.im_w_of_x.im_w_of_x.
alias fadeevaIm = libcerf.im_w_of_x.im_w_of_x;

///Alias for libcerf.err_fcts.cerfcx and libcerf.erfcx_.erfcx.
alias erfcx = libcerf.err_fcts.cerfcx;
///ditto
alias erfcx = libcerf.erfcx_.erfcx;

///Alias for libcerf.err_fcts.cerfi and libcerf.err_fcts.erfi.
alias erfi = libcerf.err_fcts.cerfi;
///ditto
alias erfi = libcerf.err_fcts.erfi;

///Alias for libcerf.err_fcts.cdawson and libcerf.err_fcts.dawson.
alias dawson = libcerf.err_fcts.cdawson;
///ditto
alias dawson = libcerf.err_fcts.dawson;

///Alias for libcerf.err_fcts.voigt
alias voigt = libcerf.err_fcts.voigt;

///Alias for libcerf.err_fcts.cerf.
alias erf = libcerf.err_fcts.cerf;

///Alias for libcerf.err_fcts.cerfc.
alias erfc = libcerf.err_fcts.cerf;
