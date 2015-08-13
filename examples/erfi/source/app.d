import std.complex, std.math;
import libcerf;

void main()
{
	auto c = erfi(complex(1.0, 0.0));
	auto f = erfi(1.0);

	assert(fabs((c.re - f)/f) < 1e-13);
}
