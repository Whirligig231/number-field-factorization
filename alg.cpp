#include "alg.h"

std::pair<Z_X, Z_X> hensel_lift(Z p, Z q, Z_X a, Z_X b, Z_X c, Z_X u, Z_X v) {
	// Algorithm 3.5.5

	Z r = gcd(p, q);
	ZN_X f = ((c - a*b)/q).convert(to_mod(r));

	// The book isn't clear on what to do here except to hint that we need to
	// do polynomial division in order to find t such that:
	// (v*f - a*t).degree() < a.degree()
	// Note that if we divide v*f by a, we get t satisfying this criterion.
	ZN_X t = (v.convert(to_mod(r))*f) / a.convert(to_mod(r));

	Z_X a0 = static_cast<Z_X>((v.convert(to_mod(r))*f) - (a.convert(to_mod(r))*t));
	Z_X b0 = static_cast<Z_X>((u.convert(to_mod(r))*f) + (b.convert(to_mod(r))*t));

	Z_X a1 = a + q*a0;
	Z_X b1 = b + q*b0;

	return std::make_pair(a1, b1);
}