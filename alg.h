#include <gmp.h>
#include <gmpxx.h>
#include <utility>
#include <tuple>
#include "polyring.h"
#include "modring.h"
#include "typedefs.h"

template <typename T>
std::tuple<poly<T>, poly<T>, poly<T>> extended_gcd(poly<T> a, poly<T> b) {
	// Algorithm 3.2.2

	// The following is a weird way to make 1, but it ensures that any
	// "extranumerary" info from a and b (such as modular base) is retained.
	// Note that b has to be nonzero for GCD to make sense, so we know b has
	// a nonzero leading coefficient.
	poly<T> u = poly<T>(b[b.degree()]/b[b.degree()]);
	poly<T> d = a;
	poly<T> v1 = poly<T>();
	poly<T> v3 = b;

	while (v3 != poly<T>()) {
		qr_pair<poly<T>> qr = d.divide(v3);
		poly<T> q = qr.quotient;
		poly<T> r = qr.remainder;
		poly<T> t = u - v1*q;
		u = v1;
		d = v3;
		v1 = t;
		v3 = r;
	}
		
	poly<T> v = (d - a*u)/b;
	return std::make_tuple(u, v, d);
}

std::pair<Z_X, Z_X> hensel_lift(Z p, Z q, Z_X a, Z_X b, Z_X c, Z_X u, Z_X v);
std::pair<Z_X, Z_X> quad_hensel_lift(Z p, Z q, Z_X a1, Z_X b1, Z_X u, Z_X v);