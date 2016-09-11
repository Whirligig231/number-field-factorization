#include <gmp.h>
#include <gmpxx.h>
#include <Eigen/Dense>
#include <utility>
#include <tuple>
#include "numbers.h"
#include "polyring.h"
#include "modring.h"
#include "typedefs.h"

template <typename T>
std::vector<vec<T>> kernel(mat<T> m) {
	// Algorithm 2.3.1
	// Note that a lot of values here have changed slightly from the book;
	// this is to avoid nasty off-by-one accesses in matrices and lists.
	// Eigen 0-indexes its matrix coefficients.
	
	int r = 0;
	std::vector<int> ci;
	for (int i = 0; i < m.rows(); i++)
		ci.push_back(-1);
	std::vector<int> di;
	// The book says this isn't necessary, but I do it in order to avoid
	// out-of-bounds accesses
	for (int i = 0; i < m.cols(); i++)
		di.push_back(-1);
	
	for (int k = 0; k < m.cols(); k++) {
		int j = -1;
		for (int i = 0; i < m.rows(); i++)
			if (m(i, k) != zero<T>(m(0, 0))*m(i, k) && ci[i] == -1)
				j = i;
			
		if (j == -1) {
			r++;
			di[k] = -1;
		}
		else {
			T d = -one<T>(m(j, k))/m(j, k);
			m(j, k) = -one<T>(m(j, k));
			
			for (int s = k+1; s < m.cols(); s++)
				m(j, s) *= d;
			
			for (int i = 0; i < m.rows(); i++) {
				if (i == j)
					continue;
				
				d = m(i, k);
				m(i, k) = zero<T>(m(i, k));
				
				for (int s = k+1; s < m.cols(); s++)
					m(i, s) += d*m(j, s);
			}
			
			ci[j] = k;
			di[k] = j;
		}
	}
	
	std::vector<vec<T>> ret;
	for (int k = 0; k < m.cols(); k++) {
		if (di[k] == -1) {
			vec<T> x(m.cols());
			for (int i = 0; i < m.cols(); i++) {
				if (di[i] >= 0)
					x(i) = m(di[i], k);
				else if (i == k)
					x(i) = one<T>(m(0, 0));
				else
					x(i) = zero<T>(m(0, 0));
			}
			
			ret.push_back(x);
		}
	}
	
	return ret;
}

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