#include <gmp.h>
#include <gmpxx.h>
#include <Eigen/Dense>
#include <utility>
#include <tuple>
#include "numbers.h"
#include "polyring.h"
#include "modring.h"
#include "polymodring.h"
#include "typedefs.h"

#pragma once

Z choose(int n, int r);
int log_bound(Z base, Z pow);

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
			if (m(i, k) != util<T>::zero(m(0, 0))*m(i, k) && ci[i] == -1)
				j = i;
			
		if (j == -1) {
			r++;
			di[k] = -1;
		}
		else {
			T d = -util<T>::one(m(j, k))/m(j, k);
			m(j, k) = -util<T>::one(m(j, k));
			
			for (int s = k+1; s < m.cols(); s++)
				m(j, s) *= d;
			
			for (int i = 0; i < m.rows(); i++) {
				if (i == j)
					continue;
				
				d = m(i, k);
				m(i, k) = util<T>::zero(m(i, k));
				
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
					x(i) = util<T>::one(m(0, 0));
				else
					x(i) = util<T>::zero(m(0, 0));
			}
			
			ret.push_back(x);
		}
	}
	
	return ret;
}

template <typename T>
std::tuple<poly<T>, poly<T>, poly<T>> extended_gcd(poly<T> a, poly<T> b) {
	// Algorithm 3.2.2
	
	if (b.degree() < 0)
		return std::make_tuple(poly<T>(util<T>::one(a[a.degree()])), poly<T>(util<T>::zero(a[a.degree()])), a);
	if (a.degree() < 0)
		return std::make_tuple(poly<T>(util<T>::one(b[b.degree()])), poly<T>(util<T>::zero(b[b.degree()])), b);

	poly<T> u = poly<T>(util<T>::one(b[b.degree()]));
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

template <typename T>
poly<T> sub_resultant_gcd(poly<T> a, poly<T> b) {
	// Algorithm 3.3.1
	
	if (b.degree() > a.degree())
		return sub_resultant_gcd(b, a);
	if (b.degree() < 0)
		return a;
	
	T ac = a.content();
	T bc = b.content();
	T d = util<T>::get_gcd(ac, bc);
	a /= ac;
	b /= bc;
	T g = util<T>::one(a[a.degree()]);
	T h = g;
	int c = 0;
	int delta = a.degree() - b.degree();
	while (true) {
		poly<T> dividend = a;
		//for (int i = 0; i < delta + 1; i++) This is stupid, why did I do it?
		//	dividend *= b[b.degree()];
		poly<T> r = dividend.pseudo_divide(b).remainder;
		if (r.degree() < 0)
			break;
		if (r.degree() == 0) {
			b = util<T>::one(r[r.degree()]);
			break;
		}
		
		a = b;
		b = r;
		b /= g;
		for (int i = 0; i < delta; i++)
			b /= h;
		c++;

		if (c < 10) {
			g = a[a.degree()];
			T h_temp = h;
			for (int i = 0; i < delta; i++)
				h *= g;
			for (int i = 0; i < delta; i++)
				h /= h_temp;
		}
		else {
			a /= a.content();
			b /= b.content();
			g = util<T>::one(a[a.degree()]);
			h = g;
			c = 0;
		}
	}
	
	return (b/b.content()) * d;
}

template <typename T>
T sub_resultant(poly<T> a, poly<T> b) {
	// Algorithm 3.3.7
	
	if (a.degree() < 0)
		return util<T>::zero(b[b.degree()]);
	if (b.degree() < 0)
		return util<T>::zero(a[a.degree()]);
	
	T a_cont = a.content();
	T b_cont = b.content();
	a /= a_cont;
	b /= b_cont;
	T g = util<T>::one(a[a.degree()]);
	T h = util<T>::one(a[a.degree()]);
	T s = util<T>::one(a[a.degree()]);
	T t = get_pow(a_cont, b.degree())*get_pow(b_cont, a.degree());
	if (a.degree() < b.degree()) {
		poly<T> temp = a;
		a = b;
		b = temp;
	}
	
	if (a.degree() % 2 && b.degree() % 2)
		s = -s;
	
	do {
		int delta = a.degree() - b.degree();
		
		if (a.degree() % 2 && b.degree() % 2)
			s = -s;

		poly<T> r = a.pseudo_divide(b).remainder;
		
		a = b;
		b = r;
		b /= g;
		for (int i = 0; i < delta; i++)
			b /= h;
		
		g = a[a.degree()];
		T h2 = get_pow(g, delta);
		if (delta == 0)
			h2 *= h;
		for (int i = 1; i < delta; i++)
			h2 /= h;
		h = h2;
	} while (b.degree() > 0);
	
	T h2 = get_pow(b[b.degree()], a.degree());
	if (a.degree() == 0)
		h2 *= h;
	for (int i = 1; i < a.degree(); i++)
		h2 /= h;
	
	return s*t*h2;
}

std::vector<ZN_X> berlekamp_small_p(ZN_X a);
std::vector<ZN_X> berlekamp(ZN_X a);
std::vector<ZN_X> berlekamp_auto(ZN_X a);

Z coeff_bound(Z_X a);

std::pair<Z_X, Z_X> hensel_lift(Z p, Z q, Z_X a, Z_X b, Z_X c, Z_X u, Z_X v);
std::pair<Z_X, Z_X> quad_hensel_lift(Z p, Z q, Z_X a1, Z_X b1, Z_X u, Z_X v);
std::pair<Z_X, Z_X> multi_hensel_lift(Z p, int exp, Z_X a, Z_X b, Z_X c);
std::vector<Z_X> poly_hensel_lift(Z p, int exp, std::vector<Z_X> ai, Z_X c);

std::vector<Z_X> factor(Z_X a);
std::vector<Q_X> factor(Q_X a);

template <typename T>
std::vector<poly<polymod<T>>> factor(poly<polymod<T>> a) {
	// Algorithm 3.6.4
	// Note that 3.6.4 assumes Q as the base field;
	// other number fields will also work.
	
	if (a.degree() < 0)
		return std::vector<poly<polymod<T>>>({a});

	poly<polymod<T>> u = a / sub_resultant_gcd(a, a.derivative());

	poly<poly<T>> g = static_cast<poly<poly<T>>>(u);
	T k = util<T>::zero(a.leading().get_value().leading());
	
	
}