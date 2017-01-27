#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include <initializer_list>
#include <iostream>
#include <functional>

#include "numbers.h"

#pragma once

template <typename T>
struct qr_pair {
	public:
		T quotient;
		T remainder;
};

template <typename T>
class poly;

template <typename T>
std::ostream &operator<<(std::ostream &os, poly<T> const &p);

template <typename T>
class poly {
	private:
		std::vector<T> coeffs;
	
	public:
		poly();
		poly(T constant);
		poly(std::vector<T> coeffs);
		poly(std::initializer_list<T> coeffs);
		poly(const poly<T> &other);
		
		void simplify();
		
		int degree() const;
		T operator[](unsigned int exponent) const;
		void set(unsigned int exponent, T new_value);
		
		poly<T> &operator=(T constant);
		poly<T> &operator=(std::vector<T> coeffs);
		poly<T> &operator=(std::initializer_list<T> coeffs);
		poly<T> &operator=(const poly<T> &other);
		
		bool operator==(const poly<T> &other) const;
		bool operator!=(const poly<T> &other) const;
		
		poly<T> &operator+=(T constant);
		poly<T> &operator-=(T constant);
		poly<T> &operator*=(T constant);
		poly<T> &operator/=(T constant);
		
		poly<T> operator+(T constant) const;
		poly<T> operator-(T constant) const;
		poly<T> operator*(T constant) const;
		poly<T> operator/(T constant) const;
		
		poly<T> operator-() const;
		
		poly<T> &operator+=(const poly<T> &other);
		poly<T> &operator-=(const poly<T> &other);
		poly<T> &operator*=(const poly<T> &other);
		
		poly<T> operator+(const poly<T> &other) const;
		poly<T> operator-(const poly<T> &other) const;
		poly<T> operator*(const poly<T> &other) const;
		
		poly<T> &operator<<=(unsigned int len);
		poly<T> &operator>>=(unsigned int len);
		
		qr_pair<poly<T>> divide(const poly<T> &other) const;
		qr_pair<poly<T>> ring_exact_divide(const poly<T> &other) const;
		qr_pair<poly<T>> pseudo_divide(const poly<T> &other) const;
		poly<T> operator/(const poly<T> &other) const;
		poly<T> operator%(const poly<T> &other) const;
		poly<T> &operator/=(const poly<T> &other);
		poly<T> &operator%=(const poly<T> &other);
		
		friend std::ostream &operator<<<>(std::ostream &os, const poly<T> &p);
		
		template <typename U>
		operator poly<U>();
		
		template <typename U>
		poly<U> convert(std::function<U(T)> converter) const;
		
		poly<T> power_mod(mpz_class power);
		poly<T> power_mod(poly<T> start, mpz_class power);
		
		poly<T> derivative();
		poly<T> reverse();
		T content();
		T norm();
		T leading();
		
		poly<T> compose(poly<T> x);
		T evaluate(T x);
};

template <typename T>
class util<poly<T>> {
public:
	static poly<T> zero();
	static poly<T> zero(const poly<T> &reference);
	static poly<T> one();
	static poly<T> one(const poly<T> &reference);
	static poly<T> from_int(int n, const poly<T> &reference);
	static poly<T> get_gcd(const poly<T> &a, const poly<T> &b);
	static poly<T> get_pow(const poly<T> &a, int exp);
};

template <typename T>
poly<T> util<poly<T>>::zero(const poly<T> &reference) {
	return poly<T>();
}

template <typename T>
poly<T> util<poly<T>>::one(const poly<T> &reference) {
	return poly<T>(util<T>::one(reference[reference.degree()]));
}

template <typename T>
poly<T> util<poly<T>>::zero() {
	return poly<T>();
}

template <typename T>
poly<T> util<poly<T>>::one() {
	return poly<T>(util<T>::one());
}

template <typename T>
poly<T> util<poly<T>>::from_int(int n, const poly<T> &reference) {
	return poly<T>(util<T>::from_int(n, reference[reference.degree()]));
}

template <typename T>
poly<T> util<poly<T>>::get_pow(const poly<T> &a, int exp) {
	poly<T> result = util<poly<T>>::one(a);
	for (int i = 0; i < exp; i++)
		result *= a;
	return result;
}

template <typename T>
poly<T>::poly() {
	this->coeffs = std::vector<T>();
}

template <typename T>
poly<T>::poly(T constant) {
	this->coeffs = std::vector<T>();
	this->coeffs.push_back(constant);
	this->simplify();
}

template <typename T>
poly<T>::poly(std::vector<T> coeffs) {
	this->coeffs = coeffs;
	this->simplify();
}

template <typename T>
poly<T>::poly(std::initializer_list<T> coeffs) {
	this->coeffs.insert(this->coeffs.end(), coeffs.begin(), coeffs.end());
	this->simplify();
}

template <typename T>
poly<T>::poly(const poly<T> &other) {
	this->coeffs = other.coeffs;
}

template <typename T>
void poly<T>::simplify() {
	while (this->coeffs.size() > 0 && this->coeffs[this->coeffs.size() - 1] == util<T>::zero(this->coeffs[this->coeffs.size() - 1]))
		this->coeffs.pop_back();
}

template <typename T>
int poly<T>::degree() const {
	return this->coeffs.size() - 1;
}

template <typename T>
T poly<T>::operator[](unsigned int exponent) const {
	if (this->degree() < 0)
		return util<T>::zero();
	if (exponent > this->degree())
		return util<T>::zero(this->coeffs[this->degree()]);
	return this->coeffs[exponent];
}

template <typename T>
void poly<T>::set(unsigned int exponent, T new_value) {
	while (exponent > this->degree())
		this->coeffs.push_back(util<T>::zero(this->coeffs[this->degree()]));
	this->coeffs[exponent] = new_value;
	this->simplify();
}

template <typename T>
poly<T> &poly<T>::operator=(T constant) {
	this->coeffs = std::vector<T>();
	this->coeffs.push_back(constant);
	this->simplify();
	return *this;
}

template <typename T>
poly<T> &poly<T>::operator=(std::vector<T> coeffs) {
	this->coeffs = coeffs;
	this->simplify();
	return *this;
}

template <typename T>
poly<T> &poly<T>::operator=(std::initializer_list<T> coeffs) {
	this->coeffs = std::vector<T>();
	this->coeffs.insert(this->coeffs.end(), coeffs.begin(), coeffs.end());
	this->simplify();
	return *this;
}

template <typename T>
poly<T> &poly<T>::operator=(const poly<T> &other) {
	this->coeffs = other.coeffs;
	return *this;
}

template <typename T>
bool poly<T>::operator==(const poly<T> &other) const {
	return (this->coeffs == other.coeffs);
}

template <typename T>
bool poly<T>::operator!=(const poly<T> &other) const {
	return (this->coeffs != other.coeffs);
}

template <typename T>
poly<T> &poly<T>::operator+=(T constant) {
	this->coeffs[0] += constant;
	this->simplify();
	return *this;
}

template <typename T>
poly<T> &poly<T>::operator-=(T constant) {
	this->coeffs[0] -= constant;
	this->simplify();
	return *this;
}

template <typename T>
poly<T> &poly<T>::operator*=(T constant) {
	for (int i = 0; i < this->coeffs.size(); i++)
		this->coeffs[i] *= constant;
	this->simplify();
	return *this;
}

template <typename T>
poly<T> &poly<T>::operator/=(T constant) {
	for (int i = 0; i < this->coeffs.size(); i++)
		this->coeffs[i] /= constant;
	this->simplify();
	return *this;
}

template <typename T>
poly<T> poly<T>::operator+(T constant) const {
	return poly<T>(*this) += constant;
}

template <typename T>
poly<T> poly<T>::operator-(T constant) const {
	return poly<T>(*this) -= constant;
}

template <typename T>
poly<T> poly<T>::operator*(T constant) const {
	return poly<T>(*this) *= constant;
}

template <typename T>
poly<T> poly<T>::operator/(T constant) const {
	return poly<T>(*this) /= constant;
}

template <typename T>
poly<T> operator*(T constant, const poly<T> &p) {
	return p*constant;
}

template <typename T>
poly<T> poly<T>::operator-() const {
	poly<T> neg;
	for (int i = 0; i < this->coeffs.size(); i++)
		neg.coeffs.push_back(-this->coeffs[i]);
	return neg;
}

template <typename T>
poly<T> &poly<T>::operator+=(const poly<T> &p) {
	if (this->coeffs.size() == 0)
		return (*this = p);
	if (p.coeffs.size() == 0)
		return *this;
	int d = this->degree();
	if (d > p.degree())
		d = p.degree();
	while (this->degree() < p.degree())
		this->coeffs.push_back(p.coeffs[this->coeffs.size()]);
	for (int i = 0; i <= d; i++)
		this->coeffs[i] += p.coeffs[i];
	this->simplify();
	return *this;
}

template <typename T>
poly<T> &poly<T>::operator-=(const poly<T> &p) {
	if (this->coeffs.size() == 0)
		return (*this = (-p));
	if (p.coeffs.size() == 0)
		return *this;
	int d = this->degree();
	if (d > p.degree())
		d = p.degree();
	while (this->degree() < p.degree())
		this->coeffs.push_back(-p.coeffs[this->degree() + 1]);
	for (int i = 0; i <= d; i++)
		this->coeffs[i] -= p.coeffs[i];
	this->simplify();
	return *this;
}

template <typename T>
poly<T> poly<T>::operator+(const poly<T> &p) const {
	return poly<T>(*this) += p;
}

template <typename T>
poly<T> poly<T>::operator-(const poly<T> &p) const {
	return poly<T>(*this) -= p;
}

template <typename T>
poly<T> poly<T>::operator*(const poly<T> &p) const {
	if (this->coeffs.size() == 0)
		return *this;
	if (p.coeffs.size() == 0)
		return p;
	poly<T> ret = poly<T>();
	for (int i = 0; i <= p.degree() + this->degree(); i++)
		ret.coeffs.push_back(util<T>::zero(this->coeffs[this->coeffs.size() - 1]));
	for (int i = 0; i <= p.degree(); i++) {
		for (int j = 0; j <= this->degree(); j++) {
			ret.coeffs[i+j] += this->coeffs[j]*p.coeffs[i];
		}
	}
	ret.simplify();
	return ret;
}

template <typename T>
poly<T> &poly<T>::operator*=(const poly<T> &p) {
	*this = (*this) * p;
}

template <typename T>
poly<T> &poly<T>::operator<<=(unsigned int len) {
	if (this->coeffs.size() == 0)
		return *this;
	for (int i = 0; i < len; i++)
		this->coeffs.insert(this->coeffs.begin(), util<T>::zero(this->coeffs[0]));
	return *this;
}

template <typename T>
poly<T> &poly<T>::operator>>=(unsigned int len) {
	if (this->coeffs.size() == 0)
		return *this;
	for (int i = 0; i < len; i++)
		this->coeffs.erase(this->coeffs.begin());
	return *this;
} 

template <typename T>
qr_pair<poly<T>> poly<T>::divide(const poly<T> &other) const {
	// Algorithm 3.1.1

	poly<T> r = *this, q;
	T invlb = util<T>::one(other[other.degree()])/(other[other.degree()]);
	while (r.degree() >= other.degree()) {
		poly<T> s = poly<T>(r[r.degree()]*invlb);
		s <<= r.degree()-other.degree();
		q += s;

		poly<T> sb = other*r[r.degree()]*invlb;
		sb <<= r.degree()-other.degree();
		r -= sb;
	}

	qr_pair<poly<T>> qr;
	qr.quotient = q;
	qr.remainder = r;
	return qr;
}

template <typename T>
qr_pair<poly<T>> poly<T>::ring_exact_divide(const poly<T> &other) const {
	// Algorithm 3.1.1

	poly<T> r = *this, q;
	while (r.degree() >= other.degree()) {
		poly<T> s = poly<T>(r[r.degree()]/(other[other.degree()]));
		s <<= r.degree()-other.degree();
		q += s;

		poly<T> sb = other*r[r.degree()]/(other[other.degree()]);
		sb <<= r.degree()-other.degree();
		r -= sb;
	}

	qr_pair<poly<T>> qr;
	qr.quotient = q;
	qr.remainder = r;
	return qr;
}

template <typename T>
qr_pair<poly<T>> poly<T>::pseudo_divide(const poly<T> &other) const {
	// Algorithm 3.1.2

	poly<T> r = *this, q;
	T d = other[other.degree()];
	int e = this->degree() - other.degree() + 1;
	while (r.degree() >= other.degree()) {
		poly<T> s = poly<T>(r[r.degree()]);
		s <<= r.degree()-other.degree();
		poly<T> s_other = other * r[r.degree()];
		s_other <<= r.degree()-other.degree();
		
		q *= d;
		q += s;
		r *= d;
		r -= s_other;
		e--;
	}
	
	T de = util<T>::one(d);
	for (int i = 0; i < e; i++)
		de *= d;

	qr_pair<poly<T>> qr;
	qr.quotient = de*q;
	qr.remainder = de*r;
	return qr;
}

template <typename T>
poly<T> poly<T>::operator/(const poly<T> &p) const {
	return this->divide(p).quotient;
}

template <typename T>
poly<T> poly<T>::operator%(const poly<T> &p) const {
	return this->divide(p).remainder;
}

template <typename T>
poly<T> &poly<T>::operator/=(const poly<T> &p) {
	*this = (*this) / p;
}

template <typename T>
poly<T> &poly<T>::operator%=(const poly<T> &p) {
	*this = (*this) % p;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const poly<T> &p) {
	os << "(";
	bool first = true;
	for (int i = 0; i < p.coeffs.size(); i++) {
		if (p.coeffs[i] == util<T>::zero(p.coeffs[i]))
			continue;
		if (!first)
			os << " + ";
		first = false;
		os << p.coeffs[i];
		if (i == 1)
			os << "x";
		else if (i > 1)
			os << "x^" << i;
	}
	if (first)
		os << "0";
	os << ")";
	return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const qr_pair<T> &qr) {
	return os << "(" << qr.quotient << ", " << qr.remainder << ")";
}

template <typename T>
template <typename U>
poly<T>::operator poly<U>() {
	std::vector<U> vec;
	for (int i = 0; i < this->coeffs.size(); i++)
		vec.push_back(static_cast<U>(this->coeffs[i]));
	return poly<U>(vec);
}

template <typename T>
template <typename U>
poly<U> poly<T>::convert(std::function<U(T)> converter) const {
	std::vector<U> vec;
	// std::cout << "sanity check, this is " << this << std::endl;
	for (int i = 0; i < this->coeffs.size(); i++) {
		U coeff = converter(this->coeffs[i]);
		// std::cout << "coeffs[" << i << "] = " << coeff << std::endl;
		vec.push_back(coeff);
	}
	return poly<U>(vec);
}

template <typename T, typename U>
poly<U> global_convert(poly<T> original, std::function<U(T)> converter) {
	return original.convert(converter);
}

template <typename T, typename U>
std::function<poly<U>(poly<T>)> curried_convert(std::function<U(T)> converter) {
	return std::bind(global_convert<T, U>, std::placeholders::_1, converter);
}

template <typename T>
poly<T> constant(T con) {
	return poly<T>(con);
}

template <typename T>
std::function<poly<T>(T)> make_constant() {
	return std::bind(constant<T>, std::placeholders::_1);
}

template <typename T>
poly<T> poly<T>::power_mod(mpz_class power) {
	// N.B. This part of the code technically fails if you try to
	// use it with a number that has more than 2^31 bits.
	// Said number would also take up about 268 MB of RAM to store,
	// so if you need to work with polynomials whose coefficients
	// are that large, ask your local supercomputer instead of some
	// random undergrad who's just trying to do his thesis.
	int numbits = mpz_sizeinbase(power.get_mpz_t(), 2);
	std::vector<poly<T>> squares;
	squares.push_back(poly<T>({util<T>::zero(this->coeffs[this->degree()]), util<T>::one(this->coeffs[this->degree()])}));
	for (int i = 1; i < numbits; i++) {
		poly<T> current = squares[squares.size() - 1];
		current *= current;
		current %= *this;
		squares.push_back(current);
	}
	
	poly<T> product = poly<T>(util<T>::one(this->coeffs[this->degree()]));
	for (int i = 0; i < numbits; i++) {
		if (mpz_tstbit(power.get_mpz_t(), i)) {
			product *= squares[i];
			product %= *this;
		}
	}
	
	return product;
}

template <typename T>
poly<T> poly<T>::power_mod(poly<T> start, mpz_class power) {
	int numbits = mpz_sizeinbase(power.get_mpz_t(), 2);
	std::vector<poly<T>> squares;
	squares.push_back(start);
	for (int i = 1; i < numbits; i++) {
		poly<T> current = squares[squares.size() - 1];
		current *= current;
		current %= *this;
		squares.push_back(current);
	}
	
	poly<T> product = poly<T>(util<T>::one(this->coeffs[this->degree()]));
	for (int i = 0; i < numbits; i++) {
		if (mpz_tstbit(power.get_mpz_t(), i)) {
			product *= squares[i];
			product %= *this;
		}
	}
	
	return product;
}

template <typename T>
poly<T> poly<T>::derivative() {
	if (this->degree() < 1)
		return poly<T>();
	
	std::vector<T> coeffs;
	for (int i = 1; i <= this->degree(); i++) {
		coeffs.push_back(util<T>::from_int(i, this->coeffs[i]) * this->coeffs[i]);
	}
	
	return poly<T>(coeffs);
}

template <typename T>
poly<T> poly<T>::reverse() {
	std::vector<T> coeffs;
	for (int i = 0; i <= this->degree(); i++) {
		coeffs.insert(coeffs.begin(), this->coeffs[i]);
	}
	
	return poly<T>(coeffs);
}

template <typename T>
T poly<T>::content() {
	T g = this->coeffs[this->degree()];
	for (int i = 0; i < this->degree(); i++) {
		if (this->coeffs[i] == util<T>::zero(this->coeffs[i]))
			continue;
		g = util<T>::get_gcd(this->coeffs[i], g);
	}
	
	return g;
}

template <typename T>
T poly<T>::norm() {
	T s = util<T>::get_abs(this->coeffs[this->degree()])*util<T>::get_abs(this->coeffs[this->degree()]);
	for (int i = 0; i < this->degree(); i++) {
		s += util<T>::get_abs(this->coeffs[i])*util<T>::get_abs(this->coeffs[i]);
	}
	
	return util<T>::get_sqrt(s);
}

template <typename T>
T poly<T>::leading() {
	//if (this->degree() < 0)
	//	return util<T>::zero();
	return this->coeffs[this->degree()];
}

template <typename T>
poly<T> poly<T>::compose(poly<T> x) {
	poly<T> result = util<poly<T>>::zero(x);
	for (int i = 0; i <= this->degree(); i++) {
		poly<T> term = this->coeffs[i];
		for (int j = 0; j < i; j++)
			term *= x;
		result += term;
	}
	
	return result;
}

template <typename T>
T poly<T>::evaluate(T x) {
	T result = util<T>::zero(x);
	for (int i = 0; i <= this->degree(); i++) {
		T term = this->coeffs[i];
		for (int j = 0; j < i; j++)
			term *= x;
		result += term;
	}
	
	return result;
}

template <typename T>
poly<poly<T>> switch_variables(poly<poly<T>> orig) {
	if (orig.degree() < 0)
		return orig;
	poly<T> x = util<poly<T>>::one(orig[orig.degree()]);
	x <<= 1;
	poly<poly<poly<T>>> mid = orig.convert(curried_convert(make_constant<T>()));
	return mid.evaluate(poly<poly<T>>(x));
}
