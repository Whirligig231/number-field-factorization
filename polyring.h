#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include <initializer_list>
#include <iostream>
#include <functional>

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
		
		unsigned int degree() const;
		const T &operator[](unsigned int exponent) const;
		
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
		
		const poly<T> operator+(T constant) const;
		const poly<T> operator-(T constant) const;
		const poly<T> operator*(T constant) const;
		const poly<T> operator/(T constant) const;
		
		const poly<T> operator-() const;
		
		poly<T> &operator+=(const poly<T> &other);
		poly<T> &operator-=(const poly<T> &other);
		poly<T> &operator*=(const poly<T> &other);
		
		const poly<T> operator+(const poly<T> &other) const;
		const poly<T> operator-(const poly<T> &other) const;
		const poly<T> operator*(const poly<T> &other) const;
		
		poly<T> &operator<<=(unsigned int len);
		poly<T> &operator>>=(unsigned int len);
		
		const qr_pair<poly<T>> divide(const poly<T> &other) const;
		const poly<T> operator/(const poly<T> &other) const;
		const poly<T> operator%(const poly<T> &other) const;
		poly<T> &operator/=(const poly<T> &other);
		poly<T> &operator%=(const poly<T> &other);
		
		friend std::ostream &operator<<<>(std::ostream &os, const poly<T> &p);
		
		template <typename U>
		operator poly<U>();
		
		template <typename U>
		const poly<U> convert(std::function<U(T)> converter) const;
};

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
	while (this->coeffs[this->coeffs.size() - 1] == static_cast<T>(0))
		this->coeffs.pop_back();
}

template <typename T>
unsigned int poly<T>::degree() const {
	if (this->coeffs.size() == 0)
		return 0;
	return this->coeffs.size() - 1;
}

template <typename T>
const T &poly<T>::operator[](unsigned int exponent) const {
	return this->coeffs[exponent];
}

template <typename T>
poly<T> &poly<T>::operator=(T constant) {
	this->coeffs = std::vector<T>();
	this->coeffs.push_back(constant);
	this->simplify();
}

template <typename T>
poly<T> &poly<T>::operator=(std::vector<T> coeffs) {
	this->coeffs = coeffs;
	this->simplify();
}

template <typename T>
poly<T> &poly<T>::operator=(std::initializer_list<T> coeffs) {
	this->coeffs = std::vector<T>();
	this->coeffs.insert(this->coeffs.end(), coeffs.begin(), coeffs.end());
	this->simplify();
}

template <typename T>
poly<T> &poly<T>::operator=(const poly<T> &other) {
	this->coeffs = other.coeffs;
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
const poly<T> poly<T>::operator+(T constant) const {
	return poly<T>(*this) += constant;
}

template <typename T>
const poly<T> poly<T>::operator-(T constant) const {
	return poly<T>(*this) -= constant;
}

template <typename T>
const poly<T> poly<T>::operator*(T constant) const {
	return poly<T>(*this) *= constant;
}

template <typename T>
const poly<T> poly<T>::operator/(T constant) const {
	return poly<T>(*this) /= constant;
}

template <typename T>
const poly<T> operator*(T constant, const poly<T> &p) {
	return p*constant;
}

template <typename T>
const poly<T> poly<T>::operator-() const {
	poly<T> neg();
	for (int i = 0; i < this->coeffs.size(); i++)
		neg.coeffs.push_back(-this->coeffs[i]);
	return neg;
}

template <typename T>
poly<T> &poly<T>::operator+=(const poly<T> &p) {
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
const poly<T> poly<T>::operator+(const poly<T> &p) const {
	return poly<T>(*this) += p;
}

template <typename T>
const poly<T> poly<T>::operator-(const poly<T> &p) const {
	return poly<T>(*this) -= p;
}

template <typename T>
const poly<T> poly<T>::operator*(const poly<T> &p) const {
	poly<T> ret = poly<T>();
	for (int i = 0; i <= p.degree() + this->degree(); i++)
		ret.coeffs.push_back(static_cast<T>(0));
	for (int i = 0; i <= p.degree(); i++) {
		for (int j = 0; j <= this->degree(); j++) {
			ret.coeffs[i+j] += this->coeffs[j]*p.coeffs[i];	
		}
	}
	return ret;
}

template <typename T>
poly<T> &poly<T>::operator*=(const poly<T> &p) {
	*this = (*this) * p;
}

template <typename T>
poly<T> &poly<T>::operator<<=(unsigned int len) {
	for (int i = 0; i < len; i++)
		this->coeffs.insert(this->coeffs.begin(), static_cast<T>(0));
	return *this;
}

template <typename T>
poly<T> &poly<T>::operator>>=(unsigned int len) {
	for (int i = 0; i < len; i++)
		this->coeffs.erase(this->coeffs.begin());
	return *this;
} 

template <typename T>
const qr_pair<poly<T>> poly<T>::divide(const poly<T> &other) const {
	// Algorithm 3.1.1

	poly<T> r = *this, q;
	T invlb = static_cast<T>(1)/(other[other.degree()]);
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
const poly<T> poly<T>::operator/(const poly<T> &p) const {
	return this->divide(p).quotient;
}

template <typename T>
const poly<T> poly<T>::operator%(const poly<T> &p) const {
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
	bool first = true;
	for (int i = 0; i < p.coeffs.size(); i++) {
		if (!first)
			os << " + ";
		if (p.coeffs[i] == static_cast<T>(0))
			continue;
		first = false;
		os << p.coeffs[i];
		if (i == 1)
			os << "x";
		else if (i > 1)
			os << "x^" << i;
	}
	if (first)
		os << "0";
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
const poly<U> poly<T>::convert(std::function<U(T)> converter) const {
	std::vector<U> vec;
	for (int i = 0; i < this->coeffs.size(); i++)
		vec.push_back(converter(this->coeffs[i]));
	return poly<U>(vec);
}