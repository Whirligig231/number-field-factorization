#include "complex.h"

complex::complex() {
	this->real = 0;
	this->imag = 0;
}

complex::complex(int real) {
	this->imag = 0;
	this->real = mpf_class(real);
}

complex::complex(mpf_class real) {
	this->imag = 0;
	this->real = real;
}

complex::complex(mpf_class real, mpf_class imag) {
	this->real = real;
	this->imag = imag;
}

complex::complex(const complex &other) {
	this->real = other.real;
	this->imag = other.imag;
}

mpf_class complex::get_real() const {
	return this->real;
}

mpf_class complex::get_imag() const {
	return this->imag;
}

complex &complex::operator=(mpf_class value) {
	this->real = value;
	this->imag = 0;
	return *this;
}

complex &complex::operator=(const complex &other) {
	this->real = other.real;
	this->imag = other.imag;
	return *this;
}

bool complex::operator==(const complex &other) const {
	return (this->real == other.real) && (this->imag == other.imag);
}

bool complex::operator!=(const complex &other) const {
	return (this->real != other.real) || (this->imag != other.imag);
}

complex &complex::operator+=(const complex &other) {
	this->real += other.real;
	this->imag += other.imag;
	return *this;
}

complex &complex::operator-=(const complex &other) {
	this->real -= other.real;
	this->imag -= other.imag;
	return *this;
}

complex &complex::operator*=(const complex &other) {
	mpf_class temp = this->real;
	this->real = this->real * other.real - this->imag * other.imag;
	this->imag = temp * other.imag + this->imag * other.real;
	return *this;
}

complex &complex::operator/=(const complex &other) {
	return (*this) *= other.inv();
}

complex complex::operator+(const complex &other) const {
	return complex(*this) += other;
}

complex complex::operator-(const complex &other) const {
	return complex(*this) -= other;
}

complex complex::operator*(const complex &other) const {
	return complex(*this) *= other;
}

complex complex::operator/(const complex &other) const {
	return complex(*this) /= other;
}

complex complex::operator-() const {
	return complex(-this->real, -this->imag);
}

complex complex::inv() const {
	return complex(
		this->real / (this->real*this->real + this->imag*this->imag),
		-this->imag / (this->real*this->real + this->imag*this->imag)
	);
}

complex complex::conjugate() const {
	return complex(this->real, -this->imag);
}

mpf_class complex::norm() const {
	return ((*this) * this->conjugate()).get_real();
}

std::ostream &operator<<(std::ostream &os, const complex &m) {
	return os << m.real << " + " << m.imag << "i";
}

complex util<complex>::zero() {
	return complex(0, 0);
}

complex util<complex>::zero(const complex &reference) {
	return complex(0, 0);
}

complex util<complex>::one(const complex &reference) {
	return complex(1, 0);
}

complex util<complex>::from_int(int n, const complex &reference) {
	return complex(n, 0);
}