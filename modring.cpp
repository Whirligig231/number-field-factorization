#include "modring.h"

mod::mod(mpz_class value) {
	this->base = 0;
	this->value = value;
}

mod::mod(mpz_class base, mpz_class value) {
	this->base = base;
	mpz_mod(this->value.get_mpz_t(), value.get_mpz_t(), base.get_mpz_t());
}

mod::mod(const mod &other) {
	this->base = other.base;
	this->value = other.value;
}

mod::mod(const mod &other, mpz_class value) {
	this->base = other.base;
	mpz_mod(this->value.get_mpz_t(), value.get_mpz_t(), other.base.get_mpz_t());
}

mpz_class mod::get_base() const {
	return this->base;
}

mod &mod::operator=(mpz_class value) {
	mpz_mod(this->value.get_mpz_t(), value.get_mpz_t(), this->base.get_mpz_t());
	return *this;
}

mod &mod::operator=(const mod &other) {
	this->base = other.base;
	this->value = other.value;
	return *this;
}

bool mod::operator==(const mod &other) const {
	return (this->value == other.value);
}

bool mod::operator!=(const mod &other) const {
	return (this->value != other.value);
}

mod &mod::operator+=(const mod &other) {
	this->value += other.value;
	mpz_class base = (this->base == 0) ? other.base : this->base;
	mpz_mod(this->value.get_mpz_t(), this->value.get_mpz_t(), base.get_mpz_t());
	return *this;
}

mod &mod::operator-=(const mod &other) {
	this->value -= other.value;
	mpz_class base = (this->base == 0) ? other.base : this->base;
	mpz_mod(this->value.get_mpz_t(), this->value.get_mpz_t(), base.get_mpz_t());
	return *this;
}

mod &mod::operator*=(const mod &other) {
	this->value *= other.value;
	mpz_class base = (this->base == 0) ? other.base : this->base;
	mpz_mod(this->value.get_mpz_t(), this->value.get_mpz_t(), base.get_mpz_t());
	return *this;
}

mod &mod::operator/=(const mod &other) {
	return (*this) *= other.inv();
}

mod mod::operator+(const mod &other) const {
	return mod(*this) += other;
}

mod mod::operator-(const mod &other) const {
	return mod(*this) -= other;
}

mod mod::operator*(const mod &other) const {
	return mod(*this) *= other;
}

mod mod::operator/(const mod &other) const {
	return mod(*this) /= other;
}

mod mod::operator-() const {
	return mod(this->base, this->base - this->value);
}

mod mod::inv() const {
	mpz_class g, s;
	mpz_gcdext(g.get_mpz_t(), s.get_mpz_t(), NULL, this->value.get_mpz_t(), this->base.get_mpz_t());
	if (g != 1) {
		std::cout << "ERROR: attempt to invert non-invertible element" << std::endl;
	}
	mpz_mod(s.get_mpz_t(), s.get_mpz_t(), this->base.get_mpz_t());
	return mod(this->base, s);
}

std::ostream &operator<<(std::ostream &os, const mod &m) {
	return os << m.value;
}

mod::operator mpz_class() {
	return this->value;
}

mod make_mod(mpz_class base, mpz_class value) {
	return mod(base, value);
}

std::function<mod(mpz_class)> to_mod(mpz_class base) {
	return std::bind(make_mod, base, std::placeholders::_1);
}