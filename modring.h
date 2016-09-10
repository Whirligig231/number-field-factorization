#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include <initializer_list>
#include <iostream>

class mod {
	private:
		mpz_class base, value;
	public:
		mod(mpz_class base);
		mod(mpz_class base, mpz_class value);
		mod(const mod &other);
		mod(const mod &other, mpz_class value);
		
		mpz_class getBase() const;
		
		mod &operator=(mpz_class value);
		mod &operator=(const mod &other);
		
		bool operator==(const mod &other) const;
		bool operator!=(const mod &other) const;
		
		mod &operator+=(const mod &other);
		mod &operator-=(const mod &other);
		mod &operator*=(const mod &other);
		mod &operator/=(const mod &other);
		
		const mod operator+(const mod &other) const;
		const mod operator-(const mod &other) const;
		const mod operator*(const mod &other) const;
		const mod operator/(const mod &other) const;
		
		const mod operator-() const;
		const mod inv() const;
		
		friend std::ostream &operator<<(std::ostream &os, const mod &p);
		
		operator mpz_class();
};

mod::mod(mpz_class base) {
	this->base = base;
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
	mpz_mod(this->value.get_mpz_t(), this->value.get_mpz_t(), this->base.get_mpz_t());
	return *this;
}

mod &mod::operator-=(const mod &other) {
	this->value -= other.value;
	mpz_mod(this->value.get_mpz_t(), this->value.get_mpz_t(), this->base.get_mpz_t());
	return *this;
}

mod &mod::operator*=(const mod &other) {
	this->value *= other.value;
	mpz_mod(this->value.get_mpz_t(), this->value.get_mpz_t(), this->base.get_mpz_t());
	return *this;
}

mod &mod::operator/=(const mod &other) {
	return (*this) *= other.inv();
}

const mod mod::operator+(const mod &other) const {
	return mod(*this) += other;
}

const mod mod::operator-(const mod &other) const {
	return mod(*this) -= other;
}

const mod mod::operator*(const mod &other) const {
	return mod(*this) *= other;
}

const mod mod::operator/(const mod &other) const {
	return mod(*this) /= other;
}

const mod mod::operator-() const {
	return mod(this->base, this->base - this->value);
}

const mod mod::inv() const {
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