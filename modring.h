#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include <functional>
#include <initializer_list>
#include <iostream>

#include "numbers.h"

#pragma once

class mod {
	private:
		mpz_class base, value;
		
	public:
		explicit mod(mpz_class value);
		mod(mpz_class base, mpz_class value);
		mod(const mod &other);
		mod(const mod &other, mpz_class value);
		
		mpz_class get_base() const;
		
		mod &operator=(mpz_class value);
		mod &operator=(const mod &other);
		
		bool operator==(const mod &other) const;
		bool operator!=(const mod &other) const;
		
		mod &operator+=(const mod &other);
		mod &operator-=(const mod &other);
		mod &operator*=(const mod &other);
		mod &operator/=(const mod &other);
		
		mod operator+(const mod &other) const;
		mod operator-(const mod &other) const;
		mod operator*(const mod &other) const;
		mod operator/(const mod &other) const;
		
		mod operator-() const;
		mod inv() const;
		
		friend std::ostream &operator<<(std::ostream &os, const mod &p);
		
		operator mpz_class();
};

std::function<mod(mpz_class)> to_mod(mpz_class base);

template <>
inline mod zero<mod>(const mod &reference) {
	return mod(0, reference.get_base());
}

template <>
inline mod one<mod>(const mod &reference) {
	return mod(1, reference.get_base());
}