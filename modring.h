#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <Eigen/Core>

#include "numbers.h"

#pragma once

class mod {
	private:
		mpz_class base, value;
		
	public:
		mod();
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
	return mod(reference.get_base(), 0);
}

template <>
inline mod one<mod>(const mod &reference) {
	return mod(reference.get_base(), 1);
}

template <>
inline mod from_int<mod>(int n, const mod &reference) {
	return mod(reference.get_base(), n);
}

namespace Eigen {
	
	template<>
	struct NumTraits<mod> : GenericNumTraits<mod> {
		typedef mod Real;
		typedef mod NonInteger;
		typedef mod Literal;
		typedef mod Nested;
		
		static inline Real epsilon() { return mod(1, 0); }
		static inline Real dummy_precision() { return mod(1, 0); }
		static inline Real digits10() { return mod(1, 0); }
		
		enum {
			IsComplex = 0,
			IsInteger = 1,
			ReadCost = 1,
			AddCost = 5,
			MulCost = 5,
			IsSigned = 1,
			RequireInitialization = 1
		};
	};
	
}