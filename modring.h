#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include <functional>
#include <initializer_list>
#include <iostream>

#pragma once

class mod {
	private:
		mpz_class base, value;
	public:
		explicit mod(mpz_class base);
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

std::function<mod(mpz_class)> to_mod(mpz_class base);