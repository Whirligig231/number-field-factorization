#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include <functional>
#include <initializer_list>
#include <iostream>

#include "numbers.h"

#pragma once

class complex {
	private:
		mpf_class real, imag;
		
	public:
		complex();
		complex(int real);
		complex(mpf_class real);
		complex(mpf_class real, mpf_class imag);
		complex(const complex &other);
		
		mpf_class get_real() const;
		mpf_class get_imag() const;
		
		complex &operator=(mpf_class value);
		complex &operator=(const complex &other);
		
		bool operator==(const complex &other) const;
		bool operator!=(const complex &other) const;
		
		complex &operator+=(const complex &other);
		complex &operator-=(const complex &other);
		complex &operator*=(const complex &other);
		complex &operator/=(const complex &other);
		
		complex operator+(const complex &other) const;
		complex operator-(const complex &other) const;
		complex operator*(const complex &other) const;
		complex operator/(const complex &other) const;
		
		complex operator-() const;
		complex inv() const;
		
		complex conjugate() const;
		mpf_class norm() const;
		
		friend std::ostream &operator<<(std::ostream &os, const complex &p);

};

template <>
class util<complex> {
public:
	static complex zero();
	static complex zero(const complex &reference);
	static complex one(const complex &reference);
	static complex from_int(int n, const complex &reference);
};