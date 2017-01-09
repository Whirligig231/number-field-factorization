#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <Eigen/Core>

#include "numbers.h"
#include "polyring.h"
#include "modring.h"
#include "polymodring.h"

#pragma once

class numberfield {
	private:
		std::unique_ptr<mpq_class> rational_value;
		std::unique_ptr<polymod<numberfield>> poly_value;
		unsigned int poly_levels;
		
		void lift(const numberfield &other);
		
	public:
		numberfield();
		numberfield(const mpq_class &value);
		numberfield(const polymod<numberfield> &value);
		numberfield(const numberfield &other);

		unsigned int get_poly_levels() const;
		mpq_class get_rational_value() const;
		polymod<numberfield> get_poly_value() const;

		numberfield &operator=(const mpq_class &value);
		numberfield &operator=(const polymod<numberfield> &value);
		numberfield &operator=(const numberfield &other);
		
		bool operator==(const numberfield &other) const;
		bool operator!=(const numberfield &other) const;
		
		numberfield &operator+=(const numberfield &other);
		numberfield &operator-=(const numberfield &other);
		numberfield &operator*=(const numberfield &other);
		numberfield &operator/=(const numberfield &other);
		
		numberfield operator+(const numberfield &other) const;
		numberfield operator-(const numberfield &other) const;
		numberfield operator*(const numberfield &other) const;
		numberfield operator/(const numberfield &other) const;
		
		numberfield operator-() const;
		numberfield inv() const;
		
		friend std::ostream &operator<<(std::ostream &os, const numberfield &p);
};

template <>
class util<numberfield> {
public:
	static numberfield zero();
	static numberfield zero(const numberfield &reference);
	static numberfield one();
	static numberfield one(const numberfield &reference);
	static numberfield from_int(int n, const numberfield &reference);
	static numberfield get_gcd(numberfield a, numberfield b);
};