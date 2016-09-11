#include <iostream>
#include <gmp.h>
#include <gmpxx.h>

#include "polyring.h"
#include "modring.h"
#include "alg.h"
#include "typedefs.h"

int main(int argc, char *argv[]) {
	Q_X a({12, -13, 0, 1});
	Q_X b({-1440, 564, 8, -15, 1});

	std::tuple<Q_X, Q_X, Q_X> uvd = extended_gcd(a, b);
	std::cout << "u = " << std::get<0>(uvd) << std::endl;
	std::cout << "v = " << std::get<1>(uvd) << std::endl;
	std::cout << "d = " << std::get<2>(uvd) << std::endl;
	
	std::cout << "a % d = " << a % std::get<2>(uvd) << std::endl;
	std::cout << "b % d = " << b % std::get<2>(uvd) << std::endl;
	std::cout << "au + bv = " << a*std::get<0>(uvd) + b*std::get<1>(uvd) << std::endl;
	
	return 0;
}