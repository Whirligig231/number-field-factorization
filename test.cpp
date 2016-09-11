#include <iostream>
#include <gmp.h>
#include <gmpxx.h>

#include "polyring.h"
#include "modring.h"
#include "alg.h"
#include "typedefs.h"

int main(int argc, char *argv[]) {
	Z_X a({3, 0, 1});
	Z_X b({4, -1, 1});
	Z_X c({2, 2, 2, 4, 1});
	Z_X u({0, 1});
	Z_X v({-1, -1});
	
	std::cout << "ab - c = " << (a*b - c) << std::endl;
	std::cout << "ua + vb = " << (u*a + v*b) << std::endl;
	
	std::pair<Z_X, Z_X> a1b1 = hensel_lift(5, 5, a, b, c, u, v);
	Z_X a1 = a1b1.first;
	Z_X b1 = a1b1.second;
	std::cout << "a1 = " << a1 << std::endl;
	std::cout << "b1 = " << b1 << std::endl;
	std::cout << "a1b1 - c = " << (a1*b1 - c) << std::endl;
	
	std::pair<Z_X, Z_X> u1v1 = quad_hensel_lift(5, 5, a1, b1, u, v);
	Z_X u1 = u1v1.first;
	Z_X v1 = u1v1.second;
	std::cout << "u1 = " << u1 << std::endl;
	std::cout << "v1 = " << v1 << std::endl;
	std::cout << "u1a1 + v1b1 = " << (u1*a1 + v1*b1) << std::endl;
	return 0;
}