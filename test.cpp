#include <iostream>
#include <gmp.h>
#include <gmpxx.h>

#include "polyring.h"
#include "modring.h"
#include "alg.h"
#include "typedefs.h"

int main(int argc, char *argv[]) {
	Z_X a({5, 1, 10, 6, 2, 2});
	std::cout << a << std::endl;
	ZN_X a2 = a.convert(to_mod(11));
	std::cout << a2 << std::endl;
	
	std::vector<ZN_X> bsp = berlekamp(a2);
	std::cout << "Factors:" << std::endl;
	for (int i = 0; i < bsp.size(); i++)
		std::cout << bsp[i] << std::endl;
	
	return 0;
}