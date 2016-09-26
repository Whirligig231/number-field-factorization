#include <iostream>
#include <gmp.h>
#include <gmpxx.h>

#include "polyring.h"
#include "modring.h"
#include "alg.h"
#include "typedefs.h"

int main(int argc, char *argv[]) {
	Z_X a({-12, 2, 2});
	Z_X b({36, 24, 4});
	std::cout << sub_resultant_gcd(a, b) << std::endl;
	
	return 0;
}