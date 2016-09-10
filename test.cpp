#include <iostream>
#include <gmp.h>
#include <gmpxx.h>

#include "polyring.h"
#include "modring.h"
#include "alg.h"

int main(int argc, char *argv[]) {
	poly<mod> mod5({mod(5, 3), mod(5, 4), mod(5, 19)});
	std::cout << mod5 << std::endl;
	mod5 += mod(5, 4);
	std::cout << mod5 << std::endl;
	poly<mpz_class> lift = static_cast<poly<mpz_class>>(mod5);
	std::cout << lift << std::endl;
	lift += mpz_class(4);
	std::cout << lift << std::endl;
	poly<mod> project = lift.convert(to_mod(3));
	std::cout << project << std::endl;
	return 0;
}