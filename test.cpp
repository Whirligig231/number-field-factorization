#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <ctime>
#include <cstdlib>

#include "polyring.h"
#include "modring.h"
#include "alg.h"
#include "typedefs.h"

ZN_X prand(Z p, int deg, gmp_randstate_t state) {
	std::vector<Z> coeffs;
	
	for (int i = 0; i <= deg; i++) {
		Z next;
		mpz_urandomm(next.get_mpz_t(), state, p.get_mpz_t());
		coeffs.push_back(next);
	}
	
	Z_X poly(coeffs);

	ZN_X p1 = poly.convert(to_mod(p));
	ZN_X p2 = poly.derivative().convert(to_mod(p));
	
	if (p1.degree() < deg)
		return prand(p, deg, state);

	p1 /= std::get<2>(extended_gcd(p1, p2));
	
	return p1;
}

int main(int argc, char *argv[]) {
	/*gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, time(NULL));
	
	int deg = atoi(argv[1]);
	Z p = atoi(argv[2]);
	int number = atoi(argv[3]);
	
	std::clock_t start = std::clock();
	for (int i = 0; i < number; i++) {
		ZN_X poly = prand(p, deg, state);
		berlekamp(poly);
	}
	std::cout << "Time per test: " << (std::clock() - start) * (double)1000 / (double)(CLOCKS_PER_SEC) / ((double) number) << " ms" << std::endl;*/
	
	Z_X a({0, 0, 2, 8, 12, 8, 2});
	std::vector<Z_X> factors = factor(a);
	for (int i = 0; i < factors.size(); i++)
		std::cout << factors[i] << std::endl;

	return 0;
}