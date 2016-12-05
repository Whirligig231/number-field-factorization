#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <ctime>
#include <cstdlib>

#include "polyring.h"
#include "modring.h"
#include "complex.h"
#include "polymodring.h"
#include "numberfield.h"
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
	mpf_set_default_prec(1000);
	
	Q_X p({Q("-10000000000000000000000000000000000000001/10000000000000000000000000000000000000000 "), Q("-999999999999999999970000000000000000000099999999999999999999/10000000000000000000000000000000000000000"), Q("299999999999999999996999999999999999999999999999999999999999999999999999999999999/1000000000000000000000000000000000000000000000000000000000000"), Q("-299999999999999999998999999999999999999969999999999999999999999999999999999999999/1000000000000000000000000000000000000000000000000000000000000"), Q("9999999999999999999999999999999999999997/100000000000000000000"), Q("1/100000000000000000000")});
	std::vector<C> roots = find_complex_roots(p, 50);
	
	for (int i = 0; i < roots.size(); i++)
		std::cout << roots[i] << std::endl;

	return 0;
}