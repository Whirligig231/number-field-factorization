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
#include "polyio.h"

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

C eval_numberfield(poly<numberfield> poly, std::vector<C> alphas, unsigned int level) {
	C_X poly_C;
	if (level > 0) {
		// Convert each coefficient
		std::vector<C> coeffs;
		for (int i = 0; i <= poly.degree(); i++) {
			numberfield nf = poly[i];
			coeffs.push_back(eval_numberfield(nf.get_poly_value().get_value(), alphas, level - 1));
		}
		poly_C = C_X(coeffs);
	}
	else {
		std::vector<C> coeffs;
		for (int i = 0; i <= poly.degree(); i++) {
			numberfield nf = poly[i];
			coeffs.push_back(mpf_class(nf.get_rational_value()));
		}
		poly_C = C_X(coeffs);
	}
	return poly_C.evaluate(alphas[level]);
}

int main(int argc, char *argv[]) {
	mpf_set_default_prec(1000);
	
	// Precisions
	unsigned int k1 = 10, k2 = 10;
	R k2_real = 1;
	for (int i = 0; i < k2; i++)
		k2_real *= 10;
	
	std::vector<poly<numberfield>> min_polys;
	std::vector<C_X> orig_polys;
	std::vector<C> alphas;
	
	int BLAH = 1;
	while (true) {
		// std::cout << "step 1: getting polynomial" << std::endl;
		// Get a polynomial from the user
		std::vector<Q> cpc;
		cpc.push_back(-2);
		for (unsigned int i = 0; i < BLAH; i++)
			cpc.push_back(0);
		cpc.push_back(1);
		Q_X current_poly(cpc);
		BLAH++;
		C_X orig_poly(current_poly);
		orig_polys.push_back(orig_poly);
		// Compute the roots
		// std::cout << "step 2: finding roots" << std::endl;
		std::vector<C> roots = find_complex_roots(current_poly, k1);
		// Display the roots
		for (unsigned int i = 0; i < roots.size(); i++)
			std::cout << roots[i] << std::endl;
		// Pick a root
		unsigned int chosen_i = 0;
		alphas.push_back(roots[chosen_i]);
		// std::cout << "step 3: converting polynomial to extended field" << std::endl;
		// Convert the polynomial in terms of the previous polynomials
		poly<numberfield> current_poly_conv = current_poly;
		for (unsigned int i = 0; i < min_polys.size(); i++) {
			std::vector<numberfield> coeffs;
			for (unsigned int j = 0; j <= current_poly_conv.degree(); j++) {
				numberfield this_coeff = current_poly_conv[j];
				this_coeff = numberfield(polymod<numberfield>(min_polys[i], poly<numberfield>(this_coeff)));
				coeffs.push_back(this_coeff);
			}
			current_poly_conv = poly<numberfield>(coeffs);
		}
		// std::cout << "step 4: factoring the polynomial" << std::endl;
		// Factor the current polynomial to find the minimal polynomial
		std::vector<poly<numberfield>> factored = factor(current_poly_conv);
		
		// Test the factors to see which is closest to 0
		if (factored.size() == 1) {
			min_polys.push_back(factored[0]);
		}
		else {
			while (true) {
				// std::cout << "step 5: evaluating each factor" << std::endl;
				// Evaluate each factor
				std::vector<C> values;
				for (unsigned int i = 0; i < factored.size(); i++) {
					values.push_back(eval_numberfield(factored[i], alphas, alphas.size() - 1));
				}
				
				// std::cout << "step 6: testing for the right factor" << std::endl;
				R min_nm = -1, max_nm = -1;
				unsigned int min_i = 0;
				for (unsigned int i = 0; i < values.size(); i++) {
					R nm = values[i].norm();
					if (nm < min_nm || min_nm < 0) {
						min_i = i;
						min_nm = nm;
					}
					if (nm > max_nm)
						max_nm = nm;
				}
				
				if (min_nm == 0 || max_nm / min_nm > k2_real) {
					// Factor i is the right factor!
					min_polys.push_back(factored[min_i]);
					break;
				}
				else {
					// std::cout << "step 7: running Newton's method" << std::endl;
					// Increase precision and re-run
					for (unsigned int i = 0; i < alphas.size(); i++) {
						// Newton's method
						alphas[i] -= orig_polys[i].evaluate(alphas[i])/orig_polys[i].derivative().evaluate(alphas[i]);
					}
				}
			}
			// std::cout << "step 8: printing output" << std::endl;
			// std::cout << min_polys[min_polys.size() - 1] << std::endl;
			print_polyterm_list(&std::cout, get_polyterm_list(min_polys[min_polys.size() - 1]));
			std::cout << std::endl;
		}
		
		// break;
	}

	return 0;
}