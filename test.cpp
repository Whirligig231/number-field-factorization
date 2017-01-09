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

C eval_numberfield(poly<numberfield> poly, std::vector<C> alphas, unsigned int level) {
	C_X poly_C;
	if (level > 0) {
		// Convert each coefficient
		std::vector<C> coeffs;
		for (unsigned int i = 0; i <= poly.degree(); i++) {
			numberfield nf = poly[i];
			coeffs.push_back(eval_numberfield(nf.get_poly_value().get_value(), alphas, level - 1));
		}
		poly_C = C_X(coeffs);
	}
	else {
		std::vector<C> coeffs;
		for (unsigned int i = 0; i <= poly.degree(); i++) {
			numberfield nf = poly[i];
			coeffs.push_back(mpf_class(nf.get_rational_value()));
		}
		poly_C = C_X(coeffs);
	}
	
	return poly_C.evaluate(alphas[level]);
}

int main(int argc, char *argv[]) {
	mpf_set_default_prec(1000);
	
	
	
	
	// numberfield n1(3);
	// poly<numberfield> base({numberfield(1), numberfield(0), numberfield(1)});
	// poly<numberfield> value({numberfield(2), numberfield(4)});
	// polymod<numberfield> pm(base, value);
	// numberfield n2(pm);
	// numberfield n3 = n1 + n2;
	
	
	
	
	
	
	
	
	
	
	
	// return 0;
	
	// Precisions
	unsigned int k1 = 10, k2 = 10;
	R k2_real = 1;
	for (int i = 0; i < k2; i++)
		k2_real *= 10;
	
	std::vector<poly<numberfield>> min_polys;
	std::vector<C_X> orig_polys;
	std::vector<C> alphas;
	
	while (true) {
		// Get a polynomial from the user
		Q_X current_poly({2, 0, 1});
		C_X orig_poly(current_poly);
		orig_polys.push_back(orig_poly);
		// Compute the roots
		std::vector<C> roots = find_complex_roots(current_poly, k1);
		// Display the roots
		for (unsigned int i = 0; i < roots.size(); i++)
			std::cout << roots[i] << std::endl;
		// Pick a root
		unsigned int chosen_i = 0;
		alphas.push_back(roots[chosen_i]);
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
		// Factor the current polynomial to find the minimal polynomial
		std::vector<poly<numberfield>> factored = factor(current_poly_conv);
		
		// Test the factors to see which is closest to 0
		if (factored.size() == 1) {
			min_polys.push_back(factored[0]);
		}
		else {
			while (true) {
				// Evaluate each factor
				std::vector<C> values;
				for (unsigned int i = 0; i < factored.size(); i++) {
					values.push_back(eval_numberfield(factored[i], alphas, alphas.size() - 1));
				}
				
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
					// Increase precision and re-run
					for (unsigned int i = 0; i < alphas.size(); i++) {
						// Newton's method
						alphas[i] -= orig_polys[i].evaluate(alphas[i])/orig_polys[i].derivative().evaluate(alphas[i]);
					}
				}
			}
			std::cout << min_polys[min_polys.size() - 1] << std::endl;
		}
		
		// break;
	}

	return 0;
}