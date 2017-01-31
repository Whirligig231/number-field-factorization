#include <iostream>
#include <string>
#include <sstream>
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

struct option {
	bool valid;
	std::string name;
	unsigned int value;
};

option parse_option(std::string opt) {
	option ret;
	ret.valid = false;
	if (opt[0] != '-')
		return ret;
	unsigned int equ = opt.find("=");
	if (equ == std::string::npos)
		return ret;
	ret.name = opt.substr(1, equ - 1);
	std::stringstream ss(opt.substr(equ + 1, opt.length() - equ - 1));
	ss >> ret.value;
	if (ss.fail())
		return ret;
	ret.valid = true;
	return ret;
}

int main(int argc, char *argv[]) {
	mpf_set_default_prec(1000);
	
	// Precisions
	unsigned int k1 = 10, k2 = 10;
	
	for (unsigned int i = 1; i < argc; i++) {
		std::string this_arg(argv[i]);
		if (this_arg.length() == 0)
			continue;
		option opt = parse_option(this_arg);
		if (!opt.valid) {
			std::cout << "Unknown command line option " << this_arg << std::endl;
		}
		else if (opt.name == "k1") {
			std::cout << "Setting k1 to " << opt.value << std::endl;
			k1 = opt.value;
		}
		else if (opt.name == "k2") {
			std::cout << "Setting k2 to " << opt.value << std::endl;
			k2 = opt.value;
		}
		else {
			std::cout << "Unknown command line argument " << opt.name << std::endl;
		}
	}
	
	R k2_real = 1;
	for (int i = 0; i < k2; i++)
		k2_real *= 10;
	
	std::vector<poly<numberfield>> min_polys;
	std::vector<C_X> orig_polys;
	std::vector<C> alphas;
	
	std::string field_str = "Q";
	
	while (true) {
		char alpha_char = (char)('a' + alphas.size());
		// std::cout << "step 1: getting polynomial" << std::endl;
		// Get a polynomial from the user
		std::cout << "Enter a rational polynomial in " << alpha_char << ", or nothing to stop adjoining roots:" << std::endl;
		std::string input;
		std::getline(std::cin, input);
		if (input.size() == 0)
			break;
		Q_X current_poly = get_rational_poly(standardize(scan_polyterm_list(input)));
		C_X orig_poly(current_poly);
		orig_polys.push_back(orig_poly);
		// Compute the roots
		// std::cout << "step 2: finding roots" << std::endl;
		std::vector<C> roots = find_complex_roots(current_poly, k1);
		// Display the roots
		for (unsigned int i = 0; i < roots.size(); i++)
			std::cout << "Root " << (i + 1) << " ~= " << roots[i] << std::endl;
		// Pick a root
		int chosen_i = -1;
		std::cout << "Which root is the value of " << alpha_char << "?" << std::endl;
		do {
			std::getline(std::cin, input);
			std::stringstream(input) >> chosen_i;
			if (chosen_i <= 0 || chosen_i > roots.size()) {
				chosen_i = -1;
				std::cout << "Please enter a valid value from 1 to " << roots.size() << "." << std::endl;
			}
			else
				chosen_i--;
		} while (chosen_i < 0);
		alphas.push_back(roots[chosen_i]);
		std::cout << "Converting to a polynomial in " << field_str << " ... " << std::endl;
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
		std::cout << "Factoring in " << field_str << " ... " << std::endl;
		// Factor the current polynomial to find the minimal polynomial
		std::vector<poly<numberfield>> factored = factor(current_poly_conv);
		
		// Test the factors to see which is closest to 0
		if (factored.size() == 1) {
			min_polys.push_back(factored[0]);
		}
		else {
			while (true) {
				std::cout << "Evaluating each factor in C ... " << std::endl;
				// Evaluate each factor
				std::vector<C> values;
				for (unsigned int i = 0; i < factored.size(); i++) {
					values.push_back(eval_numberfield(factored[i], alphas, alphas.size() - 1));
				}
				
				std::cout << "Comparing factor values ... " << std::endl;
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
					std::cout << "Refining the roots ... " << std::endl;
					// Increase precision and re-run
					for (unsigned int i = 0; i < alphas.size(); i++) {
						// Newton's method
						alphas[i] -= orig_polys[i].evaluate(alphas[i])/orig_polys[i].derivative().evaluate(alphas[i]);
					}
				}
			}
			// std::cout << "step 8: printing output" << std::endl;
			// std::cout << min_polys[min_polys.size() - 1] << std::endl;
		}
		
		std::cout << "The minimal polynomial of " << alpha_char << " over " << field_str << " is ";
		print_polyterm_list(&std::cout, get_polyterm_list(min_polys[min_polys.size() - 1]));
		std::cout << "." << std::endl;
		
		std::stringstream o;
		o << "[";
		o << alpha_char;
		o << "]";
		field_str += o.str();
	}

	char alpha_char = (char)('a' + alphas.size());
	std::cout << "Enter a polynomial in " << alpha_char << " and all previous roots:" << std::endl;
	std::string input;
	std::getline(std::cin, input);
	poly<numberfield> poly_to_factor = get_numberfield_poly(standardize(scan_polyterm_list(input)), min_polys);
	std::cout << "Factoring ... " << std::endl;
	// Factor the current polynomial to find the minimal polynomial
	std::vector<poly<numberfield>> factored = factor(poly_to_factor);
	print_polyterm_list(&std::cout, get_polyterm_list(poly_to_factor));
	std::cout << " = ";
	for (unsigned int i = 0; i < factored.size(); i++) {
		std::cout << '(';
		print_polyterm_list(&std::cout, get_polyterm_list(factored[i]));
		std::cout << ')';
	}
	std::cout << std::endl;
	
	return 0;
}