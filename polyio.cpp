#include "polyio.h"

polyterm::polyterm() {
	this->coeff = 0;
}

polyterm::polyterm(std::string in) {
	unsigned int end_of_num;
	for (end_of_num = 0; end_of_num < in.size(); end_of_num++) {
		if (in[end_of_num] >= 'a')
			break;
	}
	std::string coeff_str = in.substr(0, end_of_num);
	if (coeff_str.size() == 0)
		this->coeff = 1;
	else if (coeff_str.size() == 1 && coeff_str[0] == '+')
		this->coeff = 1;
	else if (coeff_str.size() == 1 && coeff_str[0] == '-')
		this->coeff = -1;
	else if (coeff_str[0] == '+')
		this->coeff = Q(coeff_str.substr(1));
	else
		this->coeff = Q(coeff_str);
	
	unsigned int f_start = end_of_num, f_end = f_start;
	while (f_end < in.size()) {
		f_start = f_end;
		for (f_end++; f_end < in.size(); f_end++)
			if (in[f_end] >= 'a')
				break;

		unsigned int powind = (unsigned int)(in[f_start] - 'a');
		if (powind >= this->powers.size())
			this->powers.resize(powind + 1, 0);
		
		if (f_start + 1 == f_end) {
			this->powers[powind] += 1;
		}
		else {
			unsigned int k;
			std::stringstream(in.substr(f_start + 2, f_end - f_start - 2)) >> k;
			this->powers[powind] += k;
		}
	}
}

bool polyterm::operator<(const polyterm &other) const {
	unsigned int sum1 = 0, sum2 = 0;
	for (unsigned int i = std::max(this->powers.size(), other.powers.size()) - 1; i < (1 << 31); i--) {
		unsigned int p1 = 0, p2 = 0;
		if (i < this->powers.size())
			p1 = this->powers[i];
		if (i < other.powers.size())
			p2 = other.powers[i];
		
		if (p1 > p2)
			return 1;
		if (p1 < p2)
			return 0;
	}

	return this->coeff > other.coeff;
}

void polyterm::print(std::ostream *out) {
	bool has_power = 0;
	for (unsigned int i = 0; i < this->powers.size(); i++)
		if (this->powers[i] > 0)
			has_power = 1;
	if ((this->coeff != 1 && this->coeff != -1) || !has_power)
		*out << this->coeff;
	if (this->coeff == -1 && has_power)
		*out << '-';
	bool first_power = 1;
	for (unsigned int i = 0; i < this->powers.size(); i++) {
		if (this->powers[i] == 0)
			continue;
		if (!first_power)
			*out << ' ';
		first_power = 0;
		*out << (char)('a' + i);
		if (this->powers[i] > 1)
			*out << '^' << this->powers[i];
	}
}

void polyterm::print_added(std::ostream *out) {
	if (this->coeff == 0)
		return;
	if (this->coeff > 0)
		*out << " +";
	else
		*out << " -";
	bool has_power = 0;
	for (unsigned int i = 0; i < this->powers.size(); i++)
		if (this->powers[i] > 0)
			has_power = 1;
	if ((this->coeff != 1 && this->coeff != -1) || !has_power)
		*out << ' ' << abs(this->coeff);
	bool first_power = 1;
	for (unsigned int i = 0; i < this->powers.size(); i++) {
		if (this->powers[i] == 0)
			continue;
		if (!first_power || this->coeff == 1 || this->coeff == -1)
			*out << ' ';
		first_power = 0;
		*out << (char)('a' + i);
		if (this->powers[i] > 1)
			*out << '^' << this->powers[i];
	}
}

std::vector<polyterm> scan_polyterm_list(std::string input) {
	// Remove spaces
	input.erase(std::remove_if(input.begin(), input.end(), [](char x){return std::isspace(x);}), input.end());
	
	unsigned int start_of_term = 0, end_of_term = 0;
	std::vector<polyterm> ret;
	while (true) {
		start_of_term = end_of_term;
		end_of_term++;
		while (end_of_term < input.size() && input[end_of_term] != '+' && input[end_of_term] != '-')
			end_of_term++;
		ret.push_back(polyterm(input.substr(start_of_term, end_of_term - start_of_term)));
		if (end_of_term >= input.size())
			return ret;
	}
}

void print_polyterm_list(std::ostream *out, std::vector<polyterm> list) {
	if (list.size() == 0) {
		*out << '0';
		return;
	}
	
	list[0].print(out);
	for (unsigned int i = 1; i < list.size(); i++)
		list[i].print_added(out);
}

std::vector<polyterm> standardize(std::vector<polyterm> terms) {
	std::sort(terms.begin(), terms.end());
	std::vector<polyterm> output;
	if (terms.size() == 0)
		return output;
	output.push_back(terms[0]);
	for (unsigned int i = 1; i < terms.size(); i++) {
		bool equal = 1;
		for (unsigned int j = 0; j < output[output.size() - 1].powers.size() || j < terms[i].powers.size(); j++) {
			unsigned int p1, p2;
			
			if (j >= output[output.size() - 1].powers.size())
				p1 = 0;
			else
				p1 = output[output.size() - 1].powers[j];
			
			if (j >= terms[i].powers.size())
				p2 = 0;
			else
				p2 = terms[i].powers[j];
			
			if (p1 != p2) {
				equal = 0;
				break;
			}
		}
		
		if (equal) {
			output[output.size() - 1].coeff += terms[i].coeff;
		}
		else
			output.push_back(terms[i]);
	}
	
	return output;
}

std::vector<polyterm> get_polyterm_list(poly<numberfield> input) {
	std::vector<polyterm> terms;
	for (int i = 0; i <= input.degree(); i++) {
		numberfield n = input[i];
		if (n.get_poly_levels() == 0) {
			polyterm term;
			term.coeff = n.get_rational_value();
			term.powers = std::vector<unsigned int>();
			term.powers.push_back(i);
			terms.push_back(term);
		}
		else {
			std::vector<polyterm> new_terms = get_polyterm_list(n.get_poly_value().get_value());
			for (unsigned int j = 0; j < new_terms.size(); j++) {
				new_terms[j].powers.push_back(i);
				terms.push_back(new_terms[j]);
			}
		}
	}
	
	return standardize(terms);
}

poly<Q> get_rational_poly(std::vector<polyterm> terms) {
	std::vector<Q> coeffs;
	for (unsigned int i = 0; i < terms.size(); i++) {
		unsigned int powind = 0;
		for (unsigned int j = 0; j < terms[i].powers.size(); j++)
			powind += terms[i].powers[j];
		if (powind >= coeffs.size())
			coeffs.resize(powind + 1, 0);
		coeffs[powind] += terms[i].coeff;
	}
	return poly<Q>(coeffs);
}

poly<numberfield> get_numberfield_poly(std::vector<polyterm> terms, std::vector<poly<numberfield>> min_polys, int level) {
	if (level < 0)
		level = min_polys.size();
	if (level == 0) {
		// Polynomial over Q
		return poly<numberfield>(get_rational_poly(terms));
	}

	// Find max exponent of highest variable
	unsigned int max_exp = 0;
	for (unsigned int i = 0; i < terms.size(); i++) {
		unsigned int this_exp = 0;
		if (level < terms[i].powers.size())
			this_exp = terms[i].powers[level];
		if (this_exp > max_exp)
			max_exp = this_exp;
	}
	
	// Sort into bins for the various exponents
	std::vector<std::vector<polyterm>> bins;
	for (unsigned int i = 0; i <= max_exp; i++)
		bins.push_back(std::vector<polyterm>());
	
	for (unsigned int i = 0; i < terms.size(); i++) {
		unsigned int this_exp = 0;
		if (level < terms[i].powers.size())
			this_exp = terms[i].powers[level];
		if (level < terms[i].powers.size())
			terms[i].powers[level] = 0;
		bins[this_exp].push_back(terms[i]);
	}

	// Convert each bin using recursive call
	std::vector<poly<numberfield>> bin_results;
	for (unsigned int i = 0; i <= max_exp; i++)
		bin_results.push_back(get_numberfield_poly(bins[i], min_polys, level - 1));

	// Add min polynomials to make numberfield coefficients
	std::vector<numberfield> coeffs;
	for (unsigned int i = 0; i <= max_exp; i++)
		coeffs.push_back(numberfield(polymod<numberfield>(min_polys[level - 1], bin_results[i])));

	// Combine into one polynomial
	return poly<numberfield>(coeffs);
}