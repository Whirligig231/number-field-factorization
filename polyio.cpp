#include "polyio.h"

bool polyterm::operator<(const polyterm &other) const {
	unsigned int sum1 = 0, sum2 = 0;
	for (unsigned int i = 0; i < this->powers.size(); i++)
		sum1 += this->powers[i];
	for (unsigned int i = 0; i < other.powers.size(); i++)
		sum2 += other.powers[i];
	
	if (sum1 > sum2)
		return 1;
	if (sum1 < sum2)
		return 0;
	
	for (unsigned int i = std::max(this->powers.size(), other.powers.size()) - 1; i >= 0; i--) {
		unsigned int p1, p2;
		if (i > this->powers.size())
			p1 = 0;
		else
			p1 = this->powers[i];
		
		if (i > other.powers.size())
			p2 = 0;
		else
			p2 = other.powers[i];
		
		if (p1 > p2)
			return 1;
		if (p2 < p1)
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
	for (unsigned int i = 0; i < this->powers.size(); i++) {
		if (this->powers[i] == 0)
			continue;
		*out << ' ' << (char)('a' + i);
		if (this->powers[i] > 1)
			*out << '^' << this->powers[i];
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