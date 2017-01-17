#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "polyring.h"
#include "modring.h"
#include "complex.h"
#include "polymodring.h"
#include "numberfield.h"
#include "alg.h"
#include "typedefs.h"

#pragma once

class polyterm {
public:
	Q coeff;
	std::vector<unsigned int> powers;
	
	polyterm();
	polyterm(std::string in);
	bool operator<(const polyterm &other) const;
	void print(std::ostream *out);
	void print_added(std::ostream *out);
};

std::vector<polyterm> scan_polyterm_list(std::string input);
void print_polyterm_list(std::ostream *out, std::vector<polyterm> list);
std::vector<polyterm> standardize(std::vector<polyterm> terms);
std::vector<polyterm> get_polyterm_list(poly<numberfield> input);