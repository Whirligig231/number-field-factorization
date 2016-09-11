#include <iostream>
#include <gmp.h>
#include <gmpxx.h>

#include "polyring.h"
#include "modring.h"
#include "alg.h"
#include "typedefs.h"

int main(int argc, char *argv[]) {
	mat<double> A(3, 3);
	A << 1.0, -1.0, 2.0, -1.0, 1.0, -2.0, 2.0, -2.0, 4.0;
	std::cout << A << std::endl;
	
	std::vector<vec<double>> basis = kernel(A);
	for (int i = 0; i < basis.size(); i++)
		std::cout << basis[i] << std::endl;
	
	std::cout << A << std::endl;
	
	return 0;
}