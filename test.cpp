#include <iostream>
#include <gmp.h>
#include <gmpxx.h>

#include "polyring.h"
#include "modring.h"
#include "alg.h"
#include "typedefs.h"

int main(int argc, char *argv[]) {
	Z_X a({1, 4, -3, 5});
	std::cout << a.derivative() << std::endl;
	
	return 0;
}