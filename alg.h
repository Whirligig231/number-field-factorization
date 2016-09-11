#include <gmp.h>
#include <gmpxx.h>
#include <utility>
#include "polyring.h"
#include "modring.h"
#include "typedefs.h"

std::pair<Z_X, Z_X> hensel_lift(Z p, Z q, Z_X a, Z_X b, Z_X c, Z_X u, Z_X v);
std::pair<Z_X, Z_X> quad_hensel_lift(Z p, Z q, Z_X a1, Z_X b1, Z_X u, Z_X v);