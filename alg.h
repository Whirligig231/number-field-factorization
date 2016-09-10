#include <gmp.h>
#include <gmpxx.h>
#include <utility>
#include "polyring.h"
#include "modring.h"
#include "macros.h"

std::pair<Z_X, Z_X> hensel_lift(Z p, Z q, Z_X c, Z_X a, Z_X b, Z_X u, Z_X v);