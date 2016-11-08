#include "numbers.h"

mpq_class util<mpq_class>::zero() {
	return mpq_class(0);
}

mpq_class util<mpq_class>::zero(const mpq_class &reference) {
	return mpq_class(0);
}

mpq_class util<mpq_class>::one() {
	return mpq_class(1);
}

mpq_class util<mpq_class>::one(const mpq_class &reference) {
	return mpq_class(1);
}

mpq_class util<mpq_class>::from_int(int n, const mpq_class &reference) {
	return mpq_class(n);
}

mpq_class util<mpq_class>::get_gcd(mpq_class a, mpq_class b) {
	return mpq_class(1);
}

mpq_class util<mpq_class>::get_abs(mpq_class a) {
	return abs(a);
}

mpq_class util<mpq_class>::get_pow(mpq_class a, int exp) {
	mpq_class result(1);
	for (int i = 0; i < exp; i++)
		result *= a;
	for (int i = 0; i < -exp; i++)
		result /= a;
	return result;
}