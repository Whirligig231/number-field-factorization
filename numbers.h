#pragma once

template <typename T>
T zero(const T &reference) {
	return static_cast<T>(0);
}

template <typename T>
T one(const T &reference) {
	return static_cast<T>(1);
}

template <typename T>
T from_int(int n, const T &reference) {
	return static_cast<T>(n);
}

/*template <>
inline mpz_class get_gcd<mpz_class>(int n, const mod &reference) {
	return mod(reference.get_base(), n);
}*/