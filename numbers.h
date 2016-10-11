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

template <typename T>
T get_gcd(T a, T b) {
	return gcd(a, b);
}

template <typename T>
T get_abs(T a) {
	return abs(a);
}

template <typename T>
T get_sqrt(T a) {
	return sqrt(a);
}

template <typename T>
T get_pow(T a, int exp) {
	T result = one<T>(a);
	for (int i = 0; i < exp; i++)
		result *= a;
	for (int i = 0; i < -exp; i++)
		result /= a;
	return result;
}