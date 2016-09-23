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