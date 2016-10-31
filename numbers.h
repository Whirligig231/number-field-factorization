#pragma once

template <typename T>
class util {
public:
	static T zero();
	static T zero(const T &reference);
	static T one();
	static T one(const T &reference);
	static T from_int(int n, const T &reference);
	static T get_gcd(T a, T b);
	static T get_abs(T a);
	static T get_sqrt(T a);
	static T get_pow(T a, int exp);
};

template <typename T>
T util<T>::zero() {
	return static_cast<T>(0);
}

template <typename T>
T util<T>::one() {
	return static_cast<T>(1);
}

template <typename T>
T util<T>::zero(const T &reference) {
	return static_cast<T>(0);
}

template <typename T>
T util<T>::one(const T &reference) {
	return static_cast<T>(1);
}

template <typename T>
T util<T>::from_int(int n, const T &reference) {
	return static_cast<T>(n);
}

template <typename T>
T util<T>::get_gcd(T a, T b) {
	return gcd(a, b);
}

template <typename T>
T util<T>::get_abs(T a) {
	return abs(a);
}

template <typename T>
T util<T>::get_sqrt(T a) {
	return sqrt(a);
}

template <typename T>
T util<T>::get_pow(T a, int exp) {
	T result = util<T>::one(a);
	for (int i = 0; i < exp; i++)
		result *= a;
	for (int i = 0; i < -exp; i++)
		result /= a;
	return result;
}