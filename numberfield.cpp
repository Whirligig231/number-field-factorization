#include "numberfield.h"

void numberfield_error() {
	std::cout << "Error: crossed the levels of nesting" << std::endl;
	int k = 1;
	k /= (k - k);
}

numberfield::numberfield() {
	this->is_valid = 0;
}

numberfield::numberfield(const mpq_class &value) {
	this->is_valid = 1;
	this->is_poly_value = 0;
	this->rational_value = std::unique_ptr<mpq_class>(new mpq_class(value));
}

numberfield::numberfield(const polymod<numberfield> &value) {
	this->is_valid = 1;
	this->is_poly_value = 1;
	this->poly_value = std::unique_ptr<polymod<numberfield>>(new polymod<numberfield>(value));
}

numberfield::numberfield(const numberfield &other) {
	this->is_valid = other.is_valid;
	this->is_poly_value = other.is_poly_value;
	if (this->is_poly_value) {
		this->poly_value = std::unique_ptr<polymod<numberfield>>(new polymod<numberfield>(*(other.poly_value)));
	}
	else {
		this->rational_value = std::unique_ptr<mpq_class>(new mpq_class(*(other.rational_value)));
	}
}

numberfield &numberfield::operator=(const mpq_class &value) {
	this->is_valid = 1;
	this->is_poly_value = 0;
	this->rational_value = std::unique_ptr<mpq_class>(new mpq_class(value));
	return *this;
}

numberfield &numberfield::operator=(const polymod<numberfield> &value) {
	this->is_valid = 1;
	this->is_poly_value = 1;
	this->poly_value = std::unique_ptr<polymod<numberfield>>(new polymod<numberfield>(value));
	return *this;
}

numberfield &numberfield::operator=(const numberfield &other) {
	this->is_valid = other.is_valid;
	this->is_poly_value = other.is_poly_value;
	if (this->is_poly_value) {
		this->poly_value = std::unique_ptr<polymod<numberfield>>(new polymod<numberfield>(*(other.poly_value)));
	}
	else {
		this->rational_value = std::unique_ptr<mpq_class>(new mpq_class(*(other.rational_value)));
	}
	return *this;
}

bool numberfield::is_poly() const {
	return this->is_poly_value;
}

mpq_class numberfield::get_rational_value() const {
	return (*(this->rational_value));
}

polymod<numberfield> numberfield::get_poly_value() const {
	return (*(this->poly_value));
}

bool numberfield::operator==(const numberfield &other) const {
	if (!this->is_valid) {
		return (other == *this);
	}
	if (!other.is_valid) {
		if (this->is_poly_value)
			return (*(this->poly_value) == util<polymod<numberfield>>::zero(*(this->poly_value)));
		else
			return (*(this->rational_value) == util<mpq_class>::zero(*(this->rational_value)));
	}
	if (this->is_poly_value && !other.is_poly_value) {
		numberfield_error();
		return false;
	}
	if (!this->is_poly_value && other.is_poly_value) {
		numberfield_error();
		return false;
	}
	if (this->is_poly_value) {
		return (*(this->poly_value) == *(other.poly_value));
	}
	
	return (*(this->rational_value) == *(other.rational_value));
}

bool numberfield::operator!=(const numberfield &other) const {
	return !(*this == other);
}

numberfield &numberfield::operator+=(const numberfield &other) {
	if (!this->is_valid) {
		*this = other;
		return *this;
	}
	if (!other.is_valid) {
		return *this;
	}
	if (this->is_poly_value && !other.is_poly_value) {
		numberfield_error();
		return *this;
	}
	if (!this->is_poly_value && other.is_poly_value) {
		numberfield_error();
		return *this;
	}
	if (this->is_poly_value) {
		*(this->poly_value) += *(other.poly_value);
		return *this;
	}
	
	*(this->rational_value) += *(other.rational_value);
	return *this;
}

numberfield &numberfield::operator-=(const numberfield &other) {
	if (!this->is_valid) {
		*this = (-other);
		return *this;
	}
	if (!other.is_valid) {
		return *this;
	}
	if (this->is_poly_value && !other.is_poly_value) {
		numberfield_error();
		return *this;
	}
	if (!this->is_poly_value && other.is_poly_value) {
		numberfield_error();
		return *this;
	}
	if (this->is_poly_value) {
		*(this->poly_value) -= *(other.poly_value);
		return *this;
	}
	
	*(this->rational_value) -= *(other.rational_value);
	return *this;
}

numberfield &numberfield::operator*=(const numberfield &other) {
	if (!this->is_valid) {
		this->is_poly_value = other.is_poly_value;
		if (this->is_poly_value)
			*(this->poly_value) = util<polymod<numberfield>>::zero(*(other.poly_value));
		else
			*(this->rational_value) = util<mpq_class>::zero(*(other.rational_value));
		return *this;
	}
	if (!other.is_valid) {
		if (this->is_poly_value)
			*(this->poly_value) = util<polymod<numberfield>>::zero(*(this->poly_value));
		else
			*(this->rational_value) = util<mpq_class>::zero(*(this->rational_value));
		return *this;
	}
	if (this->is_poly_value && !other.is_poly_value) {
		numberfield_error();
		return *this;
	}
	if (!this->is_poly_value && other.is_poly_value) {
		numberfield_error();
		return *this;
	}
	if (this->is_poly_value) {
		*(this->poly_value) *= *(other.poly_value);
		return *this;
	}
	
	*(this->rational_value) *= *(other.rational_value);
	return *this;
}

numberfield &numberfield::operator/=(const numberfield &other) {
	return (*this) *= other.inv();
}

numberfield numberfield::operator+(const numberfield &other) const {
	return numberfield(*this) += other;
}

numberfield numberfield::operator-(const numberfield &other) const {
	return numberfield(*this) -= other;
}

numberfield numberfield::operator*(const numberfield &other) const {
	return numberfield(*this) *= other;
}

numberfield numberfield::operator/(const numberfield &other) const {
	return numberfield(*this) /= other;
}

numberfield numberfield::operator-() const {
	if (!this->is_valid)
		return numberfield(*this);
	if (this->is_poly_value) {
		return numberfield(-(*(this->poly_value)));
	}
	
	return numberfield(-(*(this->rational_value)));
}

numberfield numberfield::inv() const {
	if (!this->is_valid) {
		std::cout << "Error: divide by zero in numberfield" << std::endl;
		int a = 1;
		a /= (a - a);
		return numberfield(*this);
	}
	if (this->is_poly_value) {
		return numberfield((*(this->poly_value)).inv());
	}
	
	return numberfield(1 / (*(this->rational_value)));
}

std::ostream &operator<<(std::ostream &os, const numberfield &m) {
	if (!m.is_valid) {
		return os << "0";
	}
	if (m.is_poly_value)
		return os << *m.poly_value;
	return os << *m.rational_value;
}

numberfield util<numberfield>::zero() {
	return numberfield();
}

numberfield util<numberfield>::zero(const numberfield &reference) {
	return numberfield() * reference;
}

numberfield util<numberfield>::one(const numberfield &reference) {
	if (reference.is_poly()) {
		return numberfield(util<polymod<numberfield>>::one(reference.get_poly_value()));
	}
	return numberfield(util<mpq_class>::one(reference.get_rational_value()));
}

numberfield util<numberfield>::from_int(int n, const numberfield &reference) {
	if (reference.is_poly()) {
		return numberfield(util<polymod<numberfield>>::from_int(n, reference.get_poly_value()));
	}
	return numberfield(util<mpq_class>::from_int(n, reference.get_rational_value()));
}

numberfield util<numberfield>::get_gcd(numberfield a, numberfield b) {
	return a;
}