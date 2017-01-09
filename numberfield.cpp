#include "numberfield.h"

numberfield::numberfield() {
	this->poly_levels = 0;
	this->rational_value = std::unique_ptr<mpq_class>(new mpq_class(0));
}

numberfield::numberfield(const mpq_class &value) {
	this->poly_levels = 0;
	this->rational_value = std::unique_ptr<mpq_class>(new mpq_class(value));
}

numberfield::numberfield(const polymod<numberfield> &value) {
	this->poly_levels = value.get_base().leading().poly_levels + 1;
	this->poly_value = std::unique_ptr<polymod<numberfield>>(new polymod<numberfield>(value));
}

numberfield::numberfield(const numberfield &other) {
	this->poly_levels = other.poly_levels;
	if (this->poly_levels > 0) {
		this->poly_value = std::unique_ptr<polymod<numberfield>>(new polymod<numberfield>(*(other.poly_value)));
	}
	else {
		mpq_class *orv = new mpq_class(*other.rational_value);
		this->rational_value = std::unique_ptr<mpq_class>(orv);
	}
}

void numberfield::lift(const numberfield &other) {
	if (other.poly_levels <= this->poly_levels)
		return;
	// std::cout << "my levels: " << this->poly_levels << std::endl;
	// std::cout << "other levels: " << other.poly_levels << std::endl;
	// std::cout << "other base: " << other.poly_value->get_base() << std::endl;
	// std::cout << "we are not returning" << std::endl;
	if (other.poly_levels - this->poly_levels > 1) {
		// std::cout << "we take the inner branch" << std::endl;
		numberfield to = other.poly_value->get_base().leading();
		this->lift(to);
	}
	
	// std::cout << "poly of this: " << poly<numberfield>({*this}) << std::endl;
	// std::cout << "polymod: " << polymod<numberfield>(other.poly_value->get_base(), poly<numberfield>({*this})) << std::endl;
	if (this->poly_levels == 0) {
		// Change from rational to poly
		this->poly_value = std::unique_ptr<polymod<numberfield>>(new polymod<numberfield>(other.poly_value->get_base(), poly<numberfield>({*this})));
		this->rational_value = nullptr;
	}
	else {
		*(this->poly_value) = polymod<numberfield>(other.poly_value->get_base(), poly<numberfield>({*this}));
	}
	// std::cout << "we did it" << std::endl;
	this->poly_levels++;
}

numberfield &numberfield::operator=(const mpq_class &value) {
	this->poly_levels = 0;
	this->rational_value = std::unique_ptr<mpq_class>(new mpq_class(value));
	return *this;
}

numberfield &numberfield::operator=(const polymod<numberfield> &value) {
	this->poly_levels = value.get_base().leading().poly_levels + 1;
	this->poly_value = std::unique_ptr<polymod<numberfield>>(new polymod<numberfield>(value));
	return *this;
}

numberfield &numberfield::operator=(const numberfield &other) {
	this->poly_levels = other.poly_levels;
	if (this->poly_levels > 0) {
		this->poly_value = std::unique_ptr<polymod<numberfield>>(new polymod<numberfield>(*(other.poly_value)));
	}
	else {
		this->rational_value = std::unique_ptr<mpq_class>(new mpq_class(*(other.rational_value)));
	}
	return *this;
}

unsigned int numberfield::get_poly_levels() const {
	return this->poly_levels;
}

mpq_class numberfield::get_rational_value() const {
	return (*(this->rational_value));
}

polymod<numberfield> numberfield::get_poly_value() const {
	return (*(this->poly_value));
}

bool numberfield::operator==(const numberfield &other) const {
	numberfield this2(*this);
	numberfield other2(other);
	this2.lift(other2);
	other2.lift(this2);
	if (this2.poly_levels > 0) {
		return *(this2.poly_value) == *(other2.poly_value);
	}
	else {
		return *(this2.rational_value) == *(other2.rational_value);
	}
}

bool numberfield::operator!=(const numberfield &other) const {
	return !(*this == other);
}

numberfield &numberfield::operator+=(const numberfield &other) {
	numberfield other2(other);
	this->lift(other2);
	other2.lift(*this);
	
	if (this->poly_levels > 0) {
		*(this->poly_value) += *(other2.poly_value);
		return *this;
	}
	
	*(this->rational_value) += *(other2.rational_value);
	return *this;
}

numberfield &numberfield::operator-=(const numberfield &other) {
	numberfield other2(other);
	this->lift(other2);
	other2.lift(*this);
	
	if (this->poly_levels > 0) {
		*(this->poly_value) -= *(other2.poly_value);
		return *this;
	}
	
	*(this->rational_value) -= *(other2.rational_value);
	return *this;
}

numberfield &numberfield::operator*=(const numberfield &other) {
	numberfield other2(other);
	this->lift(other2);
	other2.lift(*this);
	
	if (this->poly_levels > 0) {
		*(this->poly_value) *= *(other2.poly_value);
		return *this;
	}
	
	*(this->rational_value) *= *(other2.rational_value);
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
	if (this->poly_levels > 0) {
		return numberfield(-(*(this->poly_value)));
	}
	
	return numberfield(-(*(this->rational_value)));
}

numberfield numberfield::inv() const {
	if (this->poly_levels > 0) {
		return numberfield((*(this->poly_value)).inv());
	}
	
	return numberfield(1 / (*(this->rational_value)));
}

std::ostream &operator<<(std::ostream &os, const numberfield &m) {
	if (m.poly_levels > 0)
		return os << *m.poly_value;
	return os << *m.rational_value;
}

numberfield util<numberfield>::zero() {
	return numberfield();
}

numberfield util<numberfield>::zero(const numberfield &reference) {
	// std::cout << "reference is: " << reference << std::endl;
	// std::cout << "reference levels: " << reference.get_poly_levels() << std::endl;
	// if (reference.get_poly_levels() > 0)
		// std::cout << "reference base: " << reference.get_poly_value().get_base() << std::endl;
	numberfield ret = numberfield() * reference;
	// if (reference.get_poly_levels() > 0)
		// std::cout << "done specifically for reference base: " << reference.get_poly_value().get_base() << std::endl;
	return ret;
}

numberfield util<numberfield>::one(const numberfield &reference) {
	if (reference.get_poly_levels() > 0) {
		return numberfield(util<polymod<numberfield>>::one(reference.get_poly_value()));
	}
	return numberfield(util<mpq_class>::one(reference.get_rational_value()));
}

numberfield util<numberfield>::from_int(int n, const numberfield &reference) {
	if (reference.get_poly_levels() > 0) {
		return numberfield(util<polymod<numberfield>>::from_int(n, reference.get_poly_value()));
	}
	return numberfield(util<mpq_class>::from_int(n, reference.get_rational_value()));
}

numberfield util<numberfield>::get_gcd(numberfield a, numberfield b) {
	return a;
}