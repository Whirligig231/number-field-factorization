#include "alg.h"

Z choose(int n, int r) {
	if (r <= 0 || r >= n)
		return 1;
	
	return choose(n - 1, r) + choose(n - 1, r - 1);
}

int log_bound(Z base, Z pow) {
	// This uses exponentiation by squaring and a binary search
	// to find the smallest exponent e such that base^e > pow.
	// This returns an int because in practice the exponent will
	// never exceed INT_MAX; the value of pow in that case would
	// occupy nearly the entire 32-bit address space.
	
	std::vector<Z> twos; // twos[i] = base^(2^i)
	twos.push_back(base);

	int high = 1;
	
	// Populate the twos array with as much as is needed
	while (twos[twos.size()-1] <= pow) {
		twos.push_back(twos[twos.size()-1]*twos[twos.size()-1]);
		high *= 2;
	}
	
	int low = high/2;
	
	// Binary search
	while (low < high - 1) {
		int mid = (low + high) / 2;
		
		// Compute base^mid
		Z base_mid = 1;
		for (int i = 0; i < twos.size(); i++) {
			if (mid & (1 << i))
				base_mid *= twos[i];
		}
		
		if (base_mid > pow)
			high = mid;
		else
			low = mid;
	}
	
	return high;
}

std::vector<ZN_X> berlekamp_small_p(ZN_X a) {
	// Algorithm 3.4.10
	
	Z p = a[a.degree()].get_base();
	
	// The first step here is actually very tricky.
	// We could simply take the polynomial x^(pk) % a directly,
	// but pk might be too big to store in an int, and 3.1.1 has
	// Omega(pk/n) performance in the worst case.
	
	// So instead we will compute the polynomial power x^p
	// using exponentiation by squaring, keeping the polynomial
	// mod a at every step.
	std::vector<ZN_X> q_polys;
	q_polys.push_back(ZN_X(ZN(p, 1)));
	
	ZN_X xp = a.power_mod(p);
	
	for (int k = 1; k < a.degree(); k++) {
		ZN_X xpk = q_polys[q_polys.size() - 1];
		xpk *= xp;
		xpk %= a;
		q_polys.push_back(xpk);
	}
	
	mat<ZN> q(a.degree(), a.degree());
	for (int k = 0; k < a.degree(); k++)
		for (int i = 0; i < a.degree(); i++)
				q(i, k) = (q_polys[k])[i];

	mat<ZN> q2 = q - mat<ZN>::Identity(a.degree(), a.degree());
	
	std::vector<vec<ZN>> v = kernel(q2);
	
	std::vector<ZN_X> e;
	e.push_back(a);
	
	int k = 1;
	int j = 0;

	while (k < v.size()) {
		
		j++;
		ZN_X t(std::vector<ZN>(v[j].data(), v[j].data() + v[j].size()));
	
		int e_current_size = e.size();
		for (int i = 0; i < e_current_size; i++) {
	
			if (e[i].degree() <= 1)
				continue;
			
			std::vector<ZN_X> f;
			bool first_time = true;
			for (ZN s = ZN(p, 0); (s != ZN(p, 0)) || first_time; s += ZN(p, 1)) {
				first_time = false;
				ZN_X gcd = std::get<2>(extended_gcd(e[i], t-s));
				if (gcd.degree() < 1)
					continue;
				bool put_in = true;
				for (int i2 = 0; i2 < f.size(); i2++)
					if (f[i2] == gcd)
						put_in = false;
				if (put_in)
					f.push_back(gcd);
			}

			if (f.size() > 1) {
				e.erase(e.begin() + i);
				i--;
				k--;
				
				e.insert(e.end(), f.begin(), f.end());
				k += f.size();
			}
			
			if (k == v.size())
				break;

		}
	}
	
	// Berlekamp only gives us the factorization up to a constant factor,
	// because GCD is only defined up to units.
	// To correct for this, we use the leading coefficients.
	
	ZN leading_product = ZN(p, 1);
	for (int i = 0; i < e.size(); i++)
		leading_product *= e[i][e[i].degree()];
	ZN deviation = a[a.degree()] / leading_product;
	e[0] *= deviation;
	
	return e;
}

std::vector<ZN_X> berlekamp(ZN_X a) {
	// Algorithm 3.4.11
	
	gmp_randstate_t state;
	gmp_randinit_default(state);
	
	Z p = a[a.degree()].get_base();

	std::vector<ZN_X> q_polys;
	q_polys.push_back(ZN_X(ZN(p, 1)));
	
	ZN_X xp = a.power_mod(p);
	
	for (int k = 1; k < a.degree(); k++) {
		ZN_X xpk = q_polys[q_polys.size() - 1];
		xpk *= xp;
		xpk %= a;
		q_polys.push_back(xpk);
	}
	
	mat<ZN> q(a.degree(), a.degree());
	for (int k = 0; k < a.degree(); k++)
		for (int i = 0; i < a.degree(); i++)
				q(i, k) = (q_polys[k])[i];

	mat<ZN> q2 = q - mat<ZN>::Identity(a.degree(), a.degree());
	
	std::vector<vec<ZN>> v = kernel(q2);
	
	std::vector<ZN_X> e;
	e.push_back(a);
	
	int k = 1;

	while (k < v.size()) {

		// Compute t
		ZN_X t = ZN_X();
		for (int i = 0; i < v.size(); i++) {
			Z ai_z;
			mpz_urandomm(ai_z.get_mpz_t(), state, p.get_mpz_t());
			ZN ai(p, ai_z);
			ZN_X ti(std::vector<ZN>(v[i].data(), v[i].data() + v[i].size()));
			t += ai*ti;
		}
		
		int e_current_size = e.size();
		for (int i = 0; i < e_current_size; i++) {

			if (e[i].degree() <= 1)
				continue;
			
			// We can't compute t^((p-1)/2) directly due to space constraints.
			// Note that (X, Y) = (X, Y % X), so we use power_mod.
			ZN_X d = t;
			d = e[i].power_mod(d, (p - 1) / 2);
			d -= ZN_X(util<ZN>::one(a[a.degree()]));
			d = std::get<2>(extended_gcd(e[i], d));
			
			if (d.degree() < 1)
				continue;
			
			if (d.degree() >= e[i].degree())
				continue;
			
			e.push_back(d);
			e.push_back(e[i] / d);
			e.erase(e.begin() + i);
			i--;
			k++;
			
			if (k == v.size())
				break;
		}
	}
	
	// Berlekamp only gives us the factorization up to a constant factor,
	// because GCD is only defined up to units.
	// To correct for this, we use the leading coefficients.
	
	ZN leading_product = ZN(p, 1);
	for (int i = 0; i < e.size(); i++)
		leading_product *= e[i][e[i].degree()];
	ZN deviation = a[a.degree()] / leading_product;
	e[0] *= deviation;
	
	return e;
}

std::vector<ZN_X> berlekamp_auto(ZN_X a) {
	Z p = a[a.degree()].get_base();
	if (p < 3)
		return berlekamp_small_p(a);
	return berlekamp(a);
}

Z coeff_bound(Z_X a) {
	// Based on Theorem 3.5.1
	// Given a polynomial a, this returns an upper bound on the absolute values
	// of the coefficients of factors of a with degree at most deg(a)/2.
	
	// Note that we only care about the value of n = deg(b) that maximizes this.
	// From the formula in Theorem 3.5.1, it is clear that this is the highest
	// value we can give n, i.e. n = floor(deg(a)/2).
	
	// Note that we only care about the value of j that gives the highest result.
	// If n is even, this value is n/2, which maximizes the values of the two
	// binomial coefficients in the formula. If n is odd, this value is either
	// (n+1)/2 or (n-1)/2. In the former case, the second coefficient will be
	// larger, while in the latter case, the first coefficient will be larger.
	// The values of the two coefficients switch in these cases; call them c_1
	// and c_2, with c_1 <= c_2.
	
	// Now we know that |a| >= |a_m|, and all of these quantities are going to
	// be nonnegative integers. If |a_m| <= |a| and c_1 <= c_2, then
	// |a_m| c_2 + |a| c_1 <= |a_m| c_1 + |a| c_2.
	// So we want to multiply |a| by c_2, the larger coefficient; thus the
	// first coefficient should be larger, and we use j = (n-1)/2.
	
	// Summing up this information, no pun intended, we will set j = floor(n/2).
	
	int n = a.degree()/2;
	int j = n/2;
	
	// We add one because the square root may be rounded down; we want an upper bound
	Z a_norm = a.norm() + 1;
	Z a_m = util<Z>::get_abs(a[a.degree()]);
	return choose(n-1, j)*a_norm + choose(n-1, j-1)*a_m;
}

std::pair<Z_X, Z_X> hensel_lift(Z p, Z q, Z_X a, Z_X b, Z_X c, Z_X u, Z_X v) {
	// Algorithm 3.5.5

	Z r = gcd(p, q);
	ZN_X f = ((c - a*b)/q).convert(to_mod(r));

	// The book isn't clear on what to do here except to hint that we need to
	// do polynomial division in order to find t such that:
	// (v*f - a*t).degree() < a.degree()
	// Note that if we divide v*f by a, we get t satisfying this criterion.
	ZN_X t = (v.convert(to_mod(r))*f) / a.convert(to_mod(r));

	Z_X a0 = static_cast<Z_X>((v.convert(to_mod(r))*f) - (a.convert(to_mod(r))*t));
	Z_X b0 = static_cast<Z_X>((u.convert(to_mod(r))*f) + (b.convert(to_mod(r))*t));

	Z_X a1 = a + q*a0;
	Z_X b1 = b + q*b0;

	return std::make_pair(a1, b1);
}

std::pair<Z_X, Z_X> quad_hensel_lift(Z p, Z q, Z_X a1, Z_X b1, Z_X u, Z_X v) {
	// Algorithm 3.5.6
	
	Z r = p;
	ZN_X g = ((Z_X(1) - (u*a1) - (v*b1))/p).convert(to_mod(r));
	
	// As in the above function, here we divide v*g by a1.
	ZN_X t = (v.convert(to_mod(r))*g) / a1.convert(to_mod(r));
	
	Z_X u0 = static_cast<Z_X>(u.convert(to_mod(r))*g + b1.convert(to_mod(r))*t);
	Z_X v0 = static_cast<Z_X>(v.convert(to_mod(r))*g - a1.convert(to_mod(r))*t);
	
	Z_X u1 = u + p*u0;
	Z_X v1 = v + p*v0;
	
	return std::make_pair(u1, v1);
}

std::pair<Z_X, Z_X> multi_hensel_lift(Z p, int exp, Z_X a, Z_X b, Z_X c) {
	// Let c = ab (mod p); then this lifts this factorization to a factorization
	// c = a'b' (mod p^exp).
	
	Z q = p;
	Z_X u, v;
	
	// Compute u, v using Euclid
	std::tuple<ZN_X, ZN_X, ZN_X> uvr = extended_gcd(a.convert(to_mod(p)), b.convert(to_mod(p)));
	
	ZN_X un, vn;
	un = std::get<0>(uvr);
	vn = std::get<1>(uvr);
	
	// Note that this will give us r = a constant, but not necessarily 1.
	// So we adjust for this.
	
	ZN r = (un*a.convert(to_mod(p)) + vn*b.convert(to_mod(p)))[0];
	un /= r;
	vn /= r;
	
	u = static_cast<Z_X>(un);
	v = static_cast<Z_X>(vn);

	int exp_current = 1;
	while (exp_current < exp) {
		// Run 3.5.5 and 3.5.6, squaring p and q
		std::pair<Z_X, Z_X> hensel1 = hensel_lift(p, q, a, b, c, u, v);
		Z_X a1 = hensel1.first;
		Z_X b1 = hensel1.second;
		std::pair<Z_X, Z_X> hensel2 = quad_hensel_lift(p, q, a1, b1, u, v);
		Z_X u1 = hensel2.first;
		Z_X v1 = hensel2.second;
		
		// Note that p = q = r.
		// So now c = a1*b1 mod q^2 and u1*a1 + v1*b1 = 1 mod p^2.
		// We can thus replace a, b, u, v, p, q by a1, b1, u1, v1, p^2, q^2.
		
		a = a1;
		b = b1;
		u = u1;
		v = v1;
		p *= p;
		q *= q;
		exp_current <<= 1;
	}
	
	// We may have lifted too far (i.e. to p^exp2 where exp2 > exp).
	// So to reduce explosion issues, we take the coefficients mod p^exp.
	Z pexp = 1;
	for (int i = 0; i < exp; i++)
		pexp *= p;
	a = static_cast<Z_X>(a.convert(to_mod(pexp)));
	b = static_cast<Z_X>(b.convert(to_mod(pexp)));
	
	return std::make_pair(a, b);
}

std::vector<Z_X> poly_hensel_lift(Z p, int exp, std::vector<Z_X> ai, Z_X c) {
	// This is a generalization of the above algorithm to more than two factors.
	// We do this inductively.
	
	std::vector<Z_X> ai_new;
	
	// The inductive process is as follows: first, compute f, the product of all
	// but the first ai. Then c = f*ai[0] mod p.
	// Lift this factorization to a factorization f'*ai_new[0] mod p^exp.
	// Note that f' = f mod p, so f' satisfies the same condition as f on the
	// remaining factors. We repeat the process.
	
	std::vector<Z_X> fi;
	Z_X current_tail(1);
	for (int i = 1; i < ai.size(); i++) {
		current_tail *= ai[ai.size() - i];
		fi.insert(fi.begin(), current_tail);
	}
	
	Z_X sub_product = c;
	for (int i = 0; i < ai.size() - 1; i++) {
		// Now sub_product = ai[i]*fi[i] mod p.
		// Apply Hensel lifting to lift this factorization to mod p^exp.
		
		std::pair<Z_X, Z_X> hensel = multi_hensel_lift(p, exp, ai[i], fi[i], sub_product);
		
		// Append the new ai[i] to the factor list ...
		ai_new.push_back(hensel.first);
		
		// ... and use the new fi[i] as the next sub_product.
		// Note that it's still congruent to the old fi[i] mod p.
		sub_product = hensel.second;
	}
	
	// Push the final sub_product, as it will only have one factor in it.
	ai_new.push_back(sub_product);
	return ai_new;
}

std::vector<Z_X> factor(Z_X a) {
	// Algorithm 3.5.7
	
	if (a.degree() < 0)
		return std::vector<Z_X>({a});

	Z c = a.content();
	a /= c;
	Z_X u = a;
	u = u.ring_exact_divide(sub_resultant_gcd(a, a.derivative())).quotient;
	if (u[u.degree()] < 0)
		u = -u;
	
	// Cast out any factors of x
	int factors_of_x = 0;
	if (u[0] == 0) {
		u >>= 1;
		Z_X a_temp = a;
		while (a_temp[0] == 0) {
			a_temp >>= 1;
			factors_of_x++;
		}
	}
	
	// Check if u is constant
	if (u.degree() < 1) {
		std::vector<Z_X> result;
		
		result.push_back(Z_X(c));
	
		for (int i = 0; i < factors_of_x; i++)
			result.insert(result.begin(), Z_X({0, 1}));
		
		return result;
	}
	
	// If |u_0| < |u_n|, reverse U and note this down for later
	bool is_reversed = false;
	if (util<Z>::get_abs(u[0]) < util<Z>::get_abs(u[u.degree()])) {
		is_reversed = true;
		u = u.reverse();
	}
	
	// std::cout << "u = " << u << std::endl;

	Z p = 1;
	do
		mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
	while (u[u.degree()] % p == 0 || std::get<2>(extended_gcd(u.convert(to_mod(p)), u.derivative().convert(to_mod(p)))).degree() != 0);

	std::vector<ZN_X> u_factors = berlekamp_auto(u.convert(to_mod(p)));
	
	Z bound = coeff_bound(u);
	int exp = log_bound(p, 2*u[u.degree()]*bound);
	
	Z pexp = 1;
	for (int i = 0; i < exp; i++)
		pexp *= p;
	
	// std::cout << "p^e = " << p << "^" << exp << " = " << pexp << std::endl;
	
	std::vector<Z_X> ui;
	for (int i = 0; i < u_factors.size(); i++)
		ui.push_back(static_cast<Z_X>(u_factors[i]));
	ui = poly_hensel_lift(p, exp, ui, u);
	
	// Convert to monic (poly_hensel_lift doesn't do this)
	for (int i = 0; i < ui.size(); i++) {
		ZN_X ui_n = ui[i].convert(to_mod(pexp));
		ui_n /= ui_n[ui_n.degree()];
		ui[i] = static_cast<Z_X>(ui_n);
		// std::cout << "u" << i << " = " << ui[i] << std::endl;
	}
	
	std::vector<Z_X> result;
	
	int d = 1;
	while (2*d <= ui.size()) {
		
		// Initialize a combination of length d
		std::vector<int> combination;
		for (int i = 0; i < d; i++)
			combination.push_back(i);
		
		bool terminate_early = false;
		bool reset_this_d = false;
		
		while (true) {
			// We want to include ui[0] if d = 1/2 r
			if (2*d == ui.size() && combination[0] > 0)
				break;
			
			Z_X v_bar(1);
			for (int i = 0; i < d; i++)
				v_bar *= ui[combination[i]];
			
			// std::cout << "v_bar = " << v_bar << std::endl;

			Z_X v;
			
			if (v_bar.degree()*2 <= u.degree()) {
				v = static_cast<Z_X>((v_bar * u[u.degree()]).convert(to_mod(pexp)));
			}
			else {
				v = static_cast<Z_X>(u.convert(to_mod(pexp)) / v_bar.convert(to_mod(pexp)));
			}
			
			// The coefficients will be in [0, p^e - 1];
			// let's fix this!
			for (int i = 0; i <= v.degree(); i++)
				if (v[i]*2 >= pexp)
					v.set(i, v[i] - pexp);
				
			// std::cout << "v = " << v << std::endl;
			
			// Cohen recommends checking for divisibility of the constant terms first.
			if (u[u.degree()]*u[0] % v[0] == 0) {

				qr_pair<Z_X> test_qr = (u*u[u.degree()]).pseudo_divide(v);
				Z modbase = 1;
				for (int i = 0; i < u.degree() - v.degree() + 1; i++)
					modbase *= v[v.degree()];
				if (test_qr.remainder.degree() < 0 && test_qr.quotient.content() % modbase == 0) {
					// We did it! We found a factor!
					Z_X f = v / v.content();
					
					if (is_reversed)
						f = f.reverse();

					Z_X a_temp = a;
					while (true) {
						qr_pair<Z_X> qr = a_temp.pseudo_divide(f);
						modbase = 1;
						for (int i = 0; i < a_temp.degree() - f.degree() + 1; i++)
							modbase *= f[f.degree()];
						if (qr.remainder.degree() >= 0)
							break;
						if (qr.quotient.content() % modbase != 0)
							break;
						a_temp = a_temp.ring_exact_divide(f).quotient;
						result.push_back(f);
					}
					
					if (is_reversed)
						f = f.reverse();

					u = u.ring_exact_divide(f).quotient;
					if (2*d <= ui.size()) {
						for (int i = 0; i < combination.size(); i++)
							ui.erase(ui.begin() + combination[i] - i);
					}
					else {
						std::vector<Z_X> new_ui;
						for (int i = 0; i < combination.size(); i++)
							new_ui.push_back(ui[combination[i]]);
						ui = new_ui;
					}

					if (2*d > ui.size())
						terminate_early = true;
					
					d--;
					break;
				}
				
			}

			// Increment combination
			int start_point = d - 1;
			combination[d - 1]++;
			while (combination[d - 1] >= ui.size()) {
				start_point--;
				if (start_point < 0)
					break;
				combination[start_point]++;
				for (int i = start_point + 1; i < ui.size(); i++)
					combination[i] = combination[i-1] + 1;
			}
			if (start_point < 0)
				break;
		}
		
		if (terminate_early)
			break;
		
		d++;
		
	}
	
	Z_X f = u / u.content();
	Z_X a_temp = a;
	if (is_reversed)
		f = f.reverse();
	while (true) {
		qr_pair<Z_X> qr = a_temp.pseudo_divide(f);
		Z modbase = 1;
		for (int i = 0; i < a_temp.degree() - f.degree() + 1; i++)
			modbase *= f[f.degree()];
		if (qr.remainder.degree() >= 0)
			break;
		if (qr.quotient.content() % modbase != 0)
			break;
		a_temp = a_temp.ring_exact_divide(f).quotient;
		result.push_back(f);
	}
	result.push_back(Z_X(c));
	
	for (int i = 0; i < factors_of_x; i++)
		result.insert(result.begin(), Z_X({0, 1}));
	
	return result;
}

std::vector<Q_X> factor(Q_X a) {
	// Apply Gauss's lemma
	
	Z common_denominator = 1;
	for (int i = 0; i <= a.degree(); i++) {
		if (a[i] != 0)
			common_denominator = lcm(common_denominator, a[i].get_den());
	}
	
	Z_X a2 = static_cast<Z_X>(a * static_cast<Q>(common_denominator));
	std::vector<Z_X> result2 = factor(a2);
	std::vector<Q_X> result;
	for (int i = 0; i < result2.size(); i++)
		result.push_back(static_cast<Q_X>(result2[i]));
	result[result.size() - 1] /= common_denominator;
	return result;
}

std::vector<C> find_complex_roots(Q_X p_q, int precision) {
	// Algorithm 3.6.6
	
	R prec_limit = 1;
	for (int i = 0; i < precision; i++)
		prec_limit /= 10;
	
	// Note that unlike 3.6.6, here we do not assume p to be squarefree
	p_q = p_q / sub_resultant_gcd(p_q, p_q.derivative());
	std::vector<C> p_coeffs;
	for (int i = 0; i <= p_q.degree(); i++)
		p_coeffs.push_back((C)p_q[i]);
	
	C_X p(p_coeffs);
	C_X q = p;
	C_X p_prime = p.derivative();
	C_X q_prime = p_prime;
	int n = p.degree();
	
	std::vector<C> roots;
	
	while (n > 0) {
		C x(1.3, 0.314159);
		C v = q.evaluate(x);
		R m = v.norm();
		
		while (true) {
			int c = 0;
			C dx = v / q_prime.evaluate(x);
			if (sqrt(dx.norm()) < prec_limit)
				break;
			
			while (true) {
				C y = x - dx;
				C v1 = q.evaluate(y);
				R m1 = v1.norm();
				
				if (m1 < m) {
					x = y;
					v = v1;
					m = m1;
					break;
				}
				
				c++;
				dx /= 4;
				if (c >= 20) {
					// Failure
					std::cout << "Failure: this polynomial is nasty!" << std::endl;
					return roots;
				}
			}
		}
		
		x -= p.evaluate(x)/p_prime.evaluate(x);
		x -= p.evaluate(x)/p_prime.evaluate(x);
		
		if (x.get_imag() < prec_limit) {
			x = x.get_real();
			roots.push_back(x);
			q /= C_X({-x, 1});
			q_prime = q.derivative();
			n--;
		}
		else {
			roots.push_back(x);
			roots.push_back(x.conjugate());
			q /= C_X({x.norm(), R(-2*x.get_real()), 1});
			q_prime = q.derivative();
			n -= 2;
		}
	}
	
	return roots;
}