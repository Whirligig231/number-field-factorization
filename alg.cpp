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
		
		// TODO: make this actually do step 4
		
		int e_current_size = e.size();
		for (int i = 0; i < e_current_size; i++) {

			if (e[i].degree() <= 1)
				continue;
			
			// We can't compute t^((p-1)/2) directly due to space constraints.
			// Note that (X, Y) = (X, Y % X), so we use power_mod.
			ZN_X d = t;
			d = e[i].power_mod(d, (p - 1) / 2);
			d -= ZN_X(one<ZN>(a[a.degree()]));
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
	Z a_m = get_abs(a[a.degree()]);
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

std::vector<Z_X> factor(Z_X a) {
	// Algorithm 3.5.7
	
	Z c = a.content();
	a /= c;
	Z_X u = a;
	u /= sub_resultant_gcd(a, a.derivative());
	
	Z p = 1;
	do
		mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
	while (std::get<2>(extended_gcd(u.convert(to_mod(p)), u.derivative().convert(to_mod(p)))).degree() != 0);
	
	std::vector<ZN_X> u_factors = berlekamp_auto(u.convert(to_mod(p)));
	
	Z bound = coeff_bound(u);
	
	return std::vector<Z_X>();
}