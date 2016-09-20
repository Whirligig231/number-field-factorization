#include "alg.h"

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
		ZN_X t = ZN_X(zero<ZN>(a[a.degree()]));
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