#include <cstdlib>	// atoi
#include <vector>	// vector
#include <iostream>	// cout
#include <cmath>	// sqrt
#include "factorsmall.h"
#include "intpoly.h"

#define MAXP 3512
#define MAXS 32767

using std::vector;
using std::cout;
using std::endl;

#ifdef PRIMITIVEMAIN
int main(int argc, char** argv)
{
	if (argc == 2) {
		int p = atoi(argv[1]);
		
		// sieve of Eratosthenes
		int* primes = new int[MAXP]; 	// MAXP primes up to 32767
		int nump = Eratosthenes(primes);

		int r = primitiveroot(p, primes);

		cout << r << " is a primitive root mod " << p << endl;
	}
	return 0;
}
#endif /* PRIMITIVEMAIN */


#ifdef POLROOTSMAIN
int main(int argc, char** argv)
{
	if (argc == 2) {
		int f[13] = { 923478, 92387, 23476, 23476, 923746, 23473, 239874, 23874, 182734, 293874, 293847, 349857, 430958 };

		int p = atoi(argv[1]);
		// reduce f mod p
		for (int i = 0; i <= 12; i++) f[i] = mod(f[i], p);

		/*
		// sieve of Eratosthenes
		int* primes = new int[MAXP]; 	// MAXP primes up to 32767
		int nump = Eratosthenes(primes);
		*/

		int r[12];
		int degf = 12; while (f[degf] == 0) degf--;
		//int k = chienrootsmod(f, 12, p, r, primes);
		int k = polrootsmod(f, degf, r, p);
		for (int i = 0; i < k; i++) cout << r[i] << endl;
		cout << k << " roots mod " << p << endl;
	}
	return 0;
}
#endif /* POLROOTSMAIN */


// https://math.stackexchange.com/questions/124408/finding-a-primitive-root-of-a-prime-number
int primitiveroot(int p, int* q)
{
	if (p == 2) return 1;

	vector<int> phi;
	vector<int> e;
	
	int s = p-1;
	int k = factorsmall(s, q, phi, e);

	int r = 2;
	bool found = false;
	while (!found) {
		found = true;
		for (int i = 0; i < k; i++) {
			int pow = s/phi[i];
			if (modpow(r, pow, p) == 1) {
				found = false;
				break;
			}
		}
		if (found) break;
		r++;
	}
	return r;
}


int mod(int64_t r, int64_t p)
{
	return ( (r % p) + p ) % p;
}


int modpow(int64_t r, int64_t e, int64_t p)
{
	int64_t t = 1;
	while (e) {
		if (e & 1) t = (r * t) % p;
		r = (r * r) % p;
		e >>= 1;
	}
	return t;
}

// use Chien search
int chienrootsmod(int* f, int degf, int p, int* r, int* q)
{
	int g = primitiveroot(p, q);
	int64_t* gpow = new int64_t[degf+1];
	for (int j = 0; j <= degf; j++)
		gpow[j] = modpow(g, j, p);

	int64_t* t = new int64_t[degf+1];
	for (int j = 0; j <= degf; j++)
		t[j] = f[j] % p;

	int k = 0;
	if (f[0] % p == 0) { r[k++] = 0; }	// trivial case
	for (int i = 1; i < p; i++) {
		int sum = 0;
		for (int j = 0; j <= degf; j++) {
			t[j] = (t[j] * gpow[j]) % p;
			sum += t[j];
		}
		if (sum % p == 0) r[k++] = modpow(g, i, p);
		if (k == degf) break;	// there are at most degf roots
	}

	// clean up
	delete[] t;
	delete[] gpow;

	return k;
}


inline int modinv(int x, int m)
{
	int m0 = m, t, q;
	int y0 = 0, y = 1;
	if (m == 1) return 1;
	while (x > 1) {
		q = x / m;
		t = m, m = x % m, x = t;
		t = y0, y0 = y - q * y0, y = t;
	}
	if (y < 0) y += m0;
	return y;
}


inline int poldegree(int64_t* f, int maxd)
{
	int d; 
	for (d = maxd; d >= 0; d--) if (f[d] != 0) break;
	return d;
}


inline void polsetzero(int64_t* f, int maxd)
{
	for (int i = 0; i <= maxd; i++) f[i] = 0;
}


inline bool poliszero(int64_t* f, int maxd)
{
	return poldegree(f, maxd) == -1;
}


int polrootsmod(int* f, int degf, int* roots, int p)
{
	// make f monic
	int64_t fd = f[degf];
	int64_t fdinv = modinv(fd, p);
	for (int i = 0; i <= degf; i++) f[i] = mod(f[i]*fdinv, p);

	// compute r = x^p - x (mod f, p)
	int64_t* r = new int64_t[degf*2+2]();
	int64_t* r2 = new int64_t[degf*2+2]();
	r[0] = 1;
	int64_t mask = 1;
	while (mask << 1 <= p) mask <<= 1;
	while (mask) {
		// clear r2
		polsetzero(r2, degf*2);
		// square
		for (int j = 0; j <= degf; j++) {
			for (int i = 0; i <= degf; i++) {
				r2[i + j] = mod(r2[i + j] + r[i] * r[j], p);
			}
		}
		if (p & mask) {
			// multiply (in this case by x)
			int degr2 = poldegree(r2, degf*2);
			for (int i = degr2+1; i > 0; i--)
				r2[i] = r2[i-1];
			r2[0] = 0;
		}
		// now reduce r mod f using 'subtraction-based' reduction
		int degr2 = poldegree(r2, degf*2);
		while (degr2 >= degf) {
			int64_t r2d = r2[degr2];
			// subtract r2d * f * x^(degr2-degf) from r2.  Note f has been made monic
			int off = degr2 - degf;
			for (int i = 0; i <= degf; i++)
				r2[off + i] = mod(r2[off + i] - r2d * f[i], p);
			degr2 = poldegree(r2, degf*2);	// get actual degree of r2
		}
		// copy reduced r2 back to r
		polsetzero(r, degf*2+1);
		for (int i = 0; i <= degr2; i++) r[i] = r2[i];
		mask = mask >> 1;
	}
	// subtract x yielding r = x^p-x (mod f, p)
	r[1] = mod(r[1] - 1, p);

	// compute g = gcd(r, f) (mod p)
	int64_t* f1 = new int64_t[degf+1];
	int64_t* f2 = new int64_t[degf+1];
	for (int i = 0; i <= degf; i++) { f1[i] = f[i]; f2[i] = r[i]; }
	int degf1 = poldegree(f1, degf);	// get actual degree of f1
	int degf2 = poldegree(f2, degf);	// get actual degree of f2
	while (degf1 > 0) {
		degf2 = poldegree(f2, degf);	// get actual degree of f2
		if (degf2 > degf1) {	// make f1 'largest'
			int64_t* t = f1; f1 = f2; f2 = t;
			int degt = degf1; degf1 = degf2; degf2 = degt;
		}
		// make f2 monic
		int64_t f2d = f2[degf2];
		int64_t f2dinv = modinv(f2d, p);
		for (int i = 0; i <= degf2; i++) f2[i] = mod(f2[i]*f2dinv, p);
		// now reduce f1 mod f2 using 'subtraction-based' reduction
		degf1 = poldegree(f1, degf);	// get actual degree of f1
		while (degf1 >= degf2) {
			int64_t f1d = f1[degf1];
			// subtract f1d * f2 * x^(degf1-degf2) from f1.  Note f2 has been made monic
			int off = degf1 - degf2;
			for (int i = 0; i <= degf2; i++)
				f1[off + i] = mod(f1[off + i] - f1d * f2[i], p);
			degf1 = poldegree(f1, degf);	// get actual degree of f1
		}
	}
	// assign g = f1 if f1 != 0, otherwise assign g = f2
	int64_t* g = new int64_t[degf+1]();
	if (!poliszero(f1, degf1)) {
		for (int i = 0; i <= degf1; i++) g[i] = f1[i];
	}
	else {
		for (int i = 0; i <= degf2; i++) g[i] = f2[i];
	}
	// --------------------------------------------------------------- //
	int degg = poldegree(g, degf);
	// ----------- now factor g using Cantor-Zassenhaus -------------- //
	int64_t* gsplit = new int64_t[degg*2 + 2]();
	for (int i = 0; i <= degg; i++) gsplit[i] = g[i];
	int* gdegs = new int[degg]();
	int k = degg + 1;
	int kmax = 2*degg;
	if (p == 2) kmax--;
	int i0 = 0;
	gdegs[i0] = degg;
	int gidx = 0;
	int64_t* h1 = new int64_t[degg+1]();
	int64_t* h2 = new int64_t[degg+1]();
	// make g monic
	int64_t gd = g[degg];
	int64_t gdinv = modinv(gd, p);
	for (int i = 0; i <= degg; i++) g[i] = mod(g[i]*gdinv, p);
	int64_t a = 0;
	int e = (p-1)>>1;
	// in a list of k linear factors of g, there are 2*degg coefficients, or 2*degg-1 if x is a factor of original g
	while (k < kmax) { 
		// compute (x + a)^((p-1)/2) (mod g, p)
		polsetzero(r, degg);
		r[0] = 1;
		mask = 1;
		while (mask << 1 <= e) mask <<= 1;
		while (mask) {
			// clear r2
			polsetzero(r2, degf*2);
			// square
			for (int j = 0; j <= degg; j++) {
				for (int i = 0; i <= degg; i++) {
					r2[i + j] = mod(r2[i + j] + r[i] * r[j], p);
				}
			}
			if (e & mask) {
				// multiply (in this case by x + a)
				int degr2 = poldegree(r2, degg*2);
				for (int i = degr2+1; i > 0; i--)
					r2[i] = r2[i-1];
				r2[0] = 0;
				for (int i = 0; i <= degr2; i++) {
					r2[i] = mod(r2[i] + a * r2[i+1], p);
				}
			}
			// now reduce r2 mod g using 'subtraction-based' reduction
			int degr2 = poldegree(r2, degg*2);
			while (degr2 >= degg) {
				int64_t r2d = r2[degr2];
				// subtract r2d * g * x^(degr2-degg) from r2.  Note g has been made monic
				int off = degr2 - degg;
				for (int i = 0; i <= degg; i++)
					r2[off + i] = mod(r2[off + i] - r2d * g[i], p);
				degr2 = poldegree(r2, degg*2);	// get actual degree of r2
			}
			// copy reduced r2 back to r
			polsetzero(r, degg*2);
			for (int i = 0; i <= degg; i++) r[i] = r2[i];
			mask = mask >> 1;
		}
		// subtract 1
		r[0] = mod(r[0] - 1, p);
		// compute gcd(r, g) (mod p) in the hope that it is nontrivial
		for (int i = 0; i <= degg; i++) { f1[i] = g[i]; f2[i] = r[i]; }
		degf1 = poldegree(f1, degg);	// get actual degree of f1
		degf2 = poldegree(f2, degg);	// get actual degree of f2
		while (degf1 > 0 && degf2 > 0) {
			degf2 = poldegree(f2, degf);	// get actual degree of f2
			if (degf2 > degf1) {	// make f1 'largest'
				int64_t* t = f1; f1 = f2; f2 = t;
				int degt = degf1; degf1 = degf2; degf2 = degt;
			}
			// make f2 monic
			int64_t f2d = f2[degf2];
			int64_t f2dinv = modinv(f2d, p);
			for (int i = 0; i <= degf2; i++) f2[i] = mod(f2[i]*f2dinv, p);
			// now reduce f1 mod f2 using 'subtraction-based' reduction
			degf1 = poldegree(f1, degf);	// get actual degree of f1
			while (degf1 >= degf2) {
				int64_t f1d = f1[degf1];
				// subtract f1d * f2 * x^(degf1-degf2) from f1.  Note f2 has been made monic
				int off = degf1 - degf2;
				for (int i = 0; i <= degf2; i++)
					f1[off + i] = mod(f1[off + i] - f1d * f2[i], p);
				degf1 = poldegree(f1, degf);	// get actual degree of f1
			}
		}
		int degh1 = 0;
		// assign h1 = f1 if f1 != 0, otherwise assign h1 = f2
		if (degf1 > 0) {
			for (int i = 0; i <= degf1; i++) h1[i] = f1[i];
			degh1 = poldegree(h1, degf1);
		}
		else if (degf2 > 0) {
			for (int i = 0; i <= degf2; i++) h1[i] = f2[i];
			degh1 = poldegree(h1, degf2);
		}
		if (degh1 > 0 && degh1 < degg) {	// found a factor!
			// divide g by h1 to get other factor
			polsetzero(h2,degg); polsetzero(r, degf*2+1);
			for (int i = 0; i <= degg; i++) r[i] = g[i];
			int64_t c = modinv(h1[degh1], p);
			int degr = poldegree(r, degg);
			while (degr >= degh1) {
				int64_t s = mod(r[degr] * c, p);
				h2[degr - degh1] = mod(h2[degr - degh1] + s, p);
				for (int i = 0; i <= degh1; i++) {
					r[degr - degh1 + i] = mod(r[degr - degh1 + i] - s * h1[i], p);
				}
				degr = poldegree(r, degg);
			}
			// expand list of known factors
			k++;
			// shift degrees of factors after gdegs[i0] on by one
			for (int i = (k>>1) - 1; i > i0; i--) gdegs[i] = gdegs[i-1];
			// insert new gdegs
			gdegs[i0] = degh1; gdegs[i0+1] = degg - degh1;
			// shift every coeff after position gidx + degg on by one (to make room for new factors)
			for (int i = k; i > gidx + degg; i--) gsplit[i] = gsplit[i-1];
			// replace coeffs of g by those of the two new non-trivial factors
			for (int i = 0; i <= degh1; i++) gsplit[gidx + i] = h1[i];
			// add cofactor
			int degh2 = degg - degh1;
			for (int i = 0; i <= degh2; i++) gsplit[gidx + degh1 + 1 + i] = h2[i];
		}
		if (k < kmax) {
			// select first nonlinear factor
			while (gdegs[i0] == 1) gidx += gdegs[i0++] + 1;
			// store in g
			degg = gdegs[i0];
			for (int i = 0; i <= degg; i++) g[i] = gsplit[gidx + i];
			// advance to next a
			a++; if (a == p) a = 0;
		}
	}
	// enumerate roots
	for (int i = 0; i < k >> 1; i++) {
		roots[i] = mod(-gsplit[2*i] * modinv(gsplit[2*i+1], p), p);
	}

	// clean up
	delete[] h2;
	delete[] h1;
	delete[] gsplit;
	delete[] g;
	delete[] f2;
	delete[] f1;
	delete[] r2;
	delete[] r;

	return k >> 1;
}

