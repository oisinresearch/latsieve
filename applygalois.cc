#include <stdint.h>	// int64_t
#include <iostream> // cout
#include <gmpxx.h>
#include "intpoly.h"
#include <math.h>	// sqrt
#include <fstream>	// file
#include <ctime>	// clock_t
#include "mpz_poly.h"
#include <sstream>	// stringstream
#include <stack>	// stack
#include <string>
#include "factorsmall.h"
#include <assert.h>

using std::cout;
using std::endl;
using std::flush;
using std::ifstream;
using std::string;
using std::to_string;
using std::stringstream;
auto npos = std::string::npos;

inline int gcd(int a, int b);
bool known_good_prime(int pt, int* pcache, int pmin, int* sievep, int* sievenum_smodp, int goodnum, int kmin, int k);

int main (int argc, char** argv)
{
	if (argc == 1) {
		cout << endl << "Usage: ./presort factorbase relations" << endl << flush;
		return 0;
	}

	// load factor base
	ifstream fbfile(argv[1]);
	// read fbb
	getline(fbfile, line);
	int fbb = atoi(line.c_str());
	// compute small primes
	int fbmax = fbb;
	char* sieve = new char[fbmax+1]();
	int* primes = new int[1077871];	// hardcoded for the moment.  Allows fbb of 2^23 easily
	for (int i = 2; i <= sqrt(fbmax); i++)
		if(!sieve[i])
			for (int j = i*i; j <= fbmax; j += i)
				if(!sieve[j]) sieve[j] = 1;
	int nump = 0;
	for (int i = 2; i <= fbmax-1; i++)
		if (!sieve[i])
			primes[nump++] = i;

	// set up constants
	std::clock_t start; double timetaken = 0;
	mpz_t r0; mpz_init(r0);
	int* sieves0 = new int[degf * nump]();
	int* sievep0 = new int[nump]();
	int* sievenum_s0modp = new int[nump]();
	int* sieves1 = new int[degg * nump]();
	int* sievep1 = new int[nump]();
	int* sievenum_s1modp = new int[nump]();
	// read k0
	getline(fbfile, line);
	int k0 = atoi(line.c_str());
	for (int i = 0; i < k0; i++) {
		getline(fbfile, line);
		stringstream ss(line);
		string substr;
		getline(ss, substr, ',');
		sievep0[i] = atoi(substr.c_str());
		int j = 0;
		while( ss.good() ) {
			getline( ss, substr, ',' );
			sieves0[i*degf + j++] = atoi(substr.c_str());
		}
		sievenum_s0modp[i] = j;
	}
	// read k1
	getline(fbfile, line);
	int k1 = atoi(line.c_str());
	for (int i = 0; i < k1; i++) {
		getline(fbfile, line);
		stringstream ss(line);
		string substr;
		getline(ss, substr, ',');
		sievep1[i] = atoi(substr.c_str());
		int j = 0;
		while( ss.good() ) {
			getline( ss, substr, ',' );
			sieves1[i*degf + j++] = atoi(substr.c_str());
		}
		sievenum_s1modp[i] = j;
	}
	timetaken += ( clock() - start ) / (double) CLOCKS_PER_SEC;
	fbfile.close();

	// construct p0 cache and p1 cache
	int p0min = 1000;
	int* p0cache = new int[p0min]();
	int c0 = 0; int k0min = 0; 
	for (int p = 2; p < p0min; p++) {
		int t = c0;
		while (sievep0[t] < p) t++;
		if (sievep0[t] == p) { c0 = t + 1; p0cache[p] = 1; k0min++; }	// factor base prime
	}
	int p1min = 1000;
	int* p1cache = new int[p1min]();
	int c1 = 0; int k1min = 0;
	for (int p = 2; p < p1min; p++) {
		int t = c1;
		while (sievep1[t] < p) t++;
		if (sievep1[t] == p) { c1 = t + 1; p1cache[p] = 1; k1min++; }	// factor base prime
	}

	// read relations file
	ifstream rels(argv[2]);
	string separator1 = ":";
	string separator2 = ",";
	while (rels >> line) {
		string line0 = line;
		bool isrel = true;

		string Astr = line.substr(0, line.find(separator1));
		line.erase(0, Astr.length() + 1);
		int a = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		int b = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		int c = atoi(Astr.substr(0, Astr.find(separator2)).c_str());

		// side 0
		string side0 = line.substr(0, line.find(separator1));
		while (side0.length()) {
			string pstr = side0.substr(0, side0.find(separator2));
			if (pstr.length()) {
				mpz_set_str(p, pstr.c_str(), BASE);
				assert(mpz_divisible_p(N0, p));
				mpz_divexact(N0, N0, p);
				if (mpz_cmp_ui(p, fbmax) < 0) {
					int pt = mpz_get_ui(p);
					if (!known_good_prime(pt, p0cache, p0min, sievep0, k0min, k0)) { // relation is bad
						isrel = false;
						break;
					}
				}
				else { // make sure large prime is good
					int pt = mpz_get_ui(p);
					for (int i = 0; i <= degf; i++) fmodp[i] = mod(fi64[i], pt);
					int nr = polrootsmod(fmodp, degf, froots, pt);
					if (nr == 0) {
						isrel = false;
						break;
					}
				}
				int pos = side0.find(separator2);
				if (pos == npos) pos = pstr.length() - 1;
				side0.erase(0, pos + 1);
			}
		}
		if (!isrel) {
			cout << "#" << line0 << endl << flush;
			continue;
		}
		// side 1
		line.erase(0, line.find(separator1) + 1);
		string side1 = line.substr(0, line.find(separator1));
		while (side1.length()) {
			string pstr = side1.substr(0, side1.find(separator2));
			if (pstr.length()) {
				mpz_set_str(p, pstr.c_str(), BASE);
				assert(mpz_divisible_p(N1, p));
				mpz_divexact(N1, N1, p);
				if (mpz_cmp_ui(p, fbmax) < 0) {
					int pt = mpz_get_ui(p);
					if (!known_good_prime(pt, p1cache, p1min, sievep1, k1min, k1)) { // relation is bad
						isrel = false;
						break;
					}
				}
				else { // make sure large prime is good
					int pt = mpz_get_ui(p);
					for (int i = 0; i <= degg; i++) gmodp[i] = mod(gi64[i], pt);
					int nr = polrootsmod(gmodp, degg, groots, pt);
					if (nr == 0) {
						isrel = false;
						break;
					}
				}
				int pos = side1.find(separator2);
				if (pos == npos) pos = pstr.length() - 1;
				side1.erase(0, pos + 1);
			}
		}
		if (!isrel) {
			cout << "#" << line0 << endl << flush;
			continue;
		}

		// print original relation
		cout << line0 << endl << flush;
		
		// compute Galois conjugates
		for (int g = 1; g < 6; g++) {
			int ag = 4*a -2*b +c;
			int bg = 4*a + b -2*c;
			int cg = a + b + c;
			a = ag; b = bg; c = cg;

			int64_t D = b*(int64_t)b - 4*a*(int64_t)c;
			if (floor(sqrt(D)+0.5)*floor(sqrt(D)+0.5) == D) {
				cout << "a + b*x + c*x^2 not irreducible!  Skipping..." << endl << flush;
				continue;	// a+b*x+c*x^2 is not irreducible over Z
			}
			
			// restore line
			line = line0;
			line.erase(0, line.find(separator1) + 1);

			int content = gcd(a, b); content = gcd(content, c);
			a = a/content; b = b/content; c = c/content;
			
			//cout << "[a, b, c] = [" << a << ", " << b << ", " << c << "]" << endl << flush;
			mpz_poly_setcoeff_si(A, 0, a);
			mpz_poly_setcoeff_si(A, 1, b);
			mpz_poly_setcoeff_si(A, 2, c);
			mpz_poly_resultant(N0, f0, A);
			mpz_poly_resultant(N1, f1, A);
			mpz_abs(N0, N0);
			mpz_abs(N1, N1);

			// side 0
			N0primes.clear();
			string side0 = line.substr(0, line.find(separator1));
			while (side0.length()) {
				string pstr = side0.substr(0, side0.find(separator2));
				if (pstr.length()) {
					mpz_set_str(p, pstr.c_str(), BASE);
					if (mpz_divisible_p(N0, p)) {
						mpz_divexact(N0, N0, p);
						N0primes.push_back(pstr);
					}
					int pos = side0.find(separator2);
					if (pos == npos) pos = pstr.length() - 1;
					side0.erase(0, pos + 1);
				}
			}
			if (mpz_sizeinbase(N0, 2) > 63) { cout << "Error - Galois conjugate not smooth!" << endl << flush; return 0; }
			int64_t n = mpz_get_ui(N0);
			int k = factorsmall(n, primes, q, e);
			for (int i = 0; i < k; i++) {
				for (int j = 0; j < e[i]; j++) {
					N0primes.push_back(to_string(q[i]));
				}
			}
			//sort(N0primes.begin(), N0primes.end());
			side0 = ""; int l = N0primes.size();
			for (int i = 0; i < l; i++)
				side0 += N0primes[i] + (i < l - 1 ? "," : "");
			
			// side 1
			N1primes.clear();
			line.erase(0, line.find(separator1) + 1);
			string side1 = line.substr(0, line.find(separator1));
			while (side1.length()) {
				string pstr = side1.substr(0, side1.find(separator2));
				if (pstr.length()) {
					mpz_set_str(p, pstr.c_str(), BASE);
					if (mpz_divisible_p(N1, p)) {
						mpz_divexact(N1, N1, p);
						N1primes.push_back(pstr);
					}
					int pos = side1.find(separator2);
					if (pos == npos) pos = pstr.length() - 1;
					side1.erase(0, pos + 1);
				}
			}
			if (mpz_sizeinbase(N1, 2) > 63) { cout << "Error - Galois conjugate not smooth!" << endl << flush; return 0; }
			n = mpz_get_ui(N1);
			k = factorsmall(n, primes, q, e);
			for (int i = 0; i < k; i++) {
				for (int j = 0; j < e[i]; j++) {
					N1primes.push_back(to_string(q[i]));
				}
			}
			//sort(N1primes.begin(), N1primes.end());
			side1 = ""; l = N1primes.size();
			for (int i = 0; i < l; i++)
				side1 += N1primes[i] + (i < l - 1 ? "," : "");

			string str = to_string(a) + "," + to_string(b) + "," + to_string(c) + ":" + side0 + ":" + side1;

			cout << str << endl << flush;
		}
	}

	mpz_clear(p);
	mpz_clear(N1); mpz_clear(N0);
    mpz_poly_clear(A); mpz_poly_clear(f1); mpz_poly_clear(f0);
	delete[] primes;
	delete[] sieve;
	delete[] groots;
	delete[] froots;
	delete[] gmodp;
	delete[] fmodp;
	delete[] gi64;
	delete[] fi64;
	for (int i = 0; i < 20; i++) {
		mpz_clear(hpoly[i]);
		mpz_clear(gpoly[i]);
		mpz_clear(fpoly[i]);
	}
	delete[] hpoly;
	delete[] gpoly;
	delete[] fpoly;

	return 0;
}


inline int gcd(int a, int b)
{
	a = abs(a);
	b = abs(b);
	int c;
	while (b != 0) {
		c = b;
		b = a % c;
		a = c;
	}
	return a;
}


bool known_good_prime(int pt, int* pcache, int pmin, int* sievep, int* sievenum_smodp, int goodnum, int kmin, int k)
{
	if (pt < pmin) return pcache[pt];

	// else binary subdivision search
	bool known_good = false;
	int max = k;
	int min = kmin;
	int i = min + (max - min) / 2;
	while (true) {
		int i_bak = i;
		if (i != min) {
			if (pt == sievep[i]) { known_good = true; break; }
			else if (pt < sievep[i]) max = i;
			else min = i;
		}
		else {
			i = i_bak + 1;
			if (i != max) {
				if (pt == sievep[i]) { known_good = true; break; }
				else if (pt < sievep[i]) max = i;
				else min = i;
			}
			else {
				// prime pt is less than factor base max but is not good
				break;
			}
		}
		// binary subdivision
		i = min + (max - min) / 2;
	}
	if (known_good == true) if (sievenum_smodp[i] != goodnum) known_good = false;
	return known_good;
}


