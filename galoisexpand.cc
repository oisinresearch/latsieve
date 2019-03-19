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
bool known_good_prime(int pt, int* pcache, int pmin, int* sievep, int kmin, int k);

int main (int argc, char** argv)
{
	if (argc == 1) {
		cout << endl << "Usage: ./galoisexpand inputpoly factorbase relations" << endl << flush;
		return 0;
	}

	// read input polynomial
	mpz_t* fpoly = new mpz_t[20];	// max degree of 20.  Not the neatest
	mpz_t* gpoly = new mpz_t[20];	// max degree of 20.  Not the neatest
	mpz_t* hpoly = new mpz_t[20];	// max degree of 20.  Not the neatest
	for (int i = 0; i < 20; i++) {
		mpz_init(fpoly[i]);
		mpz_init(gpoly[i]);
		mpz_init(hpoly[i]);
	}
	int64_t* fi64 = new int64_t[20]();
	int64_t* gi64 = new int64_t[20]();
	int64_t* fmodp = new int64_t[20]();
	int64_t* gmodp = new int64_t[20]();
	string line;
	char linebuffer[100];
	ifstream file(argv[1]);
	getline(file, line);	// first line contains number n to factor
	// read nonlinear poly
	int degf = -1;
	while (getline(file, line) && line.substr(0,1) == "c" ) {
		line = line.substr(line.find_first_of(" ")+1);
		mpz_set_str(fpoly[++degf], line.c_str(), 10);
		fi64[degf] = mpz_get_si(fpoly[degf]);
	}
	// read other poly
	int degg = -1;
	bool read = true;
	while (read && line.substr(0,1) == "Y" ) {
		line = line.substr(line.find_first_of(" ")+1);
		mpz_set_str(gpoly[++degg], line.c_str(), 10);
		gi64[degg] = mpz_get_si(gpoly[degg]);
		read = static_cast<bool>(getline(file, line));
	}
	file.close();
	int64_t* froots = new int64_t[degf+1];
	int64_t* groots = new int64_t[degg+1];

	// load factor base
	ifstream fbfile(argv[2]);
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
	mpz_poly f0; mpz_poly f1; mpz_poly A;
	mpz_poly_init(f0, degf); mpz_poly_init(f1, degg); mpz_poly_init(A, 3);
	mpz_poly_set_mpz(f0, fpoly, degf);
	mpz_poly_set_mpz(f1, gpoly, degg);
	mpz_t N0; mpz_init(N0); mpz_t N1; mpz_init(N1);
	mpz_t p; mpz_init(p); int BASE = 16;
	vector<int> q; vector<int> e;
	vector<string> N0primes;
	vector<string> N1primes;
	ifstream rels(argv[3]);
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

		// check if relation is valid, i.e. all primes are factor base primes
		mpz_poly_setcoeff_si(A, 0, a);
		mpz_poly_setcoeff_si(A, 1, b);
		mpz_poly_setcoeff_si(A, 2, c);
		mpz_poly_resultant(N0, f0, A);
		mpz_poly_resultant(N1, f1, A);
		mpz_abs(N0, N0);
		mpz_abs(N1, N1);
		// side 0
		string N0str = line.substr(0, line.find(separator1));
		while (N0str.length()) {
			string pstr = N0str.substr(0, N0str.find(separator2));
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
				int pos = N0str.find(separator2);
				if (pos == npos) pos = pstr.length() - 1;
				N0str.erase(0, pos + 1);
			}
		}
		if (!isrel) {
			cout << "#" << line0 << endl << flush;
			continue;
		}
		// side 1
		line.erase(0, line.find(separator1) + 1);
		string N1str = line.substr(0, line.find(separator1));
		while (N1str.length()) {
			string pstr = N1str.substr(0, N1str.find(separator2));
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
				int pos = N1str.find(separator2);
				if (pos == npos) pos = pstr.length() - 1;
				N1str.erase(0, pos + 1);
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
			string N0str = line.substr(0, line.find(separator1));
			while (N0str.length()) {
				string pstr = N0str.substr(0, N0str.find(separator2));
				if (pstr.length()) {
					mpz_set_str(p, pstr.c_str(), BASE);
					if (mpz_divisible_p(N0, p)) {
						mpz_divexact(N0, N0, p);
						N0primes.push_back(pstr);
					}
					int pos = N0str.find(separator2);
					if (pos == npos) pos = pstr.length() - 1;
					N0str.erase(0, pos + 1);
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
			N0str = ""; int l = N0primes.size();
			for (int i = 0; i < l; i++)
				N0str += N0primes[i] + (i < l - 1 ? "," : "");
			
			// side 1
			N1primes.clear();
			line.erase(0, line.find(separator1) + 1);
			string N1str = line.substr(0, line.find(separator1));
			while (N1str.length()) {
				string pstr = N1str.substr(0, N1str.find(separator2));
				if (pstr.length()) {
					mpz_set_str(p, pstr.c_str(), BASE);
					if (mpz_divisible_p(N1, p)) {
						mpz_divexact(N1, N1, p);
						N1primes.push_back(pstr);
					}
					int pos = N1str.find(separator2);
					if (pos == npos) pos = pstr.length() - 1;
					N1str.erase(0, pos + 1);
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
			N1str = ""; l = N1primes.size();
			for (int i = 0; i < l; i++)
				N1str += N1primes[i] + (i < l - 1 ? "," : "");

			string str = to_string(a) + "," + to_string(b) + "," + to_string(c) + ":" + N0str + ":" + N1str;

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


bool known_good_prime(int pt, int* pcache, int pmin, int* sievep, int kmin, int k)
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
	return known_good;
}


