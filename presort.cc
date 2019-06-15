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
#include <vector>

using std::cout;
using std::endl;
using std::flush;
using std::ifstream;
using std::string;
using std::to_string;
using std::stringstream;
using std::vector;
using std::sort;
auto npos = std::string::npos;

inline int gcd(int a, int b);
bool known_good_prime(int pt, int* pcache, int pmin, int* sievep, int kmin, int k);

typedef struct relstruct {
	bool valid;
	int a;
	int b;
	int c;
	vector<int> side0p;
	vector<int> side1p;
} relation;

int main (int argc, char** argv)
{
	if (argc == 1) {
		cout << endl << "Usage: ./presort inputpoly factorbase relations" << endl << flush;
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
	int* p0cache = new int[p0min]();	// set to 0 (false)
	int c0 = 0; int k0min = 0; 
	for (int p = 2; p < p0min; p++) {
		int t = c0;
		while (sievep0[t] < p) t++;
		if (sievep0[t] == p) { c0 = t + 1; p0cache[p] = 1; k0min++; }	// factor base prime
	}
	int p1min = 1000;
	int* p1cache = new int[p1min]();	// set to 0 (false);
	int c1 = 0; int k1min = 0;
	for (int p = 2; p < p1min; p++) {
		int t = c1;
		while (sievep1[t] < p) t++;
		if (sievep1[t] == p) { c1 = t + 1; p1cache[p] = 1; k1min++; }	// factor base prime
	}

	// read relations file
	mpz_poly f0; mpz_poly f1;
	mpz_poly_init(f0, degf); mpz_poly_init(f1, degg);
	mpz_poly_set_mpz(f0, fpoly, degf);
	mpz_poly_set_mpz(f1, gpoly, degg);
	mpz_t p; mpz_init(p); int BASE = 16;
	ifstream rels(argv[3]);
	string separator1 = ":";
	string separator2 = ",";
	vector<relation> relations;
	int n = 0;
	while (rels >> line) {
		string line0 = line;
		n++;
		if (n % 100000 == 0) cout << n << " lines read" << endl;
		bool isrel = true;

		string Astr = line.substr(0, line.find(separator1));
		line.erase(0, Astr.length() + 1);
		int a = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		int b = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		int c = atoi(Astr.substr(0, Astr.find(separator2)).c_str());

		relation rel;
		rel.valid = true; rel.a = a; rel.b = b; rel.c = c;

		// side 0
		string side0 = line.substr(0, line.find(separator1));
		while (side0.length()) {
			string pstr = side0.substr(0, side0.find(separator2));
			if (pstr.length()) {
				mpz_set_str(p, pstr.c_str(), BASE);
				int pt = mpz_get_ui(p);
				/*if (mpz_cmp_ui(p, fbmax) < 0) {
					if (!known_good_prime(pt, p0cache, p0min, sievep0, k0min, k0)) { // relation is bad
						isrel = false;
						break;
					}
				}
				else { // make sure large prime is good
					for (int i = 0; i <= degf; i++) fmodp[i] = mod(fi64[i], pt);
					// the following line need not be called very often
					int nr = polrootsmod(fmodp, degf, froots, pt);
					if (nr == 0) {
						isrel = false;
						break;
					}
				}*/
				rel.side0p.push_back(pt);
				int pos = side0.find(separator2);
				if (pos == npos) pos = pstr.length() - 1;
				side0.erase(0, pos + 1);
			}
		}
		if (!isrel) {
			rel.valid = false;
			//cout << "#" << line0 << endl << flush;
			//continue;
		}
		sort(rel.side0p.begin(), rel.side0p.end());
		
		// side 1
		line.erase(0, line.find(separator1) + 1);
		string side1 = line.substr(0, line.find(separator1));
		while (side1.length()) {
			string pstr = side1.substr(0, side1.find(separator2));
			if (pstr.length()) {
				mpz_set_str(p, pstr.c_str(), BASE);
				int pt = mpz_get_ui(p);
				/*if (mpz_cmp_ui(p, fbmax) < 0) {
					if (!known_good_prime(pt, p1cache, p1min, sievep1, k1min, k1)) { // relation is bad
						isrel = false;
						break;
					}
				}
				else { // make sure large prime is good
					for (int i = 0; i <= degg; i++) gmodp[i] = mod(gi64[i], pt);
					// the following line need not be called very often
					int nr = polrootsmod(gmodp, degg, groots, pt);
					if (nr == 0) {
						isrel = false;
						break;
					}
				}*/
				rel.side1p.push_back(pt);
				int pos = side1.find(separator2);
				if (pos == npos) pos = pstr.length() - 1;
				side1.erase(0, pos + 1);
			}
		}
		if (!isrel) {
			rel.valid = false;
			//cout << "#" << line0 << endl << flush;
			//continue;
		}
		sort(rel.side1p.begin(), rel.side1p.end());

		// store relation
		relations.push_back(rel);
	}

	mpz_clear(p);
    mpz_poly_clear(f1); mpz_poly_clear(f0);
	delete[] primes;
	delete[] sieve;
	delete[] groots;
	delete[] froots;
	delete[] gmodp;
	delete[] fmodp;
	delete[] gi64;
	delete[] fi64;

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


