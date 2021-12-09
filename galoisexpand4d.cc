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
#include "mpz_poly_bivariate.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::to_string;
using std::stringstream;
auto npos = std::string::npos;

inline int gcd(int a, int b);
bool known_good_prime(int pt, int* pcache, int pmin, int* sievep, int kmin, int k);

int main (int argc, char** argv)
{
	if ((argc+1)/2 != 3 ) {	// argc must be 5 or 6
		cout << endl << "Usage: ./galoisexpand4d inputpoly factorbase relations genbadmax {nocheck}" << endl;
		return 0;
	}

	bool verbose = false;

	// read input polynomials
	mpz_t* fhtpoly = new mpz_t[20];	// max degree of 20.  Not the neatest
	mpz_t* ghtpoly = new mpz_t[20];	// max degree of 20.  Not the neatest
	mpz_t* hpoly = new mpz_t[20];	// max degree of 20.  Not the neatest
	for (int i = 0; i < 20; i++) {
		mpz_init(fhtpoly[i]);
		mpz_init(ghtpoly[i]);
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
	int degfht = -1;
	if (verbose) cout << endl << "Side 0 polynomial fh_t (ascending coefficients)" << endl;
	while (getline(file, line) && line.substr(0,3) == "fht" ) {
		line = line.substr(line.find_first_of(" ")+1);
		//mpz_set_str(c, line.c_str(), 10);
		mpz_set_str(fhtpoly[++degfht], line.c_str(), 10);
		//mpz_get_str(linebuffer, 10, fhtpoly[degf-1]);
		if (verbose) cout << line << endl;
	}
	//int degf = fhtpoly.size();
	// read other poly
	int degght = -1;
	bool read = true;
	if (verbose) cout << endl << "Side 1 polynomial gh_t: (ascending coefficients)" << endl;
	while (read && line.substr(0,3) == "ght" ) {
		line = line.substr(line.find_first_of(" ")+1);
		//mpz_set_str(c, line.c_str(), 10);
		mpz_set_str(ghtpoly[++degght], line.c_str(), 10);
		//mpz_get_str(linebuffer, 10, ghtpoly[degg-1]);
		if (verbose) cout << line << endl;
		read = static_cast<bool>(getline(file, line));
	}
	//int degg = ghtpoly.size();
	// read other poly
	int degh = -1;
	read = true;
	if (verbose) cout << endl << "Tower polynomial h: (ascending coefficients)" << endl;
	while (read && line.substr(0,1) == "h" ) {
		line = line.substr(line.find_first_of(" ")+1);
		//mpz_set_str(c, line.c_str(), 10);
		mpz_set_str(hpoly[++degh], line.c_str(), 10);
		//mpz_get_str(linebuffer, 10, ghtpoly[degg-1]);
		if (verbose) cout << line << endl;
		read = static_cast<bool>(getline(file, line));
	}
	// read bivariate F0 poly
	mpz_poly_bivariate F0;
	mpz_poly_bivariate_init(F0, 0);	// init to deg 0 (constant)
	mpz_poly F0i;
	mpz_poly_init(F0i, 0); // init to deg 0 (constant)
	mpz_t F0ij; mpz_init(F0ij);
	read = true;
	if (verbose) cout << endl << "Bivariate polynomial F0: (ascending coefficients)" << endl;
	int inow = 0;
	while (read && line.substr(0,1) == "f" ) {
		int u = line.find_first_of("_");
		string ch = line.substr(1, u-1);
		int inew = atoi(ch.c_str());
		ch = line.substr(u+1, line.find_first_of(":")-u-1);
		int  j = atoi(ch.c_str());
		if (inew == inow) {
			line = line.substr(line.find_first_of(" ")+1);
			mpz_set_str(F0ij, line.c_str(), 10);
			mpz_poly_setcoeff(F0i, j, F0ij);
		}
		else {
			mpz_poly_bivariate_setcoeff(F0, inow, F0i);
			inow = inew;
			F0i->deg = 0;
			line = line.substr(line.find_first_of(" ")+1);
			mpz_set_str(F0ij, line.c_str(), 10);
			mpz_poly_setcoeff(F0i, j, F0ij);
		}
		if (verbose) cout << line << endl;
		read = static_cast<bool>(getline(file, line));
	}
	mpz_poly_bivariate_setcoeff(F0, inow, F0i);
	// read bivariate F1 poly
	mpz_poly_bivariate F1;
	mpz_poly_bivariate_init(F1, 0);	// init to deg 0 (constant)
	mpz_poly F1i;
	mpz_poly_init(F1i, 0); // init to deg 0 (constant)
	mpz_t F1ij; mpz_init(F1ij);
	read = true;
	if (verbose) cout << endl << "Bivariate polynomial F1: (ascending coefficients)" << endl;
	inow = 0;
	while (read && line.substr(0,1) == "g" ) {
		int u = line.find_first_of("_");
		string ch = line.substr(1, u-1);
		int inew = atoi(ch.c_str());
		ch = line.substr(u+1, line.find_first_of(":")-u-1);
		int  j = atoi(ch.c_str());
		if (inew == inow) {
			line = line.substr(line.find_first_of(" ")+1);
			mpz_set_str(F1ij, line.c_str(), 10);
			mpz_poly_setcoeff(F1i, j, F1ij);
		}
		else {
			mpz_poly_bivariate_setcoeff(F1, inow, F1i);
			inow = inew;
			F1i->deg = 0;
			line = line.substr(line.find_first_of(" ")+1);
			mpz_set_str(F1ij, line.c_str(), 10);
			mpz_poly_setcoeff(F1i, j, F1ij);
		}
		if (verbose) cout << line << endl;
		read = static_cast<bool>(getline(file, line));
	}
	mpz_poly_bivariate_setcoeff(F1, inow, F1i);
	file.close();
	if (verbose) cout << endl << "Complete.  Degree fh_t = " << degfht << ", degree gh_t = " << degght << "." << endl;

	// load factor base
	std::clock_t start; double timetaken = 0;
	if (verbose) cout << endl << "Loading factor base..." << endl;
	start = clock();
	ifstream fbfile(argv[2]);
	start = clock();
	// read fbb
	getline(fbfile, line);
	int fbb = atoi(line.c_str());
	// compute small primes
	int fbmax = fbb;
	char* sieve = new char[fbmax+1]();
	int* primes = new int[2000000];	// hardcoded for the moment.  Allows fbb of 2^23 easily
	for (int i = 2; i <= sqrt(fbmax); i++)
		if(!sieve[i])
			for (int j = i*i; j <= fbmax; j += i)
				if(!sieve[j]) sieve[j] = 1;
	int nump = 0;
	for (int i = 2; i <= fbmax-1; i++)
		if (!sieve[i])
			primes[nump++] = i;

	// set up constants
	mpz_t r0; mpz_init(r0);
	int* sieves0 = new int[degfht * nump]();
	int* sievep0 = new int[nump]();
	int* sieves1 = new int[degght * nump]();
	int* sievep1 = new int[nump]();
	int* sieveP0 = new int[nump]();
	int* sieveP1 = new int[nump]();
	int* sieveS0 = new int[degfht * nump]();
	int* sieveS1 = new int[degght * nump]();
	//int* sievenum_s0modp = new int[nump]();
	//int* sievenum_s1modp = new int[nump]();
	int* sievenum_S0modp = new int[nump]();
	int* sievenum_S1modp = new int[nump]();

	// read k0
	getline(fbfile, line);
	int k0 = atoi(line.c_str());
	for (int i = 0; i < k0; i++) {
		getline(fbfile, line);
		stringstream ss(line);
		string substr;
		getline(ss, substr, ',');
		sieveP0[i] = atoi(substr.c_str());
		int j = 0;
		while( ss.good() ) {
			getline( ss, substr, ',' );
			sieveS0[i*degfht + j++] = atoi(substr.c_str());
		}
		sievenum_S0modp[i] = j;
	}
	// read kh0
	getline(fbfile, line);
	int kh0 = atoi(line.c_str());
	for (int i = 0; i < kh0; i++) {
		getline(fbfile, line);
		stringstream ss(line);
		string substr;
		getline(ss, substr, ',');
		sievep0[i] = atoi(substr.c_str());
		int j = 0;
		while( ss.good() ) {
			getline( ss, substr, ',' );
			sieves0[i*degfht + j++] = atoi(substr.c_str());
		}
		//sievenum_s0modp[i] = j;  // commented out since same as sievenum_S0modp[i]
	}
	// read k1
	getline(fbfile, line);
	int k1 = atoi(line.c_str());
	for (int i = 0; i < k1; i++) {
		getline(fbfile, line);
		stringstream ss(line);
		string substr;
		getline(ss, substr, ',');
		sieveP1[i] = atoi(substr.c_str());
		int j = 0;
		while( ss.good() ) {
			getline( ss, substr, ',' );
			sieveS1[i*degght + j++] = atoi(substr.c_str());
		}
		sievenum_S1modp[i] = j;
	}
	// read kh1
	getline(fbfile, line);
	int kh1 = atoi(line.c_str());
	for (int i = 0; i < kh1; i++) {
		getline(fbfile, line);
		stringstream ss(line);
		string substr;
		getline(ss, substr, ',');
		sievep1[i] = atoi(substr.c_str());
		int j = 0;
		while( ss.good() ) {
			getline( ss, substr, ',' );
			sieves1[i*degght + j++] = atoi(substr.c_str());
		}
		//sievenum_s1modp[i] = j;  // commented out since same as sievenum_S1modp[i]
	}
	timetaken += ( clock() - start ) / (double) CLOCKS_PER_SEC;
	if (verbose) cout << "Complete.  Time taken: " << timetaken << "s" << endl;
	if (verbose) cout << "There are " << k0 << " factor base primes on side 0." << endl;
	if (verbose) cout << "There are " << k1 << " factor base primes on side 1." << endl;
	fbfile.close();

	int64_t* froots = new int64_t[degfht];
	int64_t* groots = new int64_t[degght];

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

    int64_t genbadmax = 1000;
    genbadmax = atoi(argv[4]);

	bool nocheck = false;
	string argv5("");
    if (argc == 6) argv5.assign(argv[5], strlen(argv[5]));
	if (argv5 == "nocheck") nocheck = true;

	// read relations file
	mpz_poly f0; mpz_poly f1; mpz_poly A;
	mpz_poly_init(f0, degfht); mpz_poly_init(f1, degght); mpz_poly_init(A, 3);
	mpz_poly_set_mpz(f0, fhtpoly, degfht);
	mpz_poly_set_mpz(f1, ghtpoly, degght);
	mpz_poly h0; mpz_poly_init(h0, degh);
	mpz_poly Fqh_x; mpz_poly_init(Fqh_x, 0);
	mpz_poly_bivariate Aq; mpz_poly_bivariate_init(Aq, 0);
	mpz_poly Aq0; mpz_poly_init(Aq0, 0);
	for (int i = 0; i <= degh; i++) mpz_poly_setcoeff(h0, i, hpoly[i]);
	mpz_t N0; mpz_init(N0); mpz_t N1; mpz_init(N1);
	mpz_t p; mpz_init(p); int BASE = 16;
	vector<int> q; vector<int> e;
	vector<string> N0primes;
	vector<string> N1primes;
	ifstream rels(argv[3]);
	string separator1 = ":";
	string separator2 = ",";
	while (std::getline(rels, line)) { //rels >> line) {
		string line0 = line;
		bool isrel = true;

		if (line.substr(0, 1) == "#") continue;

		string Astr = line.substr(0, line.find(separator1));
		line.erase(0, Astr.length() + 1);
		int a = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		int b = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		int c = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		int d = atoi(Astr.substr(0, Astr.find(separator2)).c_str());

		// check if relation is valid, i.e. all primes are factor base primes
		mpz_poly_setcoeff_si(Aq0, 0, a);
		mpz_poly_setcoeff_si(Aq0, 1, b);	// Aq0 = x - R[ll]
		mpz_poly_bivariate_setcoeff(Aq, 0, Aq0);
		mpz_poly_setcoeff_si(Aq0, 0, c);
		mpz_poly_setcoeff_si(Aq0, 1, d);	// Aq0 = x - R[ll]
		mpz_poly_bivariate_setcoeff(Aq, 1, Aq0);
		// reset Fqh_x, which might need to go from e.g. deg 8 to deg 4
		//Fqh_x->deg = 0;
		Fqh_x->deg = -1; mpz_poly_cleandeg(Fqh_x, 0);
		mpz_poly_bivariate_resultant_y(Fqh_x, F0, Aq);
		mpz_poly_resultant(N0, Fqh_x, h0);
		mpz_poly_bivariate_resultant_y(Fqh_x, F1, Aq);
		mpz_poly_resultant(N1, Fqh_x, h0);
		mpz_abs(N0, N0);
		mpz_abs(N1, N1);
		// side 0
		string N0str = line.substr(0, line.find(separator1));
		while (N0str.length()) {
			string pstr = N0str.substr(0, N0str.find(separator2));
			if (pstr.length() == 0) {
				N0str.erase(0, 1);
			}
			else {
				mpz_set_str(p, pstr.c_str(), BASE);
				assert(mpz_divisible_p(N0, p));
				mpz_divexact(N0, N0, p);
				if (!nocheck) {
					if (mpz_cmp_ui(p, fbmax) < 0) {
						int pt = mpz_get_ui(p);
						if (!known_good_prime(pt, p0cache, p0min, sievep0, k0min, k0)) { // && pt > genbadmax) { // relation is bad
							isrel = false;
							break;
						}
					}
					else { // make sure large prime is good
						int pt = mpz_get_ui(p);
						//for (int i = 0; i <= degfht; i++) fmodp[i] = mod(fi64[i], pt);
						for (int i = 0; i <= degfht; i++)
							fmodp[i] = mpz_mod_ui(r0, fhtpoly[i], pt);
						int nr = polrootsmod(fmodp, degfht, froots, pt);
						if (nr == 0) {
							isrel = false;
							break;
						}
					}
				}
				int pos = N0str.find(separator2);
				if (pos == npos) pos = pstr.length() - 1;
				N0str.erase(0, pos + 1);
			}
		}
		if (!isrel) {
			cout << "#" << line0 << endl;
			continue;
		}
		// side 1
		line.erase(0, line.find(separator1) + 1);
		string N1str = line.substr(0, line.find(separator1));
		while (N1str.length()) {
			string pstr = N1str.substr(0, N1str.find(separator2));
			if (pstr.length() == 0) {
				N1str.erase(0, 1);
			}
			else {
				mpz_set_str(p, pstr.c_str(), BASE);
				assert(mpz_divisible_p(N1, p));
				mpz_divexact(N1, N1, p);
				if (!nocheck) {
					if (mpz_cmp_ui(p, fbmax) < 0) {
						int pt = mpz_get_ui(p);
						if (!known_good_prime(pt, p1cache, p1min, sievep1, k1min, k1)) { // && pt > genbadmax) { // relation is bad
							isrel = false;
							break;
						}
					}
					else { // make sure large prime is good
						int pt = mpz_get_ui(p);
						//for (int i = 0; i <= degght; i++) gmodp[i] = mod(gi64[i], pt);
						for (int i = 0; i <= degght; i++)
							gmodp[i] = mpz_mod_ui(r0, ghtpoly[i], pt);
						int nr = polrootsmod(gmodp, degght, groots, pt);
						if (nr == 0) {
							isrel = false;
							break;
						}
					}
				}
				int pos = N1str.find(separator2);
				if (pos == npos) pos = pstr.length() - 1;
				N1str.erase(0, pos + 1);
			}
		}
		if (!isrel) {
			cout << "#" << line0 << endl;
			continue;
		}

		// print original relation
		cout << line0 << endl;
		
		// compute Galois conjugates
		int a0 = a; int b0 = b; int c0 = c; int d0 = d;
		for (int g = 1; g <= 1; g++) {
			switch(g) {
				case 1:
					a = c0; b = d0; c = a0; d = b0; break;
				//case 2:
				//	a = -3*a0+b0; b = -a0; c = -3*c0+d0; d = -c0; break;
				//case 3:
				//	a = -3*c0+d0; b = -c0; c = -3*a0+b0; d = -a0; break;
			}

			// restore line
			line = line0;
			line.erase(0, line.find(separator1) + 1);

			// compute norms
			mpz_poly_setcoeff_si(Aq0, 0, a);
			mpz_poly_setcoeff_si(Aq0, 1, b);	// Aq0 = x - R[ll]
			mpz_poly_bivariate_setcoeff(Aq, 0, Aq0);
			mpz_poly_setcoeff_si(Aq0, 0, c);
			mpz_poly_setcoeff_si(Aq0, 1, d);	// Aq0 = x - R[ll]
			mpz_poly_bivariate_setcoeff(Aq, 1, Aq0);
			// reset Fqh_x, which might need to go from e.g. deg 8 to deg 4
			//Fqh_x->deg = 0;
			Fqh_x->deg = -1; mpz_poly_cleandeg(Fqh_x, 0);
			mpz_poly_bivariate_resultant_y(Fqh_x, F0, Aq);
			mpz_poly_resultant(N0, Fqh_x, h0);
			mpz_poly_bivariate_resultant_y(Fqh_x, F1, Aq);
			mpz_poly_resultant(N1, Fqh_x, h0);
			mpz_abs(N0, N0);
			mpz_abs(N1, N1);
			//cout << "[a, b, c] = [" << a << ", " << b << ", " << c << "]" << endl;

			// side 0
			N0primes.clear();
			string N0str = line.substr(0, line.find(separator1));
			while (N0str.length()) {
				string pstr = N0str.substr(0, N0str.find(separator2));
				if (pstr.length() == 0) {
					N0str.erase(0, 1);
				}
				else {
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
			if (mpz_sizeinbase(N0, 2) > 63) { cout << "Error - Galois conjugate not smooth!" << endl; return 0; }
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
				if (pstr.length() == 0) {
					N1str.erase(0, 1);
				}
				else {
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
			if (mpz_sizeinbase(N1, 2) > 63) { cout << "Error - Galois conjugate not smooth!" << endl; return 0; }
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

			string str = to_string(a) + "," + to_string(b) + ","
					   + to_string(c) + "," + to_string(d) + ":" + N0str + ":" + N1str;

			cout << str << endl;
		}
	}

	mpz_clear(p);
	mpz_clear(N1); mpz_clear(N0);
    mpz_poly_clear(A); mpz_poly_clear(f1); mpz_poly_clear(f0);
	delete[] sievenum_S1modp;
	delete[] sievenum_S0modp;
	//delete[] sievenum_s1modp;
	//delete[] sievenum_s0modp;
	delete[] sieveS1;
	delete[] sieveS0;
	delete[] sieveP1;
	delete[] sieveP0;
	delete[] sievep1;
	delete[] sieves1;
	delete[] sievep0;
	delete[] sieves0;
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
		mpz_clear(ghtpoly[i]);
		mpz_clear(fhtpoly[i]);
	}
	delete[] hpoly;
	delete[] ghtpoly;
	delete[] fhtpoly;

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
	for (;;) {
		int i_bak = i;
		if (i >= min) {
			if (pt == sievep[i]) { known_good = true; break; }
			else if (pt < sievep[i]) max = i;
            else if (min == i) break;
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
        if (min == max) break;  // couldn't find pt in sievep if this is true
	}
	return known_good;
}


