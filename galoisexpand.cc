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

using std::cout;
using std::endl;
using std::flush;
using std::ifstream;
using std::string;
using std::to_string;
auto npos = std::string::npos;

inline int gcd(int a, int b);

int main (int argc, char** argv)
{
	if (argc == 1) {
		cout << endl << "Usage: ./galoisexpand inputpoly fbb relations" << endl << flush;
		return 0;
	}

	bool verbose = false;

	// read input polynomial
	mpz_t* fpoly = new mpz_t[20];	// max degree of 20.  Not the neatest
	mpz_t* gpoly = new mpz_t[20];	// max degree of 20.  Not the neatest
	mpz_t* hpoly = new mpz_t[20];	// max degree of 20.  Not the neatest
	for (int i = 0; i < 20; i++) {
		mpz_init(fpoly[i]);
		mpz_init(gpoly[i]);
		mpz_init(hpoly[i]);
	}
	string line;
	char linebuffer[100];
	ifstream file(argv[1]);
	getline(file, line);	// first line contains number n to factor
	// read nonlinear poly
	int degf = -1;
	if (verbose) cout << endl << "Side 0 polynomial f0 (ascending coefficients)" << endl;
	while (getline(file, line) && line.substr(0,1) == "c" ) {
		line = line.substr(line.find_first_of(" ")+1);
		//mpz_set_str(c, line.c_str(), 10);
		mpz_set_str(fpoly[++degf], line.c_str(), 10);
		//mpz_get_str(linebuffer, 10, fpoly[degf-1]);
		if (verbose) cout << line << endl << flush;
	}
	//int degf = fpoly.size();
	// read other poly
	int degg = -1;
	bool read = true;
	if (verbose) cout << endl << "Side 1 polynomial f1: (ascending coefficients)" << endl;
	while (read && line.substr(0,1) == "Y" ) {
		line = line.substr(line.find_first_of(" ")+1);
		//mpz_set_str(c, line.c_str(), 10);
		mpz_set_str(gpoly[++degg], line.c_str(), 10);
		//mpz_get_str(linebuffer, 10, gpoly[degg-1]);
		if (verbose) cout << line << endl << flush;
		read = static_cast<bool>(getline(file, line));
	}
	//int degg = gpoly.size();
	file.close();
	if (verbose) cout << endl << "Complete.  Degree f0 = " << degf << ", degree f1 = " << degg << "." << endl;

	// compute small primes
	if (verbose) cout << endl << "Starting sieve of Eratosthenes for small primes..." << endl << flush;
	int fbbits = atoi(argv[2]);
	int max = 1<<fbbits;
	char* sieve = new char[max+1]();
	int* primes = new int[1077871];	// hardcoded for the moment.  Allows fbb of 2^23 easily
	for (int i = 2; i <= sqrt(max); i++)
		if(!sieve[i])
			for (int j = i*i; j <= max; j += i)
				if(!sieve[j]) sieve[j] = 1;
	int nump = 0;
	for (int i = 2; i <= max-1; i++)
		if (!sieve[i])
			primes[nump++] = i;
	if (verbose) cout << "Complete." << endl;

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

		string Astr = line.substr(0, line.find(separator1));
		line.erase(0, Astr.length() + 1);
		int a = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		int b = atoi(Astr.substr(0, Astr.find(separator2)).c_str());  Astr.erase(0, Astr.find(separator2) + 1);
		int c = atoi(Astr.substr(0, Astr.find(separator2)).c_str());

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


