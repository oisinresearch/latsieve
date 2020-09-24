#include <cstdlib>
#include <stdint.h>	// int64_t
#include <iostream> // cout
#include <iomanip> // setprecision
#include "L2lu64.h"
#include <gmpxx.h>
#include "intpoly.h"
#include <cmath>	// sqrt
#include <fstream>	// file
#include <ctime>	// clock_t
#include <cstring>	// memset
#include <omp.h>
#include "mpz_poly.h"
#include <sstream>	// stringstream
#include <stack>	// stack
#include "mpz_poly_bivariate.h"

using std::cout;
using std::endl;
using std::flush;
//using std::vector;
using std::string;
using std::ifstream;
using std::fixed;
using std::scientific;
using std::setprecision;
using std::sort;
using std::to_string;
using std::hex;
using std::stringstream;
using std::stack;
using std::abs;

struct keyval {
	int id;
	uint8_t logp;
};

__int128 MASK64;

int latsieve4d(mpz_poly_bivariate F, int64_t* h, int degh, int64_t* fh_t, int degfh_t,
		int64_t q, int l, int ll, int* allp, int nump, int* s, int* S, int* num_smodp,
		int* num_Smodp, keyval* M, int Mlen, int* B);
void histogram(keyval*M, uint8_t* H, int len);
bool lattice_sorter(keyval const& kv1, keyval const& kv2);
void csort(keyval* M, keyval* L, int* H, int len);
inline int floordiv(int a, int b);
inline int64_t modinv(int64_t x, int64_t m);
inline int nonzerolcm(int u1, int u2, int u3);
inline int nonzerolcm4d(int u1, int u2, int u3, int u4);
inline int gcd(int a, int b);
inline __int128 gcd128(__int128 a, __int128 b);
inline void getab(int u1, int u2, int u3, int v1, int v2, int v3, int x, int y, int z, 
				int B1, int B2, int B3, int* a, int* b);
inline void getab4d(int u1, int u2, int u3, int u4, int v1, int v2, int v3, int v4,
				int x, int y, int z, int t, int B1, int B2, int B3, int B4, int* a, int* b);
inline bool planeintersectsbox(int nx, int ny, int nz, int x, int y, int z, int B1, int B2, int B3);
inline bool spaceintersectsbox4d(int nx, int ny, int nz, int nt, int x, int y, int z, int t,
							int B1, int B2, int B3, int B4);
bool planeintersectsbox4d(int u1, int u2, int u3, int u4, int v1, int v2, int v3, int v4,
		int w1, int w2, int w3, int w4, int B1, int B2, int B3, int B4);
inline int minnonneg(int u, int v, int w);
inline int minnonneg4d(int u, int v, int w, int t);
inline int minabs(int u, int v, int w);
inline int maxabs(int u, int v, int w);
inline int min(int u, int v, int w);
inline int max(int u, int v, int w);
void GetlcmScalar(int B, mpz_t S, int* primes, int nump);
inline __int128 make_int128(uint64_t lo, uint64_t hi);
bool PollardPm1(mpz_t N, mpz_t S, mpz_t factor);
bool PollardPm1_mpz(mpz_t N, mpz_t S, mpz_t factor);
bool PollardPm1_int128(__int128 N, mpz_t S, __int128 &factor, int64_t mulo, int64_t muhi);
bool EECM(mpz_t N, mpz_t S, mpz_t factor, int d, int a, int X0, int Y0, int Z0);
bool EECM_mpz(mpz_t N, mpz_t S, mpz_t factor, int d, int a, int X0, int Y0, int Z0);
bool EECM_int128(__int128 N, mpz_t S, __int128 &factor, int d, int a, int X0, int Y0, int Z0, int64_t mulo, int64_t muhi);


int main(int argc, char** argv)
{
	// set constant
	MASK64 = 1L;
	MASK64 = MASK64 << 64;
	MASK64 = MASK64 - 1L;

	//cout << (uint64_t)(MASK64) << " " << (uint64_t)(MASK64 >> 64) << endl;

	if (argc != 15) {
		cout << endl << "Usage: ./latsieve4d inputpoly fbbits factorbasefile B1 B2 B3 B4 qmin qmax th0 th1 lpbbits cofacscalar qside" << endl << endl;
		return 0;
	}

	cout << "# ";
	for (int i = 0; i < argc; i++) cout << argv[i] << " ";
	cout << endl << flush;

	bool verbose = false; verbose = true;
		
	if (verbose) cout << endl << "Reading input polynomial in file " << argv[1] << "..." << flush;
	//vector<mpz_class> fpoly;
	//vector<mpz_class> gpoly;
	mpz_t* fhtpoly = new mpz_t[20];	// max degree of 20.  Not the neatest
	mpz_t* ghtpoly = new mpz_t[20];	// max degree of 20.  Not the neatest
	mpz_t* hpoly = new mpz_t[20];	// max degree of 20.  Not the neatest
	for (int i = 0; i < 20; i++) {
		mpz_init(fhtpoly[i]);
		mpz_init(ghtpoly[i]);
		mpz_init(hpoly[i]);
	}
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
		//mpz_get_str(linebuffer, 10, fpoly[degf-1]);
		if (verbose) cout << line << endl << flush;
	}
	//int degf = fpoly.size();
	// read other poly
	int degght = -1;
	bool read = true;
	if (verbose) cout << endl << "Side 1 polynomial gh_t: (ascending coefficients)" << endl;
	while (read && line.substr(0,3) == "ght" ) {
		line = line.substr(line.find_first_of(" ")+1);
		//mpz_set_str(c, line.c_str(), 10);
		mpz_set_str(ghtpoly[++degght], line.c_str(), 10);
		//mpz_get_str(linebuffer, 10, gpoly[degg-1]);
		if (verbose) cout << line << endl << flush;
		read = static_cast<bool>(getline(file, line));
	}
	//int degg = gpoly.size();
	// read other poly
	int degh = -1;
	read = true;
	if (verbose) cout << endl << "Tower polynomial h: (ascending coefficients)" << endl;
	while (read && line.substr(0,1) == "h" ) {
		line = line.substr(line.find_first_of(" ")+1);
		//mpz_set_str(c, line.c_str(), 10);
		mpz_set_str(hpoly[++degh], line.c_str(), 10);
		//mpz_get_str(linebuffer, 10, gpoly[degg-1]);
		if (verbose) cout << line << endl << flush;
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
		if (verbose) cout << line << endl << flush;
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
		if (verbose) cout << line << endl << flush;
		read = static_cast<bool>(getline(file, line));
	}
	mpz_poly_bivariate_setcoeff(F1, inow, F1i);
	file.close();
	//mpz_clear(c);
	if (verbose) cout << endl << "Complete.  Degree fh_t = " << degfht << ", degree gh_t = " << degght << "." << endl;

	if (verbose) cout << endl << "Starting sieve of Eratosthenes for small primes..." << endl << flush;
	int fbbits = 21;
	if (argc >=3) fbbits = atoi(argv[2]);
	int max = 1<<fbbits; // 10000000;// 65536;
	char* sieve = new char[max+1]();
	int* primes = new int[2097152]; //int[1077871]; // int[155611]; //new int[809228];	//new int[6542]; 	// 2039 is the 309th prime, largest below 2048
	for (int i = 2; i <= sqrt(max); i++)
		if(!sieve[i])
			for (int j = i*i; j <= max; j += i)
				if(!sieve[j]) sieve[j] = 1;
	int nump = 0;
	for (int i = 2; i <= max-1; i++)
		if (!sieve[i])
			primes[nump++] = i;
	if (verbose) cout << "Complete." << endl;

	// set up constants
	std::clock_t start; double timetaken = 0;
	mpz_t r0; mpz_init(r0);
	int* sieves0 = new int[degh * nump]();
	int* sievep0 = new int[nump]();
	int* sieves1 = new int[degh * nump]();
	int* sievep1 = new int[nump]();
	int* sieveP0 = new int[nump]();
	int* sieveP1 = new int[nump]();
	int* sieveS0 = new int[degfht * nump]();
	int* sieveS1 = new int[degght * nump]();
	int* sievenum_s0modp = new int[nump]();
	int* sievenum_s1modp = new int[nump]();
	int* sievenum_S0modp = new int[nump]();
	int* sievenum_S1modp = new int[nump]();
	// load factor base
	if (verbose) cout << endl << "Loading factor base..." << endl << flush;
	start = clock();
	ifstream fbfile(argv[3]);
	start = clock();
	// read fbb
	getline(fbfile, line);
	int fbb = atoi(line.c_str());
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
			sieves0[i*degh + j++] = atoi(substr.c_str());
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
			sieves1[i*degh + j++] = atoi(substr.c_str());
		}
		sievenum_s1modp[i] = j;
	}
	timetaken += ( clock() - start ) / (double) CLOCKS_PER_SEC;
	if (verbose) cout << "Complete.  Time taken: " << timetaken << "s" << endl << flush;
	if (verbose) cout << "There are " << k0 << " factor base primes on side 0." << endl << flush;
	if (verbose) cout << "There are " << k1 << " factor base primes on side 1." << endl << flush;
	fbfile.close();
	
	int B[4] = { 7, 7, 7, 7 };
	if (argc >= 5) B[0] = atoi(argv[4]);
	if (argc >= 6) B[1] = atoi(argv[5]);
	if (argc >= 7) B[2] = atoi(argv[6]);
	if (argc >= 8) B[3] = atoi(argv[7]);
	int B1bits = B[0]; int B2bits = B[1]; int B3bits = B[2]; int B4bits = B[3];
	int B1 = 1<<B1bits; int B2 = 1<<B2bits; int B3 = 1<<B3bits; int B4 = 1<<B4bits;
	int Mlen = (B1*2*B2*2*B3*2*B4);	// require positive z coordinate
	Mlen = (int)(3.3f * Mlen);	// upper bound on number of vectors in sieve box
	keyval* M = new keyval[Mlen];	// lattice { id, logp } pairs
	//keyval* L = new keyval[Mlen];	// copy of M
	uint8_t* H = new uint8_t[Mlen];	// histogram
	vector<int> rel;
	// clear M
	//cout << "Clearing memory..." << endl << flush;
	//start = clock();
	//for (int j = 0; j < Mlen; j++) M[j] = (keyval){ 0, 0 };
	//memset(M, 0, Mlen * sizeof(keyval));
	//timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
	//cout << "Memory cleared. Time taken: " << timetaken << "s" << endl << flush;

	int64_t qmin; int64_t qmax; mpz_t qmpz; mpz_init(qmpz);
	mpz_t* pi = new mpz_t[8]; for (int i = 0; i < 8; i++) mpz_init(pi[i]);
	int64_t* r = new int64_t[degh]();
	int64_t* R = new int64_t[degfht]();
	int64_t* fhqt = new int64_t[degfht+1](); // careful, what if degght > degfht?
	int64_t* h = new int64_t[degh+1]();
	mpz_poly f0; mpz_poly f1; mpz_poly A;
	mpz_poly_init(f0, degfht); mpz_poly_init(f1, degght); mpz_poly_init(A, 4);
	mpz_poly_set_mpz(f0, fhtpoly, degfht);
	mpz_poly_set_mpz(f1, ghtpoly, degght);
	mpz_t N0; mpz_t N1;
	mpz_init(N0); mpz_init(N1);
	stringstream stream;
	mpz_t lpb; mpz_init(lpb);
	mpz_t factor; mpz_init(factor); mpz_t p1; mpz_t p2; mpz_init(p1); mpz_init(p2);
	mpz_t tmp1; mpz_init(tmp1); 
	if (argc >= 9) qmin = strtoll(argv[8], NULL, 10);	// atoi(argv[7]);
	if (argc >= 10) qmax = strtoll(argv[9], NULL, 10);	// atoi(argv[8]);
	uint8_t th0 = 70;
	if (argc >= 11) th0 = atoi(argv[10]);
	uint8_t th1 = 70;
	if (argc >= 12) th1 = atoi(argv[11]);
	int lpbits = 29;
	if (argc >= 13) lpbits = atoi(argv[12]);
	int cofacS = 1000;
	if (argc >= 14) cofacS = atoi(argv[13]);
	mpz_t S; mpz_init(S); GetlcmScalar(cofacS, S, primes, 669);	// max S = 5000
	char* str2 = (char*)malloc(20*sizeof(char));
	int qside = atoi(argv[14]);
	mpz_poly Fqh_x; mpz_poly_init(Fqh_x, 0);
	mpz_poly_bivariate Aq; mpz_poly_bivariate_init(Aq, 0);
	mpz_poly Aq0; mpz_poly_init(Aq0, 0);
	mpz_t res; mpz_init(res);

	int64_t q = qmin;
	while (q < qmax) {
		mpz_set_ui(qmpz, q);
		mpz_nextprime(qmpz, qmpz);
		q = mpz_get_ui(qmpz);

		for (int i = 0; i <= degh; i++) h[i] = mpz_mod_ui(r0, hpoly[i], q);

		mpz_poly_bivariate* Fq;
		mpz_t* fhtpolyside = fhtpoly;
		int degfhqt = degfht;
		Fq = &F0;
		if (qside == 1) { fhtpolyside = ghtpoly; degfhqt = degght; Fq = &F1; }
		for (int i = 0; i <= degfht; i++) fhqt[i] = mpz_mod_ui(r0, fhtpolyside[i], q);
		int numl = polrootsmod(h, degh, r, q);
		int numll = polrootsmod(fhqt, degfhqt, R, q);
		if (numl == 0 || numll == 0 || q > qmax) continue;
		
		// sieve side 0
		cout << "# Starting sieve on side 0";
		if (qside == 0) cout << " for special-q " << q;
		cout << "..." << endl << flush;
		start = clock();
		//int m = latsieve3d(fq, degf, q, 0, sievep0, k0, sieves0, sievenum_s0modp, M, Mlen, B);
		int m = 0;
		// compute which r[l] are valid for given R[ll]
		int ll = 0; int l = 0;
		for (int i = 0; i < numl; i++) {
			//Fqh_x->deg = 0; mpz_poly_setcoeff_ui(Fqh_x, 0, 0);
			mpz_poly_setcoeff_ui(Aq0, 1, 1);	// Aq0 = x - R[ll]
			mpz_poly_setcoeff_si(Aq0, 0, -R[ll]);
			mpz_poly_bivariate_setcoeff(Aq, 0, Aq0);
			mpz_poly_bivariate_resultant_y(Fqh_x, *Fq, Aq);
			mpz_poly_eval_ui(res, Fqh_x, r[l]);
			if (mpz_mod_ui(r0, res, q) == 0) { l = i; break; } 
		}
		m = latsieve4d(F0, h, degh, fhqt, degfhqt, q, l, ll, sievep0, k0, sieves0, sieveS0, 
								sievenum_s0modp, sievenum_S0modp, M, Mlen, B);
		timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
		cout << "# Finished! Time taken: " << timetaken << "s" << endl << flush;
		cout << "# Size of lattice point list is " << m << "." << endl << flush;
		cout << "# Constructing histogram..." << endl << flush;
		start = clock();
		//std::stable_sort(M, M + m, &lattice_sorter);
		histogram(M, H, m);
		timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
		cout << "# Finished! Time taken: " << timetaken << "s" << endl << flush;
		int R0 = 0;
		int B1xB2 = B1*B2;
		int B2x2bits = B2bits + 1;
		int B2x2 = 2*B2;
		int B3x2 = 2*B3;
		int B1xB2x2bits = B1bits + B2bits + 1;
		int B1xB2x2xB3x2bits = B1bits + B2bits + 1 + B3bits + 1;
		rel.clear();
		for (int i = 0; i < m; i++) {
			if (H[i] > th0) {
				rel.push_back(i);
				R0++;
			}
		}
		cout << "# " << R0 << " candidates on side 0." << endl << flush;
		// sieve side 1
		cout << "# Starting sieve on side 1";
		if (qside == 1) cout << " for special-q " << q;
		cout << "..." << endl << flush;
		start = clock();
		//m = latsieve3d(fq, degg, q, 0, sievep1, k1, sieves1, sievenum_s1modp, M, Mlen, B);
		m = latsieve4d(F1, h, degh, fhqt, degfhqt, q, l, ll, sievep1, k1, sieves1, sieveS1, 
								sievenum_s1modp, sievenum_S1modp, M, Mlen, B);
		timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
		cout << "# Finished! Time taken: " << timetaken << "s" << endl << flush;
		cout << "# Size of lattice point list is " << m << "." << endl << flush;
		cout << "# Constructing histogram..." << endl << flush;
		start = clock();
		//std::stable_sort(M, M + m, &lattice_sorter);
		histogram(M, H, m);
		timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
		cout << "# Finished! Time taken: " << timetaken << "s" << endl << flush;
		int R1 = 0;
		for (int i = 0; i < m; i++) {
			if (H[i] > th1) {
				rel.push_back(i);
				R1++;
			}
		}
		cout << "# " << R1 << " candidates on side 1." << endl << flush;
		// sort hits
		sort(rel.begin(), rel.end());
		
		// compute special-q lattice L 
		int64_t L[16];
		numl = polrootsmod(h, degh, r, q);
		numll = polrootsmod(fhqt, degfhqt, R, q);
		//int l = 0; int ll = 0;
		L[0]  = q; L[1]  = -r[l]; L[2]  = -R[ll]; L[3]  = 0;
		L[4]  = 0; L[5]  = 1;     L[6]  = 0;      L[7]  = -R[ll];
		L[8]  = 0; L[9]  = 0;     L[10] = 1;      L[11] = 0;
		L[12] = 0; L[13] = 0;     L[14] = 0;      L[15] = 1;
		int64L2(L, 4);	// LLL reduce L, time log(q)^2
		
		// print list of potential relations
		int R = 0;
		for (int i = 0; i < rel.size()-1; i++)
		{
			if (rel[i] == rel[i+1] && rel[i] != 0) {
				int x = rel[i] % B1;
				int y = ((rel[i] >> B1bits) % B2x2) - B2;
				int z = ((rel[i] >> B1xB2x2bits) % B3x2) - B3;
				int t = (rel[i] >> B1xB2x2xB3x2bits) - B4;
				if (x != 0 || y != 0 || z != 0) {
					// compute [a,b,c,d]
					int a = L[0]*x+L[1]*y+L[2]*z+L[3]*t;
					int b = L[4]*x+L[5]*y+L[6]*z+L[7]*t;
					int c = L[8]*x+L[9]*y+L[10]*z+L[11]*t;
					int d = L[12]*x+L[13]*y+L[14]*z+L[15]*t;
					if (R >= 100000 & R < 100030)
						cout << rel[i] << ": " << a << "," << b << "," << c << "," << d << endl;
					R++;
				}
			}
		}
		cout << "# " << R << " potential relations found." << endl << flush;
	   
		// compute and factor resultants as much as possible, leave large gcd computation till later.
		mpz_ui_pow_ui(lpb, 2, lpbits);
		int BASE = 16;
		stack<mpz_t*> QN; stack<int> Q; int algarr[3]; mpz_t* N;
		start = clock();
		R = 0;
		if (verbose) cout << "Starting cofactorizaion..." << endl << flush;
		for (int i = 0; i < rel.size()-1; i++)
		{
			if (rel[i] == rel[i+1] && rel[i] != 0) {
				int x = rel[i] % B1;
				int y = ((rel[i] >> B1bits) % B2x2) - B2;
				int z = ((rel[i] >> B1xB2x2bits) % B3x2) - B3;
				int t = (rel[i] >> B1xB2x2xB3x2bits) - B4;
				if (x != 0 || y != 0 || z != 0) {
					// compute [a,b,c]
					int a = L[0]*x+L[1]*y+L[2]*z+L[3]*t;
					int b = L[4]*x+L[5]*y+L[6]*z+L[7]*t;
					int c = L[8]*x+L[9]*y+L[10]*z+L[11]*t;
					int d = L[12]*x+L[13]*y+L[14]*z+L[15]*t;
	
					//int64_t D = b*(int64_t)b - 4*a*(int64_t)c;
					//int64_t D = b*(int64_t)b*c*(int64_t)c - 4*b*(int64_t)b*b*d
					//    -4*(int64_t)a*c*(int64_t)c*c + 18*(int64_t)a*b*c*d
					//	-27*(int64_t)a*a*(int64_t)d*d;
					//if (floor(sqrt(D)+0.5)*floor(sqrt(D)+0.5) == D) {
					//	continue;	// a+b*x+c*x^2 is not irreducible over Z
					//}

					int content = gcd(a, b); content = gcd(content, c);
					content = gcd(content, d);
					a = a/content; b = b/content; c = c/content; d = d/content;

					//cout << "[a, b, c] = [" << a << ", " << b << ", " << c << "]" << endl << flush;
					mpz_poly_setcoeff_si(A, 0, a);
					mpz_poly_setcoeff_si(A, 1, b);
					mpz_poly_setcoeff_si(A, 2, c);
					mpz_poly_setcoeff_si(A, 2, d);
					mpz_poly_resultant(N0, f0, A);
					mpz_poly_resultant(N1, f1, A);
					mpz_abs(N0, N0);
					mpz_abs(N1, N1);
					//cout << mpz_get_str(NULL, 10, N0) << endl << flush;
					//cout << mpz_get_str(NULL, 10, N1) << endl << flush;
					string str = to_string(a) + "," + to_string(b) + "," + to_string(c) + 
						to_string(d) + ":";
					
					// trial division on side 0
					int p = primes[0]; int k = 0; 
					while (p < sieveP0[k0-1]) {
						int valp = 0;
						while (mpz_fdiv_ui(N0, p) == 0) {
							mpz_divexact_ui(N0, N0, p);
							valp++;
							stream.str("");
							stream << hex << p;
							str += stream.str() + ",";
						}
						if (p < 1000) {
							p = primes[++k];
							if (p > 1000) {
								k = 0;
								while (sieveP0[k] < 1000) k++;
							}
						}
						else {
							p = sieveP0[++k];
						}
					}
					if (mpz_fdiv_ui(N0, q) == 0 && qside == 0) {
						mpz_divexact_ui(N0, N0, q);
						stream.str("");
						stream << hex << q;
						str += stream.str();
					}
					bool isrel = true;
					bool cofactor = true;
					if (mpz_cmp_ui(N0, 1) == 0) { cofactor = false; }
					str += (qside == 0 && cofactor ? "," : "");
					// cofactorization on side 0
					int n = 0; while (!Q.empty()) Q.pop(); while (!QN.empty()) QN.pop();
					if (cofactor) {
						if (mpz_probab_prime_p(N0, 30) == 0) {  // cofactor definitely composite
							
							QN.push(&N0); Q.push(2); Q.push(1); Q.push(0); Q.push(3);
							while (!QN.empty()) {
								mpz_t* N = QN.top(); QN.pop();
								int l = Q.top(); Q.pop();
								int j = 0;
								bool factored = false;
								while (!factored) {
									int alg = Q.top(); Q.pop(); j++;
									switch (alg) {
										case 0: factored = PollardPm1(*N, S, factor);
												break;	
										case 1: factored = EECM(*N, S, factor, 25921, 83521, 19, 9537, 2737);
												break;	
										case 2: factored = EECM(*N, S, factor, 1681, 707281, 3, 19642, 19803);
												break;	
									}
									if ( !factored ) {
										if ( j >= l ) { isrel = false; break; }
									}
									else {
										mpz_set(p1, factor);
										mpz_divexact(p2, *N, factor);
										if (mpz_cmpabs(p1, p2) > 0) {
											 mpz_set(tmp1, p1); mpz_set(p1, p2); mpz_set(p2, tmp1);	// sort
										}
										// save remaining algs to array
										int lnext = l - j; int lt = lnext;
										while (lt--) { algarr[lt] = Q.top(); Q.pop(); }
										lt = lnext; if (lt) { while (lt--) Q.push(algarr[lnext-1-lt]); Q.push(lnext); }
										if (mpz_probab_prime_p(p1, 30)) {
											if (mpz_cmpabs(p1, lpb) > 0) { isrel = false; break; }
											else { mpz_get_str(str2, BASE, p1); str += str2; str += ","; }
										}
										else {
											if (!lnext) { isrel = false; break; }
											mpz_set(pi[n], p1);
											QN.push(&pi[n++]);
											lt = lnext; if (lt) { while (lt--) Q.push(algarr[lnext-1-lt]); Q.push(lnext); }
										}
										if (mpz_probab_prime_p(p2, 30)) {
											if (mpz_cmpabs(p2, lpb) > 0) { isrel = false; break; }
											else { mpz_get_str(str2, BASE, p2); str += str2; str += QN.empty() ? "" : ","; }
										}
										else {
											if (!lnext) { isrel = false; break; }
											mpz_set(pi[n], p2);
											QN.push(&pi[n++]);
											lt = lnext; if (lt) { while (lt--) Q.push(algarr[lnext-1-lt]); Q.push(lnext); }
										}
									}
								}
								if (!isrel) break;
							}
						}
						else {	// cofactor prime but is it < lpb?
							if (mpz_cmpabs(N0, lpb) > 0) isrel = false;
							else { mpz_get_str(str2, BASE, N0); str += str2; }
						}
					}
					
					str += ":";

					// trial division on side 1
					if (isrel) {
						p = primes[0]; k = 0;
						while (p < sieveP1[k1-1]) {
							int valp = 0;
							while (mpz_fdiv_ui(N1, p) == 0) {
								mpz_divexact_ui(N1, N1, p);
								valp++;
								stream.str("");
								stream << hex << p;
								str += stream.str() + ",";
							}
							if (p < 1000) {
								p = primes[++k];
								if (p > 1000) {
									k = 0;
									while (sieveP1[k] < 1000) k++;
								}
							}
							else {
								p = sieveP1[++k];
							}
						}
						if (mpz_fdiv_ui(N1, q) == 0 && qside == 1)  {
							mpz_divexact_ui(N1, N1, q);
							stream.str("");
							stream << hex << q;
							str += stream.str();
						}
						bool cofactor = true;
						if (mpz_cmp_ui(N1, 1) == 0) { cofactor = false; }
						str += (qside == 1 && cofactor ? "," : "");
						// cofactorization on side 1
						n = 0; while (!Q.empty()) Q.pop(); while (!QN.empty()) QN.pop();
						if (cofactor) {
							if (mpz_probab_prime_p(N1, 30) == 0) {  // cofactor definitely composite
								
								QN.push(&N1); Q.push(2); Q.push(1); Q.push(0); Q.push(3);
								while (!QN.empty()) {
									mpz_t* N = QN.top(); QN.pop();
									int l = Q.top(); Q.pop();
									int j = 0;
									bool factored = false;
									while (!factored) {
										int alg = Q.top(); Q.pop(); j++;
										switch (alg) {
											case 0: factored = PollardPm1(*N, S, factor);
													break;	
											case 1: factored = EECM(*N, S, factor, 25921, 83521, 19, 9537, 2737);
													break;	
											case 2: factored = EECM(*N, S, factor, 1681, 707281, 3, 19642, 19803);
													break;	
										}
										if ( !factored ) {
											if ( j >= l ) { isrel = false; break; }
										}
										else {
											mpz_set(p1, factor);
											mpz_divexact(p2, *N, factor);
											if (mpz_cmpabs(p1, p2) > 0) {
												 mpz_set(tmp1, p1); mpz_set(p1, p2); mpz_set(p2, tmp1);	// sort
											}
											// save remaining algs to array
											int lnext = l - j; int lt = lnext;
											while (lt--) { algarr[lt] = Q.top(); Q.pop(); }
											lt = lnext; if (lt) { while (lt--) Q.push(algarr[lnext-1-lt]); Q.push(lnext); }
											if (mpz_probab_prime_p(p1, 30)) {
												if (mpz_cmpabs(p1, lpb) > 0) { isrel = false; break; }
												else { mpz_get_str(str2, BASE, p1); str += str2; str += ","; }
											}
											else {
												if (!lnext) { isrel = false; break; }
												mpz_set(pi[n], p1);
												QN.push(&pi[n++]);
												lt = lnext; if (lt) { while (lt--) Q.push(algarr[lnext-1-lt]); Q.push(lnext); }
											}
											if (mpz_probab_prime_p(p2, 30)) {
												if (mpz_cmpabs(p2, lpb) > 0) { isrel = false; break; }
												else { mpz_get_str(str2, BASE, p2); str += str2; str += QN.empty() ? "" : ","; }
											}
											else {
												if (!lnext) { isrel = false; break; }
												mpz_set(pi[n], p2);
												QN.push(&pi[n++]);
												lt = lnext; if (lt) { while (lt--) Q.push(algarr[lnext-1-lt]); Q.push(lnext); }
											}
										}
									}
									if (!isrel) break;
								}
							}
							else {	// cofactor prime but is it < lpb?
								if (mpz_cmpabs(N1, lpb) > 0) isrel = false;
								else { mpz_get_str(str2, BASE, N1); str += str2; }
							}
						}

						if (isrel) { cout << str << endl << flush; R++; }
					}
				}
			}
		}
		timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
		cout << "# Finished! Cofactorization took " << timetaken << "s" << endl << flush;
		cout << "# " << R << " actual relations found." << endl << flush;
		//cout << "lpb = " << mpz_get_str(NULL, 10, lpb) << endl << flush;
	}

	free(str2);
	mpz_clear(res);
	mpz_poly_clear(Aq0);
	mpz_poly_bivariate_clear(Aq);
	mpz_poly_clear(Fqh_x);
	mpz_clear(S);
	mpz_clear(tmp1); mpz_clear(p2); mpz_clear(p1);
	mpz_clear(factor);
	mpz_clear(lpb);
    mpz_clear(N1); mpz_clear(N0);
    mpz_poly_clear(A); mpz_poly_clear(f1); mpz_poly_clear(f0);
	delete[] fhqt;
	delete[] r;
	for (int i = 0; i < 8; i++) mpz_clear(pi[i]); delete[] pi;
	mpz_clear(qmpz);
	delete[] H;
	//delete[] L;
	delete[] M;
	delete[] sievenum_S1modp;
	delete[] sieveP1;
	delete[] sieveS1;
	delete[] sievenum_S0modp;
	delete[] sieveP0;
	delete[] sieveS0;
	delete[] sievenum_s1modp;
	delete[] sievenum_s0modp;
	delete[] sievep1;
	delete[] sieves1;
	delete[] sievep0;
	delete[] sieves0;
	mpz_clear(r0);
	delete[] primes;
	delete[] sieve;
	mpz_clear(F1ij);
	mpz_poly_clear(F1i);
	mpz_poly_bivariate_clear(F1);
	mpz_clear(F0ij);
	mpz_poly_clear(F0i);
	mpz_poly_bivariate_clear(F0);
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


void histogram(keyval*M, uint8_t* H, int len)
{
	// clear H
	memset(H, 0, len * sizeof(uint8_t));
	// fill H
	for (int i = 0; i < len; i++) {
		H[M[i].id] += M[i].logp;
	}
}


void csort(keyval* M, keyval* L, int* H, int len)
{
	// copy original M
	memcpy(L, M, len * sizeof(keyval));
	// clear H
	memset(H, 0, len * sizeof(int));
	// fill H
	for (int i = 0; i < len; i++) H[M[i].id]++;
	// calculate starting index for each key
	int m = 0;
	cout << "calculating start indices..." << endl << flush;
	for (int i = 0; i < len; i++) {
		int oldH = H[i];
		H[i] = m;
		m += oldH;
	}
	cout << "writing output..." << endl << flush;
	// write output
	for (int i = 0; i < len; i++) {
		memcpy(&M[H[L[i].id]], &L[i], sizeof(keyval));
		H[L[i].id]++;
	}
}


inline void mat3x3prod(int64_t* L1, int64_t* L2, int64_t* L3)
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			L3[i*3 + j] = 0;
			for (int k = 0; k < 3; k++) {
				L3[i*3 + j] += L1[i*3 + k] * L2[k*3 + j];
			}
		}
	}
}


inline void mat4x4prod(int64_t* L1, int64_t* L2, int64_t* L3)
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			L3[i*4 + j] = 0;
			for (int k = 0; k < 4; k++) {
				L3[i*4 + j] += L1[i*4 + k] * L2[k*4 + j];
			}
		}
	}
}


bool lattice_sorter(keyval const& kv1, keyval const& kv2)
{
	return kv2.id < kv1.id;
}


int latsieve4d(mpz_poly_bivariate F, int64_t* h, int degh, int64_t* fh_t, int degfh_t,
		int64_t q, int l, int ll, int* allp, int nump, int* s, int* S, int* num_smodp,
		int* num_Smodp, keyval* M, int Mlen, int* B)
{
	int64_t L[16];
	int64_t L2[16];
	int64_t L3[16];
	int64_t* r = new int64_t[degh]();
	int64_t* R = new int64_t[degfh_t]();
	int numl = polrootsmod(h, degh, r, q);
	int numll = polrootsmod(fh_t, degfh_t, R, q);
	int B1bits = B[0];
	int B2bits = B[1];
	int B3bits = B[2];
	int B4bits = B[3];
	int B1 = 1<<B1bits;
	int B2 = 1<<B2bits;
	int B3 = 1<<B3bits;
	int B4 = 1<<B4bits;
	int B1max = B1 - 1;
	int B2max = B2 - 1;
	int B3max = B3 - 1;
	int B4max = B4 - 1;
	int B1xB2x2bits = B1bits + B2bits + 1;
	int B1xB2x2xB3x2bits = B1bits + B2bits + 1 + B3bits + 1;

	// print (q,r) ideal
	//cout << "special-q (" << q << "," << r[l] << ")" << endl;

	// compute modpinvq and modqinvp arrays
	int64_t* modpinvq = new int64_t[nump];
	int64_t* modqinvp = new int64_t[nump];
	for (int i = 0; i < nump; i++) {
		int64_t p = allp[i];
		if (p == q) continue;
		modpinvq[i] = modinv(q, p);
		modqinvp[i] = modinv(p, q);
	}

	L[0]  = q; L[1]  = -r[l]; L[2]  = -R[ll]; L[3]  = 0;
	L[4]  = 0; L[5]  = 1;     L[6]  = 0;      L[7]  = -R[ll];
	L[8]  = 0; L[9]  = 0;     L[10] = 1;      L[11] = 0;
	L[12] = 0; L[13] = 0;     L[14] = 0;      L[15] = 1;

	int64L2(L, 4);	// LLL reduce L, time log(q)^2

	// print special-q lattice
	//cout << "special-q lattice [" << L[0] << "," << L[1] << "," << L[2] << ";";
	//cout << L[3] << "," << L[4] << "," << L[5] << ";";
	//cout << L[6] << "," << L[7] << "," << L[8] << "]" << endl << flush;

	int64_t qLinv[16];
	// compute adjugate of L^-1
	qLinv[0]  = (L[15]*L[10] - L[14]*L[11])*L[5] + ((-L[15]*L[9] + L[13]*L[11])*L[6] + (L[14]*L[9] - L[13]*L[10])*L[7]);
	qLinv[1]  = (-L[15]*L[10] + L[14]*L[11])*L[1] + ((-L[14]*L[9] + L[13]*L[10])*L[3] + (L[2]*L[15]*L[9] - L[2]*L[13]*L[11]));
	qLinv[2]  = (L[15]*L[6] - L[14]*L[7])*L[1] + ((L[14]*L[5] - L[13]*L[6])*L[3] + (-L[2]*L[15]*L[5] + L[2]*L[13]*L[7]));
	qLinv[3]  = (-L[11]*L[6] + L[10]*L[7])*L[1] + ((-L[10]*L[5] + L[9]*L[6])*L[3] + (L[2]*L[11]*L[5] - L[2]*L[9]*L[7]));
	qLinv[4]  = (-L[15]*L[10] + L[14]*L[11])*L[4] + ((L[15]*L[8] - L[12]*L[11])*L[6] + (-L[14]*L[8] + L[12]*L[10])*L[7]);
	qLinv[5]  = (L[15]*L[10] - L[14]*L[11])*L[0] + ((L[14]*L[8] - L[12]*L[10])*L[3] + (-L[2]*L[15]*L[8] + L[2]*L[12]*L[11]));
	qLinv[6]  = (-L[15]*L[6] + L[14]*L[7])*L[0] + ((-L[14]*L[4] + L[12]*L[6])*L[3] + (L[2]*L[15]*L[4] - L[2]*L[12]*L[7]));
	qLinv[7]  = (L[11]*L[6] - L[10]*L[7])*L[0] + ((L[10]*L[4] - L[8]*L[6])*L[3] + (-L[2]*L[11]*L[4] + L[2]*L[8]*L[7]));
	qLinv[8]  = (L[15]*L[9] - L[13]*L[11])*L[4] + ((-L[15]*L[8] + L[12]*L[11])*L[5] + (L[13]*L[8] - L[12]*L[9])*L[7]);
	qLinv[9]  = (-L[15]*L[9] + L[13]*L[11])*L[0] + ((L[15]*L[8] - L[12]*L[11])*L[1] + (-L[13]*L[8] + L[12]*L[9])*L[3]);
	qLinv[10] = (L[15]*L[5] - L[13]*L[7])*L[0] + ((-L[15]*L[4] + L[12]*L[7])*L[1] + (L[13]*L[4] - L[12]*L[5])*L[3]);
	qLinv[11] = (-L[11]*L[5] + L[9]*L[7])*L[0] + ((L[11]*L[4] - L[8]*L[7])*L[1] + (-L[9]*L[4] + L[8]*L[5])*L[3]);
	qLinv[12] = (-L[14]*L[9] + L[13]*L[10])*L[4] + ((L[14]*L[8] - L[12]*L[10])*L[5] + (-L[13]*L[8] + L[12]*L[9])*L[6]);
	qLinv[13] = (L[14]*L[9] - L[13]*L[10])*L[0] + ((-L[14]*L[8] + L[12]*L[10])*L[1] + (L[2]*L[13]*L[8] - L[2]*L[12]*L[9]));
	qLinv[14] = (-L[14]*L[5] + L[13]*L[6])*L[0] + ((L[14]*L[4] - L[12]*L[6])*L[1] + (-L[2]*L[13]*L[4] + L[2]*L[12]*L[5]));
	qLinv[15] = (L[10]*L[5] - L[9]*L[6])*L[0] + ((-L[10]*L[4] + L[8]*L[6])*L[1] + (L[2]*L[9]*L[4] - L[2]*L[8]*L[5]));

	int i = 40; int m = 0;
	while (i < nump) {
		int64_t p = allp[i];
		if (p == q) { i++; continue; }
		//cout << p << "," << m << endl << flush;
		int64_t n = p*q;
		uint8_t logp = log2f(p);
		__int128 rl = mod(-r[l], q);
		__int128 Rll = mod(-R[ll], q);
		for (int k = 0; k < num_smodp[i]; k++) {
			__int128 sk = mod(-s[k+i*degh],p);
			int64_t t_ = q*( (sk * modpinvq[i]) % p ) + p*( (rl * modqinvp[i]) % q ); // CRT
			if (t_ >= n) t_ -= n;
			for (int K = 0; K < num_Smodp[i]; K++) {
				__int128 SK = mod(-S[K+i*degfh_t],p);
				int64_t T = q*( (SK * modpinvq[i]) % p ) + p*( (Rll * modqinvp[i]) % q ); // CRT
				if (T >= n) T -= n;
				L2[0]  = n; L2[1]  = t_; L2[2]  = T; L2[3]  = 0;
				L2[4]  = 0; L2[5]  = 1;  L2[6]  = 0; L2[7]  = T;
				L2[8]  = 0; L2[9]  = 0;  L2[10] = 1; L2[11] = 0;
				L2[12] = 0; L2[13] = 0;  L2[14] = 0; L2[15] = 1;
				int64L2(L2,4);	// LLL reduce L, time log(n)^2
				mat4x4prod(qLinv, L2, L3);	// L3 =  qLinv*L2
				for (int ii = 0; ii < 16; ii++) L3[ii] /= q;
				int64L2(L3,4);
				int u1 = L3[0]; int u2 = L3[4]; int u3 = L3[8];  int u4 = L3[12];
				int v1 = L3[1]; int v2 = L3[5]; int v3 = L3[9];  int v4 = L3[13];
				int w1 = L3[2]; int w2 = L3[6]; int w3 = L3[10]; int w4 = L3[14];
				int t1 = L3[3]; int t2 = L3[7]; int t3 = L3[11]; int t4 = L3[15];

				// compute normal (cross product)
				int nx = (-w4*v3 + w3*v4)*u2 + ( w4*v2 - w2*v4)*u3 + (-w3*v2 + w2*v3)*u4;
				int ny = ( w4*v3 - w3*v4)*u1 + (-w4*v1 + w1*v4)*u3 + ( w3*v1 - w1*v3)*u4;
				int nz = (-w4*v2 + w2*v4)*u1 + ( w4*v1 - w1*v4)*u2 + (-w2*v1 + w1*v2)*u4;
				int nt = ( w3*v2 - w2*v3)*u1 + (-w3*v1 + w1*v3)*u2 + ( w2*v1 - w1*v2)*u3;

				// enumerate lattice vectors (x,y,z,t) in box [0,B[x[-B,B[x[-B,B[x[-B,B[
				int x = 0; int y = 0; int z = 0; int t = 0;
				while (spaceintersectsbox4d(nx, ny, nz, nt, x-t1, y-t2, z-t3, t-t4, B1, B2, B3, B4)) {
					x -= t1; y -= t2; z -= t3; t -= t4;
				}
				int ts1 = x; int ts2 = y; int ts3 = z; int ts4 = t;
				while (spaceintersectsbox4d(nx, ny, nz, nt, x, y, z, t, B1, B2, B3, B4)) { 
					int j1, j2, j3, j4, jmin;
					// pin vector to start of row
					int w1B = w1 < 0 ? B1max : 0; int w2B = w2 < 0 ? B2max : -B2;
					int w3B = w3 < 0 ? B3max : -B3; int w4B = w4 < 0 ? B4max : -B4;
					j1 = w1 == 0 ? -1 : (x - w1B) / w1;
					j2 = w2 == 0 ? -1 : (y - w2B) / w2;
					j3 = w3 == 0 ? -1 : (z - w3B) / w3;
					j4 = w4 == 0 ? -1 : (t - w4B) / w4;
					jmin = minnonneg4d(j1, j2, j3, j4);
					x -= jmin*w1; y -= jmin*w2; z -= jmin*w3; t -= jmin*w4;
					int ws1 = x; int ws2 = y; int ws3 = z; int ws4 = t;
					bool inspace = true;
					while (inspace) {
						if (x >= 0 && x < B1 && y >= -B2 && y < B2
							  && z >= -B3 && z < B3 && t >= -B4 && t < B4) {
							int a = 0; int b = 0;
							getab4d(u1, u2, u3, u4, v1, v2, v3, v4, x, y, z, t, B1, B2, B3, B4, &a, &b);
							x += a*u1 + b*v1;
							y += a*u2 + b*v2;
							z += a*u3 + b*v3;
							t += a*u4 + b*v4;
							int s1 = x; int s2 = y; int s3 = z; int s4 = t;
							// move 'forward' (dir == -1) or 'backward' (dir == 1) in plane
							for (int dir = -1; dir <= 1; dir += 2) {
								bool inplane = true;
								while (inplane) {
									// pVin vector to start of row
									int u1B = u1 < 0 ? B1max : 0; int u2B = u2 < 0 ? B2max : -B2;
									int u3B = u3 < 0 ? B3max : -B3; int u4B = u4 < 0 ? B4max : -B4;
									j1 = u1 == 0 ? -1 : (x - u1B) / u1;
									j2 = u2 == 0 ? -1 : (y - u2B) / u2;
									j3 = u3 == 0 ? -1 : (z - u3B) / u3;
									j4 = u4 == 0 ? -1 : (t - u4B) / u4;
									jmin = minnonneg4d(j1, j2, j3, j4);
									x -= jmin*u1; y -= jmin*u2; z -= jmin*u3; t -= jmin*u4;
									if (x >= 0 && x < B1 && y >= -B2 && y < B2
										  && z >= -B3 && z < B3 && t >= -B4 && t < B4) {
										int u1B = u1 < 0 ? 0 : B1max; int u2B = u2 < 0 ? -B2 : B2max;
										int u3B = u3 < 0 ? -B3 : B3max; int u4B = u4 < 0 ? -B4 : B4max;
										j1 = u1 == 0 ? -1 : (u1B - x) / u1;
										j2 = u2 == 0 ? -1 : (u2B - y) / u2;
										j3 = u3 == 0 ? -1 : (u3B - z) / u3;
										j4 = u4 == 0 ? -1 : (u4B - t) / u4;
										jmin = minnonneg4d(j1, j2, j3, j4);
										for (int j = 0; j <= jmin; j++) {
											int id = x + ((y + B2) << B1bits) + ((z + B3) << B1xB2x2bits)
											  + ((t + B4) << B1xB2x2xB3x2bits);
											M[m++] = (keyval){ id, logp };
											if (id == 140329463) {
												cout << id << ": " << p << endl;
												cout << x << "," << y << "," << z << "," << t << endl;
											}
											x += u1; y += u2; z += u3; t += u4;
										}
										// move by '1-transition' vector
										if (dir == -1) { x += v1; y += v2; z += v3; t += v4; }
										else		   { x -= v1; y -= v2; z -= v3; t -= v4; }
									}
									else {
										inplane = false;	
									}
								}
								x = s1 - v1; y = s2 - v2; z = s3 - v3; t = s4 - v4;
							}
							// advance to next w grid plane
							ws1 += w1; ws2 += w2; ws3 += w3; ws4 += w4; x = ws1; y = ws2; z = ws3; t = ws4;
						}
						else {
							inspace = false;
						}
					}
					// advance to next t grid plane
					ts1 += t1; ts2 += t2; ts3 += t3; ts4 += t4; x = ts1; y = ts2; z = ts3; t = ts4;
				}
			}
		}
		// advance to next p
		i++;
	}

	// clear memory
	delete[] modqinvp;
	delete[] modpinvq;
	delete[] r;

	return m;
}

bool planeintersectsbox4d(int u1, int u2, int u3, int u4, int v1, int v2, int v3, int v4,
		int w1, int w2, int w3, int w4, int B1, int B2, int B3, int B4)
{
	int64_t A[16]; int64_t M[16];

	// vectors
	int64_t U1 = u1-w1; int64_t U2 = u2-w2; int64_t U3 = u3-w3; int64_t U4 = u4-w4;
	int64_t V1 = v1-w1; int64_t V2 = v2-w2; int64_t V3 = v3-w3; int64_t V4 = v4-w4;
	int64_t F1, F2, F3, F4, G1, G2, G3, G4;

	// 3 points determine a face in 4d, updated one at a time, giving a new face
	int64_t T[3][4] = { { 0, -B2, -B3, -B4 }, { B1, -B2, -B3, -B4},	{ B1, B2, -B3, -B4} };
	int pos = 2;	// position in cycle of 3 points (FIFO, mod 3)
	int xyzt = 2;	// component position mod 4

	// compute shortest normal vector (N1,N2,N3,N4) between plane (U,V,w) and 24 faces of 4d
	// boundary volume, then determine whether 3-space defined by (U x V x N, w) intersects
	// 4d boundary volume.
	int64_t N1 = 0; int64_t N2 = 0; int64_t N3 = 0; int64_t N4 = 0; int64_t dbest = 1l<<63;
	for (int i = 0; i < 24; i++) {
		// first construct i-th face of 4d boundary volume
		F1 = T[1][0]-T[0][0]; F2 = T[1][1]-T[0][1]; F3 = T[1][2]-T[0][2]; F4 = T[1][3]-T[0][3];
		G1 = T[2][0]-T[0][0]; G2 = T[2][1]-T[0][1]; G3 = T[2][2]-T[0][2]; G4 = T[2][3]-T[0][3];
		int64_t h1 = T[0][0]; int64_t h2 = T[0][1]; int64_t h3 = T[0][2]; int64_t h4 = T[0][3];
		// find shortest vector between planes (U,V,w) and (F,G,h) in 4d
		// d^2 = ( (w + x*U + y*V) - (h + r*F + s*G) )^2
		//  U^2*x + V*U*y - F*U*r - G*U*s + U*w - U*h = 0
		//  V*U*x + V^2*y - F*V*r - G*V*s + V*w - V*h = 0
		// -F*U*x - F*V*y + F^2*r + G*F*s - F*w + F*h = 0
		// -G*U*x - G*V*y + G*F*r + G^2*s - G*w + G*h = 0
		A[0] = U1*U1 + U2*U2 + U3*U3 + U4*U4;
		A[1] = V1*U1 + V2*U2 + V3*U3 + V4*U4;
		A[2] = -(F1*U1 + F2*U2 + F3*U3 + F4*U4);
		A[3] = -(G1*U1 + G2*U2 + G3*U3 + G4*U4);
		A[4] = A[1];
		A[5] = V1*V1 + V2*V2 + V3*V3 + V4*V4;
		A[6] = -(F1*V1 + F2*V2 + F3*V3 + F4*V4);
		A[7] = -(G1*V1 + G2*V2 + G3*V3 + G4*V4);
		A[8] = A[2]; A[9] = A[6];
		A[10] = F1*F1 + F2*F2 + F3*F3 + F4*F4;
		A[11] = G1*F1 + G2*F2 + G3*F3 + G4*F4;
		A[12] = A[3]; A[13] = A[7]; A[14] = A[11];
		A[15] = G1*G1 + G2*G2 + G3*G3 + G4*G4;
		int64_t D = fadlev4d(A, M);	// get adjugate M and determinant D of A
		// compute RHS vector
		int64_t R1 = U1*(h1-w1) + U2*(h2-w2) + U3*(h3-w3) + U4*(h4-w4);
		int64_t R2 = V1*(h1-w1) + V2*(h2-w2) + V3*(h3-w3) + V4*(h4-w4);
		int64_t R3 = F1*(w1-h1) + F2*(w2-h2) + F3*(w3-h3) + F4*(w4-h4);
		int64_t R4 = G1*(w1-h1) + G2*(w2-h2) + G3*(w3-h3) + G4*(w4-h4);
		if (D) {
			// solve linear system
			int64_t xx = M[0]*R1 + M[1]*R2 + M[2]*R3 + M[3]*R4; // actual x = xx/D
			int64_t yy = M[4]*R1 + M[5]*R2 + M[6]*R3 + M[7]*R4; // actual y = yy/D
			int64_t rr = M[8]*R1 + M[9]*R2 + M[10]*R3 + M[11]*R4; // actual r = rr/D
			int64_t ss = M[12]*R1 + M[13]*R2 + M[14]*R3 + M[15]*R4; // actual s = ss/D
			int64_t nn1 = D*w1 + xx*U1 + yy*V1; int64_t nn2 = D*w2 + xx*U2 + yy*V2;
			int64_t nn3 = D*w3 + xx*U3 + yy*V3; int64_t nn4 = D*w4 + xx*U4 + yy*V4;
			int64_t nn5 = D*h1 + rr*F1 + ss*G1; int64_t nn6 = D*h2 + rr*F2 + ss*G2;
			int64_t nn7 = D*h3 + rr*F3 + ss*G3; int64_t nn8 = D*h4 + rr*F4 + ss*G4;
			int64_t N1new = (nn1 - nn5)/D; int64_t N2new = (nn2 - nn6)/D;
			int64_t N3new = (nn3 - nn7)/D; int64_t N4new = (nn4 - nn8)/D;
			int64_t d = N1new*N1new + N2new*N2new + N3new*N3new + N4new*N4new;
			if (d < dbest) {
				N1 = N1new; N2 = N2new; N3 = N3new; N4 = N4new;
			}
		}
		// update list of points
		int old = pos;
		pos = (pos + 1) % 3;
		T[pos][0] = T[old][0]; T[pos][1] = T[old][1];
		T[pos][2] = T[old][2]; T[pos][3] = T[old][3];
		if (xyzt == 0) T[pos][xyzt] = B1 - T[pos][xyzt];
		else T[pos][xyzt] = -T[pos][xyzt];
		xyzt = (xyzt + 1) % 4;
	}
	// compute normal (cross product) W = U x V x N
	int W1 = (-N4*V3 + N3*V4)*U2 + ( N4*V2 - N2*V4)*U3 + (-N3*V2 + N2*V3)*U4;
	int W2 = ( N4*V3 - N3*V4)*U1 + (-N4*V1 + N1*V4)*U3 + ( N3*V1 - N1*V3)*U4;
	int W3 = (-N4*V2 + N2*V4)*U1 + ( N4*V1 - N1*V4)*U2 + (-N2*V1 + N1*V2)*U4;
	int W4 = ( N3*V2 - N2*V3)*U1 + (-N3*V1 + N1*V3)*U2 + ( N2*V1 - N1*V2)*U3;

	// determine whether 3-space defined by (W, w) intersects 4d boundary volume.
	return spaceintersectsbox4d(W1, W2, W3, W4, w1, w2, w3, w4, B1, B2, B3, B4);
}


/* This was a nice idea but is much more complicated than necessary.
   There exists a miraculously simple algorithm to traverse all the faces
   of a tesseract.  That is what we actually use.
bool planeintersectsbox4d(int u1, int u2, int u3, int u4, int v1, int v2, int v3, int v4,
		int w1, int w2, int w3, int w4, int B1, int B2, int B3, int B4)
{
	int64_t A[16]; int64_t M[16];

	// vectors
	U1 = u1-w1; U2 = u2-w2; U3 = u3-w3; U4 = u4-w4;
	V1 = v1-w1; V2 = v2-w2; V3 = v3-w3; V4 = v4-w4;

	// boundary point arrays.  We consider 6 groups of 4 faces, each group requires 6 points
	int64_t T1[8] = {   0,  B1,  B1,  B1,  0,    0,  };
	int64_t T2[8] = { -B2, -B2,  B2,  B2,  B2, -B2,  };
	int64_t T3[8] = { -B3, -B3, -B3, -B3, -B3, -B3,  };
	int64_t T4[8] = { -B4, -B4, -B4,  B4,  B4,  B4,  };

	// compute shortest normal vector (N1,N2,N3,N4) between plane (U,V,w) and 24 faces of 4d
	// boundary volume, then determine whether 3-space defined by (U x V x N, w) intersects
	// 4d boundary volume.
	N1 = 0; N2 = 0; N3 = 0; N4 = 0; int64_t dbest = 1<<63l;
	for (i = 0; i < 6; i++) {
		for (j = 0; j < 4; j++) {
			// first construct (i*j)-th face of 4d boundary volume
			F1 = T1[j+1]-T1[j]; F2 = T2[j+1]-T2[j]; F3 = T3[j+1]-T3[j]; F4 = T4[j+1]-T4[j];
			G1 = T1[j+2]-T1[j]; G2 = T2[j+2]-T2[j]; G3 = T3[j+2]-T3[j]; G4 = T4[j+2]-T4[j];
			h1 = T1[j]; h2 = T2[j]; h3 = T3[j]; h4 = T4[j];
			// find shortest vector between planes (U,V,w) and (F,G,h) in 4d
			// d^2 = ( (w + x*U + y*V) - (h + r*F + s*G) )^2
			//  U^2*x + V*U*y - F*U*r - G*U*s + U*w - h*U = 0
			//  V*U*x + V^2*y - F*V*r - G*V*s + V*w - h*V = 0
			// -F*U*x - F*V*y + F^2*r + G*F*s - F*w + F*h = 0
			// -G*U*x - G*V*y + G*F*r + G^2*s - G*w + G*h = 0
			A[0] = U1*U1 + U2*U2 + U3*U4 + U4*U4;
			A[1] = V1*U1 + V2*U2 + V3*U3 + V4*U4;
			A[2] = -(F1*U1 + F2*U2 + F3*U3 + F4*U4);
			A[3] = -(G1*U1 + G2*U2 + G3*U3 + G4*U4);
			A[4] = A[1];
			A[5] = V1*V1 + V2*V2 + V3*V4 + V4*V4;
			A[6] = -(F1*V1 + F2*V2 + F3*V3 + F4*V4);
			A[7] = -(G1*V1 + G2*V2 + G3*V3 + G4*V4);
			A[8] = A[2]; A[9] = A[6];
			A[10] = F1*F1 + F2*F2 + F3*F4 + F4*F4;
			A[11] = G1*F1 + G2*F2 + G3*F3 + G4*F4;
			A[12] = A[3]; A[13] = A[7]; A[14] = A[11];
			A[15] = G1*G1 + G2*G2 + G3*G3 + G4*G4;
			int64_t D = fadlev4d(A, M);	// get adjugate M and determinant D
			// compute RHS vector
			int64_t R1 = U1*(h1-w1) + U2*(h2-w2) + U3*(h3-w3) + U4*(h4-w4);
			int64_t R2 = V1*(h1-w1) + V2*(h2-w2) + V3*(h3-w3) + V4*(h4-w4);
			int64_t R3 = F1*(w1-h1) + F2*(w2-h2) + F3*(w3-h3) + F4*(w4-h4);
			int64_t R4 = G1*(w1-h1) + G2*(w2-h2) + G3*(w3-h3) + G4*(w4-h4);
			if (D) {
				// solve linear system
				int64_t xx = M[0]*R1 + M[1]*R2 + M[2]*R3 + M[3]*R4; // actual x = xx/D
				int64_t yy = M[4]*R1 + M[5]*R2 + M[6]*R3 + M[7]*R4; // actual y = yy/D
				int64_t rr = M[8]*R1 + M[9]*R2 + M[10]*R3 + M[11]*R4; // actual r = rr/D
				int64_t ss = M[12]*R1 + M[13]*R2 + M[14]*R3 + M[15]*R4; // actual s = ss/D
				int64_t nn1 = D*w1 + xx*U1 + yy*V1; int64_t nn2 = D*w2 + xx*U2 + yy*V2;
				int64_t nn3 = D*w3 + xx*U3 + yy*V3; int64_t nn4 = D*w4 + xx*U4 + yy*V4;
				int64_t nn5 = D*h1 + rr*F1 + ss*G1; int64_t nn6 = D*h2 + rr*F2 + ss*G2;
				int64_t nn7 = D*h3 + rr*F3 + ss*G3; int64_t nn8 = D*h4 + rr*F4 + ss*G4;
				int64_t N1new = (nn1 - nn5)/D; int64_t N2new = (nn2 - nn6)/D;
				int64_t N3new = (nn3 - nn7)/D; int64_t N3new = (nn4 - nn8)/D;
				int64_t d = N1new*N1new + N2new*N2new + N3new*N3new + N4new*N4new;
				if (d < dbest) {
					N1 = N1new; N2 = N2new; N3 = N3new; N4 = N4new;
				}
			}
		}
	}
	// compute normal (cross product) W = U x V x N
	int W1 = (-N4*V3 + N3*V4)*U2 + ( N4*V2 - N2*V4)*U3 + (-N3*V2 + N2*V3)*U4;
	int W2 = ( N4*V3 - N3*V4)*U1 + (-N4*V1 + N1*V4)*U3 + ( N3*V1 - N1*V3)*U4;
	int W3 = (-N4*V2 + N2*V4)*U1 + ( N4*V1 - N1*V4)*U2 + (-N2*V1 + N1*V2)*U4;
	int W4 = ( N3*V2 - N2*V3)*U1 + (-N3*V1 + N1*V3)*U2 + ( N2*V1 - N1*V2)*U3;

	// determine whether 3-space defined by (W, w) intersects 4d boundary volume.
	return spaceintersectsbox4d(W1, W2, W3, W4, w1, w2, w3, w4, B1, B2, B3, B4);
}
*/

inline int floordiv(int a, int b)
{
    int d = a / b;
    return d * b == a ? d : d - ((a < 0) ^ (b < 0));
}


// Integer programming to get (a,b) such that (x,y,z) + a*(u1,u2,u3) + b*(v1,v2,v3) within [0,B[x[-B,B[x[-B,B[
inline void getab(int u1, int u2, int u3, int v1, int v2, int v3, int x, int y, int z, 
				int B1, int B2, int B3, int* a, int* b)
{
	int U[6] = { u1, u2, u3, -u1, -u2, -u3 };
	int V[6] = { v1, v2, v3, -v1, -v2, -v3 };
	int C[6] = { B1-x-1, B2-y-1, B3-z-1, x, B2+y, B3+z };
	*a = B1; *b = B1;

	int L = abs(nonzerolcm(u1, u2, u3));

	for (int i = 0; i < 6; i++) {
		int t = abs(U[i]);
		if (t != 0) {
			V[i] *= L / t;
			C[i] *= L / t;
		}
	}
	
	for (int i = 0; i < 6; i++) {
		if (U[i] == 0) {
			if (V[i] > 0) {
				int bnew = floordiv(C[i], V[i]);
				if (bnew < *b) *b = bnew;
			}
		}
		else if (U[i] < 0) {
			for (int j = 0; j < 6; j++) {
				if (U[j] > 0 && abs(i-j) != 3) {
					int D = V[i] + V[j];
					if (D > 0) {
						int bnew = floordiv(C[i] + C[j], D);
						if (bnew < *b) *b = bnew;
					}
				}
			}
		}
	}

	for (int i = 0; i < 6; i++) {
		if (U[i] > 0) {
			int anew = floordiv(C[i] - V[i] * (*b), L);
			if (anew < *a) *a = anew;
		}
	}
}


// Integer programming to get (a,b) such that (x,y,z,t) + a*(u1,u2,u3,u4) + b*(v1,v2,3v,v4)
// within [0,B[x[-B,B[x[-B,B[x[-B,B[
inline void getab4d(int u1, int u2, int u3, int u4, int v1, int v2, int v3, int v4,
				int x, int y, int z, int t, int B1, int B2, int B3, int B4, int* a, int* b)
{
	int U[8] = { u1, u2, u3, u4, -u1, -u2, -u3, -u4 };
	int V[8] = { v1, v2, v3, v4, -v1, -v2, -v3, -v4 };
	int C[8] = { B1-x-1, B2-y-1, B3-z-1, B4-t-1, x, B2+y, B3+z, B4+t };
	*a = B1; *b = B1;

	int L = abs(nonzerolcm4d(u1, u2, u3, u4));

	for (int i = 0; i < 8; i++) {
		int s = abs(U[i]);
		if (s != 0) {
			V[i] *= L / s;
			C[i] *= L / s;
		}
	}
	
	for (int i = 0; i < 8; i++) {
		if (U[i] == 0) {
			if (V[i] > 0) {
				int bnew = floordiv(C[i], V[i]);
				if (bnew < *b) *b = bnew;
			}
		}
		else if (U[i] < 0) {
			for (int j = 0; j < 8; j++) {
				if (U[j] > 0 && abs(i-j) != 4) {
					int D = V[i] + V[j];
					if (D > 0) {
						int bnew = floordiv(C[i] + C[j], D);
						if (bnew < *b) *b = bnew;
					}
				}
			}
		}
	}

	for (int i = 0; i < 8; i++) {
		if (U[i] > 0) {
			int anew = floordiv(C[i] - V[i] * (*b), L);
			if (anew < *a) *a = anew;
		}
	}
}


inline int nonzerolcm(int u1, int u2, int u3)
{
	if (u1 == 0) u1 = 1;
	if (u2 == 0) u2 = 1;
	if (u3 == 0) u3 = 1;
	int g1 = gcd(u1, u2);
	int g2 = gcd(u2, u3);
	int g3 = gcd(g1, g2);
	return (u1 * u2 * u3) / g3;
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


inline int nonzerolcm4d(int u1, int u2, int u3, int u4)
{
	if (u1 == 0) u1 = 1;
	if (u2 == 0) u2 = 1;
	if (u3 == 0) u3 = 1;
	if (u4 == 0) u4 = 1;
	int g1 = gcd(u1, u2);
	int g2 = gcd(u2, u3);
	int g3 = gcd(g1, g2);
	int g4 = gcd(g3, u4);
	return (u1 * u2 * u3 * u4) / g4;
}


// determine if 3-space with normal (nx, ny, nz, nt) containing point (x,y,z,t) intersects
// box [0,B[x[-B,B[x[-B,B[x[-B,B[
inline bool spaceintersectsbox4d(int nx, int ny, int nz, int nt, int x, int y, int z, int t,
							int B1, int B2, int B3, int B4)
{
	int d = nx*x + ny*y + nz*z + nt*t;

	/* nxB0 = 0 */     int nxB1 = nx*(B1-1);
	int nyB0 = -ny*B2; int nyB1 = ny*(B2-1);
	int nzB0 = -nz*B3; int nzB1 = nz*(B3-1);
	int ntB0 = -nt*B4; int ntB1 = nt*(B4-1);

	int s0 = ( /*nxB0 +*/ nyB0 + nzB0 + ntB0 - d ) > 0;
	int s1 = ( nxB1   +   nyB0 + nzB0 + ntB0 - d ) > 0;
	int s2 = ( /*nxB0 +*/ nyB1 + nzB0 + ntB0 - d ) > 0;
	int s3 = ( nxB1   +   nyB1 + nzB0 + ntB0 - d ) > 0;
	int s4 = ( /*nxB0 +*/ nyB0 + nzB1 + ntB0 - d ) > 0;
	int s5 = ( nxB1   +   nyB0 + nzB1 + ntB0 - d ) > 0;
	int s6 = ( /*nxB0 +*/ nyB1 + nzB1 + ntB0 - d ) > 0;
	int s7 = ( nxB1   +   nyB1 + nzB1 + ntB0 - d ) > 0;
	int s8 = ( /*nxB0 +*/ nyB0 + nzB0 + ntB1 - d ) > 0;
	int s9 = ( nxB1   +   nyB0 + nzB0 + ntB1 - d ) > 0;
	int sa = ( /*nxB0 +*/ nyB1 + nzB0 + ntB1 - d ) > 0;
	int sb = ( nxB1   +   nyB1 + nzB0 + ntB1 - d ) > 0;
	int sc = ( /*nxB0 +*/ nyB0 + nzB1 + ntB1 - d ) > 0;
	int sd = ( nxB1   +   nyB0 + nzB1 + ntB1 - d ) > 0;
	int se = ( /*nxB0 +*/ nyB1 + nzB1 + ntB1 - d ) > 0;
	int sf = ( nxB1   +   nyB1 + nzB1 + ntB1 - d ) > 0;

	int mask = s0 + (s1<<1) + (s2<<2) + (s3<<3) + (s4<<4) + (s5<<5) + (s6<<6) + (s7<<7) +
	           (s8<<8) + (s9<<9) + (sa<<10) + (sb<<11) +
			   (sc<<12) + (sd<<13) + (se<<14) + (sf<<15);

	return ( (mask != 0) && (mask != 65535) );
}


/*
inline void getab(int u1, int u2, int u3, int v1, int v2, int v3, int x, int y, int z, int B, int* a, int* b)
{
	int n1 = u3 * (u3 < 0 ? -B-x :  B-x) - u1 * (u1 < 0 ?   -z :  B-z);
	int n2 = u1 * (u1 < 0 ? -B-y :  B-y) - u2 * (u2 < 0 ? -B-x :  B-x);
	int n3 = u2 * (u2 < 0 ?   -z :  B-z) - u3 * (u3 < 0 ? -B-y :  B-y);
	int m1 = u3 * (u3 < 0 ?  B-x : -B-x) - u1 * (u1 < 0 ?  B-z :   -z);
	int m2 = u1 * (u1 < 0 ?  B-y : -B-y) - u2 * (u2 < 0 ?  B-x : -B-x);
	int m3 = u2 * (u2 < 0 ?  B-z :   -z) - u3 * (u3 < 0 ?  B-y : -B-y);
	int d1 = u3*v1 - u1*v3;
	int d2 = u1*v2 - u2*v1;
	int d3 = u2*v3 - u3*v2;
	int b1 = d1 == 0 ? B : d1 < 0 ? n1/d1 : m1/d1;
	int b2 = d2 == 0 ? B : d2 < 0 ? n2/d2 : m2/d2;
	int b3 = d3 == 0 ? B : d3 < 0 ? n3/d3 : m3/d3;
	*b = minabs(b1, b2, b3);
	int a1 = u1 == 0 ? B : u1 < 0 ? ( -B - x - v1 * (*b) ) / u1 : ( B - x - v1 * (*b) ) / u1;
	int a2 = u2 == 0 ? B : u2 < 0 ? ( -B - y - v2 * (*b) ) / u2 : ( B - y - v2 * (*b) ) / u2;
	int a3 = u3 == 0 ? B : u3 < 0 ? (  0 - z - v3 * (*b) ) / u3 : ( B - z - v3 * (*b) ) / u3;
	*a = minabs(a1, a2, a3);
}
*/


// determine if plane with normal (nx, ny, nz) containing point (x,y,z) intersects box [0,B[x[-B,B[x[-B,B[
inline bool planeintersectsbox(int nx, int ny, int nz, int x, int y, int z, int B1, int B2, int B3)
{
	int d = nx*x + ny*y + nz*z;

	/* nxB0 = 0 */     int nxB1 = nx*(B1-1);
	int nyB0 = -ny*B2; int nyB1 = ny*(B2-1);
	int nzB0 = -nz*B3; int nzB1 = nz*(B3-1);

	int s0 = ( /*nxB0 +*/ nyB0 + nzB0 - d ) > 0;
	int s1 = ( nxB1   +   nyB0 + nzB0 - d ) > 0;
	int s2 = ( /*nxB0 +*/ nyB1 + nzB0 - d ) > 0;
	int s3 = ( nxB1   +   nyB1 + nzB0 - d ) > 0;
	int s4 = ( /*nxB0 +*/ nyB0 + nzB1 - d ) > 0;
	int s5 = ( nxB1   +   nyB0 + nzB1 - d ) > 0;
	int s6 = ( /*nxB0 +*/ nyB1 + nzB1 - d ) > 0;
	int s7 = ( nxB1   +   nyB1 + nzB1 - d ) > 0;

	int mask = s0 + (s1<<1) + (s2<<2) + (s3<<3) + (s4<<4) + (s5<<5) + (s6<<6) + (s7<<7);

	return ( (mask != 0) && (mask != 255) );
}


inline int minabs(int u, int v, int w)
{
	int m = u;
	if (abs(v) < abs(m)) m = v;
	if (abs(w) < abs(m)) m = w;
	return m;
}


inline int minnonneg(int u, int v, int w)
{
	int m = u;
	if (m < 0) m = v; if (m < 0) m = w;
	if (v >= 0 && v < m) m = v;
	if (w >= 0 && w < m) m = w;
	return m;
}


inline int minnonneg4d(int u, int v, int w, int t)
{
	int m = u;
	if (m < 0) m = v; if (m < 0) m = w; if (m < 0) m = t;
	if (v >= 0 && v < m) m = v;
	if (w >= 0 && w < m) m = w;
	if (t >= 0 && t < m) m = t;
	return m;
}


inline int maxabs(int u, int v, int w)
{
	int m = u;
	if (abs(v) > abs(m)) m = v;
	if (abs(w) > abs(m)) m = w;
	return m;
}


inline int min(int u, int v, int w)
{
	int m = u;
	if (v < m) m = v;
	if (w < m) m = w;
	return m;
}


inline int max(int u, int v, int w)
{
	int m = u;
	if (v > m) m = v;
	if (w > m) m = w;
	return m;
}


inline int64_t modinv(int64_t x, int64_t m)
{
    int64_t m0 = m, t, q;
    int64_t y0 = 0, y = 1;
    if (m == 1) return 1;
    while (x > 1) {
        q = x / m;
        t = m, m = x % m, x = t;
        t = y0, y0 = y - q * y0, y = t;
    }
    if (y < 0) y += m0;
    return y;
}

void GetlcmScalar(int B, mpz_t S, int* primes, int nump)
{
    mpz_t* tree = new mpz_t[nump];
    
    // Construct product tree
    int n = 0;
    mpz_t pe; mpz_init(pe); mpz_t pe1; mpz_init(pe1);
	int p = 2;
	while (p < B) {
		// Set tree[n] = p^e, such that p^e < B < p^(e+1)
		mpz_set_ui(pe, p);
		mpz_mul_ui(pe1, pe, p);
		while (mpz_cmp_ui(pe1, B) < 0) { mpz_set(pe, pe1); mpz_mul_ui(pe1, pe, p); }
		mpz_init(tree[n]);
		mpz_set(tree[n], pe);
		n++;
		p = primes[n];
    }
    mpz_clear(pe); mpz_clear(pe1);
    
    // Coalesce product tree
    uint64_t treepos = n - 1;
    while (treepos > 0) {
        for (int i = 0; i <= treepos; i += 2) {
            if(i < treepos)
				mpz_lcm(tree[i/2], tree[i],tree[i + 1]);
            else
				mpz_set(tree[i/2], tree[i]);
        }
        for (int i = (treepos >> 1); i < treepos - 1; i++) mpz_set_ui(tree[i + 1], 1);
        treepos = treepos >> 1;
    }
    // tree[0] is the lcm of all primes with powers bounded by B
    mpz_set(S, tree[0]);
        
    for (int i = 0; i < n; i++) mpz_clear(tree[i]);
    delete[] tree;
}


bool PollardPm1(mpz_t N, mpz_t S, mpz_t factor)
{
	int bitlen = mpz_sizeinbase(N, 2);
	if (0) { //bitlen < 64) {
		// convert N to __int128
		__int128 N128 = make_int128(mpz_getlimbn(N,0), mpz_getlimbn(N, 1));
		__int128 factor128 = 1;
		int64_t mulo = 0; int64_t muhi = 0;
		bool factored = PollardPm1_int128(N128, S, factor128, mulo, muhi);
		// convert factor128 to mpz_t
		mp_limb_t* factor_limbs = mpz_limbs_modify(factor, 2);
		factor_limbs[0] = factor128 & MASK64;
		factor_limbs[1] = factor128 >> 64;
		//if (factored) cout << mpz_get_str(NULL,10,N) << " = c * " << mpz_get_str(NULL,10,factor128) << " factored!" << endl; 
		return factored;
	}
	else if (0) { //bitlen < 128) {
		// convert N to __int128
		__int128 N128 = make_int128(mpz_getlimbn(N,0), mpz_getlimbn(N, 1));
		__int128 factor128 = 1;
		int64_t mulo = 0; int64_t muhi = 0;
		bool factored = PollardPm1_int128(N128, S, factor128, mulo, muhi);
		// convert factor128 to mpz_t
		mp_limb_t* factor_limbs = mpz_limbs_modify(factor, 2);
		factor_limbs[0] = factor128 & MASK64;
		factor_limbs[1] = factor128 >> 64;
		return factored;
	}
	else {
		return PollardPm1_mpz(N, S, factor);
	}
}


bool PollardPm1_mpz(mpz_t N, mpz_t S, mpz_t factor)
{
    int L = mpz_sizeinbase(S, 2); // exact for base = power of 2
	mpz_t g; mpz_init_set_ui(g, 2);
      
    // Scalar multiplication using square and multiply
    for (int i = 2; i <= L; i++) {
        // square
		mpz_mul(g, g, g);
		mpz_mod(g, g, N);
        if (mpz_tstbit(S, L - i) == 1) {
            // multiply
			mpz_mul_2exp(g, g, 1);
			if (mpz_cmpabs(g, N) >= 0) mpz_sub(g, g, N);
        }
    }
	// subtract 1
	mpz_sub_ui(g, g, 1);
	// compute gcd
	mpz_gcd(factor, N, g);
	bool result = mpz_cmpabs_ui(factor, 1) > 0 && mpz_cmpabs(factor, N) < 0;
	//if (result) cout << endl << endl << "\t\t\tP-1 worked!!!!" << endl << endl << flush;
	mpz_clear(g);
	return result;
}


inline __int128 gcd128(__int128 a, __int128 b)
{
	a = a < 0 ? -a : a;
	b = b < 0 ? -b : b;
	__int128 c;
	while (b != 0) {
		c = b;
		b = a % c;
		a = c;
	}
	return a;
}


bool PollardPm1_int128(__int128 N, mpz_t S, __int128 &factor, int64_t mulo, int64_t muhi)
{
    int L = mpz_sizeinbase(S, 2); // exact for base = power of 2
	__int128 g = 2;
      
    // Scalar multiplication using square and multiply
    for (int i = 2; i <= L; i++) {
        // square
		g = g * g;
		g = g % N;
        if (mpz_tstbit(S, L - i) == 1) {
            // multiply
			g = g * 2;
			if (g >= N) g -= N; 
        }
    }
	// subtract 1
	g = g - 1;
	// compute gcd
	factor = gcd128(N, g);
	bool result = (factor > 1) && (factor < N);
	//if (result) cout << endl << endl << "\t\t\tP-1 worked!!!!" << endl << endl << flush;
	return result;
}


bool EECM(mpz_t N, mpz_t S, mpz_t factor, int d, int a, int X0, int Y0, int Z0)
{
	int bitlen = mpz_sizeinbase(N, 2);
	if (0) { //bitlen < 44) {
		//cout << mpz_get_str(NULL,10,N) << endl; 
		// convert N to __int128
		__int128 N128 = make_int128(mpz_getlimbn(N,0), mpz_getlimbn(N, 1));
		__int128 factor128 = 1;
		int64_t mulo = 0; int64_t muhi = 0;
		bool factored = EECM_int128(N128, S, factor128, d, a, X0, Y0, Z0, mulo, muhi);
		// convert factor128 to mpz_t
		mp_limb_t* factor_limbs = mpz_limbs_modify(factor, 2);
		factor_limbs[0] = factor128 & MASK64;
		factor_limbs[1] = factor128 >> 64; //if (factored) cout << mpz_get_str(NULL,10,N) << " factored!" << endl; 
		return factored;
	}
	else if (0) { //bitlen < 128) {
		// convert N to __int128
		__int128 N128 = make_int128(mpz_getlimbn(N,0), mpz_getlimbn(N, 1));
		__int128 factor128 = 1;
		int64_t mulo = 0; int64_t muhi = 0;
		bool factored = EECM_int128(N128, S, factor128, d, a, X0, Y0, Z0, mulo, muhi);
		// convert factor128 to mpz_t
		mp_limb_t* factor_limbs = mpz_limbs_modify(factor, 2);
		factor_limbs[0] = factor128 & MASK64;
		factor_limbs[1] = factor128 >> 64;
		return factored;
	}
	else {
		return EECM_mpz(N, S, factor, d, a, X0, Y0, Z0);
	}
}


/* ScalarMultiplyEdwards
 * 
 * Multiply a point [X0:Y0:Z0] on a twisted edwards curve by a scalar multiple
 * d	d parameter of twisted Edwards curve
 * a	a parameter of twisted Edwards curve
 * X0,Y0,Z0	point on curve to multiply, in projective coordinates
 * N	we work modulo N
 * S	scalar multiple
 * L	length of S in bits");
*/
bool EECM_mpz(mpz_t N, mpz_t S, mpz_t factor, int d, int a, int X0, int Y0, int Z0)
{
	mpz_t SX, SY, SZ;
    mpz_t A, B, B2, B3, C, dC, B2mC, D, CaD, B2mCmD, E, EmD, F, AF, G, AG, aC, DmaC, H, Hx2, J;
    mpz_t X0aY0, X0aY0xB, X0aY0xB_mCmD, X, Y, Z, mulmod;
    
    mpz_init(A); mpz_init(B); mpz_init(B2); mpz_init(B3); mpz_init(C); mpz_init(dC); mpz_init(B2mC); 
    mpz_init(D); mpz_init(CaD); mpz_init(B2mCmD); mpz_init(E); mpz_init(EmD); mpz_init(F); mpz_init(AF);
    mpz_init(G); mpz_init(AG); mpz_init(aC); mpz_init(DmaC); mpz_init(H); mpz_init(Hx2); mpz_init(J);
    mpz_init(X0aY0xB); mpz_init(X0aY0xB_mCmD); mpz_init(X); mpz_init(Y); mpz_init(Z);
    
	mpz_init_set_ui(SX, X0);
	mpz_init_set_ui(SY, Y0);
	mpz_init_set_ui(SZ, Z0);
    
	mpz_init(mulmod);
//int len=mpz_sizeinbase(N,2);if(len > 64 && len < 128) cout << mpz_get_str(NULL,10,N) << endl; 
    int L = mpz_sizeinbase(S, 2); // exact for base = power of 2
    
    // Scalar multiplication using double & add algorithm
    // doubling formula: [2](x:y:z) = ((B-C-D)*J:F*(E-D):F*J)
    for(int i = 2; i <= L; i++) {
        // double
        mpz_add(B, SX, SY);
        mpz_mul(mulmod, B, B); mpz_mod(B2, mulmod, N);
        mpz_mul(mulmod, SX, SX); mpz_mod(C, mulmod, N);
        mpz_mul(mulmod, SY, SY); mpz_mod(D, mulmod, N);
        mpz_mul_ui(E, C, a);
        mpz_add(F, E, D);
        mpz_mul(mulmod, SZ, SZ); mpz_mod(H, mulmod, N);
        mpz_mul_2exp(Hx2, H, 1);
        mpz_sub(J, F, Hx2);
        mpz_add(CaD, C, D);
        mpz_sub(B2mCmD, B2, CaD);
        mpz_mul(X, B2mCmD, J);
        mpz_sub(EmD, E, D);
        mpz_mul(Y, F, EmD);
        mpz_mul(Z, F, J);
        mpz_mod(SX, X, N);
        mpz_mod(SY, Y, N);
        mpz_mod(SZ, Z, N);
        if(mpz_tstbit(S, L - i) == 1) {
            // add
            mpz_mul_ui(A, SZ, Z0);
            mpz_add(B, SX, SY);
            mpz_mul(mulmod, A, A); mpz_mod(B3, mulmod, N);
            mpz_mul_ui(C, SX, X0);
            mpz_mul_ui(D, SY, Y0);
            mpz_mul_ui(dC, C, d);
            mpz_add(CaD, C, D);
            mpz_mul(mulmod, dC, D); mpz_mod(E, mulmod, N);
            mpz_sub(F, B3, E);
            mpz_add(G, B3, E);
            mpz_mul_ui(mulmod, B, X0+Y0); mpz_mod(X0aY0xB, mulmod, N);
            mpz_sub(X0aY0xB_mCmD, X0aY0xB, CaD);
            mpz_mul(mulmod, A, F); mpz_mod(AF, mulmod, N);
            mpz_mul(X, AF, X0aY0xB_mCmD);
            mpz_mul(mulmod, A, G); mpz_mod(AG, mulmod, N);
            mpz_mul_ui(aC, C, a);
            mpz_sub(DmaC, D, aC);
            mpz_mul(Y, AG, DmaC);
            mpz_mul(Z, F, G);
            mpz_mod(SX, X, N);
            mpz_mod(SY, Y, N);
            mpz_mod(SZ, Z, N);
        }
    }
	// try to retrieve factor
	mpz_gcd(factor, N, SX);

	mpz_clear(mulmod); mpz_clear(SZ); mpz_clear(SY); mpz_clear(SX);
    mpz_clear(A); mpz_clear(B); mpz_clear(B2); mpz_clear(B3); mpz_clear(C); mpz_clear(dC); mpz_clear(B2mC); 
    mpz_clear(D); mpz_clear(CaD); mpz_clear(B2mCmD); mpz_clear(E); mpz_clear(EmD); mpz_clear(F); mpz_clear(AF);
    mpz_clear(G); mpz_clear(AG); mpz_clear(aC); mpz_clear(DmaC); mpz_clear(H); mpz_clear(Hx2); mpz_clear(J);
    mpz_clear(X0aY0xB); mpz_clear(X0aY0xB_mCmD); mpz_clear(X); mpz_clear(Y); mpz_clear(Z);
	
	bool result = mpz_cmpabs_ui(factor, 1) > 0 && mpz_cmpabs(factor, N) < 0;
	//if (result) cout << endl << endl << "\t\t\tEECM worked!!!!" << endl << endl << flush;
	return result;
}



/* ScalarMultiplyEdwards (__int128 version)
 * 
 * Multiply a point [X0:Y0:Z0] on a twisted edwards curve by a scalar multiple
 * d	d parameter of twisted Edwards curve
 * a	a parameter of twisted Edwards curve
 * X0,Y0,Z0	point on curve to multiply, in projective coordinates
 * N	we work modulo N
 * S	scalar multiple
 * L	length of S in bits");
*/
bool EECM_int128(__int128 N, mpz_t S, __int128 &factor, int d, int a, int X0, int Y0, int Z0, int64_t mulo, int64_t muhi)
{
	__int128 SX, SY, SZ;
    __int128 A, B, B2, B3, C, dC, B2mC, D, CaD, B2mCmD, E, EmD, F, AF, G, AG, aC, DmaC, H, Hx2, J;
    __int128 X0aY0, X0aY0xB, X0aY0xB_mCmD, X, Y, Z, mulmod;
   
   	SX = X0;
	SY = Y-1;
	SZ = Z0;
    int L = mpz_sizeinbase(S, 2); // exact for base = power of 2
    
    // Scalar multiplication using double & add algorithm
    // doubling formula: [2](x:y:z) = ((B-C-D)*J:F*(E-D):F*J)
    for(int i = 2; i <= L; i++) {
        // double
        B = SX + SY;
        B2 = (B * B) % N;
		C = (SX * SX) % N;
		D = (SY * SY) % N;
		E = C * a;
		F = E + D;
		H = (SZ * SZ) % N;
		Hx2 = H << 1;
		J = F - Hx2;
		CaD = C + D;
		B2mCmD = B2 - CaD;
		X = B2mCmD * J;
		EmD = E - D;
		Y = F * EmD;
		Z = F * J;
		SX = X % N;
		SY = Y % N;
		SZ = Z % N;
        if(mpz_tstbit(S, L - i) == 1) {
            // add
			A = SZ * Z0;
			B = SX + SY;
			B3 = (A * A) % N;
			C = SX * X0;
			D = SY * Y0;
			dC = C * d;
			CaD = C + D;
			E = (dC * D) % N;
			F = B3 - E;
			G = B3 + E;
			X0aY0xB = (B * (X0 + Y0)) % N;
			X0aY0xB_mCmD = X0aY0xB - CaD;
			AF = (A * F) % N;
			X = AF * X0aY0xB_mCmD;
			AG = (A * G) % N;
			aC = C * a;
			DmaC = D - aC;
			Y = AG * DmaC;
			Z = F * G;
			SX = X % N;
			SY = Y % N;
			SZ = Z % N;
        }
    }
	// try to retrieve factor
	factor = gcd128(N, SX);

	bool result = (factor > 0) && (factor < N);
	//if (result) cout << endl << endl << "\t\t\tEECM worked!!!!" << endl << endl << flush;
	return result;
}


inline __int128 make_int128(uint64_t lo, uint64_t hi)
{
	__int128 N = hi;
	N = N << 64;
	N += lo;
	return N;
}


