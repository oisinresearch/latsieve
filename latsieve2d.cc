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

int latsieve2d(int64_t* f, int degfq, int degfp, int64_t q, int l, int* allp, int nump, int* s, int* num_smodp,
		 keyval* M, int Mlen, int* B);
void histogram(keyval*M, uint8_t* H, int len);
bool lattice_sorter(keyval const& kv1, keyval const& kv2);
void csort(keyval* M, keyval* L, int* H, int len);
inline int floordiv(int a, int b);
inline int64_t modinv(int64_t x, int64_t m);
inline int nonzerolcm(int u1, int u2);
inline int gcd(int a, int b);
inline __int128 gcd128(__int128 a, __int128 b);
inline int minnonneg(int u, int v);
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

	if (argc != 13) {
		cout << endl << "Usage: ./latsieve3d inputpoly fbbits factorbasefile B1 B2 qmin qmax th0 th1 lpbbits cofmaxbits qside" << endl << endl;
		return 0;
	}

	cout << "# ";
	for (int i = 0; i < argc; i++) cout << argv[i] << " ";
	cout << endl;

	bool verbose = false;
		
	if (verbose) cout << endl << "Reading input polynomial in file " << argv[1] << "..." << flush;
	//vector<mpz_class> fpoly;
	//vector<mpz_class> gpoly;
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
		if (verbose) cout << line << endl;
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
		if (verbose) cout << line << endl;
		read = static_cast<bool>(getline(file, line));
	}
	//int degg = gpoly.size();
	file.close();
	//mpz_clear(c);
	if (verbose) cout << endl << "Complete.  Degree f0 = " << degf << ", degree f1 = " << degg << "." << endl;

	if (verbose) cout << endl << "Starting sieve of Eratosthenes for small primes..." << endl;
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
	int* sieves0 = new int[degf * nump]();
	int* sievep0 = new int[nump]();
	int* sievenum_s0modp = new int[nump]();
	int* sieves1 = new int[degg * nump]();
	int* sievep1 = new int[nump]();
	int* sievenum_s1modp = new int[nump]();
	// load factor base
	if (verbose) cout << endl << "Loading factor base..." << endl;
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
			sieves1[i*degg + j++] = atoi(substr.c_str());
		}
		sievenum_s1modp[i] = j;
	}
	timetaken += ( clock() - start ) / (double) CLOCKS_PER_SEC;
	if (verbose) cout << "Complete.  Time taken: " << timetaken << "s" << endl;
	if (verbose) cout << "There are " << k0 << " factor base primes on side 0." << endl;
	if (verbose) cout << "There are " << k1 << " factor base primes on side 1." << endl;
	fbfile.close();

	int64_t p0max = sievep0[k0-1];
	int64_t p1max = sievep1[k1-1];
	
	int B[2] = { 12, 12 };
	if (argc >= 5) B[0] = atoi(argv[4]);
	if (argc >= 6) B[1] = atoi(argv[5]);
	int B1bits = B[0]; int B2bits = B[1];
	int B1 = 1<<B1bits; int B2 = 1<<B2bits;
	size_t Mlen = (B1*2l*B2);	// require positive x coordinate
	Mlen = 1500000000;
	//Mlen = (size_t)(2.3f * Mlen);	// upper bound on number of vectors in sieve box
	keyval* M = new keyval[Mlen];	// lattice { id, logp } pairs
	//keyval* L = new keyval[Mlen];	// copy of M
	uint8_t* H = new uint8_t[Mlen];	// histogram
	vector<int> rel;
	// clear M
	//cout << "Clearing memory..." << endl;
	//start = clock();
	//for (int j = 0; j < Mlen; j++) M[j] = (keyval){ 0, 0 };
	//memset(M, 0, Mlen * sizeof(keyval));
	//timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
	//cout << "Memory cleared. Time taken: " << timetaken << "s" << endl;

	int64_t qmin; int64_t qmax; mpz_t qmpz; mpz_init(qmpz);
	mpz_t* pi = new mpz_t[8]; for (int i = 0; i < 8; i++) mpz_init(pi[i]);
	int64_t* r = new int64_t[degf]();
	int64_t* fq = new int64_t[degf+1]();
	mpz_poly f0; mpz_poly f1; mpz_poly A;
	mpz_poly_init(f0, degf); mpz_poly_init(f1, degg); mpz_poly_init(A, 3);
	mpz_poly_set_mpz(f0, fpoly, degf);
	mpz_poly_set_mpz(f1, gpoly, degg);
	mpz_t N0; mpz_t N1;
	mpz_init(N0); mpz_init(N1);
	stringstream stream;
	mpz_t lpb; mpz_init(lpb);
	mpz_t factor; mpz_init(factor); mpz_t p1; mpz_t p2; mpz_init(p1); mpz_init(p2); mpz_t t; mpz_init(t); 
	if (argc >= 7) qmin = strtoll(argv[6], NULL, 10);	// atoi(argv[7]);
	if (argc >= 8) qmax = strtoll(argv[7], NULL, 10);	// atoi(argv[8]);
	uint8_t th0 = 70;
	if (argc >= 9) th0 = atoi(argv[8]);
	uint8_t th1 = 70;
	if (argc >= 10) th1 = atoi(argv[9]);
	int lpbits = 29;
	if (argc >= 11) lpbits = atoi(argv[10]);
	int cofmaxbits = 52;
	if (argc >= 12) cofmaxbits = atoi(argv[11]);
	int64_t cofmax = 1l << cofmaxbits;
	//mpz_t S; mpz_init(S); GetlcmScalar(cofmax, S, primes, 669);	// max S = 5000
	char* str2 = (char*)malloc(20*sizeof(char));
	int qside = atoi(argv[12]);

	int64_t q = qmin;
	while (q < qmax) {
		mpz_set_ui(qmpz, q);
		mpz_nextprime(qmpz, qmpz);
		q = mpz_get_ui(qmpz);
		mpz_t* fpolyside = fpoly;
		int degfside = degf;
		if (qside == 1) { fpolyside = gpoly; degfside = degg; }
		for (int i = 0; i <= degfside; i++) fq[i] = mpz_mod_ui(r0, fpolyside[i], q);
		int numl = polrootsmod(fq, degfside, r, q);
		if (numl == 0 || q > qmax) continue;
		
		// sieve side 0
		cout << "# Starting sieve on side 0";
		if (qside == 0) cout << " for special-q " << q;
		cout << "..." << endl;
		start = clock();
		int m = latsieve2d(fq, degfside, degf, q, 0, sievep0, k0, sieves0, sievenum_s0modp, M, Mlen, B);
		timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
		cout << "# Finished! Time taken: " << timetaken << "s" << endl;
		cout << "# Size of lattice point list is " << m << "." << endl;
		cout << "# Constructing histogram..." << endl;
		start = clock();
		//std::stable_sort(M, M + m, &lattice_sorter);
		histogram(M, H, m);
		timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
		cout << "# Finished! Time taken: " << timetaken << "s" << endl;
		int R0 = 0;
		int B1xB2 = B1*B2;
		int B2x2bits = B2bits + 1;
		int B2x2 = 2*B2;
		rel.clear();
		for (int i = 0; i < m; i++) {
			if (H[i] > th0) {
				rel.push_back(i);
				R0++;
			}
		}
		cout << "# " << R0 << " candidates on side 0." << endl;
		// sieve side 1
		cout << "# Starting sieve on side 1";
		if (qside == 1) cout << " for special-q " << q;
		cout << "..." << endl;
		start = clock();
		m = latsieve2d(fq, degfside, degg, q, 0, sievep1, k1, sieves1, sievenum_s1modp, M, Mlen, B);
		timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
		cout << "# Finished! Time taken: " << timetaken << "s" << endl;
		cout << "# Size of lattice point list is " << m << "." << endl;
		cout << "# Constructing histogram..." << endl;
		start = clock();
		//std::stable_sort(M, M + m, &lattice_sorter);
		histogram(M, H, m);
		timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
		cout << "# Finished! Time taken: " << timetaken << "s" << endl;
		int R1 = 0;
		for (int i = 0; i < m; i++) {
			if (H[i] > th1) {
				rel.push_back(i);
				R1++;
			}
		}
		cout << "# " << R1 << " candidates on side 1." << endl;
		// sort hits
		sort(rel.begin(), rel.end());
		
		// compute special-q lattice L 
		int64_t L[4];
		numl = polrootsmod(fq, degfside, r, q);
		int l = 0;
		L[0] = q; L[1] = -r[l];
		L[2] = 0; L[3] = 1;
		int64L2(L, 2);	// LLL reduce L, time log(q)^2
		
		// print list of potential relations
		int R = 0;
		for (int i = 0; i < (int)(rel.size()-1); i++)
		{
			if (rel[i] == rel[i+1] && rel[i] != 0) {
				int x = rel[i] % B1;
				int y = ((rel[i] >> B1bits) % B2x2) - B2;
				if (x != 0 || y != 0) {
					if (R < 10) {
						// compute [a,b,c]
						int a = L[0]*x+L[1]*y;
						int b = L[2]*x+L[3]*y;
						cout << rel[i] << ": " << a << "," << b << endl;
					}
					R++;
				}
			}
		}
		cout << "# " << R << " potential relations found." << endl;
	   
		// compute and factor resultants as much as possible, leave large gcd computation till later.
		mpz_ui_pow_ui(lpb, 2, lpbits);
		int BASE = 16;
		stack<mpz_t*> QN; stack<int> Q; int algarr[3]; mpz_t* N;
		start = clock();
		R = 0;
		if (verbose) cout << "Starting cofactorizaion..." << endl;
		for (int i = 0; i < (int)(rel.size()-1); i++)
		{
			if (rel[i] == rel[i+1] && rel[i] != 0) {
				int x = rel[i] % B1;
				int y = ((rel[i] >> B1bits) % B2x2) - B2;
				if (x != 0 || y != 0) {
					// compute [a,b,c]
					int a = L[0]*x+L[1]*y;
					int b = L[2]*x+L[3]*y;
	
					int content = gcd(a, b);
					a = a/content; b = b/content;

					//cout << "[a, b] = [" << a << ", " << b << "]" << endl;
					mpz_poly_setcoeff_si(A, 0, a);
					mpz_poly_setcoeff_si(A, 1, b);
					mpz_poly_resultant(N0, f0, A);
					mpz_poly_resultant(N1, f1, A);
					mpz_abs(N0, N0);
					mpz_abs(N1, N1);
					//cout << mpz_get_str(NULL, 10, N0) << endl;
					//cout << mpz_get_str(NULL, 10, N1) << endl;
					string str = to_string(a) + "," + to_string(b) + ":";
					
					// trial division on side 0
					int p = primes[0]; int k = 0; 
					while (p < sievep0[k0-1]) {
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
								while (sievep0[k] < 1000) k++;
							}
						}
						else {
							p = sievep0[++k];
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
					str += (qside == 0 && p0max < q && cofactor ? "," : "");
					// check cofactor on side 0
					if (cofactor) {
						if (mpz_cmpabs_ui(N0, cofmax) > 0) { isrel = false; continue; }
						if (mpz_probab_prime_p(N0, 30) == 0) {  // cofactor definitely composite
							str += "pq";
						}
						else {	// cofactor prime
							str += "Q";
						}
					}
					
					str += ":";

					// trial division on side 1
					p = primes[0]; k = 0;
					while (p < sievep1[k1-1]) {
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
								while (sievep1[k] < 1000) k++;
							}
						}
						else {
							p = sievep1[++k];
						}
					}
					if (mpz_fdiv_ui(N1, q) == 0 && qside == 1)  {
						mpz_divexact_ui(N1, N1, q);
						stream.str("");
						stream << hex << q;
						str += stream.str();
					}
					cofactor = true;
					if (mpz_cmp_ui(N1, 1) == 0) { cofactor = false; }
					str += (qside == 1 && p1max < q && cofactor ? "," : "");
					// check cofactor on side 1
					if (cofactor) {
						if (mpz_cmpabs_ui(N1, cofmax) > 0) { isrel = false; continue; }
						if (mpz_probab_prime_p(N0, 30) == 0) {  // cofactor definitely composite
							str += "pq";
						}
						else {	// cofactor prime
							str += "Q";
						}
					}

					if (isrel) { cout << str << endl; R++; }
				}
			}
		}
		timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
		cout << "# Finished! Cofactorization took " << timetaken << "s" << endl;
		cout << "# " << R << " actual relations found." << endl;
		//cout << "lpb = " << mpz_get_str(NULL, 10, lpb) << endl;
	}

	free(str2);
	mpz_clear(t); mpz_clear(p2); mpz_clear(p1);
	mpz_clear(factor);
	mpz_clear(lpb);
    mpz_clear(N1); mpz_clear(N0);
    mpz_poly_clear(A); mpz_poly_clear(f1); mpz_poly_clear(f0);
	delete[] fq;
	delete[] r;
	for (int i = 0; i < 8; i++) mpz_clear(pi[i]); delete[] pi;
	mpz_clear(qmpz);
	delete[] H;
	//delete[] L;
	delete[] M;
	delete[] sievenum_s1modp;
	delete[] sievep1;
	delete[] sieves1;
	delete[] sievenum_s0modp;
	delete[] sievep0;
	delete[] sieves0;
	mpz_clear(r0);
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
	cout << "calculating start indices..." << endl;
	for (int i = 0; i < len; i++) {
		int oldH = H[i];
		H[i] = m;
		m += oldH;
	}
	cout << "writing output..." << endl;
	// write output
	for (int i = 0; i < len; i++) {
		memcpy(&M[H[L[i].id]], &L[i], sizeof(keyval));
		H[L[i].id]++;
	}
}


inline void mat2x2prod(int64_t* L1, int64_t* L2, int64_t* L3)
{
	L3[0] = L1[0]*L2[0] + L1[1]*L2[2];
	L3[1] = L1[0]*L2[1] + L1[1]*L2[3];
	L3[2] = L1[2]*L2[0] + L1[3]*L2[2];
	L3[3] = L1[2]*L2[1] + L1[3]*L2[3];
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


bool lattice_sorter(keyval const& kv1, keyval const& kv2)
{
	return kv2.id < kv1.id;
}


int latsieve2d(int64_t* f, int degfq, int degfp, int64_t q, int l, int* allp, int nump, int* s, int* num_smodp,
			 keyval* M, int Mlen, int* B)
{
	int64_t L[4];
	int64_t L2[4];
	int64_t L3[4];
	int64_t* r = new int64_t[degfq]();
	int numl = polrootsmod(f, degfq, r, q);
	int B1bits = B[0];
	int B2bits = B[1];
	int B1 = 1<<B1bits;
	int B2 = 1<<B2bits;
	int B1max = B1 - 1;
	int B2max = B2 - 1;

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

	L[0] = q; L[1] = -r[l];
	L[2] = 0; L[3] = 1;

	int64L2(L, 2);	// LLL reduce L, time log(q)^2

	// print special-q lattice
	//cout << "special-q lattice [" << L[0] << "," << L[1] << "," << L[2] << ";";
	//cout << L[3] << "," << L[4] << "," << L[5] << ";";
	//cout << L[6] << "," << L[7] << "," << L[8] << "]" << endl;

	int64_t qLinv[4];
	// compute adjugate of L^-1
	qLinv[0] = L[3];
	qLinv[1] = -L[1];
	qLinv[2] = -L[2];
	qLinv[3] = L[0];

	int i = 30; int m = 0;
	cout << "# Starting sieve at prime " << allp[i] << endl;
	while (i < nump) {
		int64_t p = allp[i];
		if (p == q) { i++; continue; }
		//cout << p << "," << m << endl;
		uint8_t logp = log2f(p);
		__int128 rl = mod(-r[l], q);
		for (int k = 0; k < num_smodp[i]; k++) {
			int64_t n = p*q;
			__int128 sk = mod(-s[i*degfp+k],p);
			int64_t t = q*( (sk * modpinvq[i]) % p ) + p*( (rl * modqinvp[i]) % q ); // CRT
			if (t >= n) t -= n;
			L2[0] = n; L2[1] = t;
			L2[2] = 0; L2[3] = 1;
			int64L2(L2,2);	// LLL reduce L, time log(n)^2
			mat2x2prod(qLinv, L2, L3);	// L3 =  qLinv*L2
            for (int ii = 0; ii < 4; ii++) L3[ii] /= q;
            int64L2(L3,2);
			int u1 = L3[0]; int u2 = L3[2];
			int v1 = L3[1]; int v2 = L3[3];

			// enumerate lattice vectors (x,y) in box [0,B[x[-B,B[
			int x = 0; int y = 0;
			int j1, j2, jmin;
			int s1 = x; int s2 = y;
			// move 'forward' in plane
			bool inplane = true;
			while (inplane) {
				// pin vector to start of row
				int u1B = u1 < 0 ? B1max : 0; int u2B = u2 < 0 ? B2max : -B2;
				j1 = u1 == 0 ? -1 : (x - u1B) / u1;
				j2 = u2 == 0 ? -1 : (y - u2B) / u2;
				jmin = minnonneg(j1, j2);
				x -= jmin*u1; y -= jmin*u2;
				if (x >= 0 && x < B1 && y >= -B2 && y < B2) {
					int u1B = u1 < 0 ? 0 : B1max; int u2B = u2 < 0 ? -B2 : B2max;
					j1 = u1 == 0 ? -1 : (u1B - x) / u1;
					j2 = u2 == 0 ? -1 : (u2B - y) / u2;
					jmin = minnonneg(j1, j2);
					for (int j = 0; j <= jmin; j++) {
						int id = x + ((y + B2) << B1bits);
						M[m++] = (keyval){ id, logp };
						x += u1; y += u2;
					}
					// move by '1-transition' vector
					x += v1; y += v2;
				}
				else {
					inplane = false;	
				}
			}
			x = s1 - v1; y = s2 - v2;
			// move 'backward' in plane
			inplane = true;
			while (inplane) {
				// pin vector to start of row
				int u1B = u1 < 0 ? B1max : 0; int u2B = u2 < 0 ? B2max : -B2;
				j1 = u1 == 0 ? -1 : (x - u1B) / u1;
				j2 = u2 == 0 ? -1 : (y - u2B) / u2;
				jmin = minnonneg(j1, j2);
				x -= jmin*u1; y -= jmin*u2;
				if (x >= 0 && x < B1 && y >= -B2 && y < B2) {
					int u1B = u1 < 0 ? 0 : B1max; int u2B = u2 < 0 ? -B2 : B2max;
					j1 = u1 == 0 ? -1 : (u1B - x) / u1;
					j2 = u2 == 0 ? -1 : (u2B - y) / u2;
					jmin = minnonneg(j1, j2);
					for (int j = 0; j <= jmin; j++) {
						int id = x + ((y + B2) << B1bits);
						M[m++] = (keyval){ id, logp };
						x += u1; y += u2;
					}
					// move by '1-transition' vector
					x -= v1; y -= v2;
				}
				else {
					inplane = false;
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


inline int floordiv(int a, int b)
{
    int d = a / b;
    return d * b == a ? d : d - ((a < 0) ^ (b < 0));
}


inline int nonzerolcm(int u1, int u2)
{
	if (u1 == 0) u1 = 1;
	if (u2 == 0) u2 = 1;
	int g1 = gcd(u1, u2);
	return (u1 * u2) / g1;
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


inline int minabs(int u, int v, int w)
{
	int m = u;
	if (abs(v) < abs(m)) m = v;
	if (abs(w) < abs(m)) m = w;
	return m;
}


inline int minnonneg(int u, int v)
{
	int m = u;
	if (m < 0) m = v;
	if (v >= 0 && v < m) m = v;
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
	//if (result) cout << endl << endl << "\t\t\tP-1 worked!!!!" << endl << endl;
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
	//if (result) cout << endl << endl << "\t\t\tP-1 worked!!!!" << endl << endl;
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
	//if (result) cout << endl << endl << "\t\t\tEECM worked!!!!" << endl << endl;
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
	//if (result) cout << endl << endl << "\t\t\tEECM worked!!!!" << endl << endl;
	return result;
}


inline __int128 make_int128(uint64_t lo, uint64_t hi)
{
	__int128 N = hi;
	N = N << 64;
	N += lo;
	return N;
}


