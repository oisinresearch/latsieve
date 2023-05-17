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
#include "mpz_poly_Fq.h"

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

struct sievedata {
	int fbb;
	int k[2];
	int64_t* q0sieve0;
	int64_t* q1sieve0;
	int64_t* q2sieve0;
	int64_t* q4sieve0;
	int64_t* q0sieve1;
	int64_t* q1sieve1;
	int64_t* q2sieve1;
	int64_t* q4sieve1;
	// special q info
	int64_t q; int qtype[8];
	int64_t r[8], R[8], m[4], a0[4], a1[4], b0[4], b1[4];
};

//int latsieve4d(int64_t* h, int degh, int degfh_t, int64_t q, int64_t r, int64_t R,
//		int* allp, int nump, int* s, int* S, int* num_Smodp, keyval* M, int Mlen, int* B);
int latsieve6d(int n, sievedata info, int side, int* allp, int nump,
	keyval* M, int Mlen, int* B);
int populate_q(sievedata* info, int side, mpz_poly h0, mpz_poly f0, mpz_poly f1,
	mpz_poly_bivariate F0, mpz_poly_bivariate F1);
void histogram(keyval*M, uint8_t* H, int len);
bool lattice_sorter(keyval const& kv1, keyval const& kv2);
void csort(keyval* M, keyval* L, int* H, int len);
inline int64_t floordiv(int64_t a, int64_t b);
inline int64_t modinv(int64_t x, int64_t m);
inline __int128 gcd128(__int128 a, __int128 b);
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
inline void getab5d(int u1, int u2, int u3, int u4, int u5, int u6, int v1, int v2, int v3,
				int v4,	int v5, int v6, int x, int y, int z, int r, int s, int t,
				int B1, int B2, int B3, int B4, int B5, int B6, int* a, int* b);
int64_t subdet(int64_t M[36], int d, int j);
int64_t matdet(int64_t M[25], int d);
inline int64_t floordiv(int64_t a, int64_t b);
inline void getab(int u1, int u2, int u3, int v1, int v2, int v3, int x, int y, int z, 
				int B1, int B2, int B3, int* a, int* b);
inline void getab4d(int u1, int u2, int u3, int u4, int v1, int v2, int v3, int v4,
				int x, int y, int z, int t, int B1, int B2, int B3, int B4, int* a, int* b);
inline void getab5d(int u1, int u2, int u3, int u4, int u5, int v1, int v2, int v3, int v4,
				int v5, int x, int y, int z, int s, int t, int B1, int B2, int B3, int B4,
				int B5, int* a, int* b);
inline void getab6d(int u1, int u2, int u3, int u4, int u5, int u6, int v1, int v2, int v3,
				int v4,	int v5, int v6, int x, int y, int z, int r, int s, int t,
				int B1, int B2, int B3, int B4, int B5, int B6, int* a, int* b);
inline int nonzerolcm(int u1, int u2, int u3);
inline int gcd(int a, int b);
inline int64_t nonzerolcm4d(int u1, int u2, int u3, int u4);
int lcm(int a, int b);
inline int64_t nonzerolcm5d(int u1, int u2, int u3, int u4, int u5);
inline int64_t nonzerolcm6d(int u1, int u2, int u3, int u4, int u5, int u6);
inline bool threespaceintersectsbox6d(int64_t nx, int64_t ny, int64_t nz, int64_t nr, 
	int64_t ns, int64_t nt, int x, int y, int z, int r, int s, int t, int B1, int B2,
	int B3, int B4, int B5, int B6);
inline bool fourspaceintersectsbox6d(int64_t nx, int64_t ny, int64_t nz, int64_t nr, 
	int64_t ns, int64_t nt, int x, int y, int z, int r, int s, int t, int B1, int B2,
	int B3, int B4, int B5, int B6);
inline bool fivespaceintersectsbox6d(int64_t nx, int64_t ny, int64_t nz, int64_t nr,
	int64_t ns, int64_t nt, int x, int y, int z, int r, int s, int t, int B1, int B2,
	int B3, int B4, int B5, int B6);
inline int minnonneg6d(int u, int v, int w, int r, int s, int t);
inline void matmul6x6(int64_t C[36], int64_t A[36], int64_t B[36]);
void matadj6x6(int64_t M[36], int64_t Madj[36]);
void printmat6x6(int64_t M[36]);


int main(int argc, char** argv)
{
	// set constant
	MASK64 = 1L;
	MASK64 = MASK64 << 64;
	MASK64 = MASK64 - 1L;

	//cout << (uint64_t)(MASK64) << " " << (uint64_t)(MASK64 >> 64) << endl;

	if (argc != 17) {
		cout << endl << "Usage: ./latsieve6d inputpoly sieve_bound factorbasefile "
			"B1 B2 B3 B4 B5 B6 qmin qmax th0 th1 lp_bound cofacscalar qside" << endl << endl;
		return 0;
	}

	cout << "# ";
	for (int i = 0; i < argc; i++) cout << argv[i] << " ";
	cout << endl;

	bool verbose = true;// false;
		
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
	if (verbose) cout << endl << "Side 0 polynomial fh_y (ascending coefficients)" << endl;
	while (getline(file, line) && line.substr(0,3) == "fhy" ) {
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
	if (verbose) cout << endl << "Side 1 polynomial gh_y: (ascending coefficients)" << endl;
	while (read && line.substr(0,3) == "ghy" ) {
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
	int fbb = 0;
	if (argc >=3) fbb = atoi(argv[2]);
	int max = fbb; // 10000000;// 65536;
	char* sieve = new char[max+1]();
	int* allp = new int[2097152]; //int[1077871]; // int[155611]; //new int[809228];	//new int[6542]; 	// 2039 is the 309th prime, largest below 2048
	for (int i = 2; i <= sqrt(max); i++)
		if(!sieve[i])
			for (int j = i*i; j <= max; j += i)
				if(!sieve[j]) sieve[j] = 1;
	int nump = 0;
	for (int i = 2; i <= max-1; i++)
		if (!sieve[i])
			allp[nump++] = i;
	if (verbose) cout << "Complete." << endl;

	// load factor base
	sievedata info;
	std::clock_t start; double timetaken = 0;
	if (verbose) cout << endl << "Loading factor base..." << endl << flush;
	start = clock();
	ifstream fbfile(argv[3]);
	start = clock();
	// read fbb
	getline(fbfile, line);
	info.fbb = atoi(line.c_str());
	// read k[0]
	getline(fbfile, line);
	info.k[0] = atoi(line.c_str());
	int imax0 = info.k[0];
	// read q2sieve0
	info.q2sieve0 = new int64_t[imax0*3]();
	for (int i = 0; i < info.k[0]; i++) {
		getline(fbfile, line);
		stringstream ss(line);
		string substr;
		getline(ss, substr, ',');
		info.q2sieve0[3*i] = atoi(substr.c_str());
		int j = 0;
		while( ss.good() ) {
			getline( ss, substr, ',' );
			info.q2sieve0[3*i + (++j)] = atoi(substr.c_str());
		}
	}
	// read k[1]
	getline(fbfile, line);
	info.k[1] = atoi(line.c_str());
	int imax1 = info.k[1];
	// read q2sieve1
	info.q2sieve1 = new int64_t[imax1*3]();
	for (int i = 0; i < info.k[1]; i++) {
		getline(fbfile, line);
		stringstream ss(line);
		string substr;
		getline(ss, substr, ',');
		info.q2sieve1[3*i] = atoi(substr.c_str());
		int j = 0;
		while( ss.good() ) {
			getline( ss, substr, ',' );
			info.q2sieve1[3*i + (++j)] = atoi(substr.c_str());
		}
	}

	timetaken += ( clock() - start ) / (double) CLOCKS_PER_SEC;
	if (verbose) cout << "Complete.  Time taken: " << timetaken << "s" << endl;
	if (verbose) cout << "There are " << info.k[0] << " factor base prime ideals on side 0." << endl;
	if (verbose) cout << "There are " << info.k[1] << " factor base prime ideals on side 1." << endl;
	fbfile.close();

	mpz_t r0; mpz_init(r0);
	int B[6] = { 3, 3, 3, 3, 3, 3 };
	if (argc >= 5) B[0] = atoi(argv[4]);
	if (argc >= 6) B[1] = atoi(argv[5]);
	if (argc >= 7) B[2] = atoi(argv[6]);
	if (argc >= 8) B[3] = atoi(argv[7]);
	if (argc >= 9) B[4] = atoi(argv[8]);
	if (argc >= 10) B[5] = atoi(argv[9]);
	int B1bits = B[0]; int B2bits = B[1]; int B3bits = B[2];
	int B4bits = B[3]; int B5bits = B[4]; int B6bits = B[5];
	int B1 = 1<<B1bits; int B2 = 1<<B2bits; int B3 = 1<<B3bits;
	int B4 = 1<<B4bits; int B5 = 1<<B5bits; int B6 = 1<<B6bits;
	uint32_t Mlen = (B1*2*B2*2*B3*2*B4*2*B5*2*B6*2);	// DO NOT! require positive z coordinate
	Mlen = (uint32_t)(2.3f * Mlen);	// upper bound on number of vectors in sieve box
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
	int64_t* h = new int64_t[degh+1]();
	mpz_poly h0; mpz_poly_init(h0, degh);
	mpz_poly f0; mpz_poly f1;
	mpz_poly_init(f0, degfht); mpz_poly_init(f1, degght);
	mpz_poly_set_mpz(f0, fhtpoly, degfht);
	mpz_poly_set_mpz(f1, ghtpoly, degght);
	mpz_t N0; mpz_t N1;
	mpz_init(N0); mpz_init(N1);
	stringstream stream;
	mpz_t lpb; mpz_init(lpb);
	mpz_t factor; mpz_init(factor); mpz_t p1; mpz_t p2; mpz_init(p1); mpz_init(p2);
	mpz_t tmp1; mpz_init(tmp1); 
	if (argc >= 11) qmin = strtoll(argv[10], NULL, 10);	// atoi(argv[7]);
	if (argc >= 12) qmax = strtoll(argv[11], NULL, 10);	// atoi(argv[8]);
	uint8_t th0 = 70;
	if (argc >= 13) th0 = atoi(argv[12]);
	uint8_t th1 = 70;
	if (argc >= 14) th1 = atoi(argv[13]);
	int64_t lpb0 = 0;
	if (argc >= 15) lpb0 = strtoll(argv[14], NULL, 10);
	int cofacS = 1000;
	if (argc >= 16) cofacS = atoi(argv[15]);
	mpz_t S; mpz_init(S); GetlcmScalar(cofacS, S, allp, 669);	// max S = 5000
	char* str2 = (char*)malloc(20*sizeof(char));
	int qside = atoi(argv[16]);
	mpz_poly Fqh_x; mpz_poly_init(Fqh_x, 0);
	mpz_poly_bivariate Aq; mpz_poly_bivariate_init(Aq, 0);
	mpz_poly Aq0; mpz_poly_init(Aq0, 0);
	mpz_t res; mpz_init(res);
	for (int i = 0; i <= degh; i++) mpz_poly_setcoeff(h0, i, hpoly[i]);

	int64_t q = qmin;
	while (q < qmax) {
		mpz_set_ui(qmpz, q);
		mpz_nextprime(qmpz, qmpz);
		q = mpz_get_ui(qmpz);
		info.q = q;

		// calculate special-q ideals for this q
		int nmax = populate_q(&info, qside, h0, f0, f1, F0, F1);
		
		cout << "# (nmax = " << nmax << ")" << endl;

		// note that we skip inert special-q, the coefficients are too big
		if (nmax == 0 || q > qmax) continue;

		int m = 0;
		for (int n = 0; n < nmax; n++) {
			// sieve side 0
			cout << "# Starting sieve on side 0" << flush;
			start = clock();
			int64_t rn = info.r[n]; int64_t Rn = info.R[n];
			if (qside == 0) {
				cout << " for special-q ";
				cout << "(q,r,R) = (" << q << "," << rn << "," << Rn << ")";
			}
			cout << "..." << endl;
			
			// we only allow degree-1 special-q ideals for the moment

			m = latsieve6d(n, info, 0, allp, nump, M, Mlen, B);
			timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
			cout << "# Finished! Time taken: " << timetaken << "s" << endl;
			cout << "# Size of lattice point list is " << m << "." << endl;
			cout << "# Constructing histogram..." << endl;
			start = clock();
			//std::stable_sort(M, M + m, &lattice_sorter);
			histogram(M, H, m);
			timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
			cout << "# Finished! Time taken: " << timetaken << "s" << endl << flush;
			int R0 = 0;
			int B1x2 = 2*B1;
			int B2x2 = 2*B2;
			int B3x2 = 2*B3;
			int B4x2 = 2*B4;
			int B5x2 = 2*B5;
			int B1x2bits = B1bits + 1;
			int B2x2bits = B2bits + 1;
			int B1x2xB2x2bits = B1bits + 1 + B2bits + 1;
			int B1x2xB2x2xB3x2bits = B1bits + 1 + B2bits + 1 + B3bits + 1;
			int B1x2xB2x2xB3x2xB4x2bits = B1bits + 1 + B2bits + 1 + B3bits + 1 + B4bits + 1;
			int B1x2xB2x2xB3x2xB4x2xB5x2bits = B1x2xB2x2xB3x2xB4x2bits + B5bits + 1;
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
			if (qside == 1) {
				cout << " for special-q ";
				cout << "(q,r,R) = (" << q << "," << rn << "," << Rn << ")";
			}
			cout << "..." << endl;
			start = clock();
			//m = latsieve3d(fq, degg, q, 0, sievep1, k1, sieves1, sievenum_s1modp, M, Mlen, B);
			m = latsieve6d(n, info, 1, allp, nump, M, Mlen, B);
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
			L[0]  = q; L[1]  = rn; L[2]  = 0;  L[3]  = Rn; L[4]  = 0;  L[5]  = 0;
			L[6]  = 0; L[7]  = 1;  L[8]  = rn; L[9]  = 0;  L[10] = Rn; L[11] = 0;
			L[12] = 0; L[13] = 0;  L[14] = 1;  L[15] = 0;  L[16] = 0;  L[17] = Rn;
			L[18] = 0; L[19] = 0;  L[20] = 0;  L[21] = 1;  L[22] = 0;  L[23] = 0;
			L[24] = 0; L[25] = 0;  L[26] = 0;  L[27] = 0;  L[28] = 1;  L[29] = 0;
			L[30] = 0; L[31] = 0;  L[32] = 0;  L[33] = 0;  L[34] = 0;  L[35] = 1;
			int64L2(L, 4);	// LLL reduce L, time log(q)^2
			
			// print list of potential relations
			int R = 0;
			for (int i = 0; i < rel.size()-1; i++)
			{
				if (rel[i] == rel[i+1] && rel[i] != 0) {
					int x = (rel[i] % B1x2) - B1;
					int y = ((rel[i] >> B1x2bits) % B2x2) - B2;
					int z = ((rel[i] >> B1x2xB2x2bits) % B3x2) - B3;
					int r = ((rel[i] >> B1x2xB2x2xB3x2bits) % B4x2) - B4;
					int s = ((rel[i] >> B1x2xB2x2xB3x2xB4x2bits) % B5x2) - B5;
					int t = (rel[i] >> B1x2xB2x2xB3x2xB4x2xB5x2bits) - B6;
					if (x != 0 || y != 0 || z != 0 || r != 0 || s != 0 || t != 0) {
						// compute [a,b,c,d]
						//int a = L[0]*x+L[1]*y+L[2]*z+L[3]*t;
						//int b = L[4]*x+L[5]*y+L[6]*z+L[7]*t;
						//int c = L[8]*x+L[9]*y+L[10]*z+L[11]*t;
						//int d = L[12]*x+L[13]*y+L[14]*z+L[15]*t;
						//if (R >= 100 & R < 110)
						//	cout << rel[i] << ": " << a << "," << b << "," << c << "," << d << endl;
						//	cout << "(x,y,z,t): " << x << "," << y << "," << z << "," << t << endl;
						R++;
					}
				}
			}
			cout << "# " << R << " potential relations found." << endl << flush;
		   
			// compute and factor resultants as much as possible, leave large gcd computation till later.
			mpz_set_ui(lpb, lpb0);
			//mpz_ui_pow_ui(lpb, 2, lpbits);
			int BASE = 16;
			stack<mpz_t*> QN; stack<int> Q; int algarr[3]; mpz_t* N;
			start = clock();
			R = 0;
			if (verbose) cout << "Starting cofactorization..." << endl << flush;
			for (int i = 0; i < rel.size()-1; i++)
			{
				if (rel[i] == rel[i+1] && rel[i] != 0) {
					int x = (rel[i] % B1x2) - B1;
					int y = ((rel[i] >> B1x2bits) % B2x2) - B2;
					int z = ((rel[i] >> B1x2xB2x2bits) % B3x2) - B3;
					int r = ((rel[i] >> B1x2xB2x2xB3x2bits) % B4x2) - B4;
					int s = ((rel[i] >> B1x2xB2x2xB3x2xB4x2bits) % B5x2) - B5;
					int t = (rel[i] >> B1x2xB2x2xB3x2xB4x2xB5x2bits) - B6;
					if (x != 0 || y != 0 || z != 0 || t != 0) {
						// compute [a,b,c]
						int a = L[0]*x+L[1]*y+L[2]*z+L[3]*r+L[4]*s+L[5]*t;
						int b = L[6]*x+L[7]*y+L[8]*z+L[9]*r+L[10]*s+L[11]*t;
						int c = L[12]*x+L[13]*y+L[14]*z+L[15]*r+L[16]*s+L[17]*t;
						int d = L[18]*x+L[19]*y+L[20]*z+L[21]*r+L[22]*s+L[23]*t;
						int e = L[24]*x+L[25]*y+L[26]*z+L[27]*r+L[28]*s+L[29]*t;
						int f = L[30]*x+L[31]*y+L[32]*z+L[33]*r+L[34]*s+L[35]*t;
		
						//int64_t D = b*(int64_t)b - 4*a*(int64_t)c;
						//int64_t D = b*(int64_t)b*c*(int64_t)c - 4*b*(int64_t)b*b*d
						//    -4*(int64_t)a*c*(int64_t)c*c + 18*(int64_t)a*b*c*d
						//	-27*(int64_t)a*a*(int64_t)d*d;
						//if (floor(sqrt(D)+0.5)*floor(sqrt(D)+0.5) == D) {
						//	continue;	// a+b*x+c*x^2 is not irreducible over Z
						//}

						int content = gcd(a, b); content = gcd(content, c);
						content = gcd(content, d); content = gcd(content, e);
						content = gcd(content, f);
						a = a/content; b = b/content; c = c/content; d = d/content;
						e = e/content; f = f/content;

						//cout << "[a, b, c] = [" << a << ", " << b << ", " << c << "]" << endl << flush;
						
						mpz_poly_setcoeff_si(Aq0, 0, a);
						mpz_poly_setcoeff_si(Aq0, 1, b);
						mpz_poly_setcoeff_si(Aq0, 2, c);
						mpz_poly_bivariate_setcoeff(Aq, 0, Aq0);
						mpz_poly_setcoeff_si(Aq0, 0, d);
						mpz_poly_setcoeff_si(Aq0, 1, e);
						mpz_poly_setcoeff_si(Aq0, 2, f);
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
						//cout << mpz_get_str(NULL, 10, N0) << endl;
						//cout << mpz_get_str(NULL, 10, N1) << endl;
						string str = to_string(a) + "," + to_string(b) + "," + 
							to_string(c) + "," + to_string(d) + to_string(e) +
							to_string(f) + ":";
						
						// trial division on side 0
						int p = allp[0]; int k = 0; 
						while (p < allp[nump-1]) {
							int valp = 0;
							while (mpz_fdiv_ui(N0, p) == 0) {
								mpz_divexact_ui(N0, N0, p);
								valp++;
								stream.str("");
								stream << hex << p;
								str += stream.str() + ",";
							}
							p = allp[++k];
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
							p = allp[0]; k = 0;
							while (p < allp[nump-1]) {
								int valp = 0;
								while (mpz_fdiv_ui(N1, p) == 0) {
									mpz_divexact_ui(N1, N1, p);
									valp++;
									stream.str("");
									stream << hex << p;
									str += stream.str() + ",";
								}
								p = allp[++k];
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
    mpz_poly_clear(f1); mpz_poly_clear(f0);
	mpz_poly_clear(h0);
	delete[] h;
	for (int i = 0; i < 8; i++) mpz_clear(pi[i]); delete[] pi;
	mpz_clear(qmpz);
	delete[] H;
	//delete[] L;
	delete[] M;
	mpz_clear(r0);
	delete[] allp;
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


int latsieve6d(int n, sievedata info, int side, int* allp, int nump,
	keyval* M, int Mlen, int* B)
{
	int64_t L[36];
	int64_t L0[36];
	int64_t L2[36];
	int64_t L3[36];
	int64_t L4[36];
	int64_t L5[25];
	int B1bits = B[0];
	int B2bits = B[1];
	int B3bits = B[2];
	int B4bits = B[3];
	int B5bits = B[4];
	int B6bits = B[5];
	int B1 = 1<<B1bits;
	int B2 = 1<<B2bits;
	int B3 = 1<<B3bits;
	int B4 = 1<<B4bits;
	int B5 = 1<<B5bits;
	int B6 = 1<<B6bits;
	int B1max = B1 - 1;
	int B2max = B2 - 1;
	int B3max = B3 - 1;
	int B4max = B4 - 1;
	int B5max = B5 - 1;
	int B6max = B6 - 1;
	int B1x2bits = B1bits + 1;
	int B1x2xB2x2bits = B1bits + 1 + B2bits + 1;
	int B1x2xB2x2xB3x2bits = B1bits + 1 + B2bits + 1 + B3bits + 1;
	int B1x2xB2x2xB3x2xB4x2bits = B1x2xB2x2xB3x2bits + B4bits + 1;
	int B1x2xB2x2xB3x2xB4x2xB5x2bits = B1x2xB2x2xB3x2xB4x2bits + B5bits + 1;

	int64_t q = info.q; int64_t r0 = info.r[n]; int64_t R0 = info.R[n];
	// print (q,r) ideal
	//cout << "special-q (" << q << "," << rl << ")" << endl;

	// compute qinvmodp and pinvmodq arrays
	int imax = info.k[side];
	int64_t* qinvmodp = new int64_t[imax]();
	int64_t* pinvmodq = new int64_t[imax]();
	int k2 = 0; int i2 = 0;
	while (true) {
		int64_t p2 = allp[i2];
		int64_t q2 = side == 0 ? info.q2sieve0[3*k2] : info.q2sieve1[3*k2];
		if (k2 < info.k[side] && p2 == q2) {
			if (p2 != q) {
				qinvmodp[k2] = modinv(q, p2);
				pinvmodq[k2] = modinv(p2, q);
			}
			k2++;
		}
		if (q2 > p2) i2++;
		if (k2 == info.k[side]) break;
	}

	L[0]  = q; L[1]  = r0; L[2]  = 0;  L[3]  = R0; L[4]  = 0;  L[5]  = 0;
	L[6]  = 0; L[7]  = 1;  L[8]  = r0; L[9]  = 0;  L[10] = R0; L[11] = 0;
	L[12] = 0; L[13] = 0;  L[14] = 1;  L[15] = 0;  L[16] = 0;  L[17] = R0;
	L[18] = 0; L[19] = 0;  L[20] = 0;  L[21] = 1;  L[22] = 0;  L[23] = 0;
	L[24] = 0; L[25] = 0;  L[26] = 0;  L[27] = 0;  L[28] = 1;  L[29] = 0;
	L[30] = 0; L[31] = 0;  L[32] = 0;  L[33] = 0;  L[34] = 0;  L[35] = 1;

	for (int jj = 0; jj < 36; jj++) L0[jj] = L[jj];
	int64L2(L, 6);	// LLL reduce L, time log(q)^2

	// compute adjugate of L^-1
	int64_t qLinv[36];
	matadj6x6(L, qLinv);

	int mm = 0;
    int i = 40;
    while (i < info.k[side]) {
        int64_t p;
        if (side == 0) p = info.q2sieve0[3*i];
        if (side == 1) p = info.q2sieve1[3*i];
        if (p == q) { i++; continue; }
        //cout << p << "," << mm << endl << flush;
        uint8_t logp = log2f(p);
        __int128 _r1 = side == 0 ? info.q2sieve0[3*i+1] : info.q2sieve1[3*i+1];
        __int128 _R1 = side == 0 ? info.q2sieve0[3*i+2] : info.q2sieve1[3*i+2];
        int64_t h0, h1;
        h0 = q*( (_r1 * qinvmodp[i]) % p ) + p*( (r0 * pinvmodq[i]) % q );
        h1 = q*( (_R1 * qinvmodp[i]) % p ) + p*( (R0 * pinvmodq[i]) % q );	
        L2[0]  = p*q; L2[1]  = h0; L2[2]  = 0;  L2[3]  = h1; L2[4]  = 0;  L2[5]  = 0;
        L2[6]  = 0;   L2[7]  = 1;  L2[8]  = h0; L2[9]  = 0;  L2[10] = h1; L2[11] = 0;
        L2[12] = 0;   L2[13] = 0;  L2[14] = 1;  L2[15] = 0;  L2[16] = 0;  L2[17] = h1;
        L2[18] = 0;   L2[19] = 0;  L2[20] = 0;  L2[21] = 1;  L2[22] = 0;  L2[23] = 0;
        L2[24] = 0;   L2[25] = 0;  L2[26] = 0;  L2[27] = 0;  L2[28] = 1;  L2[29] = 0;
        L2[30] = 0;   L2[31] = 0;  L2[32] = 0;  L2[33] = 0;  L2[34] = 0;  L2[35] = 1;
        int64L2(L2, 6);	// LLL reduce L, time log(n)^2
        matmul6x6(L3, qLinv, L2);	// L3 =  qLinv*L2
        for (int jj = 0; jj < 36; jj++) L3[jj] /= q;
        int64L2(L3, 6);
        int u1, u2, u3, u4, u5, u6;
        int v1, v2, v3, v4, v5, v6;
        int w1, w2, w3, w4, w5, w6;
        int r1, r2, r3, r4, r5, r6;
        int s1, s2, s3, s4, s5, s6;
        int t1, t2, t3, t4, t5, t6;
        u1 = L3[0]; u2 = L3[6];  u3 = L3[12]; u4 = L3[18]; u5 = L3[24]; u6 = L3[30];
        v1 = L3[1]; v2 = L3[7];  v3 = L3[13]; v4 = L3[19]; v5 = L3[25]; v6 = L3[31];
        w1 = L3[2]; w2 = L3[8];  w3 = L3[14]; w4 = L3[20]; w5 = L3[26]; w6 = L3[32];
        r1 = L3[3]; r2 = L3[9];  r3 = L3[15]; r4 = L3[21]; r5 = L3[27]; r6 = L3[33];
        s1 = L3[4]; s2 = L3[10]; s3 = L3[16]; s4 = L3[22]; s5 = L3[28]; s6 = L3[34];
        t1 = L3[5]; t2 = L3[11]; t3 = L3[17]; t4 = L3[23]; t5 = L3[29]; t6 = L3[35];

        // compute normal (cross product)
        L4[0]  = 1;  L4[1]  = 1;  L4[2]  = 1;  L4[3]  = 1;  L4[4]  = 1;  L4[5]  = 1;
        L4[6]  = u1; L4[7]  = u2; L4[8]  = u3; L4[9]  = u4; L4[10] = u5; L4[11] = u6;
        L4[12] = v1; L4[13] = v2; L4[14] = v3; L4[15] = v4; L4[16] = v5; L4[17] = v6;
        L4[18] = w1; L4[19] = w2; L4[20] = w3; L4[21] = w4; L4[22] = w5; L4[23] = w6;
        L4[24] = r1; L4[25] = r2; L4[26] = r3; L4[27] = r4; L4[28] = r5; L4[29] = r6;
        L4[30] = s1; L4[31] = s2; L4[32] = s3; L4[33] = s4; L4[34] = s5; L4[35] = s6;
        int64_t nx = subdet(L4, 6, 0);
        int64_t ny = subdet(L4, 6, 1);
        int64_t nz = subdet(L4, 6, 2);
        int64_t nr = subdet(L4, 6, 3);
        int64_t ns = subdet(L4, 6, 4);
        int64_t nt = subdet(L4, 6, 5);

        // enumerate lattice vectors (x,y,z,r,s,t) in box
        // [-B,B[x[-B,B[x[-B,B[x[-B,B[x[-B,B[x[-B,B[
        int x = 0; int y = 0; int z = 0; int r = 0; int s = 0; int t = 0;
        while (fivespaceintersectsbox6d(nx, ny, nz, nr, ns, nt, x-t1,
            y-t2, z-t3, r-t4, s-t5, t-t6, B1, B2, B3, B4, B5, B6)) {
            x -= t1; y -= t2; z -= t3; r -= t4; s -= t5; t -= t6;
        }
        int ts1 = x; int ts2 = y; int ts3 = z; int ts4 = r; int ts5 = s; int ts6 = t;
        while (fivespaceintersectsbox6d(nx, ny, nz, nr, ns, nt, x, y, z, r, s, t, B1, B2, B3, B4, B5, B6)) { 
            while (fourspaceintersectsbox6d(nx, ny, nz, nr, ns, nt, x-s1,
                y-s2, z-s3, r-s4, s-s5, t-s6, B1, B2, B3, B4, B5, B6)) {
                x -= s1; y -= s2; z -= s3; r -= s4; s -= s5; t -= s6;
            }
            int ss1 = x; int ss2 = y; int ss3 = z; int ss4 = r; int ss5 = s; int ss6 = t;
            while (fourspaceintersectsbox6d(nx, ny, nz, nr, ns, nt, x, y, z, r, s, t, B1, B2, B3, B4, B5, B6)) { 
                while (threespaceintersectsbox6d(nx, ny, nz, nr, ns, nt, x-r1,
                    y-r2, z-r3, r-r4, s-r5, t-r6, B1, B2, B3, B4, B5, B6)) {
                    x -= r1; y -= r2; z -= r3; r -= r4; s -= r5; t -= r6;
                }
                int rs1 = x; int rs2 = y; int rs3 = z; int rs4 = r; int rs5 = s; int rs6 = t;
                while (threespaceintersectsbox6d(nx, ny, nz, nr, ns, nt, x, y, z, r, s, t, B1, B2, B3, B4, B5, B6)) { 
                    int j1, j2, j3, j4, j5, j6, jmin;
                    // pin vector to start of row
                    int w1B = w1 < 0 ? B1max : -B1; int w2B = w2 < 0 ? B2max : -B2;
                    int w3B = w3 < 0 ? B3max : -B3; int w4B = w4 < 0 ? B4max : -B4;
                    int w5B = w5 < 0 ? B5max : -B5; int w6B = w6 < 0 ? B6max : -B6;
                    j1 = w1 == 0 ? -1 : (x - w1B) / w1;
                    j2 = w2 == 0 ? -1 : (y - w2B) / w2;
                    j3 = w3 == 0 ? -1 : (z - w3B) / w3;
                    j4 = w4 == 0 ? -1 : (t - w4B) / w4;
                    j5 = w5 == 0 ? -1 : (t - w5B) / w5;
                    j6 = w6 == 0 ? -1 : (t - w6B) / w6;
                    jmin = minnonneg6d(j1, j2, j3, j4, j5, j6);
                    x -= jmin*w1; y -= jmin*w2; z -= jmin*w3;
                    r -= jmin*w4; s -= jmin*w5; t -= jmin*w6;
                    int ws1 = x; int ws2 = y; int ws3 = z; int ws4 = r; int ws5 = s; int ws6 = t;
                    bool inspace = true;
                    while (inspace) {
                        if (x >= -B1 && x < B1 && y >= -B2 && y < B2 && z >= -B3 && z < B3
                         && r >= -B4 && r < B4 && s >= -B5 && s < B5 && t >= -B6 && t < B6) {
                            int a = 0; int b = 0;
                            getab6d(u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6,
                                x, y, z, r, s, t, B1, B2, B3, B4, B5, B6, &a, &b);
                            x += a*u1 + b*v1;
                            y += a*u2 + b*v2;
                            z += a*u3 + b*v3;
                            r += a*u4 + b*v4;
                            s += a*u5 + b*v5;
                            t += a*u6 + b*v6;
                            int s1 = x; int s2 = y; int s3 = z; int s4 = r; int s5 = s; int s6 = t;
                            // move 'forward' (dir == -1) or 'backward' (dir == 1) in plane
                            for (int dir = -1; dir <= 1; dir += 2) {
                                bool inplane = true;
                                while (inplane) {
                                    // pin vector to start of row
                                    int u1B = u1 < 0 ? B1max : -B1; int u2B = u2 < 0 ? B2max : -B2;
                                    int u3B = u3 < 0 ? B3max : -B3; int u4B = u4 < 0 ? B4max : -B4;
                                    int u5B = u3 < 0 ? B3max : -B5; int u6B = u4 < 0 ? B6max : -B6;
                                    j1 = u1 == 0 ? -1 : (x - u1B) / u1;
                                    j2 = u2 == 0 ? -1 : (y - u2B) / u2;
                                    j3 = u3 == 0 ? -1 : (z - u3B) / u3;
                                    j4 = u4 == 0 ? -1 : (r - u4B) / u4;
                                    j5 = u5 == 0 ? -1 : (s - u5B) / u5;
                                    j6 = u6 == 0 ? -1 : (t - u6B) / u6;
                                    jmin = minnonneg6d(j1, j2, j3, j4, j5, j6);
                                    x -= jmin*u1; y -= jmin*u2; z -= jmin*u3;
									r -= jmin*u4; s -= jmin*u5; t -= jmin*u4;
                                    if (x >= -B1 && x < B1 && y >= -B2 && y < B2
                                          && z >= -B3 && z < B3 && r >= -B4 && r < B4
										  	&& s >= -B5 && s < B5 && t >= -B6 && t < B6) {
                                        int u1B = u1 < 0 ? -B1 : B1max; int u2B = u2 < 0 ? -B2 : B2max;
                                        int u3B = u3 < 0 ? -B3 : B3max; int u4B = u4 < 0 ? -B4 : B4max;
                                        int u5B = u5 < 0 ? -B5 : B5max; int u6B = u6 < 0 ? -B6 : B6max;
                                        j1 = u1 == 0 ? -1 : (u1B - x) / u1;
                                        j2 = u2 == 0 ? -1 : (u2B - y) / u2;
                                        j3 = u3 == 0 ? -1 : (u3B - z) / u3;
                                        j4 = u4 == 0 ? -1 : (u4B - r) / u4;
                                        j5 = u5 == 0 ? -1 : (u5B - s) / u5;
                                        j6 = u6 == 0 ? -1 : (u6B - t) / u6;
                                        jmin = minnonneg6d(j1, j2, j3, j4, j5, j6);
                                        for (int j = 0; j <= jmin; j++) {
                                            int id = (x + B1) + ((y + B2) << B1x2bits) + ((z + B3) << B1x2xB2x2bits)
                                              + ((r + B4) << B1x2xB2x2xB3x2bits)
											  + ((s + B5) << B1x2xB2x2xB3x2xB4x2bits)
											  + ((t + B6) << B1x2xB2x2xB3x2xB4x2xB5x2bits);
                                            M[mm++] = (keyval){ id, logp };
                                            x += u1; y += u2; z += u3; r += u4; s += u5; t += u6;
                                        }
                                        // move by '1-transition' vector
                                        if (dir == -1) { x += v1; y += v2; z += v3; r += v4; s += v5; t += v6; }
                                        else		   { x -= v1; y -= v2; z -= v3; r -= v4; s -= v5; t -= v6; }
                                    }
                                    else {
                                        inplane = false;	
                                    }
                                }
                                x = s1 - v1; y = s2 - v2; z = s3 - v3; r = s4 - v4; s = s5 - v5; t = s6 - v6;
                            }
                            // advance to next w grid plane
                            ws1 += w1; ws2 += w2; ws3 += w3; ws4 += w4; ws5 += w5; ws6 += w6;
							x = ws1; y = ws2; z = ws3; r = ws4; s = ws5; t = ws6;
                        }
                        else {
                            inspace = false;
                        }
                    }
                    // advance to next t grid plane
                    rs1 += r1; rs2 += r2; rs3 += r3; rs4 += r4; rs5 += r5; rs6 += r6;
                    x = rs1; y = rs2; z = rs3; r = rs4; s = rs5; t = rs6;
                }
                // advance to next t grid plane
                ss1 += s1; ss2 += s2; ss3 += s3; ss4 += s4; ss5 += s5; ss6 += s6;
                x = ss1; y = ss2; z = ss3; r = ss4; s = ss5; t = ss6;
            }
            // advance to next t grid plane
            ts1 += t1; ts2 += t2; ts3 += t3; ts4 += t4; ts5 += t5; ts6 += t6;
            x = ts1; y = ts2; z = ts3; r = ts4; s = ts5; t = ts6;
        }
        // advance to next p
        i++;
    }

	// clear memory
	delete[] qinvmodp;
	delete[] pinvmodq;

	return mm;
}

void printmat6x6(int64_t M[36])
{
	cout << "[";
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			cout << M[6*i+j];
			if (j < 5) cout << ",";
			else if (i < 5) cout << ";\\\n";
			else cout << "]" << endl;
		}
	}
}

inline void matmul6x6(int64_t C[36], int64_t A[36], int64_t B[36])
{
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            C[6*i + j] = 0;
            for (int k = 0; k < 6; k++) {
                C[6*i + j] += A[6*i + k] * B[6*k + j];
            }
        }
    }
}

void matadj6x6(int64_t M[36], int64_t Madj[36])
{
    int i = 0;
    int64_t C[36] = {0};
    C[0] = 1; C[7] = 1; C[14] = 1; C[21] = 1; C[28] = 1; C[35] = 1;
    int64_t a0 = 1;
    int64_t MC[36] = {0};
    int64_t ai;
    while (true) {
        i++;
        matmul6x6(MC, M, C);
        if (i == 6) {
            ai = -(MC[0] + MC[7] + MC[14] + MC[21] + MC[28] + MC[35])/6;
            for (int k = 0; k < 36; k++) Madj[k] = -C[k];
            break;
        }
        for (int k = 0; k < 36; k++) C[k] = MC[k];
        ai = -(C[0] + C[7] + C[14] + C[21] + C[28] + C[35])/i;
        C[0] += ai; C[7] += ai; C[14] += ai; C[21] += ai; C[28] += ai; C[35] += ai;
    }
}

int64_t subdet(int64_t M[36], int d, int j)
{
	int64_t subM[36] = {0};
	int sub_i = 0, sub_j = 0;
	for (int i = 1; i < d; i++) {
		for (int k = 0; k < d; k++) {
			if (k != j) {
				subM[6*sub_i + sub_j] = M[6*i + k];
				sub_j++;
			}
		}
		sub_i++;
		sub_j = 0;
	}
	return (matdet(subM, 5));
}

// Compute the determinant of a square matrix M of dimension at most 6
int64_t matdet(int64_t M[36], int d)
{
    int64_t det = 0;

    int k = 0; int64_t c = 1; int64_t s = 1;
    while (true) {
        k++;
        if (k==d) {
            det = s*M[7*d - 7];
            break;
        }
        int64_t p = M[7*k - 7];
        if (p == 0) {
            det = 0;
            for (int i = k+1; i <= d; i++) {
                det = M[6*(i-1) + k - 1];
                if (det != 0) {
                    for (int j = k; j <= d; j++) {
                        int64_t t = M[6*(i-1) + j - 1];
                        M[6*(i-1) + j - 1] = M[6*(k-1) + j - 1];
                        M[6*(k-1) + j - 1] = t;
                    }
                    s = -s;
                    p = M[7*k - 7];
                    break;
                }
            }
            if (det == 0) break;
        }
        for (int i = k+1; i <= d; i++) {
            for (int j = k+1; j <= d; j++) {
                int64_t t = p*M[6*(i-1) + j - 1] - M[6*(i-1) + k - 1]*M[6*(k-1) + j - 1];
                M[6*(i-1) + j - 1] = t/c;
            }
        }
        c = p;
    }
    return (det);
}


inline int64_t floordiv(int64_t a, int64_t b)
{
    int64_t d = a / b;
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
// within [-B,B[x[-B,B[x[-B,B[x[-B,B[
inline void getab4d(int u1, int u2, int u3, int u4, int v1, int v2, int v3, int v4,
				int x, int y, int z, int t, int B1, int B2, int B3, int B4, int* a, int* b)
{
	int U[8] = { u1, u2, u3, u4, -u1, -u2, -u3, -u4 };
	int V[8] = { v1, v2, v3, v4, -v1, -v2, -v3, -v4 };
	int C[8] = { B1-x-1, B2-y-1, B3-z-1, B4-t-1, B1+x, B2+y, B3+z, B4+t };
	*a = B1; *b = B1;

	int64_t L = abs(nonzerolcm4d(u1, u2, u3, u4));

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


// Integer programming to get (a,b) such that
// (x,y,z,s,t) + a*(u1,u2,u3,u4,u5) + b*(v1,v2,3v,v4,v5)
// within [-B,B[x[-B,B[x[-B,B[x[-B,B[x[-B,B[
inline void getab5d(int u1, int u2, int u3, int u4, int u5, int v1, int v2, int v3, int v4,
				int v5, int x, int y, int z, int s, int t, int B1, int B2, int B3, int B4,
				int B5, int* a, int* b)
{
	int U[10] = { u1, u2, u3, u4, u5, -u1, -u2, -u3, -u4, -u5 };
	int V[10] = { v1, v2, v3, v4, v5, -v1, -v2, -v3, -v4, -v5 };
	int C[10] = { B1-x-1, B2-y-1, B3-z-1, B4-s-1, B5-t-1, B1+x, B2+y, B3+z, B4+s, B5+t };
	*a = B1; *b = B1;

	int64_t L = abs(nonzerolcm5d(u1, u2, u3, u4, u5));

	for (int i = 0; i < 10; i++) {
		int s = abs(U[i]);
		if (s != 0) {
			V[i] *= L / s;
			C[i] *= L / s;
		}
	}
	
	for (int i = 0; i < 10; i++) {
		if (U[i] == 0) {
			if (V[i] > 0) {
				int bnew = floordiv(C[i], V[i]);
				if (bnew < *b) *b = bnew;
			}
		}
		else if (U[i] < 0) {
			for (int j = 0; j < 10; j++) {
				if (U[j] > 0 && abs(i-j) != 5) {
					int D = V[i] + V[j];
					if (D > 0) {
						int bnew = floordiv(C[i] + C[j], D);
						if (bnew < *b) *b = bnew;
					}
				}
			}
		}
	}

	for (int i = 0; i < 10; i++) {
		if (U[i] > 0) {
			int anew = floordiv(C[i] - V[i] * (*b), L);
			if (anew < *a) *a = anew;
		}
	}
}


// Integer programming to get (a,b) such that
// (x,y,z,r,s,t) + a*(u1,u2,u3,u4,u5,u6) + b*(v1,v2,3v,v4,v5,v6)
// within [-B,B[x[-B,B[x[-B,B[x[-B,B[x[-B,B[x[-B,B[
inline void getab6d(int u1, int u2, int u3, int u4, int u5, int u6, int v1, int v2, int v3,
				int v4,	int v5, int v6, int x, int y, int z, int r, int s, int t,
				int B1, int B2, int B3, int B4, int B5, int B6, int* a, int* b)
{
	int U[12] = { u1, u2, u3, u4, u5, u6, -u1, -u2, -u3, -u4, -u5, -u6 };
	int V[12] = { v1, v2, v3, v4, v5, v6, -v1, -v2, -v3, -v4, -v5, -v6 };
	int C[12] = { B1-x-1, B2-y-1, B3-z-1, B4-r-1, B5-s-1, B6-t-1,
					B1+x, B2+y, B3+z, B4+4, B5+s, B6+t };
	*a = B1; *b = B1;

	int64_t L = abs(nonzerolcm6d(u1, u2, u3, u4, u5, u6));

	for (int i = 0; i < 12; i++) {
		int s = abs(U[i]);
		if (s != 0) {
			V[i] *= L / s;
			C[i] *= L / s;
		}
	}
	
	for (int i = 0; i < 12; i++) {
		if (U[i] == 0) {
			if (V[i] > 0) {
				int bnew = floordiv(C[i], V[i]);
				if (bnew < *b) *b = bnew;
			}
		}
		else if (U[i] < 0) {
			for (int j = 0; j < 12; j++) {
				if (U[j] > 0 && abs(i-j) != 6) {
					int D = V[i] + V[j];
					if (D > 0) {
						int bnew = floordiv(C[i] + C[j], D);
						if (bnew < *b) *b = bnew;
					}
				}
			}
		}
	}

	for (int i = 0; i < 12; i++) {
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

inline int64_t nonzerolcm4d(int u1, int u2, int u3, int u4)
{
    int64_t temp = lcm(u1, u2);
    temp = lcm(temp, u3);
    temp = lcm(temp, u4);
    return abs(temp);
}

// Function to compute the least common multiple of two integers
int lcm(int a, int b) {
    return (a * b) / gcd(a, b);
}

// Function to compute the least common multiple of five integers
inline int64_t nonzerolcm5d(int u1, int u2, int u3, int u4, int u5) {
    int64_t temp = lcm(u1, u2);
    temp = lcm(temp, u3);
    temp = lcm(temp, u4);
    temp = lcm(temp, u5);
    return abs(temp);
}

// Function to compute the least common multiple of six integers
inline int64_t nonzerolcm6d(int u1, int u2, int u3, int u4, int u5, int u6) {
    int64_t temp = lcm(u1, u2);
    temp = lcm(temp, u3);
    temp = lcm(temp, u4);
    temp = lcm(temp, u5);
    temp = lcm(temp, u6);
    return abs(temp);
}

// determine if 3-space with normal (nx, ny, nz, nr, ns, nt) containing point
// (x,y,z,r,s,t) intersects box [-B,B[x[-B,B[x[-B,B[x[-B,B[x[-B,B[x[-B,B[
inline bool threespaceintersectsbox6d(int64_t nx, int64_t ny, int64_t nz, int64_t nr, 
	int64_t ns, int64_t nt, int x, int y, int z, int r, int s, int t, int B1, int B2,
	int B3, int B4, int B5, int B6)
{
	int64_t d = nx*x + ny*y + nz*z + nr*r + ns*s + nt*t;

	int64_t nxB0 = -nx*B1; int64_t nxB1 = nx*(B1-1);
	int64_t nyB0 = -ny*B2; int64_t nyB1 = ny*(B2-1);
	int64_t nzB0 = -nz*B3; int64_t nzB1 = nz*(B3-1);
	int64_t nrB0 = -nr*B4; int64_t nrB1 = nr*(B4-1);
	int64_t nsB0 = -ns*B5; int64_t nsB1 = ns*(B5-1);
	int64_t ntB0 = -nt*B6; int64_t ntB1 = nt*(B6-1);

	int s0 = ( nxB0 + nyB0 + nzB0 + nrB0 + nsB0 + ntB0 - d ) > 0;
	int s1 = ( nxB1 + nyB0 + nzB0 + nrB0 + nsB0 + ntB0 - d ) > 0;
	int s2 = ( nxB0 + nyB1 + nzB0 + nrB0 + nsB0 + ntB0 - d ) > 0;
	int s3 = ( nxB1 + nyB1 + nzB0 + nrB0 + nsB0 + ntB0 - d ) > 0;
	int s4 = ( nxB0 + nyB0 + nzB1 + nrB0 + nsB0 + ntB0 - d ) > 0;
	int s5 = ( nxB1 + nyB0 + nzB1 + nrB0 + nsB0 + ntB0 - d ) > 0;
	int s6 = ( nxB0 + nyB1 + nzB1 + nrB0 + nsB0 + ntB0 - d ) > 0;
	int s7 = ( nxB1 + nyB1 + nzB1 + nrB0 + nsB0 + ntB0 - d ) > 0;
	int s8 = ( nxB0 + nyB0 + nzB0 + nrB1 + nsB0 + ntB0 - d ) > 0;
	int s9 = ( nxB1 + nyB0 + nzB0 + nrB1 + nsB0 + ntB0 - d ) > 0;
	int sa = ( nxB0 + nyB1 + nzB0 + nrB1 + nsB0 + ntB0 - d ) > 0;
	int sb = ( nxB1 + nyB1 + nzB0 + nrB1 + nsB0 + ntB0 - d ) > 0;
	int sc = ( nxB0 + nyB0 + nzB1 + nrB1 + nsB0 + ntB0 - d ) > 0;
	int sd = ( nxB1 + nyB0 + nzB1 + nrB1 + nsB0 + ntB0 - d ) > 0;
	int se = ( nxB0 + nyB1 + nzB1 + nrB1 + nsB0 + ntB0 - d ) > 0;
	int sf = ( nxB1 + nyB1 + nzB1 + nrB1 + nsB0 + ntB0 - d ) > 0;

	int mask = s0 + (s1<<1) + (s2<<2) + (s3<<3) + (s4<<4) + (s5<<5) + (s6<<6) + (s7<<7) +
	           (s8<<8) + (s9<<9) + (sa<<10) + (sb<<11) +
			   (sc<<12) + (sd<<13) + (se<<14) + (sf<<15);

	return ( (mask != 0) && (mask != 65535) );
}


// determine if 4-space with normal (nx, ny, nz, ns, nt) containing point (x,y,z,s,t)
// intersects box [-B,B[x[-B,B[x[-B,B[x[-B,B[x[-B,B[
inline bool fourspaceintersectsbox6d(int64_t nx, int64_t ny, int64_t nz, int64_t nr, 
	int64_t ns, int64_t nt, int x, int y, int z, int r, int s, int t, int B1, int B2,
	int B3, int B4, int B5, int B6)
{
	int64_t d = nx*x + ny*y + nz*z + nr*r + ns*s + nt*t;

	long nxB0 = -nx*B1; long nxB1 = nx*(B1-1);
	long nyB0 = -ny*B2; long nyB1 = ny*(B2-1);
	long nzB0 = -nz*B3; long nzB1 = nz*(B3-1);
	long nrB0 = -nr*B4; long nrB1 = nr*(B4-1);
	long nsB0 = -ns*B5; long nsB1 = ns*(B5-1);
	long ntB0 = -nt*B6; long ntB1 = nt*(B6-1);
	
	long s00 = ( nxB0 + nyB0 + nzB0 + nrB0 + nsB0 + ntB0 - d ) > 0;
	long s01 = ( nxB1 + nyB0 + nzB0 + nrB0 + nsB0 + ntB0 - d ) > 0;
	long s02 = ( nxB0 + nyB1 + nzB0 + nrB0 + nsB0 + ntB0 - d ) > 0;
	long s03 = ( nxB1 + nyB1 + nzB0 + nrB0 + nsB0 + ntB0 - d ) > 0;
	long s04 = ( nxB0 + nyB0 + nzB1 + nrB0 + nsB0 + ntB0 - d ) > 0;
	long s05 = ( nxB1 + nyB0 + nzB1 + nrB0 + nsB0 + ntB0 - d ) > 0;
	long s06 = ( nxB0 + nyB1 + nzB1 + nrB0 + nsB0 + ntB0 - d ) > 0;
	long s07 = ( nxB1 + nyB1 + nzB1 + nrB0 + nsB0 + ntB0 - d ) > 0;
	long s08 = ( nxB0 + nyB0 + nzB0 + nrB1 + nsB0 + ntB0 - d ) > 0;
	long s09 = ( nxB1 + nyB0 + nzB0 + nrB1 + nsB0 + ntB0 - d ) > 0;
	long s0a = ( nxB0 + nyB1 + nzB0 + nrB1 + nsB0 + ntB0 - d ) > 0;
	long s0b = ( nxB1 + nyB1 + nzB0 + nrB1 + nsB0 + ntB0 - d ) > 0;
	long s0c = ( nxB0 + nyB0 + nzB1 + nrB1 + nsB0 + ntB0 - d ) > 0;
	long s0d = ( nxB1 + nyB0 + nzB1 + nrB1 + nsB0 + ntB0 - d ) > 0;
	long s0e = ( nxB0 + nyB1 + nzB1 + nrB1 + nsB0 + ntB0 - d ) > 0;
	long s0f = ( nxB1 + nyB1 + nzB1 + nrB1 + nsB0 + ntB0 - d ) > 0;
	long s10 = ( nxB0 + nyB0 + nzB0 + nrB0 + nsB1 + ntB0 - d ) > 0;
	long s11 = ( nxB1 + nyB0 + nzB0 + nrB0 + nsB1 + ntB0 - d ) > 0;
	long s12 = ( nxB0 + nyB1 + nzB0 + nrB0 + nsB1 + ntB0 - d ) > 0;
	long s13 = ( nxB1 + nyB1 + nzB0 + nrB0 + nsB1 + ntB0 - d ) > 0;
	long s14 = ( nxB0 + nyB0 + nzB1 + nrB0 + nsB1 + ntB0 - d ) > 0;
	long s15 = ( nxB1 + nyB0 + nzB1 + nrB0 + nsB1 + ntB0 - d ) > 0;
	long s16 = ( nxB0 + nyB1 + nzB1 + nrB0 + nsB1 + ntB0 - d ) > 0;
	long s17 = ( nxB1 + nyB1 + nzB1 + nrB0 + nsB1 + ntB0 - d ) > 0;
	long s18 = ( nxB0 + nyB0 + nzB0 + nrB1 + nsB1 + ntB0 - d ) > 0;
	long s19 = ( nxB1 + nyB0 + nzB0 + nrB1 + nsB1 + ntB0 - d ) > 0;
	long s1a = ( nxB0 + nyB1 + nzB0 + nrB1 + nsB1 + ntB0 - d ) > 0;
	long s1b = ( nxB1 + nyB1 + nzB0 + nrB1 + nsB1 + ntB0 - d ) > 0;
	long s1c = ( nxB0 + nyB0 + nzB1 + nrB1 + nsB1 + ntB0 - d ) > 0;
	long s1d = ( nxB1 + nyB0 + nzB1 + nrB1 + nsB1 + ntB0 - d ) > 0;
	long s1e = ( nxB0 + nyB1 + nzB1 + nrB1 + nsB1 + ntB0 - d ) > 0;
	long s1f = ( nxB1 + nyB1 + nzB1 + nrB1 + nsB1 + ntB0 - d ) > 0;

	long mask = (s00)     + (s01<<1)  + (s02<<2)  + (s03<<3)  + (s04<<4)  + (s05<<5)  +
	            (s06<<6)  + (s07<<7)  + (s08<<8)  + (s09<<9)  + (s0a<<10) + (s0b<<11) +
				(s0c<<12) + (s0d<<13) + (s0e<<14) + (s0f<<15) +
			    (s10<<16) + (s11<<17) + (s12<<18) + (s13<<19) + (s14<<20) + (s15<<21) +
			    (s16<<22) + (s17<<23) + (s18<<24) + (s19<<25) + (s1a<<26) + (s1b<<27) +
			    (s1c<<28) + (s1d<<29) + (s1e<<30) + (s1f<<31);

	return ( (mask != 0) && (mask != 4294967295L) );
}


// determine if 5-space with normal (nx, ny, nz, nr, ns, nt) containing point (x,y,z,r,s,t)
// intersects box [-B,B[x[-B,B[x[-B,B[x[-B,B[x[-B,B[x[-B,B[
inline bool fivespaceintersectsbox6d(int64_t nx, int64_t ny, int64_t nz, int64_t nr,
	int64_t ns, int64_t nt, int x, int y, int z, int r, int s, int t, int B1, int B2,
	int B3, int B4, int B5, int B6)
{
	int64_t d = nx*x + ny*y + nz*z + nr*r + ns*s + nt*t;

	int64_t nxB0 = -nx*B1; int64_t nxB1 = nx*(B1-1);
	int64_t nyB0 = -ny*B2; int64_t nyB1 = ny*(B2-1);
	int64_t nzB0 = -nz*B3; int64_t nzB1 = nz*(B3-1);
	int64_t nrB0 = -nr*B4; int64_t nrB1 = nr*(B4-1);
	int64_t nsB0 = -ns*B5; int64_t nsB1 = ns*(B5-1);
	int64_t ntB0 = -nt*B6; int64_t ntB1 = nt*(B6-1);

	uint64_t s00 = ( nxB0 + nyB0 + nzB0 + nrB0 + nsB0 + ntB0 - d ) > 0;
	uint64_t s01 = ( nxB1 + nyB0 + nzB0 + nrB0 + nsB0 + ntB0 - d ) > 0;
	uint64_t s02 = ( nxB0 + nyB1 + nzB0 + nrB0 + nsB0 + ntB0 - d ) > 0;
	uint64_t s03 = ( nxB1 + nyB1 + nzB0 + nrB0 + nsB0 + ntB0 - d ) > 0;
	uint64_t s04 = ( nxB0 + nyB0 + nzB1 + nrB0 + nsB0 + ntB0 - d ) > 0;
	uint64_t s05 = ( nxB1 + nyB0 + nzB1 + nrB0 + nsB0 + ntB0 - d ) > 0;
	uint64_t s06 = ( nxB0 + nyB1 + nzB1 + nrB0 + nsB0 + ntB0 - d ) > 0;
	uint64_t s07 = ( nxB1 + nyB1 + nzB1 + nrB0 + nsB0 + ntB0 - d ) > 0;
	uint64_t s08 = ( nxB0 + nyB0 + nzB0 + nrB1 + nsB0 + ntB0 - d ) > 0;
	uint64_t s09 = ( nxB1 + nyB0 + nzB0 + nrB1 + nsB0 + ntB0 - d ) > 0;
	uint64_t s0a = ( nxB0 + nyB1 + nzB0 + nrB1 + nsB0 + ntB0 - d ) > 0;
	uint64_t s0b = ( nxB1 + nyB1 + nzB0 + nrB1 + nsB0 + ntB0 - d ) > 0;
	uint64_t s0c = ( nxB0 + nyB0 + nzB1 + nrB1 + nsB0 + ntB0 - d ) > 0;
	uint64_t s0d = ( nxB1 + nyB0 + nzB1 + nrB1 + nsB0 + ntB0 - d ) > 0;
	uint64_t s0e = ( nxB0 + nyB1 + nzB1 + nrB1 + nsB0 + ntB0 - d ) > 0;
	uint64_t s0f = ( nxB1 + nyB1 + nzB1 + nrB1 + nsB0 + ntB0 - d ) > 0;
	uint64_t s10 = ( nxB0 + nyB0 + nzB0 + nrB0 + nsB1 + ntB0 - d ) > 0;
	uint64_t s11 = ( nxB1 + nyB0 + nzB0 + nrB0 + nsB1 + ntB0 - d ) > 0;
	uint64_t s12 = ( nxB0 + nyB1 + nzB0 + nrB0 + nsB1 + ntB0 - d ) > 0;
	uint64_t s13 = ( nxB1 + nyB1 + nzB0 + nrB0 + nsB1 + ntB0 - d ) > 0;
	uint64_t s14 = ( nxB0 + nyB0 + nzB1 + nrB0 + nsB1 + ntB0 - d ) > 0;
	uint64_t s15 = ( nxB1 + nyB0 + nzB1 + nrB0 + nsB1 + ntB0 - d ) > 0;
	uint64_t s16 = ( nxB0 + nyB1 + nzB1 + nrB0 + nsB1 + ntB0 - d ) > 0;
	uint64_t s17 = ( nxB1 + nyB1 + nzB1 + nrB0 + nsB1 + ntB0 - d ) > 0;
	uint64_t s18 = ( nxB0 + nyB0 + nzB0 + nrB1 + nsB1 + ntB0 - d ) > 0;
	uint64_t s19 = ( nxB1 + nyB0 + nzB0 + nrB1 + nsB1 + ntB0 - d ) > 0;
	uint64_t s1a = ( nxB0 + nyB1 + nzB0 + nrB1 + nsB1 + ntB0 - d ) > 0;
	uint64_t s1b = ( nxB1 + nyB1 + nzB0 + nrB1 + nsB1 + ntB0 - d ) > 0;
	uint64_t s1c = ( nxB0 + nyB0 + nzB1 + nrB1 + nsB1 + ntB0 - d ) > 0;
	uint64_t s1d = ( nxB1 + nyB0 + nzB1 + nrB1 + nsB1 + ntB0 - d ) > 0;
	uint64_t s1e = ( nxB0 + nyB1 + nzB1 + nrB1 + nsB1 + ntB0 - d ) > 0;
	uint64_t s1f = ( nxB1 + nyB1 + nzB1 + nrB1 + nsB1 + ntB0 - d ) > 0;
	uint64_t s20 = ( nxB0 + nyB0 + nzB0 + nrB0 + nsB0 + ntB1 - d ) > 0;
	uint64_t s21 = ( nxB1 + nyB0 + nzB0 + nrB0 + nsB0 + ntB1 - d ) > 0;
	uint64_t s22 = ( nxB0 + nyB1 + nzB0 + nrB0 + nsB0 + ntB1 - d ) > 0;
	uint64_t s23 = ( nxB1 + nyB1 + nzB0 + nrB0 + nsB0 + ntB1 - d ) > 0;
	uint64_t s24 = ( nxB0 + nyB0 + nzB1 + nrB0 + nsB0 + ntB1 - d ) > 0;
	uint64_t s25 = ( nxB1 + nyB0 + nzB1 + nrB0 + nsB0 + ntB1 - d ) > 0;
	uint64_t s26 = ( nxB0 + nyB1 + nzB1 + nrB0 + nsB0 + ntB1 - d ) > 0;
	uint64_t s27 = ( nxB1 + nyB1 + nzB1 + nrB0 + nsB0 + ntB1 - d ) > 0;
	uint64_t s28 = ( nxB0 + nyB0 + nzB0 + nrB1 + nsB0 + ntB1 - d ) > 0;
	uint64_t s29 = ( nxB1 + nyB0 + nzB0 + nrB1 + nsB0 + ntB1 - d ) > 0;
	uint64_t s2a = ( nxB0 + nyB1 + nzB0 + nrB1 + nsB0 + ntB1 - d ) > 0;
	uint64_t s2b = ( nxB1 + nyB1 + nzB0 + nrB1 + nsB0 + ntB1 - d ) > 0;
	uint64_t s2c = ( nxB0 + nyB0 + nzB1 + nrB1 + nsB0 + ntB1 - d ) > 0;
	uint64_t s2d = ( nxB1 + nyB0 + nzB1 + nrB1 + nsB0 + ntB1 - d ) > 0;
	uint64_t s2e = ( nxB0 + nyB1 + nzB1 + nrB1 + nsB0 + ntB1 - d ) > 0;
	uint64_t s2f = ( nxB1 + nyB1 + nzB1 + nrB1 + nsB0 + ntB1 - d ) > 0;
	uint64_t s30 = ( nxB0 + nyB0 + nzB0 + nrB0 + nsB1 + ntB1 - d ) > 0;
	uint64_t s31 = ( nxB1 + nyB0 + nzB0 + nrB0 + nsB1 + ntB1 - d ) > 0;
	uint64_t s32 = ( nxB0 + nyB1 + nzB0 + nrB0 + nsB1 + ntB1 - d ) > 0;
	uint64_t s33 = ( nxB1 + nyB1 + nzB0 + nrB0 + nsB1 + ntB1 - d ) > 0;
	uint64_t s34 = ( nxB0 + nyB0 + nzB1 + nrB0 + nsB1 + ntB1 - d ) > 0;
	uint64_t s35 = ( nxB1 + nyB0 + nzB1 + nrB0 + nsB1 + ntB1 - d ) > 0;
	uint64_t s36 = ( nxB0 + nyB1 + nzB1 + nrB0 + nsB1 + ntB1 - d ) > 0;
	uint64_t s37 = ( nxB1 + nyB1 + nzB1 + nrB0 + nsB1 + ntB1 - d ) > 0;
	uint64_t s38 = ( nxB0 + nyB0 + nzB0 + nrB1 + nsB1 + ntB1 - d ) > 0;
	uint64_t s39 = ( nxB1 + nyB0 + nzB0 + nrB1 + nsB1 + ntB1 - d ) > 0;
	uint64_t s3a = ( nxB0 + nyB1 + nzB0 + nrB1 + nsB1 + ntB1 - d ) > 0;
	uint64_t s3b = ( nxB1 + nyB1 + nzB0 + nrB1 + nsB1 + ntB1 - d ) > 0;
	uint64_t s3c = ( nxB0 + nyB0 + nzB1 + nrB1 + nsB1 + ntB1 - d ) > 0;
	uint64_t s3d = ( nxB1 + nyB0 + nzB1 + nrB1 + nsB1 + ntB1 - d ) > 0;
	uint64_t s3e = ( nxB0 + nyB1 + nzB1 + nrB1 + nsB1 + ntB1 - d ) > 0;
	uint64_t s3f = ( nxB1 + nyB1 + nzB1 + nrB1 + nsB1 + ntB1 - d ) > 0;

	uint64_t mask = (s00)     + (s01<<1)  + (s02<<2)  + (s03<<3)  + (s04<<4)  + (s05<<5)  +
	                (s06<<6)  + (s07<<7)  + (s08<<8)  + (s09<<9)  + (s0a<<10) + (s0b<<11) +
				    (s0c<<12) + (s0d<<13) + (s0e<<14) + (s0f<<15) +
			        (s10<<16) + (s11<<17) + (s12<<18) + (s13<<19) + (s14<<20) + (s15<<21) +
			        (s16<<22) + (s17<<23) + (s18<<24) + (s19<<25) + (s1a<<26) + (s1b<<27) +
			        (s1c<<28) + (s1d<<29) + (s1e<<30) + (s1f<<31) +
					(s20<<32) + (s21<<33) + (s22<<34) + (s23<<35) + (s24<<36) + (s25<<37) +
	                (s26<<38) + (s27<<39) + (s28<<40) + (s29<<41) + (s2a<<42) + (s2b<<43) +
				    (s2c<<44) + (s2d<<45) + (s2e<<46) + (s2f<<47) +
			        (s30<<48) + (s31<<49) + (s32<<50) + (s33<<51) + (s34<<52) + (s35<<53) +
			        (s36<<54) + (s37<<55) + (s38<<56) + (s39<<57) + (s3a<<58) + (s3b<<59) +
			        (s3c<<60) + (s3d<<61) + (s3e<<62) + (s3f<<63);

	return ( (mask != 0) && (mask != 18446744073709551615UL) );
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


inline int minnonneg6d(int u, int v, int w, int r, int s, int t)
{
	int m = u;
	if (m < 0) m = v; if (m < 0) m = w; if (m < 0) m = r; if (m < 0) m = s; if (m < 0) m = t;
	if (v >= 0 && v < m) m = v;
	if (w >= 0 && w < m) m = w;
	if (r >= 0 && r < m) m = r;
	if (s >= 0 && s < m) m = s;
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

int populate_q(sievedata* info, int side, mpz_poly h0, mpz_poly f0, mpz_poly f1,
	mpz_poly_bivariate F0, mpz_poly_bivariate F1)
{
	gmp_randstate_t rstate;
	gmp_randinit_mt(rstate);
    gmp_randseed_ui(rstate, 123);

	mpz_poly_factor_list lf;
	mpz_poly_factor_list_init(lf);
	mpz_poly_factor_list lh;
	mpz_poly_factor_list_init(lh);
	mpz_poly_factor_list l3;
	mpz_poly_factor_list_init(l3);
	mpz_poly hlin; mpz_poly_init(hlin, 0);
	mpz_poly_bivariate Hlin; mpz_poly_bivariate_init(Hlin, 0);
	mpz_poly fhlin; mpz_poly_init(fhlin, 0);
	mpz_poly_bivariate* factors0 = new mpz_poly_bivariate[f0->deg / 2];
	mpz_poly_bivariate* factors1 = new mpz_poly_bivariate[f1->deg / 2];
	for (int i = 0; i < f0->deg / 2; i++) mpz_poly_bivariate_init(factors0[i], 0);
	for (int i = 0; i < f1->deg / 2; i++) mpz_poly_bivariate_init(factors1[i], 0);
	// compute discriminants of f0 and f1
	mpz_t Df0; mpz_init(Df0); mpz_t Df1; mpz_init(Df1);
	mpz_t dummy; mpz_init(dummy);
	mpz_poly_discriminant(Df0, f0);
	mpz_poly_discriminant(Df1, f1);
	bool allh1, allh3, allf1, allf2, allf3, allf6;
	int n = 0;	// number of special-q's found

	int64_t q0 = info->q;
	//cout << q0 << endl;
	mpz_t q; mpz_init_set_ui(q, q0);
	if (side == 0) {
		// skip q0 if it is ramified in K[x]/<f0>
		if (mpz_mod_ui(dummy, Df0, q0) == 0) return 0;

		// determine which type of special-q we have, first factor f mod q
		mpz_poly_factor(lf, f0, q, rstate);
		allf6 = true;
		for (int j = 0; j < lf->size; j++)
			if (lf->factors[j]->f->deg != 6) {
				allf6 = false;
				break;
			}
		allf3 = true;
		if (allf6) allf3 = false;
		for (int j = 0; j < lf->size; j++)
			if (lf->factors[j]->f->deg != 3) {
				allf3 = false;
				break;
			}
		
		allf2 = true;
		if (allf6 || allf3) allf2 = false;
		else
			for (int j = 0; j < lf->size; j++)
				if (lf->factors[j]->f->deg != 2) {
					allf2 = false;
					break;
				}
		allf1 = true;
		if (allf6 || allf3 || allf2) allf1 = false;
		else
			for (int j = 0; j < lf->size; j++)
				if (lf->factors[j]->f->deg != 1) {
					allf1 = false;
					break;
				}
		mpz_poly_factor(lh, h0, q, rstate);
		allh3 = true;
		for (int j = 0; j < lh->size; j++)
			if (lh->factors[j]->f->deg != 3) {
				allh3 = false;
				break;
			}
		allh1 = true;
		if (allh3) allh1 = false;
		for (int j = 0; j < lh->size; j++)
			if (lh->factors[j]->f->deg != 1) {
				allh1 = false;
				break;
			}

		// categorize special-q by allf4, allf2, allf1, allh1
		if (!allf6 && !allf3 && !allf2 && !allh3) {
			for (int j = 0; j < lf->size; j++) { // lf sorted by deg increasing
				if (lf->factors[j]->f->deg == 1) {
					int64_t R = mpz_get_ui(lf->factors[j]->f->coeff[0]);
					// find corresponding root r of h
					for (int k = 0; k < lh->size; k++) {
						int64_t r = mpz_get_ui(lh->factors[k]->f->coeff[0]);
						mpz_poly_setcoeff_ui(hlin, 0, r);
						mpz_poly_setcoeff_ui(hlin, 1, 1);
						mpz_poly_bivariate_setcoeff(Hlin, 0, hlin);
						mpz_poly_set_zero(fhlin);
						mpz_poly_bivariate_resultant_x(fhlin, F0, Hlin);
						// now factor fhlin mod q
						mpz_poly_factor(l3, fhlin, q, rstate);
						bool foundR = false;
						for (int l = 0; l < l3->size; l++) {
							int64_t Rtest = mpz_get_ui(l3->factors[l]->f->coeff[0]);
							if (R == Rtest) {
								foundR = true;
								// now we have R, r
								info->qtype[n] = 2;
								info->r[n] = r;
								info->R[n++] = R;
								break;
							}
						}
						mpz_poly_factor_list_flush(l3);
						if (foundR) break;
					}
				}
			}
		}
	}
	else {	// side == 1
		// skip q0 if it is ramified in K[x]/<f1>
		if (mpz_mod_ui(dummy, Df1, q0) == 0) return 0;

		// determine which type of special-q we have, first factor f mod q
		mpz_poly_factor(lf, f1, q, rstate);
		allf6 = true;
		for (int j = 0; j < lf->size; j++)
			if (lf->factors[j]->f->deg != 6) {
				allf6 = false;
				break;
			}
		allf3 = true;
		if (allf6) allf3 = false;
		for (int j = 0; j < lf->size; j++)
			if (lf->factors[j]->f->deg != 3) {
				allf3 = false;
				break;
			}
		
		allf2 = true;
		if (allf6 || allf3) allf2 = false;
		else
			for (int j = 0; j < lf->size; j++)
				if (lf->factors[j]->f->deg != 2) {
					allf2 = false;
					break;
				}
		allf1 = true;
		if (allf6 || allf3 || allf2) allf1 = false;
		else
			for (int j = 0; j < lf->size; j++)
				if (lf->factors[j]->f->deg != 1) {
					allf1 = false;
					break;
				}
		mpz_poly_factor(lh, h0, q, rstate);
		allh3 = true;
		for (int j = 0; j < lh->size; j++)
			if (lh->factors[j]->f->deg != 3) {
				allh3 = false;
				break;
			}
		allh1 = true;
		if (allh3) allh1 = false;
		for (int j = 0; j < lh->size; j++)
			if (lh->factors[j]->f->deg != 1) {
				allh1 = false;
				break;
			}

		// categorize special-q by allf4, allf2, allf1, allh1
		if (!allf6 && !allf3 && !allf2 && !allh3) {
			for (int j = 0; j < lf->size; j++) { // lf sorted by deg increasing
				if (lf->factors[j]->f->deg == 1) {
					int64_t R = mpz_get_ui(lf->factors[j]->f->coeff[0]);
					// find corresponding root r of h
					for (int k = 0; k < lh->size; k++) {
						int64_t r = mpz_get_ui(lh->factors[k]->f->coeff[0]);
						mpz_poly_setcoeff_ui(hlin, 0, r);
						mpz_poly_setcoeff_ui(hlin, 1, 1);
						mpz_poly_bivariate_setcoeff(Hlin, 0, hlin);
						mpz_poly_set_zero(fhlin);
						mpz_poly_bivariate_resultant_x(fhlin, F1, Hlin);
						// now factor fhlin mod q
						mpz_poly_factor(l3, fhlin, q, rstate);
						bool foundR = false;
						for (int l = 0; l < l3->size; l++) {
							int64_t Rtest = mpz_get_ui(l3->factors[l]->f->coeff[0]);
							if (R == Rtest) {
								foundR = true;
								// now we have R, r
								info->qtype[n] = 2;
								info->r[n] = r;
								info->R[n++] = R;
								break;
							}
						}
						mpz_poly_factor_list_flush(l3);
						if (foundR) break;
					}
				}
			}
		}
	}

	return n;
}
