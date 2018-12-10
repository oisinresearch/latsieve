#include <stdint.h>	// int64_t
#include <iostream> // cout
#include <iomanip> // setprecision
#include "L2lu64.h"
#include <gmpxx.h>
#include "intpoly.h"
#include <math.h>	// sqrt
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

struct keyval {
	int id;
	float logp;
};

int latsieve3d(int* f, int degf, int64_t q, int l, int* allp, int nump, int* s, int* num_smodp, keyval* M, int Mlen, int B);
void histogram(keyval*M, float* H, int len);
bool lattice_sorter(keyval const& kv1, keyval const& kv2);
void csort(keyval* M, keyval* L, int* H, int len);
inline int floordiv(int a, int b);
inline int64_t modinv(int64_t x, int64_t m);
inline int nonzerolcm(int u1, int u2, int u3);
inline int gcd(int a, int b);
inline void getab(int u1, int u2, int u3, int v1, int v2, int v3, int x, int y, int z, int B, int* a, int* b);
inline bool planeintersectsbox(int nx, int ny, int nz, int x, int y, int z, int B);
inline int minnonneg(int u, int v, int w);
inline int minabs(int u, int v, int w);
inline int maxabs(int u, int v, int w);
inline int min(int u, int v, int w);
inline int max(int u, int v, int w);
void GetlcmScalar(int B, mpz_t S, int* primes, int nump);
bool PollardPm1(mpz_t N, mpz_t S, mpz_t factor);
bool EECM(mpz_t N, mpz_t S, mpz_t factor, int d, int a, int X0, int Y0, int Z0);


int main(int argc, char** argv)
{
	if (argc == 1) {
		cout << endl << "Usage: ./latsieve3d inputpoly fbbits factorbasefile qmin qmax th0 th1 lpbbits" << endl << endl;
		return 0;
	}

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
	//mpz_clear(c);
	if (verbose) cout << endl << "Complete.  Degree f0 = " << degf << ", degree f1 = " << degg << "." << endl;

	if (verbose) cout << endl << "Starting sieve of Eratosthenes for small primes..." << endl << flush;
	int fbbits = 21;
	if (argc >=3) fbbits = atoi(argv[2]);
	int max = 1<<fbbits; // 10000000;// 65536;
	char* sieve = new char[max+1]();
	int* primes = new int[1077871]; // int[155611]; //new int[809228];	//new int[6542]; 	// 2039 is the 309th prime, largest below 2048
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
	if (verbose) cout << "Complete.  Time taken: " << timetaken << "s" << endl << flush;
	if (verbose) cout << "There are " << k0 << " factor base primes on side 0." << endl << flush;
	if (verbose) cout << "There are " << k1 << " factor base primes on side 1." << endl << flush;
	fbfile.close();
	
	int B = 256;	// 512;
	int Mlen = 512*512*256*2 + (1<<27);	//1024*1024*512*2 + (1<<27);
	keyval* M = new keyval[Mlen];	// lattice { id, logp } pairs
	//keyval* L = new keyval[Mlen];	// copy of M
	float* H = new float[Mlen];	// histogram
	vector<int> rel;
	// clear M
	//cout << "Clearing memory..." << endl << flush;
	//start = clock();
	//for (int j = 0; j < Mlen; j++) M[j] = (keyval){ 0, 0 };
	//memset(M, 0, Mlen * sizeof(keyval));
	//timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
	//cout << "Memory cleared. Time taken: " << timetaken << "s" << endl << flush;

	int qmin; int qmax; mpz_t qmpz; mpz_init(qmpz);
	mpz_t* pi = new mpz_t[8]; for (int i = 0; i < 8; i++) mpz_init(pi[i]);
	int* r = new int[degf]();
	int* fq = new int[degf+1]();
	mpz_poly f0; mpz_poly f1; mpz_poly A;
	mpz_poly_init(f0, degf); mpz_poly_init(f1, degg); mpz_poly_init(A, 3);
	mpz_poly_set_mpz(f0, fpoly, degf);
	mpz_poly_set_mpz(f1, gpoly, degg);
	mpz_t N0; mpz_t N1;
	mpz_init(N0); mpz_init(N1);
	stringstream stream;
	mpz_t lpb; mpz_init(lpb);
	mpz_t factor; mpz_init(factor); mpz_t p1; mpz_t p2; mpz_init(p1); mpz_init(p2); mpz_t t; mpz_init(t); 
	mpz_t S; mpz_init(S); GetlcmScalar(1000, S, primes, 669);	//1229);	// 900);// 669);	// for Pollard p-1
	if (argc >= 5) qmin = atoi(argv[4]);
	if (argc >= 6) qmax = atoi(argv[5]);
	float th0 = 70.0f;
	if (argc >= 7) th0 = atof(argv[6]);
	float th1 = 70.0f;
	if (argc >= 8) th1 = atof(argv[7]);
	int lpbits = 29;
	if (argc >= 9) lpbits = atoi(argv[8]);

	int q = qmin;
	while (q < qmax) {
		mpz_set_ui(qmpz, q);
		mpz_nextprime(qmpz, qmpz);
		q = mpz_get_ui(qmpz);
		for (int i = 0; i <= degf; i++) fq[i] = mpz_mod_ui(r0, fpoly[i], q);
		int numl = polrootsmod(fq, degf, r, q);
		if (numl == 0 || q > qmax) continue;
		
		// sieve side 0
		cout << "# Starting sieve on side 0 for special-q " << q << "..." << endl << flush;
		start = clock();
		int m = latsieve3d(fq, degf, q, 0, sievep0, k0, sieves0, sievenum_s0modp, M, Mlen, B);
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
		int B2 = 2*B; int BB = B*B;
		int B2bits = 1; while (1 << B2bits < B2) B2bits++;
		int BB4bits = 1; while (1 << BB4bits < BB*4) BB4bits++;
		rel.clear();
		for (int i = 0; i < m; i++) {
			if (H[i] > th0) {
				rel.push_back(i);
				R0++;
			}
		}
		cout << "# " << R0 << " candidates on side 0." << endl << flush;
		// sieve side 1
		cout << "# Starting sieve on side 1..." << endl << flush;
		start = clock();
		m = latsieve3d(fq, degf, q, 0, sievep1, k1, sieves1, sievenum_s1modp, M, Mlen, B);
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
		// print list of potential relations
		int R = 0;
		for (int i = 0; i < rel.size()-1; i++)
		{
			if (rel[i] == rel[i+1] && rel[i] != 0) {
				int x = (rel[i] % B2) - B;
				int y = ((rel[i] >> B2bits) % B2) - B;
				int z = rel[i] >> BB4bits;
				if (x != 0 || y != 0 || z != 0) {
					//cout << x << "," << y << "," << z << endl;
					R++;
				}
			}
		}
		cout << "# " << R << " potential relations found." << endl << flush;
	   
		// compute special-q lattice L 
		int64_t L[9];
		numl = polrootsmod(fq, degf, r, q);
		int l = 0;
		L[0] = q; L[1] = -r[l]; L[2] = 0;
		L[3] = 0; L[4] = 1;     L[5] = -r[l];
		L[6] = 0; L[7] = 0;     L[8] = 1;
		int64L2(L, 3);	// LLL reduce L, time log(q)^2
		
		// compute and factor resultants as much as possible, leave large gcd computation till later.
		start = clock();
		if (verbose) cout << "Computing Pollard p-1 scalar multiple..." << endl << flush;
		timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
		if (verbose) cout << "Finished! Time taken: " << timetaken << "s" << endl << flush;
		mpz_ui_pow_ui(lpb, 2, lpbits);
		int BASE = 16;
		int qside = 0;
		stack<mpz_t*> QN; stack<int> Q; int algarr[3]; mpz_t* N;
		start = clock();
		R = 0;
		if (verbose) cout << "Starting cofactorizaion..." << endl << flush;
		for (int i = 0; i < rel.size()-1; i++)
		{
			if (rel[i] == rel[i+1] && rel[i] != 0) {
				int x = (rel[i] % B2) - B;
				int y = ((rel[i] >> B2bits) % B2) - B;
				int z = rel[i] >> BB4bits;
				if (x != 0 || y != 0 || z != 0) {
					// compute [a,b,c]
					int a = L[0]*x+L[1]*y+L[2]*z;
					int b = L[3]*x+L[4]*y+L[5]*z;
					int c = L[6]*x+L[7]*y+L[8]*z;
	
					int64_t D = b*(int64_t)b - 4*a*(int64_t)c;
					if (floor(sqrt(D)+0.5)*floor(sqrt(D)+0.5) == D) {
						continue;	// a+b*x+c*x^2 is not irreducible over Z
					}

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
					//cout << mpz_get_str(NULL, 10, N0) << endl << flush;
					//cout << mpz_get_str(NULL, 10, N1) << endl << flush;
					string str = to_string(a) + "," + to_string(b) + "," + to_string(c) + ":";
					
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
					if (qside == 0) {
						mpz_divexact_ui(N0, N0, q);
						stream.str("");
						stream << hex << q;
						str += stream.str();
					}
					bool isrel = true;
					bool cofactor = true;
					if (mpz_cmp_ui(N0, 1) == 0) { cofactor = false; }
					str += (cofactor ? "," : "");
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
											 mpz_set(t, p1); mpz_set(p1, p2); mpz_set(p2, t);	// sort
										}
										// save remaining algs to array
										int lnext = l - j; int lt = lnext;
										while (lt--) { algarr[lt] = Q.top(); Q.pop(); }
										lt = lnext; if (lt) { while (lt--) Q.push(algarr[lnext-1-lt]); Q.push(lnext); }
										if (mpz_probab_prime_p(p1, 30)) {
											if (mpz_cmpabs(p1, lpb) > 0) { isrel = false; break; }
											else { str += mpz_get_str(NULL, BASE, p1); str += ","; }
										}
										else {
											if (!lnext) { isrel = false; break; }
											mpz_set(pi[n], p1);
											QN.push(&pi[n++]);
											lt = lnext; if (lt) { while (lt--) Q.push(algarr[lnext-1-lt]); Q.push(lnext); }
										}
										if (mpz_probab_prime_p(p2, 30)) {
											if (mpz_cmpabs(p2, lpb) > 0) { isrel = false; break; }
											else { str += mpz_get_str(NULL, BASE, p2); }
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
							else { str += mpz_get_str(NULL, BASE, N0); }
						}
					}
					
					str += ":";

					// trial division on side 1
					if (isrel) {
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
						if (qside == 1)  {
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
												 mpz_set(t, p1); mpz_set(p1, p2); mpz_set(p2, t);	// sort
											}
											// save remaining algs to array
											int lnext = l - j; int lt = lnext;
											while (lt--) { algarr[lt] = Q.top(); Q.pop(); }
											lt = lnext; if (lt) { while (lt--) Q.push(algarr[lnext-1-lt]); Q.push(lnext); }
											if (mpz_probab_prime_p(p1, 30)) {
												if (mpz_cmpabs(p1, lpb) > 0) { isrel = false; break; }
												else { str += mpz_get_str(NULL, BASE, p1); str += ","; }
											}
											else {
												if (!lnext) { isrel = false; break; }
												mpz_set(pi[n], p1);
												QN.push(&pi[n++]);
												lt = lnext; if (lt) { while (lt--) Q.push(algarr[lnext-1-lt]); Q.push(lnext); }
											}
											if (mpz_probab_prime_p(p2, 30)) {
												if (mpz_cmpabs(p2, lpb) > 0) { isrel = false; break; }
												else { str += mpz_get_str(NULL, BASE, p2); }
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
								else { str += mpz_get_str(NULL, BASE, N1); }
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

	mpz_clear(S);
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


void histogram(keyval*M, float* H, int len)
{
	// clear H
	memset(H, 0, len * sizeof(float));
	// fill H
	for (int i = 0; i < len; i++) H[M[i].id] += M[i].logp;
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


bool lattice_sorter(keyval const& kv1, keyval const& kv2)
{
	return kv2.id < kv1.id;
}


int latsieve3d(int* f, int degf, int64_t q, int l, int* allp, int nump, int* s, int* num_smodp, keyval* M, int Mlen, int B)
{
	int64_t L[9];
	int64_t L2[9];
	int64_t L3[9];
	int* r = new int[degf]();
	int numl = polrootsmod(f, degf, r, q);
	int B1 = B-1;
	int B2 = 2*B; int BB = B*B;
	int B2bits = 1; while (1 << B2bits < B2) B2bits++;
	int BB4bits = 1; while (1 << BB4bits < BB*4) BB4bits++;

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

	L[0] = q; L[1] = -r[l]; L[2] = 0;
	L[3] = 0; L[4] = 1;     L[5] = -r[l];
	L[6] = 0; L[7] = 0;     L[8] = 1;

	int64L2(L, 3);	// LLL reduce L, time log(q)^2

	// print special-q lattice
	//cout << "special-q lattice [" << L[0] << "," << L[1] << "," << L[2] << ";";
	//cout << L[3] << "," << L[4] << "," << L[5] << ";";
	//cout << L[6] << "," << L[7] << "," << L[8] << "]" << endl << flush;

	int64_t qLinv[9];
	// compute adjugate of L^-1
	qLinv[0] = L[4]*L[8]-L[5]*L[7];
	qLinv[1] = L[2]*L[7]-L[8]*L[1];
	qLinv[2] = L[5]*L[1]-L[2]*L[4];
	qLinv[3] = L[5]*L[6]-L[3]*L[8];
	qLinv[4] = L[8]*L[0]-L[2]*L[6];
	qLinv[5] = L[2]*L[3]-L[0]*L[5];
	qLinv[6] = L[3]*L[7]-L[4]*L[6];
	qLinv[7] = L[1]*L[6]-L[0]*L[7];
	qLinv[8] = L[0]*L[4]-L[1]*L[3];

	int i = 8; int m = 0;
	while (i < nump) {
		int64_t p = allp[i];
		if (p == q) { i++; continue; }
		//cout << p << "," << m << endl << flush;
		float logp = log2f(p);
		int64_t rl = mod(-r[l], q);
		for (int k = 0; k < num_smodp[i]; k++) {
			int64_t n = p*q;
			int64_t sk = mod(-s[i*degf+k],p);
			int64_t t = q*( (sk * modpinvq[i]) % p ) + p*( (rl * modqinvp[i]) % q ); // CRT
			if (t >= n) t -= n;
			L2[0] = n; L2[1] = t; L2[2] = 0;
			L2[3] = 0; L2[4] = 1; L2[5] = t;
			L2[6] = 0; L2[7] = 0; L2[8] = 1;
			int64L2(L2,3);	// LLL reduce L, time log(n)^2
			mat3x3prod(qLinv, L2, L3);	// L3 =  qLinv*L2
			int u1 = L3[0]/q; int u2 = L3[3]/q; int u3 = L3[6]/q;
			int v1 = L3[1]/q; int v2 = L3[4]/q; int v3 = L3[7]/q;
			int w1 = L3[2]/q; int w2 = L3[5]/q; int w3 = L3[8]/q;

			// compute normal (cross product)
			int nx = u2*v3 - u3*v2;
			int ny = u3*v1 - u1*v3;
			int nz = u1*v2 - u2*v1;

			// enumerate lattice vectors (x,y,z) in box [-B,B]x[-B,B]x[0,B]
			int x = w1; int y = w2; int z = w3;
			while (planeintersectsbox(nx, ny, nz, x-w1, y-w2, z-w3, B)) {
				x -= w1; y -= w2; z -= w3;
			}
			int j1, j2, j3, jmin;
			while (planeintersectsbox(nx, ny, nz, x, y, z, B)) {
				int ws1 = x; int ws2 = y; int ws3 = z;
				int a = 0; int b = 0;
				getab(u1, u2, u3, v1, v2, v3, x, y, z, B, &a, &b);
				x += a*u1 + b*v1;
				y += a*u2 + b*v2;
				z += a*u3 + b*v3;
				int s1 = x; int s2 = y; int s3 = z;
				// move 'forward' in plane
				bool inplane = true;
				while (inplane) {
					// pin vector to start of row
					int u1B = u1 < 0 ? B1 : -B; int u2B = u2 < 0 ? B1 : -B; int u3B = u3 < 0 ? B1 : 0;
					j1 = u1 == 0 ? -1 : (x - u1B) / u1;
					j2 = u2 == 0 ? -1 : (y - u2B) / u2;
					j3 = u3 == 0 ? -1 : (z - u3B) / u3;
					jmin = minnonneg(j1, j2, j3);
					x -= jmin*u1; y -= jmin*u2; z -= jmin*u3;
					if (x >= -B && x < B && y >= -B && y < B && z >= 0 && z < B) {
						int u1B = u1 < 0 ? -B : B1; int u2B = u2 < 0 ? -B : B1; int u3B = u3 < 0 ? 0 : B1;
						j1 = u1 == 0 ? -1 : (u1B - x) / u1;
						j2 = u2 == 0 ? -1 : (u2B - y) / u2;
						j3 = u3 == 0 ? -1 : (u3B - z) / u3;
						jmin = minnonneg(j1, j2, j3);
						for (int j = 0; j <= jmin; j++) {
							int id = (x + B) + ((y + B) << B2bits) + (z << BB4bits);
							M[m++] = (keyval){ id, logp };
							x += u1; y += u2; z += u3;
						}
						// move by '1-transition' vector
						x += v1; y += v2; z += v3;
					}
					else {
						inplane = false;	
					}
				}
				x = s1 - v1; y = s2 - v2; z = s3 - v3;
				// move 'backward' in plane
				inplane = true;
				while (inplane) {
					// pin vector to start of row
					int u1B = u1 < 0 ? B1 : -B; int u2B = u2 < 0 ? B1 : -B; int u3B = u3 < 0 ? B1 : 0;
					j1 = u1 == 0 ? -1 : (x - u1B) / u1;
					j2 = u2 == 0 ? -1 : (y - u2B) / u2;
					j3 = u3 == 0 ? -1 : (z - u3B) / u3;
					jmin = minnonneg(j1, j2, j3);
					x -= jmin*u1; y -= jmin*u2; z -= jmin*u3;
					if (x >= -B && x < B && y >= -B && y < B && z >= 0 && z < B) {
						int u1B = u1 < 0 ? -B : B1; int u2B = u2 < 0 ? -B : B1; int u3B = u3 < 0 ? 0 : B1;
						j1 = u1 == 0 ? -1 : (u1B - x) / u1;
						j2 = u2 == 0 ? -1 : (u2B - y) / u2;
						j3 = u3 == 0 ? -1 : (u3B - z) / u3;
						jmin = minnonneg(j1, j2, j3);
						for (int j = 0; j <= jmin; j++) {
							int id = (x + B) + ((y + B) << B2bits) + (z << BB4bits);
							M[m++] = (keyval){ id, logp };
							x += u1; y += u2; z += u3;
						}
						// move by '1-transition' vector
						x -= v1; y -= v2; z -= v3;
					}
					else {
						inplane = false;
					}
				}
				// advance to next plane
				x = ws1 + w1; y = ws2 + w2; z = ws3 + w3;
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


// Integer programming to get (a,b) such that (x,y,z) + a*(u1,u2,u3) + b*(v1,v2,v3) within [-B,B[x[-B,B[x[0,B[
inline void getab(int u1, int u2, int u3, int v1, int v2, int v3, int x, int y, int z, int B, int* a, int* b)
{
	int U[6] = { u1, u2, u3, -u1, -u2, -u3 };
	int V[6] = { v1, v2, v3, -v1, -v2, -v3 };
	int C[6] = { B-x-1, B-y-1, B-z-1, B+x, B+y, z };
	*a = B; *b = B;

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


// determine if plane with normal (nx, ny, nz) containing point (x,y,z) intersects box [-B,B]x[-B,B]x[0,B]
inline bool planeintersectsbox(int nx, int ny, int nz, int x, int y, int z, int B)
{
	int d = nx*x + ny*y + nz*z;

	int nxB0 = -nx*B; int nxB1 = nx*B;
	int nyB0 = -ny*B; int nyB1 = ny*B;
	int nzB1 = nz*B;

	int s0 = ( nxB0 + nyB0 + /*nz*0*/ - d ) > 0;
	int s1 = ( nxB1 + nyB0 + /*nz*0*/ - d ) > 0;
	int s2 = ( nxB0 + nyB1 + /*nz*0*/ - d ) > 0;
	int s3 = ( nxB1 + nyB1 + /*nz*0*/ - d ) > 0;
	int s4 = ( nxB0 + nyB0 +   nzB1   - d ) > 0;
	int s5 = ( nxB1 + nyB0 +   nzB1   - d ) > 0;
	int s6 = ( nxB0 + nyB1 +   nzB1   - d ) > 0;
	int s7 = ( nxB1 + nyB1 +   nzB1   - d ) > 0;

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
	return result;
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
bool EECM(mpz_t N, mpz_t S, mpz_t factor, int d, int a, int X0, int Y0, int Z0)
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


