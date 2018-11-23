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

struct keyval {
	int id;
	float logp;
};

int latsieve3d(int* f, int degf, int64_t q, int l, int* allp, int nump, int* s, int* num_smodp, keyval* M, int Mlen, int B);
void histogram(keyval*M, float* H, int len);
bool lattice_sorter(keyval const& kv1, keyval const& kv2);
void csort(keyval* M, keyval* L, int* H, int len);
inline int floordiv(int a, int b);
inline int modinv(int x, int m);
inline int nonzerolcm(int u1, int u2, int u3);
inline int gcd(int a, int b);
inline void getab(int u1, int u2, int u3, int v1, int v2, int v3, int x, int y, int z, int B, int* a, int* b);
inline bool planeintersectsbox(int nx, int ny, int nz, int x, int y, int z, int B);
inline int minnonneg(int u, int v, int w);
inline int minabs(int u, int v, int w);
inline int maxabs(int u, int v, int w);
inline int min(int u, int v, int w);
inline int max(int u, int v, int w);


int main(int argc, char** argv)
{
	bool verbose = true;
		
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
	mpz_t t;
	mpz_init(t);
	//mpz_t c;
	//mpz_init(c);
	string line;
	char linebuffer[100];
	ifstream file(argv[1]);
	getline(file, line);	// first line contains number n to factor
	// read nonlinear poly
	int degf = -1;
	if (verbose) cout << endl << "nonlinear polynomial f(x): (ascending coefficients)" << endl;
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
	if (verbose) cout << endl << "linear polynomial g(x): (ascending coefficients)" << endl;
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
	if (verbose) cout << endl << "Complete." << endl;

	if (verbose) cout << endl << "Starting sieve of Eratosthenes for small primes..." << endl << flush;
	int max = 1<<25; // 10000000;// 65536;
	char* sieve = new char[max+1]();
	int* primes = new int[2063689]; //new int[809228];	//new int[6542]; 	// 2039 is the 309th prime, largest below 2048
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
	int K = 0;
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		if (id == 0) K = omp_get_num_threads();
	}
	mpz_t r; mpz_init(r);
	int* s0 = new int[degf * nump]();
	int* sieves0 = new int[degf * nump]();
	int* num_s0modp = new int[nump]();
	int* sievep0 = new int[nump]();
	int* sievenum_s0modp = new int[nump]();
	int* s1 = new int[degg * nump]();
	int* sieves1 = new int[degg * nump]();
	int* num_s1modp = new int[nump]();
	int* sievep1 = new int[nump]();
	int* sievenum_s1modp = new int[nump]();
	int itenpc0 = nump / 10;
	int itotal = 15;
	// compute factor base
	if (verbose) cout << endl << "Constructing factor base with " << K << " threads." << endl << flush;
	//if (verbose) cout << endl << "[0%]   constructing factor base..." << endl << flush;
	start = clock();
	#pragma omp parallel
	{
		mpz_t rt; mpz_init(rt); 
		int id = omp_get_thread_num();
		int* stemp0 = new int[degf];
		int* fp = new int[degf+1]();
		int* stemp1 = new int[degg];
		int* gp = new int[degg+1]();

	#pragma omp for
		for (int i = 15; i < nump; i++) {
			int p = primes[i];
			for (int j = 0; j <= degf; j++) fp[j] = mpz_mod_ui(rt, fpoly[j], p);
			int degfp = degf; while (fp[degfp] == 0 || degfp == 0) degfp--;
			int nums0 = polrootsmod(fp, degfp, stemp0, p);
			num_s0modp[i] = nums0;
			for (int j = 0; j < nums0; j++) s0[i*degf + j] = stemp0[j];
			for (int j = 0; j <= degg; j++) gp[j] = mpz_mod_ui(rt, gpoly[j], p);
			int deggp = degg; while (gp[deggp] == 0 || deggp == 0) deggp--;
			int nums1 = polrootsmod(gp, deggp, stemp1, p);
			num_s1modp[i] = nums1;
			for (int j = 0; j < nums1; j++) s1[i*degg + j] = stemp1[j];

	#pragma omp atomic
			itotal++;
			if (itotal % itenpc0 == 0) {
	#pragma omp critical
				cout << "[" << 100 * itotal / nump + 1 << "%]\tConstructing factor base..." << endl << flush;
			}				 
		}

		delete[] gp;
		delete[] stemp1;
		delete[] fp;
		delete[] stemp0;
		mpz_clear(rt);
	}
	timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC / K;
	start = clock();
	int k0 = 0; int k1 = 0;
	for (int i = 15; i < nump; i++) {
		int nums0 = num_s0modp[i];
		if (nums0 > 0) {
			sievep0[k0] = primes[i];
			for (int j = 0; j < nums0; j++) sieves0[k0*degf + j] = s0[i*degf + j];
			sievenum_s0modp[k0++] = nums0;
		}
		int nums1 = num_s1modp[i];
		if (nums1 > 0) {
			sievep1[k1] = primes[i];
			for (int j = 0; j < nums1; j++) sieves1[k1*degg + j] = s1[i*degg + j];
			sievenum_s1modp[k1++] = nums1;
		}
	}
	timetaken += ( clock() - start ) / (double) CLOCKS_PER_SEC;
	if (verbose) cout << "Complete.  Time taken: " << timetaken << "s" << endl << flush;
	if (verbose) cout << "There are " << k0 << " factor base primes on side 0." << endl << flush;
	if (verbose) cout << "There are " << k1 << " factor base primes on side 1." << endl << flush;
	
	int B = 256;//512;
	int Mlen = 512*512*256*2;//1024*1024*512*2;
	keyval* M = new keyval[Mlen];	// lattice { id, p } pairs
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
	int q = 12345701;// 65537;
	if (argc >= 3) q = atoi(argv[2]);
	int* fq = new int[degf+1]();
	for (int i = 0; i <= degf; i++) fq[i] = mpz_mod_ui(r, fpoly[i], q);
	// sieve side 0
	cout << "Starting sieve on side 0 for special-q " << q << "..." << endl << flush;
	start = clock();
	int m = latsieve3d(fq, degf, q, 0, sievep0, k0, sieves0, sievenum_s0modp, M, Mlen, B);
	timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
	cout << "Finished! Time taken: " << timetaken << "s" << endl << flush;
	cout << "Size of lattice point list is " << m << "." << endl << flush;
	cout << "Constructing histogram..." << endl << flush;
	start = clock();
	//std::stable_sort(M, M + m, &lattice_sorter);
	histogram(M, H, m);
	timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
	cout << "Finished! Time taken: " << timetaken << "s" << endl << flush;
	float th = 100.0f;
	if (argc >= 4) th = atof(argv[3]);
	int R0 = 0;
	int B2 = 2*B; int BB = B*B;
	int B2bits = 1; while (1 << B2bits < B2) B2bits++;
	int BB4bits = 1; while (1 << BB4bits < BB*4) BB4bits++;
	for (int i = 0; i < m; i++) {
		if (H[i] > th) {
			rel.push_back(i);
			R0++;
		}
	}
	cout << R0 << " candidates on side 0." << endl << flush;
	// sieve side 1
	cout << "Starting sieve on side 1..." << endl << flush;
	start = clock();
	m = latsieve3d(fq, degf, q, 0, sievep1, k1, sieves1, sievenum_s1modp, M, Mlen, B);
	timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
	cout << "Finished! Time taken: " << timetaken << "s" << endl << flush;
	cout << "Size of lattice point list is " << m << "." << endl << flush;
	cout << "Constructing histogram..." << endl << flush;
	start = clock();
	//std::stable_sort(M, M + m, &lattice_sorter);
	histogram(M, H, m);
	timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
	cout << "Finished! Time taken: " << timetaken << "s" << endl << flush;
	th = 100.0f;
	if (argc >= 4) th = atof(argv[3]);
	int R1 = 0;
	for (int i = 0; i < m; i++) {
		if (H[i] > th) {
			rel.push_back(i);
			R1++;
		}
	}
	cout << R1 << " candidates on side 1." << endl << flush;
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
				cout << x << "," << y << "," << z << endl;
				R++;
			}
		}
	}
	cout << R << " potential relations found." << endl << flush;
	
	delete[] fq;
	delete[] H;
	//delete[] L;
	delete[] M;
	delete[] sievenum_s1modp;
	delete[] sievep1;
	delete[] num_s1modp;
	delete[] s1;
	delete[] sievenum_s0modp;
	delete[] sievep0;
	delete[] num_s0modp;
	delete[] s0;
	mpz_clear(r);
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
	cout << "special-q (" << q << "," << r[l] << ")" << endl;

	// compute modpinvq and modqinvp arrays
	int* modpinvq = new int[nump];
	int* modqinvp = new int[nump];
	for (int i = 0; i < nump; i++) {
		int p = allp[i];
		if (p == q) continue;
		modpinvq[i] = modinv(q, p);
		modqinvp[i] = modinv(p, q);
	}

	L[0] = q; L[1] = -r[l]; L[2] = 0;
	L[3] = 0; L[4] = 1;     L[5] = -r[l];
	L[6] = 0; L[7] = 0;     L[8] = 1;

	int64L2(L, 3);	// LLL reduce L, time log(q)^2

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

	int i = 0; int m = 0;
	while (i < nump) {
		int64_t p = allp[i];
		if (p == q) { i++; continue; }
		//cout << p << "," << m << endl << flush;
		float logp = logf(p);
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

