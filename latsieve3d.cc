#include <stdint.h>	// int64_t
#include <iostream> // cout
#include <iomanip> // setprecision
#include "L2lu64.h"
#include <gmpxx.h>
#include "intpoly.h"
#include <math.h>	// sqrt
#include <fstream>	// file
#include <ctime>	// clock_t

using std::cout;
using std::endl;
using std::flush;
//using std::vector;
using std::string;
using std::ifstream;
using std::fixed;
using std::scientific;
using std::setprecision;

struct keyval {
	int id;
	float logp;
};

int latsieve3d(int* f, int degf, int64_t q, int* allp, int nump, int* s, int* num_smodp, keyval* M, int Mlen, int B);
inline int modinv(int x, int m);
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
		mpz_get_str(linebuffer, 10, fpoly[degf-1]);
		if (verbose) cout << line << endl << flush;
	}
	//int degf = fpoly.size();
	// read other poly
	int degg = 0;
	bool read = true;
	if (verbose) cout << endl << "linear polynomial g(x): (ascending coefficients)" << endl;
	while (read && line.substr(0,1) == "Y" ) {
		line = line.substr(line.find_first_of(" ")+1);
		//mpz_set_str(c, line.c_str(), 10);
		mpz_set_str(gpoly[degg++], line.c_str(), 10);
		mpz_get_str(linebuffer, 10, gpoly[degg-1]);
		if (verbose) cout << line << endl << flush;
		read = getline(file, line); 
	}
	//int degg = gpoly.size();
	file.close();
	//mpz_clear(c);
	if (verbose) cout << endl << "Complete." << endl;

	if (verbose) cout << endl << "Starting sieve of Eratosthenes for small primes..." << endl << flush;
	int max = 10000000;// 65536;
	char* sieve = new char[max+1]();
	int* primes = new int[809228];	//new int[6542]; 	// 2039 is the 309th prime, largest below 2048
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
	mpz_t r; mpz_init(r);
	int* s = new int[degf * nump]();
	int* stemp = new int[degf];
	int* num_smodp = new int[nump]();
	int* fp = new int[degf+1]();
	int* sievep = new int[nump]();
	// compute factor base
	if (verbose) cout << endl << "Computing factor base..." << endl << flush;
	int k = 0;
	//# p ragma omp parallel for
	for (int i = 15; i < nump; i++) {
		int p = primes[i];
		for (int j = 0; j <= degf; j++) fp[j] = mpz_mod_ui(r, fpoly[j], p);
		int degfp = degf; while (fp[degfp] == 0 || degfp == 0) degfp--;
		int nums = polrootsmod(fp, degfp, stemp, p);
		if (nums) {
			for (int j = 0; j < nums; j++) s[k*degf + j] = stemp[j];
			num_smodp[k] = nums;
			sievep[k++] = p;
		}
		if (i % 1000 == 0) cout << i << "," << flush;
	}
	if (verbose) cout << endl << "Complete. There are " << k << " factor base primes." << endl;
	
	int B = 512;
	int Mlen = 1024*1024*512*2;// 512*512*256*10;
	keyval* M = new keyval[Mlen];
	int q = 12345701;// 65537;
	int* fq = new int[degf+1]();
	for (int i = 0; i <= degf; i++) fq[i] = mpz_mod_ui(r, fpoly[i], q);
	cout << "Starting sieve for special-q " << q << "..." << endl << flush;
	std::clock_t start = clock();
	int m = latsieve3d(fq, degf, q, sievep, k, s, num_smodp, M, Mlen, B);
	double timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC;
	cout << "Finished! Time taken: " << timetaken << "s" << endl << flush;
	cout << "Number of lattice points is " << m << "." << endl << flush;

	delete[] fq;
	delete[] M;
	delete[] fp;
	delete[] num_smodp;
	delete[] s;
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


int latsieve3d(int* f, int degf, int64_t q, int* allp, int nump, int* s, int* num_smodp, keyval* M, int Mlen, int B)
{
	int64_t L[9];
	int64_t L2[9];
	int64_t L3[9];
	int* r = new int[degf]();
	int numl = polrootsmod(f, degf, r, q);
	int B2 = 2*B; int BB = B*B;
	int B2bits = 1; while (1 << B2bits < B2) B2bits++;
	int BB4bits = 1; while (1 << BB4bits < BB*4) BB4bits++;

	// compute modpinvq and modqinvp arrays
	int* modpinvq = new int[nump];
	int* modqinvp = new int[nump];
	for (int i = 0; i < nump; i++) {
		int p = allp[i];
		modpinvq[i] = modinv(q, p);
		modqinvp[i] = modinv(p, q);
	}

	// clear M
	for (int j = 0; j < Mlen; j++) M[j] = (keyval){ 0, 0 };
	cout << "Memory cleared." << endl << flush;
	int m = 0;

	for (int l = 0; l <= 0 /* < numl*/; l++) {
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

		// clear M
		//for (int j = 0; j < Mlen; j++) M[j] = (keyval){ 0, 0 };
		int i = 0; //int m = 0;
		while (i < nump) {
			if (num_smodp[i] == 0) continue;	// skip p if no roots mod p
			int64_t p = allp[i];
			//cout << p << "," << flush;
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
						int u1B = u1 < 0 ? B : -B; int u2B = u2 < 0 ? B : -B; int u3B = u3 < 0 ? B : 0;
						j1 = u1 == 0 ? -1 : (x - u1B) / u1;
						j2 = u2 == 0 ? -1 : (y - u2B) / u2;
						j3 = u3 == 0 ? -1 : (z - u3B) / u3;
						jmin = minnonneg(j1, j2, j3);
						x -= jmin*u1; y -= jmin*u2; z -= jmin*u3;
						if (x >= -B && x <= B && y >= -B && y <= B && z >= 0 && z <= B) {
							int u1B = u1 < 0 ? -B : B; int u2B = u2 < 0 ? -B : B; int u3B = u3 < 0 ? 0 : B;
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
						int u1B = u1 < 0 ? B : -B; int u2B = u2 < 0 ? B : -B; int u3B = u3 < 0 ? B : 0;
						j1 = u1 == 0 ? -1 : (x - u1B) / u1;
						j2 = u2 == 0 ? -1 : (y - u2B) / u2;
						j3 = u3 == 0 ? -1 : (z - u3B) / u3;
						jmin = minnonneg(j1, j2, j3);
						x -= jmin*u1; y -= jmin*u2; z -= jmin*u3;
						if (x >= -B && x <= B && y >= -B && y <= B && z >= 0 && z <= B) {
							int u1B = u1 < 0 ? -B : B; int u2B = u2 < 0 ? -B : B; int u3B = u3 < 0 ? 0 : B;
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
	}
	// clear memory
	delete[] modqinvp;
	delete[] modpinvq;
	delete[] r;

	return m;
}


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
	int b1 = d1 == 0 ? -1 : d1 < 0 ? n1/d1 : m1/d1;
	int b2 = d2 == 0 ? -1 : d2 < 0 ? n2/d2 : m2/d2;
	int b3 = d3 == 0 ? -1 : d3 < 0 ? n3/d3 : m3/d3;
	*b = minnonneg(b1, b2, b3);
	int a1 = u1 == 0 ? -1 : u1 < 0 ? ( -B - x - v1 * (*b) ) / u1 : ( B - x - v1 * (*b) ) / u1;
	int a2 = u2 == 0 ? -1 : u2 < 0 ? ( -B - y - v2 * (*b) ) / u2 : ( B - y - v2 * (*b) ) / u2;
	int a3 = u3 == 0 ? -1 : u3 < 0 ? (  0 - z - v3 * (*b) ) / u3 : ( B - z - v3 * (*b) ) / u3;
	*a = minnonneg(a1, a2, a3);
}


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
	if (u < 0) m = v; if (v < 0) m = w;
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

