#include <stdint.h>	// int64_t
#include <iostream> // cout
#include <iomanip> // setprecision
#include "L2lu64.h"
#include <gmpxx.h>

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

void latsieve3d(mpz_t* f, int degf, int q, int* allp, int nump, int* s, int* num_smodp, keyval* M, int Mlen, int B);
int modinv(int x, int m);
inline bool planeintersectsbox(int u1, int u2, int u3, int v1, int v2, int v3, int w1, int w2, int w3, int B);
int minabs(int u, int v, int w);


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
	int degf = 12;
	int deg = 0;
	if (verbose) cout << endl << "nonlinear polynomial f(x): (ascending coefficients)" << endl;
	while (getline(file, line) && line.substr(0,1) == "c" ) {
		line = line.substr(line.find_first_of(" ")+1);
		//mpz_set_str(c, line.c_str(), 10);
		mpz_set_str(fpoly[deg++], line.c_str(), 10);
		mpz_get_str(linebuffer, 10, fpoly[deg-1]);
		if (verbose) cout << line << endl << flush;
	}
	//int degf = fpoly.size();
	// read other poly
	int degg = 1;
	deg = 0;
	bool read = true;
	if (verbose) cout << endl << "linear polynomial g(x): (ascending coefficients)" << endl;
	while (read && line.substr(0,1) == "Y" ) {
		line = line.substr(line.find_first_of(" ")+1);
		//mpz_set_str(c, line.c_str(), 10);
		mpz_set_str(gpoly[deg++], line.c_str(), 10);
		mpz_get_str(linebuffer, 10, gpoly[deg-1]);
		if (verbose) cout << line << endl << flush;
		read = getline(file, line); 
	}
	//int degg = gpoly.size();
	file.close();
	//mpz_clear(c);
	if (verbose) cout << endl << "Complete." << endl;

	if (verbose) cout << endl << "Starting sieve of Eratosthenes for small primes..." << endl << flush;
	int max = 2048;
	char* sieve = new char[max+1]();
	int* primes = new int[309]; 	// 2039 is the 309th prime, largest below 2048
	for (int i = 2; i <= sqrt(max); i++)
		if(!sieve[i])
			for (int j = i*i; j <= max; j += i)
				if(!sieve[j]) sieve[j] = 1;
	int nump = 0;
	for (int i = 2; i <= 2047; i++)
		if (!sieve[i])
			primes[nump++] = i;
	if (verbose) cout << "Complete." << endl;

	// set up constants and call sieve
	int* s = new int[nump]();
	int* num_smodp = new int[nump];
	for (int i = 0; i < nump; i++) num_smodp[i] = 1;
	int B = 100;
	int Mlen = 200*200*100*10;
	keyval* M = new keyval(Mlen);
	int q = 12373;
	latsieve3d(f, degf, q, primes, s, num_smodp, nump, M, Mlen, B);

	return 0;
}


void latsieve3d(mpz_t* f, int degf, int q, int* allp, int nump, int* s, int* num_smodp, keyval* M, int Mlen, int B)
{
	int64_t[9] L;
	int64_t[9] L2;
	int64_t[9] L3;
	int* r = new int[degf];
	int numl = polrootsmod(f, degf, q, r);
	int B2 = 2*B; int BB = B*B;
	int B2bits = 1; while (1 << B2bits < B2) B2bits++;
	int BBbits = 1; while (1 << BBbits < BB) BBbits++;

	// compute modpinvq and modqinvp arrays
	int* modpinvq = new int[nump];
	int* modqinvp = new int[nump];
	for (i = 0; i < nump; i++) {
		modpinvq[i] = modinv(q, p);
		modqinvp[i] = modinv(p, q);
	}

	for (l = 0; l < numl; l++) {
		L[0] = q; L[1] = r[l]; L[2] = 0;
		L[3] = 0; L[4] = 1l;   L[5] = r[l];
		L[6] = 0; L[7] = 0;    L[8] = 1l;

		int64L2(L, 3);	// LLL reduce L, time log(q)^2

		int64_t[9] qLinv;
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
			p = allp[i];
			float logp = logf(p);

			for (int k = 0; k < num_smodp[i]; k++) {
				int n = p*q;
				int t = q*( (r[l]*modpinvq[i]) % p ) + p*( (s[k]*modqinvp[i]) % q ); // CRT
				if (t >= n) t -= n;
				L2[0] = n; L2[1] = t;  L2[2] = 0;
				L2[3] = 0; L2[4] = 1l; L2[5] = t;
				L2[6] = 0; L2[7] = 0;  L2[8] = 1l;
				int64L2(L2,3);	// LLL reduce L, time log(n)^2
				mat3x3prod(L, L2, L3);	// L3 =  L*L2
				int u1 = L3[0]/q; int u2 = L3[3]/q; int u3 = L3[6]/q;
				int v1 = L3[1]/q; int v2 = L3[4]/q; int v3 = L3[7]/q;
				int w1 = L3[2]/q; int w2 = L3[5]/q; int w3 = L3[8]/q;

				u1 = -9; u2 =  -9; u3 = 11;
				v1 = 10; v2 = -12; v3 = -7;
				w1 = 25; w2 =  25; w3 = 25;
				
				// enumerate lattice vectors (x,y,z) in box [-B,B]x[-B,B]x[0,B]
				int x = w1; int y = w2; int z = w3;
				while (planeintersectsbox(u1, u2, u3, v1, v2, v3, x, y, z, B)) {
					int s1 = x; int s2 = y; int s3 = z;
					// move 'forward' in plane
					bool inplane = true;
					while (inplane) {
						int u1B = u1 < 0 ? B : -B; int u2B = u2 < 0 ? B : -B; int u3B = u3 < 0 ? B : 0;
						int jright = (u1B + x) / u1;
						int jtop   = (u2B + y) / u2;
						int jfront = (u3B + z) / u3;
						int jmax = minabs(jright, jtop, jfront);
						for (int j = 0; j < jmax; j++) {
							x += u1; y += u2; z += u3;
							int id = (x + B) + ((y + B) << B2bits) + (z << BBbits);
							M[m++] = (keyval){ id, logp };
						}
						// move by '1-transition' vector
						x += v1; y += v2; z += v3;
						u1B = -u1B; u2B = -u2B; u3B = u3 < 0 ? B : 0;
						int jleft  = (u1B + x) / u1;
						int jfloor = (u2B + y) / u2;
						int jback  = (u3B + z) / u3;
						int jmin = minabs(jleft, jfloor, jback);
						x -= jmin*u1; y -= jmin*u2; z -= jmin*u3;
						if (abs(x) > B || abs(y) > B || z < 0 || z > B)
							inplane = false;							
					}
					x = s1; y = s2; z = s3;
					// move 'backward' in plane
					inplane = true;
					while (inplane) {
						int u1B = u1 < 0 ? -B : B; int u2B = u2 < 0 ? -B : B; int u3B = u3 < 0 ? 0 : B;
						int jleft  = (u1B + x) / u1;
						int jfloor = (u2B + y) / u2;
						int jback  = (u3B + z) / u3;
						int jmin = minabs(jleft, jfloor, jback);
						for (int j = 0; j <= jmin; j++) {
							int id = (x + B) + ((y + B) << B2bits) + (z << BBbits);
							M[m++] = (keyval){ id, logp };
							x -= u1; y -= u2; z -= u3;
						}
						// move by '1-transition' vector
						x -= v1; y -= v2; z -= v3;
						u1B = -u1B; u2B = -u2B; u3B = u3 < 0 ? 0 : B;
						int jright = (u1B + x) / u1;
						int jtop   = (u2B + y) / u2;
						int jfront = (u3B + z) / u3;
						int jmax = minabs(jright, jtop, jfront);
						x += jmax*u1; y += jmax*u2; z += jmax*u3;
						if (abs(x) > B || abs(y) > B || z < 0 || z > B)
							inplane = false;							
					}
					// advance to next plane
					x = s1 + w1; y = s2 + w2; z = s3 + w3;
				}
			}
		}
	}
	// clear memory
	delete[] modqinvp;
	delete[] modpinvq;
	delete[] r;
}


int modinv(int x, int m)
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


inline bool planeintersectsbox(int u1, int u2, int u3, int v1, int v2, int v3, int w1, int w2, int w3, int B)
{
	// cross product
	int nx = u2*v3 - u3*v2;
	int ny = u3*v1 - u1*v3;
	int nz = u1*v2 - u2*v1;

	// bounds
	int sx = (nx >= 0 ? -B : B); int tx = (nx >= 0 ? B : -B);
	int sy = (ny >= 0 ? -B : B); int ty = (ny >= 0 ? B : -B);
	int sz = (nz >= 0 ?  0 : B); int tz = (nz >= 0 ? B :  0);

	float n2n = nx*nx + ny*ny + nz*nz;
	int wdotnx = w1*nx; int wdotny = w2*ny; int wdotnz = w3*nz;
	float d = sqrt( (wdotnx*wdotnx + wdotny*wdotny + wdotnz*wdotnz) / n2n );
	float pos_side = nx*tx + ny*ty + nz*tz + d;
	float neg_side = nx*sx + ny*xy + nz*sz + d;

	// result
	bool intersects = true;
	if (pos_side > 0.0f || neg_side < 0.0f) intersects = false;
	return intersects;
}


int minabs(int u, int v, int w)
{
	int m = u;
	if (abs(v) < m) m = v;
	if (abs(w) < m) m = w;
	return m;
}


