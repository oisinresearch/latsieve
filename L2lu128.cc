#include "L2lu128.h"
#include <math.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstdint>	// int64_t
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/lexical_cast.hpp>
#include <string>

using namespace boost::multiprecision;

using boost::lexical_cast;
using std::string;
using std::cout;
using std::endl;
using std::fixed;
using std::setprecision;
using std::setw;
using std::flush;

//float n = 0.51f;
//float delta = 0.99f;
//float nn = (n+0.5f) / 2.0f;
//float dd = (delta+1.0f) / 2.0f;
extern float n;
extern float delta;
extern float nn;
extern float dd;

void int128L2(int128_t* borig, int d)
{
	int d2 = d*d;
	int256_t* b = new int256_t[d2];
	int256_t* G = new int256_t[d2](); // subtle bug unless G initialized to zero
	float128* rr = new float128[d2];
	float128* uu = new float128[d2];

	// copy into b
	for (int i = 0; i < d2; i++) b[i] = borig[i].convert_to<int256_t>();

	// compute Gram matrix exactly
	for (int j = 0; j < d; j++) {
		for (int i = 0; i <= j; i++) {
			for (int l = 0; l < d; l++)
				G[j*d+i] += b[l*d+j] * b[l*d+i];
			G[i*d+j] = G[j*d+i];
		}
	}

	rr[0] = G[0].convert_to<float128>();
	int k = 1;

	while (k < d) {
		// n size-reduce b[k]

		// compute rr[k][j]'s and uu[k][j]'s from G etc
		for (int i = 0; i <= k; i++) {
			rr[i*d+k] = G[i*d+k].convert_to<float128>();
			for (int j = 0; j < i; j++)
				rr[i*d+k] -= rr[j*d+k] * uu[j*d+i];
			uu[i*d+k] = rr[i*d+k] / rr[i*d+i];
		}

		// compute max |uu[k*d+j]| for j < k
		float128 max = uu[0*d+k];
		for (int j = 0; j < k; j++) {
			float128 fabsuu = uu[j*d+k];
			if (fabsuu < 0) fabsuu *= -1;
			if (fabsuu > max) max = fabsuu;
		}
		if (max > nn) {
			for (int j = k-1; j >= 0; j--) {
				int256_t X = floor(uu[j*d+k] + 0.5).convert_to<int256_t>();
				//int256_t X = static_cast<int256_t>(floorl(uu[j*d+k] + 0.5f));
                float128 Xf = X.convert_to<float128>();
				for (int i = 0; i < d; i++)
					b[i*d+k] -= X * b[i*d+j];
				// update G
				for (int m = 0; m < d; m++) {
					for (int i = 0; i <= m; i++) {
						G[m*d+i] = 0;
						for (int l = 0; l < d; l++)
							G[m*d+i] += b[l*d+i] * b[l*d+m];
						G[i*d+m] = G[m*d+i];
					}
				}
				// update uu
				for (int i = 0; i < j; i++)
					uu[i*d+k] -= Xf * uu[i*d+j];
			}
			continue;
		}
		
		if (dd * rr[(k-1)*d+k-1] < rr[k*d+k] + uu[(k-1)*d+k]*uu[(k-1)*d+k] * rr[(k-1)*d+k-1]) {
			k++;
		}
		else {
			for (int i = 0; i < d; i++) {
				int256_t t = b[i*d+k-1];
				b[i*d+k-1] = b[i*d+k];
				b[i*d+k] = t;
			}
			// update G
			for (int j = 0; j < d; j++) {
				for (int i = 0; i <= j; i++) {
					G[j*d+i] = 0;
					for (int l = 0; l < d; l++)
						G[j*d+i] += b[l*d+i] * b[l*d+j];
					G[i*d+j] = G[j*d+i];
				}
			}
			// update rr[(k-1)*d+i] & uu[(k-1)*d+i]
			for (int i = 0; i <= k-1; i++) {
				rr[i*d+k-1] = G[i*d+k-1].convert_to<float128>();
				for (int j = 0; j < i; j++)
					rr[i*d+k-1] -= rr[j*d+k-1] * uu[j*d+i];
				uu[i*d+k-1] = rr[i*d+k-1] / rr[i*d+i];
			}
			// update rr[k*d+i] & uu[k*d+i]
			for (int i = 0; i <= k; i++) {
				rr[i*d+k] = G[i*d+k].convert_to<float128>();
				for (int j = 0; j < i; j++)
					rr[i*d+k] -= rr[j*d+k] * uu[j*d+i];
				uu[i*d+k] = rr[i*d+k] / rr[i*d+i];
			}
			
			k--;
			if (k < 1) k = 1;
		}
		// finished
	}

	// write output
	for (int i = 0; i < d2; i++) borig[i] = b[i].convert_to<int128_t>();

	// deallocate
	delete[] uu;
	delete[] rr;
	delete[] G;
	delete[] b;
}


void matprint(int d, int128_t* M)
{
	int maxw = 0;
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < d; j++) {
			string coeff = lexical_cast<string>(M[i*d+j]);
			if (coeff.length() > maxw) maxw = coeff.length();
		}
	}
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < d; j++) {
			string coeff = lexical_cast<string>(M[i*d+j]);
			cout << setw(maxw) << coeff << " ";
		}
		cout << endl << endl;
	}
}

