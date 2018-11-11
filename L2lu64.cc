#include "L2lu64.h"
#include <math.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstdint>	// int64_t

using std::cout;
using std::endl;
using std::fixed;
using std::setprecision;
using std::setw;
using std::flush;

float n = 0.51f;
float delta = 0.99f;
float nn = (n+0.5f)/2.0f;
float dd = (delta+1.0f)/2.0f;

void printbasis(int64_t* b, int d, int pr, int w)
{
	for (int j = 0; j < d; j++) {
		for (int i = 0; i < d; i++) {
			cout << fixed << setprecision(pr) << setw(w) << b[j*d+i] << (i<d-1?",":"") << "\t";
		}
		cout << ";\\" << endl << flush;
	}
}

void int64L2(int64_t* b, int d)
{
	int d2 = d*d;
	int64_t* G = new int64_t[d2](); // subtle bug unless G initialized to zero
	float* rr = new float[d2];
	float* uu = new float[d2];

	// compute Gram matrix exactly
	for (int j = 0; j < d; j++) {
		for (int i = 0; i <= j; i++) {
			for (int l = 0; l < d; l++)
				G[j*d+i] += b[l*d+j]*b[l*d+i];
			G[i*d+j] = G[j*d+i];
		}
	}

	rr[0] = G[0];
	int k = 1;

	while (k < d) {
		// n size-reduce b[k]

		// compute rr[k][j]'s and uu[k][j]'s from G etc
		for (int i = 0; i <= k; i++) {
			rr[i*d+k] = G[i*d+k];
			for (int j = 0; j < i; j++)
				rr[i*d+k] -= rr[j*d+k]*uu[j*d+i];
			uu[i*d+k] = rr[i*d+k]/rr[i*d+i];
		}

		// compute max |uu[k*d+j]| for j < k
		float max = uu[0*d+k];
		for (int j = 0; j < k; j++)
			if (fabsf(uu[j*d+k]) > max)
				max = fabs(uu[j*d+k]);
		if (max > nn) {
			for (int j = k-1; j >= 0; j--) {
				int64_t X = floorf(uu[j*d+k] + 0.5f);
				for (int i = 0; i < d; i++)
					b[i*d+k] -= X*b[i*d+j];
				// update G
				for (int m = 0; m < d; m++) {
					for (int i = 0; i <= m; i++) {
						G[m*d+i] = 0;
						for (int l = 0; l < d; l++)
							G[m*d+i] += b[l*d+i]*b[l*d+m];
						G[i*d+m] = G[m*d+i];
					}
				}
				// update uu
				for (int i = 0; i < j; i++)
					uu[i*d+k] -= X*uu[i*d+j];
			}
			continue;
		}
		
		if (dd*rr[(k-1)*d+k-1] < rr[k*d+k] + uu[(k-1)*d+k]*uu[(k-1)*d+k]*rr[(k-1)*d+k-1]) {
			k++;
		}
		else {
			for (int i = 0; i < d; i++) {
				int64_t t = b[i*d+k-1];
				b[i*d+k-1] = b[i*d+k];
				b[i*d+k] = t;
			}
			// update G
			for (int j = 0; j < d; j++) {
				for (int i = 0; i <= j; i++) {
					G[j*d+i] = 0;
					for (int l = 0; l < d; l++)
						G[j*d+i] += b[l*d+i]*b[l*d+j];
					G[i*d+j] = G[j*d+i];
				}
			}
			// update rr[(k-1)*d+i] & uu[(k-1)*d+i]
			for (int i = 0; i <= k-1; i++) {
				rr[i*d+k-1] = G[i*d+k-1];
				for (int j = 0; j < i; j++)
					rr[i*d+k-1] -= rr[j*d+k-1]*uu[j*d+i];
				uu[i*d+k-1] = rr[i*d+k-1]/rr[i*d+i];
			}
			// update rr[k*d+i] & uu[k*d+i]
			for (int i = 0; i <= k; i++) {
				rr[i*d+k] = G[i*d+k];
				for (int j = 0; j < i; j++)
					rr[i*d+k] -= rr[j*d+k]*uu[j*d+i];
				uu[i*d+k] = rr[i*d+k]/rr[i*d+i];
			}
			
			k--;
			if (k < 1) k = 1;
		}
		// finished
	}
	
	// deallocate
	delete[] uu;
	delete[] rr;
	delete[] G;
}

void copysquareint64array(int64_t* src, int64_t* dest, int d)
{
    std::copy(src, src + d*d, dest);
}

void luinv(int64_t* b, int64_t* M, int d)
{
    // P*b = L*U
    // L*U*M = P*I
    // L*Y = P
    // U*M = Y
    // M = b^-1

	int d2 = d*d;
    int* P = new int[d2]();
    int64_t* L = new int64_t[d2];
    int64_t* U = new int64_t[d2];
    int64_t* Y = new int64_t[d2];
	for (int i = 0; i < d; i++) { 
		P[i] = i;
		L[i*d+i] = 1;
	}
    copysquareint64array(b, U, d);

    for (int i = 0; i < d-1; i++) {
        // find max entry in column i
        int jmax = i;
        for (int j = i; j < d; j++)
            if (abs(U[j*d+i]) > abs(U[jmax*d+i])) jmax = j;

        // swap row i and row jmax of U
        for (int j = 0; j < d; j++) {
            int64_t t = U[jmax*d+j];
            U[jmax*d+j] = U[i*d+j];
            U[i*d+j] = t;
        }
        int s = P[i];
		P[i] = P[jmax];
		P[jmax] = s;

        // swap row i and jmax of L up to diagonal
        for (int j = 0; j < i; j++) {
            int64_t t = L[jmax*d+j];
            L[jmax*d+j] = L[i*d+j];
            L[i*d+j] = t;
        }

        // reduce rows i+1...d
        for (int j = i+1; j < d; j++) {
            int64_t r = U[j*d+i]/U[i*d+i];
            L[j*d+i] = r;
            for (int k = i; k < d; k++)
                U[j*d+k] -= U[i*d+k] * r;
        }
    }
    
    // solve for Y in L*Y = P
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            Y[j*d+i] = (i == P[j]);
            for (int k = 0; k < j; k++) {
                Y[j*d+i] -= Y[k*d+i]*L[j*d+k];
            }
        }
    }

    // solve for M in U*M = Y
    for (int i = 0; i < d; i++) {
        for (int j = d-1; j >= 0; j--) {
            M[j*d+i] = Y[j*d+i]/U[j*d+j];
            for (int k = d-1; k > j; k--) {
            	M[j*d+i] -= M[k*d+i]*U[j*d+k]/U[j*d+j];
            }
        }
    }

    delete[] Y;
    delete[] U;
    delete[] L;
    delete[] P;
}

