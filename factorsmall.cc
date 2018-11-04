#include <cstdlib>	// atoi
#include <vector>	// vector
#include <iostream>	// cout
#include <cmath>	// sqrt
#include "factorsmall.h"

#define MAXP 3512
#define MAXS 32767

using std::vector;
using std::cout;
using std::endl;

#ifdef FACTORMAIN
int main(int argc, char** argv)
{
	if (argc == 2) {
		int n = atoi(argv[1]);

		// sieve of Eratosthenes
		int* primes = new int[MAXP]; 	// MAXP primes up to 32767
		int nump = Eratosthenes(primes);

		// store result
		vector<int> p;
		vector<int> e;

		// factor n
		int k = factorsmall(n, primes, p, e);

		// print result
		for (int i = 0; i < k; i++) {
			cout << p[i];
			if (e[i] > 1) cout << "^" << e[i];
			cout << (i < k-1 ? " x " : "");
		}
		cout << endl;

		// clean up
		delete[] primes;
	}
	return 0;
}
#endif /* FACTORMAIN */

int Eratosthenes(int* primes)
{
	char* sieve = new char[MAXS]();
	int imax = (int)(sqrt(MAXS - sqrt(MAXS)) + 0.5f);
	for (int i = 2; i <= imax; i++)
		if(!sieve[i])
			for (int j = i*i; j <= MAXS; j += i)
				if(!sieve[j]) sieve[j] = 1;
	int nump = 0;
	for (int i = 2; i < MAXS; i++)
		if (!sieve[i])
			primes[nump++] = i;

	delete[] sieve;
	return nump;
}

// factor n into p1^e1 x p2^e2 x ... x pk^ek using trial division by qi up to floor(sqrt(n))
// n up to about 30 bits
int factorsmall(int n, int *q, vector<int>& p, vector<int>& e)
{
	int nq = n;
	int i = 0; int k = 0;
	while (nq > 1) {
		int qi = q[i++];
		if (nq % qi == 0) {
			k++;
			p.push_back(qi);
			e.push_back(0);
			do {
				e.back()++;
				nq /= qi;
			}
			while (nq % qi == 0);
		}
		if (i >= MAXP) break;
	}
	if (i >= MAXP && n == nq) { // n is prime
		k = 1;
		p.push_back(n);
		e.push_back(1);
	}
	if (i >= MAXP && nq < n && nq > 1) {	// cofactor is prime
		k++;
		p.push_back(nq);
		e.push_back(1);
	}
	return k;
}

