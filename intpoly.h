#ifndef INTPOLY_H
#define INTPOLY_H

#include <vector>	// vector

using std::vector;

int primitiveroot(int p, int* q);
int modpow(int64_t r, int64_t e, int64_t p);
int chienrootsmod(int* f, int degf, int p, int* r, int* q);
int mod(int r, int p);
inline int modinv(int x, int m);
inline int poldegree(int64_t* f, int maxd);
inline void polsetzero(int64_t* f, int maxd);
inline bool poliszero(int64_t* f, int maxd);
int polrootsmod(int* f, int degf, int* roots, int p);
#endif	/* INTPOLY_H */

