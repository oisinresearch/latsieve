#ifndef INTPOLY_H
#define INTPOLY_H

#include <vector>	// vector

using std::vector;

int primitiveroot(int p, int* q);
int modpow(int64_t r, int64_t e, int64_t p);
int chienrootsmod(int* f, int degf, int p, int* r, int* q);
int64_t mod(int64_t r, int64_t p);
int64_t mod128(__int128 r, int64_t p);
inline int64_t modinv(int64_t x, int64_t m);
inline int poldegree(__int128* f, int maxd);
inline void polsetzero(__int128* f, int maxd);
inline bool poliszero(__int128* f, int maxd);
int polrootsmod(int64_t* f, int degf, int64_t* roots, int64_t p);
#endif	/* INTPOLY_H */

