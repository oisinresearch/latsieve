#include <cstdint>  // int64_t

#ifndef L2lu_h
#define L2lu_h

void int64L2(int64_t* b, int d);
void printbasis(int64_t* b, int d, int pr, int w);
void copysquarefloatarray(int64_t* src, int64_t* dest, int d);
void luinv(int64_t* b, int64_t* M, int d);
int64_t fadlev4d(int64_t* A, int64_t *M);

#endif /* L2lu_h */
