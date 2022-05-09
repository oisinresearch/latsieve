#include <cstdint>  // int64_t
#include <boost/multiprecision/cpp_int.hpp>

#ifndef L2lu128_h
#define L2lu128_h

using namespace boost::multiprecision;

void int128L2(int128_t* b, int d);
void matprint(int d, int128_t* M);

#endif /* L2lu128_h */
