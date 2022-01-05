#include <gmp.h>

#ifndef MPZ_POLY_FQ_H
#define MPZ_POLY_FQ_H

void mpz_poly_Fq_factor_edf(int d, mpz_poly_bivariate f0, int64_t q, mpz_poly h,
	mpz_poly_bivariate* factors);
void mpz_poly_Fq_divrem(mpz_poly_bivariate Q, mpz_poly_bivariate R,
	mpz_poly_bivariate A, mpz_poly_bivariate B, int64_t q, mpz_poly h);
void mpz_poly_Fq_mul_Fq_shift(mpz_poly_bivariate T, mpz_poly_bivariate B, mpz_poly s, int d,
	int64_t q, mpz_poly h);
void mpz_poly_Fq_mod(mpz_poly_bivariate R, mpz_poly_bivariate A, mpz_poly_bivariate B,
    int64_t q, mpz_poly h);
int mpz_poly_Fq_mul(mpz_poly_bivariate f, mpz_poly_bivariate u, mpz_poly_bivariate v,
    int64_t q, mpz_poly h);
void mpz_poly_Fq_add(mpz_poly_bivariate f, mpz_poly_bivariate u, mpz_poly_bivariate v,
    int64_t q, mpz_poly h);
void mpz_poly_Fq_sub(mpz_poly_bivariate f, mpz_poly_bivariate u, mpz_poly_bivariate v,
    int64_t q, mpz_poly h);
void mpz_poly_Fq_makemonic(mpz_poly_bivariate G, int64_t q, mpz_poly h);
void mpz_poly_Fq_inv(mpz_poly T, mpz_poly A, int64_t q, mpz_poly h);
void mpz_poly_inv_Fq (mpz_poly f, mpz_poly u, mpz_poly h, int64_t p0);
void mpz_poly_Fq_gcd(mpz_poly_bivariate G, mpz_poly_bivariate u, mpz_poly_bivariate v,
	int64_t q, mpz_poly h);
int64_t min(int64_t a, int64_t b);
int64_t max(int64_t a, int64_t b);
void polyprintf(mpz_poly_bivariate f);
#endif /* MPZ_POLY_FQ_H */

