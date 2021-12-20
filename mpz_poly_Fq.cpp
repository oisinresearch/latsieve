#include "mpz_poly_Fq.h"

void mpz_poly_Fq_factor_edf(int d, mpz_poly_bivariate f0, int df0,
	mpz_poly h, int dh, int64_t q)
{
	// compute r = (q^2 - 1)/2
	int64_t r = (q*q - 1)/2;

	// compute random poly of degree df0 - 1
	mpz_poly_bivariate A;
	mpz_poly_bivariate_init(A, 0);
	mpz_poly A_i;
	mpz_poly_init(A_i, 0);
	for (int i = 0; i < df0-1; i++) {
		for (int j = 0; j < dh; j++) {
			mpz_poly_setcoeff_ui(A_i, j, rand() % q);
		}
	}
	int dA = df0 - 1;

	int L = 0;
	int64_t s = r;
	while (s>>=1) L++;	// number of bits in r 

	mpz_poly_bivariate B;
	mpz_poly_bivariate_init(B);
	mpz_poly_bivariate_set(B, A, dA);

	// using square and multiply, compute B = A^r - 1 mod f0, q, h
	for (int i = 2; i <= L; i++) {	
		// square
		mpz_poly_bivariate_mul(B, B, B);
		// reduce mod f0, q, h
		mpz_poly_Fq_mod(B, B, f0, q, h);
		if (r & (1<<(L - i)) == 1) {
			// multiply
			mpz_poly_bivariate_mul(B, B, A);
			// reduce mod f0, q, h
			mpz_poly_Fq_mod(B, B, f0, q, h);
		}
	}
	mzp_poly B_0; mpz_poly_init(B_0);
	for (i = 0; i < dh; i++) {
		mpz_poly_setcoeff(B_0, i, B->coeff[0]->coeff[i]);
	}
	int64_t B_0_0 = (mpz_get_ui(B_0->coeff[0]) - 1) % q;
	mpz_poly_setcoeff(B_0, 0, B_0_0);
	mpz_poly_bivariate_setcoeff(B, 0, B_0);

	// compute gcd(B, f0) mod q, h

	mpz_poly_clear(B_0);
	mpz_poly_bivariate_clear(B);
	mpz_poly_clear(A_i);
	mpz_poly_bivariate_clear(A);
}

