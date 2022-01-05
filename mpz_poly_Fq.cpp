#include "mpz_poly_bivariate.h"
#include <gmp.h>
#include "mpz_poly_Fq.h"
#include "mpz_poly_bivariate.h"
#include <cstdint>	// int64_t
#include <stack>
#include <stdlib.h>	// rand()
#include <iostream>	// cout

using std::stack;
using std::cout;
using std::endl;

void mpz_poly_Fq_factor_edf(int d, mpz_poly_bivariate f0, int64_t q, mpz_poly h,
	mpz_poly_bivariate* factors)
{
	// compute r = (q^2 - 1)/2
	int64_t r = (q*q - 1)/2;

	mpz_poly_bivariate A;
	mpz_poly_bivariate_init(A, 0);
	mpz_poly_bivariate Q;
	mpz_poly_bivariate_init(Q, 0);
	mpz_poly_bivariate R;
	mpz_poly_bivariate_init(R, 0);
	mpz_poly_bivariate G;
	mpz_poly_bivariate_init(G, 0);
	mpz_poly A_i;
	mpz_poly_init(A_i, 0);
	mpz_poly_bivariate B;
	mpz_poly_bivariate_init(B, 0);
	mpz_poly_bivariate_set(B, A);
	mpz_poly B_0; mpz_poly_init(B_0, 0);

	int m = 0;
	stack<mpz_poly_bivariate_struct_t*> composites;
	composites.push(f0);

	while (!composites.empty()) {
		// std::stack doesn't like mpz_poly_bivariate so we give it mpz_poly_bivariate_struct_t
		mpz_poly_bivariate_struct_t* f = composites.top();

		// compute random poly of degree f->deg - 1
		mpz_poly_bivariate_setzero(A);
		mpz_poly_bivariate_realloc(A, f->deg_y + 1);
		for (int i = 0; i <= f->deg_y - 1; i++) {
			mpz_poly_realloc(A_i, h->deg + 1);
			for (int j = 0; j < h->deg; j++) {
				int64_t rr = rand() % q;
				mpz_poly_setcoeff_ui(A_i, j, rr);
			}
			mpz_poly_bivariate_setcoeff(A, i, A_i);
		}

		int L = 1;
		int64_t s = r;
		while (s>>=1) L++;	// number of bits in r 
	
		mpz_poly_bivariate_set(B, A);
		// using square and multiply, compute B = A^r - 1 mod f, q, h
		for (int i = 2; i <= L; i++) {	
			// square
			mpz_poly_Fq_mul(B, B, B, q, h);
			// reduce mod f, q, h
			mpz_poly_Fq_divrem(Q, R, B, f, q, h);
			mpz_poly_bivariate_set(B, R);
			//polyprintf(B);
			if (r & (1<<(L - i))) {
				// multiply
				mpz_poly_Fq_mul(B, B, A, q, h);
				// reduce mod f, q, h
				mpz_poly_Fq_divrem(Q, R, B, f, q, h);
				mpz_poly_bivariate_set(B, R);
				//cout << "\t"; polyprintf(B);
			}
		}
		for (int i = 0; i <= B->coeff[0]->deg; i++) {
			mpz_poly_setcoeff(B_0, i, B->coeff[0]->coeff[i]);
		}
		int64_t B_0_0 = (mpz_get_ui(B_0->coeff[0]) - 1) % q;
		mpz_poly_setcoeff_ui(B_0, 0, B_0_0);
		mpz_poly_bivariate_setcoeff(B, 0, B_0);

		// compute gcd(B, f) mod q, h
		mpz_poly_Fq_gcd(G, B, f, q, h);
		
		if (0 < G->deg_y && G->deg_y < f->deg_y) {
			composites.pop();
			// we have a successful split
			mpz_poly_Fq_divrem(Q, R, f, G, q, h);
			// check if we have target degree d
			if (G->deg_y == d) {
				mpz_poly_Fq_makemonic(G, q, h);
				mpz_poly_bivariate_set(factors[m++], G);
			}
			else if (G->deg_y > 0) {
				composites.push(G);
			}
			if (Q->deg_y == d) {
				mpz_poly_Fq_makemonic(Q, q, h);
				mpz_poly_bivariate_set(factors[m++], Q);
			}
			else if (Q->deg_y > 0) {
				composites.push(Q);
			}			
		}
	}

	mpz_poly_clear(B_0);
	mpz_poly_bivariate_clear(B);
	mpz_poly_clear(A_i);
	mpz_poly_bivariate_clear(G);
	mpz_poly_bivariate_clear(R);
	mpz_poly_bivariate_clear(Q);
	mpz_poly_bivariate_clear(A);
}

// Euclidean division in F_q[x] where F_q = F_p[y]/<h>.
// As described in Cohen 3.13 Division of Polynomials
void mpz_poly_Fq_divrem(mpz_poly_bivariate Q, mpz_poly_bivariate R,
	mpz_poly_bivariate A, mpz_poly_bivariate B, int64_t q, mpz_poly h)
{
	mpz_poly LR; mpz_poly_init(LR, 0);
	mpz_poly LB; mpz_poly_init(LB, 0);
	mpz_poly t; mpz_poly_init(t, 0);
	mpz_poly_bivariate S; mpz_poly_bivariate_init(S, 0);
	mpz_poly_bivariate T; mpz_poly_bivariate_init(T, 0);
	mpz_poly_set(LB, B->coeff[B->deg_y]);
	mpz_poly ILB; mpz_poly_init(ILB, 0);
	mpz_poly_Fq_inv(ILB, LB, q, h); // ILB = 1/LB mod h, q

	if (mpz_poly_bivariate_iszero(A))
		mpz_poly_bivariate_setzero(R);
	else
		mpz_poly_bivariate_set(R, A);
	mpz_poly_bivariate_setzero(Q);
	mpz_poly_bivariate_setzero(S);

	while (R->deg_y >= B->deg_y) {
		mpz_poly_set(LR, R->coeff[R->deg_y]);
		mpz_poly_mul_mod_f_mod_ui(t, LR, ILB, h, q);
		mpz_poly_bivariate_setzero(S);
		int d = R->deg_y - B->deg_y;
		mpz_poly_bivariate_setcoeff(S, d, t);
		mpz_poly_Fq_add(Q, Q, S, q, h);
		//mpz_poly_Fq_mul(T, S, B, q, h);	// not fully necessary
		mpz_poly_Fq_mul_Fq_shift(T, B, S->coeff[d], d, q, h);
		mpz_poly_Fq_sub(R, R, T, q, h);
	}

	mpz_poly_clear(ILB);
	mpz_poly_bivariate_clear(T);
	mpz_poly_bivariate_clear(S);
	mpz_poly_clear(t);
	mpz_poly_clear(LB);
	mpz_poly_clear(LR);
}

void mpz_poly_Fq_mul_Fq_shift(mpz_poly_bivariate T, mpz_poly_bivariate B, mpz_poly s, int d,
	int64_t q, mpz_poly h)
{
	mpz_poly B_i; mpz_poly_init(B_i, 0);
	mpz_poly_bivariate_setzero(T);
	for (int i = 0; i <= B->deg_y; i++) {
		mpz_poly_set(B_i, B->coeff[i]);
		mpz_poly_cleandeg(B_i, B->coeff[i]->deg);
		mpz_poly_mul_mod_f_mod_ui(B_i, B_i, s, h, q);
		mpz_poly_bivariate_setcoeff(T, i + d, B_i);		
	}
	mpz_poly_clear(B_i);
}

void mpz_poly_Fq_mod(mpz_poly_bivariate R, mpz_poly_bivariate A, mpz_poly_bivariate B,
	int64_t q, mpz_poly h)
{
	mpz_poly_bivariate Q; mpz_poly_bivariate_init(Q, 0);
	mpz_poly_Fq_divrem(Q, R, A, B, q, h);
	mpz_poly_bivariate_clear(Q);
}

// See Cohen Algorithm 3.2.1 (Polynomial GCD)
void mpz_poly_Fq_gcd(mpz_poly_bivariate G, mpz_poly_bivariate u, mpz_poly_bivariate v,
	int64_t q, mpz_poly h)
{
	mpz_poly_bivariate A; mpz_poly_bivariate_init(A, 0);
	mpz_poly_bivariate B; mpz_poly_bivariate_init(B, 0);
	mpz_poly_bivariate Q; mpz_poly_bivariate_init(Q, 0);
	mpz_poly_bivariate R; mpz_poly_bivariate_init(R, 0);

	if (mpz_poly_bivariate_iszero(v)) {
		mpz_poly_bivariate_set(G, v);
	}
	else {
		mpz_poly_bivariate_set(B, v);		// this doesn't work if v == 0
		mpz_poly_bivariate_setzero(G);
		if (!mpz_poly_bivariate_iszero(u))
			mpz_poly_bivariate_set(A, u);
		
		while (!mpz_poly_bivariate_iszero(B)) {
			mpz_poly_Fq_divrem(Q, R, A, B, q, h);
			mpz_poly_bivariate_set(A, B);
			if (!mpz_poly_bivariate_iszero(R))
				mpz_poly_bivariate_set(B, R);
			else
				mpz_poly_bivariate_setzero(B);
		}
		if (!mpz_poly_bivariate_iszero(A))
			mpz_poly_bivariate_set(G, A);
	}

	mpz_poly_bivariate_clear(R);
	mpz_poly_bivariate_clear(Q);
	mpz_poly_bivariate_clear(B);
	mpz_poly_bivariate_clear(A);
}

/* Compute f=u*v mod q, h, where u has degree r, and v has degree s. */
int mpz_poly_Fq_mul(mpz_poly_bivariate f, mpz_poly_bivariate u, mpz_poly_bivariate v,
	int64_t q, mpz_poly h)
{
	mpz_poly T; mpz_poly_init(T, 0);
	mpz_poly f_t; mpz_poly_init(f_t, 0);
	mpz_poly_bivariate A; mpz_poly_bivariate_init(A, 0);
	mpz_poly_bivariate B; mpz_poly_bivariate_init(B, 0);
	mpz_poly_bivariate_set(A, u);
	mpz_poly_bivariate_set(B, v);
	int r = A->deg_y; int s = B->deg_y; int t = h->deg;

	mpz_poly_bivariate_setzero(f);
	mpz_poly_bivariate_realloc(f, r + s + 1);
	for (int i = 0; i <= r + s; i++)
		mpz_poly_set_zero(f->coeff[i]);
	for (int i = 0; i <= r; ++i) {
		for (int j = 0; j <= s; ++j) {
			mpz_poly_mul_mod_f_mod_ui(T, A->coeff[i], B->coeff[j], h, q);
			mpz_poly_add(f_t, f->coeff[i+j], T);
			// now reduce f_t mod q
			for (int k = 0; k <= f_t->deg; k++)
				mpz_poly_setcoeff_ui(f_t, k, (mpz_poly_getcoeff_si(f_t, k) + q) % q);				
			mpz_poly_bivariate_setcoeff(f, i+j, f_t);
		}
	}
	mpz_poly_bivariate_clear(B);
	mpz_poly_bivariate_clear(A);
	mpz_poly_clear(T);
	return r + s;
}

// Compute f = u + v mod q, h
void mpz_poly_Fq_add(mpz_poly_bivariate f, mpz_poly_bivariate u, mpz_poly_bivariate v,
	int64_t q, mpz_poly h)
{
	int r = u->deg_y; int s = v->deg_y; int t = max(r, s);
	mpz_poly f_i; mpz_poly_init(f_i, 0);
	mpz_poly_bivariate_realloc(f, t + 1);
	mpz_poly_bivariate_realloc(u, t + 1);
	mpz_poly_bivariate_realloc(v, t + 1);

	for (int i = 0; i <= t; i++) {
		mpz_poly_set_zero(f_i);
		mpz_poly_add(f_i, u->coeff[i], v->coeff[i]);
		mpz_poly_mod_ui(f->coeff[i], f_i, q);
	}

	//mpz_poly_bivariate_fixifzero(f);
	//if (!mpz_poly_bivariate_iszero)
	mpz_poly_bivariate_cleandeg(f, t);
	mpz_poly_clear(f_i);
}

// Compute f = u - v mod q, h
void mpz_poly_Fq_sub(mpz_poly_bivariate f, mpz_poly_bivariate u, mpz_poly_bivariate v,
	int64_t q, mpz_poly h)
{
	int r = u->deg_y; int s = v->deg_y; int t = max(r, s);
	mpz_poly f_i; mpz_poly_init(f_i, 0);
	mpz_poly_bivariate_realloc(f, t + 1);
	mpz_poly_bivariate_realloc(u, t + 1);
	mpz_poly_bivariate_realloc(v, t + 1);

	for (int i = 0; i <= t; i++) {
		mpz_poly_set_zero(f_i);
		mpz_poly_sub(f_i, u->coeff[i], v->coeff[i]);
		mpz_poly_mod_ui(f->coeff[i], f_i, q);
	}

	//mpz_poly_bivariate_fixifzero(f);
	//if (!mpz_poly_bivariate_iszero)
	mpz_poly_bivariate_cleandeg(f, t);
	mpz_poly_clear(f_i);
}

void mpz_poly_Fq_makemonic(mpz_poly_bivariate G, int64_t q, mpz_poly h)
{
	mpz_poly T; mpz_poly_init(T, 0);
	mpz_poly_Fq_inv(T, G->coeff[G->deg_y], q, h);
	for (int i = 0; i <= G->deg_y; i++) {
		mpz_poly_mul_mod_f_mod_ui(G->coeff[i], G->coeff[i], T, h, q);
	}
	mpz_poly_clear(T);
}

// Note A should be nonzero
void mpz_poly_Fq_inv(mpz_poly T, mpz_poly A, int64_t q, mpz_poly h)
{
	mpz_poly d; mpz_poly_init(d, 0);
	mpz_t p; mpz_init_set_ui(p, q);
	mpz_poly v; mpz_poly_init(v, 0);
	/* computes 1 = gcd(A, h) = T*A + v*h mod p, with p in mpz_t */
	mpz_poly_xgcd_mpz(d, A, h, T, v, p);
	mpz_poly_clear(v);
	mpz_clear(p);
	mpz_poly_clear(d);
}

int64_t min(int64_t a, int64_t b)
{
	return (a < b ? a : b);
}

int64_t max(int64_t a, int64_t b)
{
	return (a > b ? a : b);
}

void polyprintf(mpz_poly_bivariate f)
{
  if (f->deg_y == -1) {
    cout << "0" << endl;
    return;
  }
  for (int i = 0, printed = 0; i <= f->deg_y; ++i) {
    if (f->coeff[i]->deg == -1) {
      continue;
    }

    if (printed++) {
      cout << "+";
    }

    cout << "(";
	for (int j = 0, p2 = 0; j <= f->coeff[i]->deg; ++j) {
		if (p2++) { cout << "+"; }
		int k = mpz_cmp_ui(f->coeff[i]->coeff[j], 0);
		if (k != 0)
			cout << mpz_get_str(NULL, 10, f->coeff[i]->coeff[j]);
		if (j) { cout << "*y"; }
		if (j > 1) { cout << "^" << j; }
	}
    //mpz_poly_fprintf_endl(fp, f->coeff[i], 0);
    cout << ")";
    if (i) {
      cout << "*x";
    }

    if (i > 1) {
      cout << "^" << i;
    }
  }
  cout << endl;
}

