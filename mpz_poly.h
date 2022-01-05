#include <gmp.h>
#include <cstdint>

#ifndef MPZ_POLY_H
#define MPZ_POLY_H
/* Note, deg = -1 means P=0; otherwise, one should have coeff[deg] != 0.
   Warning: a polynomial of degree d needs d+1 allocation. */
struct mpz_poly_s {
	int alloc;
	int deg;
	mpz_t* coeff;
};

typedef struct mpz_poly_s mpz_poly[1];

/* this structure represents the interval [a/2^ka, b/2^kb] */
typedef struct {
    mpz_t a;
    int ka;
    mpz_t b;
    int kb;
} usp_root_data;

void mpz_poly_resultant(mpz_t res, mpz_poly f, mpz_poly g);
void mpz_poly_init(mpz_poly f, int d);
void mpz_poly_clear(mpz_poly f);
mpz_t* mpz_poly_lc(mpz_poly f);
void mpz_poly_set(mpz_poly f, mpz_poly g);
void mpz_poly_getcoeff(mpz_t res, int i, mpz_poly f);
int64_t mpz_poly_getcoeff_si(mpz_poly f, int i);
void mpz_poly_set_xi(mpz_poly f, int i);
void mpz_poly_set_mpz(mpz_poly f, mpz_t* g, int degg);
void mpz_poly_content(mpz_t* c, mpz_poly F);
void mpz_poly_divexact(mpz_poly Q, mpz_poly P, mpz_t a);
static int mpz_poly_divexact_clobber (mpz_poly q, mpz_poly h, mpz_poly f,
                   mpz_t p);
int mpz_poly_divexact_modp (mpz_poly q, mpz_poly h, mpz_poly f, mpz_t p);
void mpz_poly_swap(mpz_poly f, mpz_poly g);
void printmpz(mpz_t x);
void mpz_poly_pseudo_remainder(mpz_poly r, mpz_poly a, mpz_poly b);
void mpz_poly_mul_mpz(mpz_poly Q, mpz_poly P, mpz_t a);
void mpz_poly_mul(mpz_poly f, mpz_poly g, mpz_poly h);
void mpz_poly_realloc(mpz_poly f, int nc);
void mpz_poly_set_zero(mpz_poly f);
void mpz_poly_setcoeff(mpz_poly f, int i, mpz_t z);
void mpz_poly_setcoeff_ui(mpz_poly f, int i, unsigned int z);
void mpz_poly_setcoeff_si(mpz_poly f, int i, int64_t z);
void mpz_poly_cleandeg(mpz_poly f, int d);
static inline int mpz_poly_coeff_sgn(mpz_poly f, int i);
void mpz_poly_add(mpz_poly f, mpz_poly g, mpz_poly h);
void mpz_poly_sub(mpz_poly f, mpz_poly g, mpz_poly h);
void mpz_poly_setcoeff_double(mpz_poly f, int i, double z);
int numberOfRealRoots(mpz_poly f, usp_root_data* Roots);
void divide(mpz_t a, int k, mpz_poly f);
int usp(mpz_t a, mpz_t b, int m, int up, int va, int vb, int *nroots,
     mpz_poly f, mpz_poly r, usp_root_data *R);
double ln2(mpz_t a);
void getroot(mpz_t a, int ka, mpz_t b, int kb, int *nroots, usp_root_data *R);
int signValue(mpz_t a, int k, mpz_poly f);
int sign (mpz_t a);
void mpz_poly_discriminant(mpz_t disc, mpz_poly f);
void mpz_poly_derivative(mpz_poly df, mpz_poly f);
void mpz_poly_eval(mpz_t res, mpz_poly f, mpz_t x);
void mpz_poly_eval_ui (mpz_t res, mpz_poly f, unsigned long x);

void mpz_poly_makemonic_mod_mpz (mpz_poly Q, mpz_poly P, mpz_t m);
int mpz_poly_mod_mpz (mpz_poly R, mpz_poly A, mpz_t m);
int mpz_poly_mod_ui(mpz_poly R, mpz_poly A, uint64_t m);
static void mpz_poly_gcd_mpz_clobber (mpz_poly f, mpz_poly g, mpz_t p);
void mpz_poly_gcd_mpz (mpz_poly f, mpz_poly a, mpz_poly b, mpz_t p);

struct mpz_poly_with_m_s {
    mpz_poly f;
    int m;
};
typedef struct mpz_poly_with_m_s mpz_poly_with_m[1];
typedef struct mpz_poly_with_m_s * mpz_poly_with_m_ptr;
typedef const struct mpz_poly_with_m_s * mpz_poly_with_m_srcptr;

struct mpz_poly_factor_list_s {
    mpz_poly_with_m * factors;
    int alloc;
    int size;
};
typedef struct mpz_poly_factor_list_s mpz_poly_factor_list[1];
typedef struct mpz_poly_factor_list_s * mpz_poly_factor_list_ptr;
typedef const struct mpz_poly_factor_list_s * mpz_poly_factor_list_srcptr;

void mpz_poly_factor_list_init(mpz_poly_factor_list_ptr l);
void mpz_poly_factor_list_clear(mpz_poly_factor_list_ptr l);
void mpz_poly_factor_list_flush(mpz_poly_factor_list_ptr l);
void mpz_poly_factor_list_push(mpz_poly_factor_list_ptr l, mpz_poly f, int m);
//void mpz_poly_factor_list_fprintf(FILE* fp, mpz_poly_factor_list_srcptr l);
int mpz_poly_factor_sqf(mpz_poly_factor_list_ptr lf, mpz_poly f, mpz_t p);
int mpz_poly_factor_ddf(mpz_poly_factor_list_ptr lf, mpz_poly f0, mpz_t p);
int mpz_poly_factor_edf(mpz_poly_factor_list_ptr lf, mpz_poly f, int k, mpz_t p, gmp_randstate_t rstate);

/* output is sorted by degree and lexicographically */
int mpz_poly_factor(mpz_poly_factor_list lf, mpz_poly f, mpz_t p, gmp_randstate_t rstate);
int mpz_poly_is_irreducible(mpz_poly f, mpz_t p);

void mpz_ndiv_r (mpz_t a, mpz_t b, mpz_t c);

void mpz_poly_pow_mod_f_mod_mpz (mpz_poly Q, mpz_poly P, mpz_poly f, mpz_t a, mpz_t p);
void mpz_poly_mul_mod_f_mod_ui (mpz_poly Q, mpz_poly P1, mpz_poly P2,
			    mpz_poly f, int64_t m0);
void mpz_poly_mul_mod_f_mod_mpz (mpz_poly Q, mpz_poly P1, mpz_poly P2, mpz_poly f,
	mpz_t m, mpz_t invf);
void mpz_poly_mul_mod_f (mpz_poly Q, mpz_poly P1, mpz_poly P2, mpz_poly f);
void mpz_poly_sqr_mod_f_mod_mpz (mpz_poly Q, mpz_poly P, mpz_poly f,
			    mpz_t m, mpz_t invf);
void mpz_poly_pow_ui_mod_f_mod_mpz (mpz_poly Q, mpz_poly P, mpz_poly f,
                          unsigned long a, mpz_t p);
int mpz_poly_mod_f_mod_mpz (mpz_poly R, mpz_poly f, mpz_t m, mpz_t invf);
static int mpz_poly_sqr (mpz_t *f, mpz_t *g, int r);
int mpz_poly_div_qr_z (mpz_poly q, mpz_poly r, mpz_poly f, mpz_poly g);
int mpz_poly_div_r_z (mpz_poly r, mpz_poly f, mpz_poly g);
void mpz_poly_add_ui (mpz_poly g, mpz_poly f, unsigned long a);
void mpz_poly_sub_ui (mpz_poly g, mpz_poly f, unsigned long a);
int mpz_poly_cmp (mpz_poly a, mpz_poly b);
size_t mpz_poly_size (mpz_poly f);
void mpz_poly_div_xi(mpz_poly g, mpz_poly f, int i);
void mpz_poly_mul_xi (mpz_poly g, mpz_poly f, int i);
int mpz_poly_valuation(mpz_poly f);
void mpz_poly_xgcd_mpz (mpz_poly d, mpz_poly f, mpz_poly g, mpz_poly u, mpz_poly v, mpz_t p);
void mpz_poly_sub_mod_mpz(mpz_poly f, mpz_poly g, mpz_poly h, mpz_t m);
#endif /* MPZ_POLY_H */

