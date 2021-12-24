#ifndef MPZ_POLY_BIVARIATE_H
#define MPZ_POLY_BIVARIATE_H 

#include <stdio.h>
#include <gmp.h>
#include "mpz_poly.h"

/*
 * Bivariate polynomial are written as polynomial in y with coefficients in
 * mpz_poly in x.
 */

typedef struct {
  int alloc;
  int deg_y;
  int deg_x;
  mpz_poly * coeff;
} mpz_poly_bivariate_struct_t;

typedef mpz_poly_bivariate_struct_t mpz_poly_bivariate[1];

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Initialize a mpz_poly_bivariate with alloc = d + 1.
 */
void mpz_poly_bivariate_init(mpz_poly_bivariate f, int d);
/*
 * Initialize a mpz_poly_bivariate with alloc = d + 1 and alloc of all mpz_poly
 * equals to dx + 1.
 */
void mpz_poly_bivariate_init_y_x(mpz_poly_bivariate f, int dy, int dx);
/*
 * Increase alloc.
 */
void mpz_poly_bivariate_realloc(mpz_poly_bivariate f, int nc);
/*
 * Increase alloc and set alloc of all mpz_poly equals dx + 1.
 */
void mpz_poly_bivariate_realloc_x(mpz_poly_bivariate f, int nc, int dx);
/*
 * Clear a mpz_poly_bivariate.
 */
void mpz_poly_bivariate_clear(mpz_poly_bivariate f);
/*
 * Find deg_y of f.
 */
void mpz_poly_bivariate_cleandeg(mpz_poly_bivariate f, int deg_y);
/*
 * Set the ith coefficient of f.
 */
void mpz_poly_bivariate_setcoeff(mpz_poly_bivariate f, int i,
    mpz_poly z);
/*
 * Set the ith coefficient of f to an unisigned int constant polynomial.
 */
void mpz_poly_bivariate_setcoeff_ui(mpz_poly_bivariate f, int i,
    uint64_t z)
/*
 * Print f in a file.
 */
void mpz_poly_bivariate_fprintf(FILE * fp, mpz_poly_bivariate f);
/*
 * Compute res(x) = f(x, y).
 */
void mpz_poly_bivariate_eval_y(mpz_poly res,
    mpz_poly_bivariate f, mpz_t y);
/*
 * Compute res(y) = f(x, y).
 */
void mpz_poly_bivariate_eval_x(mpz_poly res,
    mpz_poly_bivariate f, mpz_t x);

/*
 * Compute resultant(x) = res(f, g).
 */
void mpz_poly_bivariate_resultant_y(mpz_poly resultant,
    mpz_poly_bivariate f, mpz_poly_bivariate g);

/*
 * Compute resultant(y) = res(f, g).
 */
void mpz_poly_bivariate_resultant_x(mpz_poly resultant,
    mpz_poly_bivariate f, mpz_poly_bivariate g);


#ifdef __cplusplus
}
#endif

#endif /* MPZ_POLY_BIVARIATE_H */
