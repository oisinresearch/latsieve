#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mpz_poly.h"
#include <algorithm>
#include <iostream> // cout
#include <cassert> // assert
#include <cstdint>

using std::cout;
using std::endl;
using std::flush;
using std::max;

/*
 * Compute the resultant of f and g and set the resultat in res.
 *  See Henri Cohen, "A Course in Computational Algebraic Number Theory",
 *  for more information.
 *
 * Assume that the polynomials are normalized.
 */
void mpz_poly_resultant(mpz_t res, mpz_poly f, mpz_poly g)
{
	if (f->deg == -1 || g->deg == -1) {
		mpz_set_ui(res, 0);
		return;
	}

	//ASSERT(mpz_cmp_ui(f->coeff[f->deg], 0) != 0);
	//ASSERT(mpz_cmp_ui(g->coeff[g->deg], 0) != 0);

	long int s = 1;
	mpz_t k;
	mpz_t h;
	mpz_t t;
	mpz_t tmp;
	int d;
	mpz_poly r;
	mpz_poly a;
	mpz_poly b;

	mpz_init(k);
	mpz_init(h);
	mpz_init(t);
	mpz_init(tmp);
	mpz_poly_init(r, -1);
	mpz_poly_init(a, f->deg);
	mpz_poly_init(b, g->deg);

	mpz_poly_set(a, f);
	mpz_poly_set(b, g);

	mpz_poly_content(&k, a);
	mpz_poly_content(&h, b);

	mpz_poly_divexact(a, a, k);
	mpz_poly_divexact(b, b, h);

	#ifndef NDEBUG
	mpz_poly_content(&tmp, a);
	//ASSERT(mpz_cmp_ui(tmp, 1) == 0);
	mpz_poly_content(&tmp, b);
	//ASSERT(mpz_cmp_ui(tmp, 1) == 0);
	#endif // NDEBUG

	mpz_pow_ui(t, k, (unsigned long int) b->deg);
	mpz_pow_ui(tmp, h, (unsigned long int) a->deg);
	mpz_mul(t, t, tmp);

	mpz_set_ui(k, 1);
	mpz_set_ui(h, 1);

	if (a->deg < b->deg) {
		mpz_poly_swap(a, b);

		if ((a->deg % 2) == 1 && (b->deg % 2) == 1) {
			s = -1;
		}
	}

	while (b->deg > 0) {
		d = a->deg - b->deg;

		if ((a->deg % 2) == 1 && (b->deg % 2) == 1) {
			s = -s;
		}

		mpz_poly_pseudo_remainder(r, a, b);
		mpz_poly_set(a, b);

		//ASSERT(d >= 0);

		mpz_pow_ui(tmp, h, (unsigned long int) d);
		mpz_mul(tmp, k, tmp);

		mpz_poly_divexact(b, r, tmp);

		mpz_set(k, *mpz_poly_lc(a));

		#ifdef NDEBUG
		if (d == 0) {
			//ASSERT(mpz_cmp_ui(h, 1) == 0);
		}
		#endif // NDEBUG
		mpz_pow_ui(h, h, (unsigned long int) (d - 1));
		mpz_pow_ui(tmp, k, (unsigned long int) d);
		mpz_divexact(h, tmp, h);
	}

	//Prevent an error if b = 0.
	if (b->deg == -1) {
		mpz_set_ui(res, 0);
	}
	else {
		//ASSERT(a->deg > 0);
		//ASSERT(b->deg == 0);

		mpz_pow_ui(h, h, (unsigned long int) (a->deg - 1));

		//ASSERT(a->deg >= 0);

		mpz_pow_ui(tmp, b->coeff[0], (unsigned long int) a->deg);
		mpz_divexact(h, tmp, h);

		mpz_mul_si(t, t, s);
		mpz_mul(h, h, t);
		mpz_set(res, h);
	}

	mpz_clear(k);
	mpz_clear(h);
	mpz_clear(t);
	mpz_clear(tmp);
	mpz_poly_clear(a);
	mpz_poly_clear(b);
	mpz_poly_clear(r);
}

/* Allocate a polynomial that can contain 'd+1' coefficients and set to zero.
   We allow d < 0, which is equivalent to d = -1.
 */
void mpz_poly_init(mpz_poly f, int d)
{
	f->deg = -1;
	if (d < 0) {
		f->alloc = 0;
		f->coeff = (mpz_t*) NULL;
	}
	else {
		int i;
		f->alloc = d+1;
		f->coeff = (mpz_t*)malloc((d+1)*sizeof(mpz_t));
		//FATAL_ERROR_CHECK ((*f)->coeff == NULL, "not enough memory");
		for (i = 0; i <= d; ++i)
			mpz_init(f->coeff[i]);
	}
}

/* Free polynomial f in mpz_poly. */
void mpz_poly_clear(mpz_poly f) 
{
	int i;
	for (i = 0; i < f->alloc; ++i)
		mpz_clear(f->coeff[i]);
	if (f->coeff != NULL)
		free(f->coeff);
	f->coeff = NULL; /* to avoid a double-free */
	memset(f, 0, sizeof(mpz_poly));
	f->deg = -1;
	f->alloc = 0; /* to avoid a double-free */
}

mpz_t* mpz_poly_lc(mpz_poly f)
{
    //assert(f->deg >= 0);
    return &(f->coeff[f->deg]);
}

/* Copy g to f, where f must be initialized (but there is no need it has
   enough allocated coefficients, since mpz_poly_setcoeff reallocates if
   needed). */
void mpz_poly_set(mpz_poly f, mpz_poly g) 
{
	if (f == g)
		return; /* nothing to do */

	f->deg = g->deg;
	for (int i = g->deg; i >= 0; --i)
		mpz_poly_setcoeff(f, i, g->coeff[i]);
}

/* Get coefficient for the i-th term. */
void mpz_poly_getcoeff(mpz_t res, int i, mpz_poly f)
{
    if (i > f->deg)
        mpz_set_ui (res, 0);
    else
        mpz_set (res, f->coeff[i]);
}

/* Get coefficient for the i-th term.  Note: f->coeff[i] must be an int64_t or equivalent. */
int64_t mpz_poly_getcoeff_si(mpz_poly f, int i)
{
    if (i > f->deg)
    	return 0;
    else
        return mpz_get_ui(f->coeff[i]);
}

/* x^i is often useful */
void mpz_poly_set_xi(mpz_poly f, int i)
{
    mpz_poly_realloc (f, i + 1);
    for(int j = 0 ; j <= i ; j++) {
        mpz_set_ui(f->coeff[j], j == i);
    }
    f->deg = i;
}

void mpz_poly_set_mpz(mpz_poly f, mpz_t* g, int degg) 
{
	f->deg = degg;
	for (int i = degg; i >= 0; --i)
		mpz_poly_setcoeff(f, i, g[i]);
}

/* put in c the content of f */
void mpz_poly_content(mpz_t* c, mpz_poly F)
{
	mpz_t* f = F->coeff;
	int d = F->deg;

	mpz_set(*c, f[0]);
	for (int i = 1; i <= d; i++)
		mpz_gcd(*c, *c, f[i]);
	mpz_abs(*c, *c);
}


/* Set Q=P/a, where a is an mpz_t. Assume a divides the content of P (the
   division is done with mpz_divexact). Otherwise the result is not correct. */
void mpz_poly_divexact(mpz_poly Q, mpz_poly P, mpz_t a)
{
	mpz_t aux;
	mpz_init (aux);
	Q->deg = P->deg;
	for (int i = 0; i <= P->deg; ++i)
	{
		mpz_divexact(aux, P->coeff[i], a);
		mpz_set(Q->coeff[i], aux);
	}
	mpz_clear (aux);
}

/* q=divexact(h, f) mod p, f not necessarily monic.
   Assumes lc(h) <> 0 mod p.
   Clobbers h. */
/* Coefficients of f must be reduced mod p
 * Coefficients of h need not be reduced mod p
 * Coefficients of q are reduced mod p
 */
static int mpz_poly_divexact_clobber (mpz_poly q, mpz_poly h, mpz_poly f,
                   mpz_t p)
{
  int i, d = f->deg, dh = h->deg;
  mpz_t t, aux;

  mpz_init (t);
  mpz_init (aux);
  assert (d >= 0);
  assert (dh >= 0);
  assert (dh >= d);
  assert (mpz_divisible_p (h->coeff[dh], p) == 0);

  mpz_poly_realloc (q, dh + 1 - d);
  q->deg = dh - d;
  /* t is 1/f[d] mod p */
  mpz_set (aux, f->coeff[d]);
  mpz_mod (aux, aux, p);
  if (!mpz_invert (t, aux, p)) {
      mpz_clear(t);
      mpz_clear(aux);
      return 0;
  }


  while (dh >= d) {

    /* subtract h[dh]/f[d]*x^(dh-d)*f from h */
    if (mpz_cmp_ui(t, 1) != 0) {
      mpz_mul (h->coeff[dh], h->coeff[dh], t);
      mpz_mod (h->coeff[dh], h->coeff[dh], p);
    }
    mpz_set (q->coeff[dh-d], h->coeff[dh]);
    mpz_mod (q->coeff[dh-d], q->coeff[dh - d], p);

    /* we only need to update the coefficients of degree >= d of h,
       i.e., we want i >= 2d - dh */
    for (i = (2 * d > dh) ? 2 * d - dh : 0; i < d; i++) {
      mpz_mul (aux, h->coeff[dh], f->coeff[i]);
      mpz_sub (h->coeff[dh - d + i], h->coeff[dh - d + i], aux);
      mpz_mod (h->coeff[dh - d + i], h->coeff[dh - d + i], p);
    }
    dh --;
  }
  /* since lc(h) <> 0 mod p, q is normalized */

  mpz_clear (t);
  mpz_clear (aux);
  return 1;
}

/* q <- divexact(h, f) mod p, f not necessarily monic. */
/* coefficients of f must be reduced mod p
 * coefficients of h need not be reduced mod p
 * coefficients of q are reduced mod p
 */
int mpz_poly_divexact_modp (mpz_poly q, mpz_poly h, mpz_poly f, mpz_t p)
{
    mpz_poly hh;
    mpz_poly_init(hh, h->deg);
    mpz_poly_set(hh, h);
    int r = mpz_poly_divexact_clobber(q, hh, f, p);
    mpz_poly_clear(hh);
    return r;
}


/* swap f and g */
void mpz_poly_swap(mpz_poly f, mpz_poly g)
{
	int i;
	mpz_t *t;

	i = f->alloc;
	f->alloc = g->alloc;
	g->alloc = i;
	i = f->deg;
	f->deg = g->deg;
	g->deg = i;
	t = f->coeff;
	f->coeff = g->coeff;
	g->coeff = t;
}


void printmpz(mpz_t x)
{
    cout << mpz_get_str(NULL, 10, x);
    cout << endl;
}


/*
 * Compute the pseudo remainder of a and b such that
 *  lc(b)^(deg(a) - deg(b) + 1) * a = b * q + r with deg(r) < deg(b).
 *  See Henri Cohen, "A Course in Computational Algebraic Number Theory",
 *  for more information.
 *
 * Assume that deg(a) >= deg(b) and b is not the zero polynomial.
 */
void mpz_poly_pseudo_remainder(mpz_poly r, mpz_poly a, mpz_poly b)
{
	//ASSERT(a->deg >= b->deg);
	//ASSERT(b->deg != -1);

	int m = a->deg;
	int n = b->deg;
	mpz_t d;
	int e;
	mpz_poly s;

	#ifndef NDEBUG
	/*MAYBE_UNUSED*/ mpz_poly q_tmp;
	#endif // NDEBUG

	mpz_init(d);
	mpz_set(d, *mpz_poly_lc(b));

	#ifndef NDEBUG
	mpz_poly_init(q_tmp, 0);
	#endif // NDEBUG

	mpz_poly_set(r, a);

	e = m - n + 1;

	while (r->deg >= n) {
		mpz_poly_init(s, r->deg - n);
		mpz_poly_setcoeff(s, r->deg - n, *mpz_poly_lc(r));

		#ifndef NDEBUG
        mpz_poly_mul_mpz(q_tmp, q_tmp, d);
        mpz_poly_add(q_tmp, s, q_tmp);
		#endif // NDEBUG

		mpz_poly_mul_mpz(r, r, d);
		mpz_poly_mul(s, b, s);
		mpz_poly_sub(r, r, s);
		mpz_poly_clear(s);
		e--;
	}

	//ASSERT(e >= 0);

	mpz_pow_ui(d, d, (unsigned long int) e);

	#ifndef NDEBUG
	mpz_poly_mul_mpz(q_tmp, q_tmp, d);
	#endif // NDEBUG

	mpz_poly_mul_mpz(r, r, d);

	#ifndef NDEBUG
	mpz_poly f, g;
	mpz_poly_init(f, a->deg);
	mpz_poly_init(g, b->deg);
	mpz_poly_set(f, a);
	mpz_poly_set(g, b);

	mpz_set(d, *mpz_poly_lc(g));

	//ASSERT(m - n + 1 >= 0);

	mpz_pow_ui(d, d, (unsigned long int) (m - n + 1));
	mpz_poly_mul_mpz(f, f, d);

	mpz_poly_mul(g, g, q_tmp);
	mpz_poly_add(g, g, r);

	//ASSERT(mpz_poly_cmp(f, g) == 0);

	mpz_poly_clear(f);
	mpz_poly_clear(g);
	mpz_poly_clear(q_tmp);
	#endif // NDEBUG

	mpz_clear(d);
}


/* Set Q=a*P, where a is an mpz_t */
void mpz_poly_mul_mpz(mpz_poly Q, mpz_poly P, mpz_t a)
{
	int i;
	mpz_t aux;

	mpz_init (aux);
	Q->deg = P->deg;
	for (i = 0; i <= P->deg; ++i)
	{
		mpz_mul(aux, P->coeff[i], a);
		mpz_poly_setcoeff(Q, i, aux);
	}
	mpz_clear (aux);
}


/* Set f=g*h. Note: f might equal g or h.
   Assumes the g[g->deg] and h[h->deg[] are not zero. */
void mpz_poly_mul(mpz_poly f, mpz_poly g, mpz_poly h)
{
	if (f == h || f == g)
	{
		mpz_poly aux;
		mpz_poly_init(aux, -1);
		mpz_poly_mul(aux, g, h);
		mpz_poly_set(f, aux);
		mpz_poly_clear(aux);
		return;
	}

	if ((g->deg == -1) || (h->deg == -1)) {
		f->deg = -1;
		return;
	}

	mpz_poly_realloc(f, g->deg + h->deg + 1);

	if (g->deg == 0) {
		mpz_poly_mul_mpz(f, h, g->coeff[0]);
		return;
	}

	if (h->deg == 0) {
		mpz_poly_mul_mpz(f, g, h->coeff[0]);
		return;
	}

	//ASSERT(mpz_cmp_ui(g->coeff[g->deg], 0) != 0);
	//ASSERT(mpz_cmp_ui(h->coeff[h->deg], 0) != 0);
	//ASSERT(f != g);
	//ASSERT(f != h);

	mpz_t t;
	mpz_init(t);

	mpz_poly_realloc(f, (g->deg) + (h->deg) + 1);
	for (int i = 0; i <= f->deg; i++) mpz_set_ui(f->coeff[i], 0);
	for (int i = 0; i <= g->deg; i++) {
		for (int j = 0; j <= h->deg; j++) {
			mpz_mul(t, g->coeff[i], h->coeff[j]);
			mpz_add(f->coeff[i+j], f->coeff[i+j], t);
		}
	}
	f->deg = g->deg + h->deg;

	/* there is no need to run mpz_poly_cleandeg since g[g->deg] <> 0
	 and h[h->deg] <> 0 */
	//ASSERT(mpz_cmp_ui(f->coeff[f->deg], 0) != 0);

	mpz_clear(t);
}


/* realloc f to (at least) nc coefficients */
void mpz_poly_realloc(mpz_poly f, int nc)
{
	if (f->alloc < nc) {
		f->coeff = (mpz_t*)realloc(f->coeff, nc*sizeof(mpz_t));
		//FATAL_ERROR_CHECK(f->coeff == NULL, "not enough memory");
		for (int i = f->alloc; i < nc; i++)
			mpz_init(f->coeff[i]);
		f->alloc = nc;
	}
}

/* Set a zero polynomial. */
void mpz_poly_set_zero(mpz_poly f)
{
  f->deg = -1;
}

/* Set mpz_t coefficient for the i-th term. */
void mpz_poly_setcoeff(mpz_poly f, int i, mpz_t z)
{
	mpz_poly_realloc(f, i + 1);
	mpz_set(f->coeff[i], z);
	if (i >= f->deg)
		mpz_poly_cleandeg(f, i);
}

/* Set mpz_t coefficient for the i-th term. */
void mpz_poly_setcoeff_ui(mpz_poly f, int i, unsigned int z)
{
	mpz_poly_realloc(f, i + 1);
	mpz_set_ui(f->coeff[i], z);
	if (i >= f->deg)
		mpz_poly_cleandeg(f, i);
}

/* Set mpz_t coefficient for the i-th term. */
void mpz_poly_setcoeff_si(mpz_poly f, int i, int64_t z)
{
	mpz_poly_realloc(f, i + 1);
	mpz_set_si(f->coeff[i], z);
	if (i >= f->deg)
		mpz_poly_cleandeg(f, i);
}

/* Find polynomial degree. */
void mpz_poly_cleandeg(mpz_poly f, int d)
{
	//ASSERT(d >= -1);
	while ((d >= 0) && (mpz_poly_coeff_sgn(f, d)==0))
		d--;
	f->deg = d;
}

/* Return 0 if f[i] is zero, -1 is f[i] is negative and +1 if f[i] is positive,
   like mpz_sgn function. */
inline int mpz_poly_coeff_sgn(mpz_poly f, int i)
{
	if (i >= f->alloc)
		return 0;
	else
		return mpz_sgn(f->coeff[i]);
}

/* Set f=g+h.
   Note: f can be the same as g or h;
         g can be the same as h. */
void mpz_poly_add(mpz_poly f, mpz_poly g, mpz_poly h) {
    int i, maxdeg;
    mpz_t z;
    mpz_init(z);
    maxdeg = max(g->deg, h->deg);
    mpz_poly_realloc(f, maxdeg + 1);
    for (i = 0 ; i <= maxdeg ; i++) {
        if (i <= g->deg)
            mpz_set(z, g->coeff[i]);
        else
            mpz_set_ui(z, 0);
        if (i <= h->deg)
            mpz_add(z, z, h->coeff[i]);
        mpz_set(f->coeff[i], z);
    }
    f->deg = maxdeg;
    mpz_clear(z);
    mpz_poly_cleandeg(f, maxdeg);
}

/* Set f=g-h.
   Note: f can be the same as g or h;
         g can be the same as h. */
void mpz_poly_sub(mpz_poly f, mpz_poly g, mpz_poly h)
{
    int i, maxdeg;
    mpz_t z;
    mpz_init(z);
    maxdeg = max(g->deg, h->deg);
    mpz_poly_realloc(f, maxdeg + 1);
    for (i = 0 ; i <= maxdeg ; i++) {
        if (i <= g->deg)
            mpz_set(z, g->coeff[i]);
        else
            mpz_set_ui(z, 0);
        if (i <= h->deg)
            mpz_sub(z, z, h->coeff[i]);
        mpz_set(f->coeff[i], z);
    }
    f->deg = maxdeg;
    mpz_clear(z);
    mpz_poly_cleandeg(f, maxdeg);
}

/* f[i] <- z */
void mpz_poly_setcoeff_double(mpz_poly f, int i, double z)
{
    mpz_poly_realloc (f, i + 1);
    mpz_set_d (f->coeff[i], z);
    if (i >= f->deg)
        mpz_poly_cleandeg(f, i);
}

/* return the number of real roots of the polynomial f[0]+f[1]*x+...+f[n]*x^d
   Assume f[d] is not zero.
*/
int numberOfRealRoots(mpz_poly f, usp_root_data* Roots)
{
    int i, nroots;
    mpz_t a, R, R1;
    double C, fd, x;
    mpz_poly r;
    mpz_poly_init(r, f->deg);

    mpz_init(a);
    //assert(mpz_cmp_ui (f[d], 0) != 0);
    nroots = 0; /* initialize number of roots found */
    if (mpz_cmp_ui(f->coeff[0], 0) == 0) /* root at 0 */
    {
        mpz_set_ui(a, 0);
        while (mpz_cmp_ui(f->coeff[0], 0) == 0)
        {
            divide(a, 0, f);
        }
    }
    double T = 0.0;
    fd = ln2(*mpz_poly_lc(f)); /* leading coefficient */
    C = (mpz_cmp_ui(f->coeff[f->deg - 1], 0) == 0) ? 0.0 : ln2(f->coeff[f->deg - 1]) - fd - log(1.0 * f->deg) / log(2.0);
    T = 0.0;
    for (i = 1; i <= f->deg; i++)
    {
        if (mpz_cmp_ui(f->coeff[f->deg - i], 0))
        {
            /* for i=1, we get ln2(f[d-1]) - fd, which is always
             larger than C = ln2(f[d-1]) - fd - log2(d) */
            x = (ln2(f->coeff[f->deg-i]) - fd) / (double) i;
            if (x > T) T = x;
        }
    }
    T += 1.0;
    T = T + log (1 + exp((C - T) / log(2.0))) / log(2.0);
  
    i = 1 + (int)T;
    mpz_set_ui(a, 1);
    mpz_mul_2exp(a, a, i);
    mpz_init(R);
    mpz_set(R, a);

    mpz_init(R1);
    mpz_neg(R1, R);
    i = usp(R1, R, 0, f->deg, signValue(R1, 0, f), signValue(R, 0, f),
           &nroots, f, r, Roots);

    mpz_clear(a);
    mpz_clear(R);
    mpz_clear(R1);
    mpz_poly_clear(r);

    return nroots;
}

/* divides in-place the input polynomial by 2^k*x-a */
void divide(mpz_t a, int k, mpz_poly f)
{
    int i;
    mpz_t u, v, w;

    mpz_init(u);
    mpz_init(v);
    mpz_init(w);
    mpz_set(u, f->coeff[f->deg]);
    for (i = f->deg - 1; i >= 0; i--)
    {
        mpz_tdiv_q_2exp(w, u, k); /* p[i] <- u/2^k */
        mpz_mul(v, a, w);
        mpz_add(u, f->coeff[i], v);
        mpz_set(f->coeff[i], w);
    }
    f->deg--; /* reduces degree by 1 */
    mpz_clear(u);
    mpz_clear(v);
    mpz_clear(w);
}

/* returns number of real roots (isolated) in a/2^m..b/2^m of polynomial
   f[0]+f[1]*x+...+f[d]*x^d
   r[0..d] is an auxiliary array.
   If R is not NULL, puts the roots in R.
*/
int usp(mpz_t a, mpz_t b, int m, int up, int va, int vb, int *nroots,
     mpz_poly f, mpz_poly r, usp_root_data *R)
{
    int lmi, i, k, c, s, smi, last, d;
    mpz_t mi, u, v, w;

    if (va * vb == 2 * (up % 2) - 1)
        up--;
    if (up == 0)
        return 0;
    else if (up == 1) {
        getroot(a, m, b, m, nroots, R);
        return 1;
    }
    mpz_init(mi);
    mpz_add(mi, a, b);
    lmi = m;
    if (mpz_fdiv_ui(mi, 2) == 0)
        mpz_tdiv_q_2exp(mi, mi, 1);
    else
        lmi++;
    smi = signValue(mi, lmi, f);
    if (smi == 0) { /* rational root at mi */
        int le, ri, i, d0 = f->deg;
        mpz_poly q;
        /* we cannot divide in-place, otherwise we will modify the input
        polynomial for the rest of the algorithm */
        mpz_poly_init(q, f->deg);
        mpz_poly_set(q, f);
        while (smi == 0) {
            getroot(mi, lmi, mi, lmi, nroots, R);
            divide(mi, lmi, q);
            smi = signValue(mi, lmi, q);
        }
        if (lmi > m) {
            mpz_mul_2exp(a, a, 1);
            mpz_mul_2exp(b, b, 1);
        }
        le = usp(a, mi, lmi, f->deg, signValue(a, lmi, q),
        signValue(mi, lmi, q), nroots, q, r, R);
        ri = usp(mi, b, lmi, f->deg, signValue(mi, lmi, q),
        signValue(b, lmi, q), nroots, q, r, R);
        mpz_clear(mi);
        mpz_poly_clear(q);
        return 1 + le + ri;
    }
    if (va * smi < 0) {
        if (up == 2) {
            getroot(a, m, mi, lmi, nroots, R);
            getroot(mi, lmi, b, m, nroots, R);
            mpz_clear(mi);
            return 2;
        }
        else if (vb * smi < 0) {
            mpz_t aa;

            mpz_init(aa);
            if (lmi > m)
            mpz_mul_2exp(aa, a, 1);
            else
            mpz_set(aa, a);
            c = usp(aa, mi, lmi, up, va, smi, nroots, f, r, R);
            if (c < up) {
                if (lmi > m)
                    mpz_mul_2exp(aa, b, 1);
                else
                    mpz_set(aa, b);
                c += usp(mi, aa, lmi, up - c, smi, vb, nroots, f, r, R);
            }
            mpz_clear(mi);
            mpz_clear(aa);
            return c;
        }
    }
    mpz_set(r->coeff[f->deg], f->coeff[f->deg]);
    mpz_init(w);
    for (i = f->deg-1; i >= 0; i--) {
        mpz_mul(r->coeff[i], b, r->coeff[i+1]); 
        mpz_mul_2exp(w, f->coeff[i], (f->deg-i) * m);
        mpz_add(r->coeff[i], r->coeff[i], w);
    }
    mpz_init(v);
    mpz_sub(v, a, b);
    mpz_init(u);
    mpz_set(u, v);
    for (k = 1; k < f->deg; k++) {
        for (i = f->deg-1; i >= k; i--) {
            mpz_mul(w, b, r->coeff[i+1]);
            mpz_add(r->coeff[i], r->coeff[i], w);
        }
        mpz_mul(r->coeff[k], r->coeff[k], u);
        mpz_mul(u, u, v);
    }
    mpz_clear(v);
    mpz_clear(w);
    mpz_mul(r->coeff[f->deg], r->coeff[f->deg], u);
    mpz_clear(u);
    last = sign(r->coeff[0]);
    d = f->deg-1;
    for (c = k = s = 0; k <= f->deg && c < 2; k++) {
        /* invariant: all signs in r[0]..r[f->deg-(d+1)] are equal */
        while (d > k && sign(r->coeff[f->deg - d]) * last >= 0)
        d--;
        if (d < k) {
            /* d+1 <= k, thus all signs in r[0]..r[f->deg-k] are equal,
            thus only one more sign change is possible */
            c += (sign(r->coeff[f->deg - k]) * s < 0);
            k = f->deg;
        }
        else {
            for (i = f->deg - 1; i >= k; i--)
                mpz_add(r->coeff[f->deg - i], r->coeff[f->deg - i], r->coeff[f->deg - i - 1]);
            i = mpz_cmp_ui(r->coeff[f->deg-k], 0);
            if (s * i < 0) {
                c++;
                if (va * vb > 0)
                c = 2;
            }
            if (i != 0)
                s = i; /* s is the last sign */
        }
        /* when k=f->deg-1 here and c=1, necessarily va * vb < 0, otherwise
        we would have c>=2 already, thus when we exit we cannot have
        c = 2 and k=f->deg+1 */
    }
    if (c == 1)
        getroot(a, m, b, m, nroots, R);
    else if (c > 1) {
        mpz_t aa;

        mpz_init(aa);
        //ASSERT(k <= f->deg);
        if (lmi > m)
            mpz_mul_2exp(aa, a, 1);
        else
            mpz_set(aa, a);
        c = usp(aa, mi, lmi, up, va, smi, nroots, f, r, R);
        if (c < up) {
            if (lmi > m)
                mpz_mul_2exp(aa, b, 1);
            else
                mpz_set(aa, b);
            c += usp(mi, aa, lmi, up-c, smi, vb, nroots, f, r, R);
        }
        mpz_clear(aa);
    }
    mpz_clear(mi);
    return c;
}

double ln2(mpz_t a)
{
    int l;
    double r;

    l = mpz_sizeinbase(a, 2);
    if (l <= 1024) /* a fits in a double */
        r = log(fabs (mpz_get_d (a))) / log (2.0);
    else {
        mpz_t b;

        mpz_init(b);
        mpz_tdiv_q_2exp(b, a, l - 900);
        r = mpz_get_d(b);
        r = (double) l - 900.0 + log(fabs (r)) / log(2.0);
        mpz_clear(b);
    }
    return r;
}

/* isolating interval is [a/2^ka, b/2^kb] */
void getroot(mpz_t a, int ka, mpz_t b, int kb, int *nroots, usp_root_data *R)
{
    if (R != NULL)
    {
        mpz_set(R[*nroots].a, a);
        R[*nroots].ka = ka;
        mpz_set(R[*nroots].b, b);
        R[*nroots].kb = kb;
    }
    *nroots += 1;
}

/* returns the sign of p(a/2^k) */
int signValue(mpz_t a, int k, mpz_poly f)
{
    int i, ret;
    mpz_t v, w;

    mpz_init(v);
    mpz_init(w);
    mpz_set(w, f->coeff[f->deg]);
    for (i = f->deg - 1; i >= 0; i--)
    {
        mpz_mul(w, w, a);
        mpz_mul_2exp(v, f->coeff[i], k * (f->deg - i));
        mpz_add(w, w, v);
    }
    ret = sign(w);
    mpz_clear(v);
    mpz_clear(w);
    return ret;
}

int sign (mpz_t a)
{
	int i;

	i = mpz_cmp_ui (a, 0);
	if (i > 0)
		return 1;
	else if (i < 0)
		return -1;
	else
	  	return 0;
}

void mpz_poly_discriminant(mpz_t disc, mpz_poly f)
{
    mpz_poly df;
    mpz_poly_init(df, f->deg);
    mpz_poly_derivative(df, f);
    mpz_poly_resultant(disc, f, df);
    //ASSERT(mpz_divisible_p(*disc, mpz_poly_lc(f)));
    mpz_divexact(disc, disc, *mpz_poly_lc(f));
    mpz_poly_clear(df);
}

/* Affects the derivative of f to df. Assumes df different from f.
   Assumes df has been initialized with degree at least f->deg-1. */
void mpz_poly_derivative(mpz_poly df, mpz_poly f)
{
  int i;

  df->deg = (f->deg <= 0) ? -1 : f->deg - 1;
  for (i = 0; i <= f->deg - 1; i++)
    mpz_mul_si (df->coeff[i], f->coeff[i + 1], i + 1);
}

/* Set res=f(x). Assumes res and x are different variables. */
void mpz_poly_eval(mpz_t res, mpz_poly f, mpz_t x) {
  int i, d;
  d = f->deg;
  if (d == -1) {
    mpz_set_ui(res, 0);
    return;
  }
  mpz_set(res, f->coeff[d]);
  for (i = d-1; i>=0; --i) {
    mpz_mul(res, res, x);
    mpz_add(res, res, f->coeff[i]);
  }
}

/* Set res=f(x) where x is an unsigned long. */
void mpz_poly_eval_ui (mpz_t res, mpz_poly f, unsigned long x)
{
  int d = f->deg;

  mpz_set (res, f->coeff[d]);
  for (int i = d - 1; i >= 0; i--)
  {
    mpz_mul_ui (res, res, x);
    mpz_add (res, res, f->coeff[i]);
  }
}

/* h=rem(h, f) mod N, f not necessarily monic, N not necessarily prime */
/* Coefficients of f must be reduced mod N
 * Coefficients of h need not be reduced mod N on input, but are reduced
 * on output.
 */
static int
mpz_poly_pseudodiv_r (mpz_poly h, mpz_poly f, mpz_t N, mpz_t factor)
{
  int i, d = f->deg, dh = h->deg;
  mpz_t tmp, inv;

  mpz_init_set_ui (inv, 1);
  mpz_init_set_ui (tmp, 1);

  mpz_set (tmp, f->coeff[d]);
  if (mpz_cmp_ui(tmp, 1) != 0) {
    /* inv is 1/f[d] mod N */
    if (!mpz_invert (inv, tmp, N)) {
      if (factor != NULL)
        mpz_gcd(factor, tmp, N);
      mpz_clear(inv);
      mpz_clear(tmp);
      return 0;
    }
  }

  while (dh >= d)
  {
    /* subtract h[dh]/f[d]*x^(dh-d)*f from h */
    if (mpz_cmp_ui(inv, 1) != 0) {
      mpz_mul (h->coeff[dh], h->coeff[dh], inv);
      mpz_ndiv_r (h->coeff[dh], h->coeff[dh], N);
    }

    for (i = 0; i < d; i++) {
      mpz_mul (tmp, h->coeff[dh], f->coeff[i]);
      mpz_mod (tmp, tmp, N);
      mpz_sub (h->coeff[dh - d + i], h->coeff[dh - d + i], tmp);
      mpz_ndiv_r (h->coeff[dh - d + i], h->coeff[dh - d + i], N);
    }

    do {
      dh --;
    }
    while (dh >= 0 && mpz_divisible_p (h->coeff[dh], N));

    h->deg = dh;
  }

  mpz_clear (inv);
  mpz_clear (tmp);
  return 1;
}

/* h=rem(h, f) mod p, f not necessarily monic. */
/* returns 0 if an inverse of the leading coeff could not be found. */
int
mpz_poly_div_r (mpz_poly h, mpz_poly f, mpz_t p)
{
  return mpz_poly_pseudodiv_r (h, f, p, NULL);
}

/*
   computes q, r such that f = q*g + r mod p, with deg(r) < deg(g)
   and p in mpz_t
   q and r must be allocated!
   the only aliasing allowed is f==r.
*/
int mpz_poly_div_qr (mpz_poly q, mpz_poly r, mpz_poly f, mpz_poly g, mpz_t p)
{
  int k, j, df = f->deg, dg = g->deg, dq = df - dg;
  mpz_t lg, invlg;

  if (df < dg) /* f is already reduced mod g */
    {
      mpz_poly_set_zero (q);
      mpz_poly_set (r, f);
      return 1;
    }

  /* now df >= dg */
  mpz_poly_realloc(q, dq + 1);

  mpz_init(lg);
  mpz_init_set_ui(invlg, 1);

  mpz_poly_set(r, f);
  q->deg = dq;

  mpz_set (lg, g->coeff[dg]);
  mpz_mod (lg, lg, p);
  /* invlg = 1/g[dg] mod p */
  if (mpz_cmp_ui(lg, 1) != 0)
    if (!mpz_invert(invlg, lg, p)) {
        mpz_clear(lg);
        mpz_clear(invlg);
        return 0;
    }


  for (k = df-dg ; k >=0 ; k--) {
    mpz_mul(q->coeff[k], r->coeff[k+dg], invlg);
    mpz_mod(q->coeff[k], q->coeff[k], p);
    for (j = dg+k ; j >= k ; j--) {
      mpz_submul(r->coeff[j], q->coeff[k], g->coeff[j-k]);
      mpz_mod(r->coeff[j], r->coeff[j], p);
    }
  }
  mpz_poly_cleandeg(r, r->deg);

  mpz_clear(invlg);
  mpz_clear(lg);
  return 1;
}

/* Set Q = P/lc(P) (mod m). Q and P might be identical. */
/* Coefficients of P need not be reduced mod m
 * Coefficients of Q are reduced mod m
 */
void mpz_poly_makemonic_mod_mpz (mpz_poly Q, mpz_poly P, mpz_t m)
{
  int i;
  mpz_t aux;
  mpz_init(aux);
  for(i = P->deg; i>=0; i--) {
      mpz_mod(aux, P->coeff[i], m);
      if (mpz_cmp_ui(aux, 0) != 0) break;
  }
  /* i is the degree of the leading monomial */
  Q->deg = i;
  if (i < 0) {
      /* if i == -1, then Q is the zero polynomial, there's nothing to do */
      mpz_clear(aux);
      return;
  }

  mpz_t aux2;
  mpz_init(aux2);
  mpz_invert(aux2, aux, m);
  for (i = 0; i < Q->deg; ++i) {
      mpz_mul(aux, aux2, P->coeff[i]);
      mpz_mod(aux, aux, m);
      mpz_poly_setcoeff(Q, i, aux);
  }
  mpz_clear(aux2);

  /* we can directly set the leading coefficient to 1 */
  mpz_poly_setcoeff_si (Q, Q->deg, 1);
  mpz_clear(aux);
}

/* Coefficients of A need not be reduced mod m
 * Coefficients of R are reduced mod m
 */
int mpz_poly_mod_mpz (mpz_poly R, mpz_poly A, mpz_t m)
{
  /* reduce lower coefficients */
  mpz_poly_realloc(R, A->deg + 1);
  for (int i = 0; i <= A->deg; ++i)
    mpz_mod (R->coeff[i], A->coeff[i], m);

  mpz_poly_cleandeg(R, A->deg);
  return R->deg;
}

/* f=gcd(f, g) mod p, with p in mpz_t */
/* clobbers g */
/* Coefficients of f and g need not be reduced mod p on input.
 * Coefficients of f are reduced mod p on output */
static void
mpz_poly_gcd_mpz_clobber (mpz_poly f, mpz_poly g, mpz_t p)
{
  /* First reduce mod p */
  mpz_poly_mod_mpz (f, f, p);
  mpz_poly_mod_mpz (g, g, p);
  while (g->deg >= 0)
    {
      mpz_poly_div_r (f, g, p);
      /* now deg(f) < deg(g): swap f and g */
      mpz_poly_swap (f, g);
    }
}

/* f <- gcd(a, b) mod p. */
/* Coefficients of a and b need not be reduced mod p
 * Coefficients of f are reduced mod p */
void
mpz_poly_gcd_mpz (mpz_poly f, mpz_poly a, mpz_poly b, mpz_t p)
{
    mpz_poly hh;
    if (f == b) {
        mpz_poly_init(hh, a->deg);
        mpz_poly_set(hh, a);
    } else {
        mpz_poly_init(hh, b->deg);
        mpz_poly_set(hh, b);
        mpz_poly_set(f, a); /* will do nothing if f = a */
    }
    mpz_poly_gcd_mpz_clobber(f, hh, p);
    mpz_poly_clear(hh);
}

/* factoring polynomials */

void mpz_poly_factor_list_init(mpz_poly_factor_list_ptr l)
{
    memset(l, 0, sizeof(*l));
}

void mpz_poly_factor_list_clear(mpz_poly_factor_list_ptr l)
{
    for(int i = 0 ; i < l->alloc ; i++)
        mpz_poly_clear(l->factors[i]->f);
    free(l->factors);
    memset(l, 0, sizeof(*l));
}

void mpz_poly_factor_list_flush(mpz_poly_factor_list_ptr l)
{
    /* There's a design choice here. We may elect to free everything.
     * Instead, we'll simply mark everything as zero, but keep all
     * allocated space.
     */
    for(int i = 0 ; i < l->alloc ; i++)
        l->factors[i]->f->deg = -1;
    l->size = 0;
}


void mpz_poly_factor_list_prepare_write(mpz_poly_factor_list_ptr l, int index)
{
    if (index >= l->alloc) {
        l->alloc = index + 4 + l->alloc / 4;
        l->factors = (mpz_poly_with_m*) realloc(l->factors, l->alloc * sizeof(mpz_poly_with_m));
        /* We need to set something. A zero polynomial has NULL storage
         * area, so that will do (realloc behaves as needed).  */
        for(int i = l->size ; i < l->alloc ; i++) {
            mpz_poly_init(l->factors[i]->f, -1);
            l->factors[i]->m = 0;
        }
    }
    if (l->size <= index)
        l->size = index + 1;
}

void mpz_poly_factor_list_push(mpz_poly_factor_list_ptr l, mpz_poly f, int m)
{
    mpz_poly_factor_list_prepare_write(l, l->size);
    mpz_poly_set(l->factors[l->size - 1]->f, f);
    l->factors[l->size - 1]->m = m;
}

/*
void mpz_poly_factor_list_fprintf(FILE* fp, mpz_poly_factor_list_srcptr l)
{
    for (int i = 0 ; i < l->size ; i++){
        char * res;
        int rc = mpz_poly_asprintf(&res,l->factors[i]->f);
        assert(rc >= 0);
        if (i) fprintf(fp, "*");
        fprintf(fp, "(%s)^%d", res, l->factors[i]->m);
        free(res);
    }
    fprintf(fp, "\n");
}
*/

/* Squarefree factorization */

/* This auxiliary function almost does the sqf. It fills
 * lf->factors[stride*i] (i from 1 to deg(f)) with the factors with
 * multiplicity i in f.  lf->factors[0] is filled with the product whose
 * multiplicity is a multiple of the field characteristic.  returns max
 * multiplicity stored (multiplied by the stride value).
 *
 * (stride has an importance for the recursive call, for p small. E.g. on
 * f=g^p, we get called with g=f^(1/p) and stride=p).
 *
 * We make no effort to check that lf is clean on input, which is to be
 * guaranteed by the caller (e.g. sufficiently many polynomials, all
 * equal to 1 -- or 0, which can easily be identified as something
 * unset).
 *
 * Coefficients of f need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
static int mpz_poly_factor_sqf_inner(mpz_poly_factor_list_ptr lf, mpz_poly f, int stride, mpz_t p)
{
    int r = 0;

    mpz_poly g, mi, mi1;
    mpz_poly t0,t1, T, tmp;
    mpz_poly_init(g, f->deg);
    mpz_poly_init(mi, f->deg);
    mpz_poly_init(mi1, f->deg);
    mpz_poly_init(t0, f->deg);
    mpz_poly_init(t1, f->deg);
    mpz_poly_init(T, f->deg);
    mpz_poly_init(tmp, f->deg);

    mpz_poly_derivative(t0, f);
    mpz_poly_gcd_mpz(g, f, t0, p);
    mpz_poly_divexact_modp(mi, f, g, p);
    /* mi is f/gcd(f,f') == all repeated prime factors of f whose
     * multiplicity isn't a multiple of the field characteristic.
     */

    mpz_poly_set_xi(T, 0);
    for(int i = 1 ; mi->deg > 0 ; i++) {
        /* note the weird argument ordering */
        mpz_poly_pow_ui_mod_f_mod_mpz(t0, mi, NULL, i, p);
        mpz_poly_divexact_modp(t1, f, t0, p);
        /* t1 = all polynomials in mi taken out from f with multiplicity i */
        mpz_poly_gcd_mpz(mi1, t1, mi, p);
        /* mi1 = almost like mi, but since factors with multiplicity i
         * are no longer in t1, there's absent from mi1 too. Whence
         * mi/mi1 is exactly the product of factors of multiplicity 1.
         */
        mpz_poly_factor_list_prepare_write(lf, i * stride);
        /* Use tmp so that we don't absurdly keep storage within
         * lf->factors */
        mpz_poly_divexact_modp(tmp, mi, mi1, p);
        mpz_poly_set(lf->factors[i * stride]->f, tmp);
        /* multiplicity field still unused at this point */
        mpz_poly_pow_ui_mod_f_mod_mpz(t0, lf->factors[i * stride]->f, NULL, i, p);
        mpz_poly_mul(T, T, t0);
        mpz_poly_mod_mpz(T, T, p);
        mpz_poly_swap(mi, mi1);
        r = i * stride;
    }

    mpz_poly_factor_list_prepare_write(lf, 0);
    mpz_poly_divexact_modp(lf->factors[0]->f, f, T, p);

    mpz_poly_clear(g);
    mpz_poly_clear(tmp);
    mpz_poly_clear(mi);
    mpz_poly_clear(mi1);
    mpz_poly_clear(t0);
    mpz_poly_clear(t1);
    mpz_poly_clear(T);
    return r;
}

/* Fills lf->factors[i] (i from 1 to deg(f)) with product of the factors
 * with multiplicity i in f, and to 1 if there are none. lf->factors[0]
 * is set to 1.  return the largest multiplicity stored. */
/* Coefficients of f0 need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
int mpz_poly_factor_sqf(mpz_poly_factor_list_ptr lf, mpz_poly f0,
			mpz_t p)
{
    /* factoring 0 doesn't make sense, really */
    assert(f0->deg >= 0);

    /* We'll call mpz_poly_factor_sqf_inner, possibly several times if
     * we are in small characteristic.
     */
    mpz_poly f;
    mpz_poly_init(f, f0->deg);
    mpz_poly_makemonic_mod_mpz(f, f0, p);
    assert(mpz_cmp_ui(*mpz_poly_lc(f), 1) == 0);

    int m = 0;
    int pu = mpz_get_ui(p);  // see below
    /* reset the factor list completely */
    mpz_poly_factor_list_flush(lf);

    for(int stride = 1 ; ; stride *= pu) {
        int r = mpz_poly_factor_sqf_inner(lf, f, stride, p);
        if (r > m) m = r;
        if (lf->factors[0]->f->deg == 0) {
            // if p is LAAAARGE, then of course we'll never have a linear
            // polynomial out of sqf_inner, thus we'll break early here.
            break;
        }
        /* divide coefficients */
        for(int i = 0 ; i <= lf->factors[0]->f->deg ; i++) {
            if (i % pu == 0) {
                mpz_set(f->coeff[i / pu], lf->factors[0]->f->coeff[i]);
            } else {
                assert (mpz_cmp_ui(lf->factors[0]->f->coeff[i], 0) == 0);
            }
        }
        f->deg = lf->factors[0]->f->deg / pu;
        mpz_poly_set_xi(lf->factors[0]->f, 0);
    }
    /* Now make sure that all factors in the factor list are non-zero */
    for(int i = 0 ; i < lf->size ; i++) {
        if (lf->factors[i]->f->deg < 0) {
            mpz_poly_set_xi(lf->factors[i]->f, 0);
        }
    }
    mpz_poly_clear(f);
    return m;
}




/* This performs distinct degree factorization */
/* Input polynomial must be squarefree -- otherwise repeated factors
 * probably won't show up in the factor list, or maybe at the wrong place
 * as parasites. */
/* Returns max degree of factors found (i.e. if the largest factors we
 * have are two factors of degree 7, we return 7, not 14). */
/* Coefficients of f0 need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
static int mpz_poly_factor_ddf_inner(mpz_poly_factor_list_ptr lf, mpz_poly f0, mpz_t p, int only_check_irreducible)
{
    mpz_poly g, gmx, x, tmp;
    mpz_poly f;
    mpz_poly_init(f, f0->deg);
    int i;

    /* factoring 0 doesn't make sense, really */
    assert(f0->deg >= 0);

    mpz_poly_makemonic_mod_mpz(f, f0, p);
    assert(mpz_cmp_ui(*mpz_poly_lc(f), 1) == 0);

    /* reset the factor list completely */
    mpz_poly_factor_list_flush(lf);

    mpz_poly_init (g, 2 * f->deg - 1);
    mpz_poly_init (gmx, 2 * f->deg - 1);
    mpz_poly_init (x, 1);
    mpz_poly_init (tmp, f->deg);

    mpz_poly_set_xi(x, 1);
    mpz_poly_set_xi(g, 1);

    for (i = 1; i <= f->deg ; ++i) {
        if (2 * i > f->deg) {
            /* Then we know that the remaining f is irreducible.  */
            mpz_poly_factor_list_prepare_write(lf, f->deg);
            for( ; i < f->deg ; i++) {
                /* multiplicity field still unused at this point */
                mpz_poly_set_xi(lf->factors[i]->f, 0);
            }
            /* multiplicity field still unused at this point */
            mpz_poly_swap(lf->factors[f->deg]->f, f);
            break;
        }

        /* g <- g^p mod fp */
        mpz_poly_pow_mod_f_mod_mpz (g, g, f, p, p);

        /* subtract x */
        mpz_poly_sub (gmx, g, x);
        mpz_poly_mod_mpz(gmx, gmx, p);

        /* lf[i] <- gcd (f, x^(p^i)-x) */
        mpz_poly_factor_list_prepare_write(lf, i);
        /* multiplicity field still unused at this point */

        /* see remark in _sqf regarding the relevance of tmp for storage */
        mpz_poly_gcd_mpz(tmp, f, gmx, p);
        mpz_poly_divexact_modp(f, f, tmp, p);
        mpz_poly_set(lf->factors[i]->f, tmp);

        /* Note for a mere irreducibility test: the length of the loop in
         * the irreducible case would still be deg(f)/2, and the penalty
         * caused by storing factors can be neglected.
         */
        if (only_check_irreducible && lf->factors[i]->f->deg > 0)
            break;

        if (f->deg == 0)
            break;
    }

    mpz_poly_clear (g);
    mpz_poly_clear (x);
    mpz_poly_clear (gmx);
    mpz_poly_clear (tmp);
    mpz_poly_clear (f);

    return i;
}

int mpz_poly_factor_ddf(mpz_poly_factor_list_ptr lf, mpz_poly f0,
			mpz_t p)
{
  return mpz_poly_factor_ddf_inner(lf, f0, p, 0);
}

/* Note that this also works for non squarefree polynomials -- the factor
 * list returned by mpz_poly_factor_ddf will be rubbish, but the m ==
 * f->deg test will tell the truth. */
int mpz_poly_is_irreducible(mpz_poly f, mpz_t p)
{
    mpz_poly_factor_list lf;
    mpz_poly_factor_list_init(lf);
    int m = mpz_poly_factor_ddf_inner(lf, f, p, 1);
    mpz_poly_factor_list_clear(lf);
    return m == f->deg;
}

/* Naive factorization of polynomials over GF(2). We assume we're always
 * dealing with polynomials of degree at most 10, so we're beter off
 * simply enumerating potential factors...
 */

/*
 * Add 1 to f. If the constant term is equal to 1, set this term to 0 and
 *  propagate the addition of 1 to the next coefficient, and so on.
 *
 * f: the polynomial on which the addition is computed, the modifications are
 *  made on f.
 */
static void mpz_poly_add_one_in_F2(mpz_poly f)
{
  assert(f->deg >= 1);

  int i = 0;
  while (mpz_cmp_ui(f->coeff[i], 1) == 0) {
    mpz_poly_setcoeff_si(f, i, 0);
    i++;
    if (i > f->deg) {
      break;
    }
  }
  mpz_poly_setcoeff_si(f, i, 1);
}

/*
 * Factorize naively a mpz_poly mod 2.
 *
 * Return the number of factors found.
 */
static int mpz_poly_factor2(mpz_poly_factor_list_ptr list, mpz_poly f,
    mpz_t p)
{
  assert(mpz_cmp_ui(p, 2) == 0);

  //make a copy of f.
  mpz_poly fcopy;
  mpz_poly_init(fcopy, f->deg);
  mpz_poly_set(fcopy, f);

  //reduce all the coefficient mod p.
  mpz_t coeff;
  mpz_init(coeff);
  for (int i = 0; i <= f->deg; i++) {
    mpz_poly_getcoeff(coeff, i, f);
    mpz_mod(coeff, coeff, p);
    mpz_poly_setcoeff(fcopy, i, coeff);
  }

  //Purge list.
  mpz_poly_factor_list_flush(list);

  //If deg(f) in F2 is less than 1, we have the factor.
  if (fcopy->deg < 1) {
    mpz_clear(coeff);
    mpz_poly_clear(fcopy);
    assert(list->size == 0);
    return list->size;
  } else if (fcopy->deg == 1) {
    mpz_clear(coeff);
    mpz_poly_factor_list_push(list, fcopy, 1);
    mpz_poly_clear(fcopy);
    assert(list->size == 1);
    return list->size;
  }

  if (mpz_poly_is_irreducible(f, p)) {
    //If f is irreducible mod 2, fcopy is the factorisation of f mod 2.
    mpz_poly_factor_list_push(list, fcopy, 1);
  } else {
    //Create the first possible factor.
    mpz_poly tmp;
    mpz_poly_init(tmp, 1);
    mpz_poly_setcoeff_si(tmp, 1, 1);

    //Enumerate all the possible factor of f mod 2.
    while (tmp->deg <= fcopy->deg) {
      //tmp is a possible factor.
      if (mpz_poly_is_irreducible(tmp, p)) {
        mpz_poly q;
        mpz_poly_init(q, 0);
        mpz_poly r;
        mpz_poly_init(r, 0);
        //Euclidean division of fcopy
        mpz_poly_div_qr(q, r, fcopy, tmp, p);
        //Power of the possible factor.
        unsigned int m = 0;
        //While fcopy is divisible by tmp.
        while (r->deg == -1) {
          //Increase the power of tmp.
          m++;
          mpz_poly_set(fcopy, q);
          if (fcopy->deg == 0 || fcopy->deg == -1) {
            //No other possible factor.
            break;
          }
          mpz_poly_div_qr(q, r, fcopy, tmp, p);
        }
        if (m != 0) {
          //Push tmp^m as a factor of f mod 2.
          mpz_poly_factor_list_push(list, tmp, m);
        }
        mpz_poly_clear(q);
        mpz_poly_clear(r);
      }
      //Go to the next possible polynomial in F2.
      mpz_poly_add_one_in_F2(tmp);
    }
    mpz_poly_clear(tmp);
  }

#ifndef NDBEBUG
  //Verify if the factorisation is good.
  mpz_poly_cleandeg(fcopy, -1);
  for (int i = 0; i <= f->deg; i++) {
    mpz_poly_getcoeff(coeff, i, f);
    mpz_mod(coeff, coeff, p);
    mpz_poly_setcoeff(fcopy, i, coeff);
  }

  mpz_poly fmul;
  mpz_poly_init(fmul, -1);
  mpz_poly_set(fmul, list->factors[0]->f);
  for (int j = 1; j < list->factors[0]->m; j++) {
    mpz_poly_mul(fmul, fmul, list->factors[0]->f);
  }
  for (int i = 1; i < list->size ; i++) {
    for (int j = 0; j < list->factors[i]->m; j++) {
      mpz_poly_mul(fmul, fmul, list->factors[i]->f);
    }
  }
  for (int i = 0; i <= fmul->deg; i++) {
    mpz_poly_getcoeff(coeff, i, fmul);
    mpz_mod(coeff, coeff, p);
    mpz_poly_setcoeff(fmul, i, coeff);
  }

  assert(mpz_poly_cmp(fcopy, fmul) == 0);

  mpz_poly_clear(fmul);
#endif // NDBEBUG

  mpz_clear(coeff);
  mpz_poly_clear(fcopy);

  return list->size;
}

/*
 * this tries to split between squares and non-squares -- it's the most
 * basic split of course, and the other splits merely randomize on top of
 * this. Note that this building block must be changed for characteristic
 * two
 *
 * returns non-zero if a non-trivial or split is obtained.
 *
 * k is the degree of factors we are looking for.
 *
 * We split in two parts:
 *
 *  - factors whose roots are squares in GF(p^k).
 *  - factors whose roots are non squares in GF(p^k).
 *
 * Note that we do not find the factor X this way; this is to be done by
 * the caller.
 *
 * Obviously, we include some shift, and hope that eventually there is a
 * shift that works. Large characteristic is generally happy with some
 * translation shift. Small characteristic may need more general shifts.
 *
 * Coefficients of f0 need not be reduced mod p.
 * Coefficients of g[0] and g[1] are reduced mod p.
 */
static void mpz_poly_factor_edf_pre(mpz_poly g[2], mpz_poly f, int k,
				    mpz_t p)
{
    int nontrivial = 0;
    mpz_poly_set_xi(g[0], 0);
    mpz_poly_set_xi(g[1], 0);

    assert (f->deg > k);

    mpz_poly xplusa;
    mpz_poly_init(xplusa, 1);

    mpz_t half_pk;
    mpz_init(half_pk);
    mpz_pow_ui(half_pk, p, k);
    mpz_fdiv_q_ui(half_pk, half_pk, 2); /* (p^k-1)/2 */

    for(unsigned long a = 0 ; ! nontrivial ; a++) {
        /* we want to iterate on monic polynomials of degree <= k-1. */
        /* In order to bear in mind what happens in large enough
         * characteristic, we'll name these polynomials xplusa, although
         * it does not have to be x+a (and it can't be restricted to only
         * x+a if the characteristic is small -- that does not give
         * enough legroom).
         */
        if (mpz_fits_ulong_p(p)) {
            unsigned long pz = mpz_get_ui(p);
            if (a == 0) {
                /* special case, really */
                mpz_poly_set_xi(xplusa, 1);
            } else {
                /* write a in base p, and add 1 */
                int i = 0;
                for(unsigned long tmp = a ; tmp ; i++, tmp /= pz) {
                    mpz_poly_setcoeff_ui(xplusa, i, tmp % pz);
                }
                mpz_poly_setcoeff_ui(xplusa, i, 1);
            }
        } else {
            /* take the polynomial x+a */
            mpz_poly_set_xi(xplusa, 1);
            mpz_poly_setcoeff_ui(xplusa, 0, a);
        }

        mpz_poly_pow_mod_f_mod_mpz(g[0], xplusa, f, half_pk, p);

        mpz_poly_add_ui(g[1], g[0], 1);     /* (x+a)^((p^k-1)/2) + 1 */
        mpz_poly_sub_ui(g[0], g[0], 1);     /* (x+a)^((p^k-1)/2) - 1 */

        mpz_poly_mod_mpz(g[0], g[0], p);
        mpz_poly_mod_mpz(g[1], g[1], p);

        mpz_poly_gcd_mpz(g[0], g[0], f, p);
        mpz_poly_gcd_mpz(g[1], g[1], f, p);

        if (g[0]->deg + g[1]->deg < f->deg) {
            /* oh, we're lucky. x+a is a factor ! */
            int s = g[0]->deg > g[1]->deg;
            /* multiply g[s] by (x+a) */
            mpz_poly_mul_mod_f_mod_mpz(g[s], g[s], xplusa, f, p, NULL);
        }
        assert(g[0]->deg + g[1]->deg == f->deg);
        assert(g[0]->deg % k == 0);
        assert(g[1]->deg % k == 0);

        nontrivial += g[0]->deg != 0 && g[0]->deg != f->deg;
        nontrivial += g[1]->deg != 0 && g[1]->deg != f->deg;
    }
    mpz_clear(half_pk);
    mpz_poly_clear(xplusa);
}

/* This factors f, and for each factor q found, store q in lf.
 * Return the number of distinct factors found. */
/* Coefficients of f need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
static int mpz_poly_factor_edf_inner(mpz_poly_factor_list_ptr lf, mpz_poly f, int k, mpz_t p, gmp_randstate_t rstate)
{
    if (f->deg == k) {
        mpz_poly_factor_list_push(lf, f, 1);
        return 1;
    }
    if (f->deg == 0) {
        return 0;
    }

    mpz_poly h[2];

    mpz_poly_init(h[0], f->deg);
    mpz_poly_init(h[1], f->deg);

    mpz_poly_factor_edf_pre(h, f, k, p);

    int n = 0;

    n += mpz_poly_factor_edf_inner(lf, h[0], k, p, rstate);
    mpz_poly_clear(h[0]);

    n += mpz_poly_factor_edf_inner(lf, h[1], k, p, rstate);
    mpz_poly_clear(h[1]);

    return n;
}

// returns f0->deg / d
/* Coefficients of f need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
int mpz_poly_factor_edf(mpz_poly_factor_list_ptr lf, mpz_poly f, int k, mpz_t p, gmp_randstate_t rstate)
{
    mpz_poly_factor_list_flush(lf);

    if (mpz_cmp_ui(p, 2) == 0) {
        /* we need some other code for edf. Currently we have very naive
         * code, but that's good enough for small degree. */
        return mpz_poly_factor2(lf, f, p);
    }

    int v = mpz_poly_valuation(f);
    if (v) {
        /* Since our input is square-free, then we expect v==1.
         * Furthermore, k prescribes the extension field where the
         * expected roots live, thus for 0 it should really be 1. */
        assert(v == 1 && k == 1);
        mpz_poly_factor_list_prepare_write(lf, lf->size);
        mpz_poly_set_xi(lf->factors[lf->size-1]->f, 1);

        mpz_poly f1;
        mpz_poly_init(f1, f->deg - 1);
        mpz_poly_div_xi(f1, f, v);
        int n = 1 + mpz_poly_factor_edf_inner(lf, f1, k, p, rstate);
        mpz_poly_clear(f1);
        return n;
    }

    int ret = mpz_poly_factor_edf_inner(lf, f, k, p, rstate);
    return ret;
}

typedef int (*sortfunc_t) (const void *, const void *);

static int mpz_poly_with_m_cmp(mpz_poly_with_m * a, mpz_poly_with_m * b)
{
    int r = mpz_poly_cmp((*a)->f, (*b)->f);
    if (r) return r;
    return ((*a)->m > (*b)->m) - ((*b)->m > (*a)->m);
}


/* putting it all together */
/* Coefficients of f need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
int mpz_poly_factor(mpz_poly_factor_list lf, mpz_poly f, mpz_t p, gmp_randstate_t rstate)
{
    mpz_poly_factor_list sqfs, ddfs, edfs;
    mpz_poly_factor_list_init(sqfs);
    mpz_poly_factor_list_init(ddfs);
    mpz_poly_factor_list_init(edfs);

    mpz_poly_factor_list_flush(lf);

    int maxmult = mpz_poly_factor_sqf (sqfs, f, p);
    for(int m = 1 ; m <= maxmult ; m++) {
        assert(sqfs->factors[m]->f->deg >= 0);
        if (sqfs->factors[m]->f->deg == 0) continue;
        int maxdeg = mpz_poly_factor_ddf(ddfs, sqfs->factors[m]->f, p);
        for(int k = 1 ; k <= maxdeg ; k++) {
            assert(ddfs->factors[k]->f->deg >= 0);
            if(ddfs->factors[k]->f->deg == 0) continue;
            mpz_poly_factor_edf(edfs, ddfs->factors[k]->f, k, p, rstate);
            for(int j = 0 ; j < edfs->size ; j++) {
                /* register this factor, with multiplicity m */
                /* cheat a bit... */
                mpz_poly_factor_list_prepare_write(lf, lf->size);
                mpz_poly_swap(lf->factors[lf->size-1]->f, edfs->factors[j]->f);
                mpz_poly_makemonic_mod_mpz(lf->factors[lf->size-1]->f, lf->factors[lf->size-1]->f, p);
                lf->factors[lf->size-1]->m = m;
            }
        }
    }
    mpz_poly_factor_list_clear(edfs);
    mpz_poly_factor_list_clear(ddfs);
    mpz_poly_factor_list_clear(sqfs);

    /* sort factors by degree and lexicographically */
    qsort(lf->factors, lf->size, sizeof(mpz_poly_with_m), (sortfunc_t) &mpz_poly_with_m_cmp);

    return lf->size;
}

/* a <- b cmod c with -c/2 <= a < c/2 */
void mpz_ndiv_r (mpz_t a, mpz_t b, mpz_t c)
{
  mpz_mod (a, b, c); /* now 0 <= a < c */

  size_t n = mpz_size (c);
  mp_limb_t aj, cj;
  int sub = 0, sh = GMP_NUMB_BITS - 1;

  if (mpz_getlimbn (a, n-1) >= (mp_limb_t) 1 << sh)
    sub = 1;
  else
    {
      while (n-- > 0)
	{
	  cj = mpz_getlimbn (c, n);
	  aj = mpz_getlimbn (a, n) << 1;
	  if (n > 0)
	    aj |= mpz_getlimbn (a, n-1) >> sh;
	  if (aj > cj)
	    {
	      sub = 1;
	      break;
	    }
	  else if (aj < cj)
	    break;
	}
    }
  if (sub)
    mpz_sub (a, a, c);
}

/* Coefficients of f must be reduced mod m on input
 * Coefficients of P need not be reduced mod p.
 * Coefficients of Q are reduced mod p.
 */
void mpz_poly_pow_ui_mod_f_mod_mpz (mpz_poly Q, mpz_poly P, mpz_poly f,
                          unsigned long a, mpz_t p)
{
    mpz_t az;
    mpz_init_set_ui(az, a);
    mpz_poly_pow_mod_f_mod_mpz(Q, P, f, az, p);
    mpz_clear(az);
}

/* Q = P^a mod f, mod p. Note, p is mpz_t */
/* f may be NULL, in case there is only reduction mod p */
/* Coefficients of f must be reduced mod p on input
 * Coefficients of P need not be reduced mod p.
 * Coefficients of Q are reduced mod p
 */
void mpz_poly_pow_mod_f_mod_mpz (mpz_poly Q, mpz_poly P,
			    mpz_poly f, mpz_t a, mpz_t p)
{
  int k = mpz_sizeinbase(a, 2), l, L = 0, j;
  mpz_poly R, *T = NULL;
  mpz_t invf;

  if (mpz_cmp_ui(a, 0) == 0) {
    mpz_poly_set_xi (Q, 0);
    return;
  }

  mpz_poly_init(R, f ? 2*f->deg : -1);

  // Initialize R to P
  mpz_poly_set(R, P);

  /* compute invf = 1/p mod lc(f) */
  if (f != NULL)
    {
      mpz_init (invf);
      mpz_invert (invf, p, f->coeff[f->deg]);
    }

  /* We use base-2^l exponentiation with sliding window,
     thus we need to precompute P^2, P^3, ..., P^(2^l-1).
     The expected average cost k squarings plus k/(l+1) + 2^(l-1) multiplies.
  */

  for (l = 1; k / (l + 1) + (1 << (l - 1)) > k / (l + 2) + (1 << l); l++);
  /* this gives (roughly) l=1 for k < 8, l=2 for k < 27, l=3 for k < 84,
     l=4 for k < 245 */
  if (l > 1)
    {
      L = 1 << (l - 1);
      T = (mpz_poly*) malloc (L * sizeof (mpz_poly));
      /* we store P^2 in T[0], P^3 in T[1], ..., P^(2^l-1) in T[L-1] */
      for (j = 0; j < L; j++)
	mpz_poly_init (T[j], f ? 2*f->deg : -1);
      mpz_poly_sqr_mod_f_mod_mpz (T[0], R, f, p, invf);             /* P^2 */
      mpz_poly_mul_mod_f_mod_mpz (T[1], T[0], R, f, p, invf);       /* P^3 */
      for (j = 2; j < L; j++)
	mpz_poly_mul_mod_f_mod_mpz (T[j], T[j-1], T[0], f, p, invf);
    }

  // Horner
  for (k -= 2; k >= 0;)
  {
    while (k >= 0 && mpz_tstbit (a, k) == 0)
      {
	mpz_poly_sqr_mod_f_mod_mpz (R, R, f, p, invf);
	k --;
      }
    if (k < 0)
      break;
    j = mpz_scan1 (a, (k >= l) ? k - (l - 1) : 0);
    /* new window starts at bit k, and ends at bit j <= k */
    int e = 0;
    while (k >= j)
      {
	mpz_poly_sqr_mod_f_mod_mpz (R, R, f, p, invf);
	e = 2 * e + mpz_tstbit (a, k);
	k --;
      }
    mpz_poly_mul_mod_f_mod_mpz (R, R, (e == 1) ? P : T[e/2], f, p, invf);
  }

  mpz_poly_swap(Q, R);
  mpz_poly_clear(R);
  if (f != NULL)
    mpz_clear (invf);
  if (l > 1)
    {
      for (k = 0; k < L; k++)
	mpz_poly_clear (T[k]);
      free (T);
    }
}


/* Q = P1*P2 mod f, mod m
   f is the original algebraic polynomial (non monic but small coefficients)
   Coefficients of P1 and P2 need not be reduced mod m.
   Coefficients of Q are reduced mod m. */
void mpz_poly_mul_mod_f_mod_ui (mpz_poly Q, mpz_poly P1, mpz_poly P2,
			    mpz_poly f, int64_t m0)
{
	mpz_t m; mpz_init_set_ui(m, m0);
	
	mpz_poly_mul_mod_f_mod_mpz(Q, P1, P2, f, m, NULL);

	mpz_clear(m);
}

/* Q = P1*P2 mod f, mod m
   f is the original algebraic polynomial (non monic but small coefficients)
   Coefficients of P1 and P2 need not be reduced mod m.
   Coefficients of Q are reduced mod m.
   If invf is not NULL, it is 1/m mod lc(f). */
void mpz_poly_mul_mod_f_mod_mpz (mpz_poly Q, mpz_poly P1, mpz_poly P2,
			    mpz_poly f, mpz_t m, mpz_t invf)
{
  int d1 = P1->deg;
  int d2 = P2->deg;
  int d = d1+d2;
  mpz_poly R;

  mpz_poly_init(R, d);
#ifdef MPZ_POLY_TIMINGS
  if (mpz_fits_sint_p (f->coeff[0])){
    START_TIMER;
    //d = mpz_poly_mul_tc (R->coeff, P1->coeff, d1, P2->coeff, d2);
	mpz_poly_mul(R, P1, P2);
	d = R->deg;
    mpz_poly_cleandeg(R, d);
    END_TIMER (TIMER_MUL);
    // reduce mod f
    RESTART_TIMER;
    mpz_poly_mod_f_mod_mpz (R, f, m, invf);
    END_TIMER (TIMER_RED);
  }else
#endif
  {
    //d = mpz_poly_mul_tc (R->coeff, P1->coeff, d1, P2->coeff, d2);
	mpz_poly_mul(R, P1, P2);
	d = R->deg;
    mpz_poly_cleandeg(R, d);
    // reduce mod f
    mpz_poly_mod_f_mod_mpz (R, f, m, invf);
  }
  mpz_poly_set(Q, R);
  mpz_poly_clear(R);
}

/* Q = P1*P2 mod f, assuming f is monic */
void mpz_poly_mul_mod_f (mpz_poly Q, mpz_poly P1, mpz_poly P2,
                        mpz_poly f)
{
    mpz_poly_mul(Q,P1,P2);
    mpz_poly_div_r_z(Q,Q,f);
}

/* Q = P^2 mod f, mod m
   f is the original algebraic polynomial (non monic but small coefficients)
   Coefficients of P need not be reduced mod m.
   Coefficients of Q are reduced mod m.
   If not NULL, invf = 1/m mod lc(f). */
void mpz_poly_sqr_mod_f_mod_mpz (mpz_poly Q, mpz_poly P, mpz_poly f,
			    mpz_t m, mpz_t invf)
{
  int d1 = P->deg;
  int d = d1 + d1;
  mpz_poly R;

  mpz_poly_init(R, d);

  /* Fast squaring in 2d1+1 squares, i.e., 2d-1 squares.
     For d=5, this gives 9 squares. */
  // compute timing only if f has short coefficients
#ifdef MPZ_POLY_TIMINGS
  if (mpz_fits_sint_p (f->coeff[0])){
    START_TIMER;
    d = mpz_poly_sqr (R->coeff, P->coeff, d1);
    mpz_poly_cleandeg(R, d);
    END_TIMER (TIMER_SQR);
    // reduce mod f
    RESTART_TIMER;
    mpz_poly_mod_f_mod_mpz (R, f, m, invf);
    END_TIMER (TIMER_RED);
  }else
#endif
  {
    d = mpz_poly_sqr (R->coeff, P->coeff, d1);
    mpz_poly_cleandeg(R, d);
    // reduce mod f
    mpz_poly_mod_f_mod_mpz (R, f, m, invf);
  }

  mpz_poly_set(Q, R);
  mpz_poly_clear(R);
}

/* Reduce R[d]*x^d + ... + R[0] mod f[df]*x^df + ... + f[0] modulo m.
   Return the degree of the remainder.
   Coefficients of f must be reduced mod m on input.
   Coefficients of R need not be reduced mod m on input, but are reduced
   on output.
   If invf is not NULL, it should be 1/m mod lc(f). */
int mpz_poly_mod_f_mod_mpz (mpz_poly R, mpz_poly f, mpz_t m,
			mpz_t invf)
{
  mpz_t aux, c;
  size_t size_m, size_f;

  if (f == NULL)
    goto reduce_R;

  if (invf == NULL)
    {
      mpz_init (aux);
      /* aux = 1/m mod lc(f) */
      mpz_invert (aux, m, f->coeff[f->deg]);
    }

  size_m = mpz_size (m);
  size_f = mpz_poly_size (f);

  mpz_init (c);
  // FIXME: write a subquadratic variant
  while (R->deg >= f->deg)
    {
      /* Here m is large (thousand to million bits) and lc(f) is small
       * (typically one word). We first subtract lambda * m * x^(dR-df) ---
       * which is zero mod m --- to R such that the new coefficient of degree
       * dR is divisible by lc(f), i.e., lambda = lc(R)/m mod lc(f). Then if
       * c = (lc(R) - lambda * m) / lc(f), we subtract c * x^(dR-df) * f. */
      mpz_mod (c, R->coeff[R->deg], f->coeff[f->deg]); /* lc(R) mod lc(f) */
      mpz_mul (c, c, (invf == NULL) ? aux : invf);
      mpz_mod (c, c, f->coeff[f->deg]);    /* lc(R)/m mod lc(f) */
      mpz_submul (R->coeff[R->deg], m, c);  /* lc(R) - m * (lc(R) / m mod lc(f)) */
      assert (mpz_divisible_p (R->coeff[R->deg], f->coeff[f->deg]));
      mpz_divexact (c, R->coeff[R->deg], f->coeff[f->deg]);
      /* If R[deg] has initially size 2n, and f[deg] = O(1), then c has size
	 2n here. However, in the equal-degree factorization, even if f[deg]
	 = O(1), the lower coefficients of f might have n bits. Thus we decide
	 to reduce whenever the total size exceeds 2n. */
      if (mpz_size (c) + size_f > 2 * size_m)
	mpz_mod (c, c, m);
      for (int i = R->deg - 1; i >= R->deg - f->deg; --i)
	mpz_submul (R->coeff[i], c, f->coeff[f->deg - R->deg + i]);
      R->deg--;
    }

  mpz_clear (c);
  if (invf == NULL)
    mpz_clear (aux);

 reduce_R:
  mpz_poly_mod_mpz (R, R, m);

  return R->deg;
}

/* f <- g^2 where g has degree r */
static int mpz_poly_sqr (mpz_t *f, mpz_t *g, int r)
{
  int i, j;
  assert(f != g);
  for (i = 0; i <= 2 * r; i++)
    mpz_set_ui (f[i], 0);
  for (i = 0; i <= r; ++i)
    for (j = 0; j < i; ++j)
      mpz_addmul(f[i+j], g[i], g[j]);
  for (i = 0; i < 2*r ; i++)
      mpz_mul_2exp(f[i], f[i], 1);
  for (i = 0; i < r ; i++) {
      mpz_mul(f[2*r], g[i], g[i]);
      mpz_add(f[2*i], f[2*i], f[2*r]);
  }
  mpz_mul(f[2*r], g[r], g[r]);
  return 2 * r;
}

/* This also computes q and r such that f = q * g + r, but over Z, not
 * modulo a prime. Also, we do not assume that g is monic. Of course, if
 * it is not, then most often the result will be undefined (over Z). We
 * return something well-defined if q and r happen to be integer
 * polynomials.
 * We return 0 if this is not the case (in which case q and r are
 * undefined).
 * r==f is allowed.
 */
int mpz_poly_div_qr_z (mpz_poly q, mpz_poly r, mpz_poly f, mpz_poly g)
{
  int k, j, df = f->deg, dg = g->deg, dq = df - dg;

  if (df < dg) /* f is already reduced mod g */
    {
      mpz_poly_set_zero (q);
      mpz_poly_set (r, f);
      return 1;
    }

  /* now df >= dg */
  mpz_poly_realloc(q, dq + 1);

  mpz_poly_set(r, f);
  q->deg = dq;

  mpz_srcptr lg = g->coeff[dg];

  for (k = df-dg ; k >=0 ; k--) {
    if (!mpz_divisible_p(r->coeff[k+dg], lg))
        return 0;
    mpz_divexact(q->coeff[k], r->coeff[k+dg], lg);
    for (j = dg+k ; j >= k ; j--) {
      mpz_submul(r->coeff[j], q->coeff[k], g->coeff[j-k]);
    }
  }
  mpz_poly_cleandeg(r, r->deg);
  return 1;
}
int mpz_poly_div_r_z (mpz_poly r, mpz_poly f, mpz_poly g)
{
    mpz_poly quo;
    mpz_poly_init(quo,-1);
    int ret = mpz_poly_div_qr_z(quo,r,f,g);
    mpz_poly_clear(quo);
    return ret;
}

/* Set g=f+a where a is unsigned long. */
void mpz_poly_add_ui (mpz_poly g, mpz_poly f, unsigned long a)
{
    mpz_poly_set(g, f);
    mpz_poly_realloc(g, 1);
    if (f->deg >= 0) {
        mpz_add_ui(g->coeff[0], f->coeff[0], a);
        g->deg = f->deg;
        if (g->deg == 0)
            mpz_poly_cleandeg(g, g->deg);
    } else {
        mpz_set_ui(g->coeff[0], a);
        g->deg = 0;
    }
}

/* Set g=f-a where a is unsigned long. */
void mpz_poly_sub_ui (mpz_poly g, mpz_poly f, unsigned long a)
{
    mpz_poly_set(g, f);
    mpz_poly_realloc(g, 1);
    if (f->deg >= 0) {
        mpz_sub_ui(g->coeff[0], f->coeff[0], a);
        g->deg = f->deg;
        if (g->deg == 0)
            mpz_poly_cleandeg(g, g->deg);
    } else {
        mpz_set_ui(g->coeff[0], a);
        mpz_neg(g->coeff[0], g->coeff[0]);
        g->deg = 0;
    }
}

/* return 0 if f and g are equal,
 * -1 if f is "smaller" and 1 if f is "bigger", for some arbitrary
 * ordering (lowest degree first, then lex order).
 *
 * Assumes f and g are normalized */
int mpz_poly_cmp (mpz_poly a, mpz_poly b)
{
    int r = (a->deg > b->deg) - (b->deg > a->deg);
    if (r) return r;
    for(int d = a->deg; d >= 0 ; d--) {
        r = mpz_cmp(a->coeff[d], b->coeff[d]);
        if (r) return r;
    }
    return 0;
}

/* return the maximal limb-size of the coefficients of f */
size_t mpz_poly_size (mpz_poly f)
{
  size_t S = 0, s;
  int i;
  int d = f->deg;

  for (i = 0; i <= d; i++)
  {
    s = mpz_size (f->coeff[i]);
    if (s > S)
      S = s;
  }
  return S;
}

/* g <- quo (f, x^i) */
void mpz_poly_div_xi(mpz_poly g, mpz_poly f, int i)
{
    if (f->deg < i) {
        mpz_poly_set_zero(g);
        return;
    }
    if (i == 0) {
        mpz_poly_set(g, f);
        return;
    }
    if (g == f) {
        /* rotate the coefficients, don't do any freeing */
        mpz_t * temp = (mpz_t*) malloc(i * sizeof(mpz_t));
        memcpy(temp, g->coeff, i * sizeof(mpz_t));
        memmove(g->coeff, g->coeff + i, (g->deg + 1 - i) * sizeof(mpz_t));
        memcpy(g->coeff + (g->deg + 1 - i), temp, i * sizeof(mpz_t));
        g->deg -= i;
        free(temp);
        return;
    }

    mpz_poly_realloc(g, 1 + f->deg - i);
    for(int j = i ; j <= f->deg ; j++) {
        mpz_set(g->coeff[j - i], f->coeff[j]);
    }
    g->deg = f->deg - i;
}

/* g <- f * x^i */
void mpz_poly_mul_xi (mpz_poly g, mpz_poly f, int i)
{
    if (i == 0) {
        mpz_poly_set(g, f);
        return;
    }

    mpz_poly_realloc(g, 1 + f->deg + i);

    if (g == f) {
        for(int j = g->deg ; j >= 0 ; j--) {
            mpz_swap(g->coeff[j + i], g->coeff[j]);
        }
    } else {
        for(int j = g->deg ; j >= 0 ; j--) {
            mpz_set(g->coeff[j + i], g->coeff[j]);
        }
    }
    for(int j = 0 ; j < i ; j++) {
        mpz_set_ui(g->coeff[j], 0);
    }
    g->deg = f->deg + i;
}

/* return the valuation of f */
int mpz_poly_valuation(mpz_poly f)
{
    int n = 0;
    assert(f->deg >= 0);
    for( ; n < f->deg  && mpz_cmp_ui(f->coeff[n], 0) == 0 ; n++) ;
    return n;
}

/* computes d = gcd(f, g) = u*f + v*g mod p, with p in mpz_t */
/* Coefficients of f and g need not be reduced mod p.
 * Coefficients of d, u, v are reduced mod p */
void mpz_poly_xgcd_mpz (mpz_poly d, mpz_poly f, mpz_poly g, mpz_poly u, mpz_poly v, mpz_t p)
{
  mpz_poly q, tmp;
  mpz_poly gg;

  if (f->deg < g->deg) {
      mpz_poly_xgcd_mpz(d, g, f, v, u, p);
      return;
  }
  mpz_poly_init(gg, g->alloc);
  mpz_poly_set(d, f);
  mpz_poly_set(gg, g);
  mpz_poly_mod_mpz(d, d, p);
  mpz_poly_mod_mpz(gg, gg, p);

  mpz_poly uu, vv;
  mpz_poly_init (uu, 0);
  mpz_poly_init (vv, 0);

  mpz_poly_set_xi(u, 0);
  mpz_poly_set_zero(uu);

  mpz_poly_set_xi(vv, 0);
  mpz_poly_set_zero(v);

  mpz_poly_init(q, d->deg);
  mpz_poly_init(tmp, d->deg + gg->deg);

  while (gg->deg >= 0)
    {

      /* q, r := f div g mod p */
      mpz_poly_div_qr (q, d, d, gg, p);

      /* u := u - q * uu mod p */
      mpz_poly_mul(tmp, q, uu);
      mpz_poly_sub_mod_mpz(u, u, tmp, p);
      mpz_poly_swap (u, uu);

      /* v := v - q * vv mod p */
      mpz_poly_mul(tmp, q, vv);
      mpz_poly_sub_mod_mpz(v, v, tmp, p);
      mpz_poly_swap (v, vv);

      /* now deg(f) < deg(g): swap f and g */
      mpz_poly_swap (d, gg);
    }

  /* make monic */
  mpz_t inv;
  if (mpz_cmp_ui(d->coeff[d->deg], 1) != 0)
    {
      mpz_init(inv);
      mpz_invert(inv, d->coeff[0], p);
      mpz_poly_mul_mpz(d, d, inv);
      mpz_poly_mod_mpz(d, d, p);
      mpz_poly_mul_mpz(u, u, inv);
      mpz_poly_mod_mpz(u, u, p);
      mpz_poly_mul_mpz(v, v, inv);
      mpz_poly_mod_mpz(v, v, p);
      mpz_clear(inv);
    }

  mpz_poly_clear(gg);
  mpz_poly_clear(uu);
  mpz_poly_clear(vv);
  mpz_poly_clear(q);
  mpz_poly_clear(tmp);
}

/* Set f=g-h (mod m). */
void mpz_poly_sub_mod_mpz(mpz_poly f, mpz_poly g, mpz_poly h, mpz_t m)
{
    mpz_poly_sub(f, g, h);
    mpz_poly_mod_mpz(f, f, m);
}

