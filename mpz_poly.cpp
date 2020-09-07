#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mpz_poly.h"
#include <algorithm>
#include <iostream> // cout

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
    //ASSERT_ALWAYS(mpz_cmp_ui (f[d], 0) != 0);
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

