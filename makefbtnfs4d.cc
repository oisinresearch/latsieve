#include <stdint.h>	// int64_t
#include <iostream> // cout
#include <iomanip> // setprecision
#include <gmpxx.h>
#include "intpoly.h"
#include <math.h>	// sqrt
#include <fstream>	// file
#include <ctime>	// clock_t
#include <cstring>	// memset
#include <omp.h>
#include "mpz_poly.h"
#include "mpz_poly_bivariate.h"
#include "mpz_poly_Fq.h"
#include <sstream>	// stringstream

using std::cout;
using std::endl;
using std::flush;
//using std::vector;
using std::string;
using std::ifstream;
using std::fixed;
using std::scientific;
using std::setprecision;
using std::sort;
using std::to_string;
using std::hex;
using std::stringstream;


int main(int argc, char** argv)
{
	if (argc == 1) {
		cout << "usage: ./makefbtnfs polyfile.poly fbb factorbasefile" << endl;
		return 0;
	}

	bool verbose = true;
		
	if (verbose) cout << endl << "Reading input polynomial in file " << argv[1] << "..." << flush;
	//vector<mpz_class> fpoly;
	//vector<mpz_class> gpoly;
	mpz_t* fpoly = new mpz_t[20];	// max degree of 20.  Not the neatest
	mpz_t* gpoly = new mpz_t[20];	// max degree of 20.  Not the neatest
	mpz_t* hpoly = new mpz_t[20];	// max degree of 20.  Not the neatest
	for (int i = 0; i < 20; i++) {
		mpz_init(fpoly[i]);
		mpz_init(gpoly[i]);
		mpz_init(hpoly[i]);
	}
	string line;
	char linebuffer[100];
	ifstream file(argv[1]);
	getline(file, line);	// first line contains number n to factor
	// read nonlinear poly
	int degfhy = -1;
	if (verbose) cout << endl << "Side 0 polynomial fh_y (ascending coefficients)" << endl;
	while (getline(file, line) && line.substr(0,3) == "fhy" ) {
		line = line.substr(line.find_first_of(" ")+1);
		//mpz_set_str(c, line.c_str(), 10);
		mpz_set_str(fpoly[++degfhy], line.c_str(), 10);
		//mpz_get_str(linebuffer, 10, fpoly[degfhy-1]);
		if (verbose) cout << line << endl << flush;
	}
	//int degfhy = fpoly.size();
	// read other poly
	int degghy = -1;
	bool read = true;
	if (verbose) cout << endl << "Side 1 polynomial gh_y: (ascending coefficients)" << endl;
	while (read && line.substr(0,3) == "ghy" ) {
		line = line.substr(line.find_first_of(" ")+1);
		//mpz_set_str(c, line.c_str(), 10);
		mpz_set_str(gpoly[++degghy], line.c_str(), 10);
		//mpz_get_str(linebuffer, 10, gpoly[degghy-1]);
		if (verbose) cout << line << endl << flush;
		read = static_cast<bool>(getline(file, line));
	}
	//int degghy = gpoly.size();
	// read other poly
	int degh = -1;
	read = true;
	if (verbose) cout << endl << "Tower polynomial h: (ascending coefficients)" << endl;
	while (read && line.substr(0,1) == "h" ) {
		line = line.substr(line.find_first_of(" ")+1);
		//mpz_set_str(c, line.c_str(), 10);
		mpz_set_str(hpoly[++degh], line.c_str(), 10);
		//mpz_get_str(linebuffer, 10, gpoly[degg-1]);
		if (verbose) cout << line << endl << flush;
		read = static_cast<bool>(getline(file, line));
	}
	// read bivariate F0 poly
	mpz_poly_bivariate F0;
	mpz_poly_bivariate_init(F0, 0);	// init to deg 0 (constant)
	mpz_poly F0i;
	mpz_poly_init(F0i, 0); // init to deg 0 (constant)
	mpz_t F0ij; mpz_init(F0ij);
	read = true;
	if (verbose) cout << endl << "Bivariate polynomial F0: (ascending coefficients)" << endl;
	int inow = 0;
	while (read && line.substr(0,1) == "f" ) {
		int u = line.find_first_of("_");
		string ch = line.substr(1, u-1);
		int inew = atoi(ch.c_str());
		ch = line.substr(u+1, line.find_first_of(":")-u-1);
		int  j = atoi(ch.c_str());
		if (inew == inow) {
			line = line.substr(line.find_first_of(" ")+1);
			mpz_set_str(F0ij, line.c_str(), 10);
			mpz_poly_setcoeff(F0i, j, F0ij);
		}
		else {
			mpz_poly_bivariate_setcoeff(F0, inow, F0i);
			inow = inew;
			F0i->deg = 0;
			line = line.substr(line.find_first_of(" ")+1);
			mpz_set_str(F0ij, line.c_str(), 10);
			mpz_poly_setcoeff(F0i, j, F0ij);
		}
		if (verbose) cout << line << endl << flush;
		read = static_cast<bool>(getline(file, line));
	}
	mpz_poly_bivariate_setcoeff(F0, inow, F0i);
	// read bivariate F1 poly
	mpz_poly_bivariate F1;
	mpz_poly_bivariate_init(F1, 0);	// init to deg 0 (constant)
	mpz_poly F1i;
	mpz_poly_init(F1i, 0); // init to deg 0 (constant)
	mpz_t F1ij; mpz_init(F1ij);
	read = true;
	if (verbose) cout << endl << "Bivariate polynomial F1: (ascending coefficients)" << endl;
	inow = 0;
	while (read && line.substr(0,1) == "g" ) {
		int u = line.find_first_of("_");
		string ch = line.substr(1, u-1);
		int inew = atoi(ch.c_str());
		ch = line.substr(u+1, line.find_first_of(":")-u-1);
		int  j = atoi(ch.c_str());
		if (inew == inow) {
			line = line.substr(line.find_first_of(" ")+1);
			mpz_set_str(F1ij, line.c_str(), 10);
			mpz_poly_setcoeff(F1i, j, F1ij);
		}
		else {
			mpz_poly_bivariate_setcoeff(F1, inow, F1i);
			inow = inew;
			F1i->deg = 0;
			line = line.substr(line.find_first_of(" ")+1);
			mpz_set_str(F1ij, line.c_str(), 10);
			mpz_poly_setcoeff(F1i, j, F1ij);
		}
		if (verbose) cout << line << endl << flush;
		read = static_cast<bool>(getline(file, line));
	}
	mpz_poly_bivariate_setcoeff(F1, inow, F1i);
	file.close();
	//mpz_clear(c);
	mpz_poly h0; mpz_poly_init(h0, degh);
	mpz_poly f0; mpz_poly f1;
	mpz_poly_init(f0, degfhy); mpz_poly_init(f1, degghy);
	mpz_poly_set_mpz(h0, hpoly, degh);
	mpz_poly_set_mpz(f0, fpoly, degfhy);
	mpz_poly_set_mpz(f1, gpoly, degghy);
	if (verbose) cout << endl << "Complete.  Degree fh_y = " << degfhy << ", degree gh_y = " << degghy << "." << endl;

	if (verbose) cout << endl << "Starting sieve of Eratosthenes for small primes..." << endl << flush;
	int64_t fbb = 1<<21;
	if (argc >=3) fbb = atoi(argv[2]);
	int64_t max = fbb; // 10000000;// 65536;
	char* sieve = new char[max+1]();
	int64_t* primes = new int64_t[max]; // int64_t[155611]; //new int64_t[809228];	//new int64_t[6542]; 	// 2039 is the 309th prime, largest below 2048
	for (int64_t i = 2; i <= sqrt(max); i++)
		if(!sieve[i])
			for (int64_t j = i*i; j <= max; j += i)
				if(!sieve[j]) sieve[j] = 1;
	int64_t nump = 0;
	for (int64_t i = 2; i <= max-1; i++)
		if (!sieve[i])
			primes[nump++] = i;
	if (verbose) cout << "Complete." << endl;

	// set up constants
	std::clock_t start; double timetaken = 0;
	int K = 0;
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		if (id == 0) K = omp_get_num_threads();
	}
 	
	gmp_randstate_t rstate;
	gmp_randinit_mt(rstate);
    gmp_randseed_ui(rstate, 123);

	vector<int64_t> q0sieve0;
	vector<int64_t> q1sieve0;
	vector<int64_t> q2sieve0;
	vector<int64_t> q4sieve0;
	vector<int64_t> q0sieve1;
	vector<int64_t> q1sieve1;
	vector<int64_t> q2sieve1;
	vector<int64_t> q4sieve1;

	int64_t itenpc0 = nump / 10;
	int64_t itotal = 0;
	// compute factor base
	if (verbose) cout << endl << "Constructing factor base with " << K << " threads." << endl << flush;
	//if (verbose) cout << endl << "[0%]   constructing factor base..." << endl << flush;
	start = clock();
	#pragma omp parallel
	{
		int id = omp_get_thread_num();

		mpz_poly_factor_list lf;
		mpz_poly_factor_list_init(lf);
		mpz_poly_factor_list lh;
		mpz_poly_factor_list_init(lh);
		mpz_poly_factor_list l3;
		mpz_poly_factor_list_init(l3);
		mpz_poly hlin; mpz_poly_init(hlin, 0);
		mpz_poly_bivariate Hlin; mpz_poly_bivariate_init(Hlin, 0);
		mpz_poly f0hlin; mpz_poly_init(f0hlin, 0);
		mpz_poly_bivariate* factors0 = new mpz_poly_bivariate[f0->deg / 2];
		mpz_poly_bivariate* factors1 = new mpz_poly_bivariate[f1->deg / 2];
		for (int i = 0; i < f0->deg / 2; i++) mpz_poly_bivariate_init(factors0[i], 0);
		for (int i = 0; i < f1->deg / 2; i++) mpz_poly_bivariate_init(factors1[i], 0);
		// compute discriminants of f0 and f1
		mpz_t Df0; mpz_init(Df0); mpz_t Df1; mpz_init(Df1);
		mpz_poly_discriminant(Df0, f0);
		mpz_poly_discriminant(Df1, f1);

	#pragma omp for
		for (int64_t i = 0; i < nump; i++) {
			int64_t q0 = primes[i];
			// skip q0 if it is ramified in Q[x]/<f0>
			if (mpz_mod_ui(NULL, Df0, q0) == 0)
				continue;
			mpz_t q; mpz_init_set_ui(q, q0);
			// determine which type of special-q we have, first factor f mod q
			mpz_poly_factor(lf, f0, q, rstate);
			bool allf4 = true;
			for (int j = 0; j < lf->size; j++)
				if (lf->factors[j]->f->deg != 4) {
					allf4 = false;
					break;
				}
			bool allf2 = true;
			if (allf4) allf2 = false;
			else
				for (int j = 0; j < lf->size; j++)
					if (lf->factors[j]->f->deg != 2) {
						allf2 = false;
						break;
					}
			bool allf1 = true;
			if (allf4 || allf2) allf1 = false;
			else
				for (int j = 0; j < lf->size; j++)
					if (lf->factors[j]->f->deg != 1) {
						allf1 = false;
						break;
					}
			mpz_poly_factor(lh, h0, q, rstate);
			bool allh1 = true;
			for (int j = 0; j < lh->size; j++)
				if (lh->factors[j]->f->deg != 1) {
					allh1 = false;
					break;
				}

			// categorize special-q by allf4, allf2, allh1
			if (allf4) // q inert, 0 roots
				q0sieve0.push_back(q0);
			else if ((allf2 && allh1) || (!allf4 && !allf2 && !allf1)) { // 1 root 
				if (allf2 && allh1) {
					for (int j = 0; j < lh->size; j++) {
						int64_t m = mpz_get_ui(lh->factors[j]->f->coeff[0]);
						q1sieve0.push_back(q0);
						q1sieve0.push_back(m);
					}
				}
				else {	// some degree 1 ideals, some degree 2 ideal(s)
					for (int j = 0; j < lh->size; j++) {
						int64_t m = mpz_get_ui(lh->factors[j]->f->coeff[0]);
						mpz_poly_setcoeff_ui(hlin, 0, m);
						mpz_poly_setcoeff_ui(hlin, 1, 1);
						mpz_poly_bivariate_setcoeff(Hlin, 0, hlin);
						mpz_poly_bivariate_resultant_y(f0hlin, F0, Hlin);
						// now factor f0hlin mod q
						mpz_poly_factor(l3, f0hlin, q, rstate);
						if (l3->factors[0]->f->deg == 2) {	// smallest degree is 2
							q1sieve0.push_back(q0);
							q1sieve0.push_back(m);
						}
						mpz_poly_factor_list_flush(l3);
					}
				}
			}
			else if (allf1 || (!allf4 && !allf2 && !allf1)) { // 1 root of f and 1 root of h
				for (int j = 0; j < lf->size; j++) {
					if (lf->factors[j]->f->deg == 1) {
						int64_t R = mpz_get_ui(lf->factors[j]->f->coeff[0]);
						// find corresponding root r of h
						for (int k = 0; k < lh->size; k++) {
							int64_t r = mpz_get_ui(lh->factors[k]->f->coeff[0]);
							mpz_poly_setcoeff_ui(hlin, 0, r);
							mpz_poly_setcoeff_ui(hlin, 1, 1);
							mpz_poly_bivariate_setcoeff(Hlin, 0, hlin);
							mpz_poly_bivariate_resultant_y(f0hlin, F0, Hlin);
							// now factor f0hlin mod q
							mpz_poly_factor(l3, f0hlin, q, rstate);
							bool foundR = false;
							for (int l = 0; l < l3->size; l++) {
								int64_t Rtest = mpz_get_ui(l3->factors[l]->f->coeff[0]);
								if (R == Rtest) {
									foundR = true;
									// now we have R, r
									q2sieve0.push_back(q0);
									q2sieve0.push_back(r);
									q2sieve0.push_back(R);
									break;
								}
							}
							mpz_poly_factor_list_flush(l3);
							if (foundR) break;
						}
					}
				}
			}
			else { // mpz_poly_Fq_factor, 4 values
				mpz_poly_Fq_factor_edf(2, F0, q0, h0, factors0);
				for (int j = 0; j < f0->deg / 2; j++) {
					int64_t a0 = mpz_get_ui(factors0[0]->coeff[0]->coeff[0]);
					int64_t a1 = mpz_get_ui(factors0[0]->coeff[0]->coeff[1]);
					// rel1 = x + a1*y + a0
					// rel2 = rel1*y + rel1
					int64_t h0_0 = mpz_get_ui(h0->coeff[0]);
					int64_t h0_1 = mpz_get_ui(h0->coeff[1]);
					int64_t b0 = a0 + (-h0_0*a1);
					int64_t b1 = a0 + (-h0_1 + 1);
					// Note:  The above only works for monic, quadratic h0
					q4sieve0.push_back(q0);
					q4sieve0.push_back(a0);
					q4sieve0.push_back(a1);
					q4sieve0.push_back(b0);
					q4sieve0.push_back(b1);
				}
			}

			mpz_poly_factor_list_flush(lf);
			mpz_poly_factor_list_flush(lh);
	#pragma omp atomic
			itotal++;
			if (itotal % itenpc0 == 0) {
	#pragma omp critical
				if (verbose) cout << "[" << 100 * itotal / nump + 1 << "%]\tConstructing factor base..." << endl;
			}				 
		}
	}
	timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC / K;
	start = clock();
	timetaken += ( clock() - start ) / (double) CLOCKS_PER_SEC;
	int k0 = q0sieve0.size();
	int k1 = q1sieve0.size();
	int k2 = q2sieve0.size();
	int k4 = q4sieve0.size();
	int l0 = q0sieve1.size();
	int l1 = q1sieve1.size();
	int l2 = q2sieve1.size();
	int l4 = q4sieve1.size();
	if (verbose) cout << "Complete.  Time taken: " << timetaken << "s" << endl << flush;
	if (verbose) cout << "There are " << k0 + (k1/2) + (k2/3) + (k4/5) << " factor base primes on side 0." << endl << flush;
	if (verbose) cout << "There are " << l0 + (l1/2) + (l2/3) + (l4/5) << " factor base primes on side 1." << endl << flush;
	//if (verbose) cout << "There are " << kh << " factor base primes in tower." << endl << flush;

	// write output file
	FILE* out;
	out = fopen(argv[3], "w+");
	fprintf(out, "%d\n", fbb);
	// side 0
	fprintf(out, "%d\n", k0);
	fprintf(out, "%d\n", k1);
	fprintf(out, "%d\n", k2);
	fprintf(out, "%d\n", k4);
	for (int i = 0; i < k0; i++) {
		fprintf(out, "%d\n", q0sieve0[i]);
	}
	for (int i = 0; i < k1/2; i++) {
		fprintf(out, "%d", q1sieve0[i*2]);
        fprintf(out, ",%d", q1sieve0[i*2+1]);
		fprintf(out, "\n");
	}
	for (int i = 0; i < k2/3; i++) {
		fprintf(out, "%d", q2sieve0[i*3]);
		for (int j = 1; j < 3; j++)
			fprintf(out, ",%d", q2sieve0[i*3+j]);
		fprintf(out, "\n");
	}
	for (int i = 0; i < k4/5; i++) {
		fprintf(out, "%d", q4sieve0[i*5]);
		for (int j = 1; j < 5; j++)
			fprintf(out, ",%d", q4sieve0[i*5+j]);
		fprintf(out, "\n");
	}
	// side 1
	fprintf(out, "%d\n", l0);
	fprintf(out, "%d\n", l1);
	fprintf(out, "%d\n", l2);
	fprintf(out, "%d\n", l4);
	for (int i = 0; i < l0; i++) {
		fprintf(out, "%d\n", q0sieve1[i]);
	}
	for (int i = 0; i < l1/2; i++) {
		fprintf(out, "%d", q1sieve1[i*2]);
        fprintf(out, ",%d", q1sieve1[i*2+1]);
		fprintf(out, "\n");
	}
	for (int i = 0; i < l2/3; i++) {
		fprintf(out, "%d", q2sieve1[i*3]);
		for (int j = 1; j < 3; j++)
			fprintf(out, ",%d", q2sieve1[i*3+j]);
		fprintf(out, "\n");
	}
	for (int i = 0; i < l4/5; i++) {
		fprintf(out, "%d", q4sieve1[i*5]);
		for (int j = 1; j < 5; j++)
			fprintf(out, ",%d", q4sieve1[i*5+j]);
		fprintf(out, "\n");
	}
	fclose(out);

	// free memory
	gmp_randclear(rstate);
	delete[] primes;
	delete[] sieve;
	mpz_clear(F1ij);
	mpz_poly_clear(F1i);
	mpz_poly_bivariate_clear(F1);
	mpz_clear(F0ij);
	mpz_poly_clear(F0i);
	mpz_poly_bivariate_clear(F0);
	for (int i = 0; i < 20; i++) {
		mpz_clear(hpoly[i]);
		mpz_clear(gpoly[i]);
		mpz_clear(fpoly[i]);
	}
	delete[] hpoly;
	delete[] gpoly;
	delete[] fpoly;

	return 0;
}


