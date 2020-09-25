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
	int degfht = -1;
	if (verbose) cout << endl << "Side 0 polynomial fh_t (ascending coefficients)" << endl;
	while (getline(file, line) && line.substr(0,3) == "fht" ) {
		line = line.substr(line.find_first_of(" ")+1);
		//mpz_set_str(c, line.c_str(), 10);
		mpz_set_str(fpoly[++degfht], line.c_str(), 10);
		//mpz_get_str(linebuffer, 10, fpoly[degfht-1]);
		if (verbose) cout << line << endl << flush;
	}
	//int degfht = fpoly.size();
	// read other poly
	int degght = -1;
	bool read = true;
	if (verbose) cout << endl << "Side 1 polynomial gh_t: (ascending coefficients)" << endl;
	while (read && line.substr(0,3) == "ght" ) {
		line = line.substr(line.find_first_of(" ")+1);
		//mpz_set_str(c, line.c_str(), 10);
		mpz_set_str(gpoly[++degght], line.c_str(), 10);
		//mpz_get_str(linebuffer, 10, gpoly[degght-1]);
		if (verbose) cout << line << endl << flush;
		read = static_cast<bool>(getline(file, line));
	}
	//int degght = gpoly.size();
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
	if (verbose) cout << endl << "Complete.  Degree fh_t = " << degfht << ", degree gh_t = " << degght << "." << endl;

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
	mpz_t r0; mpz_init(r0);
	int64_t* S0 = new int64_t[degfht * nump]();
	int64_t* sieveS0 = new int64_t[degfht * nump]();
	int* num_S0modp = new int[nump]();
	int64_t* sieveP0 = new int64_t[nump]();
	int* sievenum_S0modp = new int[nump]();
	int64_t* S1 = new int64_t[degght * nump]();
	int64_t* sieveS1 = new int64_t[degght * nump]();
	int* num_S1modp = new int[nump]();
	int64_t* sieveP1 = new int64_t[nump]();
	int* sievenum_S1modp = new int[nump]();
	int64_t* s = new int64_t[degh * nump]();
	int* num_smodp = new int[nump]();
	int64_t* sieves0 = new int64_t[degh * nump]();
	int64_t* sievep0 = new int64_t[nump]();
	int64_t* sieves1 = new int64_t[degh * nump]();
	int64_t* sievep1 = new int64_t[nump]();
	int* sievenum_s0modp = new int[nump]();
	int* sievenum_s1modp = new int[nump]();
	int64_t* sj_f = new int64_t[degfht * nump]();
	int64_t* sj_g = new int64_t[degght * nump]();
	int64_t itenpc0 = nump / 10;
	int64_t itotal = 0;
	mpz_poly Fh_x; mpz_poly_init(Fh_x, 0);
	mpz_poly_bivariate Ap; mpz_poly_bivariate_init(Ap, 0);
	mpz_poly Ap0; mpz_poly_init(Ap0, 0);
	mpz_t res; mpz_init(res);
	// compute factor base
	if (verbose) cout << endl << "Constructing factor base with " << K << " threads." << endl << flush;
	//if (verbose) cout << endl << "[0%]   constructing factor base..." << endl << flush;
	start = clock();
	#pragma omp parallel
	{
		mpz_t rt; mpz_init(rt); 
		int id = omp_get_thread_num();
		int64_t* S0temp = new int64_t[degfht];
		int64_t* fp = new int64_t[degfht+1]();
		int64_t* S1temp = new int64_t[degght];
		int64_t* gp = new int64_t[degght+1]();
		int64_t* stemp = new int64_t[degh];
		int64_t* hp = new int64_t[degh+1]();

	#pragma omp for
		for (int64_t i = 0; i < nump; i++) {
			int64_t p = primes[i];
			for (int j = 0; j <= degfht; j++) fp[j] = mpz_mod_ui(rt, fpoly[j], p);
			int degfp = degfht; while (fp[degfp] == 0 || degfp == 0) degfp--;
			int numS0 = polrootsmod(fp, degfp, S0temp, p);
			num_S0modp[i] = numS0;
			for (int j = 0; j < numS0; j++) S0[i*degfht + j] = S0temp[j];
			for (int j = 0; j <= degght; j++) gp[j] = mpz_mod_ui(rt, gpoly[j], p);
			int deggp = degght; while (gp[deggp] == 0 || deggp == 0) deggp--;
			int numS1 = polrootsmod(gp, deggp, S1temp, p);
			num_S1modp[i] = numS1;
			for (int j = 0; j < numS1; j++) S1[i*degght + j] = S1temp[j];
			for (int j = 0; j <= degh; j++) hp[j] = mpz_mod_ui(rt, hpoly[j], p);
			int deghp = degh; while (hp[deghp] == 0 || deghp == 0) deghp--;
			int nums = polrootsmod(hp, deghp, stemp, p);
			num_smodp[i] = nums;
			for (int j = 0; j < nums; j++) s[i*degh + j] = stemp[j];
			// compute roots of h corresponding to roots of f for this p
			for (int j = 0; j < numS0; j++) {
				mpz_poly_setcoeff_si(Ap0, 0, -S0temp[j]);
				mpz_poly_bivariate_setcoeff(Ap, 0, Ap0);
				mpz_poly_setcoeff_ui(Ap0, 0, 1);	// Aq = x - SK
				mpz_poly_bivariate_setcoeff(Ap, 1, Ap0);
				mpz_poly_bivariate_resultant_y(Fh_x, F0, Ap);
				for (int k = 0; k < nums; k++) {
					mpz_poly_eval_ui(res, Fh_x, stemp[k]);
					if (mpz_mod_ui(r0, res, p) == 0) {  // is s[k] valid?
						sj_f[i*degfht + j] = stemp[k];
						break;
					}
				}
			}
			// compute roots of h corresponding to roots of g for this p
			for (int j = 0; j < numS1; j++) {
				mpz_poly_setcoeff_si(Ap0, 0, -S1temp[j]);
				mpz_poly_bivariate_setcoeff(Ap, 0, Ap0);
				mpz_poly_setcoeff_ui(Ap0, 0, 1);	// Aq = x - SK
				mpz_poly_bivariate_setcoeff(Ap, 1, Ap0);
				mpz_poly_bivariate_resultant_y(Fh_x, F1, Ap);
				for (int k = 0; k < nums; k++) {
					mpz_poly_eval_ui(res, Fh_x, stemp[k]);
					if (mpz_mod_ui(r0, res, p) == 0) {  // is s[k] valid?
						sj_g[i*degfht + j] = stemp[k];
						break;
					}
				}
			}
	#pragma omp atomic
			itotal++;
			if (itotal % itenpc0 == 0) {
	#pragma omp critical
				if (verbose) cout << "[" << 100 * itotal / nump + 1 << "%]\tConstructing factor base..." << endl;
			}				 
		}

		delete[] hp;
		delete[] stemp;
		delete[] gp;
		delete[] S1temp;
		delete[] fp;
		delete[] S0temp;
		mpz_clear(rt);
	}
	timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC / K;
	start = clock();
	int k0 = 0; int k1 = 0;
	for (int i = 0; i < nump; i++) {
		int nums = num_smodp[i];
		int numS0 = num_S0modp[i];
		if (numS0 > 0 && nums > 0) {
			sieveP0[k0] = primes[i];
			for (int j = 0; j < numS0; j++) sieveS0[k0*degfht + j] = S0[i*degfht + j];
			sievenum_S0modp[k0] = numS0;
			sievep0[k0] = primes[i];
			for (int j = 0; j < nums; j++) sieves0[k0*degh + j] = s[i*degh + j];
			sievenum_s0modp[k0++] = nums;
		}
		int numS1 = num_S1modp[i];
		if (numS1 > 0 && nums > 0) {
			sieveP1[k1] = primes[i];
			for (int j = 0; j < numS1; j++) sieveS1[k1*degght + j] = S1[i*degght + j];
			sievenum_S1modp[k1] = numS1;
			sievep1[k1] = primes[i];
			for (int j = 0; j < nums; j++) sieves1[k1*degh + j] = s[i*degh + j];
			sievenum_s1modp[k1++] = nums;
		}
	}
	timetaken += ( clock() - start ) / (double) CLOCKS_PER_SEC;
	if (verbose) cout << "Complete.  Time taken: " << timetaken << "s" << endl << flush;
	if (verbose) cout << "There are " << k0 << " factor base primes on side 0." << endl << flush;
	if (verbose) cout << "There are " << k1 << " factor base primes on side 1." << endl << flush;
	//if (verbose) cout << "There are " << kh << " factor base primes in tower." << endl << flush;

	// write output file
	FILE* out;
	out = fopen(argv[3], "w+");
	fprintf(out, "%d\n", fbb);
	fprintf(out, "%d\n", k0);
	for (int i = 0; i < k0; i++) {
		fprintf(out, "%d", sieveP0[i]);
		for (int j = 0; j < sievenum_S0modp[i]; j++)
            fprintf(out, ",%d", sieveS0[i*degfht + j]);
		fprintf(out, "\n");
	}
	fprintf(out, "%d\n", k0);
	for (int i = 0; i < k0; i++) {
		fprintf(out, "%d", sievep0[i]);
		for (int j = 0; j < sievenum_S0modp[i]; j++)
            fprintf(out, ",%d", sj_f[i*degfht + j]);
		fprintf(out, "\n");
	}
	fprintf(out, "%d\n", k1);
	for (int i = 0; i < k1; i++) {
		fprintf(out, "%d", sieveP1[i]);
		for (int j = 0; j < sievenum_S1modp[i]; j++)
            fprintf(out, ",%d", sieveS1[i*degght + j]);
		fprintf(out, "\n");
	}
	fprintf(out, "%d\n", k1);
	for (int i = 0; i < k1; i++) {
		fprintf(out, "%d", sievep1[i]);
		for (int j = 0; j < sievenum_S1modp[i]; j++)
            fprintf(out, ",%d", sj_g[i*degh + j]);
		fprintf(out, "\n");
	}
	fclose(out);

	// free memory
	mpz_clear(res);
	mpz_poly_clear(Ap0);
	mpz_poly_bivariate_clear(Ap);
	mpz_poly_clear(Fh_x);
	delete[] sj_g;
	delete[] sj_f;
	delete[] sievenum_s1modp;
	delete[] sievenum_s0modp;
	delete[] sievep1;
	delete[] sievep0;
	delete[] num_smodp;
	delete[] s;
	delete[] sievenum_S1modp;
	delete[] sieveP1;
	delete[] num_S1modp;
	delete[] S1;
	delete[] sievenum_S0modp;
	delete[] sieveP0;
	delete[] num_S0modp;
	delete[] S0;
	mpz_clear(r0);
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


