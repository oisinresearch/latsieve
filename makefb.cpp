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
		cout << "usage: ./makefb polyfile.poly fbb factorbasefile" << endl;
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
	int degf = -1;
	if (verbose) cout << endl << "Side 0 polynomial f0 (ascending coefficients)" << endl;
	while (getline(file, line) && line.substr(0,1) == "c" ) {
		line = line.substr(line.find_first_of(" ")+1);
		//mpz_set_str(c, line.c_str(), 10);
		mpz_set_str(fpoly[++degf], line.c_str(), 10);
		//mpz_get_str(linebuffer, 10, fpoly[degf-1]);
		if (verbose) cout << line << endl << flush;
	}
	//int degf = fpoly.size();
	// read other poly
	int degg = -1;
	bool read = true;
	if (verbose) cout << endl << "Side 1 polynomial f1: (ascending coefficients)" << endl;
	while (read && line.substr(0,1) == "Y" ) {
		line = line.substr(line.find_first_of(" ")+1);
		//mpz_set_str(c, line.c_str(), 10);
		mpz_set_str(gpoly[++degg], line.c_str(), 10);
		//mpz_get_str(linebuffer, 10, gpoly[degg-1]);
		if (verbose) cout << line << endl << flush;
		read = static_cast<bool>(getline(file, line));
	}
	//int degg = gpoly.size();
	file.close();
	//mpz_clear(c);
	if (verbose) cout << endl << "Complete.  Degree f0 = " << degf << ", degree f1 = " << degg << "." << endl;

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
	int64_t* s0 = new int64_t[degf * nump]();
	int64_t* sieves0 = new int64_t[degf * nump]();
	int* num_s0modp = new int[nump]();
	int64_t* sievep0 = new int64_t[nump]();
	int* sievenum_s0modp = new int[nump]();
	int64_t* s1 = new int64_t[degg * nump]();
	int64_t* sieves1 = new int64_t[degg * nump]();
	int* num_s1modp = new int[nump]();
	int64_t* sievep1 = new int64_t[nump]();
	int* sievenum_s1modp = new int[nump]();
	int64_t itenpc0 = nump / 10;
	int64_t itotal = 0;
	// compute factor base
	if (verbose) cout << endl << "Constructing factor base with " << K << " threads." << endl << flush;
	//if (verbose) cout << endl << "[0%]   constructing factor base..." << endl << flush;
	start = clock();
	#pragma omp parallel
	{
		mpz_t rt; mpz_init(rt); 
		int id = omp_get_thread_num();
		int64_t* stemp0 = new int64_t[degf];
		int64_t* fp = new int64_t[degf+1]();
		int64_t* stemp1 = new int64_t[degg];
		int64_t* gp = new int64_t[degg+1]();

	#pragma omp for
		for (int64_t i = 0; i < nump; i++) {
			int64_t p = primes[i];
			for (int j = 0; j <= degf; j++) fp[j] = mpz_mod_ui(rt, fpoly[j], p);
			int degfp = degf; while (fp[degfp] == 0 || degfp == 0) degfp--;
			int nums0 = polrootsmod(fp, degfp, stemp0, p);
			num_s0modp[i] = nums0;
			for (int j = 0; j < nums0; j++) s0[i*degf + j] = stemp0[j];
			for (int j = 0; j <= degg; j++) gp[j] = mpz_mod_ui(rt, gpoly[j], p);
			int deggp = degg; while (gp[deggp] == 0 || deggp == 0) deggp--;
			int nums1 = polrootsmod(gp, deggp, stemp1, p);
			num_s1modp[i] = nums1;
			for (int j = 0; j < nums1; j++) s1[i*degg + j] = stemp1[j];

	#pragma omp atomic
			itotal++;
			if (itotal % itenpc0 == 0) {
	#pragma omp critical
				if (verbose) cout << "[" << 100 * itotal / nump + 1 << "%]\tConstructing factor base..." << endl << flush;
			}				 
		}

		delete[] gp;
		delete[] stemp1;
		delete[] fp;
		delete[] stemp0;
		mpz_clear(rt);
	}
	timetaken = ( clock() - start ) / (double) CLOCKS_PER_SEC / K;
	start = clock();
	int k0 = 0; int k1 = 0;
	for (int i = 0; i < nump; i++) {
		int nums0 = num_s0modp[i];
		if (nums0 > 0) {
			sievep0[k0] = primes[i];
			for (int j = 0; j < nums0; j++) sieves0[k0*degf + j] = s0[i*degf + j];
			sievenum_s0modp[k0++] = nums0;
		}
		int nums1 = num_s1modp[i];
		if (nums1 > 0) {
			sievep1[k1] = primes[i];
			for (int j = 0; j < nums1; j++) sieves1[k1*degg + j] = s1[i*degg + j];
			sievenum_s1modp[k1++] = nums1;
		}
	}
	timetaken += ( clock() - start ) / (double) CLOCKS_PER_SEC;
	if (verbose) cout << "Complete.  Time taken: " << timetaken << "s" << endl << flush;
	if (verbose) cout << "There are " << k0 << " factor base primes on side 0." << endl << flush;
	if (verbose) cout << "There are " << k1 << " factor base primes on side 1." << endl << flush;

	// write output file
	FILE* out;
	out = fopen(argv[3], "w+");
	fprintf(out, "%d\n", fbb);
	fprintf(out, "%d\n", k0);
	for (int i = 0; i < k0; i++) {
		fprintf(out, "%d", sievep0[i]);
		for (int j = 0; j < sievenum_s0modp[i]; j++)
            fprintf(out, ",%d", sieves0[i*degf + j]);
		fprintf(out, "\n");
	}
	fprintf(out, "%d\n", k1);
	for (int i = 0; i < k1; i++) {
		fprintf(out, "%d", sievep1[i]);
		for (int j = 0; j < sievenum_s1modp[i]; j++)
            fprintf(out, ",%d", sieves1[i*degg + j]);
		fprintf(out, "\n");
	}
	fclose(out);

	// free memory	
	delete[] sievenum_s1modp;
	delete[] sievep1;
	delete[] num_s1modp;
	delete[] s1;
	delete[] sievenum_s0modp;
	delete[] sievep0;
	delete[] num_s0modp;
	delete[] s0;
	mpz_clear(r0);
	delete[] primes;
	delete[] sieve;
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


