#include <stdio.h>
#include <gmpxx.h>
#include <iostream>	// cout
#include <vector>	// vector
#include <sstream>	// istringstream
#include <fstream>	// ifstream
#include <omp.h>
#include <iomanip>	// setprecision
#include <dirent.h>
#include <sys/stat.h>

using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::stringstream;
using std::istringstream;
using std::vector;
using std::getline;
using std::fixed;
using std::setprecision;
using std::ofstream;
using std::flush;

inline void endian_swap(unsigned int& x);

string separator0 = " ";
string separator1 = ":";
string separator2 = ",";

class spM {
	private:
	public:
		int nrows;
		int32_t** valj;	// valuations, jagged array
		uint32_t** colj;	// column position of jth entry in row i
		uint32_t* rw;	   // row weight

		spM(string relsfilename, int numrels, int* ideals) {
			ifstream rels(relsfilename);
			nrows = numrels;
			rw = new uint32_t[nrows];
			valj = new int32_t*[nrows];
			colj = new uint32_t*[nrows];
			string line;
			for (int i = 0; i < nrows; i++) {
				getline(rels, line);
				string Astr = line.substr(0, line.find(separator1));
				line.erase(0, Astr.length() + 1);

				stringstream ss(line);
				stringstream hh;
				vector<uint32_t> relideals;
				
				relideals.clear();
				uint32_t lastcol = -1;
				uint32_t rwi = 0;
				while(ss.good()) {
					string str;
					getline(ss, str, ',');
					int ideal = stoi(str, NULL, 16);
					uint32_t col = ideals[ideal];
					relideals.push_back(col);
					if (col != lastcol) {
						lastcol = col;
						rwi++;
					}
				}
				//sort(relideals.begin(), relideals.end());
				rw[i] = rwi; 
				valj[i] = new int32_t[rwi]();	//  very important to set to zero
				colj[i] = new uint32_t[rwi];
				int t = relideals.size();
				int jj = 0;
				int k = 0;
				while (jj < rwi) {
					uint32_t col = relideals[k];
					colj[i][jj] = col;
					while (relideals[k] == col) {
						valj[i][jj]++;
						k++;
						if (k >= t) break;
					}
					jj++;
				}					
			}
			rels.close();
		}
		~spM() {
			for (int i = 0; i < nrows; i++) {
				delete[] valj[i];
				delete[] colj[i];
			}
			delete[] colj;
			delete[] valj;
			delete[] rw;
		}
};

void spMv(int nrows, int ncols, int nsm, spM* M0, mpz_t* d, mpz_t* w0, mpz_t ell, mpz_t* v0,
	int t_total, int* usedcols);
void vset(int nrows, mpz_t* v, mpz_t* w);
void dot(int nrows, mpz_t* v, mpz_t* w, mpz_t t);
void vscalarprod(int nrows, mpz_t* v, mpz_t* w, mpz_t t);
void vadd(int nrows, mpz_t* w, mpz_t* u, mpz_t* v);
void vsub(int nrows, mpz_t* w, mpz_t* u, mpz_t* v);
void vmod(int nrows, mpz_t* v, mpz_t* w, mpz_t ell);

int main(int argc, char** argv)
{
	if (argc != 10) {
		cout << endl << "Usage: ./sprelvtxt rels_file nrows ideals idealbound rhs.sm "
			"nsm vector.txt num_threads ell" << endl << endl;
		return 0;
	}

	cout << "# ";
	for (int i = 0; i < argc; i++) cout << argv[i] << " ";
	cout << endl;

	bool verbose = false;

	string relsfilename = string(argv[1]);
	int nrows = atoi(argv[2]);
	string idealsfilename = string(argv[3]);
	int maxideal = atoi(argv[4]);

	string line;
	cout << "Reading ideals file..." << flush;
	ifstream idealsfile(idealsfilename);
	getline(idealsfile, line);
	int ncols = atoi(line.substr(2).c_str());
	int* ideals = new int[maxideal+1]();
	int* idealsgood = new int[ncols]();
	for (int i = 0; i < ncols; i++) {
		getline(idealsfile, line);
		int col = atoi(line.substr(0, line.find(separator0)).c_str());
		int ideal = stoi(line.substr(line.find(separator0)+1), NULL, 16);
		ideals[ideal] = col;	// easy lookup
	}
	cout << "done." << endl;
	cout << "Number of ideal columns is " << ncols << endl;
	
	string smfilename(argv[5]);
	int nsm = atoi(argv[6]);
	string vfilename(argv[7]);
	int num_threads = atoi(argv[8]);
	mpz_t ell; mpz_init(ell);
	mpz_set_str(ell, argv[9], 10);
	int numbits = mpz_sizeinbase(ell, 2);

	omp_set_num_threads(num_threads);

	// initialize scalars
	mpz_t t; mpz_init(t);

	// initialize vectors
	cout << "Initializing vectors..." << flush;
	mpz_t* v = new mpz_t[ncols+nsm];
	mpz_t* z = new mpz_t[nrows];
	for (int i = 0; i < nrows; i++) {
		mpz_init_set_ui(v[i], 0);
		mpz_init_set_ui(z[i], 0);
	}
	cout << "done." << endl;

	cout << "Reading input vector..." << flush;
	ifstream vfile(vfilename);
	string vline;
	for (int i = 0; i < ncols+nsm; i++) {
		getline(vfile, vline);
		mpz_set_str(v[i], vline.c_str(), 10);
	}
	cout << "done." << endl;
	cout << "First entry of vector is " << mpz_get_str(NULL, 10, v[0]) << endl;
	cout << "Last entry of vector is " << mpz_get_str(NULL, 10, v[ncols+nsm-1]) << endl;

	// initialize and read matrix constructed from relations
	cout << "Reading matrix into memory..." << flush;
	spM* M0 = new spM(relsfilename, nrows, ideals);
	cout << "done." << endl;

	// initialize and read rhs d (Schirokauer maps).  Note:  sm file is a text file
	cout << "Reading Schirokauer maps..." << flush;
	mpz_t* d = new mpz_t[nrows * nsm];
	ifstream smfile(smfilename);
	getline(smfile, line);	// ignore first line
	vector<string> tokens;
	string token;
	for (int i = 0; i < nrows; i++) {
		getline(smfile, line);
		istringstream iss2(line);
		tokens.clear();
		while (getline(iss2, token, ' ')) {
			if (!token.empty())
				tokens.push_back(token);
		}
		for (int j = 0; j < nsm; j++) {
			mpz_init(d[i*nsm + j]);
			mpz_set_str(d[i*nsm + j], tokens[j].c_str(), 10);
		}
	}
	smfile.close();
	cout << "done." << endl;

	// our matrix M = M0|d


	// compute z = M*v
	cout << "Computing sparse matrix vector product..." << flush;
	int* usedcols = new int[ncols]();
	spMv(nrows, ncols, nsm, M0, d, v, ell, z, num_threads, usedcols);
	// write kernel
	ofstream zfile("zz.txt");
	bool kerneltrue = true;
	for (int i = 0; i < nrows; i++) {
		zfile << mpz_get_str(NULL, 10, z[i]) << endl;
		if (mpz_cmp_ui(z[i], 0) != 0) kerneltrue = false;
	}
	zfile.close();
	cout << "done." << endl;
	if (kerneltrue) cout << "Vector is true kernel element." << endl;
	else cout << "Vector is not a kernel element." << endl;
	cout << "spMv product written to zz.txt" << endl;
	bool allused = true;
	for (int i = 0; i < ncols; i++) if (usedcols[i] == 0) allused = false;
	if (allused) cout << "All logs in kernel used." << endl;

	// release memory
	delete[] usedcols;
	for (int i = 0; i < nrows*nsm; i++) {
		mpz_clear(d[i]);
	}
	delete[] d;
	delete M0;
	mpz_clear(t);
	for (int i = 0; i < nrows; i++) {
		mpz_clear(z[i]);
		mpz_clear(v[i]);
	}
	delete[] z;
	delete[] v;
	mpz_clear(ell);
}

void spMv(int nrows, int ncols, int nsm, spM* M0, mpz_t* d, mpz_t* w0, mpz_t ell, mpz_t* v0,
	int t_total, int* usedcols)
{
	// v0 = (M0|d)*w0 (mod ell)
	#pragma omp parallel
	{
		// v0 = (M0|d)*w0
		#pragma omp for
		for (int i = 0; i < nrows; i++) {
			mpz_set_ui(v0[i], 0);
            int jmax = M0->rw[i];
			for (int j = 0; j < jmax; j++) {
                uint32_t col = M0->colj[i][j];
                int32_t val = M0->valj[i][j];
				if (val >= 0) mpz_addmul_ui(v0[i], w0[col], val);
				else mpz_submul_ui(v0[i], w0[col], -val);
				usedcols[col] = 1;
			}
			for (int j = 0; j < nsm; j++) {
				mpz_addmul(v0[i], w0[ncols + j], d[i*nsm + j]);
			}
			mpz_mod(v0[i], v0[i], ell);
		}
	}
}

void vset(int nrows, mpz_t* v, mpz_t* w)
{
	for (int i = 0; i < nrows; i++) {
		mpz_set(v[i], w[i]);
	}
}

void dot(int nrows, mpz_t* v, mpz_t* w, mpz_t t)
{
	mpz_t m;
	mpz_init(m);

	mpz_set_ui(t, 0);
	for (int i = 0; i < nrows; i++) {
		mpz_mul(m, v[i], w[i]);
		mpz_add(t, t, m);
	}
	
	mpz_clear(m);
}

void vscalarprod(int nrows, mpz_t* v, mpz_t* w, mpz_t t)
{
	mpz_t m;
	mpz_init(m);

	for (int i = 0; i < nrows; i++) {
		mpz_mul(m, w[i], t);
		mpz_set(v[i], m);
	}

	mpz_clear(m);
}

void vadd(int nrows, mpz_t* w, mpz_t* u, mpz_t* v)
{
	for (int i = 0; i < nrows; i++) {
		mpz_add(w[i], u[i], v[i]);
	}
}

void vsub(int nrows, mpz_t* w, mpz_t* u, mpz_t* v)
{
	for (int i = 0; i < nrows; i++) {
		mpz_sub(w[i], u[i], v[i]);
	}
}

void vmod(int nrows, mpz_t* v, mpz_t* w, mpz_t ell)
{
	for (int i = 0; i < nrows; i++) {
		mpz_mod(v[i], w[i], ell);
	}
}

inline void endian_swap(unsigned int& x)
{
	if (x <= 0xFF) return;
	if (x <= 0xFFFF) { x = (x>>8) | ((x & 0xFF)<<8); return; }
	if (x <= 0xFFFFFF) { x = (x>>16) | (x & 0xFF00) | ((x & 0xFF)<<16); return; }
	x = (x>>24) | ((x & 0xFF0000)>>8) | ((x & 0xFF00)<<8) | ((x & 0xFF)<<24);
}

