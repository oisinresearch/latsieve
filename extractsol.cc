#include <stdio.h>
#include <gmpxx.h>
#include <iostream>	// cout
#include <vector>	// vector
#include <sstream>	// stringstream
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
using std::vector;
using std::getline;
using std::fixed;
using std::setprecision;
using std::ofstream;
using std::flush;
using std::hex;

inline void endian_swap(unsigned int& x);
string long2hex(long n);

class spM {
	private:
	public:
		int nrows;
		int32_t** valj;	// valuations, jagged array
		uint32_t** colj;	// column position of jth entry in row i
		uint32_t* rw;	   // row weight

		spM(int nrows1, string M0filename) {
			ifstream M0file(M0filename);
			nrows = nrows1;
			rw = new uint32_t[nrows];
			valj = new int32_t*[nrows];
			colj = new uint32_t*[nrows];
			for (int i = 0; i < nrows; i++) {
				// read ith row weight
				M0file.read((char*)&rw[i], sizeof(rw[i]));
				valj[i] = new int32_t[rw[i]];
				colj[i] = new uint32_t[rw[i]];
				for (uint32_t j = 0; j < rw[i]; j++) {
					M0file.read((char*)&colj[i][j], sizeof(colj[i][j]));
					M0file.read((char*)&valj[i][j], sizeof(valj[i][j]));
				}
			}
			M0file.close();
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

void spMv(int nrows, int nsm, spM* M0, mpz_t* d, mpz_t* w0, mpz_t ell, mpz_t* v0,
	int t_total);
void vset(int nrows, mpz_t* v, mpz_t* w);
void dot(int nrows, mpz_t* v, mpz_t* w, mpz_t t);
void vscalarprod(int nrows, mpz_t* v, mpz_t* w, mpz_t t);
void vadd(int nrows, mpz_t* w, mpz_t* u, mpz_t* v);
void vsub(int nrows, mpz_t* w, mpz_t* u, mpz_t* v);
void vmod(int nrows, mpz_t* v, mpz_t* w, mpz_t ell);

int main(int argc, char** argv)
{
	if (argc != 13) {
		cout << endl << "Usage: ./extractsol rels_file nrels Kbad Z output_rels "
			"nsm input_sm output_sm Kgood idealbound ideals_file ideals_out" << endl << endl;
		return 0;
	}

	cout << "# ";
	for (int i = 0; i < argc; i++) cout << argv[i] << " ";
	cout << endl;

	bool verbose = false;

	string line;
	string separator0 = " ";
	string separator1 = ":";
	string separator2 = ",";

	string relsfilename = string(argv[1]);
	int nrows = atoi(argv[2]);
	string Kbadfilename = string(argv[3]);
	string Zfilename = string(argv[4]);
	string outrelsfilename = string(argv[5]);
	int nsm = atoi(argv[6]);
	string smfilename = string(argv[7]);
	string outsmfilename = string(argv[8]);
	string Koutfilename = string(argv[9]);
	int maxideal = atoi(argv[10]);
	string idealsfilename = string(argv[11]);
	string idealsoutfilename = string(argv[12]);
	
	cout << "Reading ideals file..." << flush;
	ifstream idealsfile(idealsfilename);
	getline(idealsfile, line);
	int ncols = atoi(line.substr(2).c_str());
	int* idlookup = new int[maxideal+1]();
	int* idealsgood = new int[ncols]();
	int* ideals = new int[ncols]();
	for (int i = 0; i < ncols; i++) {
		getline(idealsfile, line);
		int col = atoi(line.substr(0, line.find(separator0)).c_str());
		int ideal = stoi(line.substr(line.find(separator0)+1), NULL, 16);
		idlookup[ideal] = col;	// easy lookup
		ideals[i] = ideal;
	}
	cout << "done." << endl;
	
	// initialize vectors
	mpz_t* Kbad = new mpz_t[ncols+nsm];
	for (int i = 0; i < nrows; i++) {
		mpz_init_set_ui(Kbad[i], 0);
	}

	cout << "Reading Kbad vector..." << flush;
	string vline;
	ifstream vfile0(Kbadfilename);
	int numgood = 0;
	for (int i = 0; i < ncols+nsm; i++) {
		getline(vfile0, vline);
		mpz_set_str(Kbad[i], vline.c_str(), 10);
	}
	cout << "done." << endl;
	cout << "Reading Z vector..." << flush;
	ifstream vfile1(Zfilename);
	mpz_t Zi; mpz_init_set_ui(Zi, 0);
	int* rowgood = new int[nrows]();
	for (int i = 0; i < nrows; i++) {
		getline(vfile1, vline);
		mpz_set_str(Zi, vline.c_str(), 10);
		if (mpz_cmp_ui(Zi, 0) == 0) {
			numgood++;
			rowgood[i] = 1;
		}
	}
	cout << "done." << endl;

	// read relations file
	int  num = 0; int mark = 1024;
	ifstream rels(relsfilename);
	cout << "Reading relations..." << endl;
	while (std::getline(rels, line)) { //rels >> line) {
		string line0 = line;

		string Astr = line.substr(0, line.find(separator1));
		line.erase(0, Astr.length() + 1);

		stringstream ss(line);
		stringstream hh;
		vector<int> relideals;

		while(ss.good()) {
			string str;
			getline(ss, str, ',');
			int ideal = stoi(str, NULL, 16);
			relideals.push_back(ideal);
		}

		int numideals = relideals.size();
		sort(relideals.begin(), relideals.end());

		if (rowgood[num] == 1) {
			for (int i = 0; i < numideals; i++) {
				idealsgood[idlookup[relideals[i]]] = 1;
			}
		}
		
		num++;
		if (num % mark == 0) {
			cout << num << " relations read..." << endl;
			mark *= 2;
		}
	}
	cout << "all relations read." << endl;

	// count good ideals
	int numgoodideals = 0;
	for (int i = 0; i < ncols; i++) if (idealsgood[i] == 1) numgoodideals++;

	// write output relations
	num = 0; mark = 1024;
	int numwritten = 0;
	ifstream rels2(relsfilename);
	ofstream outrels(outrelsfilename);
	cout << "Writing subset of relations with good logs..." << endl;
	while (std::getline(rels2, line)) { //rels >> line) {

		if (rowgood[num] == 1) {
			outrels << line << endl;
			numwritten++;
		}
		
		num++;
		if (num % mark == 0) {
			cout << num << " relations processed..." << endl;
			mark *= 2;
		}
	}
	cout << numwritten << " relations written." << endl;
	outrels.close();

	// write kernel
	cout << "Writing output vector Kgood..." << flush;
	ofstream Koutfile(Koutfilename);
	for (int i = 0; i < ncols; i++) {
		if (idealsgood[i] == 1) {
			Koutfile << mpz_get_str(NULL, 10, Kbad[i]) << endl;
		}
	}
	// write last nsm logs regardless
	for (int i = 0; i < nsm; i++) {
		Koutfile << mpz_get_str(NULL, 10, Kbad[ncols + i]) << endl;
	}
	Koutfile.close();
	cout << "done." << endl;
	
	// write output ideals
	cout << "Writing output ideals..." << flush;
	int ii = 0;
	ofstream idealsout(idealsoutfilename);
	idealsout << "# " << numgoodideals << endl;
	for (int id = 0; id < ncols; id++) {
		if (idealsgood[id] == 1) {
			idealsout << ii << " " << long2hex(ideals[id]) << endl;
			ii++;
		}
	}
	idealsout.close();
	cout << "done." << endl;

	// initialize and read rhs (Schirokauer maps).  Note:  sm file is a text file
	cout << "Reading Schirokauer maps..." << flush;
	ifstream smfile(smfilename);
	ofstream outsm(outsmfilename);
	cout << "Writing subset of Schirokauer maps with good logs..." << endl;
	num = 0;
	// write first line
	getline(smfile, line);
	outsm << numgood << " " << line.substr(line.find(separator0)+1) << endl;
	for (int i = 0; i < nrows; i++) {
		getline(smfile, line);
		if (rowgood[num] == 1) {
			outsm << line << endl;
		}
		
		num++;
		if (num % mark == 0) {
			cout << num << " relations processed..." << endl;
			mark *= 2;
		}
	}
	smfile.close();
	outsm.close();
	cout << "done." << endl;
	
	// release memory
	delete[] idealsgood;
	delete[] ideals;
	delete[] idlookup;
	for (int i = 0; i < nrows; i++) {
		mpz_clear(Kbad[i]);
	}
	delete[] Kbad;
	mpz_clear(Zi);
	delete[] rowgood;
}

void spMv(int nrows, int nsm, spM* M0, mpz_t* d, mpz_t* w0, mpz_t ell, mpz_t* v0,
	int t_total)
{
	cout << endl;
	// v0 = (M0|d)*w0 (mod ell)
	//omp_set_num_threads(t);
	#pragma omp parallel
	{
		// v0 = (M0|d)*w0
		#pragma omp for
		for (int i = 0; i < nrows; i++) {
			mpz_set_ui(v0[i], 0);
            int jmax = M0->rw[i];
			//if (M0->colj[i][0]!=0) cout << i << ":" << M0->colj[i][0] << endl;
			for (int j = 0; j < jmax; j++) {
                uint32_t col = M0->colj[i][j];
                int32_t val = M0->valj[i][j];
				if (i==1&&val<0) cout << col << ":" << mpz_get_str(NULL, 10, w0[col]) << endl;
				//if (i<=1 && val>0) cout << col << "," << flush;
				if (val >= 0) mpz_addmul_ui(v0[i], w0[col], val);
				else mpz_submul_ui(v0[i], w0[col], -val);
			}
			for (int j = 0; j < nsm; j++) {
				mpz_addmul(v0[i], w0[nrows - nsm + j], d[i*nsm + j]);
			}
			mpz_mod(v0[i], v0[i], ell);
			if (i==1) cout << endl;
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

string long2hex(long n)
{
	string str1;
	stringstream stream;
	stream.str("");
	stream << hex << n;
	stream >> str1;
	return str1;
}

