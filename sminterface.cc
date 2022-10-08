#include <stdint.h>	// int64_t
#include <iostream> // cout
#include <gmpxx.h>
#include <math.h>	// sqrt
#include <fstream>	// file
#include <ctime>	// clock_t
#include <sstream>	// stringstream
#include <stack>	// stack
#include <string>
#include <assert.h>
#include <vector>
#include <pthread.h>

using std::cout;
using std::endl;
using std::flush;
using std::ifstream;
using std::string;
using std::pair;
using std::make_pair;
using std::stringstream;
using std::vector;
using std::hex;
using std::to_string;
using std::sort;
auto npos = std::string::npos;

typedef struct struct2 {
	void* RETURNVALUE;
	string ellstr;
	int nt;
	int numrels;
	int numrows;
	int nsm;
	int startrow;
	vector<string>* index;
	vector<string>* smin;
	vector<string>* smout;
} thread_root_data;

void *thread_root(void* context_data);
long hex2long(string str1);
string long2hex(long n);
vector<string> split(string text, char delim);

int main(int argc, char** argv)
{
	if (argc != 9) {
		cout << "Usage:  ./sminterface numrels index num_Sch_0 num_Sch_1 "
			 << "sm_in sm_out numthreads ell" << endl << endl;
		return 0;
	}

	int numrels = atoi(argv[1]);
	string indexfile = string(argv[2]);
	int nsm0 = atoi(argv[3]);
	int nsm1 = atoi(argv[4]);
	string sminfile = string(argv[5]);
	string smoutfile = string(argv[6]);
	int nt = atoi(argv[7]);
	string ellstr = argv[8];

	// read index
	vector<string> index;
	cout << "Reading index file..." << endl;
	string line;
	char linebuffer[100];
	ifstream file1(indexfile);
	getline(file1, line);	// first line contains number of rows
	int numrows = atoi(line.c_str()); 
	int mark = 1024;
	int t = 0;
	while (getline(file1, line)) {
		index.push_back(line);
		t++;
		if (t % mark == 0) {
			cout << to_string(t) << " index rows read..." << endl;
			mark *= 2;
		}
	}
	if (numrows != t) {
		cout << "number of index rows is incorrect.  Exiting..." << endl;
		return 1;
	}
	cout << t << " index rows read." << endl << endl;
	file1.close();

	// read Schirokauer maps
	vector<string> smin;
	cout << "Reading Schirokauer maps..." << endl;
	ifstream file2(sminfile);
	getline(file2, line);	// first line contains number of rows
	vector<string> line0b = split(line, ' ');
	mark = 1024;
	t = 0;
	while (getline(file2, line)) {
		smin.push_back(line);
		t++;
		if (t % mark == 0) {
			cout << to_string(t) << " Schirokauer maps read..." << endl;
			mark *= 2;
		}
	}
	if (numrels != t) {
		cout << "number of relations is incorrect.  Exiting..." << endl;
		return 1;
	}
	cout << t << " Schirokauer maps read." << endl << endl;
	file2.close();

	thread_root_data* data = new thread_root_data[nt];
	pthread_t* th = new pthread_t[nt];
	vector<string> smout;
	for (int i = 0; i < numrels; i++) {
		smout.push_back("");
	}

	// main loop
	// create nt threads of work
	for (int i = 0; i < nt; i++) {
		data[i].ellstr = ellstr;
		data[i].index = &index;
		data[i].smin = &smin;
		data[i].smout = &smout;
		data[i].nt = nt;
		data[i].numrels = numrels;
		data[i].numrows = numrows;
		data[i].nsm = nsm0 + nsm1;
		data[i].startrow = i;
		pthread_create(&th[i], NULL, &thread_root, (void*)&data[i]);
	}

	// join nt threads
	for (int i = 0; i < nt; i++) {
		pthread_join(th[i], (void**)&data[i]); // this is blocking
	}

	// write output file
	FILE* out;
	out = fopen(smoutfile.c_str(), "w+");
	t = 0;
	cout << endl << "Writing output Schirokauer maps to " << smoutfile << "..." << endl;
	line = to_string(numrows) + " " + to_string(nsm0+nsm1) + " " + ellstr;
	fprintf(out, "%s\n", line.c_str());
	for (int i = 0; i < numrows; i++) {
		fprintf(out, "%s\n", smout[i].c_str());
		t++;
		if (t % mark == 0) {
			cout << to_string(t) << " Schirokauer maps written..." << endl;
			mark *= 2;
		}
	}
	int numidealhits = t;
	cout << to_string(t) << " Schirokauer maps written." << endl;
	
	cout << endl << "Finished!" << endl;

	// tidy up
	delete[] th;
	delete[] data;

	return 0;
}

void *thread_root(void* context_data)
{
	thread_root_data* data = (thread_root_data*)context_data;
	mpz_t ell; mpz_init(ell);
	mpz_set_str(ell, data->ellstr.c_str(), 10);

	mpz_t* s = new mpz_t[data->nsm];
	for (int k = 0; k < data->nsm; k++) mpz_init(s[k]);
	mpz_t t; mpz_init(t);
	vector<string> indexi;
	vector<string> smini;
	int mark = 1024;
	int i = data->startrow;
	while (i < data->numrows) {
		indexi = split((*data->index)[i], ' ');
		int n = atoi(indexi[0].c_str());
		for (int k = 0; k < data->nsm; k++) mpz_set_ui(s[k], 0);
		for (int j = 1; j <= n; j++) {
			long r = hex2long(indexi[j].substr(0, indexi[j].find_first_of(":")));
			int v = atoi(indexi[j].substr(indexi[j].find_first_of(":")+1).c_str());
			smini = split((*data->smin)[r], ' ');
			for (int k = 0; k < data->nsm; k++) {
				mpz_set_str(t, smini[k].c_str(), 10);
				mpz_mul_si(t, t, v);
				mpz_add(s[k], s[k], t);
			}
			smini.clear();
		}
		string smout = "";
		for (int k = 0; k < data->nsm; k++) {
			mpz_mod(s[k], s[k], ell);
			smout = smout + string(mpz_get_str(NULL, 10, s[k]));
			if (k < data->nsm - 1) smout = smout + " ";
		}
		(*data->smout)[i] = smout;
		indexi.clear();
		i += data->nt;
		if (i > mark && data->startrow == 0) {
			cout << i << " Schirokauer maps computed..." << endl;
			mark *= 2;
		}
		
	}

	// tidy up
	mpz_clear(t);
	for (int k = 0; k < data->nsm; k++) mpz_clear(s[k]);
	delete[] s;
	mpz_clear(ell);

	return NULL;
}

long hex2long(string str1)
{
	stringstream stream;
	long jcol = -1;
	stream.str("");
	stream << str1;
	stream >> hex >> jcol;
	return jcol;
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

vector<string> split(string text, char delim)
{
    string line;
    vector<string> vec;
    stringstream ss(text);
    while(getline(ss, line, delim)) {
        vec.push_back(line);
    }
    return vec;
}

