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
//#include "mpz_poly.h"
#include <pari/pari.h>
#include <chrono>
#include <thread>

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
using std::stoi;
auto npos = std::string::npos;

typedef struct struct2 {
	void* RETURNVALUE;
	struct pari_thread pth;
	string ellstr;
	int nt;
	int numrels;
	int nsm0;
	int nsm1;
	int ellexp0;
	int ellexp1;
	int startrow;
	vector<string>* index;
	vector<string>* rels;
	vector<string>* smout;
} thread_root_data;

void *thread_root(void* context_data);
long hex2long(string str1);
string long2hex(long n);
vector<string> split(string text, char delim);

int main(int argc, char** argv)
{
	if (argc != 9) {
		cout << "Usage:  ./sm2dmono numrels purged_rels polyfile "
			 << "num_Sch_0 ellexp0 sm_out numthreads ell" << endl << endl;
		return 0;
	}

	int numrels = atoi(argv[1]);
	string relsfile = string(argv[2]);
	string polyfile = string(argv[3]);
	int nsm0 = atoi(argv[4]);
	int ellexp0 = atoi(argv[5]);
	string smoutfile = string(argv[6]);
	int nt = atoi(argv[7]);
	string ellstr = argv[8];

	string line;
	char linebuffer[100];
	int mark = 1024;
	int t = 0;

	// initialize pari library
	pari_init(2000000000, 65536);

	string loadnfs3dstr = "\\r ~/nfs3d.gp";
	GEN str1 = gp_read_str(loadnfs3dstr.c_str());
	geval(str1);

	// read poly file
	string loadpolystr = "loadpolymono(\"" + polyfile + "\");";
	GEN str2 = gp_read_str(loadpolystr.c_str());
	GEN vec1 = geval(str2);
	GEN p = gel(vec1,1);
	GEN f = gel(vec1,2);

	// set up number fields and J-ideals
	cout << "Constructing number fields and J-ideals in pari-gp..." << endl;
	long prec = 5;
	GEN x = pol_x(0);
	GEN L3 = nfinit0(f, 3, prec);
	GEN x3 = gel(L3,2);
	cout << GENtostr(x3) << endl;
	GEN J3 = idealinv(L3, idealhnf0(L3, gen_1, x3));

	// read relations file
	cout << "Reading relations file..." << endl;
	vector<string> rels;
	ifstream file2(relsfile);
	mark = 1024;
	t = 0;
	while (getline(file2, line)) {
		rels.push_back(line);
		t++;
		if (t % mark == 0) {
			cout << to_string(t) << " relations read..." << endl;
			mark *= 2;
		}
	}
	if (numrels != t) {
		cout << "number of relations is incorrect.  Exiting..." << endl;
		return 1;
	}
	cout << to_string(t) << " relations read." << endl << endl;
	file2.close();
	
	GEN args = cgetg(7, t_VEC);
	gel(args, 1) = f; 
	gel(args, 2) = x3; 

	// begin threading code
	thread_root_data* data = new thread_root_data[nt];

	vector<string> smout;
	for (int i = 0; i < numrels; i++) {
		smout.push_back("");
	}

	pthread_t* th = new pthread_t[nt];
	for (int i = 0; i < nt; i++) {
		pari_thread_alloc(&data[i].pth, 2000000000, args);
	}

	cout << "Computing Schirokauer maps..." << endl;
	// create nt threads of work
	for (int i = 0; i < nt; i++) {
		data[i].ellstr = ellstr;
		data[i].rels = &rels;
		data[i].smout = &smout;
		data[i].nt = nt;
		data[i].numrels = numrels;
		data[i].nsm0 = nsm0;
		data[i].ellexp0 = ellexp0;
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
	mark = 1024;
	cout << endl << "Writing output Schirokauer maps to " << smoutfile << "..." << endl;
	line = to_string(numrels) + " " + to_string(nsm0) + " " + ellstr;
	fprintf(out, "%s\n", line.c_str());
	for (int i = 0; i < numrels; i++) {
		fprintf(out, "%s\n", smout[i].c_str());
		t++;
		if (t % mark == 0) {
			cout << to_string(t) << " Schirokauer maps written..." << endl;
			mark *= 2;
		}
	}
	int nummaprowswritten = t;
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
	
	// start pari thread
	GEN args = pari_thread_start(&data->pth);

	pari_sp ltop = avma;
	long prec = 5;
	GEN f = gel(args, 1); 
	GEN x3 = gel(args, 2); 
	GEN I3, I3M;	
	GEN g3;
	GEN Ip, pj, idj, fj, sj;
	GEN item0 = cgetg(5, t_VEC);
	long i0; GEN logi0; GEN dims;

	// set up number fields and J-ideals
	GEN L3 = nfinit(f, prec);
	GEN J3 = idealinv(L3, idealhnf0(L3, gen_1, x3));
	GEN x3inv = gsubst(x3, gvar(x3), gel(L3, 2));
	//long x3var = gvar(x3);
	//GEN L3iso = gel(L3,2);
	
	// set up Schirokauer map data
	GEN ell = gp_read_str(data->ellstr.c_str());
	GEN eps0 = gsubgs(gpowgs(ell, data->ellexp0), 1);
	
	GEN a, b;
	GEN T1, A1;
	GEN feps;
	int mark = 1024;
	int i = data->startrow;
	while (i < data->numrels) {
		pari_sp ltop = avma;
		A1 = gen_1;
		string rel = (*data->rels)[i];
		if (rel[rel.length()-1] == ',') rel = rel.substr(0, rel.length()-1);
		string ABhex = rel.substr(0, rel.find(":"));
		string Ahex = ABhex.substr(0, ABhex.find(","));
		string Bhex = ABhex.substr(ABhex.find(",")+1);
		int64_t A = stol(Ahex, NULL, 16);
		int64_t B = stol(Bhex, NULL, 16);
		
		int64_t a1 = A - (1l<<46);
		int64_t b1 = B - (1l<<46);

		a = stoi(a1);
		b = stoi(b1);
		
		T1 = gadd(a, gmul(b, x3));
		//A1 = liftall(gmul(A1, gpowgs(gmodulo(gmodulo(T1, f), gsqr(ell)), v)));
		A1 = gmul(A1, T1);
		//std::this_thread::sleep_for(std::chrono::milliseconds(2));
		string smout = "";
		feps = gdiv(liftall(gsubgs(gpow(gmodulo(gmodulo(A1, f), gsqr(ell)), eps0, prec), 1)), ell);
		for (int i = 0; i < data->nsm0; i++) {
			char* sm_char = GENtostr(polcoef(feps, i+1, -1));
			string sm_exp_i = string(sm_char);
			smout = smout + sm_exp_i + " ";
			free(sm_char);
		}
		(*data->smout)[i] = smout;
		i += data->nt;
		if (i > mark && data->startrow == 0) {
			cout << i << " Schirokauer maps computed..." << endl;
			mark *= 2;
		}
		avma = ltop;
	}

	// tidy up
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

