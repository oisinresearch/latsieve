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
	int numrows;
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
	if (argc != 13) {
		cout << "Usage:  ./sminterface numrels purged_rels index polyfile isofile "
			 << "num_Sch_0 num_Sch_1 ellexp0 ellexp1 sm_out numthreads ell" << endl << endl;
		return 0;
	}

	int numrels = atoi(argv[1]);
	string relsfile = string(argv[2]);
	string indexfile = string(argv[3]);
	string polyfile = string(argv[4]);
	string isofile = string(argv[5]);
	int nsm0 = atoi(argv[6]);
	int nsm1 = atoi(argv[7]);
	int ellexp0 = atoi(argv[8]);
	int ellexp1 = atoi(argv[9]);
	string smoutfile = string(argv[10]);
	int nt = atoi(argv[11]);
	string ellstr = argv[12];

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

	// initialize pari library
	pari_init(2000000000, 65536);

	string loadnfs3dstr = "\\r ~/nfs3d.gp";
	GEN str1 = gp_read_str(loadnfs3dstr.c_str());
	geval(str1);

	// read poly file
	string loadpolystr = "loadpolytnfs(\"" + polyfile + "\");";
	GEN str2 = gp_read_str(loadpolystr.c_str());
	GEN vec1 = geval(str2);
	GEN p = gel(vec1,1);
	GEN f = gel(vec1,2);
	GEN g = gel(vec1,3);
	GEN h = gel(vec1,4);
	GEN f0 = gel(vec1,5);
	GEN g0 = gel(vec1,6);

	// read isofile
	GEN iso1 = readstr(isofile.c_str());
	GEN y3 = geval(gel(iso1,1));
	GEN x3 = geval(gel(iso1,2));
	GEN y4 = geval(gel(iso1,3));
	GEN x4 = geval(gel(iso1,4));

	// set up number fields and J-ideals
	cout << "Constructing number fields and J-ideals in pari-gp..." << endl;
	long prec = 5;
	GEN x = pol_x(0);
	GEN L3 = nfinit0(f, 3, prec);
	GEN L4 = nfinit0(g, 3, prec);
	GEN J3 = idealinv(L3, idealhnf0(L3, gen_1, gsubst(x3, gvar(x3), gel(L3, 2))));
	GEN J4 = idealinv(L4, idealhnf0(L4, gen_1, gsubst(x4, gvar(x4), gel(L4, 2))));

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
	gel(args, 2) = g; 
	gel(args, 3) = y3; 
	gel(args, 4) = x3; 
	gel(args, 5) = y4; 
	gel(args, 6) = x4; 

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
		data[i].index = &index;
		data[i].rels = &rels;
		data[i].smout = &smout;
		data[i].nt = nt;
		data[i].numrels = numrels;
		data[i].numrows = numrows;
		data[i].nsm0 = nsm0;
		data[i].nsm1 = nsm1;
		data[i].ellexp0 = ellexp0;
		data[i].ellexp1 = ellexp1;
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
	GEN g = gel(args, 2); 
	GEN y3 = gel(args, 3); 
	GEN x3 = gel(args, 4); 
	GEN y4 = gel(args, 5); 
	GEN x4 = gel(args, 6); 
	GEN I3, I4, I3M, I4M;	
	GEN g3, g4;
	GEN Ip, pj, idj, fj, sj;
	GEN item0 = cgetg(5, t_VEC);
	long i0; GEN logi0; GEN dims;

	// set up number fields and J-ideals
	GEN x = pol_x(0);
	GEN L3 = nfinit(f, prec);
	GEN L4 = nfinit(g, prec);
	GEN J3 = idealinv(L3, idealhnf0(L3, gen_1, gsubst(x3, gvar(x3), gel(L3, 2))));
	GEN x3inv = gsubst(x3, gvar(x3), gel(L3, 2));
	long x3var = gvar(x3);
	GEN L3iso = gel(L3,2);
	//GEN J4 = idealinv(L4, idealhnf0(L4, gen_1, gsubst(x4, gvar(x4), gel(L4, 2))));
	
	// set up Schirokauer map data
	GEN ell = gp_read_str(data->ellstr.c_str());
	GEN eps0 = gsubgs(gpowgs(ell, data->ellexp0), 1);
	GEN eps1 = gsubgs(gpowgs(ell, data->ellexp1), 1);
	int ly3 = glength(y3);
	int lx3 = glength(x3);
	GEN v3 = cgetg(3, t_VEC);
	GEN vy3 = cgetg(ly3+1, t_VEC);
	GEN vx3 = cgetg(lx3+1, t_VEC);
	for (int j = 1; j <= ly3; ++j)
		gel(vy3, j) = denominator(polcoef(y3, j-1, -1), NULL);
	for (int j = 1; j <= lx3; ++j)
		gel(vx3, j) = denominator(polcoef(y3, j-1, -1), NULL);
	gel(v3, 1) = vy3;
	gel(v3, 2) = vx3;
	GEN m1 = glcm0(gconcat1(v3), NULL);
	GEN v4 = cgetg(3, t_VEC);
	int ly4 = glength(y4);
	int lx4 = glength(x4);
	GEN vy4 = cgetg(ly4+1, t_VEC);
	GEN vx4 = cgetg(lx4+1, t_VEC);
	for (int j = 1; j <= ly4; ++j)
		gel(vy4, j) = denominator(polcoef(y4, j-1, -1), NULL);
	for (int j = 1; j <= lx4; ++j)
		gel(vx4, j) = denominator(polcoef(y4, j-1, -1), NULL);
	gel(v4, 1) = vy4;
	gel(v4, 2) = vx4;
	GEN m2 = glcm0(gconcat1(v4), NULL);
	
	GEN a, b, c, d;
	GEN T1, T2, A1, A2;
	GEN feps, geps;
	vector<string> indexi;
	int mark = 1024;
	int i = data->startrow;
	while (i < data->numrows) {
		pari_sp ltop = avma;
		A1 = gen_1;
		A2 = gen_1;
		indexi = split((*data->index)[i], ' ');
		int n = stoi(indexi[0]);
		for (int j = 1; j <= n; j++) {
			long r = stol(indexi[j].substr(0, indexi[j].find_first_of(":")), NULL, 16);
			int v = stoi(indexi[j].substr(indexi[j].find_first_of(":")+1));
			string rel = (*data->rels)[r];
			if (rel[rel.length()-1] == ',') rel = rel.substr(0, rel.length()-1);
			string ABhex = rel.substr(0, rel.find(":"));
			string Ahex = ABhex.substr(0, ABhex.find(","));
			string Bhex = ABhex.substr(ABhex.find(",")+1);
			int64_t A = stol(Ahex, NULL, 16);
			int64_t B = stol(Bhex, NULL, 16);
			
			int64_t a1 = (A>>24) - (1l<<23);
			int64_t b1 = (A&((1l<<24)-1)) - (1l<<23);
			int64_t c1 = (B>>24) - (1l<<23);
			int64_t d1 = (B&((1l<<24)-1)) - (1l<<23);

			a = stoi(a1);
			b = stoi(b1);
			c = stoi(c1);
			d = stoi(d1);
			
			T1 = gmul(m1, gadd(gadd(a, gmul(b, y3)), gmul(gadd(c, gmul(d, y3)), x3)));
			T2 = gmul(m2, gadd(gadd(a, gmul(b, y4)), gmul(gadd(c, gmul(d, y4)), x4)));
			A1 = liftall(gmul(A1, gpowgs(gmodulo(gmodulo(T1, f), gsqr(ell)), v)));
			A2 = liftall(gmul(A2, gpowgs(gmodulo(gmodulo(T2, g), gsqr(ell)), v)));
		}
		//std::this_thread::sleep_for(std::chrono::milliseconds(2));
		string smout = "";
		feps = gdiv(liftall(gsubgs(gpow(gmodulo(gmodulo(A1, f), gsqr(ell)), eps0, prec), 1)), ell);
		for (int i = 0; i < data->nsm0; i++) {
			char* sm_char = GENtostr(polcoef(feps, i+1, -1));
			string sm_exp_i = string(sm_char);
			smout = smout + sm_exp_i + " ";
			free(sm_char);
		}
		geps = gdiv(liftall(gsubgs(gpow(gmodulo(gmodulo(A2, g), gsqr(ell)), eps1, prec), 1)), ell);
		for (int i = 0; i < data->nsm1; i++) {
			char* sm_char = GENtostr(polcoef(geps, i+1, -1));
			string sm_exp_i = string(sm_char);
			smout = smout + sm_exp_i;
			free(sm_char);
			if (i < data->nsm1 - 1) smout = smout + " ";
		}
		(*data->smout)[i] = smout;
		indexi.clear();
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

