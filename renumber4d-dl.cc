#include <stdint.h>	// int64_t
#include <iostream> // cout
#include <gmpxx.h>
#include "intpoly.h"
#include <math.h>	// sqrt
#include <fstream>	// file
#include <ctime>	// clock_t
#include "mpz_poly.h"
#include <sstream>	// stringstream
#include <stack>	// stack
#include <string>
#include "factorsmall.h"
#include <assert.h>
#include <pthread.h>
#include <pari/pari.h>
#include <vector>

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
using std::lower_bound;
using std::distance;
auto npos = std::string::npos;

typedef struct struct2 {
	void* RETURNVALUE;
	int nt;
	struct pari_thread pth;
	vector<string>* rels;
	vector<pair<string,int>>* ideals;
	vector<string> rnrels;
	int startrow;
	int endrow;
} thread_root_data;

void *thread_root(void* context_data);
long hex2long(string str1);
string long2hex(long n);

int main(int argc, char** argv)
{
	if (argc != 7) {
		cout << "Usage:  ./renumber4d-dl polyfile isofile renumberfile relationsfile ";
		cout << "numthreads outputfile" << endl << endl;
		return 0;
	}

	string polyfile = string(argv[1]);
	string isofile = string(argv[2]);
	string renumberfile = string(argv[3]);
	string relsfile = string(argv[4]);
	int nt = atoi(argv[5]);
	string outputfile = argv[6];
	
	string loadnfs3dstr = "\\r ~/nfs3d.gp";
	string loadpolystr = "loadpolytnfs(\"" + polyfile + "\");";

	// initialize pari library
	pari_init(2000000000, 65536);

	GEN str1 = gp_read_str(loadnfs3dstr.c_str());
	geval(str1);
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

	// read factor base (renumber file)
	cout << "Reading factor base (renumber file)..." << endl;
	int nideals = 0;
	vector<pair<string,int>> ideals;
	ifstream file1(renumberfile);
	int mark = 1024;
	int t = 0;
	string line;
	while (getline(file1, line)) {
		ideals.push_back(make_pair(line, t));
		t++;
		if (t >= mark) {
			cout << to_string(t) << " ideals read..." << endl;
			mark *= 2;
		}
	}
	cout << to_string(t) << " ideals read." << endl << endl;
	long numideals = t;
	file1.close();

	cout << "Sorting factor base as vector of indexed strings..." << endl;
	sort(ideals.begin(), ideals.end());
			
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
	int nrels = 0;
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
	cout << to_string(t) << " relations read." << endl << endl;
	nrels = t;
	file2.close();
	
	// begin threading code
	thread_root_data* data = new thread_root_data[nt];
	for (int i = 0; i < nt; i++) {
		data[i].nt = nt;
		data[i].rels = &rels;
		data[i].ideals = &ideals;
	}

	GEN args = cgetg(7, t_VEC);
	gel(args, 1) = f; 
	gel(args, 2) = g; 
	gel(args, 3) = y3; 
	gel(args, 4) = x3; 
	gel(args, 5) = y4; 
	gel(args, 6) = x4; 

	// main loop
	pthread_t* th = new pthread_t[nt];
	for (int i = 0; i < nt; i++) {
		pari_thread_alloc(&data[i].pth, 2000000000, args);
	}

	int gap = (nrels/nt) + 1;
	int iters = 1000;
	for (int i = 0; i < nt; i++) {
		data[i].startrow = i*gap;
		data[i].endrow = i*gap+iters-1;
		if (data[i].endrow >= (i+1)*gap) data[i].endrow = (i+1)*gap - 1;
		if (data[i].endrow >= nrels) data[i].endrow = nrels - 1;
	}

	t = 0;
	mark = 1024;
	cout << "Renumbering relations..." << endl;
	while (t < nrels) {
		// create nt threads of work
		for (int i = 0; i < nt; i++) {
			pthread_create(&th[i], NULL, &thread_root, (void*)&data[i]);
		}

		// join nt threads
		for (int i = 0; i < nt; i++) {
			pthread_join(th[i], (void**)&data[i]); // this is blocking
		}

		// update factor base logs with new logs
		t = 0;
		for (int i = 0; i < nt; i++) {
			t += data[i].rnrels.size();
		}
		if (t > mark) {
			cout << to_string(t) << " relations renumbered..." << endl;
			mark *= 2;
		}

		// continue only if t < nrels
		if (t < nrels) {
			for (int i = 0; i < nt; i++) {
				data[i].startrow += iters;
				data[i].endrow += iters;
				if (data[i].endrow >= (i+1)*gap) data[i].endrow = (i+1)*gap - 1;
				if (data[i].endrow >= nrels) data[i].endrow = nrels - 1;
			}
		}
	}
	cout << to_string(t) << " relations renumbered." << endl;

	// write output file
	FILE* out;
	out = fopen(outputfile.c_str(), "w+");
	t = 0;
	cout << endl << "Writing relations to output file " << outputfile << "..." << endl;
	for (int i = 0; i < nt; i++) {
		for (int j = 0; j < data[i].rnrels.size(); j++) {
			fprintf(out, "%s\n", data[i].rnrels[j].c_str());
			t++;			
			if (t > mark) {
				cout << to_string(t) << " renumbered relations written..." << endl;
				mark *= 2;
			}
		}
	}
	cout << to_string(t) << " renumbered relations written." << endl;
	fclose(out);

	cout << endl << "Finished!" << endl;

	// tidy up
	for (int i = 0; i < nt; i++) pari_thread_free(&data[i].pth);
	delete[] th;
	pari_close();

	return 0;
}

void *thread_root(void* context_data)
{
	thread_root_data* data = (thread_root_data*)context_data;

	vector<pair<string,int>>* list = data->ideals;

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
	
	int mark = 1024;
	int row = data->startrow;
	while (row <= data->endrow) {
		string rel = (*data->rels)[row];
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

		GEN a = stoi(a1);
		GEN b = stoi(b1);
		GEN c = stoi(c1);
		GEN d = stoi(d1);

		vector<int> idealcols;
		idealcols.clear();
		idealcols.push_back(0);	// J ideal

		// side 0
		g3 = gadd(gadd(a, gmul(b, y3)), gmul(gadd(c, gmul(d, y3)), x3));
		I3 = idealmul(L3, idealhnf0(L3, gsubst(g3, gvar(g3), gel(L3, 2)), NULL), J3);
		I3M = idealfactor(L3, I3);
		dims = matsize(I3M);
		int nv = itos(gel(dims,1));
		for (int j = 0; j < nv; ++j) {
			Ip = gcoeff(I3M, j+1, 1);
			int v = itos(gcoeff(I3M, j+1, 2));
			pj = gel(Ip, 1);
			idj = gel(Ip, 2);
			fj = gel(Ip, 4);
			sj = gen_0;
			char* pjstr = GENtostr(pj);
			char* idjstr = GENtostr(idj);
			char* fjstr = GENtostr(fj);
			string ideal = "[[" + string(pjstr) + ", " + string(idjstr) + "], "
				+ string(fjstr) + ", 0]";
			free(fjstr); free(idjstr); free(pjstr);
			auto lb = lower_bound(list->begin(), list->end(), make_pair(ideal, 0),
        		[](pair<string,int> a, pair<string,int> b) { return a.first < b.first; });
    		if (lb->first == ideal) {
				int i0 = lb->second;
				idealcols.push_back(i0);
			}
		}

		// side 1
		g4 = gadd(gadd(a, gmul(b, y4)), gmul(gadd(c, gmul(d, y4)), x4));
		//I4 = idealmul(L4, idealhnf0(L4, gsubst(g4, gvar(g4), gel(L4, 2)), NULL), J4);
		I4 = idealhnf(L4, g4);
		I4M = idealfactor(L4, I4);
		dims = matsize(I4M);
		nv = itos(gel(dims,1));
		for (int j = 0; j < nv; ++j) {
			Ip = gcoeff(I4M, j+1, 1);
			int v = itos(gcoeff(I4M, j+1, 2));
			pj = gel(Ip, 1);
			idj = gel(Ip, 2);
			fj = gel(Ip, 4);
			sj = gen_1;
			char* pjstr = GENtostr(pj);
			char* idjstr = GENtostr(idj);
			char* fjstr = GENtostr(fj);
			string ideal = "[[" + string(pjstr) + ", " + string(idjstr) + "], "
				+ string(fjstr) + ", 1]";
			free(fjstr); free(idjstr); free(pjstr);
			auto lb = lower_bound(list->begin(), list->end(), make_pair(ideal, 0),
        		[](pair<string,int> a, pair<string,int> b) { return a.first < b.first; });
    		if (lb->first == ideal) {
				int i0 = lb->second;
				idealcols.push_back(i0);
			}
		}

		// construct renumbered relation
		sort(idealcols.begin(), idealcols.end());
		rel = ABhex + ":";
		for (int j = 0; j < idealcols.size(); j++) {
			rel += long2hex(idealcols[j]);
			if (j < idealcols.size()-1) rel += ",";
		}
		data->rnrels.push_back(rel);
		
		// advance to next relation for this thread
		row++;
	}
	avma = ltop;

	// finish
	pari_thread_close();
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
