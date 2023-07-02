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
	int nt;
	struct pari_thread pth;
	string ellstr;
	vector<string>* rels;
	vector<string>* fblogs;
	vector<string>* smlogs;
	int numrels;
	int startrow;
	int numnewlogs;
	vector<string> newlogs;
} thread_root_data;

typedef struct struct3 {
	void *RETURNVALUE;
	vector<string> rel;
	vector<string> log;
	int count;
} newrel;

string deduce_from_both_sides(int64_t a, int64_t b, int64_t c, int64_t d, int nsm0, int nsm1,
	thread_root_data* data, vector<pair<long, int>> knownlogs, vector<string> sm_exp);
void *thread_root(void* context_data);
long hex2long(string str1);
string long2hex(long n);

int main(int argc, char** argv)
{
	if (argc != 13) {
		cout << "Usage:  ./reconstructlog-dl-tnfs polyfile isofile renumberfile knownlogs ";
		cout << "ideals purgedrels relsdel num_Sch_0 num_Sch_1 numthreads ell outlogs";
		cout << endl << endl; // origrels "
//			 << "outmat outvec outideals" << endl << endl;
		return 0;
	}

	string polyfile = string(argv[1]);
	string isofile = string(argv[2]);
	string renumberfile = string(argv[3]);
	string knownlogsfile = string(argv[4]);
	string idealsfile = string(argv[5]);
	string purgedrelsfile = string(argv[6]);
	string relsdelfile = string(argv[7]);
	int nsm0 = atoi(argv[8]);
	int nsm1 = atoi(argv[9]);
	int nt = atoi(argv[10]);
	string ellstr = argv[11];
	
	string loadnfs3dstr = "\\r ~/nfs3d.gp";
	string loadpolystr = "loadpolytnfs(\"" + polyfile + "\");";

	// initialize pari library
	pari_init(20000000000, 65536);

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
	GEN ideals = readstr(renumberfile.c_str());
	long numideals = glength(ideals);
	cout << to_string(numideals) << " ideals read." << endl;

	// read known logs
	GEN logi; GEN v1;
	GEN mix = readstr(knownlogsfile.c_str());
	long numknownlogs = glength(mix);
	GEN knownlogs = cgetg(numknownlogs+1, t_VEC);
	cout << "Reading known logs..." << endl;
	int mark = 1024;
	for (long i = 1; i <= numknownlogs; ++i) {
		v1 = gel(mix, i);
		logi = geval(v1);
		gel(knownlogs, i) = gcopy(logi);
		if (i % mark == 0) {
			cout << to_string(i) << " known logs read..." << endl;
			mark *= 2;
		}
	}
	cout << to_string(numknownlogs) << " known logs read." << endl << endl;

	// read ideals
	cout << "Reading ideals file..." << endl;
	int numJ = 1;	// watch out - hard coded number of J ideals to 1
	int* col2id = new int[numknownlogs - nsm0 - nsm1]();
	string line;
	char linebuffer[100];
	ifstream file1(idealsfile);
	getline(file1, line);	// first line contains number of (col,idnum) pairs
	mark = 1024;
	int t = 0;
	while (getline(file1, line)) {
		int col = atoi(line.substr(0, line.find_first_of(" ")).c_str());
		string str1 = line.substr(line.find_first_of(" ")+1);
		long idnum = hex2long(str1);
		col2id[col] = idnum;
		t++;
		if (t % mark == 0) {
			cout << to_string(t) << " ideal column lookups read..." << endl;
			mark *= 2;
		}
	}
	int numlookups = t;
	cout << to_string(numlookups) << " ideal column lookups read." << endl << endl;
	file1.close();

	// write known logs to fblogs array
	vector<string> fblogs;
	for (int i = 0; i < numideals; i++) {
		fblogs.push_back("");
	}
	for (int i = 0; i < numknownlogs-nsm0-nsm1; i++) {
		fblogs[col2id[i]] = string(GENtostr(gel(knownlogs, i+1)));
	}
	// also record Schirokauer map logs
	vector<string> smlogs;
	for (int i = 0; i < nsm0+nsm1; i++) {
		int j = numknownlogs-nsm0-nsm1+i;
		smlogs.push_back(string(GENtostr(gel(knownlogs, j+1))));
	}

	// set up number fields and J-ideals
	long prec = 5;
	GEN x = pol_x(0);
	GEN L3 = nfinit0(f, 3, prec);
	GEN L4 = nfinit0(g, 3, prec);
	GEN J3 = idealinv(L3, idealhnf0(L3, gen_1, gsubst(x3, gvar(x3), gel(L3, 2))));
	GEN J4 = idealinv(L4, idealhnf0(L4, gen_1, gsubst(x4, gvar(x4), gel(L4, 2))));

	// set up Schirokauer map data
	GEN ff = factor(gaddgs(gsqr(p), 1));
	GEN dims;
	dims = matsize(ff);
	int nv = itos(gel(dims,1));
	GEN ell = gcoeff(ff, nv, 1);
	GEN eps0 = gsubgs(gpowgs(ell, 4), 1);	//2), 1);  // N.B. CRITICAL CHOICE OF EXPONENT!  THIS WAS MAJOR HASSLE! SCUPPERED RECONSTRUCTION STEP!
	GEN eps1 = gsubgs(gpowgs(ell, 4), 1);	//1), 1);  // N.B. CRITICAL CHOICE OF EXPONENT!
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

	// begin threading code
	thread_root_data* data = new thread_root_data[nt];
	for (int i = 0; i < nt; i++) {
		data[i].nt = nt;
		data[i].ellstr = ellstr;
		data[i].fblogs = &fblogs;
		data[i].smlogs = &smlogs;
	}

	// read purgedrels
	cout << " reading purged rels file..." << endl;
	vector<string> rels;
	ifstream file2(purgedrelsfile);
	//getline(file2, line);	// first line contains number of relations
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
	for (int i = 0; i < nt; i++) {
		data[i].rels = &rels;
	}
	file2.close();

	// read relsdel
	cout << " reading relsdel file..." << endl;
	//vector<string> relsdel;
	ifstream file3(relsdelfile);
	getline(file3, line);	// first line contains number of relations
	mark = 1024;
	t = 0;
	while (getline(file3, line)) {
		rels.push_back(line);
		t++;
		if (t % mark == 0) {
			cout << to_string(t) << " relations read..." << endl;
			mark *= 2;
		}
	}
	cout << to_string(t) << " relations read." << endl << endl;
	int numrels = rels.size();
	for (int i = 0; i < nt; i++) {
		data[i].numrels = numrels - nsm0 - nsm1;
	}
	file3.close();

	GEN args = cgetg(13, t_VEC);
	gel(args, 0) = f; 
	gel(args, 1) = g; 
	gel(args, 2) = ell;
	gel(args, 3) = eps0; 
	gel(args, 4) = eps1; 
	gel(args, 5) = y3; 
	gel(args, 6) = x3; 
	gel(args, 7) = y4; 
	gel(args, 8) = x4; 
	gel(args, 9) = m1; 
	gel(args, 10) = m2; 
	gel(args, 11) = stoi(nsm0);
	gel(args, 12) = stoi(nsm1);

	pthread_t* th = new pthread_t[nt];
	newrel* newrels = new newrel[nt];
	for (int i = 0; i < nt; i++) {
		pari_thread_alloc(&data[i].pth, 4000000000, args);
		newrels[i].count = 0;
	}

	int pass = 0;
	int newlogs = 0;
	// main loop
	do {
		cout << "starting pass #" << pass << endl;
		// reset newlogs each pass
		newlogs = 0;
		for (int i = 0; i < nt; i++) {
			data[i].newlogs.clear();
			data[i].numnewlogs = 0;
		}
		
		// create nt threads of work
		for (int i = 0; i < nt; i++) {
			data[i].startrow = i;
			pthread_create(&th[i], NULL, &thread_root, (void*)&data[i]);
		}

		// join nt threads
		for (int i = 0; i < nt; i++) {
			pthread_join(th[i], (void**)&data[i]); // this is blocking
		}

		// update factor base logs with new logs
		for (int i = 0; i < nt; i++) {
			for (int j = 0; j < data[i].numnewlogs; j++) {
				string newlog = data[i].newlogs[j];
				long jcol = strtoll(newlog.substr(0, newlog.find(" ")).c_str(), NULL, 10);
				string jlogrel = newlog.substr(newlog.find(" ")+1);
				string jlog = jlogrel.substr(0, jlogrel.find(" "));
				string jrel = jlogrel.substr(jlogrel.find(" ")+1);
				newrels[i].rel.push_back(jrel);
				newrels[i].log.push_back(jlog);
				newrels[i].count++;
				fblogs[jcol] = jlog;
			}
			newlogs += data[i].numnewlogs;
		}
		
		// update screen
		cout << "pass #" << pass << ": " << newlogs << " new logs found" << endl;

		pass++;
	}
	while (newlogs != 0);

	// write output file, ideals linked to known logs
	FILE* out;
	out = fopen(argv[12], "w+");
	t = 0;
	cout << endl << "Writing logs to output file " << argv[11] << "..." << endl;
	for (int i = 0; i < numideals; i++) {
		if (!fblogs[i].empty()) {
			string line = long2hex(i) + " " + fblogs[i] + " " 
				+ string(GENtostr(geval(gel(ideals, i+1))));
			if (i < 10) cout << line << endl;
			fprintf(out, "%s\n", line.c_str());
			t++;
		}
		if (t % mark == 0) {
			cout << to_string(t) << " ideal logs written..." << endl;
			mark *= 2;
		}
	}
	int numidealhits = t;
	cout << to_string(t) << " ideal logs written." << endl;
	
	// write Schirokauer maps
	for (int i = 0; i < nsm0; i++) {
		string line = "sm " + smlogs[i] + " " + to_string(i) + " 0";
		fprintf(out, "%s\n", line.c_str());
	}
	for (int i = 0; i < nsm1; i++) {
		string line = "sm " + smlogs[nsm0+i] + " " + to_string(nsm0+i) + " 1";
		fprintf(out, "%s\n", line.c_str());
	}
	cout << to_string(nsm0 + nsm1) << " Schirokauer map logs written." << endl;
	fclose(out);

/*
	// write output ideals linked to known logs
	FILE* outideals;
	outideals = fopen(argv[15], "w+");
	cout << endl << "Writing ideal columns to output file " << argv[15] << "..." << endl;
	string line1 = "# " + to_string(numidealhits);
	fprintf(outideals, "%s\n", line1.c_str());
	t = 0;
	for (int i = 0; i < numideals; i++) {
		if (!fblogs[i].empty()) {
			string line = to_string(t) + " " + long2hex(i);
			fprintf(outideals, "%s\n", line.c_str());
			t++;
		}
		if (t % mark == 0) {
			cout << to_string(t) << " ideal columns written..." << endl;
			mark *= 2;
		}
	}
	cout << to_string(t) << " ideal columns written." << endl;
	fclose(outideals);

	// write output relations file
	FILE* outmat;
	outmat = fopen(argv[13], "w+");
	t = 0;
	cout << endl << "Writing output relations from " << argv[12] << "..." << endl;
	// first copy original relations
	ifstream rels(argv[12]);
	t = 0;
	mark = 1024;
    while (getline(rels, line)) {
		fprintf(outmat, "%s\n", line.c_str());
		t++;
		if (t % mark == 0) {
			mark *= 2;
		}
	}
	// now write 16 threads worth of discovered relations
	t = 0;
	mark = 1024;
	for (int i = 0; i < nt; i++) {
		for (int j = 0; j < newrels[i].count; j++) {
			fprintf(outmat, "%s\n", newrels[i].rel[j].c_str());
			t++;
			if (t > mark) {
				mark *= 2;
			}
		}
	}
	fclose(outmat);

	// write output vector file
	FILE* outvec;
	outvec = fopen(argv[14], "w+");
	cout << endl << "Writing output nullspace vector " << argv[13] << "..." << endl;
	t = 0;
	for (int i = 0; i < numideals; i++) {
		if (!fblogs[i].empty()) {
			string line = fblogs[i];
			fprintf(outvec, "%s\n", line.c_str());
			t++;
		}
		if (t % mark == 0) {
			cout << to_string(t) << " vector logs written..." << endl;
			mark *= 2;
		}
	}
	cout << to_string(t) << " vector logs written." << endl;
	fclose(outvec);

	cout << endl << "Finished!" << endl;
*/
	// tidy up
	for (int i = 0; i < nt; i++) pari_thread_free(&data[i].pth);
	delete[] col2id;
	delete[] th;
	delete[] newrels;
	delete[] data;
	pari_close();

	return 0;
}

void *thread_root(void* context_data)
{
	thread_root_data* data = (thread_root_data*)context_data;

	// start pari thread
	GEN args = pari_thread_start(&data->pth);

	long prec = 5;
	GEN f = gel(args, 0); 
	GEN g = gel(args, 1); 
	GEN ell = gel(args, 2);
	GEN eps0 = gel(args, 3); 
	GEN eps1 = gel(args, 4); 
	GEN y3 = gel(args, 5); 
	GEN x3 = gel(args, 6); 
	GEN y4 = gel(args, 7); 
	GEN x4 = gel(args, 8); 
	GEN m1 = gel(args, 9); 
	GEN m2 = gel(args, 10);
	int nsm0 = itos(gel(args, 11));
	int nsm1 = itos(gel(args, 12));

	data->numnewlogs = 0;
	int mark = 1024;
	int row = data->startrow;
	while (row < data->numrels) {
		pari_sp ltop = avma;
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

		//cout << a1 << " " << b1 << " " << c1 << " " << d1 << endl;
		//cout << "2 J 1 " << GENtostr(log0) << endl;
		string side0logstr = "";
		string side1logstr = "";
		// split string of ideals
		string idealcsv = rel.substr(rel.find(":")+1);
		vector<long> idealcols;
		idealcols.clear();
		stringstream ss;
		ss.str("");
		ss << idealcsv;
		while (ss.good()) {
			string id;
			getline(ss, id, ',');
			idealcols.push_back(hex2long(id));
		}
		sort(idealcols.begin(), idealcols.end());

		int numunk = 0;
		long unklogj = 0;
		int nv = idealcols.size();
		vector<pair<long, int>> knownjlogs;
		knownjlogs.clear();
		int e = 0;
		long jlast = -1;	// no ideal has this column number
		for (int j = 0; j <= nv; ++j) {
			long jcol = -1;
			if (j <= nv) {
				jcol = idealcols[j];
			}
			if (jcol != jlast) {
				if (jlast >= 0) {
					if ((*data->fblogs)[jlast].empty()) {
						numunk++;
						unklogj = jlast;	// found an unknown log, hopefully just 1
					}
					else {
						knownjlogs.push_back(make_pair(jlast, e));
					}
				}
				jlast = jcol;
				e = 1;
			}
			else {
				e++;
			}
		}

		if (numunk == 1) {
			// first compute Schirokauer maps for this (a,b,c,d)
			vector<string> sm_exp;
			sm_exp.clear();
			// side 0
			GEN A1 = liftall(gmul(m1, gmodulo(gadd(gadd(a, gmul(b, y3)), gmul(gadd(c, gmul(d, y3)), x3)), f)));
			GEN feps = gdiv(liftall(gsubgs(gpow(gmodulo(gmodulo(A1, f), gsqr(ell)), eps0, prec), 1)), ell);
			for (int i = 0; i < nsm0; i++) {
				string sm_exp_i = string(GENtostr(polcoef(feps, i+1, -1)));
				sm_exp.push_back(sm_exp_i);
			}
			// side 1
			GEN A2 = liftall(gmul(m2, gmodulo(gadd(gadd(a, gmul(b, y4)), gmul(gadd(c, gmul(d, y4)), x4)), g)));
			GEN geps = gdiv(liftall(gsubgs(gpow(gmodulo(gmodulo(A2, g), gsqr(ell)), eps1, prec), 1)), ell);
			for (int i = 0; i < nsm1; i++) {
				string sm_exp_i = string(GENtostr(polcoef(geps, i+1, -1)));
				sm_exp.push_back(sm_exp_i);
			}
			// deduce single unknown log from the known ones
			string log = deduce_from_both_sides(a1, b1, c1, d1, nsm0, nsm1, data,
				knownjlogs, sm_exp);
			data->newlogs.push_back(to_string(unklogj) + " " + log + " " + rel);
			data->numnewlogs++;
			//cout << "	found log of ideal " << unklogj << endl;
		}

		// advance to next relation for this thread
		row += data->nt;
		if (row > mark && data->startrow == 0) {
			cout << row << " relations processed..." << endl;
			mark *= 2;
		}
		avma = ltop;
	}

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

string deduce_from_both_sides(int64_t a, int64_t b, int64_t c, int64_t d, int nsm0, int nsm1,
	thread_root_data* data, vector<pair<long, int>> knownlogs, vector<string> sm_exp)
{
	mpz_t logJ; mpz_init(logJ);
	mpz_set_str(logJ, (*data->fblogs)[0].c_str(), 10);

	// deduce target
	mpz_t tlog; mpz_init(tlog);
	mpz_t log; mpz_init(log);
	mpz_t ell; mpz_init(ell);
	mpz_t exp; mpz_init(exp);
	mpz_set_str(ell, data->ellstr.c_str(), 10);
	mpz_set_ui(tlog, 0);
	// ideal logs first
	for (int i = 0; i < knownlogs.size(); i++) {
		int jcol = knownlogs[i].first;
		int e = knownlogs[i].second;
		string val = (*data->fblogs)[jcol];
		if (val != "") {
			//	cout << " + " << e << "*" << val;
			mpz_set_str(log, val.c_str(), 10);
			mpz_mul_si(log, log, e);
			mpz_add(tlog, tlog, log);
			mpz_mod(tlog, tlog, ell);	// reduce mod ell
		}
	}
	// then Schirokauer maps
	for (int i = 0; i < nsm0+nsm1; i++) {
		string eval = sm_exp[i];
		string val = (*data->smlogs)[i];
		//	cout << " + " << eval << "*" << val;
		mpz_set_str(exp, eval.c_str(), 10);
		mpz_set_str(log, val.c_str(), 10);
		mpz_mul(log, log, exp);
		mpz_add(tlog, tlog, log);
		mpz_mod(tlog, tlog, ell);	// reduce mod ell		
	}
	
	//mpz_add(tlog, tlog, logJ);	// account for J ideal (warning - J should be for both sides)
	mpz_sub(tlog, ell, tlog);	// move known logs to other side
	mpz_mod(tlog, tlog, ell);

	string tlogstr = mpz_get_str(NULL, 10, tlog);

	//cout << "vlog(" << q << ") = " << mpz_get_str(NULL, 10, tlog) << endl;

	mpz_clear(exp);
	mpz_clear(ell);
	mpz_clear(log);
	mpz_clear(tlog);

	return tlogstr;
}

