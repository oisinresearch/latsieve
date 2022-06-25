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
auto npos = std::string::npos;

typedef struct struct1 {
	int s;
	uint64_t p;
	GEN ideal;
	int e;
	string log;
} dlog;

typedef struct struct2 {
	vector<string>* relsdel;
	int rowstart;
	vector<string> newrows;
	int numnewrows;
	int workid;
} thread_root_data;

string deduce_from_both_sides(string ellstr, string filename);

void *thread_root(void* context_data);

void init_val3d(mpz_t* poly, int d, GEN *f, GEN *K, GEN *x, GEN *J1, GEN *J2)
{
    // construct poly f
	*f = gen_0;
	for (int j = 0; j <= d; ++j) {
		char str[100];
		mpz_get_str(str, 10, poly[j]);
		GEN fij = gp_read_str(str);
		*f = gadd(*f, gmul(fij, gpowgs(*x, j)));
	}
    // construct number field K, xx primitive element
    if (cmpii(leading_coeff(*f), gen_1) == 0) {
        *K = nfinit0(*f, 0, 10);
    }
    else {
        GEN p3 = nfinit0(*f, 0, 10);    // prec = 10
        *K = gcopy(gel(p3, 1));
        *x = gcopy(gel(p3, 2));
    }
    // construct constant ideals
    *J1 = idealhnf0(*K, gen_1, *x);
    *J1 = idealinv(*K, *J1);
    *J2 = idealpow0(*K, *J1, gen_2, 0);
}

bool find_log(GEN Ip, GEN K, GEN J1, GEN x, vector<dlog> dloglist,
	int side, uint64_t p, uint64_t* r, uint64_t* k)
{
	GEN IpM = gcoeff(idealfactor(K, Ip), 1, 1);
	GEN Id;
	GEN IdM;

	// bisection search
	size_t numlogs = dloglist.size();
	size_t b = numlogs/2;
	bool foundp = false;
	int i = b;
	while (b > 0) {
		b = b >> 1;
		if (dloglist[i].p < p) {
			i += b;
		}
		else if (dloglist[i].p > p) {
			i -= b;
		}
		else {
			if (dloglist[i].s == side) {
				foundp = true;
				break;
			}
			else if (dloglist[i].s < side) {
				while (dloglist[i].s < side) i--;	// side 1 comes before side 0 in dlogs
			}
			else if (dloglist[i].s > side) {
				while (dloglist[i].s > side) i++;
			}
		}
	}
	bool foundlog = false;
	if (foundp) {
		// track back to first occurance of p on this side
		if (i > 0) {
			while (dloglist[i-1].s == side && dloglist[i-1].p == p) i--;
		}
		while (i < numlogs) {
			if (dloglist[i].p != p) {
				break;
			}
			/*
			if (dloglist[i].r == p) {	// projective case
				Id = idealhnf0(K, stoi(p), gsub(gmul(stoi(p), x), gen_1));
				Id = idealmul(K, Id, J1);
			}
			else {
				Id = idealhnf0(K, stoi(p), gsub(x, stoi(dloglist[i].r)));
				Id = idealmul(K, Id, J1);
			}
			IdM = gcoeff(idealfactor(K, Id), 1, 1);
			if (cmp_universal(IdM, IpM) == 0) {	// found a match
				foundlog = true;
				*r = dloglist[i].r;
				*k = i;
				break;
			}*/
			i++;
		}
	}
	return foundlog;
}

vector<dlog> dlogs;
vector< pair<uint64_t,int> > ideals;
vector< pair<int,string> > sm;

int main(int argc, char** argv)
{
	if (argc != 12) {
		cout << "Usage:  ./descend_tnfs polyfile isofile renumberfile knownlogs ";
		cout << "ideals relsdel num_Sch_0 num_Sch_1 numthreads ell outputfile" << endl << endl;
		return 0;
	}

	string polyfile = string(argv[1]);
	string isofile = string(argv[2]);
	string renumberfile = string(argv[3]);
	string knownlogsfile = string(argv[4]);
	string idealsfile = string(argv[5]);
	string relsdelfile = string(argv[6]);
	int nsm0 = atoi(argv[7]);
	int nsm1 = atoi(argv[8]);
	int nt = atoi(argv[9]);
	string ellstr = argv[10];
	string dlogfile = string(argv[11]);
	
	string loadnfs3dstr = "\\r ~/nfs3d.gp";
	string loadpolystr = "loadpolytnfs(\"" + polyfile + "\");";

	// initialize pari library
	pari_init(2000000000,65536);

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
	GEN ideals0 = readstr(renumberfile.c_str());
	long numideals = glength(ideals0);
	GEN ideals = cgetg(numideals+1, t_VEC);
	GEN v0 = gel(ideals0, 1);
	string row1str = "[[0,[0]~],0]";
	gel(ideals, 1) = geval(strtoGENstr(row1str.c_str()));
	GEN v1, idsi, pi, idi, si, idealrow;
	int mark = 1024;
	cout << "Reading ideals..." << endl;
	for (long i = 2; i <= numideals; ++i) {
		v1 = gel(ideals0, i);
		idsi = geval(gconcat1(vecslice0(v1, 3, glength(v1))));
		pi = gcopy(gel(gel(idsi, 1), 1));
		idi = gcopy(gel(gel(idsi, 1), 2));
		si = gcopy(gel(idsi, 2));
		idealrow = cgetg(4, t_VEC);
		gel(idealrow, 1) = pi;
		gel(idealrow, 2) = idi;
		gel(idealrow, 3) = si;
		gel(ideals, i) = idealrow;
		if (i % mark == 0) {
			cout << string(i) << " ideals read..." << endl;
			mark *= 2;
		}
	}
	cout << string(numideals) << " ideals read." << endl;

	// read known logs
	GEN knownlogs;
	GEN mix = readstr(knownlogsfile.c_str());
	long numknownlogs = glength(mix);
	GEN v1, idsi, pi, idi, si, idealrow;
	cout << "Reading known logs..." << endl;
	mark = 1024;
	for (long i = 1; i <= numknownlogs; ++i) {
		v1 = gel(mix, i);
		logi = geval(v1);
		gel(knownlogs, i) = logi;
		if (i % mark == 0) {
			cout << string(i) << " known logs read..." << endl;
			mark *= 2;
		}
	}
	cout << string(numknownlogs) << " known logs read." << endl;

	// read ideals
	cout << "Reading ideals file..." << endl;
	int numJ = 1;	// watch out - hard coded number of J ideals to 1
	int col2id = new int[numknownlogs - nsm0 - nsm1];
	string line;
	char linebuffer[100];
	ifstream file(idealsfile);
	getline(file, line);	// first line contains number of (col,idnum) pairs
	mark = 1024;
	int t = 0;
	while (getline(file, line)) {
		int col = atoi(line.substr(0, line.find_first_of(" ")));
		long idnum = strtoll(line.substr(line.find_first_of(" ")+1));
		col2id[col] = idnum - numJ;
		t++;
		if (t % mark == 0) {
			cout << string(t) << " ideal column lookups read..." << endl;
			mark *= 2;
		}
	}
	int numlookups = t;
	cout << string(numlookups) << " ideal column lookups read." << endl;

	// write known logs to fblogs array
	string* fblogs = new string[numideals]();
	for (int i = 0; i < numknownlogs; i++) {
		fblogs[col2id[i]] = string(GENtostr(gel(knownlogs, i+1)));
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
	GEN eps0 = gsubgs(gpowgs(ell, 2), 1);  // N.B. CRITICAL CHOICE OF EXPONENT!
	GEN eps1 = gsubgs(gpowgs(ell, 1), 1);  // N.B. CRITICAL CHOICE OF EXPONENT!
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
		data[i].newrows = new vector<string>[nt];
		data[i].numnewrows = new int[nt];
	}

	// read relsdel
	ifstream file(relsdelfile);
	getline(file, line);	// first line contains number of relations
	mark = 1024;
	t = 0;
	while (getline(file, line)) {
		data.relsdel.push_back(line);
		t++;
		if (t % mark == 0) {
			cout << string(t) << " relations read..." << endl;
			mark *= 2;
		}
	}
	int numrels = t;
	cout << string(numrels) << " relations read." << endl;

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

	// main loop
	struct pari_thread* pth = new pari_thread[nt];
	pthread_t* th = new pthread_t[nt];
	
	int pass = 0;
	int newlogs = 0;
	stringstream stream;
	do {
		// reset newlogs each pass
		newlogs = 0;
		for (int i = 0; i < nt; i++) {
			data[i].newrows.clear();
			date[i].numnewrows = 0;
		}
		
		// create nt threads of work
		for (int i = 0; i < nt; i++) {
			datai[i].rowstart = i;
			pari_thread_alloc(&pth[i], 4000000000, args);
			pthread_create(&th[i], NULL, &thread_root, (void*)&data[i]);
		}

		// join nt threads
		for (int i = 0; i < nt; i++) {
			pthread_join(th[i], (void*)&data[i]); // this is blocking
		}

		// update factor base logs with new logs
		for (int i = 0; i < nt; i++) {
			for (int j = 0; j < data[i].numnewrows; j++) {
				string newrow = data[i].newrows[j];
				stream.str("");
				string colstr = string(atoi(newrow.substr(0, newrow.find(" "))));
				stream << colstr;
				int col;
				stream >> hex >> col;
				string log = newrow.substr(newrow.find(" ")+1);
				fblogs[col] = log;
			}
			newlogs += data[i].numnewrows;
		}
		
		// update screen
		cout << "pass " << pass << ": " << newlogs << " new logs found" << endl;

		pass++;
	}
	while (newlogs != 0);

	// write output file, ideals linked to known logs
	FILE* out;
	out = fopen(dlogfile, "w+");
	for (int i = 0; i < numideals; i++) {
		if (!fblogs[i].empty()) {
			string line = string(i) + " " + fblogs[i] + " " 
				+ string(GENtostr(gel(ideals, i+1)));
			fprintf(out, "%s\n", line);
			t++;
		}
		if (t % mark == 0) {
			cout << string(t) << " ideal logs written..." << endl;
			mark *= 2;
		}
	}
	cout << string(t) << " ideal logs written." << endl;
	
	// write Schirokauer maps
	for (int i = 0; i < nsm0; i++) {
		fprintf(out, "%s\n", "sm " + fblogs[numideals+i] + " " + string(i) + " 0");
	}
	for (int i = 0; i < nsm1; i++) {
		fprintf(out, "%s\n", "sm " + fblogs[numideals+nsm0+i] + " " + string(nsm0+i) + " 1");
	}
	cout << string(nsm0 + nsm1) << " Schirokauer map logs written." << endl;
	fclose(out);

	// tidy up
	for (int i = 0; i < nt; i++) {
		delete[] data[i].numnewrows;
		delete[] data[i].newrows;
	}
	delete[] thread_root_data;
	delete[] col2id;
	for (int i = 0; i < nt; i++) pari_thread_free(&pth[i]);
	delete[] pth;
	delete[] th;
	delete[] result;
	pari_close();

	return 0;
}

void *thread_root(void* context_data)
{
	thread_root_data* data = (thread_root_data*)context_data;

	// start pari thread
	GEN args = pari_thread_start(&data->pth[workid]);

	pari_sp ltop = avma;
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

	int row = data->startrow;
	while (row < data->numrels) {
		string rel = data->relsdel[row];
	
		GEN a = stoi(inta);
		GEN b = stoi(intb);
		GEN c = stoi(intc);
		GEN d = stoi(intd);

		GEN I3, I4, I3M, I4M;	
		GEN g3, g4;
		GEN Ip, pj, idj, sj;
		GEN item0 = cgetg(4, t_VEC);
		long i0; GEN logi0;
		cout << inta << " " << intb << " " << intc << " " << intd << endl;
		cout << "2 J 1 " << GENtostr(log0) << endl;
		string side0logstr = "";
		string side1logstr = "";
		if (side == 0 || side == 2) {
			cout << "side 0:" << endl;
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
				sj = gen_0;
				gel(item0, 1) = gcopy(pj);
				gel(item0, 2) = gcopy(idj);
				gel(item0, 3) = gcopy(sj);
				i0 = gtos(stoi(vecsearch(ideals, item0, NULL)));
				if (i0 != 0) {
					logi0 = gcopy(gel(logs, i0));			 
					cout << GENtostr(sj) << " " << GENtostr(pj) << " " << GENtostr(idj) << " "
						<< v << " " << GENtostr(logi0) << endl;
					side0logstr += ("+"+to_string(v)+"*"+GENtostr(logi0));
				}
				else {
					// attempt to deduce log from file named "deduce_{pj}.txt" for prime {pj}
					string deduce_pj = "deduce_" + string(GENtostr(pj)) + ".txt";
					string strlog = deduce_from_both_sides(ellstr, deduce_pj);
					cout << GENtostr(sj) << " " << GENtostr(pj) << " " << GENtostr(idj) << " "
						<< v << " " << strlog << endl;
					side0logstr += ("+"+to_string(v)+"*"+strlog);
				}
			}
			// print Schirokauer maps for side 0
			GEN A1 = liftall(gmul(m1, gmodulo(gadd(gadd(a, gmul(b, y3)), gmul(gadd(c, gmul(d, y3)), x3)), f)));
			GEN feps = gdiv(liftall(gsubgs(gpow(gmodulo(gmodulo(A1, f), gsqr(ell)), eps0, prec), 1)), ell);
			for (int i = 0; i < nsm0; i++) {
				string sm_exp_i = string(GENtostr(polcoef(feps, i+1, -1)));
				cout << "0 sm_exp " << i+1 << " " << sm_exp_i << endl;
				v1 = strsplit(gel(mix0, numlogs+1+i), strtoGENstr(" "));
				string sm0i = string(GENtostr(geval(gel(v1, 2))));
				side0logstr += ("+"+ sm_exp_i + "*" + sm0i);
			}
		}
		if (side == 1 || side == 2) {
			cout << "side 1:" << endl;
			g4 = gadd(gadd(a, gmul(b, y4)), gmul(gadd(c, gmul(d, y4)), x4));
			I4 = idealmul(L4, idealhnf0(L4, gsubst(g4, gvar(g4), gel(L4, 2)), NULL), J4);
			I4M = idealfactor(L4, I4);
			dims = matsize(I4M);
			int nv = itos(gel(dims,1));
			for (int j = 0; j < nv; ++j) {
				Ip = gcoeff(I4M, j+1, 1);
				int v = itos(gcoeff(I4M, j+1, 2));
				pj = gel(Ip, 1);
				idj = gel(Ip, 2);
				sj = gen_1;
				gel(item0, 1) = gcopy(pj);
				gel(item0, 2) = gcopy(idj);
				gel(item0, 3) = gcopy(sj);
				i0 = gtos(stoi(vecsearch(ideals, item0, NULL)));
				if (i0 != 0) {
					logi0 = gcopy(gel(logs, i0));			 
					cout << GENtostr(sj) << " " << GENtostr(pj) << " " << GENtostr(idj) << " "
						<< v << " " << GENtostr(logi0) << endl;
					side1logstr += ("+"+to_string(v)+"*"+GENtostr(logi0));
				}
				else {
					// attempt to deduce log from file named "deduce_{pj}.txt" for prime {pj}
					string deduce_pj = "deduce_" + string(GENtostr(pj)) + ".txt";
					string strlog = deduce_from_both_sides(ellstr, deduce_pj); 
					cout << GENtostr(sj) << " " << GENtostr(pj) << " " << GENtostr(idj) << " "
						<< v << " " << strlog << endl;
					side1logstr += ("+"+to_string(v)+"*"+strlog);
				}
			}
			// print Schirokauer maps for side 1
			GEN A2 = liftall(gmul(m2, gmodulo(gadd(gadd(a, gmul(b, y4)), gmul(gadd(c, gmul(d, y4)), x4)), g)));
			GEN geps = gdiv(liftall(gsubgs(gpow(gmodulo(gmodulo(A2, g), gsqr(ell)), eps1, prec), 1)), ell);
			for (int i = 0; i < nsm1; i++) {
				string sm_exp_i = string(GENtostr(polcoef(geps, i+1, -1)));
				cout << "1 sm_exp " << i+1 << " " << sm_exp_i << endl;
				v1 = strsplit(gel(mix0, numlogs+numsch0+1+i), strtoGENstr(" "));
				string sm1i = string(GENtostr(geval(gel(v1, 2))));
				side1logstr += ("+"+ sm_exp_i + "*" + sm1i);
			}
		}
		cout << "Schirokauer Maps:" << endl;
		for (int i = 0; i < numsch0 + numsch1; i++) {
			v1 = strsplit(gel(mix0, numlogs+1+i), strtoGENstr(" "));
			cout << GENtostr(geval(gel(v1, 4))) << " SM " << i+1 << " "
				<< GENtostr(geval(gel(v1, 2))) << endl;
		}
		if (side == 0) cout << side0logstr;
		if (side == 1) cout << side1logstr;

		row += data->nt;
	}

	avma = ltop;

	// finish
	pari_thread_close();
	return NULL;
}


string deduce_from_both_sides(string ellstr, string filename)
{
	vector<pair<int,string> > sm_exp;
	vector<dlog> dlogs;
	vector<pair<int,string> > sm;

    string separator0 = " ";
    string separator1 = ":";
    string separator2 = ",";
	string separator3 = "[";
	string separator5 = "~";
	string line;
	string tlogstr = "";

	ifstream file(filename);
	if (file.good()) {
		string pfilestr = filename.substr(7, filename.find(".") - 7);
		int64_t pfile = strtoull(pfilestr.c_str(), NULL, 10);

		// read a b c d
		getline(file, line);
		int64_t a = strtoll(line.substr(0, line.find(separator0)).c_str(), NULL, 10);
		line.erase(0, line.find(separator0) + 1);
		int64_t b = strtoll(line.substr(0, line.find(separator0)).c_str(), NULL, 10);
		line.erase(0, line.find(separator0) + 1);
		int64_t c = strtoll(line.substr(0, line.find(separator0)).c_str(), NULL, 10);
		line.erase(0, line.find(separator0) + 1);
		int64_t d = strtoll(line.substr(0, line.find(separator0)).c_str(), NULL, 10);
		line.erase(0, line.find(separator0) + 1);

		mpz_t logJ; mpz_init(logJ);
		getline(file, line);
		line.erase(0, line.find(separator0) + 1);
		line.erase(0, line.find(separator0) + 1);
		line.erase(0, line.find(separator0) + 1);
		mpz_set_str(logJ, line.c_str(), 10);

		int tside = 0;

		// read "side 0:"
		getline(file, line);

		// read side 0 ideal logs/Schirokauer exponents
		while (getline(file, line)) {
			if (line.find("side 1:", 0) == 0) break;

			int side = atoi(line.substr(0, line.find(separator0)).c_str());
			line.erase(0, line.find(separator0) + 1);

			string pstr = line.substr(0, line.find(separator0));
			line.erase(0, line.find(separator0) + 1);

			GEN id = gp_read_str(line.substr(0, line.find(separator5) + 1).c_str());
			line.erase(0, line.find(separator5) + 2);

			int e = atoi(line.substr(0, line.find(separator0)).c_str());
			if (line.find(separator0) == npos) line = "";
			else line.erase(0, line.find(separator0) + 1);

			string val = line;

			if (pstr == "sm_exp") {
				sm_exp.push_back(make_pair(side, val));
			}
			else {
				uint64_t p = strtoull(pstr.c_str(), NULL, 10);
				if (val == "" && p != pfile) {
					val = deduce_from_both_sides(ellstr, "deduce_" + to_string(p) + ".txt");
					//cout << "\tvlog(" << p << ") = " << val << " on side " << 0 << endl;
				}
				dlogs.push_back((dlog){ side, p, gcopy(id), e, val });
			}
		}

		// read side 1 ideal logs/Schirokauer exponents
		while (getline(file, line)) {
			if (line.find("Schirokauer Maps:", 0) == 0) break;

			int side = atoi(line.substr(0, line.find(separator0)).c_str());
			line.erase(0, line.find(separator0) + 1);

			string pstr = line.substr(0, line.find(separator0));
			line.erase(0, line.find(separator0) + 1);

			GEN id = gp_read_str(line.substr(0, line.find(separator5) + 1).c_str());
			line.erase(0, line.find(separator5) + 2);

			int e = atoi(line.substr(0, line.find(separator0)).c_str());
			if (line.find(separator0) == npos) line = "";
			else line.erase(0, line.find(separator0) + 1);

			string val = line;

			if (pstr == "sm_exp") {
				sm_exp.push_back(make_pair(side, val));
			}
			else {
				uint64_t p = strtoull(pstr.c_str(), NULL, 10);
				if (val == "" && p != pfile) {
					val = deduce_from_both_sides(ellstr, "deduce_" + to_string(p) + ".txt");
					//cout << "\tvlog(" << p << ") = " << val << " on side " << 1 << endl;
				}
				dlogs.push_back((dlog){ side, p, gcopy(id), e, val });
			}
		}

		// read Schirokauer map logs
		while (getline(file, line)) {
			line.erase(0, line.find(separator0) + 1);

			int side = atoi(line.substr(0, line.find(separator0)).c_str());
			line.erase(0, line.find(separator0) + 1);

			uint64_t r = strtoll(line.substr(0, line.find(separator0)).c_str(), NULL, 10);
			line.erase(0, line.find(separator0) + 1);
			
			string val = line;

			sm.push_back(make_pair(side, val));
		}

		// deduce target
		int64_t q = 0;
		mpz_t tlog; mpz_init(tlog);
		mpz_t log; mpz_init(log);
		mpz_t ell; mpz_init(ell);
		mpz_t exp; mpz_init(exp);
		mpz_set_str(ell, ellstr.c_str(), 10);
		int ndlogs = dlogs.size();
		mpz_set_ui(tlog, 0);
		// ideal logs first
		for (int i = 0; i < ndlogs; i++) {
			int side = dlogs[i].s;
			int64_t p = dlogs[i].p;
			int e = dlogs[i].e;
			string val = dlogs[i].log;
			if (val != "") {
				//	cout << " + " << e << "*" << val;
				mpz_set_str(log, val.c_str(), 10);
				mpz_mul_si(log, log, e);
				mpz_add(tlog, tlog, log);
				mpz_mod(tlog, tlog, ell);	// reduce mod ell
			}
			else {
				q = p;
			}
		}
		// then Schirokauer maps
		int nsm = sm.size();
		for (int i = 0; i < nsm; i++) {
			int side = sm[i].first;
			string eval = sm_exp[i].second;
			string val = sm[i].second;
			//	cout << " + " << eval << "*" << val;
			mpz_set_str(exp, eval.c_str(), 10);
			mpz_set_str(log, val.c_str(), 10);
			mpz_mul(log, log, exp);
			mpz_add(tlog, tlog, log);
			mpz_mod(tlog, tlog, ell);	// reduce mod ell		
		}
		//if (target_side == 1) {
		mpz_add(tlog, tlog, logJ);	// account for J ideal (warning - J should be for both sides)
		//	cout << " + " << mpz_get_str(NULL, 10, logJ);
		//}

		//cout << endl << "\ttarget side = " << target_side << endl;

		//if (target_side == 0) {
		mpz_sub(tlog, ell, tlog);	// move known logs to other side
		//}

		mpz_mod(tlog, tlog, ell);

		tlogstr = mpz_get_str(NULL, 10, tlog);

		//	cout << endl;
		//cout << "vlog(" << q << ") = " << mpz_get_str(NULL, 10, tlog) << endl;

		mpz_clear(exp);
		mpz_clear(ell);
		mpz_clear(log);
		mpz_clear(tlog);
	}

	return tlogstr;
}

