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
auto npos = std::string::npos;

typedef struct struct1 {
	int s;
	uint64_t p;
	uint64_t r;
	string log;
} dlog;

typedef struct struct2 {
	vector<string>* lines;
	struct pari_thread pth;
	int workid;
} thread_root_arg;

typedef struct struct3 {
	uint64_t p;
	uint64_t r;
	int s;
	string a;
	int n;
} badideal;

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

bool is_bad(GEN Ip, vector<badideal> badlist)
{
	bool found = false;
	int nbad = badlist.size();
	for (int i = 0; i < nbad; i++) {
		string astr(GENtostr(gel(Ip, 2)));
		if (astr == badlist[i].a) {
			found = true;
			break;
		}
	}
	return found;
}

int find_bad(GEN Ip, vector<badideal> badlist)
{
	int k = 0;
	int nbad = badlist.size();
	for (int i = 0; i < nbad; i++) {
		string astr(GENtostr(gel(Ip, 2)));
		if (astr == badlist[i].a) {
			k = i;
			break;
		}
	}
	return k;
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
			}
			i++;
		}
	}
	return foundlog;
}

vector<badideal> badideals;
int badpmax = 0;
vector<dlog> dlogs;
vector< pair<uint64_t,int> > ideals;
vector< pair<int,string> > sm;

int main(int argc, char** argv)
{
	if (argc != 8) {
		cout << "Usage:  ./descend polyfile badidealfile dlogfile a b c side" << endl << endl;
		cout << "Output:  Lines of the form" << endl;
		cout << "side prime root exponent log" << endl << endl;
		return 0;
	}

	// read input polynomial
	mpz_t* fpoly = new mpz_t[20];	// max degree of 20.  Not the neatest
	mpz_t* gpoly = new mpz_t[20];	// max degree of 20.  Not the neatest
	for (int i = 0; i < 20; i++) {
		mpz_init(fpoly[i]);
		mpz_init(gpoly[i]);
	}
	int64_t* fi64 = new int64_t[20]();
	int64_t* gi64 = new int64_t[20]();
	int64_t* fmodp = new int64_t[20]();
	int64_t* gmodp = new int64_t[20]();
	string line;
	char linebuffer[100];
	ifstream file(argv[1]);
	getline(file, line);	// first line contains number n to factor
	// read nonlinear poly
	int degf = -1;
	while (getline(file, line) && line.substr(0,1) == "c" ) {
		line = line.substr(line.find_first_of(" ")+1);
		mpz_set_str(fpoly[++degf], line.c_str(), 10);
		fi64[degf] = mpz_get_si(fpoly[degf]);
	}
	// read other poly
	int degg = -1;
	bool read = true;
	while (read && line.substr(0,1) == "Y" ) {
		line = line.substr(line.find_first_of(" ")+1);
		mpz_set_str(gpoly[++degg], line.c_str(), 10);
		gi64[degg] = mpz_get_si(gpoly[degg]);
		read = static_cast<bool>(getline(file, line));
	}
	file.close();
	int64_t* froots = new int64_t[degf+1];
	int64_t* groots = new int64_t[degg+1];


	// initialize pari library
	pari_init(100000000,65536);
	
	// set up constants
	GEN f = pol_x(fetch_user_var("f"));
	GEN Ka = cgetg(DEFAULTPREC, typ_NF);
	GEN x = pol_x(fetch_user_var("x"));
	GEN J1a = pol_x(fetch_user_var("J1a"));
	GEN J2a = pol_x(fetch_user_var("J2a"));
	
	init_val3d(fpoly, degf, &f, &Ka, &x, &J1a, &J2a);
	
	GEN g = pol_x(fetch_user_var("g"));
	GEN Kb = cgetg(DEFAULTPREC, typ_NF);
	GEN y = pol_x(fetch_user_var("y"));
	GEN J1b = pol_x(fetch_user_var("J1b"));
	GEN J2b = pol_x(fetch_user_var("J2b"));
	
	init_val3d(gpoly, degg, &g, &Kb, &y, &J1b, &J2b);


	// read bad ideals file
    string separator0 = " ";
    string separator1 = ":";
    string separator2 = ",";
	string badidealsfile = argv[2];
	GEN I1a;
	GEN I1aM;
	GEN adims;
	GEN I1b;
	GEN I1bM;
	int nbad = 1;
	ifstream input(badidealsfile);
	while (getline(input, line)) {
		uint64_t p = atoi(line.substr(0, line.find(separator2)).c_str());
		line.erase(0, line.find(separator2) + 1);
		uint64_t r = atoi(line.substr(0, line.find(separator1)).c_str());
		line.erase(0, line.find(separator1) + 1);
		int s = atoi(line.substr(0, line.find(separator1)).c_str());
		line.erase(0, line.find(separator1) + 1);
		//int e = atoi(line.c_str());

		if (s == 0) {
			if (p == r) {
				I1a = idealhnf0(Ka, stoi(p), gsub(gmul(stoi(p), x), gen_1));
			}
			else {
				I1a = idealhnf0(Ka, stoi(p), gsub(x, stoi(r)));
			}
			I1a = idealmul(Ka, I1a, J1a);
			I1aM = idealfactor(Ka, I1a);
			adims = matsize(I1aM);
			int nv = itos(gel(adims, 1));
			for (int j = 0; j < nv; j++) {
				string a(GENtostr(gel(gcoeff(I1aM, j+1, 1),2)));
				badideals.push_back((badideal){ p, r, s, a, nbad++ });
			}
		}
		else {
			if (p == r) {
				I1b = idealhnf0(Kb, stoi(p), gsub(gmul(stoi(p), y), gen_1));
			}
			else {
				I1b = idealhnf0(Kb, stoi(p), gsub(y, stoi(r)));
			}
			I1b = idealmul(Kb, I1b, J1b);
			I1bM = idealfactor(Kb, I1b);
			adims = matsize(I1bM);
			int nv = itos(gel(adims, 1));
			for (int j = 0; j < nv; j++) {
				string astr(GENtostr(gel(gcoeff(I1bM, j+1, 1),2)));
				badideals.push_back((badideal){ p, r, s, astr, nbad++ });
			}
		}
		if (p > badpmax) badpmax = p;
	}
	
	// read dlogs
	ifstream dlogfile(argv[3]);
	stringstream ss;
	getline(dlogfile, line);	// 0 column
	int n = 0;
	// first read bad ideal logs
	badpmax = 0;
	for (int j = 1; j < nbad; j++) {
		getline(dlogfile, line);
		string log = line.substr(line.find("bad ideals ") + 11);
		badideal bad1 = badideals[n++];
		dlogs.push_back((dlog){ bad1.s, bad1.p, bad1.r, log });
		if (bad1.p > badpmax) badpmax = bad1.p;
	}
	// read ideal logs
	uint64_t p, r; int s; string str;
	while (getline(dlogfile, line)) {
		line.erase(0, line.find(separator0) + 1); // remove line number
		str = line.substr(0, line.find(separator0));
		if (str == "SM") break;	// we have reached the Schirokauer maps
		ss << hex << str; ss >> p; ss.clear();
		line.erase(0, line.find(separator0) + 1);
		s = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);
		str = line.substr(0, line.find(separator0));
		ss << hex << str; ss >> r; ss.clear();
		line.erase(0, line.find(separator0) + 1);
		string log = line;
		dlogs.push_back((dlog){ s, p, r, log });
		n++;
	}
	// read Schirokauer map logs
	line.erase(0, line.find(separator0) + 1);
	s = atoi(line.substr(0, line.find(separator0)).c_str());
	line.erase(0, line.find(separator0) + 1);
	line.erase(0, line.find(separator0) + 1);
	sm.push_back(make_pair(s, line));
	while (getline(dlogfile, line)) {
		line.erase(0, line.find(separator0) + 1);
		line.erase(0, line.find(separator0) + 1);
		s = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);
		line.erase(0, line.find(separator0) + 1);
		sm.push_back(make_pair(s, line));
	}


	// read a,b,c, side
	int64_t inta = atoi(argv[4]);
	int64_t intb = atoi(argv[5]);
	int64_t intc = atoi(argv[6]);
	int side = atoi(argv[7]);
	
	GEN a = stoi(inta);
	GEN b = stoi(intb);
	GEN c = stoi(intc);

	GEN Ip; GEN bdims;
	if (side == 0) {
		I1a = idealhnf0(Ka, gadd(gadd(a, gmul(b, x)), gmul(c, gsqr(x))), NULL);
		I1a = idealmul(Ka, I1a, J2a);
		I1aM = idealfactor(Ka, I1a);
		adims = matsize(I1aM);
		int nv = itos(gel(adims,1));
		for (int j = 0; j < nv; ++j) {
			Ip = gcoeff(I1aM, j+1, 1);
			uint64_t p = itos(gel(Ip, 1));
			int v = itos(gcoeff(I1aM, j+1, 2));
			ideals.push_back(make_pair(p, v));
			if (is_bad(Ip, badideals)) {
				int k = find_bad(Ip, badideals);
				cout << side << " " << p << " " << badideals[k].r << " "
					<< v << " " << dlogs[k].log << endl;
			}
			else {
				uint64_t r = 0;
				uint64_t k = 0;
				bool found = find_log(Ip, Ka, J1a, x, dlogs, side, p, &r, &k);
				string log = "";
				if (found) log = dlogs[k].log;
				cout << side << " " << p << " " << r << " " << v << " " << log << endl;
			}
		}
	}
	else if (side == 1) {
		I1b = idealhnf0(Kb, gadd(gadd(a, gmul(b, y)), gmul(c, gsqr(y))), NULL);
		I1b = idealmul(Kb, I1b, J2b);
		I1bM = idealfactor(Kb, I1b);
		bdims = matsize(I1bM);
		int nv = itos(gel(bdims,1));
		for (int j = 0; j < nv; ++j) {
			Ip = gcoeff(I1bM, j+1, 1);
			uint64_t p = itos(gel(Ip, 1));
			int v = itos(gcoeff(I1bM, j+1, 2));
			ideals.push_back(make_pair(p, v));
			if (is_bad(Ip, badideals)) {
				int k = find_bad(Ip, badideals);
				cout << side << " " << p << " " << badideals[k].r << " "
					<< v << " " << dlogs[k].log << endl;
			}
			else {
				uint64_t r = 0;
				uint64_t k = 0;
				bool found = find_log(Ip, Kb, J1b, y, dlogs, side, p, &r, &k);
				string log = "";
				if (found) log = dlogs[k].log;
				cout << side << " " << p << " " << r << " " << v << " " << log << endl;
			}
		}
	}
}


