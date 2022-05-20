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
	GEN ideal;
	string log;
} dlog;

typedef struct struct2 {
	vector<string>* lines;
	struct pari_thread pth;
	int workid;
} thread_root_arg;

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
	if (argc != 11) {
		cout << "Usage:  ./descend_tnfs polyfile isofile dlogfile num_Sch_0 num_Sch_1 a b c d side" << endl << endl;
		cout << "Output:  Lines of the form" << endl;
		cout << "side prime root exponent log" << endl;
		cout << "followed by global Schirokauer maps" << endl;
		cout << "If side is 0 or 1, descend just that side, if side is 2 descend both sides"
			 << endl << endl;
		return 0;
	}

	// initialize pari library
	pari_init(2000000000,65536);

	string polyfile = string(argv[1]);
	string isofile = string(argv[2]);
	string dlogfile = string(argv[3]);
	string loadnfs3dstr = "\\r ~/nfs3d.gp";
	string loadpolystr = "loadpolytnfs(\"" + polyfile + "\");";

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

	// read dlogfile
	int numsch0 = atoi(argv[4]);
	int numsch1 = atoi(argv[5]);
	GEN mix0 = readstr(dlogfile.c_str());
	long numlogs = glength(mix0) - numsch0 - numsch1;
	GEN mix = cgetg(numlogs+1, t_VEC);
	GEN v0 = strsplit(gel(mix0, 1), strtoGENstr(" "));
	GEN log0 = geval(gel(v0, 2));
	string row1str = "[0,[0]~,0," + string(GENtostr(log0)) + "]";
	gel(mix, 1) = geval(strtoGENstr(row1str.c_str()));
	GEN v1, logi, idsi, pi, idi, si, logrow;
	for (long i = 2; i <= numlogs; ++i) {
		v1 = strsplit(gel(mix0, i), strtoGENstr(" "));
		logi = geval(gel(v1, 2));
		idsi = geval(gconcat1(vecslice0(v1, 3, glength(v1))));
		pi = gcopy(gel(gel(idsi, 1), 1));
		idi = gcopy(gel(gel(idsi, 1), 2));
		si = gcopy(gel(idsi, 2));
		logrow = cgetg(5, t_VEC);
		gel(logrow, 1) = pi;
		gel(logrow, 2) = idi;
		gel(logrow, 3) = si;
		gel(logrow, 4) = logi;
		gel(mix, i) = logrow;
		//if (i % 1000 == 0) cout << i << endl;
	}

	// sort logs and split into two lists:  searchable ideals and corresponding logs
	//cout << "number of logs is " << numlogs << endl;
	//cout << "number of Schirokauer maps on side 0 is " << numsch0 << endl;
	//cout << "number of Schirokauer maps on side 1 is " << numsch1 << endl;
	mix = vecsort0(mix, NULL, 0);
	/*for (int i = 1; i <= 10; i++) {
		output(gel(mix, i));
	}*/

	GEN ideals = cgetg(numlogs+1, t_VEC);
	GEN logs = cgetg(numlogs+1, t_VEC);

	for (long i = 1; i <= numlogs; ++i) {
		GEN pi = gcopy(gel(gel(mix, i),1));
		GEN idi = gcopy(gel(gel(mix, i),2));
		GEN si = gcopy(gel(gel(mix, i),3));
		GEN logi = gcopy(gel(gel(mix, i),4));
		GEN idealrow = cgetg(4, t_VEC);
		gel(idealrow, 1) = pi;
		gel(idealrow, 2) = idi;
		gel(idealrow, 3) = si;
		gel(ideals, i) = idealrow;
		gel(logs, i) = logi;
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
	GEN eps0 = gsubgs(gpowgs(ell, 4), 1);
	GEN eps1 = gsubgs(gpowgs(ell, 4), 1);
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

	// read a,b,c, side
	int64_t inta = atoi(argv[6]);
	int64_t intb = atoi(argv[7]);
	int64_t intc = atoi(argv[8]);
	int64_t intd = atoi(argv[9]);
	int side = atoi(argv[10]);
	
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
	cout << "J 1 " << GENtostr(log0) << endl;
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
			}
			else {
				cout << GENtostr(sj) << " " << GENtostr(pj) << " " << GENtostr(idj) << " "
					<< v << " " << endl;
			}
		}
	}
	GEN A1 = liftall(gmul(m1, gmodulo(gadd(gadd(a, gmul(b, y3)), gmul(gadd(c, gmul(d, y3)), x3)), f)));
	GEN feps = gdiv(liftall(gsubgs(gpow(gmodulo(gmodulo(A1, f), gsqr(ell)), eps0, prec), 1)), ell);
	for (int i = 0; i < numsch0; i++) {
		cout << "0 sm_exp " << i+1 << " 1 " << GENtostr(polcoef(feps, i+1, -1)) << endl;
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
			}
			else {
				cout << GENtostr(sj) << " " << GENtostr(pj) << " " << GENtostr(idj) << " "
					<< v << " " << endl;
			}
		}
	}
	GEN A2 = liftall(gmul(m2, gmodulo(gadd(gadd(a, gmul(b, y4)), gmul(gadd(c, gmul(d, y4)), x4)), g)));
	GEN geps = gdiv(liftall(gsubgs(gpow(gmodulo(gmodulo(A2, g), gsqr(ell)), eps1, prec), 1)), ell);
	for (int i = 0; i < numsch1; i++) {
		cout << "1 sm_exp " << i+1 << " 1 " << GENtostr(polcoef(geps, i+1, -1)) << endl;
	}
	cout << "Schirokauer Maps:" << endl;
	for (int i = 0; i < numsch0 + numsch1; i++) {
		v1 = strsplit(gel(mix0, numlogs+1+i), strtoGENstr(" "));
		cout << GENtostr(geval(gel(v1, 4))) << " SM " << GENtostr(geval(gel(v1, 3))) << " "
			<< GENtostr(geval(gel(v1, 2))) << endl;
	}
}


