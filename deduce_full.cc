#include <stdint.h>	// int64_t
#include <iostream> // cout
#include <gmpxx.h>
#include <math.h>	// sqrt
#include <fstream>	// file
#include <ctime>	// clock_t
#include <sstream>	// stringstream
#include <string>
#include <assert.h>

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
	int64_t p;
	int64_t r;
	int e;
	string log;
} dlog;

typedef struct struct2 {
	int64_t p;
	int64_t r;
	int s;
} prime;

vector<pair<int,string> > sm_exp;
vector<dlog> dlogs;
vector<pair<int,string> > sm;

string deduce_on_one_side(string ellstr, string filename);
string deduce_from_both_sides(string ellstr, string filename);
void getabc(string filename, int64_t &a, int64_t &b, int64_t &c);
int num_unknown_primes_both_sides(string filename);
int num_unknown_primes_one_side(string filename);

int main(int argc, char** argv)
{
	if (argc != 7) {
		cout << "Usage:  ./deduce_full poly badidealfile dlogfile ell inputfile side" << endl << endl;
		return 0;
	}

	string poly(argv[1]);
	string bad(argv[2]);
	string dlog(argv[3]);
	string ellstr(argv[4]);
	string inputfile(argv[5]);
	string side(argv[6]);

	deduce_on_one_side(ellstr, inputfile);

	return 0;
}

string deduce_on_one_side(string ellstr, string filename)
{
    string separator0 = " ";
    string separator1 = ":";
    string separator2 = ",";
	string line;
	ifstream file(filename);

	// read a b c
	getline(file, line);
	int64_t a = atoi(line.substr(0, line.find(separator0)).c_str());
	line.erase(0, line.find(separator0) + 1);
	int64_t b = atoi(line.substr(0, line.find(separator0)).c_str());
	line.erase(0, line.find(separator0) + 1);
	int64_t c = atoi(line.substr(0, line.find(separator0)).c_str());
	line.erase(0, line.find(separator0) + 1);

	int tside = 0;

	// read "side S:"
	getline(file, line);

	// read side 1 ideal logs/Schirokauer exponents
	while (getline(file, line)) {
		if (line.find("Schirokauer Maps:", 0) == 0) break;

		int side = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		string pstr = line.substr(0, line.find(separator0));
		line.erase(0, line.find(separator0) + 1);

		int64_t r = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		int e = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		string val = line;

		if (pstr == "sm_exp") {
			sm_exp.push_back(make_pair(side, val));
		}
		else {
			int64_t p = strtoull(pstr.c_str(), NULL, 10);
			if (val == "")
				val = deduce_from_both_sides(ellstr, "deduce_" + to_string(p) + ".txt");
			dlogs.push_back((dlog){ side, p, r, e, val });
		}
	}

	// read Schirokauer map logs
	while (getline(file, line)) {
		line.erase(0, line.find(separator0) + 1);

		int side = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		uint64_t r = atoi(line.substr(0, line.find(separator0)).c_str());
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
			//cout << " + " << e << "*" << val;
			mpz_set_str(log, val.c_str(), 10);
			mpz_mul_si(log, log, e);
			mpz_add(tlog, tlog, log);
			mpz_mod(tlog, tlog, ell);	// reduce mod ell
		}
		else { // can't hit this since val has always been looked up
			string pfilename = "deduce_" + to_string(p) + ".txt";
			val = deduce_from_both_sides(ell, pfilename);
		}
	}
	// then Schirokauer maps
	int nsm_exp = sm_exp.size();
	int offset = 0;
	while (sm[offset].first != dlogs[0].s) offset++;
	for (int i = 0; i < nsm_exp; i++) {
		int side = sm[offset + i].first;
		string eval = sm_exp[i].second;
		string val = sm[offset + i].second;
		//cout << " + " << eval << "*" << val;
		mpz_set_str(exp, eval.c_str(), 10);
		mpz_set_str(log, val.c_str(), 10);
		mpz_mul(log, log, exp);
		mpz_add(tlog, tlog, log);
		mpz_mod(tlog, tlog, ell);	// reduce mod ell		
	}
	// Note that the following commented-out line is wrong - the J ideal has norm 1, therefore log zero.
	mpz_add_ui(tlog, tlog, 1);	// account for J ideal (warning - hardcoded to side 1)
	//cout << " + 1";
	//mpz_sub(tlog, ell, tlog);	// move known logs to other side
	mpz_mod(tlog, tlog, ell);	// reduce mod ell

	string tlogstr = mpz_get_str(NULL, 10, tlog);

	//cout << endl;
	cout << "vlog(" << q << ") = " << mpz_get_str(NULL, 10, tlog) << endl;

	mpz_clear(exp);
	mpz_clear(ell);
	mpz_clear(log);
	mpz_clear(tlog);

	return tlogstr;
}

string deduce_from_both_sides(string ellstr, string filename)
{
    string separator0 = " ";
    string separator1 = ":";
    string separator2 = ",";
	string line;
	ifstream file(filename);

	// read a b c
	getline(file, line);
	int64_t a = atoi(line.substr(0, line.find(separator0)).c_str());
	line.erase(0, line.find(separator0) + 1);
	int64_t b = atoi(line.substr(0, line.find(separator0)).c_str());
	line.erase(0, line.find(separator0) + 1);
	int64_t c = atoi(line.substr(0, line.find(separator0)).c_str());
	line.erase(0, line.find(separator0) + 1);

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

		int64_t r = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		int e = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		string val = line;

		if (pstr == "sm_exp") {
			sm_exp.push_back(make_pair(side, val));
		}
		else {
			int64_t p = strtoull(pstr.c_str(), NULL, 10);
			if (val == "")
				val = deduce_from_both_sides(ellstr, "deduce_" + to_string(p) + ".txt");
			dlogs.push_back((dlog){ side, p, r, e, val });
		}
	}

	// read side 1 ideal logs/Schirokauer exponents
	while (getline(file, line)) {
		if (line.find("Schirokauer Maps:", 0) == 0) break;

		int side = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		string pstr = line.substr(0, line.find(separator0));
		line.erase(0, line.find(separator0) + 1);

		int64_t r = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		int e = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		string val = line;

		if (pstr == "sm_exp") {
			sm_exp.push_back(make_pair(side, val));
		}
		else {
			int64_t p = strtoull(pstr.c_str(), NULL, 10);
			if (val == "")
				val = deduce_from_both_sides(ellstr, "deduce_" + to_string(p) + ".txt");
			dlogs.push_back((dlog){ side, p, r, e, val });
		}
	}

	// read Schirokauer map logs
	while (getline(file, line)) {
		line.erase(0, line.find(separator0) + 1);

		int side = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		uint64_t r = atoi(line.substr(0, line.find(separator0)).c_str());
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
			//cout << " + " << e << "*" << val;
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
		//cout << " + " << eval << "*" << val;
		mpz_set_str(exp, eval.c_str(), 10);
		mpz_set_str(log, val.c_str(), 10);
		mpz_mul(log, log, exp);
		mpz_add(tlog, tlog, log);
		mpz_mod(tlog, tlog, ell);	// reduce mod ell		
	}
	// Note that the following commented-out line is wrong - the J ideal has norm 1, therefore log zero.
	//cout << " + 1";
	mpz_add_ui(tlog, tlog, 1);	// account for J ideal (warning - hardcoded to side 1)
	mpz_sub(tlog, ell, tlog);	// move known logs to other side

	string tlogstr = mpz_get_str(NULL, 10, tlog);

	//cout << endl;
	//cout << "vlog(" << q << ") = " << mpz_get_str(NULL, 10, tlog) << endl;

	mpz_clear(exp);
	mpz_clear(ell);
	mpz_clear(log);
	mpz_clear(tlog);

	return tlogstr;
}

int num_unknown_primes_both_sides(string filename)
{
	vector<prime> unknownprimes;

    string separator0 = " ";
    string separator1 = ":";
    string separator2 = ",";
	string line;
	ifstream file(filename);

	// read a b c
	getline(file, line);
	int64_t a = atoi(line.substr(0, line.find(separator0)).c_str());
	line.erase(0, line.find(separator0) + 1);
	int64_t b = atoi(line.substr(0, line.find(separator0)).c_str());
	line.erase(0, line.find(separator0) + 1);
	int64_t c = atoi(line.substr(0, line.find(separator0)).c_str());
	line.erase(0, line.find(separator0) + 1);

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

		int64_t r = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		int e = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		string val = line;

		if (pstr != "sm_exp") {
			int64_t p = strtoull(pstr.c_str(), NULL, 10);
			if (val == "") unknownprimes.push_back((prime){ p, r, 0 });
		}
	}

	// read side 1 ideal logs/Schirokauer exponents
	while (getline(file, line)) {
		if (line.find("Schirokauer Maps:", 0) == 0) break;

		int side = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		string pstr = line.substr(0, line.find(separator0));
		line.erase(0, line.find(separator0) + 1);

		int64_t r = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		int e = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		string val = line;

		if (prime != "sm_exp") {
			int64_t p = strtoull(pstr.c_str(), NULL, 10);
			if (val == "") unknownprimes.push_back((prime){ p, r, 1 });
		}
	}

	return unknownprimes->size();
}

int num_unknown_primes_one_side(string filename)
{
	vector<prime> unknownprimes;

    string separator0 = " ";
    string separator1 = ":";
    string separator2 = ",";
	string line;
	ifstream file(filename);

	// read a b c
	getline(file, line);
	int64_t a = atoi(line.substr(0, line.find(separator0)).c_str());
	line.erase(0, line.find(separator0) + 1);
	int64_t b = atoi(line.substr(0, line.find(separator0)).c_str());
	line.erase(0, line.find(separator0) + 1);
	int64_t c = atoi(line.substr(0, line.find(separator0)).c_str());
	line.erase(0, line.find(separator0) + 1);

	int tside = 0;

	// read "side S:"
	getline(file, line);

	// read side 1 ideal logs/Schirokauer exponents
	while (getline(file, line)) {
		if (line.find("Schirokauer Maps:", 0) == 0) break;

		int side = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		string pstr = line.substr(0, line.find(separator0));
		line.erase(0, line.find(separator0) + 1);

		int64_t r = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		int e = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		string val = line;

		if (prime != "sm_exp") {
			int64_t p = strtoull(pstr.c_str(), NULL, 10);
			if (val == "") unknownprimes.push_back((prime){ p, r, 1 });
		}
	}

	return unknownprimes->size();
}
