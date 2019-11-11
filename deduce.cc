#include <stdint.h>	// int64_t
#include <iostream> // cout
#include <gmpxx.h>
#include <math.h>	// sqrt
#include <fstream>	// file
#include <ctime>	// clock_t
#include <sstream>	// stringstream
#include <string>
#include <assert.h>
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
auto npos = std::string::npos;

typedef struct struct1 {
	int s;
	int64_t p;
	int64_t r;
	int e;
	string log;
} dlog;

vector<pair<int,string> > sm_exp;
vector<dlog> dlogs;
vector<pair<int,string> > sm;

int main(int argc, char** argv)
{
	if (argc != 3) {
		cout << "Usage:  ./deduce ell inputfile" << endl << endl;
		return 0;
	}

    string separator0 = " ";
    string separator1 = ":";
    string separator2 = ",";
	string line;
	ifstream file(argv[2]);

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

		string prime = line.substr(0, line.find(separator0));
		line.erase(0, line.find(separator0) + 1);

		int64_t r = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		int e = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		string val = line;

		if (prime == "sm_exp") {
			sm_exp.push_back(make_pair(side, val));
		}
		else {
			int64_t p = strtoull(prime.c_str(), NULL, 10);
			dlogs.push_back((dlog){ side, p, r, e, val });
			if (val == "") tside = side;
		}
	}

	// read side 1 ideal logs/Schirokauer exponents
	while (getline(file, line)) {
		if (line.find("Schirokauer Maps:", 0) == 0) break;

		int side = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		string prime = line.substr(0, line.find(separator0));
		line.erase(0, line.find(separator0) + 1);

		int64_t r = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		int e = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		string val = line;

		if (prime == "sm_exp") {
			sm_exp.push_back(make_pair(side, val));
		}
		else {
			int64_t p = strtoull(prime.c_str(), NULL, 10);
			dlogs.push_back((dlog){ side, p, r, e, val });
			if (val == "") tside = side;
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
	mpz_set_str(ell, argv[1], 10);
	int ndlogs = dlogs.size();
	mpz_set_ui(tlog, 0);
	// ideal logs first
	for (int i = 0; i < ndlogs; i++) {
		int side = dlogs[i].s;
		int64_t p = dlogs[i].p;
		int e = dlogs[i].e;
		string val = dlogs[i].log;
		if (val != "") {
			cout << " + " << e << "*" << val;
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
		cout << " + " << eval << "*" << val;
		mpz_set_str(exp, eval.c_str(), 10);
		mpz_set_str(log, val.c_str(), 10);
		mpz_mul(log, log, exp);
		mpz_add(tlog, tlog, log);
		mpz_mod(tlog, tlog, ell);	// reduce mod ell		
	}
	// Note that the following commented-out line is wrong - the J ideal has norm 1, therefore log zero.
	cout << " + 1";
	mpz_add_ui(tlog, tlog, 1);	// account for J ideal (warning - hardcoded to side 1)
	mpz_sub(tlog, ell, tlog);	// move known logs to other side

	cout << endl;
	cout << "vlog(" << q << ") = " << mpz_get_str(NULL, 10, tlog) << endl;

	mpz_clear(exp);
	mpz_clear(ell);
	mpz_clear(log);
	mpz_clear(tlog);

	return 0;
}
