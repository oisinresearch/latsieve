#include <cstdlib>
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

string deduce_on_one_side(string ellstr, string filename);
string deduce_from_both_sides(string ellstr, string filename, int target_side);

int main(int argc, char** argv)
{
	if (argc != 5) {
		//cout << "Usage:  ./deduce_full poly badidealfile dlogfile ell inputfile side" << endl << endl;
		cout << "Usage:  ./rotate_full dlogfile ell inputfile side" << endl << endl;
		return 0;
	}

	string dlog(argv[1]);
	string ellstr(argv[2]);
	string inputfile(argv[3]);
	string side(argv[4]);

	deduce_on_one_side(ellstr, inputfile);

	return 0;
}

string deduce_on_one_side(string ellstr, string filename)
{
	vector<pair<int,string> > sm_exp;
	vector<dlog> dlogs;
	vector<pair<int,string> > sm;

    string separator0 = " ";
    string separator1 = ":";
    string separator2 = ",";
	string line;
	ifstream file(filename);

	// read a b c
	getline(file, line);
	int64_t a = strtoll(line.substr(0, line.find(separator0)).c_str(), NULL, 10);
	line.erase(0, line.find(separator0) + 1);
	int64_t b = strtoll(line.substr(0, line.find(separator0)).c_str(), NULL, 10);
	line.erase(0, line.find(separator0) + 1);
	int64_t c = strtoll(line.substr(0, line.find(separator0)).c_str(), NULL, 10);
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

		int64_t r = strtoll(line.substr(0, line.find(separator0)).c_str(), NULL, 10);
		line.erase(0, line.find(separator0) + 1);

		int e = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		string val = line;

		if (pstr == "sm_exp") {
			sm_exp.push_back(make_pair(side, val));
		}
		else {
			int64_t p = strtoull(pstr.c_str(), NULL, 10);
			if (val == "") {
				val = deduce_from_both_sides(ellstr, "deduce_" + to_string(p) + ".txt", side);
				//cout << "\tvlog(" << p << ") = " << val << " on side " << side << endl;
				cout << "[a,b,c] = " << "[" << a << "," << b << "," << c << "]; p = " << p << endl;
			}
			dlogs.push_back((dlog){ side, p, r, e, val });
		}
		tside = side;
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
			//cout << " + " << e << "*" << val;
			mpz_set_str(log, val.c_str(), 10);
			mpz_mul_si(log, log, e);
			mpz_add(tlog, tlog, log);
			mpz_mod(tlog, tlog, ell);	// reduce mod ell
		}
		else { // can't hit this since val has always been looked up
			string pfilename = "deduce_" + to_string(p) + ".txt";
			val = deduce_from_both_sides(ellstr, pfilename, 1);
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
	if (tside == 1) {
		mpz_add_ui(tlog, tlog, 1);	// account for J ideal (warning - hardcoded to side 1)
		//cout << " + 1";
	}
	
	mpz_mod(tlog, tlog, ell);	// reduce mod ell

	if (tside == 0) {
		mpz_sub(tlog, ell, tlog);	// move known logs to other side
	}

	string tlogstr = mpz_get_str(NULL, 10, tlog);

	//cout << "[a,b,c] = " << "[" << a << "," << b << "," << c << "]; p = " << q << endl;

	cout << endl;
	cout << "vlog(" << q << ") = " << mpz_get_str(NULL, 10, tlog) << endl;

	mpz_clear(exp);
	mpz_clear(ell);
	mpz_clear(log);
	mpz_clear(tlog);

	return tlogstr;
}

string deduce_from_both_sides(string ellstr, string filename, int target_side)
{
	vector<pair<int,string> > sm_exp;
	vector<dlog> dlogs;
	vector<pair<int,string> > sm;

    string separator0 = " ";
    string separator1 = ":";
    string separator2 = ",";
	string line;
	ifstream file(filename);

	string pfilestr = filename.substr(7, filename.find(".") - 7);
	int64_t pfile = strtoull(pfilestr.c_str(), NULL, 10);

	// read a b c
	getline(file, line);
	int64_t a = strtoll(line.substr(0, line.find(separator0)).c_str(), NULL, 10);
	line.erase(0, line.find(separator0) + 1);
	int64_t b = strtoll(line.substr(0, line.find(separator0)).c_str(), NULL, 10);
	line.erase(0, line.find(separator0) + 1);
	int64_t c = strtoll(line.substr(0, line.find(separator0)).c_str(), NULL, 10);
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

		int64_t r = strtoll(line.substr(0, line.find(separator0)).c_str(), NULL, 10);
		line.erase(0, line.find(separator0) + 1);

		int e = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		string val = line;

		if (pstr == "sm_exp") {
			sm_exp.push_back(make_pair(side, val));
		}
		else {
			int64_t p = strtoull(pstr.c_str(), NULL, 10);
			if (val == "" && p != pfile) {
				val = deduce_from_both_sides(ellstr, "deduce_" + to_string(p) + ".txt", 0);
				//cout << "\tvlog(" << p << ") = " << val << " on side " << 0 << endl;
				cout << "[a,b,c] = " << "[" << a << "," << b << "," << c << "]; p = " << p << endl;
			}
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

		int64_t r = strtoll(line.substr(0, line.find(separator0)).c_str(), NULL, 10);
		line.erase(0, line.find(separator0) + 1);

		int e = atoi(line.substr(0, line.find(separator0)).c_str());
		line.erase(0, line.find(separator0) + 1);

		string val = line;

		if (pstr == "sm_exp") {
			sm_exp.push_back(make_pair(side, val));
		}
		else {
			int64_t p = strtoull(pstr.c_str(), NULL, 10);
			if (val == "" && p != pfile) {
				val = deduce_from_both_sides(ellstr, "deduce_" + to_string(p) + ".txt", 1);
				//cout << "\tvlog(" << p << ") = " << val << " on side " << 1 << endl;
				cout << "[a,b,c] = " << "[" << a << "," << b << "," << c << "]; p = " << p << endl;
			}
			dlogs.push_back((dlog){ side, p, r, e, val });
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
	//if (target_side == 1) {
	mpz_add_ui(tlog, tlog, 1);	// account for J ideal (warning - hardcoded to side 1)
	//cout << " + 1";
	//}

	//cout << endl << "\ttarget side = " << target_side << endl;

	//if (target_side == 0) {
	mpz_sub(tlog, ell, tlog);	// move known logs to other side
	//}

	string tlogstr = mpz_get_str(NULL, 10, tlog);

	//cout << endl;
	//cout << "vlog(" << q << ") = " << mpz_get_str(NULL, 10, tlog) << endl;

    if (q != 0)
    	cout << "[a,b,c] = " << "[" << a << "," << b << "," << c << "]; p = " << q << endl;

	mpz_clear(exp);
	mpz_clear(ell);
	mpz_clear(log);
	mpz_clear(tlog);

	return tlogstr;
}

